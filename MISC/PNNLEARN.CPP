/******************************************************************************/
/*                                                                            */
/*  PNNLEARN - Use jackknife method to find best sigma                        */
/*                                                                            */
/* Copyright (c) 1993 by Academic Press, Inc.                                 */
/*                                                                            */
/* All rights reserved.  Permission is hereby granted, until further notice,  */
/* to make copies of this diskette, which are not for resale, provided these  */
/* copies are made from this master diskette only, and provided that the      */
/* following copyright notice appears on the diskette label:                  */
/* (c) 1993 by Academic Press, Inc.                                           */
/*                                                                            */
/* Except as previously stated, no part of the computer program embodied in   */
/* this diskette may be reproduced or transmitted in any form or by any means,*/
/* electronic or mechanical, including input into storage in any information  */
/* system for resale, without permission in writing from the publisher.       */
/*                                                                            */
/* Produced in the United States of America.                                  */
/*                                                                            */
/* ISBN 0-12-479041-0                                                         */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

static int evaluate ( int nvars , int npops , int *ntrain , double **tsets ,    // npops vector of ptrs to training sets
                      int model , double sigma , double *outputs ) ;

extern int pnn ( int nvars , int npops , int *n , double **t ,
                 int model , double sigma , double *u , double *o ) ;

int pnnlearn (
      int nvars ,         // Number of variables in x
      int npops ,         // Number of populations
      int *ntrain ,       // npops vector of number of samples of each pop
      double **tsets ,    // npops vector of ptrs to training sets
      int model ,         // Distance/Weight model (see comments, normally 1)
      int nsigs ,         // Number of sigmas for initial search
      int maxits ,        // Max number of iterations
      double sigmin ,     // Minimum possible value to try for sigma
      double sigmax ,     // And max
      double *sigma ,     // Returned best scale parameter
      double *outputs     // Work vector npops long
      )
{
   int i, score, ibest, left_score, mid_score, right_score, prev, iter ;
   int try_score ;
   double sig, sigmult, left, mid, right, itry ;

   sigmult = exp ( log ( sigmax / sigmin ) / (nsigs-1) ) ;
   sig = sigmin ;

/*
   Evaluate confusion for nsigs sigmas ranging from sigmin to sigmax, with
   their values spaced by equal ratios (sigmult).  Keep track of the best
   sigma in 'mid' and its subscript so we know if it is the first or last.
   Also keep track of the best score and its neighbors.
*/

   prev = 0 ;        // Not needed.  Shuts up compilers over use before init
   mid_score = -1 ;  // Keep track of best here

   for (i=0 ; i<nsigs ; i++) {
      score = evaluate ( nvars , npops , ntrain , tsets , model , sig ,
                         outputs ) ;
      if (score > mid_score) {
         mid = sig ;
         ibest = i ;
         mid_score = score ;
         left_score = prev ;
         }
      else if (i == (ibest+1))
         right_score = score ;
      prev = score ;
      sig *= sigmult ;
      }

   left = mid / sigmult ;
   right = mid * sigmult ;

/*
   At this point we (probably) have three sigmas (left, mid, right) with their
   scores (left_score, mid_score, right_score) such that the mid has the
   highest score.

   There is another possibility.  If ibest=0 it means the best score was had at
   the extreme low end of the user's sigma range, and left_score is unknown.
   Bail out the silly user by shrinking until either we get a three point tie
   (in which case we quit, choosing the center sigma as the best) or until
   our score drops, in which case we just go on to the refinement phase.

   Also, we have the same possibility at the high end.
*/

   if (! ibest) {   // Best at low end, so keep going down
      for (;;) {
         left_score = evaluate ( nvars , npops , ntrain , tsets , model , left ,
                                 outputs ) ;
         if ((left_score == mid_score)  &&  (right_score == mid_score))
            goto DONE ;  // Totally flat, so quit trying
         if (left_score < mid_score) // Trio with mid best, so go refine
            break ;
         right = mid ;               // Still getting better, so move on down
         right_score = mid_score ;
         mid = left ;
         mid_score = left_score ;
         left = mid / sigmult ;
         }
      }

   else if (ibest == (nsigs-1)) {   // Best at high end, so keep going up
      for (;;) {
         right_score = evaluate( nvars , npops , ntrain , tsets , model , right,
                                 outputs ) ;
         if ((left_score == mid_score)  &&  (right_score == mid_score))
            goto DONE ;  // Totally flat, so quit trying
         if (right_score < mid_score) // Trio with mid best, so go refine
            break ;
         left = mid ;               // Still getting better, so move on up
         left_score = mid_score ;
         mid = right ;
         mid_score = right_score ;
         right = mid * sigmult ;
         }
      }

/*
   At this point we definitely have three sigmas (left, mid, right) with their
   scores (left_score, mid_score, right_score) such that the mid has the
   highest score.  Refine.

*/

   for (iter=0 ; iter<maxits ; iter++) {

      if (mid / left  >  right / mid)  {  // Left interval larger, so split it
         itry = sqrt ( mid * left ) ;     // Split ratio in half
         try_score = evaluate ( nvars , npops , ntrain , tsets , model , itry ,
                                outputs ) ;
         if ((try_score > mid_score)  ||
             ((try_score == mid_score)  &&  (left_score > right_score))) {
            right = mid ;
            right_score = mid_score ;
            mid = itry ;
            mid_score = try_score ;
            }
         else {
            left = itry ;
            left_score = try_score ;
            }
         }

      else {  // Right interval larger, so split it
         itry = sqrt ( mid * right ) ;      // Split ratio in half
         try_score = evaluate ( nvars , npops , ntrain , tsets , model , itry ,
                                outputs ) ;
         if ((try_score > mid_score)  ||
             ((try_score == mid_score)  &&  (left_score < right_score))) {
            left = mid ;
            left_score = mid_score ;
            mid = itry ;
            mid_score = try_score ;
            }
         else {
            right = itry ;
            right_score = try_score ;
            }
         }

      if ((left_score == mid_score)  &&  (right_score == mid_score))
         break ;  // Totally flat, so quit trying
      }

DONE:
   *sigma = mid ;
   return mid_score ;
}


/*
--------------------------------------------------------------------------------

  evaluate - Local routine to evaluate performance for a given sigma

--------------------------------------------------------------------------------
*/

static int evaluate (
      int nvars ,         // Number of variables in x
      int npops ,         // Number of populations
      int *ntrain ,       // npops vector of number of samples of each pop
      double **tsets ,    // npops vector of ptrs to training sets
      int model ,         // Distance/Weight model (see comments, normally 1)
      double sigma ,     // Scale parameter
      double *outputs     // Work vector npops long
      )
{
   int i, tclass, ipop, score, popscore, n, exclude ;
   double *tptr, temp, *v1ptr, *v2ptr ;

   score = 0 ;
   for (ipop=0 ; ipop<npops ; ipop++) {   // Evaluate for each population
      popscore = 0 ;                      // Score each pop separately

      n = ntrain[ipop] ;                  // Cases from this population
      tptr = tsets[ipop] ;                // Point to them
      --ntrain[ipop] ;                    // Temporarily exclude test case

      exclude = n ;                       // Test with this case
      while (--exclude >= 0) {            // Exclude each case
         if (exclude < n-1) {             // Do not need to swap last case!
            v1ptr = tptr + exclude * nvars ;
            v2ptr = tptr + (n-1) * nvars ;
            for (i=0 ; i<nvars ; i++) {
               temp = *v1ptr ;
               *v1ptr++ = *v2ptr ;
               *v2ptr++ = temp ;
               }
            }
         tclass = pnn ( nvars , npops , ntrain , tsets , model , sigma ,
                        tptr+(n-1)*nvars , outputs ) ;
         if (tclass == ipop)
            ++popscore ;
         } // For all excludes

      ++ntrain[ipop] ;                    // Restore to correct value
      score += popscore ;                 // Cumulate overall score
      } // For all pops
   return score ;
}
