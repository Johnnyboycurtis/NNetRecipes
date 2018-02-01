/******************************************************************************/
/*                                                                            */
/*  ANNEAL - Minimize a function using simulated annealing                    */
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

#define RANDMAX 32767

static void shake ( int nvars , double *center , double *x , double temp ) ;

double anneal (
      double (*func)( int n , double *x ) ,  // Function to be minimized
      int ntemps ,        // Number of temperatures to use
      int niters ,        // Number of iterations at each temperature
      int setback ,       // Set back iters this much if improved
      int maxcalls ,      // Limit number of function calls (due to setback)
      double fquit ,      // Quit if function reduced to this amount
      double starttemp ,  // Starting temperature (standard deviation)
      double stoptemp ,   // Stopping temperature
      int nvars ,         // Number of variables in x
      double *x ,         // Set to starting estimate when called; returns best
      double *cent        // Work vector nvars long
      )
{
   int i, iter, seed, bestseed, improved, itemp ;
   double tempmult, temp, fval, bestfval ;

/*
   For every temperature, the center around which we will perturb is the
   best point so far.  This is kept in 'cent', so initialize cent to the
   user's starting estimate.   Also, initialize 'bestfval', the best function
   value so far, to be the function value at that starting point.
*/

   for (i=0 ; i<nvars ; i++)
      cent[i] = x[i] ;

   bestfval = (*func) ( nvars , x ) ;


/*
   This is the temperature reduction loop and the iteration within temperature
   loop.  We use a slick trick to keep track of the best point at a given
   temperature.  We certainly don't want to replace 'cent' every time an
   improvement is had, as then we would be moving our center about, compromising
   the global nature of the algorithm.  We could, of course, have a second work
   area in which we save the 'best so far for this temperature' point.  But if
   there are a lot of variables, the usual case, this wastes memory.  What we
   do is to save the seed of the random number generator which created the
   improvement.  Then later, when we need to retrieve the best, simply set
   the random seed and regenerate it.  This technique also saves a lot of
   copying time if many improvements are made for a single temperature.
*/

   temp = starttemp ;
   tempmult = exp ( log ( stoptemp / starttemp ) / (ntemps-1) ) ;

   for (itemp=0 ; itemp<ntemps ; itemp++) {   // Temp reduction loop

      improved = 0 ;                          // Flags if this temp improved

      for (iter=0 ; iter<niters ; iter++) {   // Iters per temp loop

         seed = rand () ;                     // Get a random seed
         srand ( seed ) ;                     // Brute force set it
         shake ( nvars , cent , x , temp ) ;  // Randomly perturb about cent
         fval = (*func) ( nvars , x ) ;       // evaluate function there
         if (fval < bestfval) {               // If this iteration improved
            bestfval = fval ;                 // then update the best so far
            bestseed = seed ;                 // and save seed to recreate it
            improved = 1 ;                    // Flag that we improved
            if (bestfval <= fquit)            // If we reached the user's
               break ;                        // limit, we can quit
            iter -= setback ;                 // It often pays to keep going
            if (iter < 0)                     // at this temperature if we
               iter = 0 ;                     // are still improving
            }
         if (! --maxcalls)                    // This is for safety in case
            break ;                           // setback causes excess looping
         }

      if (improved) {                         // If this temp saw improvement
         srand ( bestseed ) ;                 // set seed to what caused it
         shake ( nvars , cent , x , temp ) ;  // and recreate that point
         for (i=0 ; i<nvars ; i++)            // which will become center
            cent[i] = x[i] ;                  // for next temperature iters
         srand ( bestseed / 2 + 999 ) ;       // Jog seed away from best
         }

      if (bestfval <= fquit)                  // If we reached the user's
         break ;                              // limit, we can quit

      if (! maxcalls)                         // This is for safety in case
         break ;                              // setback causes excess looping

      temp *= tempmult ;                      // Reduce temp for next pass
      }                                       // through this temperature loop

   for (i=0 ; i<nvars ; i++)                  // Return the best point in x
      x[i] = cent[i] ;

   return bestfval ;
}

/*
--------------------------------------------------------------------------------

   SHAKE - Randomly perturb the point

--------------------------------------------------------------------------------
*/

static void shake ( int nvars , double *center , double *x , double temp )
{
   double r ;

// Recall that the variance of a uniform deviate on 0-1 is 1/12.
// Adding four such random variables multiplies the variance by 4,
// while dividing by 2 divides the variance by 4.

   temp *= 3.464101615 / (2. * (double) RANDMAX ) ;  // SQRT ( 12 ) = 3.464...

   while (nvars--) {
      r = (double) rand() + (double) rand() -
          (double) rand() - (double) rand() ;
      *x++ = *center++ + temp * r ;
      }
}
