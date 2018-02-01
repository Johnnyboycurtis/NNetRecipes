/******************************************************************************/
/*                                                                            */
/*  ANNEAL - Use simulated annealing to optimize LayerNet weights             */
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
#include "const.h"       // System and limitation constants, typedefs, structs
#include "classes.h"     // Includes all class headers
#include "funcdefs.h"    // Function prototypes

void LayerNet::anneal (
   TrainingSet *tptr ,        // Training set to use
   struct LearnParams *lptr , // User's general learning parameters
   LayerNet *bestnet ,        // Work area used to keep best network
   int init                   // Use zero suffix (initialization) anneal parms?
   )
{
   int ntemps, niters, setback, reg, nvars, key, user_quit ;
   int iter, improved, ever_improved, itemp ;
   long seed, bestseed ;
   char msg[80] ;
   double tempmult, temp, fval, bestfval, starttemp, stoptemp, fquit ;
   SingularValueDecomp *sptr ;
   struct AnnealParams *aptr ; // User's annealing parameters
                             
   aptr = lptr->ap ;

/*
   The parameter 'init' is nonzero if we are initializing
   weights for learning.  If zero we are attempting to break
   out of a local minimum.  The main effect  of this parameter
   is whether or not we use the zero suffix variables in the
   anneal parameters.
   A second effect is that regression is used only for
   initialization, not for escape.
*/

   if (init) {
      ntemps = aptr->temps0 ;
      niters = aptr->iters0 ;
      setback = aptr->setback0 ;
      starttemp = aptr->start0 ;
      stoptemp = aptr->stop0 ;
      }
   else {
      ntemps = aptr->temps ;
      niters = aptr->iters ;
      setback = aptr->setback ;
      starttemp = aptr->start ;
      stoptemp = aptr->stop ;
      }

/*
   Initialize other local parameters.  Note that there is no sense using
   regression if there are no hidden layers.  Also, regression is almost
   always counterproductive for local minimum escape.
*/

   fquit = lptr->quit_err ;
   reg = init  &&  nhid1  &&  (lptr->init != 1) ;

/*
   Allocate the singular value decomposition object for REGRESS.
   Also allocate a work area for REGRESS to preserve matrix.
*/

   if (reg) {
      if (nhid1 == 0)         // No hidden layer
         nvars = nin + 1 ;
      else if (nhid2 == 0)    // One hidden layer
         nvars = nhid1 + 1 ;
      else                    // Two hidden layers
         nvars = nhid2 + 1 ;

      MEMTEXT ( "ANNEAL: new SingularValueDecomp" ) ;
      sptr = new SingularValueDecomp ( tptr->ntrain , nvars , 1 ) ;

      if ((sptr == NULL)  || ! sptr->ok) {
         memory_message (
            "for annealing with regression. Try ANNEAL NOREGRESS.");
         if (sptr != NULL)
            delete sptr ;
         neterr = 1.0 ; // Flag failure to LayerNet::learn which called us
         return ;
         }
      }

/*
   For every temperature, the center around which we will perturb is the
   best point so far.  This is kept in 'bestnet', so initialize it to the
   user's starting estimate.   Also, initialize 'bestfval', the best
   function value so far, to be the function value at that starting point.
*/

   copy_weights ( bestnet , this ) ; // Current weights are best so far
   if (init)
      bestfval = 1.e30 ;  // Force it to accept SOMETHING
   else 
      bestfval = trial_error ( tptr ) ;

/*
   This is the temperature reduction loop and the iteration within
   temperature loop.  We use a slick trick to keep track of the
   best point at a given temperature.  We certainly don't want to
   replace the best every time an improvement is had, as then we
   would be moving our center about, compromising the global nature
   of the algorithm.  We could, of course, have a second work area
   in which we save the 'best so far for this temperature' point.
   But if there are a lot of variables, the usual case, this wastes
   memory.  What we do is to save the seed of the random number
   generator which created the improvement.  Then later, when we
   need to retrieve the best, simply set the random seed and
   regenerate it.  This technique also saves a lot of copying time
   if many improvements are made for a single temperature.
*/

   temp = starttemp ;
   tempmult = exp( log( stoptemp / starttemp ) / (ntemps-1)) ;
   ever_improved = 0 ;                       // Flags if improved at all
   user_quit = 0 ;                           // Flags user pressed ESCape

   for (itemp=0 ; itemp<ntemps ; itemp++) {  // Temp reduction loop

      improved = 0 ;                         // Flags if this temp improved

      if (init) {
         sprintf ( msg , "\nANNEAL temp=%.2lf ", temp ) ;
         progress_message ( msg ) ;
         }

      for (iter=0 ; iter<niters ; iter++) {  // Iters per temp loop

         seed = longrand () ;                // Get a random seed
         slongrand ( seed ) ;                // Brute force set it
         perturb (bestnet, this, temp, reg) ;// Randomly perturb about best

         if (reg)                            // If using regression, estimate
            fval = regress ( tptr , sptr ) ; // out weights now
         else                                // Otherwise just evaluate
            fval = trial_error ( tptr ) ;

         if (fval < bestfval) {              // If this iteration improved
            bestfval = fval ;                // then update the best so far
            bestseed = seed ;                // and save seed to recreate it
            ever_improved = improved = 1 ;   // Flag that we improved

            if (bestfval <= fquit)           // If we reached the user's
               break ;                       // limit, we can quit

            iter -= setback ;                // It often pays to keep going
            if (iter < 0)                    // at this temperature if we
               iter = 0 ;                    // are still improving
            }
         }                                   // Loop: for all iters at a temp

      if (improved) {                        // If this temp saw improvement
         slongrand ( bestseed ) ;            // set seed to what caused it
         perturb (bestnet, this, temp, reg) ;// and recreate that point
         copy_weights ( bestnet , this ) ;   // which will become next center
         slongrand ( bestseed / 2 + 999 ) ;  // Jog seed away from best

         if (init) {
            sprintf ( msg , " err=%.3lf%% ", 100.0 * bestfval ) ;
            progress_message ( msg ) ;
            }
         }

      if (bestfval <= fquit)  // If we reached the user's
         break ;              // limit, we can quit

      if (kbhit()) {          // Was a key pressed?
         key = getch () ;     // Read it if so
         while (kbhit())      // Flush key buffer in case function key
            getch () ;        // or key was held down
         if (key == 27) {     // ESCape
            user_quit = 1 ;   // Flags user that ESCape was pressed
            break ;
            }
         }

      if (user_quit)
         break ;

      temp *= tempmult ;      // Reduce temp for next pass
      }                       // through this temperature loop


/*
   The trials left this weight set and neterr in random condition.
   Make them equal to the best, which will be the original
   if we never improved.

   Also, if we improved and are using regression, recall that bestnet
   only contains the best hidden weights, as we did not bother to run
   regress when we updated bestnet.  Do that now before returning.
*/

   copy_weights ( this , bestnet ) ; // Return best weights in this net
   neterr = bestfval ;               // Trials destroyed weights, err

   if (ever_improved  &&  reg)
      neterr = regress ( tptr , sptr ) ; // regressed output weights

   if (reg) {
      MEMTEXT ( "ANNEAL: delete SingularValueDecomp" ) ;
      delete sptr ;
      }
}


/*
--------------------------------------------------------------------------------

   Local routine to perturb coefficients

--------------------------------------------------------------------------------
*/

void LayerNet::perturb ( LayerNet *cent , LayerNet *perturbed , double temp ,
                         int reg )
{
   int n ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      shake ( n , cent->out_coefs , perturbed->out_coefs , temp ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      shake ( n , cent->hid1_coefs , perturbed->hid1_coefs , temp ) ;
      if (! reg) {
         n = nout * (nhid1+1) ;
         shake ( n , cent->out_coefs , perturbed->out_coefs , temp ) ;
         }
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      shake ( n , cent->hid1_coefs , perturbed->hid1_coefs , temp ) ;
      n = nhid2 * (nhid1+1) ;
      shake ( n , cent->hid2_coefs , perturbed->hid2_coefs , temp ) ;
      if (! reg) {
         n = nout * (nhid2+1) ;
         shake ( n , cent->out_coefs , perturbed->out_coefs , temp ) ;
         }
      }
}
