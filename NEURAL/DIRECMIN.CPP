/******************************************************************************/
/*                                                                            */
/*  direcmin - Minimize along a direction                                     */
/*                                                                            */
/*  Normally this returns the mean square error, which will be 0-1.           */
/*  If the user interrupted, it returns the negative mean square error.       */
/*                                                                            */
/*  This is a two step process.  First we find three points whose center has  */
/*  the smallest function value (we bound the minimum).                       */
/*  Then we use Brent's algorithm to refine the interval.                     */
/*                                                                            */
/*  We leave 'coefs' set at the point that produced the minimum and return    */
/*  the error function at that point.  We change the direction 'dir' to       */
/*  be the actual distance moved.                                             */
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

double LayerNet::direcmin (
   TrainingSet *tptr , // Training set to use
   double start_err ,  // Error (function value) at starting coefficients
   int itmax ,         // Upper limit on number of iterations allowed
   double eps ,        // Small, but greater than machine precision
   double tol ,        // Brent's tolerance (>= sqrt machine precision)
   double *base ,      // Work area (stepping out point)
   double *direc )     // Work area (stepping out direction)
{
   int key, user_quit, iter ;
   double step, x1, x2, x3, t1, t2, numer, denom, max_step ;
   double xlow, xhigh, xbest, testdist ;
   double current_err, err, previous_err, step_err ;
   double prevdist, frecent, fthirdbest, fsecbest, fbest ;
   double tol1, tol2, xrecent, xthirdbest, xsecbest, xmid;
   double  first_step = 2.5 ; // Heuristically found best

   user_quit = 0 ;

/*
   Take one step out in the gradient direction.  First preserve
   original weights for use as departure point parameterized by STEP.
*/

   preserve ( base ) ;   // Establishes a base for stepping out
   step_out ( first_step , direc , base ) ;
   err = trial_error ( tptr ) ;

/*
   If it increased, we had numerical problems computing the direction or
   the direction itself is too large a step.
   Negate the direction and use -1, 0 and 1.618 as first three steps.
   Otherwise use 0, 1 and 2.618 as first three steps.
*/

   if (err > start_err) {
      negate_dir ( direc ) ;
      x1 = -first_step ;
      x2 = 0. ;
      previous_err = err ;
      current_err = start_err ;
      }
   else {
      x1 = 0. ;
      x2 = first_step ;
      previous_err = start_err ;
      current_err = err ;
      }

/*
   At this point we have taken a single step and the function decreased.
   Take one more step in the golden ratio.
   Also keep errors lined up as 'previous_err', 'current_err' and 'err'.
   The corresponding abscissae will be x1, x2 and x3.
*/

   if (kbhit()) {          // Was a key pressed?
      key = getch () ;     // Read it if so
      while (kbhit())      // Flush key buffer in case function key
         getch () ;        // or key was held down
      if (key == 27)       // ESCape
         return (- err) ;
      }

   x3 = x2 + 1.618034 * first_step ;
   step_out ( x3 , direc , base ) ;
   err = trial_error ( tptr ) ;

/*
   We now have three points x1, x2 and x3 with corresponding errors
   of 'previous_err', 'current_err' and 'err'.
   Endlessly loop until we bracket the minimum with the outer two.
*/

   while (err < current_err) {   // As long as we are descending...

      if (kbhit()) {          // Was a key pressed?
         key = getch () ;     // Read it if so
         while (kbhit())      // Flush key buffer in case function key
            getch () ;        // or key was held down
         if (key == 27) {     // ESCape
            user_quit = 1 ;
            break ;
            }
         }

/*
   Try a parabolic fit to estimate the location of the minimum.
*/

      t1 = (x2 - x1) * (current_err - err) ;
      t2 = (x2 - x3) * (current_err - previous_err) ;
      denom = 2. * ( t2 - t1 ) ;
      if (fabs ( denom ) < eps) {
         if (denom > 0.)
            denom = eps ;
         else
            denom = -eps ;
         }
      step = x2 + ((x2 - x1) * t1  -  (x2 - x3) * t2) / denom ;//Here if perfect
      max_step = x2 + 200. * (x3 - x2) ; // Don't jump too far

      if ((x2 - step) * (step - x3)  >  0.) {         // It's between x2 and x3
         step_out ( step , direc , base ) ;
         step_err = trial_error ( tptr ) ;

         if (step_err < err) {   // It worked!  We found min between b and c.

            x1 = x2 ;
            x2 = step ;
            previous_err = current_err ;
            current_err = step_err ;
            goto BOUNDED ;
            }
         else if (step_err > current_err) { // Slight miscalc.  Min at x2.
            x3 = step ;
            err = step_err ;
            goto BOUNDED ;
            }
         else {             // Parabolic fit was total waste.  Use default.
            step = x3 + 1.618034 * (x3 - x2) ;
            step_out ( step , direc , base ) ;
            step_err = trial_error ( tptr ) ;
            }
         }

      else if ((x3 - step) * (step - max_step) > 0.0) { // Between x3 and lim
         step_out ( step , direc , base ) ;
         step_err = trial_error ( tptr ) ;
         if (step_err < err) {  // Decreased, so advance by golden ratio
            x2 = x3 ;
            x3 = step ;
            step = x3 + 1.618034 * (x3 - x2) ;
            current_err = err ;
            err = step_err ;
            step_out ( step , direc , base ) ;
            step_err = trial_error ( tptr ) ;
            }
         }

      else if ((step - max_step) * (max_step - x3)  >= 0.) {  // Beyond limit
         step = max_step ;
         step_out ( step , direc , base ) ;
         step_err = trial_error ( tptr ) ;
         if (step_err < err) {  // Decreased, so advance by golden ratio
            x2 = x3 ;
            x3 = step ;
            step = x3 + 1.618034 * (x3 - x2) ;
            current_err = err ;
            err = step_err ;
            step_out ( step , direc , base ) ;
            step_err = trial_error ( tptr ) ;
            }
         }

      else {  // Wild!  Reject parabolic and use golden ratio.
         step = x3 + 1.618034 * (x3 - x2) ;
         step_out ( step , direc , base ) ;
         step_err = trial_error ( tptr ) ;
         }

/*
   Shift three points and continue endless loop
*/

      x1 = x2 ;
      x2 = x3 ;
      x3 = step ;
      previous_err = current_err ;
      current_err = err ;
      err = step_err ;
      } // Endless stepping out loop

BOUNDED:
   step_out ( x2 , direc , base);//Leave coefs at min

   if (x1 > x3) {  // We may have switched direction at start.
      t1 = x1 ;    // Brent's method which follows assumes ordered parameter.
      x1 = x3 ;
      x3 = t1 ;
      }

   if (user_quit) {
      update_dir ( x2 , direc ) ;// Make it be the actual dist moved
      return -current_err ;
      }

/*
--------------------------------------------------------------------------------

   At this point we have bounded the minimum between x1 and x3.

   Go to the refinement stage.  We use Brent's algorithm.

--------------------------------------------------------------------------------
*/

/*
  Initialize prevdist, the distance moved on the previous step, to 0 so that
  the 'if (fabs ( prevdist )  >  tol1)' encountered on the first iteration
  below will fail, forcing a golden section the first time.  Also initialize
  step to 0 to avoid a zealous compiler from pointing out that it was
  referenced before being set.
*/

   prevdist = step = 0.0 ;

/*
   We always keep the minimum bracketed between xlow and xhigh.
   xbest has the min function so far (or latest if tie).
   xsecbest and xthirdbest are the second and third best.
*/

   xbest = xsecbest = xthirdbest = x2 ;
   xlow = x1 ;
   xhigh = x3 ;

   fbest = fsecbest = fthirdbest = current_err ;


/*
   Main loop.  For safety we impose a limit on iterations.
*/

   for (iter=0 ; iter<itmax ; iter++) {

      xmid = 0.5 * (xlow + xhigh) ;
      tol1 = tol * (fabs ( xbest ) + eps) ;
      tol2 = 2. * tol1 ;

      if (kbhit()) {          // Was a key pressed?
         key = getch () ;     // Read it if so
         while (kbhit())      // Flush key buffer in case function key
            getch () ;        // or key was held down
         if (key == 27) {     // ESCape
            user_quit = 1 ;
            break ;
            }
         }

/*
   The following convergence test simultaneously makes sure xhigh and
   xlow are close relative to tol2, and that xbest is near the midpoint.
*/
      if (fabs ( xbest - xmid )  <=  (tol2 - 0.5 * (xhigh - xlow)))
         break ;

      if (fabs ( prevdist )  >  tol1) {  // If we moved far enough try parabolic fit
         t1 = (xbest - xsecbest) * (fbest - fthirdbest) ; // Temps for the
         t2 = (xbest - xthirdbest) * (fbest - fsecbest) ; // parabolic estimate
         numer = (xbest - xthirdbest) * t2  -  (xbest - xsecbest) * t1 ;
         denom = 2. * (t1 - t2) ;  // Estimate will be numer / denom
         testdist = prevdist ;     // Will soon verify interval is shrinking
         prevdist = step ;         // Save for next iteration
         if (denom != 0.0)         // Avoid dividing by zero
            step = numer / denom ; // This is the parabolic estimate to min
         else 
            step = 1.e30 ;         // Assures failure of next test

         if ((fabs ( step ) < fabs ( 0.5 * testdist ))// If shrinking
          && (step + xbest > xlow)             // and within known bounds
          && (step + xbest < xhigh)) {         // then we can use the
            xrecent = xbest + step ;           // parabolic estimate
            if ((xrecent - xlow  <  tol2)  ||  // If we are very close
                (xhigh - xrecent <  tol2)) {   // to known bounds
               if (xbest < xmid)               // then stabilize
                  step = tol1 ;
               else
                  step = -tol1 ;
               }
            }
         else {  // Parabolic estimate poor, so use golden section
            prevdist = (xbest >= xmid)  ?  xlow - xbest  :  xhigh - xbest ;          // Poor so use    
            step = .3819660 * prevdist ;
            }
         }
      else { // prevdist did not exceed tol1: we did not move far enough
             // to justify a parabolic fit.  Use golden section.
         prevdist = (xbest >= xmid)  ?  xlow - xbest  :  xhigh - xbest ;
         step = .3819660 * prevdist ;
         }

      if (fabs (step)  >=  tol1)     // In order to numerically justify
         xrecent = xbest + step ;    // another trial we must move a
      else {                         // decent distance.
         if (step > 0.)
            xrecent = xbest + tol1 ;
         else
            xrecent = xbest - tol1 ;
         }

/*
   At long last we have a trial point 'xrecent'.  Evaluate the function.
*/

      step_out ( xrecent , direc , base ) ;
      frecent = trial_error ( tptr ) ;

      if (frecent <= fbest) {    // If we improved...
         if (xrecent >= xbest)   // Shrink the (xlow,xhigh) interval by
            xlow = xbest ;       // replacing the appropriate endpoint
         else
            xhigh = xbest ;
         xthirdbest = xsecbest ; // Update x and f values for best,
         xsecbest = xbest ;      // second and third best
         xbest = xrecent ;
         fthirdbest = fsecbest ;
         fsecbest = fbest ;
         fbest = frecent ;
         }

      else {                  // We did not improve
         if (xrecent < xbest) // Shrink the (xlow,xhigh) interval by
            xlow = xrecent ;  // replacing the appropriate endpoint
         else
            xhigh = xrecent ;

         if ((frecent <= fsecbest)   // If we at least beat the second best
          || (xsecbest == xbest)) {  // or we had a duplication
            xthirdbest = xsecbest ;  // we can update the second and third
            xsecbest = xrecent ;     // best, though not the best.
            fthirdbest = fsecbest ;  // Recall that we started iters with
            fsecbest = frecent ;     // best, sec and third all equal.
            }
         else if ((frecent <= fthirdbest) // Oh well.  Maybe at least we can
           || (xthirdbest == xbest)       // beat the third best or rid
           || (xthirdbest == xsecbest)) { // ourselves of a duplication
            xthirdbest = xrecent ;        // (which is how we start the
            fthirdbest = frecent ;        // iterations)
            }
         }
      }
 
   step_out ( xbest , direc , base );//Leave coefs at min
   update_dir ( xbest , direc ) ;// Make it be the actual distance moved

   if (user_quit)
      return -fbest ;
   else 
      return fbest ;
}


/*
--------------------------------------------------------------------------------

   Local routine to preserve coefs in 'base'

--------------------------------------------------------------------------------
*/

void LayerNet::preserve ( double *base )
{
   int n ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      memcpy ( base , out_coefs , n * sizeof(double) ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      memcpy ( base , hid1_coefs , n * sizeof(double) ) ;
      base += n ;
      n = nout * (nhid1+1) ;
      memcpy ( base , out_coefs , n * sizeof(double) ) ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      memcpy ( base , hid1_coefs , n * sizeof(double) ) ;
      base += n ;
      n = nhid2 * (nhid1+1) ;
      memcpy ( base , hid2_coefs , n * sizeof(double) ) ;
      base += n ;
      n = nout * (nhid2+1) ;
      memcpy ( base , out_coefs , n * sizeof(double) ) ;
      }
}


/*
--------------------------------------------------------------------------------

   Local routine to step out from base

--------------------------------------------------------------------------------
*/

void LayerNet::step_out ( double step , double *direc , double *base )
{
   int i, n ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      for (i=0 ; i<n ; i++)
         out_coefs[i] = *(base++) + *(direc++) * step ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      for (i=0 ; i<n ; i++)
         hid1_coefs[i] = *(base++) + *(direc++) * step ;
      n = nout * (nhid1+1) ;
      for (i=0 ; i<n ; i++)
         out_coefs[i] = *(base++) + *(direc++) * step ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      for (i=0 ; i<n ; i++)
         hid1_coefs[i] = *(base++) + *(direc++) * step ;
      n = nhid2 * (nhid1+1) ;
      for (i=0 ; i<n ; i++)
         hid2_coefs[i] = *(base++) + *(direc++) * step ;
      n = nout * (nhid2+1) ;
      for (i=0 ; i<n ; i++)
         out_coefs[i] = *(base++) + *(direc++) * step ;
      }
}

/*
--------------------------------------------------------------------------------

   Local routine to make 'dir' be the actual distance moved

--------------------------------------------------------------------------------
*/

void LayerNet::update_dir ( double step , double *direc )
{
   int n ;

  if (nhid1 == 0)               // No hidden layer
      n = nout * (nin+1) ;
   else if (nhid2 == 0)          // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
   else                          // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;

   while (n--) {
      *direc *= step ;
      ++direc ;
      }
}

/*
--------------------------------------------------------------------------------

   Local routine to negate 'dir'

--------------------------------------------------------------------------------
*/

void LayerNet::negate_dir ( double *direc )
{
   int n ;

  if (nhid1 == 0)               // No hidden layer
      n = nout * (nin+1) ;
   else if (nhid2 == 0)          // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
   else                          // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;

   while (n--) {
      *direc = -*direc ;
      ++direc ;
      }
}
