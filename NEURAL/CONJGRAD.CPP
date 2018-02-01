/******************************************************************************/
/*                                                                            */
/*  conjgrad - Conjugate gradient learning                                    */
/*                                                                            */
/*  Normally this returns the mean square error, which will be 0-1.           */
/*  If the user interrupted, it returns the negative mean square error.       */
/*  Insufficient memory returns -2.                                           */
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

double LayerNet::conjgrad (
   TrainingSet *tptr , // Training set to use
   int maxits ,        // Maximum iterations allowed
   double reltol ,     // Relative error change tolerance
   double errtol       // Quit if error drops this low
   )
{
   int i, n, iter, key, retry, max_retry ;
   double gam, *g, *h, *outdelta, *hid2delta, *grad, *base ;
   double error ;
   double prev_err ;
                     
   max_retry = 5 ;

/*
   Allocate work memory
*/

   MEMTEXT ( "CONJGRAD work" ) ;
   if (nhid2) {
      hid2delta = (double *) MALLOC ( nhid2 * sizeof(double) ) ;
      if (hid2delta == NULL)
         return -2.0 ;
      }
   else
      hid2delta = NULL ;

   outdelta = (double *) MALLOC ( nout * sizeof(double) ) ;

   if (nhid1 == 0)               // No hidden layer
      n = nout * (nin+1) ;
   else if (nhid2 == 0)          // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
   else                          // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;

   grad = (double *) MALLOC ( n * sizeof(double) ) ;
   base = (double *) MALLOC ( n * sizeof(double) ) ;
   g = (double *) MALLOC ( n * sizeof(double) ) ;
   h = (double *) MALLOC ( n * sizeof(double) ) ;

   if ((outdelta == NULL) || (grad == NULL) ||
       (base == NULL) || (g == NULL) || (h == NULL)) {
      if (hid2delta != NULL)
         FREE ( hid2delta ) ;
      if (outdelta != NULL)
         FREE ( outdelta ) ;
      if (grad != NULL)
         FREE ( grad ) ;
      if (base != NULL)
         FREE ( base ) ;
      if (g != NULL)
         FREE ( g ) ;
      if (h != NULL)
         FREE ( h ) ;
      return -2.0 ;   // Flags error
      }

   prev_err = 1.e30 ;
   error = find_grad ( tptr , hid2delta , outdelta , grad ) ;

   memcpy ( g , grad , n * sizeof(double) ) ;
   memcpy ( h , grad , n * sizeof(double) ) ;

/*
   Main iteration loop is here
*/

   for (iter=0 ; iter<maxits ; iter++) {  // Each iter is an epoch

/*
   Check current error against user's max.  Abort if user pressed ESCape
*/

      if (error <= errtol)   // If our error is within user's limit
         break ;             // then we are done!

      if (error <= reltol)   // Generally not necessary: reltol<errtol in
         break ;             // practice, but help silly users

      if (kbhit()) {         // Was a key pressed?
         key = getch () ;    // Read it if so
         while (kbhit())     // Flush key buffer in case function key
            getch () ;       // or key was held down
         if (key == 27) {    // ESCape
            error = -error ; // Flags user that ESCape was pressed
            break ;
            }
         }

      prev_err = error ;
      error = direcmin ( tptr , error , 10 , 1.e-10 ,
                         0.5 , base , grad ) ;
      if (error < 0.0)  // Indicates user pressed ESCape
         goto CGFINISH ;

      if ((2.0 * (prev_err - error)) <=       // If this direc gave poor result
          (reltol * (prev_err + error + 1.e-10))) { // will use random direc
         prev_err = error ;                   // But first exhaust grad
         error = find_grad ( tptr , hid2delta , outdelta , grad ) ;
         error = direcmin ( tptr , error , 15 , 1.e-10 ,
                            1.e-3 , base , grad ) ;
         for (retry=0 ; retry<max_retry ; retry++) {
            for (i=0 ; i<n ; i++)
               grad[i] = (double) (rand() - RANDMAX/2) / (RANDMAX * 10.0) ;
            error = direcmin ( tptr , error , 10 , 1.e-10 ,
                               1.e-2 , base , grad ) ;
            if (error < 0.0)  // Indicates user pressed ESCape
               goto CGFINISH ;
            if (retry < max_retry/2)
               continue ;
            if ((2.0 * (prev_err - error)) >
                (reltol * (prev_err + error + 1.e-10)))
               break ;   // Get out of retry loop if we improved enough
            } // For retry
         if (retry == max_retry)   // If we exhausted all tries
            break ;                // probably hopeless
         memcpy ( g , grad , n * sizeof(double) ) ;
         memcpy ( h , grad , n * sizeof(double) ) ;
         } // If this dir gave poor result

      prev_err = error ;

/*
   Setup for next iteration
*/

      error = find_grad ( tptr , hid2delta , outdelta , grad ) ;
      gam = gamma ( g , grad ) ;
      if (gam < 0.0)
         gam = 0.0 ;
      if (gam > 1.0)
         gam = 1.0 ;

      find_new_dir ( gam , g , h , grad ) ;
      }  // This is the end of the main iteration loop

/*
   Free work memory
*/

CGFINISH:
   MEMTEXT ( "CONJGRAD work" ) ;
   if (hid2delta != NULL)
      FREE ( hid2delta ) ;
   FREE ( outdelta ) ;
   FREE ( grad ) ;
   FREE ( base ) ;
   FREE ( g ) ;
   FREE ( h ) ;

   return error ;
}

/*
--------------------------------------------------------------------------------

   Local routine to compute gradient for a trial epoch

--------------------------------------------------------------------------------
*/

double LayerNet::find_grad (
   TrainingSet *tptr ,
   double *hid2delta ,
   double *outdelta ,
   double *grad
   )
{
   int i, j, size, tset, tclass, n, nprev, nnext ;
   double error, *dptr, diff, delta, *hid1grad, *hid2grad, *outgrad ;
   double *outprev, *prevact, *nextcoefs, *nextdelta, *gradptr ;

/*
   Compute size of each training sample
*/

   if (outmod == OUTMOD_CLASSIFY)
      size = nin + 1 ;
   else if (outmod == OUTMOD_AUTO)
      size = nin ;
   else if (outmod == OUTMOD_GENERAL)
      size = nin + nout ;

/*
   Compute length of grad vector and gradient positions in it.
   Also point to layer previous to output and its size.
   Ditto for layer after hid1.
*/

   if (nhid1 == 0) {      // No hidden layer
      n = nout * (nin+1) ;
      outgrad = grad ;
      nprev = nin ;
      }
   else if (nhid2 == 0) { // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
      hid1grad = grad ;
      outgrad = grad + nhid1 * (nin+1) ;
      outprev = hid1 ;
      nprev = nhid1 ;
      nnext = nout ;
      nextcoefs = out_coefs ;
      nextdelta = outdelta ;
      }
   else {                 // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;
      hid1grad = grad ;
      hid2grad = grad + nhid1 * (nin+1) ;
      outgrad = hid2grad + nhid2 * (nhid1+1) ;
      outprev = hid2 ;
      nprev = nhid2 ;
      nnext = nhid2 ;
      nextcoefs = hid2_coefs ;
      nextdelta = hid2delta ;
      }

   for (i=0 ; i<n ; i++)  // Zero gradient for summing
      grad[i] = 0.0 ;

   error = 0.0 ;  // Will cumulate total error here
   for (tset=0 ; tset<tptr->ntrain ; tset++) { // Do all samples

      dptr = tptr->data + size * tset ;     // Point to this sample
      trial ( dptr ) ;                      // Evaluate network for it

      if (outmod == OUTMOD_AUTO) {          // If this is AUTOASSOCIATIVE
         for (i=0 ; i<nout ; i++) {         // then the expected outputs
            diff = *dptr++ - out[i] ;       // are just the inputs
            error += diff * diff ;
            outdelta[i] = diff * actderiv ( out[i] ) ;
            }
         }

      else if (outmod == OUTMOD_CLASSIFY) {  // If this is Classification
         tclass = (int) dptr[nin] - 1 ;     // class is stored after inputs
         for (i=0 ; i<nout ; i++) {         // Recall that train added a
            if (tclass == i)                // fraction so that the above
               diff = NEURON_ON - out[i] ;  // truncation to get tclass is
            else                            // always safe in any radix
               diff = NEURON_OFF - out[i] ;
            error += diff * diff ;
            outdelta[i] = diff * actderiv ( out[i] ) ;
            }
         }

      else if (outmod == OUTMOD_GENERAL) {  // If this is GENERAL output
         dptr += nin ;                      // outputs stored after inputs
         for (i=0 ; i<nout ; i++) {
            diff = *dptr++ - out[i] ;
            error += diff * diff ;
            outdelta[i] = diff * actderiv ( out[i] ) ;
            }
         }

/*
   Cumulate output gradient
*/

      if (nhid1 == 0)         // No hidden layer
         prevact = tptr->data + size * tset ;
      else
         prevact = outprev ;  // Point to previous layer
      gradptr = outgrad ;
      for (i=0 ; i<nout ; i++) {
         delta = outdelta[i] ;
         for (j=0 ; j<nprev ; j++)
            *gradptr++ += delta * prevact[j] ;
         *gradptr++ += delta ;   // Bias activation is always 1
         }

/*
   Cumulate hid2 gradient (if it exists)
*/
   
      if (nhid2) {
         gradptr = hid2grad ;
         for (i=0 ; i<nhid2 ; i++) {
            delta = 0.0 ;
            for (j=0 ; j<nout ; j++)
               delta += outdelta[j] * out_coefs[j*(nhid2+1)+i] ;
            delta *= actderiv ( hid2[i] ) ;
            hid2delta[i] = delta ;
            for (j=0 ; j<nhid1 ; j++)
               *gradptr++ += delta * hid1[j] ;
            *gradptr++ += delta ;   // Bias activation is always 1
            }
         }

/*
   Cumulate hid1 gradient (if it exists)
*/
   
      if (nhid1) {
         prevact = tptr->data + size * tset ;
         gradptr = hid1grad ;
         for (i=0 ; i<nhid1 ; i++) {
            delta = 0.0 ;
            for (j=0 ; j<nnext ; j++)
               delta += nextdelta[j] * nextcoefs[j*(nhid1+1)+i] ;
            delta *= actderiv ( hid1[i] ) ;
            for (j=0 ; j<nin ; j++)
               *gradptr++ += delta * prevact[j] ;
            *gradptr++ += delta ;   // Bias activation is always 1
            }
         }

      } // for all tsets
   
   return error / ((double) tptr->ntrain * (double) nout) ;
}

/*
--------------------------------------------------------------------------------

   Local routine to find gamma

--------------------------------------------------------------------------------
*/

double LayerNet::gamma ( double *g , double *grad )
{
   int i, n ;
   double denom, numer ;

   if (nhid1 == 0)        // No hidden layer
      n = nout * (nin+1) ;
   else if (nhid2 == 0)   // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
   else                   // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;

   numer = denom = 0. ;

   for (i=0 ; i<n ; i++) {
      denom += g[i] * g[i] ;
      numer += (grad[i] - g[i]) * grad[i] ;  // Grad is neg gradient
      }

   if (denom == 0.)   // Should never happen (means gradient is zero!)
      return 0. ;
   else
      return numer / denom ;
}

/*
--------------------------------------------------------------------------------

   Local routine to find correction for next iteration

--------------------------------------------------------------------------------
*/

void LayerNet::find_new_dir ( double gam , double *g ,
                              double *h , double *grad )
{
   int i, n ;

   if (nhid1 == 0)        // No hidden layer
      n = nout * (nin+1) ;
   else if (nhid2 == 0)   // One hidden layer
      n = nhid1 * (nin+1) + nout * (nhid1+1) ;
   else                   // Two hidden layers
      n = nhid1 * (nin+1) + nhid2 * (nhid1+1) + nout * (nhid2+1) ;

   for (i=0 ; i<n ; i++) {
      g[i] = grad[i] ;
      grad[i] = h[i] = g[i] + gam * h[i] ;
      }
}


/*
--------------------------------------------------------------------------------

   Local routine for debugging

--------------------------------------------------------------------------------
*/

void LayerNet::check_grad ( TrainingSet *tptr , double *grad )
{
   int i, j, n ;
   double f0, f1, deriv, dot, len1, len2 ;

   dot = len1 = len2 = 0.0 ;
   f0 = trial_error ( tptr ) ;

   for (i=0 ; i<nhid1 ; i++) {
      for (j=0 ; j<=nin ; j++) {
         hid1_coefs[i*(nin+1)+j] += .001 ;
         f1 = trial_error ( tptr ) ;
         hid1_coefs[i*(nin+1)+j] -= .001 ;
         deriv = 10000.0 * (f0 - f1) ;
         len1 += *grad * *grad ;
         len2 += deriv * deriv ;
         dot += *grad++ * deriv ;
         }
      }

   for (i=0 ; i<nhid2 ; i++) {
      for (j=0 ; j<=nhid1 ; j++) {
         hid2_coefs[i*(nhid1+1)+j] += .001 ;
         f1 = trial_error ( tptr ) ;
         hid2_coefs[i*(nhid1+1)+j] -= .001 ;
         deriv = 10000.0 * (f0 - f1) ;
         len1 += *grad * *grad ;
         len2 += deriv * deriv ;
         dot += *grad++ * deriv ;
         }
      }

   if (nhid1 == 0)        // No hidden layer
      n = nin ;
   else if (nhid2 == 0)   // One hidden layer
      n = nhid1 ;
   else                   // Two hidden layers
      n = nhid2 ;
   for (i=0 ; i<nout ; i++) {
      for (j=0 ; j<=n ; j++) {
         out_coefs[i*(n+1)+j] += .001 ;
         f1 = trial_error ( tptr ) ;
         out_coefs[i*(n+1)+j] -= .001 ;
         deriv = 10000.0 * (f0 - f1) ;
         len1 += *grad * *grad ;
         len2 += deriv * deriv ;
         dot += *grad++ * deriv ;
         }
      }
}
