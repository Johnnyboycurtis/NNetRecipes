/******************************************************************************/
/*                                                                            */
/*  PNN - Classify using probabilistic neural network                         */
/*                                                                            */
/*  This uses a probabilistic neural network model to classify an unknown.    */
/*  Training set i (i=0,...,npop-1) contains ntrain[i] samples, each of which */
/*  is pointed to by tsets[i].                                                */
/*  Each sample is a vector containing nvars variables.                       */
/*  The model parameter controls which of several choices of weight functions */
/*  and distance measures are used.  The usual Gaussian distribution is       */
/*  model 1.  Model 2 uses the sum of the absolute differences rather than    */
/*  the sum of squared differences (squared Euclidean distance) as the        */
/*  argument to the exponential function.  Model 3 uses 1 / (1 + d*d), where  */
/*  d is the Euclidean distance, as the weighting function.                   */
/*  Sigma is the scale parameter.                                             */
/*  It returns the number (0 to npops-1) of the most likely population to     */
/*  which the unknown belongs.  It also computes the outputs vector so the    */
/*  caller can apply prior probabilities or costs.                            */
/*                                                                            */
/*  Note that the input sigma is scaled by the number of variables.           */
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

int pnn (
      int nvars ,         // Number of variables in x
      int npops ,         // Number of populations
      int *ntrain ,       // npops vector of number of samples of each pop
      double **tsets ,    // npops vector of ptrs to training sets
      int model ,         // Distance/Weight model (see comments, normally 1)
      double sigma ,      // Scale parameter
      double *unknown ,   // nvars vector to be classified
      double *outputs     // npops vector of classification functions values
      )
{
   int i, n, ipop, icase, ibest ;
   double diff, dist, *tptr, *optr, best ;

   sigma *= (double) nvars ;              // Keep sigma meaningful
   best = -1. ;                           // Keep track of best mean function

   for (ipop=0 ; ipop<npops ; ipop++) {   // Evaluate for each population

      n = ntrain[ipop] ;                  // Number of samples of this pop
      tptr = tsets[ipop] ;                // Point to them
      optr = outputs+ipop ;               // Cumulate mean weight function here
      *optr = 0.0 ;

      for (icase=0 ; icase<n ; icase++) { // Do all cases in this population

         dist = 0.0 ;                     // Will sum distance measure here

         if (model == 2) {                // Use sum of absolute differences
            for (i=0 ; i<nvars ; i++) {   // Compute the distance measure
               diff = unknown[i] - *tptr++ ;
               dist += fabs ( diff ) ;
               }
            dist /= sigma ;
            }

         else {                           // Use squared Euclidean distance
            for (i=0 ; i<nvars ; i++) {   // Compute the distance measure
               diff = unknown[i] - *tptr++ ;
               dist += diff * diff ;
               }
            dist /= sigma * sigma ;
            }

         if (model < 3)                   // Use exponential weighting
            *optr += exp ( - dist ) ;
         else 
            *optr += 1. / (1. + dist ) ;
         }

      *optr /= (double) n ;  // Mean function value for this population

      if (*optr > best) {    // Keep track of best for classification
         best = *optr ;
         ibest = ipop ;
         }
      }

   return ibest ;
}  
