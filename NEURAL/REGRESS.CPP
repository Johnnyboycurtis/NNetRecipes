/******************************************************************************/
/*                                                                            */
/*  REGRESS - Use regression to compute LayerNet output weights               */
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

double LayerNet::regress (
   TrainingSet *tptr ,       // Training set used for regression input
   SingularValueDecomp *sptr // Work areas and needed functions
   )

{
   int i, out, tset, size, nvars ;
   double *aptr, *bptr, *dptr, err, temp, diff ;


/*
   Compute the size of each training sample in tptr->data and the number
   of independent variables (columns of matrix)
*/

   if (outmod == OUTMOD_CLASSIFY)
      size = nin + 1 ;
   else if (outmod == OUTMOD_AUTO)
      size = nin ;
   else if (outmod == OUTMOD_GENERAL)
      size = nin + nout ;

   if (nhid1 == 0)         // No hidden layer
      nvars = nin + 1 ;
   else if (nhid2 == 0)    // One hidden layer
      nvars = nhid1 + 1 ;
   else                    // Two hidden layers
      nvars = nhid2 + 1 ;

/*
   Pass through training set, building matrix, then find its singular value
   decomposition.  We keep a copy of it so we can compute the error later.
*/

   aptr = sptr->a ;                     // Will build matrix here

   for (tset=0 ; tset<tptr->ntrain ; tset++) { // Do all training samples

      dptr = tptr->data + size * tset ; // Point to this sample

      if (nhid1 == 0) {                 // No hidden layer
         for (i=0 ; i<nin ; i++)        // so matrix is just inputs
            *aptr++ = *dptr++ ;
         }

      else if (nhid2 == 0) {            // One hidden layer
         for (i=0 ; i<nhid1 ; i++)      // so matrix is hidden1 activations
            *aptr++ = activity ( dptr , hid1_coefs+i*(nin+1) , nin ) ;
         }

      else {                            // Two hidden layers
         for (i=0 ; i<nhid1 ; i++)
            hid1[i] = activity ( dptr , hid1_coefs+i*(nin+1) , nin ) ;
         for (i=0 ; i<nhid2 ; i++)
            *aptr++ = activity ( hid1 , hid2_coefs+i*(nhid1+1) , nhid1 ) ;
         }
      *aptr++ = 1.0 ;  // Bias term is last column of matrix
      } // For each training sample


/*
   Do the singular value decomposition.
   Then solve for weights for each output neuron.
   After each output weight vector is computed (using backsub),
   compute the activation of that output neuron, compare it to
   its desired value in the training set, and cumulate the error.
*/

   sptr->svdcmp () ;

   err = 0.0 ;

   for (out=0 ; out<nout ; out++) {  // For each output neuron

      bptr = sptr->b ;               // Backsub routine wants RHS here

      for (tset=0 ; tset<tptr->ntrain ; tset++) {

         dptr = tptr->data + size * tset ;    // Training sample starts here

         if (outmod == OUTMOD_AUTO) {         // If this is AUTOASSOCIATIVE
            temp = dptr[out] ;                // output is just input
            if (temp > 0.999999)              // Avoid problems in
               temp = 0.999999 ;              // inverse_act function
            if (temp < 0.000001)
               temp = 0.000001 ;
            *bptr++ = inverse_act ( temp ) ;  // Inverse activation function
            }

         else if (outmod == OUTMOD_CLASSIFY) { // If this is Classification
            if ((int) dptr[nin] == out+1)     // class identifier past inputs
               *bptr++ = inverse_act ( NEURON_ON ) ; // Inverse of NEURON_ON
            else 
               *bptr++ = inverse_act ( NEURON_OFF ) ;
            }

         else if (outmod == OUTMOD_GENERAL) { // If this is GENERAL output
            temp = dptr[nin+out] ;            // output is just past input
            if (temp > 0.999999)
               temp = 0.999999 ;
            if (temp < 0.000001)
               temp = 0.000001 ;
            *bptr++ = inverse_act ( temp ) ;    // Inverse activation function
            }
         } // For all training samples

      bptr = out_coefs + out * nvars ;    // Weight vector for this output will
      sptr->backsub ( 1.e-8 , bptr ) ;    // go here.  Find those weights.

      for (i=0 ; i<nvars ; i++) {   // Limit to reasonable values
         if (bptr[i] > 5.)
            bptr[i] = 5. ;
         if (bptr[i] < -5.)
            bptr[i] = -5. ;
         }

/*
   The weights for output neuron 'out' are now in place in out_coefs and are
   pointed to by bptr.  Pass through the training set, using the activations
   of the layer just before the output layer, still in sptr->a, to compute
   the activation of the output neuron.  Compare this attained activation to
   the desired in the training sample, and cumulate the mean square error.
   Note that we use nvars-1 in the call to 'activity' because the bias term
   is taken care of in that subroutine.
*/

      for (tset=0 ; tset<tptr->ntrain ; tset++) {// Cumulate err of this output

         dptr = tptr->data + size * tset ;       // Training sample starts here
         aptr = sptr->a + tset * nvars ;         // Inputs to output layer
         diff = activity ( aptr , bptr , nvars-1 ) ; // Find this output

         if (outmod == OUTMOD_AUTO)              // If this is AUTOASSOCIATIVE
            diff -= dptr[out] ;                  // the desired output is input

         else if (outmod == OUTMOD_CLASSIFY) {    // If this is Classification
            if ((int) dptr[nin] == out+1)        // class identifier past inputs
               diff -= NEURON_ON ;
            else 
               diff -= NEURON_OFF ;
            }

         else if (outmod == OUTMOD_GENERAL)      // If this is GENERAL output
            diff -= dptr[nin+out] ;              // output is just past input

         err += diff * diff ;
         }
      } // For each output

   err /= (double) tptr->ntrain * (double) nout ;
   neterr = err ;
   return err ;
}

