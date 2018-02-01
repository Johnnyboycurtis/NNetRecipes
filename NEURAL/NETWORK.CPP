/******************************************************************************/
/*                                                                            */
/*  NETWORK Routines specific to the Network parent class                     */
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

Network::~Network()
{
   return ;
}

/*
--------------------------------------------------------------------------------

   trial_error - Compute the mean square error for the entire training set

--------------------------------------------------------------------------------
*/

double Network::trial_error ( TrainingSet *tptr )
{
   int i, size, tset, tclass ;
   double tot_err, *dptr, diff ;

   if (outmod == OUTMOD_CLASSIFY)   // Compute size of each training sample
      size = nin + 1 ;
   else if (outmod == OUTMOD_AUTO)
      size = nin ;
   else if (outmod == OUTMOD_GENERAL)
      size = nin + nout ;

   tot_err = 0.0 ;  // Total error will be cumulated here

   for (tset=0 ; tset<tptr->ntrain ; tset++) {  // Do all samples

      dptr = tptr->data + size * tset ;     // Point to this sample
      trial ( dptr ) ;                      // Evaluate network for it

      if (outmod == OUTMOD_AUTO) {          // If this is AUTOASSOCIATIVE
         for (i=0 ; i<nout ; i++) {         // then the expected outputs
            diff = *dptr++ - out[i] ;       // are just the inputs
            tot_err += diff * diff ;
            }
         }

      else if (outmod == OUTMOD_CLASSIFY) {  // If this is Classification
         tclass = (int) dptr[nin] - 1 ;     // class is stored after inputs
         for (i=0 ; i<nout ; i++) {         // Recall that train added a
            if (tclass == i)                // fraction so that the above
               diff = NEURON_ON - out[i] ;  // truncation to get tclass is
            else                            // always safe in any radix
               diff = NEURON_OFF - out[i] ;
            tot_err += diff * diff ;
            }
         }

      else if (outmod == OUTMOD_GENERAL) {  // If this is GENERAL output
         dptr += nin ;                      // outputs stored after inputs
         for (i=0 ; i<nout ; i++) {
            diff = *dptr++ - out[i] ;
            tot_err += diff * diff ;
            }
         }

      } // for all tsets

   tot_err /= (double) tptr->ntrain * (double) nout ;
   neterr = tot_err ;
   return tot_err ;
}
