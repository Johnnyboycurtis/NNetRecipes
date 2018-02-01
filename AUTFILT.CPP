/******************************************************************************/
/*                                                                            */
/*  Autoassociative Filtering                                                 */
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

/*
   This is the declaration for a subroutine which executes a trained
   autoassociative neural network.
*/

void evaluate (    // Execute a trained autoassociative neural network
   int n ,         // Number of inputs and outputs
   double *in ,    // Input to the network
   double *out     // Outputs returned here
   ) ;

/*
   This is the actual filtering routine
*/

void aut_filt (    // Autoassociative filter routine
   int n ,         // Number of network inputs and outputs (window width)
   int m ,         // Length of the series to be filtered
   double *in ,    // Input series to be filtered
   double *out ,   // Output series of filtered input
   double *work    // Work vector n long
   )
{
   int i, j, n_presentations, lim ;

   for (i=0 ; i<m ; i++)  // Zero output in preparation for
      out[i] = 0. ;       // cumulating window results

   n_presentations = m - n + 1 ;  // This many windows will be presented

   for (i=0 ; i<n_presentations ; i++) {  // Loop does all presentations
      evaluate ( n , in+i , work ) ;      // Run network for this window
      for (j=0 ; j<n ; j++)               // Cumulate sum of
         out[i+j] += work[j] ;            // all presentations
      }

   if (n_presentations > 1) {             // Divide to get means
      lim = (n_presentations < n) ? n_presentations : n ; // Max possible
      for (i=0 ; i<m ; i++) {             // Loop divides each output
         j = m - i ;                      // If we are at the far end
         if ((i+1) < j)                   // or at the near end we must
            j = i + 1 ;                   // compute the number in sum
         if (j > lim)                     // In the interior the
            j = lim ;                     // max went into the sum
         out[i] /= (double) j ;           // Divide to get mean
         }
      }
}
