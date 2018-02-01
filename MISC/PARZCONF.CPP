/******************************************************************************/
/*                                                                            */
/*  PARZCONF.CPP - Estimate confidence with Parzen windows                    */
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

double parzconf (    // Returns confidence in alternative hypothesis
   int n0 ,          // Number of samples in null hypothesis collection
   int n1 ,          // And alternative
   double *h0 ,      // Null hypothesis samples
   double *h1 ,      // And alternative
   double sigma ,    // Scale parameter (0.02 to 0.05 best)
   double observed   // Observed value to be classified
   )
{
   int i ;
   double d, sum0, sum1 ;

   sum0 = 1.e-30 ;       // Insurance against dividing by 0 later
   for (i=0 ; i<n0 ; i++) {
      d = (observed - h0[i]) / sigma ;
      sum0 += exp ( - d * d ) ;
      }
   sum0 /= (double) n0 ;

   sum1 = 1.e-30 ;
   for (i=0 ; i<n1 ; i++) {
      d = (observed - h1[i]) / sigma ;
      sum1 += exp ( - d * d ) ;
      }
   sum1 /= (double) n1 ;

   return 100. * sum1 / (sum0 + sum1) ;
}
