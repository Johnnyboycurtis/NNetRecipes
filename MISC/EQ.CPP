/******************************************************************************/
/*                                                                            */
/* EQ - Subroutine to generate equilateral input transformation matrix        */
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

void equilat (
   int n ,       // Number of values, at least 2
   double *mat , // Output n by n-1 matrix of activations
   int range     // 0 for .1 to .9, 1 for -1 to 1
   )
{
   int i, j, k, nm1 ;
   double r, f ;

/*
   Initialize mat[0,0]=-1 and mat[1,0]=1
*/

   nm1 = n - 1 ;
   mat[0] = -1.0 ;
   mat[nm1] = 1.0 ;

/*
   Recursively build rest of matrix
*/

   for (k=2 ; k<n ; k++) {

      // Scale matrix so far
      r = (double) k ;
      f = sqrt ( r * r - 1.0 ) / r ; // Scaling factor
      for (i=0 ; i<k ; i++) {        // Rows of mat so far
         for (j=0 ; j<k-1 ; j++)     // And cols
            mat[i*nm1+j] *= f ;      // Scale it
         }

      // Append a column of all -1/k
      r = -1.0 / r ;                 // Value of new column
      for (i=0 ; i<k ; i++)          // Rows of mat so far
         mat[i*nm1+k-1] = r ;        // New column on right

      // Append new row of all 0's except 1 at end
      for (i=0 ; i<k-1 ; i++)
         mat[k*nm1+i] = 0.0 ;
      mat[k*nm1+k-1] = 1.0 ;
      }

/*
   The above algorithm left the values in -1 to 1.
   If they need to be rescaled .1 to .9, do it now.
*/

   if (range)
      return ;

   i = n * nm1 ;
   while (i--)
      mat[i] = 0.4 * mat[i] + 0.5 ;
}
