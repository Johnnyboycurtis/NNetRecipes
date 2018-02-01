/******************************************************************************/
/*                                                                            */
/*  DOTPROD - Compute dot product of two vectors                              */
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

/*
   This routine is a major eater of time, so it should be made
   as efficient as possible.  Fanatics will code it in assembler.

   The summation loop in this version is opened into groups of
   four.  This saves much overhead in loop counting.  Even more
   important, it avoids breaking the pipeline in RISC and other 
   processors which count on pipelines for speed. 

   Note that we use *(vec1+1) et cetera here.  This is for the 
   benefit of 386 family processors which can include a constant
   address offset in instructions.  Other processors may benefit
   from incrementing each time: *vec1++ et cetera. 
*/
     
double dotprod (
   int n ,         // Length of vectors
   double *vec1 ,  // One of the vectors to be dotted
   double *vec2 )  // The other vector
{
   int k, m ;
   double sum ;

   sum = 0.0 ;  // Will cumulate dot product here
   k = n / 4 ;  // Divide vector into this many groups of 4
   m = n % 4 ;  // This is the remainder of that division

   while (k--) {    // Do each group of 4
      sum += *vec1 * *vec2 ;
      sum += *(vec1+1) * *(vec2+1) ;
      sum += *(vec1+2) * *(vec2+2) ;
      sum += *(vec1+3) * *(vec2+3) ;
      vec1 += 4 ;
      vec2 += 4 ;
      }

   while (m--)      // Do the remainder
      sum += *vec1++ * *vec2++ ;

   return sum ;
}

