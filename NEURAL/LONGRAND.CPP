/******************************************************************************/
/*                                                                            */
/*  LONGRAND - Generate 32 bit period-free random numbers.                    */
/*                                                                            */
/*    void slongrand ( long iseed ) - Set the random seed                     */
/*    long longrand () - Return a 32 bit random integer in range 0-IM-1       */
/*    long longrandmax () - Return the maximum number which pf_rand may return*/
/*    double unifrand () - Return uniform random in [0,1)                     */
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

#define IM 714025L      // These constants are from Press et. al.
#define IA 1366L        // "Numerical Recipes in C"
#define IC 150889L      // Do not tamper with them unless you are an expert
#define TABLE_LENGTH 97 // Length of shuffle table

static long seed = 797 ;           // Keep the current seed here
static long randout ;              // Random output
static long table[TABLE_LENGTH] ;  // Keep shuffle table here
static int table_initialized = 0 ; // Has it been initialized?

/*
   Set the random seed
*/

void slongrand ( long iseed )
{
   seed = iseed ;          // We keep the static seed here
   table_initialized = 0 ; // Must also rebuild table with it!
}

/*
   Return the maximum random number
*/

long longrandmax ()
{
   return IM - 1L ;
}

/*
   This is the actual random number generator
*/

long longrand ()     // Return the next random number
{
   int i ;

   if (! table_initialized) {  // Initialize shuffle table before use
      table_initialized = 1 ;  // Flag to avoid more inits
      for (i=0 ; i<TABLE_LENGTH ; i++) { // Fill entire table
         seed = (IA * seed + IC) % IM ;  // Make a random number
         table[i] = seed ;               // Put it in the table
         }
      seed = (IA * seed + IC) % IM ;     // One more random number
      randout = seed ;                   // for upcoming first use
      }

   i = (int) ((double) TABLE_LENGTH * (double) randout / (double) IM) ;
   randout = table[i] ;            // This output comes from table
   seed = (IA * seed + IC) % IM ;  // Make new random number
   table[i] = seed ;               // to replace used entry
   return randout ;                // then return old entry
}


double unifrand ()
{
   return (double) longrand () / (double) IM ;
}
