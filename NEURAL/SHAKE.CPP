/******************************************************************************/
/*                                                                            */
/*  SHAKE - Randomly perturb a point                                          */
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

void shake ( int nvars , double *center , double *x , double temp )
{
   double r ;

/*
   Recall that the variance of a uniform deviate on 0-1 is 1/12.
   Adding four such random variables multiplies the variance by 4,
   while dividing by 2 divides the variance by 4.
*/

   temp *= 3.464101615 / (2. * longrandmax () ) ; // SQRT(12)=3.464...

   while (nvars--) {
      r = (double) longrand() + (double) longrand() -
          (double) longrand() - (double) longrand() ;
      *x++ = *center++ + temp * r ;
      }
}
