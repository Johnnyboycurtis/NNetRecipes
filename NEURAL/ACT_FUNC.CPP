/******************************************************************************/
/*                                                                            */
/*  ACT_FUNC - Compute the activation function                                */
/*                                                                            */
/*  This also contains the routine actderiv which computes the derivative.    */
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
   This uses a piecewise linear approximation.  The routine
   'act_func_init' must be called before the 'act_func' function
   in order to initialize the interpolation table.
*/

static double f_factor, f_f[LTAB_LENGTH], f_d[LTAB_LENGTH] ;

void act_func_init ()
{
   int i ;

   f_factor = (double) (LTAB_LENGTH - 1) / (double) LTAB_MAX ;
   for (i=0 ; i<LTAB_LENGTH ; i++) {
      f_f[i] = 1.0 / (1.0 + exp ( - ((double) i) / f_factor )) ;
      if (i)
         f_d[i-1] = f_f[i] - f_f[i-1] ;
      }
}

double act_func ( double x )
{
   int i ;
   double xd ;

   if (x >= 0.0) {                   // One symmetric half
      xd = x * f_factor ;            // Find location in table
      i = (int) xd ;                 // Subscript in table
      if (i >= (LTAB_LENGTH - 1))    // If outside table
         return f_f[LTAB_LENGTH-1] ; // hold at highest entry
      return f_f[i] + f_d[i] * (xd - i) ; // Else interpolate
      }

   else {                            // Other symmetric half
      xd = -x * f_factor ;
      i = (int) xd ;
      if (i >= (LTAB_LENGTH - 1))
         return 1.0 - f_f[LTAB_LENGTH-1] ;
      return 1.0 - (f_f[i] + f_d[i] * (xd - i)) ;
      }
}

/*
   ACTDERIV - Compute the derivative of the activation function.
   Note that it is written as a function of the output activation,
   f(net), rather than net itself.
   This is all that is needed for the logistic activation function.
   If the user tries other functions, appropriate reparameterization
   must be done.  To make this dependent on the net would
   necessitate storing the nets, which is otherwise unnecessary.
*/

double actderiv ( double f )
{
   return f * (1.0 - f) ;   // Logistic derivative
}

/*
   INVERSE_ACT - Compute the inverse of the activation function.
*/

double inverse_act ( double f )
{
   return -log ( 1.0 / f - 1.0 ) ; // Inverse logistic function
}

