/******************************************************************************/
/*                                                                            */
/*  MESSAGES - All routines for issuing messages to user                      */
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

void error_message ( char *msg )
{
   printf ( "\nERROR: %s", msg ) ;
}

void warning_message ( char *msg )
{
   printf ( "\nWARNING: %s", msg ) ;
}

void normal_message ( char *msg )
{
   printf ( "\n%s", msg ) ;
}

void progress_message ( char *msg )
{
   printf ( "%s", msg ) ;
}

void memory_message ( char *msg )
{
   printf ( "\nInsufficient memory %s", msg ) ;
}

void str_to_upr ( char *str )
{
   char *ptr = str-1;
   while (*++ptr)
      *ptr = toupper ( *ptr ) ;
}
