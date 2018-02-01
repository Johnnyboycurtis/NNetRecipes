/******************************************************************************/
/*                                                                            */
/*  EXECUTE - All routines related to AUTO and GENERAL network execution      */
/*                                                                            */
/*  CLASSIF stuff is in CONFUSE.CPP                                           */
/*  But note that EXECUTE is certainly valid in CLASSIF mode.  It outputs     */
/*  neuron activations, just like in the other modes.                         */
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

void Network::execute_from_file (
   char *dataname ,  // Input file name
   char *outname     // Output file name
   )
{
   int i, maxlin, did_any ;
   double *inputs, *iptr ;
   char msg[81], *line, *lptr ;
   FILE *fp_in, *fp_out ;

/*
   Open the file which contains the data to be classified
*/

   if ((fp_in = fopen ( dataname , "rt" )) == NULL) {
      strcpy ( msg , "Cannot open input data file " ) ;
      strcat ( msg , dataname ) ;
      error_message ( msg ) ;
      return ;
      }

/*
   Open the file to which outputs will be written.
   If it already exists, write a newline at its end.
*/

   fp_out = fopen ( outname , "rt" ) ;
   if (fp_out != NULL) {
      did_any = 1 ;
      fclose ( fp_out ) ;
      }
   else
      did_any = 0 ;

   if ((fp_out = fopen ( outname , "at" )) == NULL) {
      strcpy ( msg , "Cannot open output file " ) ;
      strcat ( msg , outname ) ;
      error_message ( msg ) ;
      fclose ( fp_in ) ;
      return ;
      }

   if (did_any)
      fprintf ( fp_out , "\n" ) ;

/*
   Allocate for the file lines as read. Conservatively guess length.
   Also allocate for network input vector.
*/

   maxlin = nin * 20 + 100 ;
   if (maxlin < 1024)
      maxlin = 1024 ;

   MEMTEXT ( "EXECUTE:line, inputs" ) ;
   line = (char *) MALLOC ( maxlin ) ;
   inputs = (double *) MALLOC ( nin * sizeof(double) ) ;

   if ((line == NULL)  ||  (inputs == NULL)) {
      memory_message ( "to execute" ) ;
      fclose ( fp_in ) ;
      fclose ( fp_out ) ;
      if (line != NULL)
         FREE ( line ) ;
      if (inputs != NULL)
         FREE ( inputs ) ;
      return ;
      }

/*
   Read and process the file.
*/

   did_any = 0 ;  /* If file runs out on first try, ERROR! */

   for (;;) {  // Endless loop reads until file exhausted

      if ((fgets ( line , maxlin , fp_in ) == NULL) || (strlen ( line ) < 2)) {
         if ((! did_any)  ||  ferror ( fp_in )) {
            strcpy ( msg , "Problem reading file " ) ;
            strcat ( msg , dataname ) ;
            error_message ( msg ) ;
            }
         break ;
         }

      lptr = line ;           // Parse the data from this line
      iptr = inputs ;         // This will be the network inputs
      for (i=0 ; i<nin ; i++)
         *iptr++ = ParseDouble ( &lptr ) ;

      if (did_any)            // Start each new case on new line
         fprintf ( fp_out , "\n" ) ;

      did_any = 1 ;           // Flag that at least one found
      trial ( inputs ) ;      // Compute network's outputs

      for (i=0 ; i<nout ; i++)
         fprintf ( fp_out , " %.4lf", out[i] ) ;

      while ((! feof ( fp_in ))  &&  (line[strlen(line)-1] != '\n'))
         fgets ( line , maxlin , fp_in ) ; // Line length may exceed maxlin

      if (feof ( fp_in ))
         break ;
      }  /* Endless loop until a file runs out */

   MEMTEXT ( "EXECUTE:line, inputs" ) ;
   fclose ( fp_in ) ;
   fclose ( fp_out ) ;
   FREE ( line ) ;
   FREE ( inputs ) ;
}
