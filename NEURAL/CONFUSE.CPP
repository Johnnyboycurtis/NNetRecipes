/******************************************************************************/
/*                                                                            */
/*  CONFUSE - All routines related to CLASSIFY confusion                      */
/*                                                                            */
/*    reset_confusion                                                         */
/*    show_confusion                                                          */
/*    save_confusion                                                          */
/*    classify_from_file                                                      */
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
--------------------------------------------------------------------------------

   reset_confusion

--------------------------------------------------------------------------------
*/

void Network::reset_confusion ()
{
   memset ( confusion , 0 , (nout+1) * sizeof(int) ) ;
}

/*
--------------------------------------------------------------------------------

   show_confusion

--------------------------------------------------------------------------------
*/

void Network::show_confusion ()
{
   int i ;
   char *msg ;

   MEMTEXT ( "CONFUSE:show msg" ) ;
   if ((msg = (char *) MALLOC ( (nout+1) * 5 + 11 )) == NULL ) {
      memory_message ( "to SHOW CONFUSION" ) ;
      return ;
      }

   strcpy ( msg , "Confusion:" ) ;
   for (i=0 ; i<nout ; i++)
      sprintf ( msg+5*i+10 , "%5d" , confusion[i] ) ;
   sprintf ( msg+5*nout+10, "%5d", confusion[nout] ) ;
   msg[5*nout+15] = 0 ;
   normal_message ( msg ) ;
   FREE ( msg ) ;
   return ;
}

/*
--------------------------------------------------------------------------------

   save_confusion

--------------------------------------------------------------------------------
*/

void Network::save_confusion ( char *name )
{
   int i ;
   char *msg ;
   FILE *fp ;

   MEMTEXT ( "CONFUSE:save msg" ) ;
   if ((msg = (char *) MALLOC ( (nout+1) * 5 + 1 )) == NULL ) {
      memory_message ( "to SAVE CONFUSION" ) ;
      return ;
      }

/*
   Open the file to which confusion will be written.
   If it already exists, write a newline at its end.
*/

   fp = fopen ( name , "rt" ) ;
   if (fp != NULL) {
      i = 1 ;
      fclose ( fp ) ;
      }
   else
      i = 0 ;

   if ((fp = fopen ( name , "at" )) == NULL) {
      error_message ( "Cannot open SAVE CONFUSION file" ) ;
      FREE ( msg ) ;
      return ;
      }

   if (i)
      fprintf ( fp , "\n" ) ;

/*
   Write confusion
*/

   for (i=0 ; i<nout ; i++)
      sprintf ( msg+5*i , "%5d" , confusion[i] ) ;
   sprintf ( msg+5*nout , "%5d", confusion[nout] ) ;
   msg[5*nout+5] = 0 ;

   fprintf ( fp , "%s", msg ) ;
   fclose ( fp ) ;
   FREE ( msg ) ;
   return ;
}

/*
--------------------------------------------------------------------------------

   classify_from_file

--------------------------------------------------------------------------------
*/

void Network::classify_from_file ( char *name , double thresh )
{
   int i, maxlin, did_any, best ;
   double *inputs, *iptr, maxact ;
   char msg[81], *line, *lptr ;
   FILE *fp ;

/*
   Open the file which contains the data to be classified
*/

   if ((fp = fopen ( name , "rt" )) == NULL) {
      strcpy ( msg , "Cannot open " ) ;
      strcat ( msg , name ) ;
      error_message ( msg ) ;
      return ;
      }

/*
   Allocate for the file lines as read. Conservatively guess length.
   Also allocate for network input vector.
*/

   maxlin = nin * 20 + 100 ;
   if (maxlin < 1024)
      maxlin = 1024 ;

   MEMTEXT ( "CONFUSE:line, inputs" ) ;
   line = (char *) MALLOC ( maxlin ) ;
   inputs = (double *) MALLOC ( nin * sizeof(double) ) ;

   if ((line == NULL)  ||  (inputs == NULL)) {
      memory_message ( "to classify" ) ;
      fclose ( fp ) ;
      if (line != NULL)
         FREE ( line ) ;
      if (inputs != NULL)
         FREE ( inputs ) ;
      return ;
      }

/*
   Read the file.
*/

   did_any = 0 ;  /* If file runs out on first try, ERROR! */

   for (;;) {  // Endless loop reads until file exhausted

      if ((fgets ( line , maxlin , fp ) == NULL) || (strlen ( line ) < 2)) {
         if ((! did_any)  ||  ferror ( fp )) {
            strcpy ( msg , "Problem reading file " ) ;
            strcat ( msg , name ) ;
            error_message ( msg ) ;
            }
         break ;
         }

      lptr = line ;              // Parse the data from this line
      iptr = inputs ;            // This will be the network inputs
      for (i=0 ; i<nin ; i++)
         *iptr++ = ParseDouble ( &lptr ) ;

      did_any = 1 ;              // Flag that at least one found
      trial ( inputs ) ;         // Compute network's outputs

      maxact = -1.e30 ;          // Will keep highest activity here
      best = 0 ;                 // Insurance only (good habit)
      for (i=0 ; i<nout ; i++) { // Find winning output
         if (out[i] > maxact) {
            maxact = out[i] ;
            best = i ;
            }
         }

      if (maxact >= thresh)   // If winner has enough activation
         ++confusion[best] ;  // count it in confusion
      else                    // If too little, throw it
         ++confusion[nout] ;  // in the reject category

      while ((! feof ( fp ))  &&  (line[strlen(line)-1] != '\n'))
         fgets ( line , maxlin , fp ) ; // Line length may exceed maxlin

      if (feof ( fp ))
         break ;
      }  /* Endless loop until a file runs out */

   fclose ( fp ) ;
   MEMTEXT ( "CONFUSE:line, inputs" ) ;
   FREE ( line ) ;
   FREE ( inputs ) ;
}
