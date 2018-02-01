/******************************************************************************/
/*                                                                            */
/*  CONTROL - Routines related to processing user's control commands          */
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

   Prototypes for local routines

--------------------------------------------------------------------------------
*/

static void remove_comments ( char *str ) ;
static void remove_leading_blanks ( char **str ) ;
static void remove_trailing_blanks ( char *str ) ;

/*
--------------------------------------------------------------------------------

   split_control_line - This splits the control line into the command part
      and the parameter part.  It replaces the colon (which terminates the
      command) with a 0, and returns a pointer to the rest of the line.
      It also returns a pointer to the command part of the line, which will
      not be the input line if it started with blanks.  The command part is
      converted to upper case.
      Note that no actual copies are made.  The returned string pointers
      simply point to parts of the input line, which is modified as needed.
      If the input line did not contain a colon, but had characters not
      attributable to comments, n_command is returned 0 and rest returns
      the line.

--------------------------------------------------------------------------------
*/

void split_control_line (
   char *line ,     // Control line input here
   char **command , // Returned pointer to command part
   int *n_command , // Length of command part
   char **rest ,    // Returned pointer to parameter part
   int *n_rest      // Length of parameter part
   )
{
   int i, n ;

   remove_leading_blanks ( &line ) ;
   remove_comments ( line ) ;

   n = strlen ( line ) ;
   for (i=0 ; i<n ; i++) {
      if (line[i] == ':') {                 // Colon terminates command part
         line[i] = 0 ;                      // so terminate this substring
         remove_trailing_blanks ( line ) ;  // Blanks may preceed colon
         *command = line ;                  // Return pointer to command part
         *n_command = strlen ( *command ) ; // and its length
         str_to_upr ( *command ) ;          // Return it upper case
         *rest = line + i + 1 ;             // This is the parameter part
         remove_leading_blanks ( rest ) ;   // which may have leading blanks
         str_to_upr ( *rest ) ;             // Return it upper case
         *n_rest = strlen ( *rest ) ;       // Return its length
         return ;                           // Normal return here
         }
      }

   *n_command = 0 ;  // Error if we get here
   *rest = *command = line ;
   *n_rest = strlen ( line ) ;
}

/*
--------------------------------------------------------------------------------

      get_control_line - Read a control line from a file or stdin

--------------------------------------------------------------------------------
*/

void get_control_line (
   char *control_line ,       // Control line returned here
   int *control_file_number , // Stack pointer for control files (-1 = none)
   FILE *control_files[]      // Stack of file pointers to control files
   )
{
   int n ;

   *control_line = 0 ;   // Always clear any old command in it

   while (*control_file_number >= 0) {  // First priority is control files

      if (fgets ( control_line , CONTROL_LINE_LENGTH ,
                  control_files[*control_file_number] ) == NULL ) { // EOF?
         fclose ( control_files[*control_file_number] ) ; // If so, close it
         --*control_file_number ;                      // and pop it from stack
         continue ;
         }

      n = strlen ( control_line ) ;    // fgets returns the newline at end
      if (control_line[n-1] == '\n')   // so we must remove it
         control_line[--n] = 0 ;       // Using an 'if' is good style

      if (n)                           // It may be that it was a null line
         return ;
      }

   printf ( "\nInput (or QUIT:) " ) ;  // No control files
   gets ( control_line ) ;             // so read from stdin
}

/*
--------------------------------------------------------------------------------

   stack_control_file - Process a CONTROL command by opening, stacking the file 

--------------------------------------------------------------------------------
*/

void stack_control_file (
   char *filename ,           // Full path name of new control file
   int *control_file_number , // Input of current stack position, output inc'd
   FILE *control_files[]      // Stack of file pointers to control files
   )
{
   FILE *fp ;

   if (*control_file_number >= MAX_CONTROL_FILES-1) {
      printf ( "\nERROR... Too many nested control files.  Max is %d.",
               MAX_CONTROL_FILES ) ;
      return ;
      }

   fp = fopen ( filename , "rt" ) ;

   if (fp == NULL) {
      printf ( "\nERROR... Cannot open %s.", filename ) ;
      return ;
      }

   control_files[++*control_file_number] = fp ;
}

/*
--------------------------------------------------------------------------------

   Local utility routines

--------------------------------------------------------------------------------
*/

static void remove_leading_blanks ( char **str )
{
   while (**str == ' ')
      ++*str ;
}

static void remove_trailing_blanks ( char *str )
{
   int n = strlen ( str ) ;
   while (n--) {
      if (str[n] == ' ')
         str[n] = 0 ;
      else
         break ;
      }
}

static void remove_comments ( char *str )
{
   int i, n ;
   n = strlen ( str ) ;
   for (i=0 ; i<n ; i++) {
      if (str[i] == ';') {
         str[i] = 0 ;
         while (i--) {
            if (str[i] == ' ')
               str[i] = 0 ;
            else
               return ;
            }
         }
      }
}


