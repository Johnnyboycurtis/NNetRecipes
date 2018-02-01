/******************************************************************************/
/*                                                                            */
/*  WT_SAVE - Save and restore learned weights to/from disk files             */
/*                                                                            */
/*  This returns 1 if error, else 0.                                          */
/*                                                                            */
/*  It handles common errors, opens and closes files, etc.                    */
/*  The actual work is done in the appropriate network virtual function.      */
/*                                                                            */
/*  The 'updating' parameter would be input 1 if this routine is called for   */
/*  automatic updating during training, else 0.                               */
/*                                                                            */
/*  The file header is as follows:                                            */
/*                                                                            */
/*    Header 1:                                                               */
/*      NEURAL WEIGHT FILE - ASCII string for identification purposes         */
/*      version - Version number of this file (currently 1)                   */
/*      model - NETMOD_LAYER or NETMOD_KOH                                    */
/*      n_in - Number of input neurons (sum of lengths)                       */
/*      n_hid1 - Number of hidden layer one neruons                           */
/*      n_hid2 - Ditto for hidden layer two                                   */
/*      n_out - Number of output neurons                                      */
/*      outmodel - AUTO / CLASSIF / GENERAL                                   */
/*      normal - Network dependent normalization method                       */
/*      h2len - Length of header 2, currently 0                               */
/*      48 bytes of unused data, reserved for the future                      */
/*                                                                            */
/*    Header 2:  currently not used                                           */
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

static int write_header ( Network *net , int net_model , FILE *fp ) ;
static Network *read_header ( FILE *fp , int *net_model ) ;

struct wt_header1 {
   char id[20]    ; // Identification string
   int version    ; // Version number of this file (currently 1)
   int model      ; // NETMOD_LAYER or NETMOD_KOH
   int n_in       ; // Number of input neurons (sum of lengths)
   int n_hid1     ; // Number of hidden layer one neruons
   int n_hid2     ; // Ditto for hidden layer two
   int n_out      ; // Number of output neurons
   int outmodel   ; // AUTO / CLASSIF / GENERAL
   int normal     ; // Network dependent normalization method
   int h2len      ; // Length of header 2, currently 0
   char dummy[48] ; // unused data, reserved for the future, saved as 0
   } ;

/*
--------------------------------------------------------------------------------

   WT_SAVE

--------------------------------------------------------------------------------
*/

void wt_save ( Network *net , int net_model , int updating , char *name )
{
   char msg[81] ;
   FILE *fp ;

   if ((fp = fopen ( name , "wb" )) == NULL) {
      if (! updating) {
         strcpy ( msg , "Cannot open WEIGHT file " ) ;
         strcat ( msg , name ) ;
         error_message ( msg ) ;
         }
      return ;
      }

   if (write_header ( net , net_model , fp ) || net->wt_save ( fp ) ||
       ferror ( fp )) {
      if (! updating) {
         strcpy ( msg , "Could not write WEIGHT file " ) ;
         strcat ( msg , name ) ;
         error_message ( msg ) ;
         }
      fclose ( fp ) ;
      return ;
      }

   fclose ( fp ) ;
}


/*
--------------------------------------------------------------------------------

   WT_RESTORE

--------------------------------------------------------------------------------
*/

Network *wt_restore ( char *name , int *net_model )
{
   char msg[81] ;
   FILE *fp ;
   Network *net ;

   if ((fp = fopen ( name , "rb" )) == NULL) {
      strcpy ( msg , "Cannot open WEIGHT file " ) ;
      strcat ( msg , name ) ;
      error_message ( msg ) ;
      return NULL ;
      }

   net = read_header ( fp , net_model ) ;
   if (net == NULL) {
      strcpy ( msg , "Could not read WEIGHT file " ) ;
      strcat ( msg , name ) ;
      error_message ( msg ) ;
      fclose ( fp ) ;
      return NULL ;
      }

   net->wt_restore ( fp ) ;
   fclose ( fp ) ;

   if (! net->ok) {  // Set to 0 if file read error
      strcpy ( msg , "Could not read WEIGHT file " ) ;
      strcat ( msg , name ) ;
      error_message ( msg ) ;
      MEMTEXT ( "WT_SAVE: file error delete network" ) ;
      net->ok = 1 ;  // Tells destructor to free vectors
      delete net ;
      net = NULL ;
      }
   return net ;
}


/*
--------------------------------------------------------------------------------

   write_header - Local routine writes weight file header

--------------------------------------------------------------------------------
*/

static int write_header ( Network *net , int net_model , FILE *fp )
{
   int i ;
   struct wt_header1 h1 ;

   strcpy ( h1.id , "NEURAL WEIGHT FILE" ) ;
   h1.version = 1 ;
   h1.model = net_model ;
   h1.n_in = net->nin ;
   h1.n_out = net->nout ;
   h1.outmodel = net->outmod ;
   h1.h2len = 0 ;
   for (i=0 ; i<48 ; i++)
      (h1.dummy)[i] = (char) 0 ;

   if (net_model == NETMOD_LAYER) {
      h1.n_hid1 = ((LayerNet *) net)->nhid1 ;
      h1.n_hid2 = ((LayerNet *) net)->nhid2 ;
      h1.normal = 0 ;  // Ignored
      }
   else if (net_model == NETMOD_KOH)
      h1.normal = ((KohNet *) net)->normalization ;

   if (! fwrite ( &h1 , sizeof(h1) , 1 , fp ))
      return 1 ;

   return 0 ;
}


/*
--------------------------------------------------------------------------------

   read_header - Local routine reads weight file header

--------------------------------------------------------------------------------
*/

static Network *read_header ( FILE *fp , int *net_model )
{
   int nin, nh1, nh2, nout, outmod ;
   struct wt_header1 h1 ;
   Network *net ;
   KohParams kp ;

   if (! fread ( &h1 , sizeof(h1) , 1 , fp ))
      return NULL ;

   if (strcmp ( h1.id , "NEURAL WEIGHT FILE" )) {
      error_message ( "This is not a NEURAL WEIGHT file." ) ;
      return NULL ;
      }

   *net_model = h1.model ;
   nin = h1.n_in ;
   nh1 = h1.n_hid1 ;
   nh2 = h1.n_hid2 ;
   nout = h1.n_out ;
   outmod = h1.outmodel ;
   kp.normalization = h1.normal ;

   MEMTEXT ( "WT_SAVE: new network for restore" ) ;
   if (*net_model == NETMOD_LAYER)
      net = new LayerNet ( outmod , nin , nh1 , nh2 , nout , 1 , 0 ) ;
   else if (*net_model == NETMOD_KOH)
      net = new KohNet ( nin , nout , &kp , 1 , 0 ) ;
   else {
      error_message ( "WEIGHT file specified illegal network model" ) ;
      return NULL ;
      }

   if ((net == NULL)  ||  (! net->ok)) {  // Malloc failure?
      memory_message ( "to create network." ) ;
      if (net != NULL)
         delete net ;
      return NULL ;
      }

   return net ;
}
