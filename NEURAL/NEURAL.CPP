/******************************************************************************/
/*                                                                            */
/*  NEURAL - Main program for implementing all neural network functions       */
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
   These are used in MEM.CPP for diagnostic memory allocation
*/

extern int mem_log ;       // Keep a log on disk?
extern long mem_used ;     // Keeps track of memory usage
extern char mem_name[80] ; // Full path name of log file

/*
   Forward declarations for static routines
*/

static int check_anneal ( char *command , char *rest ,
                          struct AnnealParams *anneal_params ) ;
static int check_genetic ( char *command , char *rest ,
                           struct GenInitParams *geninit_params ) ;
static int check_kohonen ( char *command , char *rest ,
                           struct KohParams *koh_params , Network **net ) ;
static int check_learn_params ( char *command , char *rest ,
                               struct LearnParams *learn_params , int netmod ) ;
static int get_koh_init ( char *rest ) ;
static int get_layer_init ( char *rest ) ;
static int get_yn ( char *msg ) ;
static int ok_to_clear_tset ( TrainingSet **tset ) ;
static int ok_to_clear_weights ( Network **network ) ;

/*
--------------------------------------------------------------------------------

   Main entry point

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{

/*
   Declarations of local variables
*/

/*
   User's command control line related variables are here.
   Control_file_number and control_files permit nesting of 'CONTROL' commands.
   If control_file_number equals -1, control commands are read from stdin.
   Otherwise they are read from that file in FILE *control_files.
   Up to MAX_CONTROL_FILES can be stacked.
*/

   int control_file_number = -1 ;           // Stack pointer for control files
   FILE *control_files[MAX_CONTROL_FILES] ; // This is the stack

   char *control_line ;    // User's commands here
   char *command, *rest ;  // Pointers to its command and parameter parts
   int n_command, n_rest ; // Lengths of those parts

/*
   These are network parameters which may be set by the user via commands.
   They are initialized to defaults which indicate that the user has not
   yet set them.  As they are set, their current values are placed here.
   When learning is done for a network, their values are copied from here
   into the network object.  When a network is read, the object's values
   are copied from it to here.  Otherwise, these variables are not used;
   the values in the network object itself are used.  The only purpose of
   these variables is to keep track of current values.
*/

   int net_model = -1 ;     // Network model (see NETMOD_? in CONST.H)
   int out_model = -1 ;     // Output model (see OUTMOD_? in CONST.H)
   int n_inputs = -1 ;      // Number of input neurons
   int n_outputs = -1 ;     // Number of output neurons
   int n_hidden1 = -1 ;     // Number of hidden layer one neurons
   int n_hidden2 = -1 ;     // Ditto layer 2 (0 if just one hidden layer)


   TrainingSet *tset = NULL ;            // Training set here
   Network *network = NULL ;             // Network here
   struct LearnParams learn_params ;     // General learning parameters
   struct AnnealParams anneal_params ;   // Simulated annealing parameters
   struct GenInitParams geninit_params ; // Genetic initialization parameters
   struct KohParams koh_params ;         // Kohonen parameters

   int classif_output = -1 ;  // Current class (0=reject) for classif training
   char out_file[80] = "" ;   // File for EXECUTE output
   double threshold ;         // CLASSIFY confusion reject cutoff

/*
   Miscellaneous variables
*/

   int i, n, m ;
   double p ;
   char *msg ;
   FILE *fp ;

/*
--------------------------------------------------------------------------------

   Program starts here.

   Verify that a careless user didn't fail to set the integer size
   correctly when compiling.

--------------------------------------------------------------------------------
*/

#if VERSION_16_BIT
   if (sizeof(int) > 2) {
      printf ( "\nRecompile with VERSION_16_BIT set to 0 in CONST.H" ) ;
      exit ( 1 ) ;
      }
#else
   if (sizeof(int) < 4) {
      printf ( "\nRecompile with VERSION_16_BIT set to 1 in CONST.H" ) ;
      exit ( 1 ) ;
      }
#endif

printf ( "\nNEURAL - Program to train and test neural networks" ) ;
printf("\nCopyright (c) 1993 by Academic Press, Inc.");
printf("\nAll rights reserved.  Permission is hereby granted, until further notice,");
printf("\nto make copies of this diskette, which are not for resale, provided these");
printf("\ncopies are made from this master diskette only, and provided that the");
printf("\nfollowing copyright notice appears on the diskette label:");
printf("\n(c) 1993 by Academic Press, Inc.");
printf("\nExcept as previously stated, no part of the computer program embodied in");
printf("\nthis diskette may be reproduced or transmitted in any form or by any means,");
printf("\nelectronic or mechanical, including input into storage in any information");
printf("\nsystem for resale, without permission in writing from the publisher.");
printf("\nProduced in the United States of America.");
printf("\nISBN 0-12-479041-0");

/*
   Process command line parameters
*/

   mem_name[0] = 0 ;  // Default is no memory allocation file

   for (i=1 ; i<argc ; i++) {  // Process all command line args
      str_to_upr ( argv[i] ) ; // Easier if all upper case

      if (! strcmp ( argv[i] , "/DEBUG" )) {
         sscanf ( argv[++i] , "%s" , mem_name ) ;
         if ((strlen ( mem_name ) > 1)  ||  ! isalpha ( mem_name[0] )) {
            printf ( "\nIllegal DEBUG drive (%s); must be 1 letter." ) ;
            exit ( 1 ) ;
            }
         continue ;
         }

      printf ( "\nIllegal command line parameter (%s)", argv[i] ) ;
      exit ( 1 ) ;
      }

/*
   Initialize memory allocation debugging
*/

   if (strlen ( mem_name )) {
      strcat ( mem_name , ":mem.log" ) ;
      fp = fopen ( mem_name , "wt" ) ;
      if (fp == NULL) {
         printf ( "\nCannot open debugging file %s", mem_name ) ;
         exit ( 1 ) ;
         }
      fclose ( fp ) ;
      mem_log = 1 ;
      }
   else 
      mem_log = 0 ;

   mem_used = 0 ;

/*
   Initialize defaults
*/

   learn_params.init = -1 ;
   learn_params.quit_err = 0.0 ;
   learn_params.retries = 32767 ;

   anneal_params.temps0 = 3 ;
   anneal_params.temps = 4 ;
   anneal_params.iters0 = 50 ;
   anneal_params.iters = 20 ;
   anneal_params.setback0 = 50 ;
   anneal_params.setback = 20 ;
   anneal_params.start0 = 3.0 ;
   anneal_params.start = 4.0 ;
   anneal_params.stop0 = 1.0 ;
   anneal_params.stop = 0.02 ;

   geninit_params.pool = 50 ;
   geninit_params.gens = 3 ;
   geninit_params.climb = 0 ;
   geninit_params.overinit = 1.5 ;
   geninit_params.pcross = 0.8 ;
   geninit_params.pmutate = 0.0001 ;

   koh_params.normalization = 0 ;  // 0=multiplicative, 1=Z 
   koh_params.learn_method = 1 ;   // 0=additive, 1=subtractive
   koh_params.rate = 0.4 ;         // learning rate
   koh_params.reduction = 0.99 ;   // learning rate reduction

   learn_params.ap = &anneal_params ;
   learn_params.gp = &geninit_params ;
   learn_params.kp = &koh_params ;

   act_func_init () ; // Initialize interpolation table for activation function

   MEMTEXT ( "NEURAL: control_line, msg" ) ;
   if (((control_line = (char *) MALLOC ( CONTROL_LINE_LENGTH+1 )) == NULL)
    || ((msg = (char *) MALLOC ( CONTROL_LINE_LENGTH+1 )) == NULL)) {
      printf ( "\nInsufficient memory" ) ;
      exit ( 1 ) ;
      }

/*
   Main loop processes all commands
*/

   for (;;) {

      get_control_line ( control_line , &control_file_number, control_files ) ;

      split_control_line ( control_line , &command , &n_command ,
                           &rest , &n_rest ) ;

      if (! n_command) {
         if (n_rest) {
            sprintf ( msg , "No colon after command: %s", rest ) ;
            error_message ( msg ) ;
            }
         continue ;
         }

      sprintf ( msg , "%s : %s", command, rest ) ;
      normal_message ( msg ) ;

/*
   Act on the command
*/

      if (! strcmp ( command , "QUIT" ))
         break ;

      if (! strcmp ( command , "CONTROL" )) {
         stack_control_file ( rest , &control_file_number , control_files ) ;
         continue ;
         }

      if (! strcmp ( command , "NETWORK MODEL" )) {
         if (! strcmp ( rest , "LAYER" ))
            n = NETMOD_LAYER ;
         else if (! strcmp ( rest , "KOHONEN" ))
            n = NETMOD_KOH ;
         else {
            sprintf ( msg , "Illegal NETWORK MODEL: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (net_model == n)
            continue ;
         if (ok_to_clear_weights( &network )) {
            net_model = n ;
            learn_params.init = -1 ;
            }
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "OUTPUT MODEL" )) {
         if (! strcmp ( rest , "CLASSIFY" ))
            n = OUTMOD_CLASSIFY ;
         else if (! strcmp ( rest , "AUTO" ))
            n = OUTMOD_AUTO ;
         else if (! strcmp ( rest , "GENERAL" ))
            n = OUTMOD_GENERAL ;
         else {
            sprintf ( msg , "Illegal OUTPUT MODEL: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (out_model == n)
            continue ;
         if ((ok_to_clear_tset( &tset )) && (ok_to_clear_weights( &network)))
            out_model = n ;
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "N INPUTS" )) {
         m = sscanf ( rest , "%d" , &n ) ;
         if ((m <= 0)  ||  (n <= 0)  ||  (n > MAX_INPUTS)) {
            sprintf ( msg , "Illegal N INPUTS: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (n_inputs == n)
            continue ;
         if ((ok_to_clear_tset( &tset)) && (ok_to_clear_weights(&network)))
            n_inputs = n ;
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "N OUTPUTS" )) {
         m = sscanf ( rest , "%d" , &n ) ;
         if ((m <= 0)  ||  (n <= 0)  ||  (n > MAX_OUTPUTS)) {
            sprintf ( msg , "Illegal N OUTPUTS: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (n_outputs == n)
            continue ;
         if ((ok_to_clear_tset( &tset)) && (ok_to_clear_weights(&network)))
            n_outputs = n ;
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "N HIDDEN1" )) {
         m = sscanf ( rest , "%d" , &n ) ;
         if ((m <= 0)  ||  (n < 0)  ||  (n > MAX_HIDDEN)) {
            sprintf ( msg , "Illegal N HIDDEN1: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (n_hidden1 == n)
            continue ;
         if (ok_to_clear_weights( &network ))
            n_hidden1 = n ;
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "N HIDDEN2" )) {
         m = sscanf ( rest , "%d" , &n ) ;
         if ((m <= 0)  ||  (n < 0)  ||  (n > MAX_HIDDEN)) {
            sprintf ( msg , "Illegal N HIDDEN2: %s", rest ) ;
            error_message ( msg ) ;
            continue ;
            }
         if (n  &&  ! n_hidden1) {
            error_message ( "N HIDDEN2 must be 0 if N HIDDEN1 IS 0." ) ;
            continue ;
            }
         if (n_hidden2 == n)
            continue ;
         if (ok_to_clear_weights( &network ))
            n_hidden2 = n ;
         else
            warning_message ( "Command aborted" ) ;
         continue ;
         }

      if (! strcmp ( command , "TRAIN" )) {
         if ((out_model == OUTMOD_AUTO)  &&  (n_outputs != n_inputs)) {
            warning_message ( "Setting N OUTPUTS = N INPUTS" ) ;
            n_outputs = n_inputs ;
            }
         if (out_model <= 0)
            error_message ( "TRAIN used before OUTPUT MODEL set." ) ;
         else if (n_inputs <= 0)
            error_message ( "TRAIN used before N INPUTS set." ) ;
         else if (n_outputs <= 0)
            error_message ( "TRAIN used before N OUTPUTS set." ) ;
         else if ((net_model != NETMOD_KOH) && (out_model == OUTMOD_CLASSIFY)
                  &&  (classif_output < 0))
            error_message( "CLASSIFY output mode but CLASSIFY OUTPUT not set.");
         else if ((net_model == NETMOD_KOH)  &&  (out_model != OUTMOD_CLASSIFY))
            error_message( "KOHONEN network requires CLASSIFY output mode.");
         else {
            if (tset == NULL) {
               MEMTEXT ( "NEURAL: new tset" ) ;
               tset = new TrainingSet ( out_model , n_inputs , n_outputs ) ;
               }
            tset->train ( rest , classif_output ) ;
            }
         continue ;
         }

      if (check_anneal ( command , rest , &anneal_params ))
         continue ;

      if (check_genetic ( command , rest , &geninit_params ))
         continue ;

      if (check_kohonen ( command , rest , &koh_params , &network ))
         continue ;

      if (check_learn_params ( command , rest , &learn_params , net_model ))
         continue ;

      if (! strcmp ( command , "LEARN" )) {
         if ((tset == NULL)  ||  (tset->ntrain == 0)) {
            error_message ( "Cannot LEARN; No training set exists." ) ;
            continue ;
            }
         if ((net_model == NETMOD_KOH)  &&  (out_model != OUTMOD_CLASSIFY)) {
            error_message( "KOHONEN network requires CLASSIFY output mode.");
            continue ;
            }
         if (learn_params.init < 0) {
            error_message( "Initialization method not set.");
            continue ;
            }
         if (network == NULL) {
            if (net_model == NETMOD_LAYER) {
               if (n_hidden1 < 0) {
                  error_message ( "LEARN used before N HIDDEN1 set." ) ;
                  continue ;
                  }
               else if (n_hidden2 < 0) {
                  error_message ( "LEARN used before N HIDDEN2 set." ) ;
                  continue ;
                  }
               else {
                  MEMTEXT ( "NEURAL: new LayerNet" ) ;
                  network = new LayerNet ( out_model , n_inputs , n_hidden1 ,
                                           n_hidden2 , n_outputs , 1 , 1 ) ;
                  }
               }
            else if (net_model == NETMOD_KOH) {
               MEMTEXT ( "NEURAL: new KohNet" ) ;
               network = new KohNet ( n_inputs , n_outputs ,
                                      &koh_params , 1 , 1 ) ;
               }
            }
         if ((network == NULL)  ||  (! network->ok)) {  // Malloc failure?
            memory_message ( "to create network." ) ;
            if (network != NULL) {
               delete network ;
               network = NULL ;
               }
            continue ;
            }
         network->learn ( tset , &learn_params ) ;
         if (network->neterr > 0.999999) {  // Indicates massive failure
            MEMTEXT ( "NEURAL: learn failure delete network" ) ;
            delete network ;
            network = NULL ;
            }
         else {
            sprintf ( msg , "Final error = %.4lf%% of max possible",
                      100.0 * network->neterr ) ;
            normal_message ( msg ) ;
            }
         continue ;
         }

      if (! strcmp ( command , "SAVE WEIGHTS" )) {
         if (network == NULL)
            error_message ( "There are no learned weights to save." ) ;
         else
            wt_save ( network , net_model , 0 , rest ) ;
         continue ;
         }

      if (! strcmp ( command , "RESTORE WEIGHTS" )) {
         if (network != NULL) {
            MEMTEXT ( "NEURAL: delete network for restore" ) ;
            delete network ;
            network = NULL ;
            }
         network = wt_restore ( rest , &net_model ) ;
         if (network == NULL)
            continue ;
         if (tset != NULL) {
            if ((tset->nin != network->nin)
             || (tset->nout != network->nout)
             || (tset->outmod != network->outmod)) {
               error_message ( "Network conflicts with existing training set.");
               continue ;
               }
            }
         out_model = network->outmod ;
         n_inputs = network->nin ;
         n_outputs = network->nout ;
         if (net_model == NETMOD_LAYER) {
            n_hidden1 = ((LayerNet*) network)->nhid1 ;
            n_hidden2 = ((LayerNet*) network)->nhid2 ;
            }
         if (net_model == NETMOD_KOH)
            koh_params.normalization = ((KohNet *) network)->normalization ;
         learn_params.init = -1 ;
         continue ;
         }

      if (! strcmp ( command , "CLEAR TRAINING" )) {
         if (tset != NULL) {
            MEMTEXT ( "NEURAL: delete tset" ) ;
            delete tset ;
            tset = NULL ;
            }
         continue ;
         }

      if (! strcmp ( command , "CLEAR WEIGHTS" )) {
         if (network != NULL) {
            MEMTEXT ( "NEURAL: delete network" ) ;
            delete network ;
            network = NULL ;
            }
         continue ;
         }

      if (! strcmp ( command , "CLASSIFY OUTPUT" )) {
         if (net_model == NETMOD_KOH) {
            error_message ( "Cannot specify output for KOHONEN model." ) ;
            continue ;
            }
         if (n_outputs < 0) {
            error_message ( "CLASSIFY OUTPUT used before N OUTPUTS set." ) ;
            continue ;
            }
         if (out_model != OUTMOD_CLASSIFY) {
            error_message
                  ( "CLASSIFY OUTPUT only valid when OUTPUT MODEL:CLASSIFY" ) ;
            continue ;
            }
         m = sscanf ( rest , "%d" , &n ) ;
         if ((m <= 0)  ||  (n < 0)) {
            sprintf ( msg , "Illegal CLASSIFY OUTPUT: %s", rest ) ;
            error_message ( msg ) ;
            }
         else if (n > n_outputs) {
            sprintf ( msg , "CLASSIFY OUTPUT (%d) exceeds N OUTPUTS (%d)",
                      n, n_outputs ) ;
            error_message ( msg ) ;
            }
         else
            classif_output = n ;
         continue ;
         }


      if (! strcmp ( command , "OUTPUT FILE" )) {
         strcpy ( out_file , rest ) ;
         continue ;
         }

      if (! strcmp ( command , "EXECUTE" )) {
         if (network == NULL)
            error_message ( "There is no trained network" ) ;
         else
            network->execute_from_file ( rest , out_file ) ;
         continue ;
         }

      if (! strcmp ( command , "CLASSIFY" )) {
         if (network == NULL)
            error_message ( "There is no trained network" ) ;
         else if (out_model != OUTMOD_CLASSIFY)
            error_message ( "CLASSIFY valid only in CLASSIFY output mode" ) ;
         else
            network->classify_from_file ( rest , threshold ) ;
         continue ;
         }

      if (! strcmp ( command , "RESET CONFUSION" )) {
         if (network == NULL)
            error_message ( "There is no trained network" ) ;
         else
            network->reset_confusion () ;
         continue ;
         }

      if (! strcmp ( command , "CONFUSION THRESHOLD" )) {
         p = atof ( rest ) ;
         if ((p < 0.0)  ||  (p > 100.0)) {
            sprintf ( msg , "Illegal CONFUSION THRESHOLD: %s", rest ) ;
            error_message ( msg ) ;
            }
         else
            threshold = p / 100.0 ;
         continue ;
         }

      if (! strcmp ( command , "SHOW CONFUSION" )) {
         if (network == NULL)
            error_message ( "There is no trained network" ) ;
         else if (out_model != OUTMOD_CLASSIFY)
            error_message ( "CONFUSION valid only in CLASSIFY output mode" ) ;
         else
            network->show_confusion () ;
         continue ;
         }

      if (! strcmp ( command , "SAVE CONFUSION" )) {
         if (network == NULL)
            error_message ( "There is no trained network" ) ;
         else if (out_model != OUTMOD_CLASSIFY)
            error_message ( "CONFUSION valid only in CLASSIFY output mode" ) ;
         else
            network->save_confusion ( rest ) ;
         continue ;
         }

      sprintf ( msg , "Unknown command: %s", command ) ;
      error_message ( msg ) ;

      } // Endless command loop

   MEMTEXT ( "NEURAL: control_line, msg" ) ;
   FREE ( control_line ) ;
   FREE ( msg ) ;
   MEMCLOSE () ;
   return 0 ;
}

/*
--------------------------------------------------------------------------------

   Static routines

--------------------------------------------------------------------------------
*/

static int ok_to_clear_tset( TrainingSet **tset )
{
   if (*tset == NULL)
      return 1 ;

   if (get_yn ( "Existing training set must be cleared.  OK" )) {
      MEMTEXT ( "NEURAL: delete tset" ) ;
      delete *tset ;
      *tset = NULL ;
      return 1 ;
      }
   else
      return 0 ;
}

static int ok_to_clear_weights( Network **network )
{
   if (*network == NULL)
      return 1 ;

   if (get_yn ( "Existing learned weights must be cleared.  OK" )) {
      MEMTEXT ( "NEURAL: delete network" ) ;
      delete *network ;
      *network = NULL ;
      return 1 ;
      }
   else
      return 0 ;
}

static int get_yn ( char *msg )
{
   int key ;
   printf ( "\n%s? (Y/N):", msg ) ;
   for (;;) {
      key = getch () ;
      if ((key == 'y')  ||  (key == 'Y'))
         return 1 ;
      if ((key == 'n')  ||  (key == 'N'))
         return 0 ;
      printf ( " Y or N:" ) ;
      }
}

static int get_layer_init ( char *rest )
{
   if (! strcmp ( rest , "NOINIT" ))
      return 0 ;
   else if (! strcmp ( rest , "ANNEAL NOREGRESS" ))
      return 1 ;
   else if (! strcmp ( rest , "ANNEAL" ))
      return 2 ;
   else if (! strcmp ( rest , "GENETIC" ))
      return 3 ;
   else if (! strcmp ( rest , "REGRESS" ))
      return 4 ;
   else
      return -1 ;
}

static int get_koh_init ( char *rest )
{
   if (! strcmp ( rest , "NOINIT" ))
      return 0 ;
   else if (! strcmp ( rest , "RANDOM" ))
      return 1 ;
   else
      return -1 ;
}

static int check_anneal ( char *command , char *rest ,
                          struct AnnealParams *anneal_params )
{
   int m, n ;
   double p ;
   char msg[80] ;

   if (! strcmp ( command , "ANNEAL INIT TEMPS" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal ANNEAL TEMPS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->temps0 = n ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL TEMPS" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal ANNEAL INIT TEMPS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->temps = n ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL INIT ITERS" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal ANNEAL INIT ITERS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->iters0 = n ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL ITERS" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal ANNEAL ITERS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->iters = n ;
      return 1 ;
      }
    if (! strcmp ( command , "ANNEAL INIT SETBACK" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n < 0)) {
         sprintf ( msg , "Illegal ANNEAL INIT SETBACK: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->setback0 = n ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL SETBACK" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n < 0)) {
         sprintf ( msg , "Illegal ANNEAL SETBACK: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->setback = n ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL INIT START" )) {
      p = atof ( rest ) ;
      if (p <= 0.0) {
         sprintf ( msg , "Illegal ANNEAL INIT START: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->start0 = p ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL START" )) {
      p = atof ( rest ) ;
      if (p <= 0.0) {
         sprintf ( msg , "Illegal ANNEAL START: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->start = p ;
      return 1 ;
      }

   if (! strcmp ( command , "ANNEAL INIT STOP" )) {
      p = atof ( rest ) ;
      if (p <= 0.0) {
         sprintf ( msg , "Illegal ANNEAL INIT STOP: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->stop0 = p ;
      return 1 ;
      }


   if (! strcmp ( command , "ANNEAL STOP" )) {
      p = atof ( rest ) ;
      if (p <= 0.0) {
         sprintf ( msg , "Illegal ANNEAL STOP: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         anneal_params->stop = p ;
      return 1 ;
      }
   return 0 ;
}

static int check_genetic ( char *command , char *rest ,
                           struct GenInitParams *geninit_params )
{
   int m, n ;
   double p ;
   char msg[80] ;

   if (! strcmp ( command , "GENETIC INIT POOL" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal GENETIC INIT POOL: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         geninit_params->pool = n ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT GENS" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n <= 0)) {
         sprintf ( msg , "Illegal GENETIC INIT GENS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         geninit_params->gens = n ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT CLIMB" )) {
      geninit_params->climb = 1 ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT NOCLIMB" )) {
      geninit_params->climb = 0 ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT CROSS" )) {
      p = atof ( rest ) ;
      if ((p < 0.0)  ||  (p > 1.0)) {
         sprintf ( msg , "Illegal GENETIC INIT CROSS: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         geninit_params->pcross = p ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT MUTATE" )) {
      p = atof ( rest ) ;
      if ((p < 0.0)  ||  (p > 1.0)) {
         sprintf ( msg , "Illegal GENETIC INIT MUTATE: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         geninit_params->pmutate = p ;
      return 1 ;
      }

   if (! strcmp ( command , "GENETIC INIT OVERINIT" )) {
      p = atof ( rest ) ;
      if (p < 1.0) {
         sprintf ( msg , "Illegal GENETIC INIT OVERINIT: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         geninit_params->overinit = p ;
      return 1 ;
      }

   return 0 ;
}

static int check_kohonen ( char *command , char *rest ,
                           struct KohParams *koh_params , Network **net )
{
   double p ;
   char msg[80] ;

   if (! strcmp ( command , "KOHONEN NORMALIZATION MULTIPLICATIVE" )) {
      if (koh_params->normalization == 0)
         return 1 ;
      if (ok_to_clear_weights( net ))
         koh_params->normalization = 0 ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN NORMALIZATION Z" )) {
      if (koh_params->normalization == 1)
         return 1 ;
      if (ok_to_clear_weights( net ))
         koh_params->normalization = 1 ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN LEARN ADDITIVE" )) {
      koh_params->learn_method = 0 ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN LEARN SUBTRACTIVE" )) {
      koh_params->learn_method = 1 ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN LEARNING RATE" )) {
      p = atof ( rest ) ;
      if ((p < 0.0)  ||  (p > 1.0)) {
         sprintf ( msg , "Illegal KOHONEN LEARNING RATE: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         koh_params->rate = p ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN LEARNING REDUCTION" )) {
      p = atof ( rest ) ;
      if ((p < 0.0)  ||  (p > 1.0)) {
         sprintf ( msg , "Illegal KOHONEN LEARNING REDUCTION: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         koh_params->reduction = p ;
      return 1 ;
      }

   return 0 ;
}

static int check_learn_params ( char *command , char *rest ,
                                struct LearnParams *learn_params , int netmod )
{
   int m, n ;
   double p ;
   char msg[80] ;

   if (! strcmp ( command , "QUIT ERROR" )) {
      p = atof ( rest ) / 100.0 ; // User expresses it as percent
      if ((p < 0.0)  ||  (p > 1.0)) {
         sprintf ( msg , "Illegal QUIT ERROR: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         learn_params->quit_err = p ;
      return 1 ;
      }

   if (! strcmp ( command , "QUIT RETRIES" )) {
      m = sscanf ( rest , "%d" , &n ) ;
      if ((m <= 0)  ||  (n < 0)) {
         sprintf ( msg , "Illegal QUIT RETRIES: %s", rest ) ;
         error_message ( msg ) ;
         }
      else
         learn_params->retries = n ;
      return 1 ;
      }

   if (! strcmp ( command , "LAYER INIT" )) {
      if (netmod != NETMOD_LAYER) {
         error_message ( "LAYER INIT makes no sense for this model." ) ;
         return 1 ;
         }
      if ((n = get_layer_init ( rest )) < 0) {
         error_message ( "Illegal LAYER INIT method." ) ;
         return 1 ;
         }
      learn_params->init = n ;
      return 1 ;
      }

   if (! strcmp ( command , "KOHONEN INIT" )) {
      if (netmod != NETMOD_KOH) {
         error_message ( "KOHONEN INIT makes no sense for this model." ) ;
         return 1 ;
         }
      if ((n = get_koh_init ( rest )) < 0) {
         error_message ( "Illegal KOHONEN INIT method." ) ;
         return 1 ;
         }
      learn_params->init = n ;
      return 1 ;
      }

   return 0 ;
}
