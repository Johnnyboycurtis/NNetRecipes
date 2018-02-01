/******************************************************************************/
/*                                                                            */
/*  LAYERNET - All principal routines for LayerNet processing                 */
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

static void free_non_null ( void **p ) ;

/*
--------------------------------------------------------------------------------

   Constructor

   The parameter 'executable' determines whether work areas for hidden and
   output neuron activations are also allocated.  These are needed if we
   will ever apply inputs and want to compute outputs.

   In case of malloc failure, we set 'ok' to zero so the user knows about it.
   Also, we always leave unallocated weight pointers set to NULL.  There is no
   hard reason for doing this; calling programs should always know enough not
   to reference them.  However, it is simply good style.  Most compilers are
   much better at producing code that intercepts NULL pointer references than
   just wild pointers.  An ounce of prevention...

--------------------------------------------------------------------------------
*/

LayerNet::LayerNet (
   int out_model ,
   int n_inputs ,
   int n_hidden1 ,
   int n_hidden2 ,
   int n_outputs ,
   int executable ,  // Also allocate hidden and output neurons?
   int zero          // Zero all weights?
   )
{
   int n1, n2, n3 ;

   outmod = out_model ;
   nin = n_inputs ;
   nhid1 = n_hidden1 ;
   nhid2 = n_hidden2 ;
   nout = n_outputs ;
   exe = executable ;
   neterr = 1.0 ;

   confusion = NULL ;
   hid1_coefs = hid2_coefs = out_coefs = hid1 = hid2 = out = NULL ;

   ok = 0 ;   // Indicates failure of malloc (What a pessimist!)

   if (exe  &&  (confusion = (int *) MALLOC ( (nout+1) * sizeof(int))) == NULL)
      return ;

   if (nhid1 == 0) {                // No hidden layer
      n1 = nout * (nin+1) ;
      if (((out_coefs = (double *) MALLOC ( n1 * sizeof(double) )) == NULL)
       || (exe && (out = (double *) MALLOC ( nout * sizeof(double) )) == NULL)){
         free_non_null ( (void **) &out_coefs ) ;
         free_non_null ( (void **) &confusion ) ;
         return ;
         }
      if (zero) {
         while (n1--)
            out_coefs[n1] = 0.0 ;
         }
      }

   else if (nhid2 == 0) {           // One hidden layer
      n1 = nhid1 * (nin+1) ;
      n2 = nout * (nhid1+1) ;
      if (((hid1_coefs = (double *) MALLOC ( n1 * sizeof(double) )) == NULL)
       || ((out_coefs = (double *) MALLOC ( n2 * sizeof(double) ))==NULL)
       || (exe && (hid1 = (double *) MALLOC ( nhid1 * sizeof(double) ))==NULL)
       || (exe && (out = (double *) MALLOC ( nout * sizeof(double) )) == NULL)){
         free_non_null ( (void **) &hid1_coefs ) ;
         free_non_null ( (void **) &out_coefs ) ;
         free_non_null ( (void **) &hid1 ) ;
         free_non_null ( (void **) &confusion ) ;
         return ;
         }
      if (zero) {
         while (n1--)
            hid1_coefs[n1] = 0.0 ;
         while (n2--)
            out_coefs[n2] = 0.0 ;
         }
      }

   else {                           // Two hidden layers
      n1 = nhid1 * (nin+1) ;
      n2 = nhid2 * (nhid1+1) ;
      n3 = nout * (nhid2+1) ;
      if (((hid1_coefs = (double *) MALLOC ( n1 * sizeof(double) )) == NULL)
       || ((hid2_coefs = (double *) MALLOC ( n2 * sizeof(double) )) == NULL)
       || ((out_coefs = (double *) MALLOC ( n3 * sizeof(double) ))==NULL)
       || (exe && (hid1 = (double *) MALLOC ( nhid1 * sizeof(double) ))==NULL)
       || (exe && (hid2 = (double *) MALLOC ( nhid2 * sizeof(double) ))==NULL)
       || (exe && (out = (double *) MALLOC ( nout * sizeof(double) )) == NULL)){
         free_non_null ( (void **) &hid1_coefs ) ;
         free_non_null ( (void **) &hid2_coefs ) ;
         free_non_null ( (void **) &out_coefs ) ;
         free_non_null ( (void **) &hid1 ) ;
         free_non_null ( (void **) &hid2 ) ;
         free_non_null ( (void **) &confusion ) ;
         return ;
         }
      if (zero) {
         while (n1--)
            hid1_coefs[n1] = 0.0 ;
         while (n2--)
            hid2_coefs[n2] = 0.0 ;
         while (n3--)
            out_coefs[n3] = 0.0 ;
         }
      }

   if (exe)
      memset ( confusion , 0 , (nout+1) * sizeof(int) ) ;

   ok = 1 ;            // Indicate to caller that all mallocs succeeded
}

/*
   Local routine to free non-null pointers
*/

static void free_non_null ( void **p )
{
   if (*p != NULL) {
      FREE ( *p ) ;
      *p = NULL ;
      }
}


/*
--------------------------------------------------------------------------------

   Destructor

--------------------------------------------------------------------------------
*/

LayerNet::~LayerNet()
{
   if (! ok)    // If constructor's mallocs failed
      return ;  // there is nothing to free

   FREE ( out_coefs ) ;
   if (exe) {
      FREE ( out ) ;
      FREE ( confusion ) ;
      }

   if (nhid1) {
      FREE ( hid1_coefs ) ;
      if (exe)
         FREE ( hid1 ) ;
      if (nhid2) {
         FREE ( hid2_coefs ) ;
         if (exe)
            FREE ( hid2 ) ;
         }
      }
}

/*
--------------------------------------------------------------------------------

   copy_weights - Copy the weights from one network to another
                  Note that this is NOT like a copy or assignment,
                  as it does not copy other parameters.  In fact,
                  it gets sizes from the calling instance!

--------------------------------------------------------------------------------
*/

void LayerNet::copy_weights ( LayerNet *dest , LayerNet *source )
{
   int n ;

   dest->neterr = source->neterr ;
   if (source->exe  &&  dest->exe) // These may be important too!
      memcpy ( dest->confusion , source->confusion , (nout+1) * sizeof(int) ) ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      memcpy ( dest->out_coefs , source->out_coefs , n * sizeof(double) ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      memcpy ( dest->hid1_coefs , source->hid1_coefs , n * sizeof(double) ) ;
      n = nout * (nhid1+1) ;
      memcpy ( dest->out_coefs , source->out_coefs , n * sizeof(double) ) ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      memcpy ( dest->hid1_coefs , source->hid1_coefs , n * sizeof(double) ) ;
      n = nhid2 * (nhid1+1) ;
      memcpy ( dest->hid2_coefs , source->hid2_coefs , n * sizeof(double) ) ;
      n = nout * (nhid2+1) ;
      memcpy ( dest->out_coefs , source->out_coefs , n * sizeof(double) ) ;
      }
}


/*
--------------------------------------------------------------------------------

   zero_weights - Zero all weights in a network

--------------------------------------------------------------------------------
*/

void LayerNet::zero_weights ()
{
   int n ;

   neterr = 1.0 ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      while (n--)
         out_coefs[n] = 0.0 ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      while (n--)
         hid1_coefs[n] = 0.0 ;
      n = nout * (nhid1+1) ;
      while (n--)
         out_coefs[n] = 0.0 ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      while (n--)
         hid1_coefs[n] = 0.0 ;
      n = nhid2 * (nhid1+1) ;
      while (n--)
         hid2_coefs[n] = 0.0 ;
      n = nout * (nhid2+1) ;
      while (n--)
         out_coefs[n] = 0.0 ;
      }
}


/*
--------------------------------------------------------------------------------

   trial - Compute the output for a given input by evaluating network

--------------------------------------------------------------------------------
*/

void LayerNet::trial ( double *input )
{
   int i ;
   
   if (! exe) {   // Should NEVER happen, but good style to aid debugging
      error_message ( "Internal error in LayerNet::trial" ) ;
      return ;
      }

   if (nhid1 == 0) {                // No hidden layer
      for (i=0 ; i<nout ; i++)
         out[i] = activity ( input , out_coefs+i*(nin+1) , nin ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      for (i=0 ; i<nhid1 ; i++)
         hid1[i] = activity ( input , hid1_coefs+i*(nin+1) , nin ) ;
      for (i=0 ; i<nout ; i++)
         out[i] = activity ( hid1 , out_coefs+i*(nhid1+1) , nhid1 ) ;
      }

   else {                           // Two hidden layers
      for (i=0 ; i<nhid1 ; i++)
         hid1[i] = activity ( input , hid1_coefs+i*(nin+1) , nin ) ;
      for (i=0 ; i<nhid2 ; i++)
         hid2[i] = activity ( hid1 , hid2_coefs+i*(nhid1+1) , nhid1 ) ;
      for (i=0 ; i<nout ; i++)
         out[i] = activity ( hid2 , out_coefs+i*(nhid2+1) , nhid2 ) ;
      }
}


/*
--------------------------------------------------------------------------------

   learn

--------------------------------------------------------------------------------
*/

void LayerNet::learn ( TrainingSet *tptr , struct LearnParams *lptr )
{
   int i, itry, n_escape, n_retry, bad_count ;
   double err, prev_err, best_err, start_of_loop_error ;
   char msg[80] ;
   SingularValueDecomp *sptr ;
   LayerNet *worknet, *bestnet ;

   if (! exe) {   // Should NEVER happen, but good style to aid debugging
      error_message ( "Internal error in LayerNet::learn" ) ;
      return ;
      }

   n_escape = n_retry = 0 ;

/*
   Allocate scratch memory
*/

   MEMTEXT ( "LAYERNET::learn new worknet, bestnet" ) ;
   worknet = new LayerNet ( outmod , nin , nhid1 , nhid2 , nout , 0 , 0 ) ;
   bestnet = new LayerNet ( outmod , nin , nhid1 , nhid2 , nout , 0 , 1 ) ;

   if ((worknet == NULL)  ||  (! worknet->ok)
    || (bestnet == NULL)  ||  (! bestnet->ok)) {
      memory_message ( "to learn" ) ;
      if (worknet != NULL)
         delete worknet ;
      if (bestnet != NULL)
         delete bestnet ;
      neterr = 1.0 ;
      return ;
      }

/*
   Find initial weight estimates via simulated annealing or genetics
   or simply regression if no hidden layers
*/

   if (lptr->init == 0)                             // NO INIT
      neterr = trial_error ( tptr ) ;

   if ((lptr->init == 1)  ||  (lptr->init == 2))    // ANNEAL
      anneal ( tptr , lptr , worknet , 1 ) ;

   else if (lptr->init == 3) {                      // GENETIC
      if (nhid1)
         gen_init ( tptr , lptr ) ;
      else {
         error_message ( "Genetic init pointless with no hidden layer" ) ;
         neterr = 1.0 ;
         goto FINISH ;
         }
      }

   else if (lptr->init == 4) {                      // REGRESSION
      if (nhid1) {
         warning_message ( "REGRESS illegal if hidden layer.");
         neterr = 1.0 ;
         goto FINISH ;
         }

      MEMTEXT ( "LAYERNET: new SingularValueDecomp" ) ;
      sptr = new SingularValueDecomp ( tptr->ntrain , nin+1 , 1 ) ;
      if (! sptr->ok) {    // Indicates insufficient memory
         memory_message("for regression. Try ANNEAL NOREGRESS." ) ;
         neterr = 1.0 ;  // Flag failure to LayerNet::learn which called us
         delete sptr ;
         }
      neterr = regress ( tptr , sptr ) ;
      MEMTEXT ( "LAYERNET: delete SingularValueDecomp" ) ;
      delete sptr ;
      }

   if (lptr->init  &&  neterr > 0.999999) // Memory allocation failure
      goto FINISH ;

/*
   Initialization is done.  Learning loop is here.
   First, do conjugate gradient optimization, finding local minimum.
   Then anneal to break out of it.  If successful, loop back up to
   do conjugate gradient again.  Otherwise restart totally random.
*/

   copy_weights ( bestnet , this ) ;
   best_err = neterr ;
   bad_count = 0 ;           // Handles flat local mins

   for (itry=1 ; ; itry++) {

      sprintf ( msg , "Try %d  (best=%lf):", itry, 100.0 * best_err ) ;
      normal_message ( msg ) ;

      start_of_loop_error = neterr ;
      err = conjgrad ( tptr , 32767 , 1.e-8 , lptr->quit_err ) ;
      neterr = fabs ( err ) ; // err<0 if user pressed ESCape

      sprintf ( msg , "  Gradient err=%lf", 100.0 * neterr ) ;
      progress_message ( msg ) ;

      if (neterr < best_err) {   // Keep track of best
         copy_weights ( bestnet , this ) ;
         best_err = neterr ;
         }

      if (err < lptr->quit_err) // err<0 if user pressed ESCape
         break ;

      i = itry * 97 + 101 ;   // Insure new seed for anneal
      if (i<0)
         i = -i ;
      slongrand ( (long) i ) ;

      prev_err = neterr ;  // So we can see if anneal helped
      anneal ( tptr , lptr , worknet , 0 ) ;

      sprintf ( msg , "  Anneal err=%lf", 100.0 * neterr ) ;
      progress_message ( msg ) ;

      if (neterr < best_err) {  // Keep track of best
         copy_weights ( bestnet , this ) ;
         best_err = neterr ;
         }

      if (best_err < lptr->quit_err)
         break ;

      if (neterr < prev_err) { // Did we break out of local min?
         if ((start_of_loop_error - neterr) < 1.e-3)
            ++bad_count ;  // Avoid many unprofitable iters
         else
            bad_count = 0 ;
         if (bad_count < 4) {
            ++n_escape ;          // For user interest only
            continue ;            // Escaped, so gradient learn again
            }
         }

      if (++n_retry > lptr->retries)
         break ;

      progress_message ( "  RESTART" ) ;
      zero_weights () ;  // Failed to break out, so retry random
      anneal ( tptr , lptr , worknet , 1 ) ;
      }

FINISH:
   copy_weights ( this , bestnet ) ;
   MEMTEXT ( "LAYERNET::learn delete worknet, bestnet" ) ;
   delete worknet ;
   delete bestnet ;
   sprintf ( msg , "%d successful escapes, %d retries", n_escape, n_retry ) ;
   normal_message ( msg ) ;

   return ;
}


/*
--------------------------------------------------------------------------------

   wt_save - Save weights to disk (called from WT_SAVE.CPP)
   wt_restore - Restore weights from disk (called from WT_SAVE.CPP)

--------------------------------------------------------------------------------
*/

int LayerNet::wt_save ( FILE *fp )
{
   int n ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      fwrite ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      fwrite ( hid1_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nout * (nhid1+1) ;
      fwrite ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      fwrite ( hid1_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nhid2 * (nhid1+1) ;
      fwrite ( hid2_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nout * (nhid2+1) ;
      fwrite ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   if (ferror ( fp ))
      return 1 ;
   return 0 ;
}

void LayerNet::wt_restore ( FILE *fp )
{
   int n ;

   if (nhid1 == 0) {                // No hidden layer
      n = nout * (nin+1) ;
      fread ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   else if (nhid2 == 0) {           // One hidden layer
      n = nhid1 * (nin+1) ;
      fread ( hid1_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nout * (nhid1+1) ;
      fread ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   else {                           // Two hidden layers
      n = nhid1 * (nin+1) ;
      fread ( hid1_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nhid2 * (nhid1+1) ;
      fread ( hid2_coefs , n * sizeof(double) , 1 , fp ) ;
      n = nout * (nhid2+1) ;
      fread ( out_coefs , n * sizeof(double) , 1 , fp ) ;
      }

   if (ferror ( fp ))
      ok = 0 ;
}
