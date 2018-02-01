/******************************************************************************/
/*                                                                            */
/*  KOHNET - All principal routines for KohNet processing                     */
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

   Note that some normalization methods generate an extra input.
   Therefore we always allocate nin+1 length vectors, even though
   we may not need the extra weight.

   The parameter 'executable' determines whether work areas for
   output neuron activations are also allocated. These are needed
   if we will ever apply inputs and want to compute outputs.

   In case of malloc failure, we set 'ok' to zero so the user knows about it.
   Also, we always leave unallocated pointers set to NULL.  There is no
   hard reason for doing this; calling programs should always know enough not
   to reference them.  However, it is simply good style.  Most compilers are
   much better at producing code that intercepts NULL pointer references than
   just wild pointers.  An ounce of prevention...

--------------------------------------------------------------------------------
*/

KohNet::KohNet (
   int n_inputs ,
   int n_outputs ,
   KohParams *kp ,  // Specialized parameters
   int executable , // Also allocate hidden and output neurons?
   int zero         // Zero all weights?
   )
{
   int n ;

   outmod = OUTMOD_CLASSIFY ;
   nin = n_inputs ;
   nout = n_outputs ;
   normalization = kp->normalization ;
   exe = executable ;
   neterr = 1.0 ;

   confusion = NULL ;
   out_coefs = out = NULL ;

   ok = 0 ;   // Indicates failure of malloc (What a pessimist!)

   if (exe && (confusion=(int *) MALLOC((nout+1) * sizeof(int))) == NULL)
      return ;

   n = nout * (nin+1) ; // Some normalizations generate extra input
   if (((out_coefs = (double *) MALLOC ( n * sizeof(double) )) == NULL)
    || (exe && (out = (double *) MALLOC ( nout * sizeof(double) )) == NULL)){
      free_non_null ( (void **) &out_coefs ) ;
      free_non_null ( (void **) &confusion ) ;
      return ;
      }

   if (zero) {
      while (n--)
         out_coefs[n] = 0.0 ;
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

KohNet::~KohNet()
{
   if (! ok)    // If constructor's mallocs failed
      return ;  // there is nothing to free

   FREE ( out_coefs ) ;
   if (exe) {
      FREE ( out ) ;
      FREE ( confusion ) ;
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

void KohNet::copy_weights ( KohNet *dest , KohNet *source )
{
   int n ;

   dest->neterr = source->neterr ;
   if (source->exe  &&  dest->exe) // These may be important too!
      memcpy ( dest->confusion , source->confusion , (nout+1) * sizeof(int) ) ;

   n = nout * (nin+1) ;
   memcpy ( dest->out_coefs , source->out_coefs , n * sizeof(double) ) ;
}


/*
--------------------------------------------------------------------------------

   zero_weights - Zero all weights in a network

--------------------------------------------------------------------------------
*/

void KohNet::zero_weights ()
{
   int n ;

   neterr = 1.0 ;

   n = nout * (nin+1) ;
   while (n--)
      out_coefs[n] = 0.0 ;
}

/*
--------------------------------------------------------------------------------

   KOH_NORM - Routines for normalizing Kohonen vectors to unit length

   in_norm - Normalize an input vector by computing a normalizing
             factor and the synthetic last input.
             The input vector itself is not touched.
             It is assumed that all inputs are in the range -1 to 1.
             The end result is that if the inputs are multiplied by
             normfac, that vector with synth appended has unit length.

   wt_norm - Normalize a weight vector in place.  The synthetic last
             component is NOT computed.

--------------------------------------------------------------------------------
*/

void KohNet::in_norm (
   double *input ,   // Input vector
   double *normfac , // Output: multiply input by this
   double *synth     // Output: synthetic last input
   )
   
{
   double length, d ;

   length = veclen ( nin , input ) ; // Squared length
   if (length < 1.e-30)              // Safety
      length = 1.e-30 ;

   if (normalization == 0) {      // Multiplicative
      *normfac = 1.0 / sqrt ( length ) ;
      *synth = 0.0 ;
      }

   else if (normalization == 1) { // Z
      *normfac = 1.0 / sqrt ( nin ) ;
      d = (double) nin - length ;
      if (d > 0.0)
         *synth = sqrt ( d ) * *normfac ;
      else                // If the inputs are all -1 to 1
         *synth = 0.0 ;   // this error never occurs
      }
}

void KohNet::wt_norm ( double *w )
{
   int i ;
   double len ;

   len = veclen ( nin , w ) ;     // Ignore last weight
   if (len < 1.e-30)              // Safety
      len = 1.e-30 ;

   if (normalization == 0) {      // Multiplicative
      len = 1.0 / sqrt ( len ) ;
      for (i=0 ; i<nin ; i++)
         w[i] *= len ;
      w[nin] = 0. ;
      }

   else if (normalization == 1) { // Z
      len += w[nin] * w[nin] ;
      len = 1.0 / sqrt ( len ) ;
      for (i=0 ; i<=nin ; i++)
         w[i] *= len ;
      }
}

/*
--------------------------------------------------------------------------------

   trial - Compute the outputs for a given input by evaluating the network.
           It is assumed that all inputs are from -1 to 1, but not
           necessarily normalized (that is done here).

--------------------------------------------------------------------------------
*/

void KohNet::trial ( double *input )
{
   int i ;
   double normfac, synth, *optr ;
   
   in_norm ( input , &normfac , &synth ) ;  // Normalize input

   for (i=0 ; i<nout ; i++) {
      optr = out_coefs + i * (nin+1) ;  // i'th weight vector
      out[i] = dotprod ( nin , input , optr ) * normfac
               + synth * optr[nin] ;
      out[i] = 0.5 * (out[i] + 1.0) ;   // Remap -1,1 to 0,1
      if (out[i] > 1.0)   // Only trivially happens due to rounding
         out[i] = 1.0 ;
      if (out[i] < 0.0)
         out[i] = 0.0 ;
      }
}


/*
--------------------------------------------------------------------------------

   winner - Return the subscript of the winning neuron.
            This is identical to 'trial' above except that
            it also returns the normalization info and winner.

--------------------------------------------------------------------------------
*/

int KohNet::winner (
   double *input ,   // Input vector
   double *normfac , // Output: multiply input by this
   double *synth     // Output: synthetic last input
   )
{
   int i, win ;
   double biggest, *optr ;

   in_norm ( input , normfac , synth ) ;  // Normalize input

   biggest = -1.e30 ;
   for (i=0 ; i<nout ; i++) {
      optr = out_coefs + i * (nin+1) ;  // i'th weight vector
      out[i] = dotprod ( nin , input , optr ) * *normfac
               + *synth * optr[nin] ;
      out[i] = 0.5 * (out[i] + 1.0) ;   // Remap -1,1 to 0,1
      if (out[i] > biggest) {
         biggest = out[i] ;
         win = i ;
         }
      if (out[i] > 1.0)   // Only trivially happens due to rounding
         out[i] = 1.0 ;
      if (out[i] < 0.0)
         out[i] = 0.0 ;
      }

   return win ;
}


/*
--------------------------------------------------------------------------------

   learn

--------------------------------------------------------------------------------
*/

void KohNet::learn (
   TrainingSet *tptr ,       // Training set
   struct LearnParams *lptr  // Learning parameters
   )

{
   int i, key, tset ;
   int iter ;         // Iterations (epochs)
   int n_retry ;      // Number of random retries
   int nwts ;         // Total number of weights
   int *won ;         // Counts how many times each neuron won
   int winners ;      // How many neurons won per epoch
   char msg[80] ;     // For messages to user
   double *work ;     // Scratch for additive learning
   double *correc ;   // Scratch for cumulative correction vector
   double rate ;      // Current learning rate
   double bigerr ;    // Biggest single error in epoch
   double bigcorr ;   // Biggest cumulative correction in epoch
   double best_err ;  // Minimum error so far
   double *dptr ;     // Points to a training case
   KohNet *bestnet ;  // Preserve best here
   KohParams *kp ;    // User's parameters here

   if (! exe) {   // Should NEVER happen, but good style to aid debugging
      error_message ( "Internal error in KohNet::learn" ) ;
      exit ( 1 ) ;
      }

   kp = lptr->kp ;  // Simplify pointing to parameters
   neterr = 1.0 ;

/*
   If this is multiplicative normalization, make sure all training
   cases are non-null.
*/

   if (normalization == 0) {      // Multiplicative
      for (tset=0 ; tset<tptr->ntrain ; tset++) {
         dptr = tptr->data + (nin+1) * tset ;
         if (veclen ( nin , dptr ) < 1.e-30) {
            error_message (
                      "Multiplicative normalization has null training case" ) ;
            return ;
            }
         }
      }

/*
   Allocate scratch memory, then initialize weights
*/

   MEMTEXT ( "KOHNET::learn new bestnet" ) ;
   bestnet = new KohNet ( nin , nout , kp , 0 , 1 ) ;

   if ((bestnet == NULL)  ||  (! bestnet->ok)) {
      memory_message ( "to learn" ) ;
      if (bestnet != NULL)
         delete bestnet ;
      return ;
      }

   nwts = nout * (nin+1) ;
   MEMTEXT ( "KOHNET: Learn scratch" ) ;
   won = (int *) MALLOC ( nout * sizeof(int) ) ;
   correc = (double *) MALLOC ( nwts * sizeof(double) ) ;
   if (! kp->learn_method)  // Needed only for additive method
      work = (double *) MALLOC ( (nin+1) * sizeof(double)) ;
   else 
      work = NULL ;

   if ((won == NULL)  ||  (correc == NULL)  ||
       ((! kp->learn_method)  &&  (work == NULL))) {
      if (won != NULL)
         FREE ( won ) ;
      if (correc != NULL)
         FREE ( correc ) ;
      if (work != NULL)
         FREE ( work ) ;
      delete bestnet ;
      memory_message ( "to learn" ) ;
      return ;
      }

   rate = kp->rate ;

/*
   If the user specified NOINIT, they are continuing to learn from
   existing weights.  Call epoch1 to find the error associated with
   those weights, and save the weights as best so far.
   Then don't waste that call to epoch1.  Call epoch2 to update
   the weights.
*/

   if (lptr->init == 0) {      // NOINIT (continue learning)
      epoch1 ( tptr , rate , 1 , won , &bigerr , correc , work ) ;
      best_err = neterr = bigerr ;
      copy_weights ( bestnet , this ) ;
      epoch2 ( rate , kp->learn_method , won , &bigcorr , correc ) ;
      }
   else if (lptr->init == 1) { // RANDOM Initialize weights
      initialize () ;
      best_err = 1.e30 ;
      }

/*
   Main loop is here.  Each iter is a complete epoch.
*/

   n_retry = 0 ;
   for (iter=0 ; ; iter++) {

      epoch1 ( tptr , rate , kp->learn_method , won ,
               &bigerr , correc , work ) ;

      neterr = bigerr ;

      if (neterr < best_err) {  // Keep track of best
         best_err = neterr ;
         copy_weights ( bestnet , this ) ;
         }

      winners = 0 ;     // Count how many neurons won this iter
      i = nout ;
      while (i--) {
         if (won[i])
            ++winners ;
         }

      sprintf( msg ,
         "Iter %d err=%.2lf (best=%.2lf)  %d won",
          iter, 100.0 * neterr, 100.0 * best_err, winners ) ;
      normal_message ( msg ) ;

      if (kbhit()) {            // Was a key pressed?
         key = getch () ;       // Read it if so
         while (kbhit())        // Flush key buffer in case function key
            getch () ;          // or key was held down
         if (key == 27)         // ESCape
            break ;
         }

      if (bigerr < lptr->quit_err) // Are we done?
         break ;

/*
   If one or more neurons failed to ever win, make it a winner.
   Note that this has a theoretical flaw.
   If the training set has duplication such that there are fewer
   unique values than neurons, we can get in a loop of flipping
   case values around neurons.  Thus, rather than verifying
   winners<tptr->ntrain below, we should ideally count how many
   unique values are in the training set, and use that number.
   However, that would be time consuming and protect against an
   extremely unlikely event anyway.
*/

      if ((winners < nout)  &&  (winners < tptr->ntrain)) {
         force_win ( tptr , won ) ;
         continue ;
         }

      epoch2 ( rate , kp->learn_method , won , &bigcorr , correc ) ;

      sprintf( msg , "  correction=%.2lf", 100.0 * bigcorr ) ;
      progress_message ( msg ) ;

      if (bigcorr < 1.e-5) { // Trivial improvement?
         if (++n_retry > lptr->retries) // If so, start over
            break ;          // unless at user's limit
         initialize () ;     // Use totally random weights
         iter = -1 ;         // End of loop incs this to 0
         rate = kp->rate ;   // Rate starts high again
         continue ;
         }

      if (rate > 0.01)  // Reduce learning rate each time
         rate *= kp->reduction ;

      } // Endless learning loop

/*
   We are done.  Retrieve the best weights.  Learning should have left
   them very close to normalized, but it doesn't hurt to touch them up.
   Unfortunately, this can slightly change the network error.
*/

   copy_weights ( this , bestnet ) ;

   for (i=0 ; i<nout ; i++)
      wt_norm ( out_coefs + i * (nin+1) ) ;

   MEMTEXT ( "KOHNET: Learn scratch" ) ;
   delete bestnet ;
   FREE ( won ) ;
   FREE ( correc ) ;
   if (! kp->learn_method)  // Needed only for additive method
      FREE ( work ) ;
   return ;
}

/*
--------------------------------------------------------------------------------

   initialize - Initialize weights

--------------------------------------------------------------------------------
*/

void KohNet::initialize ()
{
   int i ;
   double *optr ;

   zero_weights () ;
   shake ( nout * (nin+1) , out_coefs , out_coefs , 1.0 ) ;
   for (i=0 ; i<nout ; i++) {
      optr = out_coefs + i * (nin+1) ;  // This weight vector
      wt_norm ( optr ) ;
      }
}


/*
--------------------------------------------------------------------------------

   epoch1 - Compute the error and correction vector

--------------------------------------------------------------------------------
*/

void KohNet::epoch1 (
   TrainingSet *tptr , // Training set
   double rate ,       // Learning rate
   int learn_method ,  // 0=additive, 1=subtractive
   int *won ,          // Work vector holds times each neuron won
   double *bigerr ,    // Returns max error length across training set
   double *correc ,    // Work vector nout*(nin+1) long for corrections
   double *work        // Work vector nin+1 long for additive learning
   )   

{
   int i, best, size, nwts, tset ;
   double *dptr, normfac, synth, *cptr, *wptr, length, diff ;

   nwts = nout * (nin+1) ;
   size = nin + 1 ;   // Size of each case in training set

/*
   Zero cumulative corrections and winner counts
*/

   i = nwts ;
   while (i--)
      correc[i] = 0.0 ;

   memset ( won , 0 , nout * sizeof(int) ) ;
   *bigerr = 0.0 ;  // Length of biggest error vector

/*
   Cumulate the correction vector 'correc' across the epoch
*/

   for (tset=0 ; tset<tptr->ntrain ; tset++) {
      dptr = tptr->data + size * tset ; // Point to this case
      best = winner ( dptr , &normfac , &synth ) ; // Winning neuron
      ++won[best] ;                   // Record this win
      wptr = out_coefs+best*(nin+1) ; // Winner's weights here
      cptr = correc+best*(nin+1) ;    // Corrections summed here
      length = 0.0 ;                  // Length of error vector

      for (i=0 ; i<nin ; i++) {  // Do all inputs
         diff = dptr[i] * normfac - wptr[i] ; // Input minus weight
         length += diff * diff ; // Cumulate length of error
         if (learn_method)       // Subtractive method
            cptr[i] += diff ;    // just uses differences
         else                    // Additive more complex
            work[i] = rate * dptr[i] * normfac + wptr[i] ;
         }                       // Loop does actual inputs
      diff = synth - wptr[nin] ; // Don't forget synthetic input
      length += diff * diff ;    // It is part of input too!
      if (learn_method)          // Subtractive method
         cptr[nin] += diff ;     // Cumulate across epoch
      else                       // Additive more complex
         work[nin] = rate * synth + wptr[nin] ;

      if (length > *bigerr)      // Keep track of largest error
         *bigerr = length ;

      if (! learn_method) {      // Additive method
         wt_norm ( work ) ;
         for (i=0 ; i<=nin ; i++)
            cptr[i] += work[i] - wptr[i] ;
         }

      } // Pass through all training sets, cumulating correction vector

   *bigerr = sqrt ( *bigerr ) ;
}


/*
--------------------------------------------------------------------------------

   epoch2 - Adjust weights per corrections from epoch1

--------------------------------------------------------------------------------
*/

void KohNet::epoch2 (
   double rate ,       // Learning rate
   int learn_method ,  // 0=additive, 1=subtractive
   int *won ,          // Work vector holds times each neuron won
   double *bigcorr ,   // Returns length of largest correction vector
   double *correc      // Work vector nout*(nin+1) long for corrections
   )   

{
   int i, j ;
   double corr, *cptr, *wptr, length, f ;

   *bigcorr = 0.0 ;                // Length of largest correction

   for (i=0 ; i<nout ; i++) {      // Apply mean correction to each

      if (! won[i])                // If this neuron never won
         continue ;                // might as well skip update

      wptr = out_coefs+i*(nin+1) ; // i's weights here
      cptr = correc+i*(nin+1) ;    // Corrections were summed here

      f = 1.0 / (double) won[i] ;  // Finds mean across epoch
      if (learn_method)            // Subtractive method
         f *= rate ;               // needs learning rate included

      length = 0.0 ;               // Will sum length of correction

      for (j=0 ; j<=nin ; j++) {   // Weight vector for this neuron
         corr = f * cptr[j] ;      // Mean correction
         wptr[j] += corr ;         // Update weight vector
         length += corr * corr ;   // Sum length of this correction
         }

      if (length > *bigcorr)       // Keep track of biggest correction
         *bigcorr = length ;
      }

/*
   Scale the correction length per learning rate so that we
   are not fooled into thinking we converged when really all
   that happened is that the learning rate got small.
   Note that it can exceed 1.0 if the weights and data
   pointed in opposing directions.
*/

   *bigcorr = sqrt ( *bigcorr ) / rate ;
}


/*
--------------------------------------------------------------------------------

   force_win - Force a neuron to win.

--------------------------------------------------------------------------------
*/

void KohNet::force_win (
   TrainingSet *tptr , // Training set
   int *won            // Work vector holds times each neuron won
   )   

{
   int i, tset, best, size, which ;
   double *dptr, normfac, synth, dist, *optr ;

   size = nin + 1 ;  // Size of each training case

/*
   Find the training case which is farthest from its winning neuron.
   It is reasonable to believe that this case is not adequately
   represented by that neuron, and deserves a neuron of its very own.
*/

   dist = 1.e30 ;
   for (tset=0 ; tset<tptr->ntrain ; tset++) {
      dptr = tptr->data + size * tset ; // Point to this case
      best = winner ( dptr , &normfac , &synth ) ; // Winning neuron
      if (out[best] < dist) {  // Far indicated by low activation
         dist = out[best] ;    // Maintain record
         which = tset ;        // and which case did it
         }
      }

/*
   Now find the non-winning neuron which is most similar to
   the under-represented case found above.
*/

   dptr = tptr->data + size * which ;
   best = winner ( dptr , &normfac , &synth ) ;

   dist = -1.e30 ;
   i = nout ;
   while (i--) {           // Try all neurons
      if (won[i])          // If this one won then skip it
         continue ;        // We want a non-winner
      if (out[i] > dist) { // High activation means similar
         dist = out[i] ;   // Keep track of best
         which = i ;       // and its subscript
         }
      }

/*
   Use that training case to define the new weights.
   Strictly speaking, we should multiply the inputs by normfac,
   then append synth.  But since we normalize, it is equivalent
   (and faster) to copy the inputs, then append synth / normfac.
*/

   optr = out_coefs + which * (nin+1) ;        // Non-winner's weights
   memcpy( optr , dptr , nin*sizeof(double)) ; // become case
   optr[nin] = synth / normfac ;               // Append synth
   wt_norm ( optr ) ;                          // Keep normal
}


/*
--------------------------------------------------------------------------------

   wt_save - Save weights to disk (called from WT_SAVE.CPP)
   wt_restore - Restore weights from disk (called from WT_SAVE.CPP)

--------------------------------------------------------------------------------
*/

int KohNet::wt_save ( FILE *fp )
{
   int n ;

   n = nout * (nin+1) ;
   fwrite ( out_coefs , n * sizeof(double) , 1 , fp ) ;
   if (ferror ( fp ))
      return 1 ;
   return 0 ;
}

void KohNet::wt_restore ( FILE *fp )
{
   int n ;

   n = nout * (nin+1) ;
   fread ( out_coefs , n * sizeof(double) , 1 , fp ) ;
   if (ferror ( fp ))
      ok = 0 ;
}

