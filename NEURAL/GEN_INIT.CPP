/******************************************************************************/
/*                                                                            */
/*  GEN_INIT - Use genetic method to initialize LayerNet weights              */
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
   Declarations for local subroutines
*/

static void decode( char *popptr , int nin , int nh1 , int nh2 ,
                    double *w1 , double *w2);
static void error_to_fitness ( int popsize , double favor_best ,
            double fitfac , double *errors , double *fitness ) ;
static void fitness_to_choices ( int popsize , double *fitness , int *choices );
static void mutate ( char *child , int nvars , double pmutate );
static void pick_parents ( int *nchoices , int *choices ,
                           int *parent1 , int *parent2 ) ;
static void rand_ind ( char *popptr , int chromsize ) ;
static void reproduce ( char *p1 , char *p2 , int first_child ,
                        int nvars , char *child , int *crosspt , int *split ) ;

/*
   For translation speed we convert from gray codes to binary with a lookup
   table which is built in the first call.
*/

static unsigned char gray_code_table[256] ;  // Translation table
static int gray_initialized = 0 ;            // Has it been built yet?

/*
   Entry point is here
*/

void LayerNet::gen_init (
   TrainingSet *tptr ,        // Training set to use
   struct LearnParams *lptr   // User's general learning parameters
   )
{
   int i, istart, individual, best_individual, generation, n_cross ;
   int first_child, parent1, parent2, improved, crosspt, nchoices, *choices ;
   int initpop, popsize, gens, climb, nvars, chromsize, split, ind ;
   double pcross, pmutate, error, besterror, *errors, *fitness, worst ;
   double fquit, favor_best, fitfac, maxerr, minerr, avgerr, overinit ;
   SingularValueDecomp *sptr ;
   struct GenInitParams *gptr ;  // User's genetic initialization parameters
   char *pool1, *pool2, *oldpop, *newpop, *popptr, *temppop, *best ;
   char msg[80] ;

   gptr = lptr->gp ;

   popsize = gptr->pool ;
   gens = gptr->gens ;
   climb = gptr->climb ;
   overinit = gptr->overinit ;
   pcross = gptr->pcross ;
   pmutate = gptr->pmutate ;

   fquit = lptr->quit_err ;

   favor_best = 3.1 ;
   fitfac = -20.0 ;


/*
--------------------------------------------------------------------------------

   Do all scratch memory allocation.

--------------------------------------------------------------------------------
*/

/*
   Allocate the singular value decomposition object for REGRESS.
*/

   if (nhid2 == 0)         // One hidden layer
      nvars = nhid1 + 1 ;
   else                    // Two hidden layers
      nvars = nhid2 + 1 ;

   MEMTEXT ( "GEN_INIT: new SingularValueDecomp" ) ;
   sptr = new SingularValueDecomp ( tptr->ntrain , nvars , 1 ) ;

   if ((sptr == NULL)  || ! sptr->ok) {
      memory_message("for genetic initialization. Try ANNEAL NOREGRESS.");
      neterr = 1.0 ;    // Flag failure to LayerNet::learn which called us
      if (sptr != NULL)
         delete sptr ;
      return ;
      }

   chromsize = nhid1 * (nin+1) ;        // Length of an individual's chromosome
   if (nhid2)                           // is the number of hidden weights
      chromsize += nhid2 * (nhid1+1) ;

   errors = fitness = NULL ;
   choices = NULL ;
   pool1 = pool2 = NULL ;
   MEMTEXT ( "GEN_INIT: errors, fitness, choices, best, pool1,pool2");
   if (((errors = (double*) MALLOC ( popsize * sizeof(double))) == NULL)
    || ((fitness = (double*) MALLOC ( popsize * sizeof(double))) == NULL)
    || ((best = (char*) MALLOC( chromsize )) == NULL)
    || ((choices = (int*) MALLOC ( popsize * sizeof(int))) == NULL)
    || ((pool1 = (char*) MALLOC( popsize * chromsize )) == NULL)
    || ((pool2 = (char*) MALLOC( popsize * chromsize )) == NULL)) {
      if (errors != NULL)
         FREE ( errors ) ;
      if (fitness != NULL)
         FREE ( fitness ) ;
      if (choices != NULL)
         FREE ( choices ) ;
      if (pool1 != NULL)
         FREE ( pool1 ) ;
      if (pool2 != NULL)
         FREE ( pool2 ) ;
      delete sptr ;
      memory_message("for genetic initialization.  Try ANNEAL NOREGRESS." ) ;
      neterr = 1.0 ;  // Flag failure to LayerNet::learn which called us
      return ;
      }

/*
   Generate initial population pool.

   We also preserve the best weights across all generations,
   as this is what we will ultimately return to the user.
   Its mean square error is besterror.
*/

   besterror = 1.e30 ;     // For saving best (across all individuals and gens)
   maxerr = avgerr = 0.0 ; // For progress display only
   best_individual = 0 ;   // Safety only
   initpop = popsize * overinit ; // Overinitialization of initial population
   progress_message ( "\nGenerating initial population" ) ;

   for (ind=0 ; ind<initpop ; ind++) { // Try overinitialization times

      if (ind<popsize)                 // If still in pop size limit
         individual = ind ;            // just use next avail space
      else {                           // Else we search entire pop
         worst = -1. ;                 // for the worst member
         for (i=0 ; i<popsize ; i++) { // which we will then replace
            if (errors[i] > worst) {
               worst = errors[i] ;
               individual = i ;
               }
            }
         avgerr -= worst ;             // Exclude discards from average
         }

      popptr = pool1 + individual * chromsize ; // Build init pop in pool1
      rand_ind ( popptr , chromsize ) ;         // Randomly generate individual
      decode ( popptr , nin , nhid1 , nhid2 ,   // Convert genotype (chromosome)
               hid1_coefs , hid2_coefs );       // to phenotype (weights)
      error = regress ( tptr , sptr ) ;         // Evaluate network error
      errors[individual] = error ;              // and keep all errors

      if (error < besterror) {                  // Keep track of best
         besterror = error ;                    // as it is returned to user
         best_individual = individual ;         // This is its index in pool1
         }

      if (error > maxerr)                       // Max and average error are
         maxerr = error ;                       // for progress display only
      avgerr += error ;

      if (error <= fquit)
         break ;

      progress_message ( "." ) ;
      }

   sprintf (msg , "\nInitial pop:    Min err=%7.4lf   Max=%7.4lf   Avg=%7.4lf",
            100. * besterror, 100. * maxerr, 100.0 * avgerr / (double) popsize);
   progress_message ( msg ) ;


/*
   The initial population has been built in pool1.
   Copy its best member to 'best' in case it never gets beat (unlikely
   but possible!).
   Also, we will need best if the climb option is true.
*/

   popptr = pool1 + best_individual * chromsize ; // Point to best
   memcpy ( best , popptr , chromsize ) ;         // and save it

/*
   This is the main generation loop.  There are two areas for population pool
   storage: pool1 and pool2.  At any given time, oldpop will be set to one of
   them, and newpop to the other.  This avoids a lot of copying.
*/

   oldpop = pool1 ;       // This is the initial population
   newpop = pool2 ;       // The next generation is created here

   for (generation=0 ; generation<gens ; generation++) {

      if (error <= fquit) // We may have satisfied this in init pop
         break ;          // So we test at start of generation loop

      error_to_fitness ( popsize , favor_best , fitfac , errors , fitness ) ;

      fitness_to_choices ( popsize , fitness , choices ) ;

      nchoices = popsize ;         // Will count down as choices array emptied
      n_cross = pcross * popsize ; // Number crossing over
      first_child = 1 ;            // Generating first of parent's 2 children?
      improved = 0 ;               // Flags if we beat best

      if (climb) {                 // If we are to hill climb
         memcpy ( newpop , best , chromsize ) ; // start with best
         errors[0] = besterror ;   // Record its error
         istart = 1 ;              // and start children past it
         }
      else
         istart = 0 ;

/*
   Generate the children
*/

      maxerr = avgerr = 0.0 ;   // For progress display only
      minerr = 1.0 ;            // Ditto

      for (individual=istart ; individual<popsize ; individual++) {

         popptr = newpop + individual * chromsize ; // Will put this child here

         if (first_child)  // If this is the first of 2 children, pick parents
            pick_parents ( &nchoices , choices , &parent1 , &parent2 ) ;

         if (n_cross-- > 0)    // Do crossovers first
            reproduce ( oldpop + parent1 * chromsize , oldpop + parent2 * chromsize ,
                        first_child , chromsize , popptr , &crosspt , &split ) ;
         else if (first_child) // No more crossovers, so just copy parent
            memcpy ( popptr , oldpop + parent1 * chromsize , chromsize ) ;
         else
            memcpy ( popptr , oldpop + parent2 * chromsize , chromsize );

         if (pmutate > 0.0)
            mutate ( popptr , chromsize , pmutate ) ;

         decode ( popptr , nin , nhid1 , nhid2 , hid1_coefs , hid2_coefs ) ;
         error = regress ( tptr , sptr ) ; // Evaluate child's error
         errors[individual] = error ;      // and keep each

         if (error < besterror) {          // Keep track of best
            besterror = error ;            // It will be returned to user
            best_individual = individual ; // This is its index in newpop
            improved = 1 ;                 // Flag so we copy it later
            }

         if (error > maxerr)               // Min, max and average error
            maxerr = error ;               // for progress display only
         if (error < minerr)
            minerr = error ;
         avgerr += error ;

         if (error <= fquit)
            break ;

         first_child = ! first_child ;
         } // For all genes in population

/*
   We finished generating all children.  If we improved (one of these
   children beat the best so far) then copy that child to the best.
   Swap oldpop and newpop for the next generation.
*/

      if (improved) {
         popptr = newpop + best_individual * chromsize ; // Point to best
         memcpy ( best , popptr , chromsize ) ;          // and save it
         }

      temppop = oldpop ;   // Switch old and new pops for next generation
      oldpop = newpop ;
      newpop = temppop ;

      sprintf(msg, "\nGeneration %3d: Min err=%7.4lf   Max=%7.4lf   Avg=%7.4lf",
              generation+1, 100. * minerr, 100. * maxerr,
              100.0 * avgerr / (double) popsize ) ;
      progress_message ( msg ) ;
      }

/*
   We are all done.
*/

   decode ( best , nin , nhid1 , nhid2 , hid1_coefs , hid2_coefs ) ;
   besterror = regress ( tptr , sptr ) ;              // Evaluate network error

   MEMTEXT ( "GEN_INIT: errors, fitness, choices, best, pool1,pool2");
   FREE ( errors ) ;
   FREE ( fitness ) ;
   FREE ( choices ) ;
   FREE ( best ) ;
   FREE ( pool1 ) ;
   FREE ( pool2 ) ;
   MEMTEXT ( "GEN_INIT: delete sptr" ) ;
   delete sptr ;
}

/*
--------------------------------------------------------------------------------

   error_to_fitness - Convert the objective function value of each individual
                     to a scaled fitness value.  The scaled fitness may be
                     considered an expected frequency of choice.

--------------------------------------------------------------------------------
*/

static void error_to_fitness (
   int popsize ,       // Length of errors, fitness vectors
   double favor_best , // Factor for favoring best over average (2-3 is good)
   double fitfac ,     // Factor for converting error to raw fitness (-20 good)
   double *errors ,    // Input popsize vector of values of objective function
   double *fitness     // Output popsize vector of scaled fitnesses
   )
{
   int individual ;
   double fit, avgfitness, minfitness, maxfitness, ftemp, tmult, tconst ;

/*
   The first step is to convert the objective function value (which is to
   be minimized) into a raw (unscaled) fitness value.  The best method
   can be problem dependent.  Certainly, the mapping function must be
   decreasing, as we want smaller values of the objective function to map to
   larger values of fitness.  Also, later calculations are simplified if the
   fitness is always positive.

   The conversion function used here is f(v) = exp ( k * v ) where k is a
   negative number.  For objective functions which range from zero to one,
   as would be the case of a relative error function, a constant of about
   -20 is appropriate.  This would map .001 to .98, .01 to .82 and .1 to .14.
*/

   avgfitness = 0.0 ;
   maxfitness = -1.e30 ;
   minfitness = 1.e30 ;

   for (individual=0 ; individual<popsize ; individual++) {
      fitness[individual] = fit = exp ( fitfac * errors[individual] ) ;
      avgfitness += fit ;
      if (fit > maxfitness)
         maxfitness = fit ;
      if (fit < minfitness)
         minfitness = fit ;
      }

   avgfitness /= (double) popsize ;

/*
   The second step is to apply a linear transform to these fitnesses to prevent
   extraordinarily fit individuals from dominating early on, and at the same
   time still favor the most fit later in the run when a large number of
   individuals are very fit.
   
   This transform is:  f' = tmult * f + tconst.
   
   The coefficients are chosen so that the transformed maximum fitness is
   favor_best times the transformed average, while the average after transform
   is equal to that before.  A typical value for favor_best is 2-3.
   One problem is that late in the run, when the average is close to the max,
   very small fitnesses may map negative.  In this case, map the smallest
   to zero and do the best we can for the max.

   Note that a common alternative is to use the mapping just described, and
   truncate transformed fitnesses at zero.  However, the method shown here
   is usually superior, as it preserves genetic diversity.
*/

   ftemp = maxfitness - avgfitness ;
   if (ftemp > 1.e-20) {  // Insurance: average may equal max!
      tmult = (favor_best - 1.0) * avgfitness / ftemp ;
      tconst = avgfitness * (maxfitness - favor_best * avgfitness) / ftemp ;
      }
   else {
      tmult = 1.0 ;
      tconst = 0.0 ;
      }

/*
   The 'ideal' scaling factor was just computed.  Use it to map the minimum
   fitness.  If it comes out negative, compute an alternative scaling factor
   which will map the minimum to zero and keep the average unchanged.
*/

   if (tmult * minfitness + tconst < 0.0) { // Do not allow negative fitness
      ftemp = avgfitness - minfitness ;
      if (ftemp > 1.e-20) {
         tmult = avgfitness / ftemp ;
         tconst = -minfitness * avgfitness / ftemp ;
         }
      else {
         tmult = 1.0 ;
         tconst = 0.0 ;
         }
      }

/*
   The scaling factors have been computed.  Do the scaling now.
   The truncation at zero is theoretically unnecessary, as we avoided
   negatives when we computed the scaling factor above.  However, floating
   point problems can sometimes still produce a 'negative zero'.  In deference
   to possible user modifications which DEMAND nonnegative fitnesses, it is
   good programming practice to enforce this.
*/

   avgfitness = 0.0 ;
   for (individual=0 ; individual<popsize ; individual++) {
      fit = tmult * fitness[individual] + tconst ;
      if (fit < 0.0)
         fit = 0.0 ;
      fitness[individual] = fit ;
      avgfitness += fit ;
      }

   avgfitness /= (double) popsize ;

/*
   The final step is to normalize the fitnesses by dividing each by the
   average fitness.  The effect is that then each fitness can be interpreted
   as the expected number of times it would be chosen from the population
   pool if its probability of selection were proportional to its fitness.
*/

   for (individual=0 ; individual<popsize ; individual++)
      fitness[individual] /= avgfitness ;
}

/*
--------------------------------------------------------------------------------

   fitness_to_choices - Convert the array of fitnesses (which contain
               expected frequency of selection) into the array of parent
               choices.  This will allow random selection of parents
               without replacement later, while still insuring that
               we select (to within one) the expected number of each.


--------------------------------------------------------------------------------
*/

static void fitness_to_choices (
   int popsize ,      // Length of fitness, choices vectors
   double *fitness ,  // Input array of expected selection frequencies
   int *choices       // Output array of parents
   )
{
   int individual, expected, k ;

/*
   We build the choices array in two steps.  This, the first step, assigns
   parents according to the integer part of their expected frequencies.
*/

   k = 0 ;  // Will index choices array
   for (individual=0 ; individual<popsize ; individual++) {
      expected = (int) fitness[individual] ; // Assign this many now
      fitness[individual] -= expected ;      // Save fractional remainder
      while (expected--)                     // Forcibly use the int expected
         choices[k++] = individual ;         // quantity of this individual
      }

/*
   The second step is to take care of the remaining fractional expected
   frequencies.  Pass through the population, randomly selecting members
   with probability equal to their remaining fractional expectation.
   It is tempting to think that the algorithm below could loop excessively
   due to a very small fitness.  But recall that the sum of the fitnesses will
   be AT LEAST as large as the number remaining to be selected, and generally
   much more.  Thus, the ones with very small fitness (likely to cause trouble)
   will never become the only remaining possibilities.
*/

   while (k < popsize) {  // Select until choices is full
      individual = unifrand() * popsize ;// Randomly select individual
      if (fitness[individual] > 0.0) {   // Try members still having expectation
         if (fitness[individual] >= unifrand()) { // Selects with this prob
            choices[k++] = individual ;   // Bingo!  Select this individual
            fitness[individual] -= 1.0 ;  // and make it ineligable for future
            }
         }
      }
}

/*
--------------------------------------------------------------------------------

   pick_parents - Randomly select two parents from 'choices' candidate array

--------------------------------------------------------------------------------
*/

static void pick_parents (
   int *nchoices ,  // Number of choices (returned decremented by two)
   int *choices ,   // Array (nchoices long) of candidates for parent
   int *parent1 ,   // One parent returned here
   int *parent2     // and the other here
   )
{
   int k ;

   k = unifrand() * *nchoices ;        // Select position in choices array
   *parent1 = choices[k] ;             // Then return that parent
   choices[k] = choices[--*nchoices] ; // without replacement

   k = unifrand() * *nchoices ;
   *parent2 = choices[k] ;
   choices[k] = choices[--*nchoices] ;
}

/*
--------------------------------------------------------------------------------

   reproduce - Create a child from half of each of two parents.
      If first_child is true, randomly generate the crossover point.
      Else use the supplied crossover point and take other halves.
      If the chromosome size is at least 16, use two point crossover
      half of the time.  Signal this by returning a negative crosspt.
      We implement splitting the chromosome at a bit level (i.e. within
      a byte) with the 'split' parameter.  This is in the range 0-7,
      and is the bit number where the split takes place.  Like crosspt
      it is output for the first child and input for the second.

--------------------------------------------------------------------------------
*/

static void reproduce (
   char *p1 ,        // Pointer to one parent
   char *p2 ,        // and the other
   int first_child , // Is this the first of their 2 children?
   int chromsize ,   // Number of genes in chromosome
   char *child ,     // Output of a child
   int *crosspt ,    // If first_child, output of xover pt, else input it.
   int *split        // In/out of within byte splitting point
   )

{
   int i, n1, n2, n3, n4 ;
   char left, right, *pa, *pb ;

   if (first_child) {

      *split = longrand() % 8 ; // We will split boundary bytes here
      *crosspt = 1 + unifrand() * chromsize ;  // Randomly select cross pt

      if ((chromsize >= 16)  &&  (unifrand() < 0.33333)) // Two point?
         *crosspt = -*crosspt ; // flag this for second child

      pa = p1 ;
      pb = p2 ;
      } // If first child

   else {                       // Second child
      pa = p2 ;                 // so parents reverse roles
      pb = p1 ;
      } // If second child

/*
   Prepare for reproduction
*/

   if (*split) {              // Create left and right splitting masks
      right = 1 ;
      i = *split ;
      while (--i)
         right = (right << 1) | 1 ;
      left = 255 ^ right ;
      }

   if (*crosspt > 0) {        // Use one point crossover
      n1 = chromsize / 2 ;    // This many genes in first half of child
      n2 = chromsize - n1 ;   // and this many in second half
      n3 = n4 = 0 ;           // We are using one point crossover
      i = *crosspt - 1 ;      // We will start building child here
      }
   else {                             // Use two point crossover
      n1 = n2 = n3 = chromsize / 4 ;  // This many in first three quarters
      n4 = chromsize - n1 - n2 - n3 ; // And the last quarter gets the rest
      i = -*crosspt - 1 ;             // 2 point method was flagged by neg
      }

/*
   Do reproduction here
*/

   if (*split) {
      i = (i+1) % chromsize ;
      child[i] = (left & pa[i])  |  (right & pb[i]) ;
      --n1 ;
      }

   while (n1--) {
      i = (i+1) % chromsize ;
      child[i] = pb[i] ;
      }

   if (*split) {
      i = (i+1) % chromsize ;
      child[i] = (left & pb[i])  |  (right & pa[i]) ;
      --n2 ;
      }

   while (n2--) {
      i = (i+1) % chromsize ;
      child[i] = pa[i] ;
      }

   if (n4) {               // Two point crossover?

      if (*split) {
         i = (i+1) % chromsize ;
         child[i] = (left & pa[i])  |  (right & pb[i]) ;
         --n3 ;
         }

      while (n3--) {
         i = (i+1) % chromsize ;
         child[i] = pb[i] ;
         }

      if (*split) {
         i = (i+1) % chromsize ;
         child[i] = (left & pb[i])  |  (right & pa[i]) ;
         --n4 ;
         }

      while (n4--) {
         i = (i+1) % chromsize ;
         child[i] = pa[i] ;
         }

      } // If two point crossover
}

/*
--------------------------------------------------------------------------------

   mutate - apply the mutation operator to a single child

--------------------------------------------------------------------------------
*/
static void mutate (
   char *child ,   // Input/Output of the child
   int chromsize , // Number of variables in objective function
   double pmutate  // Probability of mutation
   )
{
   while (chromsize--) {
      if (unifrand() < pmutate)                          // Mutate this gene?
         child[chromsize] ^= (char) 1 << (longrand() % 8) ;  // Flip random bit
      }

}

/*
--------------------------------------------------------------------------------

   rand_ind - Randomly generate an individual's chromosome

--------------------------------------------------------------------------------
*/
static void rand_ind ( char *popptr , int chromsize )
{
   while (chromsize--)
      *popptr++ = 255 & longrand() ;
}

/*
--------------------------------------------------------------------------------

   decode - Decode the genes in this chromosome to the network weights.
            In genetic parlance, convert the genotype to the phenotype.

--------------------------------------------------------------------------------
*/

static void decode ( char *popptr , int nin , int nh1 , int nh2 ,
                     double *w1 , double *w2 )
{
   int n ;
   unsigned char gray, bit, parity, sum ;
   double *wptr ;

   if (! gray_initialized) {     // If the translation table has not yet
      gray_initialized = 1 ;     // been built, do it now (but just once!)
      for (n=0 ; n<256 ; n++) {  // Each gene is a one byte gray code
         gray = (unsigned char) n ;
         sum =  0 ;
         bit = 128 ;
         parity = 0 ;
         while (bit) {
            if (bit & gray)
               parity = ! parity ;
            if (parity)
               sum |= bit ;
            bit = bit >> 1 ;
            }
         gray_code_table[n] = sum ;
         }
      }

   n = nh1 * (nin+1) ; // Do hid layer 1 first.  It has this many weights.
   wptr = w1 ;         // Point to its weights

   while (n--) {
      gray = (unsigned char) *popptr++ ;  // The gene is a one byte gray code
      *wptr++ = ((double) gray_code_table[gray] - 127.5) * .0392 ; // -5 to 5
      }

   if (nh2) {              // Do second hidden layer if any
      n = nh2 * (nh1+1) ;
      wptr = w2 ;
      while (n--) {
         gray = (unsigned char) *popptr++ ;
         *wptr++ = ((double) gray_code_table[gray] - 127.5) * .0392 ;
         }
      }
}
