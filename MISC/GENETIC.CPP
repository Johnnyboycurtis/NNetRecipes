/******************************************************************************/
/*                                                                            */
/*  GENETIC - Minimize a function using genetic methods                       */
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

#define RANDMAX 32767

static void shake ( int nvars , double *center , double *x , double temp ) ;
static void fval_to_fitness ( int popsize , double favor_best ,
                            double fitfac , double *fvals , double *fitness ) ;
static void fitness_to_choices ( int popsize , double *fitness , int *choices );
static void pick_parents ( int *nchoices , int *choices ,
                           int *parent1 , int *parent2 ) ;
static void reproduce ( double *p1 , double *p2 , int first_child ,
                        int nvars , double *child , int *crosspt ) ;
static void mutate ( double *child , int nvars , double stddev ,
                     double pmutate ) ;

double genetic (
      double (*func)( int n , double *x ) ,  // Function to be minimized
      int ngens ,        // Number of complete generations
      int popsize ,      // Number of individuals in population
      int climb ,        // Do we hill climb via elitism?
      double stddev ,    // Standard deviation for initial population
      double pcross ,    // Probability of crossover (.6-.9 typical)
      double pmutate ,   // Probability of mutation (.0001 to .001 typical)
      double favor_best ,// Factor for favoring best over average (2-3 is good)
      double fitfac ,    // Factor for converting fval to raw fitness (-20 good)
      double fquit ,     // Quit if function reduced to this amount
      int nvars ,        // Number of variables in x
      double *x ,        // Set to starting estimate when called; returns best
      double *best,      // Work vector nvars long
      int *choices ,     // Work vector popsize long
      double *fvals ,    // Work vector popsize long
      double *fitness ,  // Work vector popsize long
      double *pool1,     // Work vector nvars * popsize long
      double *pool2 )    // Work vector nvars * popsize long
{
   int istart, individual, best_individual, generation, n_cross ;
   int first_child, parent1, parent2, improved, crosspt, nchoices ;
   double fval, bestfval, *oldpop, *newpop, *popptr, *temppop ;

/*
   Generate initial population pool.

   We also preserve the best point across all generations, as this is what
   we will ultimately return to the user.  Its objective function value is
   bestfval.
*/

   bestfval = 1.e30 ;     // For saving best (across all individuals)
   best_individual = 0 ;  // Safety only

   for (individual=0 ; individual<popsize ; individual++) {

      popptr = pool1 + individual * nvars ;   // Build population in pool1
      shake ( nvars , x , popptr , stddev ) ; // Randomly perturb about init x
      fval = (*func) ( nvars , popptr ) ;     // Evaluate function there
      fvals[individual] = fval ;              // and keep function value of each

      if (fval < bestfval) {                  // Keep track of best
         bestfval = fval ;                    // as it will be returned to user
         best_individual = individual ;       // This is its index in pool1
         }

      if (fval <= fquit)
         break ;
      }

/*
   The initial population has been built in pool1.
   Copy its best member to 'best' in case it never gets beat (unlikely
   but possible!).
*/

   popptr = pool1 + best_individual * nvars ;           // Point to best
   memcpy ( best , popptr , nvars * sizeof(double) ) ;  // and save it

/*
   This is the main generation loop.  There are two areas for population pool
   storage: pool1 and pool2.  At any given time, oldpop will be set to one of
   them, and newpop to the other.  This avoids a lot of copying.
*/

   oldpop = pool1 ;   // This is the initial population
   newpop = pool2 ;   // The next generation is created here

   for (generation=0 ; generation<ngens ; generation++) {

      if (fval <= fquit)   // We may have satisfied this in initial population
         break ;           // So we test at start of generation loop

      fval_to_fitness ( popsize , favor_best , fitfac , fvals , fitness ) ;

      fitness_to_choices ( popsize , fitness , choices ) ;

      nchoices = popsize ;         // Will count down as choices array emptied
      n_cross = pcross * popsize ; // Number crossing over
      first_child = 1 ;            // Generating first of parent's 2 children?
      improved = 0 ;               // Flags if we beat best

      if (climb) {                  // If we are to hill climb
         memcpy ( newpop , best , nvars * sizeof(double) ) ; // start with best
         fvals[0] = bestfval ;      // Record its error
         istart = 1 ;               // and start children past it
         }
      else
         istart = 0 ;

/*
   Generate the children
*/

      for (individual=istart ; individual<popsize ; individual++) {

         popptr = newpop + individual * nvars ; // Will put this child here

         if (first_child)  // If this is the first of 2 children, pick parents
            pick_parents ( &nchoices , choices , &parent1 , &parent2 ) ;

         if (n_cross-- > 0)    // Do crossovers first
            reproduce ( oldpop + parent1 * nvars , oldpop + parent2 * nvars ,
                        first_child , nvars , popptr , &crosspt ) ;
         else if (first_child) // No more crossovers, so just copy parent
            memcpy( popptr , oldpop + parent1 * nvars , nvars * sizeof(double));
         else
            memcpy( popptr , oldpop + parent2 * nvars , nvars * sizeof(double));

         if (pmutate > 0.0)
            mutate ( popptr , nvars , stddev , pmutate ) ;

         fval = (*func) ( nvars , popptr ) ; // Evaluate function for this child
         fvals[individual] = fval ;          // and keep function value of each

         if (fval < bestfval) {              // Keep track of best
            bestfval = fval ;                // It will be returned to user
            best_individual = individual ;   // This is its index in newpop
            improved = 1 ;                   // Flag so we copy it later
            }

         if (fval <= fquit)
            break ;

         first_child = ! first_child ;
         } // For all genes in population

/*
   We finished generating all children.  If we improved (one of these
   children beat the best so far) then copy that child to the best.
   Swap oldpop and newpop for the next generation.
*/

      if (improved) {
         popptr = newpop + best_individual * nvars ;          // Point to best
         memcpy ( best , popptr , nvars * sizeof(double) ) ;  // and save it
         }

      temppop = oldpop ;   // Switch old and new pops for next generation
      oldpop = newpop ;
      newpop = temppop ;
      }

/*
   We are all done.  Copy the best to x, as that is how we return it.
*/

   memcpy ( x , best , nvars * sizeof(double) ) ;  // Return best
   return bestfval ;
}

/*
--------------------------------------------------------------------------------

   SHAKE - Randomly perturb the point

--------------------------------------------------------------------------------
*/

static void shake ( int nvars , double *center , double *x , double temp )
{
   double r ;

// Recall that the variance of a uniform deviate on 0-1 is 1/12.
// Adding four such random variables multiplies the variance by 4,
// while dividing by 2 divides the variance by 4.

   temp *= 3.464101615 / (2. * (double) RANDMAX ) ;  // SQRT ( 12 ) = 3.464...

   while (nvars--) {
      r = (double) rand() + (double) rand() -
          (double) rand() - (double) rand() ;
      *x++ = *center++ + temp * r ;
      }
}

/*
--------------------------------------------------------------------------------

   fval_to_fitness - Convert the objective function value of each individual
                     to a scaled fitness value.  The scaled fitness may be
                     considered an expected frequency of choice.

--------------------------------------------------------------------------------
*/

static void fval_to_fitness (
   int popsize ,       // Length of fvals, fitness vectors
   double favor_best , // Factor for favoring best over average (2 is good)
   double fitfac ,     // Factor for converting fval to raw fitness (-20 good)
   double *fvals ,     // Input popsize vector of values of objective function
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
      fitness[individual] = fit = exp ( fitfac * fvals[individual] ) ;
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

   fitness_to_choices - Convert the array of fitnesses (which contain expected
                        frequency of selection) into the array of parent
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
   double rn ;

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
      rn = (double) rand () / ((double) RANDMAX + 1.0) ; // 0-1 random
      individual = rn * popsize ;        // Randomly select an individual
      if (fitness[individual] > 0.0) {   // Try members still having expectation
         rn = (double) rand () / ((double) RANDMAX + 1.0) ; // 0-1 random
         if (fitness[individual] >= rn) { // Selects with this probability
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
   double rn ;

   rn = (double) rand () / ((double) RANDMAX + 1.0) ; // Random 0-1
   k = rn * *nchoices ;                // Select position in choices array
   *parent1 = choices[k] ;             // Then return that parent
   choices[k] = choices[--*nchoices] ; // without replacement

   rn = (double) rand () / ((double) RANDMAX + 1.0) ; // Ditto parent 2
   k = rn * *nchoices ;
   *parent2 = choices[k] ;
   choices[k] = choices[--*nchoices] ;
}

/*
--------------------------------------------------------------------------------

   reproduce - Create a child from half of each of two parents.
               If first_child is true, randomly generate the crossover point.
               Else use the supplied crossover point and take other halves.

--------------------------------------------------------------------------------
*/

static void reproduce (
   double *p1 ,        // Pointer to one parent
   double *p2 ,        // and the other
   int first_child ,   // Is this the first of their 2 children?
   int nvars ,         // Number of variables in objective function
   double *child ,     // Output of a child
   int *crosspt        // If first_child, output of xover pt, else input it.
   )
{
   int i, n1, n2 ;
   double rn, *pa, *pb ;

   if (first_child) {
      rn = (double) rand () / ((double) RANDMAX + 1.0) ; // Random 0-1
      *crosspt = rn * nvars ;  // Randomly select crossover point
      pa = p1 ;
      pb = p2 ;
      }
   else {
      pa = p2 ;
      pb = p1 ;
      }

   n1 = nvars / 2 ;    // This many genes in first half of child
   n2 = nvars - n1 ;   // and this many in second half
   i = *crosspt ;

   while (n1--) {
      i = (i+1) % nvars ;
      child[i] = pa[i] ;
      }

   while (n2--) {
      i = (i+1) % nvars ;
      child[i] = pb[i] ;
      }
}

/*
--------------------------------------------------------------------------------

   mutate - apply the mutation operator to a single child

--------------------------------------------------------------------------------
*/
static void mutate (
   double *child ,     // Input/Output of the child
   int nvars ,         // Number of variables in objective function
   double stddev ,     // Standard deviation of mutation
   double pmutate      // Probability of mutation
   )
{
   double rn ;

   while (nvars--) {
      rn = (double) rand () / ((double) RANDMAX + 1.0) ; // Random 0-1
      if (rn < pmutate) {                                // Mutate this gene?
         rn = (double) rand() + (double) rand() -
              (double) rand() - (double) rand() ;
         child[nvars] += rn * stddev * 3.464101615 / (2. * (double) RANDMAX ) ;
         }
      }

}
