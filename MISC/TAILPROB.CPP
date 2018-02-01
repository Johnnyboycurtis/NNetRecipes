/******************************************************************************/
/*                                                                            */
/*  TAILPROB.CPP - Estimate tail probabilities froma sample                   */
/*                                                                            */
/*  Also includes leftconf and rightconf for interpolating to find confidence */
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

void shellsort ( int n , double *x ) ;

/*
--------------------------------------------------------------------------------

   tailprob - Given an array of sample data, this sorts it into ascending order
              and removes repeats.  It returns lprob and rprob such that
              lprob[i] is the fraction of samples less than or equal to x[i]
              (after x modified by this routine!) and rprob[i] is the number
              greater than or equal to x[i].
              It also returns ks_err, which, for given ks_alpha, is an alpha
              level confidence band around the lprob and rprob values.
              For example, if ks_alpha is input 0.05 then there is only a
              five percent chance that the true lprob and rprob differ from
              the estimate here by more than that amount.  This estimate is
              not exact; it is always conservative (ks_err returned too big).
              For n at least 35 it is very close; for n at least 100 it is
              essentially perfect.

--------------------------------------------------------------------------------
*/
void tailprob (
   int n ,             // Length of data array as input
   double *x ,         // Input of data, output sorted with repeats removed
   int *new_n ,        // Output length of x after repeats removed
   double *lprob ,     // Output new_n long of left tail probabilities
   double *rprob ,     // Ditto for right tail
   double ks_alpha ,   // Input (often .05) of type I error level for KS test
   double *ks_err      // Output of corresponding error in lprob, rprob
   )
{
   int i, j ;

   shellsort ( n , x ) ;                       // Sort ascending

   i = j = 0 ;                                 // i indexes source, j dest
   while (i < n) {                             // Pass through entire array
      x[j] = x[i] ;                            // Copy first occurence of value
      rprob[j] = (double) (n-i) / (double) n ; // Number greater or equal
      while ((++i < n)  &&  (x[i] == x[j])) ;  // Bypass repeats
      lprob[j++] = (double) i / (double) n ;   // Number less than or equal
      }

   *new_n = j ;

// The following is good if n>35, excellent if n>100, and always conservative

   *ks_err = sqrt ( -0.5 * log ( 0.5 * ks_alpha ) / (double) n ) ;
}

/*
--------------------------------------------------------------------------------

   leftconf - Given the results of tailprob and a data value, this returns a
              confidence value for rejecting the hypothesis that the
              distribution which gave rise to this data value is that from
              which the tailprob sample was obtained, in favor of an alternative
              that favors smaller data values.

--------------------------------------------------------------------------------
*/

double leftconf (
   int n ,             // Length of x and lprob arrays 
   double *x ,         // Input of data, sorted with repeats removed
   double *lprob ,     // Input of left tail probabilities from tailprob
   double observed     // Observed data value to be tested
   )
{
   int lo, mid, hi ;
   double p ;

   if (observed > x[n-1])
      return 0. ;

   if (observed <= x[0])
      return 100. * (1.0 - lprob[0]) ;

   lo = 0 ;    // At all times we will keep x[lo] strictly less than observed
   hi = n-1 ;  // and x[hi] greater or equal to observed

   for (;;) {                   // This loop cuts interval in half each time
      mid = (lo + hi) / 2 ;     // Center of interval
      if (mid == lo)            // Happens when lo and hi adjacent
         break ;                // So then we are done
      if (x[mid] < observed)    // Replace appropriate interval end with mid
         lo = mid ;
      else
         hi = mid ;
      }

   p = (observed - x[hi-1]) / (x[hi] - x[hi-1]) * (lprob[hi] - lprob[hi-1]) ;
   p += lprob[hi-1] ;
   return 100. * (1.0 - p) ;
}


/*
--------------------------------------------------------------------------------

   rightconf- Given the results of tailprob and a data value, this returns a
              confidence value for rejecting the hypothesis that the
              distribution which gave rise to this data value is that from
              which the tailprob sample was obtained, in favor of an alternative
              that favors larger data values.

--------------------------------------------------------------------------------
*/

double rightconf (
   int n ,             // Length of x and rprob arrays 
   double *x ,         // Input of data, sorted with repeats removed
   double *rprob ,     // Input of right tail probabilities from tailprob
   double observed     // Observed data value to be tested
   )
{
   int lo, mid, hi ;
   double p ;

   if (observed > x[n-1])
      return 100. * (1.0 - rprob[n-1]) ;

   if (observed <= x[0])
      return 0. ;

   lo = 0 ;    // At all times we will keep x[lo] strictly less than observed
   hi = n-1 ;  // and x[hi] greater or equal to observed

   for (;;) {                   // This loop cuts interval in half each time
      mid = (lo + hi) / 2 ;     // Center of interval
      if (mid == lo)            // Happens when lo and hi adjacent
         break ;                // So then we are done
      if (x[mid] < observed)    // Replace appropriate interval end with mid
         lo = mid ;
      else
         hi = mid ;
      }

   p = (observed - x[hi-1]) / (x[hi] - x[hi-1]) * (rprob[hi] - rprob[hi-1]) ;
   p += rprob[hi-1] ;
   return 100. * (1.0 - p) ;
}
