/******************************************************************************/
/*                                                                            */
/*  FUZZY - Fuzzy set stuff                                                   */
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

class FuzzySet {

public:
   FuzzySet () ;
   FuzzySet ( int npts , double *xpts , double *ypts ) ;
   FuzzySet ( double x0 , double y0 , double x1 , double y1 ,
              double x2 , double y2 , double x3 , double y3 ) ;
   FuzzySet ( const FuzzySet& s ) ;
   void operator= ( const FuzzySet& s ) ;
   ~FuzzySet () ;
   double centroid () ;
   void cj ( const FuzzySet& s , int conj ) ;
   int invalid() ;
   double membership ( double pt ) ;
   void negate () ;
   void print () ;
   void scale ( double fac ) ;

private:
   int n ;         // Number of points in membership function def
   double *x ;     // Domain points
   double *y ;     // Corresponding function values
} ;


/*
--------------------------------------------------------------------------------

   Constructors, destructor

--------------------------------------------------------------------------------
*/

FuzzySet::FuzzySet ()
{
   n = 0 ;
}

FuzzySet::FuzzySet ( int npts , double *xpts , double *ypts )
{
   if (npts < 1) {
      n = 0 ;    // Flag as invalid
      return ;
      }

   x = (double *) malloc ( 2 * npts * sizeof(double) ) ;
   if (x == NULL) {
      n = 0 ;    // Flag as invalid
      return ;
      }

   y = x + npts ;  // This saves a call to malloc

   memcpy ( x , xpts , npts * sizeof(double) ) ;
   memcpy ( y , ypts , npts * sizeof(double) ) ;
   n = npts ;
}


FuzzySet::FuzzySet ( double x0 , double y0 , double x1 , double y1 ,
                     double x2 , double y2 , double x3 , double y3 )
{
   double xx[4], yy[4] ;

   xx[0] = x0 ;
   yy[0] = y0 ;
   n = 1 ;

   if ((x1 != x0)  ||  (y1 != y0)) {
      xx[n] = x1 ;
      yy[n] = y1 ;
      ++n ;
      }

   if ((x2 != x1)  ||  (y2 != y1)) {
      xx[n] = x2 ;
      yy[n] = y2 ;
      ++n ;
      }

   if ((x3 != x2)  ||  (y3 != y2)) {
      xx[n] = x3 ;
      yy[n] = y3 ;
      ++n ;
      }

   x = (double *) malloc ( 2 * n * sizeof(double) ) ;
   if (x == NULL) {
      n = 0 ;    // Flag as invalid
      return ;
      }

   y = x + n ;  // This saves a call to malloc

   memcpy ( x , xx , n * sizeof(double) ) ;
   memcpy ( y , yy , n * sizeof(double) ) ;
}


FuzzySet::~FuzzySet ()
{
   if (n)
      free ( x ) ;  // Recall malloc'd only once
}

/*
--------------------------------------------------------------------------------

   Copy and Assignment Constructors

--------------------------------------------------------------------------------
*/

FuzzySet::FuzzySet ( const FuzzySet& s )
{
   n = s.n ;
   if (!n)
      return ;

   x = (double *) malloc ( 2 * n * sizeof(double) ) ;
   if (x == NULL) {
      n = 0 ;
      return ;
      }
   y = x + n ;
   memcpy ( x , s.x , n * sizeof(double) ) ;
   memcpy ( y , s.y , n * sizeof(double) ) ;
}

void FuzzySet::operator= ( const FuzzySet& s )
{
   if (this == &s)  //  Handles s=s correctly
      return ;

   if (n)
      free ( x ) ;

   n = s.n ;
   if (!n)
      return ;

   x = (double *) malloc ( 2 * n * sizeof(double) ) ;
   if (x == NULL) {
      n = 0 ;
      return ;
      }
   y = x + n ;
   memcpy ( x , s.x , n * sizeof(double) ) ;
   memcpy ( y , s.y , n * sizeof(double) ) ;
}

/*
--------------------------------------------------------------------------------

   Membership evaluator, invalid indicator, printer

--------------------------------------------------------------------------------
*/

double FuzzySet::membership ( double pt )
{
   int lo, mid, hi ;
   double yy ;

   if (! n)
      return 0.0 ;

   if (pt <= x[0])
      return y[0] ;

   if (pt >= x[n-1])
      return y[n-1] ;

   lo = 0 ;    // We will keep x[lo] strictly less than pt
   hi = n-1 ;  // and x[hi] greater or equal to pt

   for (;;) {                   // Cuts interval in half each time
      mid = (lo + hi) / 2 ;     // Center of interval
      if (mid == lo)            // Happens when lo and hi adjacent
         break ;                // So then we are done
      if (x[mid] < pt)    // Replace appropriate interval end with mid
         lo = mid ;
      else
         hi = mid ;
      }

   yy = (pt - x[hi-1]) / (x[hi] - x[hi-1]) * (y[hi] - y[hi-1]) ;
   return yy + y[hi-1] ;
}

FuzzySet::invalid()
{
   return !n ;
}

void FuzzySet::print()
{
   int i ;
   for (i=0 ; i<n ; i++)
      printf ( "\n%8.3lf %8.3lf", x[i], y[i] ) ;
}

/*
--------------------------------------------------------------------------------

   Negation, Inference by scaling

--------------------------------------------------------------------------------
*/

void FuzzySet::negate ()
{
   int i ;

   i = n ;
   while (i--)
      y[i] = 1.0 - y[i] ;
}

void FuzzySet::scale ( double fac )
{
   int i ;

   i = n ;
   while (i--)
      y[i] = fac * y[i] ;
}

/*
--------------------------------------------------------------------------------

   Conjunction, disjunction

--------------------------------------------------------------------------------
*/

static void cj_out ( double x , double y , double *xn , double *yn , int *nn ) ;
static void intsec ( double x1 , double y1 , double x2 , double y2 ,
                     double y3 , double y4 , double *xint , double *yint ) ;

void FuzzySet::cj (
   const FuzzySet& s, // Modifying set
   int conj           // conjunction if 1, disjunction if 0
   )
{
   int i, nn, vertex0, vertex1, use_func_0 ;
   double *xn, *yn, xl, xr, y0l, y0r, y1l, y1r, xint, yint ;
   double next_x0, next_y0, next_x1, next_y1, rightmost_x, frac ;

   if (!n)         // Cannot operate on an invalid set
      return ;

   if (!s.n) {     // Propagate invalidity
      free ( x ) ;
      n = 0 ;
      return ;
      }

   xn = (double *) malloc ( 4 * (n + s.n) * sizeof(double) ) ;
   if (xn == NULL) {
      n = 0 ;
      free ( x ) ;
      return ;
      }
   yn = xn + 2 * (n + s.n) ;

   conj = (conj != 0) ;  // Sloppy user may not restrict to 0/1

/*
   In the following, the suffix '0' refers to the current set,
   being modified.  '1' means the external, modifying set.
   The domain endpoints are always xl and xr.  The corresponding
   function values are y0l, y0r, y1l, y1r.
*/


/*
   Initialize the leftmost point.  If one of the domains starts
   lower than the other, imply a continued function value for the
   other so they can start together.
*/

   if (x[0] < s.x[0]) { // Start result at leftmost x
      xl = x[0] ;       // Which is this one
      vertex0 = 1 ;     // Subscript of next point
      vertex1 = 0 ;     // Other was an implicit point
      }
   else if (x[0] > s.x[0]) { // Same thing, other set
      xl = s.x[0] ;
      vertex0 = 0 ;
      vertex1 = 1 ;
      }
   else {
      xl = x[0] ;
      vertex0 = vertex1 = 1 ;
      }

   y0l = y[0] ;
   y1l = s.y[0] ;

   nn = 0 ;
   cj_out ( xl , (conj ^ (y0l>y1l)) ? y0l : y1l , xn , yn , &nn ) ;

/*
   Also initialize the rightmost point so we end together.
*/

   if (x[n-1] >= s.x[s.n-1])
      rightmost_x = x[n-1] ;
   else 
      rightmost_x = s.x[s.n-1] ;

/*
   At any time from now on, next_?? will be the next point in
   each function.  We will advance by choosing the closest.
*/

   if (vertex0 < n) {
      next_x0 = x[vertex0] ;
      next_y0 = y[vertex0] ;
      }
   else {
      next_x0 = rightmost_x ;
      next_y0 = y[vertex0-1] ;
      }

   if (vertex1 < s.n) {
      next_x1 = s.x[vertex1] ;
      next_y1 = s.y[vertex1] ;
      }
   else {
      next_x1 = rightmost_x ;
      next_y1 = s.y[vertex1-1] ;
      }

/*
   Main loop is here
*/

   while ((vertex0 < n) || (vertex1 < s.n)) {

      if (next_x0 < next_x1)          // Function 0 has closest next x
         use_func_0 = 1 ;
      else if (next_x1 < next_x0)     // Function 1 has closest next x
         use_func_0 = 0 ;
      else
         use_func_0 = (vertex0 < n) ; // Break tie by who can advance

      if (use_func_0) {

         xr = next_x0 ;
         y0r = next_y0 ;
         if (next_x1 == xl)  // Interpolate to corresponding y
            frac = 0.0 ;
         else 
            frac = (xr - xl) / (next_x1 - xl) ;
         y1r = y1l + frac * (next_y1 - y1l) ;

         if (++vertex0 < n) {         // Advance in this function
            next_x0 = x[vertex0] ;
            next_y0 = y[vertex0] ;
            }
         else 
            next_x0 = rightmost_x ;
         }

      else {                         // Function 1 has closest next x

         xr = next_x1 ;
         y1r = next_y1 ;
         if (next_x0 == xl)  // Interpolate to corresponding y
            frac = 0.0 ;
         else 
            frac = (xr - xl) / (next_x0 - xl) ;
         y0r = y0l + frac * (next_y0 - y0l) ;

         if (++vertex1 < s.n) {        // Advance in this function
            next_x1 = s.x[vertex1] ;
            next_y1 = s.y[vertex1] ;
            }
         else 
            next_x1 = rightmost_x ;
         }

/*
   If the functions cross, generate new point at intersection
*/

      if ((xr > xl)  &&  ((y0l-y1l) * (y0r-y1r) < 0.0)) {
         intsec ( xl , y0l , xr , y0r , y1l , y1r , &xint , &yint ) ;
         cj_out ( xint , yint , xn , yn , &nn ) ;
         }

      cj_out ( xr , (conj ^ (y0r>y1r)) ? y0r : y1r , xn , yn , &nn ) ;

      xl = xr ;
      y0l = y0r ;
      y1l = y1r ;

      } // Main loop

/*
   The logic is done.  Save the results.
*/

   free ( x ) ;
   for (i=0 ; i<nn ; i++)
      xn[nn+i] = yn[i] ;
   x = (double *) realloc ( xn , 2 * nn * sizeof(double) ) ;
   y = x + nn ;
   n = nn ;
}

/*
--------------------------------------------------------------------------------

   cj_out - Local routine to output a point
   intsec - Local routine to intersect two line segments

--------------------------------------------------------------------------------
*/

static void cj_out (  // save a resultant point
   double x ,         // x value of new point
   double y ,         // and y
   double *xn ,       // Array of x's
   double *yn ,       // and y's
   int *nn            // Number in above
   )
{
   double d ;

/*
   Do not store identical points
*/

   if ((*nn > 0)  &&
      (fabs(x-xn[*nn-1]) < 1.e-20)  &&  (fabs(y-yn[*nn-1]) < 1.e-20))
      return ;

/*
   Do not store colinear points
*/

   if (*nn > 1) {
      d = xn[*nn-1] * (y - yn[*nn-2]) +
          yn[*nn-1] * (xn[*nn-2] - x) +
          x * yn[*nn-2] - y * xn[*nn-2] ;
      if (fabs(d) < 1.e-10)
         --*nn ; // Colinear, so replace midpoint
      }

   xn[*nn] = x ;
   yn[*nn] = y ;
   ++*nn ;
}


static void intsec ( double x1 , double y1 , double x2 , double y2 ,
                     double y3 , double y4 , double *xint , double *yint )
{
   double den ;

   den = y1 - y2 - y3 + y4 ; // Guaranteed >0 if they cross
   *xint = (x2 * y1 - x1 * y2 - x2 * y3 + x1 * y4) / den ;
   *yint = (y1 * y4 - y2 * y3) / den ;
}

/*
--------------------------------------------------------------------------------

   centroid - Compute the centroid of a fuzzy set

--------------------------------------------------------------------------------
*/

double FuzzySet::centroid ()
{
   int i ;
   double t, numer, denom ;

   if (n < 2)
      return 0.0 ;

   numer = denom = 0.0 ;

   for (i=1 ; i<n ; i++) {
      t = y[i-1] * (2.0 * x[i-1] + x[i]) + y[i] * (x[i-1] + 2.0 * x[i]) ;
      numer += t * (x[i] - x[i-1]) ;
      denom += (x[i] - x[i-1]) * (y[i-1] + y[i]) ;
      }
   return numer / (3.0 * denom) ;
}
