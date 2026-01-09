/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_search_intervals.hh
//

#pragma once

#ifndef UTILS_SEARCH_INTERVALS_HH
#define UTILS_SEARCH_INTERVALS_HH

#include "Utils.hh"

namespace Utils
{

  // ============================================================================
  // INTERVAL SEARCH FUNCTIONS
  // ============================================================================

  /**
   * \brief Searches for the interval containing a value in a sorted array
   * \tparam T_int Integer type for indices
   * \tparam T_real Real type for array values
   * \param npts Number of points in the array (must be >= 2)
   * \param X Sorted array of real numbers
   * \param x Value to locate (may be adjusted if closed range)
   * \param last_interval Reference to last known interval index, updated with found interval
   * \param closed If true, treats array as closed periodic range
   * \param can_extend If false, x must be within array bounds
   */
  template <typename T_int, typename T_real> inline void search_interval(
    T_int        npts,
    T_real const X[],
    T_real &     x,
    T_int &      last_interval,
    bool         closed,
    bool         can_extend )
  {
    using std::lower_bound;

    // check points
    T_int n{ npts - 1 };
    UTILS_ASSERT(
      npts > 1 && last_interval >= 0 && last_interval < n,
      "In search_interval( npts={}, X, x={}, last_interval={}, "
      "closed={}, can_extend={})\n"
      "npts must be >= 2 and last_interval must be in [0,npts-2]\n",
      npts,
      x,
      last_interval,
      closed,
      can_extend );

    // check range
    T_real xl{ X[0] };
    T_real xr{ X[n] };
    if ( closed )
    {  // put x in range (change also its value)
      T_real L{ xr - xl };
      x -= xl;
      x = fmod( x, L );
      if ( x < 0 ) x += L;
      x += xl;
    }
    else
    {
      UTILS_ASSERT(
        can_extend || ( x >= xl && x <= xr ),
        "In search_interval( npts={}, X, x={}, last_interval={}, "
        "closed={}, can_extend={})\n"
        "out of range: [{},{}]\n",
        npts,
        x,
        last_interval,
        closed,
        can_extend,
        xl,
        xr );
    }

    // find the interval of the support of the B-spline
    T_real const * XL{ X + last_interval };

    // Use <= instead of < for the right boundary check to handle exact matches
    // but only if we're not at the last interval
    if ( XL[1] < x || ( x == XL[1] && last_interval < n - 1 && !Utils::is_zero( XL[1] - XL[2] ) ) )
    {  // x on the right (or exactly at the right boundary and not the last interval)
      if ( x >= X[n - 1] )
      {
        last_interval = n - 1;  // last interval
      }
      else if ( x < XL[2] )
      {  // x in (XL[1],XL[2]) or exactly at XL[1] when not degenerate
        ++last_interval;
      }
      else
      {  // x >= XL[2] search the right interval
        T_real const * XE{ X + n };
        last_interval += T_int( lower_bound( XL, XE, x ) - XL );
        T_real const * XX = X + last_interval;
        if ( x < XX[0] || Utils::is_zero( XX[0] - XX[1] ) ) --last_interval;
      }
    }
    else if ( x < XL[0] )
    {  // on the left
      if ( x <= X[1] )
      {                     // x in [X[0],X[1]]
        last_interval = 0;  // first interval
      }
      else if ( XL[-1] <= x )
      {  // x in [XL[-1],XL[0])
        --last_interval;
      }
      else
      {
        last_interval     = T_int( lower_bound( X + 1, XL, x ) - X );
        T_real const * XX = X + last_interval;
        if ( x < XX[0] || Utils::is_zero( XX[0] - XX[1] ) ) --last_interval;
      }
    }
    else
    {
      // x in the interval [ XL[0], XL[1] ] nothing to do
      // If x is exactly at XL[1] and we're at the last interval or it's a degenerate interval, stay here
    }
    // check computed interval
    UTILS_ASSERT(
      last_interval >= 0 && last_interval < n,
      "In search_interval( npts={}, X, x={}, last_interval={}, "
      "closed={}, can_extend={})\n"
      "computed last_interval of range: [{},{}]\n",
      npts,
      x,
      last_interval,
      closed,
      can_extend,
      xl,
      xr );
  }

  /**
   * \brief Legacy wrapper for search_interval function
   * \deprecated Use search_interval() instead
   */
  template <typename T_int, typename T_real> inline void searchInterval(
    T_int        npts,
    T_real const X[],
    T_real &     x,
    T_int &      last_interval,
    bool         closed,
    bool         can_extend )
  {
    search_interval( npts, X, x, last_interval, closed, can_extend );
  }

}  // namespace Utils

#endif  // UTILS_INTERVALS_HH

//
// eof: Utils_search_intervals.hh
//
