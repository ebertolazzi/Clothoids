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
// file: Utils.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils.hh"

#include "Utils_fmt.hh"

namespace Utils
{

  using std::lower_bound;
  using std::string;

  /*\
   |                           _     ___       _                       _
   |   ___  ___  __ _ _ __ ___| |__ |_ _|_ __ | |_ ___ _ ____   ____ _| |
   |  / __|/ _ \/ _` | '__/ __| '_ \ | || '_ \| __/ _ \ '__\ \ / / _` | |
   |  \__ \  __/ (_| | | | (__| | | || || | | | ||  __/ |   \ V / (_| | |
   |  |___/\___|\__,_|_|  \___|_| |_|___|_| |_|\__\___|_|    \_/ \__,_|_|
  \*/
  //!
  //!  \brief Searches for the interval that contains a given value in a sorted
  //!  array.
  //!
  //!  This function searches for the interval within which the value `x` falls
  //!  based on a provided sorted array of `npts` points. It updates the
  //!  `last_interval` to reflect the correct interval in the array.
  //!
  //!  If the `closed` parameter is true, the function will adjust the value of
  //!  `x` to fit within the range defined by the first and last elements of the
  //!  array. If `can_extend` is false and `x` is outside this range, an
  //!  assertion will fail.
  //!
  //!  \tparam T_int The type used for the integer index (e.g., `int`).
  //!  \tparam T_real The type used for the real number (e.g., `double`).
  //!
  //!  \param npts The number of points in the array `X`. Must be greater than
  //!  or equal to 2.
  //!  \param X An array of sorted real numbers defining the points.
  //!  \param x A reference to the value being searched for in the array.
  //!  \param last_interval A reference to the last known interval index. It
  //!  will be updated
  //!                       to the correct interval index after the search.
  //!  \param closed A boolean indicating if the range is closed. If true, `x`
  //!  will be
  //!                adjusted to fit within the range defined by the first and
  //!                last elements.
  //!  \param can_extend A boolean indicating if `x` is allowed to be outside
  //!  the range
  //!                    defined by the first and last elements.
  //!
  //!  \note The function makes use of assertions to ensure that input values
  //!  are valid
  //!        and that the calculated interval is within the expected range.
  //!
  //!  \throw std::runtime_error If assertions fail due to invalid input
  //!  parameters or
  //!         out-of-range conditions.
  //!
  template <typename T_int, typename T_real>
  void
  search_interval( T_int npts, T_real const X[], T_real & x, T_int & last_interval, bool closed, bool can_extend )
  {
    // check points
    T_int n{ npts - 1 };
    UTILS_ASSERT( npts > 1 && last_interval >= 0 && last_interval < n,
                  "In search_interval( npts={}, X, x={}, last_interval={}, "
                  "closed={}, can_extend={})\n"
                  "npts must be >= 2 and last_interval must be in [0,npts-2]\n",
                  npts, x, last_interval, closed, can_extend );

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
      UTILS_ASSERT( can_extend || ( x >= xl && x <= xr ),
                    "In search_interval( npts={}, X, x={}, last_interval={}, "
                    "closed={}, can_extend={})\n"
                    "out of range: [{},{}]\n",
                    npts, x, last_interval, closed, can_extend, xl, xr );
    }

    // find the interval of the support of the B-spline
    T_real const * XL{ X + last_interval };
    if ( XL[1] < x )
    {  // x on the right
      if ( x >= X[n - 1] )
      {
        last_interval = n - 1;  // last interval
      }
      else if ( x < XL[2] )
      {  // x in (XL[1],XL[2])
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
    }
    // check computed interval
    UTILS_ASSERT( last_interval >= 0 && last_interval < n,
                  "In search_interval( npts={}, X, x={}, last_interval={}, "
                  "closed={}, can_extend={})\n"
                  "computed last_interval of range: [{},{}]\n",
                  npts, x, last_interval, closed, can_extend, xl, xr );
  }

  extern template void search_interval( int32_t     npts,
                                        float const X[],
                                        float &     x,
                                        int32_t &   last_interval,
                                        bool        closed,
                                        bool        can_extend );

  template void search_interval( int32_t      npts,
                                 double const X[],
                                 double &     x,
                                 int32_t &    last_interval,
                                 bool         closed,
                                 bool         can_extend );

  template void search_interval( int64_t     npts,
                                 float const X[],
                                 float &     x,
                                 int64_t &   last_interval,
                                 bool        closed,
                                 bool        can_extend );

  template void search_interval( int64_t      npts,
                                 double const X[],
                                 double &     x,
                                 int64_t &    last_interval,
                                 bool         closed,
                                 bool         can_extend );
}  // namespace Utils

#endif

//
// eof: Utils.cc
//
