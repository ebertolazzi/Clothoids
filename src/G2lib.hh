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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: G2lib.hh
///

#ifndef G2LIB_HH
#define G2LIB_HH

#include <iostream>
#include <cmath>
#include <cfloat>
#include <sstream>
#include <stdexcept>
#include <limits>

#ifndef G2LIB_ASSERT
  #define G2LIB_ASSERT(COND,MSG)           \
    if ( !(COND) ) {                       \
      std::ostringstream ost;              \
      ost << "On line: " << __LINE__       \
          << " file: " << __FILE__         \
          << '\n' << MSG << '\n';          \
      throw std::runtime_error(ost.str()); \
    }
#endif

// select computer architecture
#if defined(__APPLE__) && defined(__MACH__)
  // osx architecture
  #define G2LIB_OS_OSX 1
  #if defined(__i386__)
    #define G2LIB_ARCH32 1
  #elif defined(__x86_64__)
    #define G2LIB_ARCH64 1
  #endif
#elif defined(__unix__)
  // linux architecture
  #define G2LIB_OS_LINUX 1
  #if defined(__i386__)
    #define G2LIB_ARCH32 1
  #elif defined(__x86_64__)
    #define G2LIB_ARCH64 1
  #endif
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define G2LIB_OS_WINDOWS 1
  #if defined(_M_X64) || defined(_M_AMD64)
    #define G2LIB_ARCH64 1
  #else
    #define G2LIB_ARCH32 1
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else
  #error "unsupported OS!"
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
  #ifndef G2LIB_DO_NOT_USE_CXX11
    #define G2LIB_USE_CXX11
  #endif
#else
  // not C++11 compiler
  #ifndef nullptr
    #define nullptr NULL
  #endif
#endif


//! Clothoid computations routine
namespace G2lib {

  typedef std::basic_ostream<char> ostream_type;

  typedef double real_type;
  typedef int    int_type;

  extern real_type const machepsi;
  extern real_type const machepsi10;
  extern real_type const machepsi100;
  extern real_type const machepsi1000;
  extern real_type const m_pi;        // pi
  extern real_type const m_pi_2;      // pi/2
  extern real_type const m_2pi;       // 2*pi
  extern real_type const m_1_pi;      // 1/pi
  extern real_type const m_1_sqrt_pi; // 1/sqrt(pi)

  static
  inline
  bool
  isZero( real_type x )
  { return FP_ZERO == std::fpclassify(x); }

  static
  inline
  bool
  isInfinite( real_type x )
  { return FP_INFINITE == std::fpclassify(x); }

  static
  inline
  bool
  isNaN( real_type x )
  { return FP_NAN == std::fpclassify(x); }

  static
  inline
  bool
  isRegular( real_type x )
  { return !( FP_INFINITE == std::fpclassify(x) ||
              FP_NAN      == std::fpclassify(x) ); }

  /*
  // sin(x)/x
  */
  real_type Sinc( real_type x );
  real_type Sinc_D( real_type x );
  real_type Sinc_DD( real_type x );
  real_type Sinc_DDD( real_type x );

  /*
  // (1-cos(x))/x
  */
  real_type Cosc( real_type x );
  real_type Cosc_D( real_type x );
  real_type Cosc_DD( real_type x );
  real_type Cosc_DDD( real_type x );

  /*
  // atan(x)/x
  */
  real_type Atanc( real_type x );
  real_type Atanc_D( real_type x );
  real_type Atanc_DD( real_type x );
  real_type Atanc_DDD( real_type x );

  //! Add or remove multiple of \f$ 2\pi \f$ to an angle  in order to put it in the range \f$ [-\pi,\pi]\f$.
  void rangeSymm( real_type & ang );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  projectPointOnLine( real_type x0,
                      real_type y0,
                      real_type c0, //!< cos(theta0)
                      real_type s0, //!< sin(theta0)
                      real_type x,
                      real_type y ) {
    real_type dx = x - x0;
    real_type dy = y - y0;
    return (s0 * dy + c0 * dx);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  real_type
  projectPointOnCircle( real_type x0,
                        real_type y0,
                        real_type c0, //!< cos(theta0)
                        real_type s0, //!< sin(theta0)
                        real_type k,
                        real_type L,
                        real_type qx,
                        real_type qy ) {
    real_type dx  = x0 - qx;
    real_type dy  = y0 - qy;
    real_type a0  = c0 * dy - s0 * dx;
    real_type b0  = s0 * dy + c0 * dx;
    real_type tmp = a0*k;

    if ( 1+2*tmp > 0 ) {

      tmp = b0/(1+tmp);
      tmp *= -Atanc(tmp*k); // lunghezza

      if ( tmp < 0 ) {
        real_type absk = std::abs(k);
        // if 2*pi*R + tmp <= L add 2*pi*R  to the solution
        if ( m_2pi <= absk*(L-tmp) ) tmp += m_2pi / absk;
      }

      return tmp;

    } else {

      real_type om = atan2( b0, a0+1/k );
      if ( k < 0 ) om += m_pi;
      real_type ss = -om/k;
      real_type t  = m_2pi/std::abs(k);
      if      ( ss < 0 ) ss += t;
      else if ( ss > t ) ss += t;
      return ss;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  inline
  bool
  pointInsideCircle( real_type x0,
                     real_type y0,
                     real_type c0, //!< cos(theta0)
                     real_type s0, //!< sin(theta0)
                     real_type k,
                     real_type qx,
                     real_type qy ) {
    real_type cx  = x0 - s0/k;
    real_type cy  = y0 + c0/k;
    real_type dst = hypot( qx - cx, qy - cy );
    return dst*k <= 1;
  }

  /*\
   |   ____        _           ____       ____
   |  / ___|  ___ | |_   _____|___ \__  _|___ \
   |  \___ \ / _ \| \ \ / / _ \ __) \ \/ / __) |
   |   ___) | (_) | |\ V /  __// __/ >  < / __/
   |  |____/ \___/|_| \_/ \___|_____/_/\_\_____|
  \*/

  class Solve2x2 {
    int_type  i[2], j[2];
    real_type LU[2][2];
    real_type epsi;
    bool      singular;

  public:
  
    Solve2x2() : epsi(1e-10) {}
    bool factorize( real_type A[2][2] );
    bool solve( real_type const b[2], real_type x[2] ) const;
  };

  /*!
  //  return +1 = CounterClockwise
  //  return -1 = Clockwise
  //  return  0 = flat
  //
  //  CounterClockwise:
  //    the path P1->P2->P3 turns Counter-Clockwise, i.e.,
  //    the point P3 is located "on the left" of the line P1-P2.
  //  Clockwise:
  //    the path turns Clockwise, i.e.,
  //    the point P3 lies "on the right" of the line P1-P2.
  //  flat:
  //    the point P3 is located on the line segment [P1 P2].
  //
  //  Algorithm from FileExchage geom2d adapated from Sedgewick's book.
  */
  int_type
  isCounterClockwise( real_type const P1[2],
                      real_type const P2[2],
                      real_type const P3[2] );
  /*!
  //  return +1 = Inside
  //  return -1 = Outsize
  //  return  0 = on border
  */
  int_type
  isPointInTriangle( real_type const pt[2],
                     real_type const P1[2],
                     real_type const P2[2],
                     real_type const P3[2] );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  updateInterval( int_type      & lastInterval,
                  real_type       x,
                  real_type const Xvec[],
                  int_type        npts );

}

#endif

///
/// eof: G2lib.hh
///
