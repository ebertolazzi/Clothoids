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

#include <vector>
#include <utility>

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

#define G2LIB_PURE_VIRTUAL = 0
#ifdef G2LIB_USE_CXX11
  #define G2LIB_OVERRIDE override
#else
  #define G2LIB_OVERRIDE
#endif

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wc++98-compat"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored  "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
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
  void
  minmax3( real_type   a,
           real_type   b,
           real_type   c,
           real_type & vmin,
           real_type & vmax ) {
    vmin = vmax = a;
    if ( b < vmin ) vmin = b;
    else            vmax = b;
    if ( c < vmin ) vmin = c;
    else if ( c > vmax ) vmax = c;
  }

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

  real_type
  projectPointOnCircle( real_type x0,
                        real_type y0,
                        real_type c0, //!< cos(theta0)
                        real_type s0, //!< sin(theta0)
                        real_type k,
                        real_type L,
                        real_type qx,
                        real_type qy );

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
   |  ____                  ____
   | | __ )  __ _ ___  ___ / ___|   _ _ ____   _____
   | |  _ \ / _` / __|/ _ \ |  | | | | '__\ \ / / _ \
   | | |_) | (_| \__ \  __/ |__| |_| | |   \ V /  __/
   | |____/ \__,_|___/\___|\____\__,_|_|    \_/ \___|
  \*/

  typedef enum {
    G2LIB_LINE,
    G2LIB_POLYLINE,
    G2LIB_CIRCLE,
    G2LIB_BIARC,
    G2LIB_CLOTHOID,
    G2LIB_CLOTHOID_LIST
  } CurveType;

  class BaseCurve {

    // block default constructor
    BaseCurve( BaseCurve const & );

  protected:
    CurveType _type;

    typedef std::pair<real_type,real_type> Ipair;
    typedef std::vector<Ipair>             IntersectList;

  public:

    BaseCurve( CurveType const & __type )
    : _type(__type)
    {}

    virtual
    ~BaseCurve() {}

    CurveType
    type() const
    { return _type; }

    virtual
    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    bbox( real_type   offs,
          real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    length() const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    length( real_type offs ) const G2LIB_PURE_VIRTUAL;

    virtual real_type xBegin() const;
    virtual real_type yBegin() const;
    virtual real_type xEnd()   const;
    virtual real_type yEnd()   const;

    virtual real_type xBegin( real_type offs ) const;
    virtual real_type yBegin( real_type offs ) const;
    virtual real_type xEnd( real_type offs )   const;
    virtual real_type yEnd( real_type offs )   const;

    virtual real_type tx_Begin() const;
    virtual real_type ty_Begin() const;
    virtual real_type tx_End()   const;
    virtual real_type ty_End()   const;

    virtual real_type nx_Begin() const;
    virtual real_type ny_Begin() const;
    virtual real_type nx_End()   const;
    virtual real_type ny_End()   const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    X( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    Y( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    X_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    Y_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    X_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    Y_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    X_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    Y_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    eval_DD( real_type   s,
             real_type & x_DD,
             real_type & y_DD ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    eval_DDD( real_type   s,
              real_type & x_DDD,
              real_type & y_DDD ) const G2LIB_PURE_VIRTUAL;

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    virtual
    real_type
    nx( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ny( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    nx_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ny_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    nx_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ny_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    nx_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ny_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    tx( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ty( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    tx_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ty_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    tx_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ty_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    tx_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    ty_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    nor( real_type s, real_type n[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    nor_D( real_type s, real_type n_D[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    nor_DD( real_type s, real_type n_DD[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    nor_DDD( real_type s, real_type n_DDD[2] ) const G2LIB_PURE_VIRTUAL;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    tan( real_type s, real_type t[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    tan_D( real_type s, real_type t_D[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    tan_DD( real_type s, real_type t_DD[2] ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    tan_DDD( real_type s, real_type t_DDD[2] ) const G2LIB_PURE_VIRTUAL;

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    virtual
    real_type
    X( real_type s, real_type offs ) const;

    virtual
    real_type
    Y( real_type s, real_type offs ) const;

    virtual
    real_type
    X_D( real_type s, real_type offs ) const;

    virtual
    real_type
    Y_D( real_type s, real_type offs ) const;

    virtual
    real_type
    X_DD( real_type s, real_type offs ) const;

    virtual
    real_type
    Y_DD( real_type s, real_type offs ) const;

    virtual
    real_type
    X_DDD( real_type s, real_type offs ) const;

    virtual
    real_type
    Y_DDD( real_type s, real_type offs ) const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    eval( real_type   s,
          real_type   offs,
          real_type & x,
          real_type & y ) const;

    virtual
    void
    eval_D( real_type   s,
            real_type   offs,
            real_type & x_D,
            real_type & y_D ) const;

    virtual
    void
    eval_DD( real_type   s,
             real_type   offs,
             real_type & x_DD,
             real_type & y_DD ) const;

    virtual
    void
    eval_DDD( real_type   s,
              real_type   offs,
              real_type & x_DDD,
              real_type & y_DDD ) const;

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    virtual
    bool
    collision( BaseCurve const & ) const G2LIB_PURE_VIRTUAL;

    virtual
    bool
    collision( real_type         offs,
               BaseCurve const & obj,
               real_type         offs_obj ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    intersect( BaseCurve const & obj,
               IntersectList   & ilist ) const G2LIB_PURE_VIRTUAL;

    virtual
    void
    intersect( real_type         offs,
               BaseCurve const & obj,
               real_type         offs_obj,
               IntersectList   & ilist ) const G2LIB_PURE_VIRTUAL;
    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    virtual
    real_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type & x,
                  real_type & y,
                  real_type & s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type   offs,
                  real_type & x,
                  real_type & y,
                  real_type & s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    distance( real_type qx, real_type qy ) const {
      real_type x, y, s;
      return closestPoint( qx, qy, x, y, s );
    }

    virtual
    real_type
    distance( real_type qx,
              real_type qy,
              real_type offs ) const {
      real_type x, y, s;
      return closestPoint( qx, qy, offs, x, y, s );
    }

    /*!
     | \param  qx  x-coordinate of the point
     | \param  qy  y-coordinate of the point
     | \param  x   x-coordinate of the projected point on the curve
     | \param  y   y-coordinate of the projected point on the curve
     | \param  s   parameter on the curve of the projection
     | \return 1  = unique orthogonal projection
     |         0  = more than one projection (first returned)
     |         -1 = projection line not othogonal to curve
     |         -2 = projection line not othogonal andnot unique
    \*/
    virtual
    int_type
    projection( real_type   qx,
                real_type   qy,
                real_type & x,
                real_type & y,
                real_type & s ) const G2LIB_PURE_VIRTUAL;

    /*!
     | \param  qx   x-coordinate of the point
     | \param  qy   y-coordinate of the point
     | \param  offs offset of the curve
     | \param  x    x-coordinate of the projected point on the curve
     | \param  y    y-coordinate of the projected point on the curve
     | \param  s    parameter on teh curve of the projection
     | \return 1  = unique orthogonal projection
     |         0  = more than one projection (first returned)
     |         -1 = projection line not othogonal to curve
     |         -2 = projection line not othogonal andnot unique
    \*/
    virtual
    int_type // true if projection is unique and orthogonal
    projection( real_type   qx,
                real_type   qy,
                real_type   offs,
                real_type & x,
                real_type & y,
                real_type & s ) const G2LIB_PURE_VIRTUAL;

    /*\
     |    __ _           _ ____ _____
     |   / _(_)_ __   __| / ___|_   _|
     |  | |_| | '_ \ / _` \___ \ | |
     |  |  _| | | | | (_| |___) || |
     |  |_| |_|_| |_|\__,_|____/ |_|
    \*/

    virtual
    bool
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const G2LIB_PURE_VIRTUAL;

    /*\
     |   _   _ _   _ ____  ____ ____
     |  | \ | | | | |  _ \| __ ) ___|
     |  |  \| | | | | |_) |  _ \___ \
     |  | |\  | |_| |  _ <| |_) |__) |
     |  |_| \_|\___/|_| \_\____/____/
    \*/

    /*!
     | \param n_knots  the number of knots
     | \param n_pnts   the number of polygon points
    \*/

    virtual
    void
    paramNURBS( int_type & n_knots,
                int_type & n_pnts ) const G2LIB_PURE_VIRTUAL;

    /*!
     | \brief Compute rational B-spline coefficients for a line segment
     |
     | \param knots  knots of the B-spline
     | \param Poly   polygon of the B-spline
    \*/

    virtual
    void
    toNURBS( real_type knots[],
             real_type Poly[][3] ) const G2LIB_PURE_VIRTUAL;


    /*!
     | \brief Compute B-spline coefficients for a line segment
     |
     | \param knots  knots of the B-spline
     | \param Poly   polygon of the B-spline
    \*/

    virtual
    void
    toBS( real_type knots[],
          real_type Poly[][2] ) const G2LIB_PURE_VIRTUAL;

  };
}

#endif

///
/// eof: G2lib.hh
///
