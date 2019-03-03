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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
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
  extern real_type const sqrtMachepsi;
  extern real_type const m_pi;        // pi
  extern real_type const m_pi_2;      // pi/2
  extern real_type const m_2pi;       // 2*pi
  extern real_type const m_1_pi;      // 1/pi
  extern real_type const m_1_sqrt_pi; // 1/sqrt(pi)
  extern bool            intersect_with_AABBtree;

  static
  inline
  void
  noAABBtree()
  { intersect_with_AABBtree = false; }

  static
  inline
  void
  yesAABBtree()
  { intersect_with_AABBtree = true; }

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

  real_type
  projectPointOnArc(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type k,
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

  /*!
   | Solve the nonlinear system
   |
   | \f[ A x + B y = C \f]
   | \f[ a x^2 + b y^2 = c \f]
   |
   | \param A first parameter of the linear equation
   | \param B second parameter of the linear equation
   | \param C third parameter of the linear equation
   | \param a first parameter of the quadratic equation
   | \param b second parameter of the quadratic equation
   | \param c third parameter of the quadratic equation
   | \return the number of solution 0, 1 or 2
   |
  \*/
  int_type
  solveLinearQuadratic( real_type A,
                        real_type B,
                        real_type C,
                        real_type a,
                        real_type b,
                        real_type c,
                        real_type x[],
                        real_type y[] );

  /*!
   | Solve the nonlinear system
   |
   | \f[ A x + B y = C \f]
   | \f[ x^2 + y^2 = 1 \f]
   |
   | \param A first parameter of the linear equation
   | \param B second parameter of the linear equation
   | \param C third parameter of the linear equation
   | \return the number of solution 0, 1 or 2
   |
  \*/
  int_type
  solveLinearQuadratic2( real_type A,
                         real_type B,
                         real_type C,
                         real_type x[],
                         real_type y[] );

  /*!
   | Intersect the parametric arc
   |
   | \f[ x = x_1+\frac{\sin(\kappa_1 s+\theta_1)-sin(\theta_1)}{\kappa_1} \f]
   | \f[ y = y_1+\frac{\cos(\theta_1)-\cos(\kappa_1 s+\theta_1)}{\kappa_1} \f]
   |
   | with the parametric arc
   | \f[ x = x_2+\frac{\sin(\kappa_2 s+\theta_2)-sin(\theta_2)}{\kappa_2} \f]
   | \f[ y = y_2+\frac{\cos(\theta_2)-\cos(\kappa_2 s+\theta_2)}{\kappa_2} \f]
   |
   | \param x1     x-origin of the first arc
   | \param y1     y-origin of the first arc
   | \param theta1 initial angle of the first arc
   | \param kappa1 curvature of the first arc
   | \param x2     x-origin of the second arc
   | \param y2     y-origin of the second arc
   | \param theta2 initial angle of the second arc
   | \param kappa2 curvature of the second arc
   | \param s1     parameter2 of intersection for the first circle arc
   | \param s2     parameter2 of intersection for the second circle arc
   |
   | \return the number of solution 0, 1 or 2
   |
  \*/

  int_type
  intersectCircleCircle( real_type x1,
                         real_type y1,
                         real_type theta1,
                         real_type kappa1,
                         real_type x2,
                         real_type y2,
                         real_type theta2,
                         real_type kappa2,
                         real_type s1[],
                         real_type s2[] );

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
    G2LIB_LINE=0,
    G2LIB_POLYLINE,
    G2LIB_CIRCLE,
    G2LIB_BIARC,
    G2LIB_CLOTHOID,
    G2LIB_CLOTHOID_LIST
  } CurveType;

  extern char const *CurveType_name[];

  typedef std::pair<real_type,real_type> Ipair;
  typedef std::vector<Ipair>             IntersectList;

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  class BaseCurve;

  bool
  collision( BaseCurve const & C1,
             BaseCurve const & C2 );

  bool
  collision( BaseCurve const & C1,
             real_type         offs_C1,
             BaseCurve const & C2,
             real_type         offs_C2 );

  void
  intersect( BaseCurve const & C1,
             BaseCurve const & C2,
             IntersectList   & ilist,
             bool              swap_s_vals );

  void
  intersect( BaseCurve const & C1,
             real_type         offs_C1,
             BaseCurve const & C2,
             real_type         offs_C2,
             IntersectList   & ilist,
             bool              swap_s_vals );

  class BaseCurve {

    // block default constructor
    BaseCurve( BaseCurve const & );

  protected:
    CurveType _type;

  public:

    BaseCurve( CurveType const & __type )
    : _type(__type)
    {}

    virtual
    ~BaseCurve() {}

    CurveType
    type() const
    { return _type; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    length() const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    length( real_type offs ) const G2LIB_PURE_VIRTUAL;

    /*\
     |   _     _
     |  | |__ | |__   _____  __
     |  | '_ \| '_ \ / _ \ \/ /
     |  | |_) | |_) | (_) >  <
     |  |_.__/|_.__/ \___/_/\_\
    \*/

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

    /*\
     |   ____             _          _______           _
     |  | __ )  ___  __ _(_)_ __    / / ____|_ __   __| |
     |  |  _ \ / _ \/ _` | | '_ \  / /|  _| | '_ \ / _` |
     |  | |_) |  __/ (_| | | | | |/ / | |___| | | | (_| |
     |  |____/ \___|\__, |_|_| |_/_/  |_____|_| |_|\__,_|
     |              |___/
    \*/

    virtual real_type thetaBegin() const;
    virtual real_type thetaEnd() const;

    virtual real_type kappaBegin() const;
    virtual real_type kappaEnd() const;

    virtual real_type xBegin() const;
    virtual real_type yBegin() const;
    virtual real_type xEnd() const;
    virtual real_type yEnd() const;

    virtual real_type xBegin( real_type offs ) const;
    virtual real_type yBegin( real_type offs ) const;
    virtual real_type xEnd( real_type offs ) const;
    virtual real_type yEnd( real_type offs ) const;

    virtual real_type tx_Begin() const;
    virtual real_type ty_Begin() const;
    virtual real_type tx_End() const;
    virtual real_type ty_End() const;

    virtual real_type nx_Begin() const;
    virtual real_type ny_Begin() const;
    virtual real_type nx_End() const;
    virtual real_type ny_End() const;

    /*\
     |  _   _          _
     | | |_| |__   ___| |_ __ _
     | | __| '_ \ / _ \ __/ _` |
     | | |_| | | |  __/ || (_| |
     |  \__|_| |_|\___|\__\__,_|
    \*/

    virtual
    real_type
    theta( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    theta_D( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    theta_DD( real_type s ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    theta_DDD( real_type s ) const G2LIB_PURE_VIRTUAL;

    /*\
     |   _
     |  | | ____ _ _ __  _ __   __ _
     |  | |/ / _` | '_ \| '_ \ / _` |
     |  |   < (_| | |_) | |_) | (_| |
     |  |_|\_\__,_| .__/| .__/ \__,_|
     |            |_|   |_|
    \*/

    real_type
    kappa( real_type s ) const
    { return theta_D(s); }

    real_type
    kappa_D( real_type s ) const
    { return theta_DD(s); }

    real_type
    kappa_DD( real_type s ) const
    { return theta_DDD(s); }

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    virtual
    real_type
    tx( real_type s ) const;

    virtual
    real_type
    ty( real_type s ) const;

    virtual
    real_type
    tx_D( real_type s ) const;

    virtual
    real_type
    ty_D( real_type s ) const;

    virtual
    real_type
    tx_DD( real_type s ) const;

    virtual
    real_type
    ty_DD( real_type s ) const;

    virtual
    real_type
    tx_DDD( real_type s ) const;

    virtual
    real_type
    ty_DDD( real_type s ) const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type
    nx( real_type s ) const
    { return ty(s); }

    real_type
    nx_D( real_type s ) const
    { return ty_D(s); }

    real_type
    nx_DD( real_type s ) const
    { return ty_DD(s); }

    real_type
    nx_DDD( real_type s ) const
    { return ty_DDD(s); }

    real_type
    ny( real_type s ) const
    { return -tx(s); }

    real_type
    ny_D( real_type s ) const
    { return -tx_D(s); }

    real_type
    ny_DD( real_type s ) const
    { return -tx_DD(s); }

    real_type
    ny_DDD( real_type s ) const
    { return -tx_DDD(s); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    tg( real_type   s,
        real_type & tg_x,
        real_type & tg_y ) const {
      tg_x = this->tx(s);
      tg_y = this->ty(s);
    }

    virtual
    void
    tg_D( real_type   s,
          real_type & tg_x_D,
          real_type & tg_y_D ) const {
      tg_x_D = this->tx_D(s);
      tg_y_D = this->ty_D(s);
    }

    virtual
    void
    tg_DD( real_type   s,
           real_type & tg_x_DD,
           real_type & tg_y_DD ) const {
      tg_x_DD = this->tx_DD(s);
      tg_y_DD = this->ty_DD(s);
    }

    virtual
    void
    tg_DDD( real_type   s,
            real_type & tg_x_DDD,
            real_type & tg_y_DDD ) const {
      tg_x_DDD = this->tx_DDD(s);
      tg_y_DDD = this->ty_DDD(s);
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    nor( real_type s, real_type & nx, real_type & ny ) const {
      tg( s, ny, nx ); ny = -ny;
    }

    void
    nor_D( real_type s, real_type & nx_D, real_type & ny_D ) const {
      tg_D( s, ny_D, nx_D ); ny_D = -ny_D;
    }

    virtual
    void
    nor_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const {
      tg_DD( s, ny_DD, nx_DD ); ny_DD = -ny_DD;
    }

    virtual
    void
    nor_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const {
      tg_DDD( s, ny_DDD, nx_DDD ); ny_DDD = -ny_DDD;
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    evaluate( real_type   s,
              real_type & th,
              real_type & k,
              real_type & x,
              real_type & y ) const {
      eval( s, x, y );
      th = theta( s );
      k  = theta_D( s );
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    evaluate( real_type   s,
              real_type   offs,
              real_type & th,
              real_type & k,
              real_type & x,
              real_type & y ) const {
      eval( s, x, y );
      th = theta( s );
      k  = theta_D( s );
      k /= 1-offs*k; // scale curvature
    }

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
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    virtual
    void
    translate( real_type tx, real_type ty ) G2LIB_PURE_VIRTUAL;

    virtual
    void
    rotate( real_type angle, real_type cx, real_type cy ) G2LIB_PURE_VIRTUAL;

    virtual
    void
    scale( real_type sc ) G2LIB_PURE_VIRTUAL;

    virtual
    void
    reverse() G2LIB_PURE_VIRTUAL;

    virtual
    void
    changeOrigin( real_type newx0, real_type newy0 ) G2LIB_PURE_VIRTUAL;

    virtual
    void
    trim( real_type s_begin, real_type s_end ) G2LIB_PURE_VIRTUAL;

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    bool
    collision( BaseCurve const & C ) const
    { return G2lib::collision( *this, C ); }

    bool
    collision( real_type         offs,
               BaseCurve const & C,
               real_type         offs_C ) const
    { return G2lib::collision( *this, offs, C, offs_C ); }

    void
    intersect( BaseCurve const & C,
               IntersectList   & ilist,
               bool              swap_s_vals ) const
    { G2lib::intersect( *this, C, ilist, swap_s_vals ); }

    void
    intersect( real_type         offs,
               BaseCurve const & C,
               real_type         offs_C,
               IntersectList   & ilist,
               bool              swap_s_vals ) const
    { G2lib::intersect( *this, offs, C, offs_C, ilist, swap_s_vals ); }

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    /*!
     | \param  qx  x-coordinate of the point
     | \param  qy  y-coordinate of the point
     | \param  x   x-coordinate of the projected point on the curve
     | \param  y   y-coordinate of the projected point on the curve
     | \param  s   parameter on the curve of the projection
     | \param  t   curvilinear coordinate of the point x,y (if orthogonal projection)
     | \param  dst distance point projected point
     | \return 1 = point is projected orthogonal
     |         0 = more than one projection (first returned)
     |        -1 = minimum point is not othogonal projection to curve
    \*/
    virtual
    int_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type & x,
                  real_type & y,
                  real_type & s,
                  real_type & t,
                  real_type & dst ) const G2LIB_PURE_VIRTUAL;

    /*!
     | \param  qx  x-coordinate of the point
     | \param  qy  y-coordinate of the point
     | \param  offs offset of the curve
     | \param  x   x-coordinate of the projected point on the curve
     | \param  y   y-coordinate of the projected point on the curve
     | \param  s   parameter on the curve of the projection
     | \param  t   curvilinear coordinate of the point x,y (if orthogonal projection)
     | \param  dst distance point projected point
     | \return 1 = point is projected orthogonal
     |         0 = more than one projection (first returned)
     |        -1 = minimum point is not othogonal projection to curve
    \*/
    virtual
    int_type // true if projection is unique and orthogonal
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type   offs,
                  real_type & x,
                  real_type & y,
                  real_type & s,
                  real_type & t,
                  real_type & dst ) const G2LIB_PURE_VIRTUAL;

    virtual
    real_type
    distance( real_type qx, real_type qy ) const {
      real_type x, y, s, t, dst;
      closestPoint( qx, qy, x, y, s, t, dst );
      return dst;
    }

    virtual
    real_type
    distance( real_type qx,
              real_type qy,
              real_type offs ) const {
      real_type x, y, s, t, dst;
      closestPoint( qx, qy, offs, x, y, s, t, dst );
      return dst;
    }

    /*\
     |    __ _           _ ____ _____
     |   / _(_)_ __   __| / ___|_   _|
     |  | |_| | '_ \ / _` \___ \ | |
     |  |  _| | | | | (_| |___) || |
     |  |_| |_|_| |_|\__,_|____/ |_|
    \*/

    bool
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const {
      real_type X, Y, dst;
      int_type icode = closestPoint( x, y, X, Y, s, t, dst );
      return icode >= 0;
    }

    virtual
    void
    info( ostream_type & stream ) const G2LIB_PURE_VIRTUAL;

  };

}

#endif

///
/// eof: G2lib.hh
///
