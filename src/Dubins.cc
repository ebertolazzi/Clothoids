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
/// file: Dubins.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "PolynomialRoots.hh"

namespace G2lib {

  using PolynomialRoots::Quadratic;

  /*\
   |   ____        _     _
   |  |  _ \ _   _| |__ (_)_ __  ___
   |  | | | | | | | '_ \| | '_ \/ __|
   |  | |_| | |_| | |_) | | | | \__ \
   |  |____/ \__,_|_.__/|_|_| |_|___/
  \*/

  static
  void
  into_0_2pi( real_type & a ) {
    while ( a < 0            ) a += Utils::m_2pi;
    while ( a > Utils::m_2pi ) a -= Utils::m_2pi;
  }

  static
  void
  minus_pi_pi( real_type & a ) {
    while ( a < -Utils::m_pi ) a += Utils::m_2pi;
    while ( a >  Utils::m_pi ) a -= Utils::m_2pi;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  to_integer( DubinsType d ) {
    switch ( d ) {
    case DubinsType::LSL: return 0;
    case DubinsType::RSR: return 1;
    case DubinsType::LSR: return 2;
    case DubinsType::RSL: return 3;
    case DubinsType::LRL: return 4;
    case DubinsType::RLR: return 5;
    default: break;
    }
    return 6;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::setup( GenericContainer const & gc ) {
    string cwhere{ fmt::format("Dubins[{}]::setup( gc ):", this->name() ) };
    char const * where{ cwhere.c_str() };
    real_type x0     = gc.get_map_number("x0",     where );
    real_type y0     = gc.get_map_number("y0",     where );
    real_type theta0 = gc.get_map_number("theta0", where );
    real_type x1     = gc.get_map_number("x1",     where );
    real_type y1     = gc.get_map_number("y1",     where );
    real_type theta1 = gc.get_map_number("theta1", where );
    real_type kmax   = gc.get_map_number("kmax",   where );
    bool ok = this->build( x0, y0, theta0, x1, y1, theta1, kmax );
    UTILS_ASSERT( ok, "Dubins[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::build( LineSegment const & )   { UTILS_ERROR("can convert from LineSegment to Dubins\n"); }
  void Dubins::build( CircleArc const & )     { UTILS_ERROR("can convert from CircleArc to Dubins\n"); }
  void Dubins::build( Biarc const & )         { UTILS_ERROR("can convert from Biarc to Dubins\n"); }
  void Dubins::build( ClothoidCurve const & ) { UTILS_ERROR("can convert from ClothoidCurve to Dubins\n"); }
  void Dubins::build( PolyLine const & )      { UTILS_ERROR("can convert from PolyLine to Dubins\n"); }
  void Dubins::build( BiarcList const & )     { UTILS_ERROR("can convert from BiarcList to Dubins\n"); }
  void Dubins::build( ClothoidList const & )  { UTILS_ERROR("can convert from ClothoidList to Dubins\n"); }
  void Dubins::build( Dubins3p const & )      { UTILS_ERROR("can convert from Dubins3p to Dubins\n"); }

  void
  Dubins::build( Dubins const & DB ) {
    m_C0 = DB.m_C0;
    m_C1 = DB.m_C1;
    m_C2 = DB.m_C2;
    m_solution_type = DB.m_solution_type;
  }

  bool
  Dubins_build(
    real_type    x0,
    real_type    y0,
    real_type    theta0,
    real_type    x3,
    real_type    y3,
    real_type    theta3,
    real_type    k_max,
    DubinsType & type,
    real_type  & L1,
    real_type  & L2,
    real_type  & L3,
    real_type    grad[2]
  ) {

    using std::abs;

    // setup in standard form
    real_type dx    { x3 - x0 };
    real_type dy    { y3 - y0 };
    real_type th    { atan2( dy, dx ) };
    real_type alpha { theta0 - th };
    real_type beta  { theta3 - th };
    real_type d     { k_max * hypot( dx, dy ) };

    // metto angolo in -Pi, Pi
    minus_pi_pi( alpha );
    minus_pi_pi( beta );

    // common constants
    real_type ca   { cos(alpha) };
    real_type sa   { sin(alpha) };
    real_type cb   { cos(beta)  };
    real_type sb   { sin(beta)  };
    real_type ab   { alpha-beta };
    real_type cab  { cos(ab)    };
    real_type sab  { sin(ab)    };
    real_type sasb { sa+sb      };
    real_type cacb { ca+cb      };
    real_type dx2  { 2*d        };
    real_type dd2  { d*d/2      };
    real_type dca  { d*ca       };
    real_type dsa  { d*sa       };
    real_type dcb  { d*cb       };
    real_type dsb  { d*sb       };

    real_type a2 { 2-cacb };
    real_type b2 { 2*sasb };
    real_type c2 { 2+cacb };

    real_type min_len{ Utils::Inf<real_type>() };

    type = DubinsType::DUBINS_ERROR;

    #define CHECK_DUBINS( TYPE ) \
    real_type ll{l1+l2+l3};      \
    if ( ll < min_len ) { min_len=ll; L1=l1; L2=l2; L3=l3; type=DubinsType::TYPE; }

    auto LSL = [&]() -> void {
      real_type A{ d + sa - sb };
      real_type B{ cb - ca };
      real_type thM{ atan2( B, A ) };
      real_type l2{ A * cos(thM) + B * sin(thM) };
      if ( l2 < 0 ) { thM += Utils::m_pi; l2 = -l2; }
      real_type l1{ thM  - alpha }; into_0_2pi( l1 );
      real_type l3{ beta - thM   }; into_0_2pi( l3 );
      CHECK_DUBINS( LSL );
    };

    auto RSR = [&]() -> void {
      real_type A{ d - sa + sb };
      real_type B{ ca - cb };
      real_type thM{ atan2( B, A ) };
      real_type l2{ A * cos(thM) + B * sin(thM) };
      if ( l2 < 0 ) { thM += Utils::m_pi; l2 = -l2; }
      real_type l1{ alpha - thM }; into_0_2pi( l1 );
      real_type l3{ thM - beta  }; into_0_2pi( l3 );
      CHECK_DUBINS( RSR );
    };

    auto LSR = [&]() -> void {
      real_type X[2];
      bool good{abs(a2) > abs(c2)};

      Quadratic Q( good ? a2 : c2, b2 + dx2, good ? c2 : a2 );
      integer nr{ Q.get_real_roots( X ) };

      for ( integer ir{0}; ir < nr; ++ir ) {
        real_type Xsol{ X[ir] };
        real_type Xsol2{ Xsol*Xsol };
        real_type th{ atan2( -2*Xsol, good ? Xsol2-1:1-Xsol2 ) };
        real_type l2{ (d + sasb)*cos(th) - cacb*sin(th) };
        if ( l2 >= 0 ) {
          real_type l1{ th - alpha }; into_0_2pi( l1 );
          real_type l3{ th - beta  }; into_0_2pi( l3 );
          CHECK_DUBINS( LSR );
        }
      }
    };

    auto RSL = [&]() -> void {
      real_type X[2];
      bool good{abs(a2) > abs(c2)};

      Quadratic Q( good ? a2 : c2, b2 - dx2, good ? c2 : a2 );
      integer nr{ Q.get_real_roots( X ) };

      for ( integer ir{0}; ir < nr; ++ir ) {
        real_type Xsol{ X[ir] };
        real_type Xsol2{ Xsol*Xsol };
        real_type th{ atan2( -2*Xsol, good ? Xsol2-1:1-Xsol2 ) };
        real_type l2{ ( d - sasb)*cos(th) + cacb*sin(th) };
        if ( l2 >= 0 ) {
          real_type l1{ alpha - th }; into_0_2pi( l1 );
          real_type l3{ beta  - th }; into_0_2pi( l3 );
          CHECK_DUBINS( RSL );
        }
      }
    };

    auto LRL = [&]() -> void {

      Quadratic Q( cab-1+dd2+d*sasb, 4*(dcb+sab), 3+dd2+dsa-3*(dsb+cab) );
      real_type Z[2];
      integer nr{ Q.get_real_roots( Z ) };

      for ( integer ir{0}; ir < nr; ++ir ) {
        real_type Zsol{ Z[ir] };
        real_type l3{ Utils::m_pi + 2*atan( Zsol ) };
        real_type th{ l3 + ab };
        real_type t{ dca+sab-2*sin(th)   };
        real_type b{ dsa+1-cab+2*cos(th) };
        real_type l1{ atan2( t, b ) }; into_0_2pi( l1 );
        real_type l2{ l1 + th };       into_0_2pi( l2 );
        CHECK_DUBINS( LRL );
      }
    };

    auto RLR = [&]() -> void {

      Quadratic Q( cab-1+dd2-d*sasb, 4*(dcb-sab), 3+dd2-dsa+3*(dsb-cab) );
      real_type Z[2];
      integer nr{ Q.get_real_roots( Z ) };

      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type Zsol{ Z[ir] };
        real_type l3{ Utils::m_pi + 2*atan( Zsol ) };
        real_type th{ ab - l3 };
        real_type t{ dca-sab+2*sin(th)   };
        real_type b{ 1-dsa-cab+2*cos(th) };
        real_type l1{ atan2( t, b ) }; into_0_2pi( l1 );
        real_type l2{ l1 - th };       into_0_2pi( l2 );
        CHECK_DUBINS( RLR );
      }
    };

    LSL(); RSR();
    LSR(); RSL();
    LRL(); RLR();

    switch ( type ) {
    case DubinsType::LSL:
    case DubinsType::RSR:
    case DubinsType::LSR:
    case DubinsType::RSL:
      {
        real_type C1{cos(L1)};
        real_type C3{cos(L3)};
        switch ( type ) {
        case DubinsType::LSL: grad[0] = C1-1; grad[1] = 1-C3; break;
        case DubinsType::RSR: grad[0] = 1-C1; grad[1] = C3-1; break;
        case DubinsType::LSR: grad[0] = C1-1; grad[1] = C3-1; break;
        case DubinsType::RSL: grad[0] = 1-C1; grad[1] = 1-C3; break;
        default:                                              break;
        }
      }
      break;
    case DubinsType::LRL:
      {
        real_type t2{2*sin(ab+(L1+L3))};
        grad[0] = (dca+sab)/t2-1;
        grad[1] = 1-(dcb+sab)/t2;
      }
      break;
    case DubinsType::RLR:
      {
        real_type t2{2*sin(ab-(L1+L3))};
        grad[0] = (dca-sab)/t2+1;
        grad[1] = (sab-dcb)/t2-1;
      }
      break;
    default:
      grad[0] = grad[1] = 0;
      break;
    }

    return type != DubinsType::DUBINS_ERROR;
  }

  bool
  Dubins::build(
    real_type x0,
    real_type y0,
    real_type theta0,

    real_type x3,
    real_type y3,
    real_type theta3,

    real_type k_max
  ) {
    real_type L1, L2, L3, grad[2];
    bool ok{ Dubins_build(
      x0, y0, theta0, x3, y3, theta3, k_max,
      m_solution_type, L1, L2, L3, grad
    ) };

    if ( ok ) {
      real_type s1{1}, s2{0}, s3{1};
      switch ( m_solution_type ) {
      case DubinsType::LSL:                       break;
      case DubinsType::RSR: s1 = s3 = -1;         break;
      case DubinsType::LSR: s3 = -1;              break;
      case DubinsType::RSL: s1 = -1;              break;
      case DubinsType::LRL: s2 = -1;              break;
      case DubinsType::RLR: s1 = s3 = -1; s2 = 1; break;
      default:                                    break;
      }

      // scale derivative
      m_length        = (L1+L2+L3)/k_max;
      m_length_Dalpha = grad[0]/k_max;
      m_length_Dbeta  = grad[1]/k_max;

      // salva e riscala
      m_C0.build( x0,           y0,           theta0,           s1*k_max, L1/k_max );
      m_C1.build( m_C0.x_end(), m_C0.y_end(), m_C0.theta_end(), s2*k_max, L2/k_max );
      m_C2.build( m_C1.x_end(), m_C1.y_end(), m_C1.theta_end(), s3*k_max, L3/k_max );
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::length_ISO( real_type offs ) const {
    return m_C0.length_ISO( offs )+m_C1.length_ISO( offs )+m_C2.length_ISO( offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT(FUN)                    \
if ( s < m_C0.length() ) return m_C0.FUN(s);  \
s -= m_C0.length();                           \
if ( s < m_C1.length() ) return m_C1.FUN(s);  \
s -= m_C1.length();                           \
return m_C2.FUN(s)

  real_type Dubins::theta  ( real_type s ) const { DUBINS_SELECT( theta ); }
  real_type Dubins::theta_D( real_type s ) const { DUBINS_SELECT( theta_D ); }
  real_type Dubins::X( real_type s )       const { DUBINS_SELECT( X ); }
  real_type Dubins::X_D( real_type s )     const { DUBINS_SELECT( X_D ); }
  real_type Dubins::X_DD( real_type s )    const { DUBINS_SELECT( X_DD ); }
  real_type Dubins::X_DDD( real_type s )   const { DUBINS_SELECT( X_DDD ); }
  real_type Dubins::Y( real_type s )       const { DUBINS_SELECT( Y ); }
  real_type Dubins::Y_D( real_type s )     const { DUBINS_SELECT( Y_D ); }
  real_type Dubins::Y_DD( real_type s )    const { DUBINS_SELECT( Y_DD ); }
  real_type Dubins::Y_DDD( real_type s )   const { DUBINS_SELECT( Y_DDD ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT_EVAL(FUN,...)    \
  if ( s < m_C0.length() ) {           \
    m_C0.FUN( s, __VA_ARGS__ );        \
  } else {                             \
    s -= m_C0.length();                \
    if ( s < m_C1.length() ) {         \
      m_C1.FUN( s, __VA_ARGS__ );      \
    } else {                           \
      s -= m_C1.length();              \
      m_C2.FUN( s, __VA_ARGS__ );      \
    }                                  \
  }

  void
  Dubins::eval(
    real_type   s,
    real_type & theta,
    real_type & kappa,
    real_type & x,
    real_type & y
  ) const {
    DUBINS_SELECT_EVAL( evaluate, theta, kappa, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    DUBINS_SELECT_EVAL( eval, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    DUBINS_SELECT_EVAL( eval_D, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    DUBINS_SELECT_EVAL( eval_DD, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    DUBINS_SELECT_EVAL( eval_DDD, x_DDD, y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void
  Dubins::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO, offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_D, offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_DD, offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_DDD, offs, x_DDD, y_DDD );
  }

  void
  Dubins::reverse() {
    CircleArc TMP(m_C0); m_C0 = m_C2; m_C2 = TMP;
    m_C0.reverse();
    m_C1.reverse();
    m_C2.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::scale( real_type scl ) {
    m_C0.scale( scl );
    m_C1.scale( scl );
    m_C1.change_origin( m_C0.x_end(), m_C0.y_end() );
  }

  void
  Dubins::change_origin( real_type newx0, real_type newy0 ) {
    m_C0.change_origin(newx0,newy0);
    m_C1.change_origin(m_C0.x_end(),m_C0.y_end());
    m_C2.change_origin(m_C1.x_end(),m_C1.y_end());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::trim( real_type, real_type ) {
    UTILS_ERROR0( "Dubins::trim not defined, convert to ClothoidList to trim the curve!");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_C0.bbox( xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_C1.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
    m_C2.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::bbox_ISO(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_C0.bbox_ISO( offs, xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_C1.bbox_ISO( offs, xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
    m_C2.bbox_ISO( offs, xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  /*\
   |             _ _ _     _
   |    ___ ___ | | (_)___(_) ___  _ __
   |   / __/ _ \| | | / __| |/ _ \| '_ \
   |  | (_| (_) | | | \__ \ | (_) | | | |
   |   \___\___/|_|_|_|___/_|\___/|_| |_|
  \*/

  bool
  Dubins::collision( Dubins const & B ) const {
    return m_C0.collision( B.m_C0 ) ||
           m_C0.collision( B.m_C1 ) ||
           m_C0.collision( B.m_C2 ) ||

           m_C1.collision( B.m_C0 ) ||
           m_C1.collision( B.m_C1 ) ||
           m_C1.collision( B.m_C2 ) ||

           m_C2.collision( B.m_C0 ) ||
           m_C2.collision( B.m_C1 ) ||
           m_C2.collision( B.m_C2 );
  }

  bool
  Dubins::collision_ISO(
    real_type      offs,
    Dubins const & B,
    real_type      offs_B
  ) const {
    return m_C0.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C0.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C0.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C1.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C1.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C2.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C2, offs_B );
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Dubins::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      return this->collision( C );
    } else {
      return G2lib::collision( this, pC );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Dubins::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      return G2lib::collision_ISO( this, offs, pC, offs_C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::intersect(
    Dubins const  & B,
    IntersectList & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist02,
                  ilist10, ilist11, ilist12,
                  ilist20, ilist21, ilist22;

    m_C0.intersect( B.m_C0, ilist00 );
    m_C0.intersect( B.m_C1, ilist01 );
    m_C0.intersect( B.m_C2, ilist02 );

    m_C1.intersect( B.m_C0, ilist10 );
    m_C1.intersect( B.m_C1, ilist11 );
    m_C1.intersect( B.m_C2, ilist12 );

    m_C2.intersect( B.m_C0, ilist20 );
    m_C2.intersect( B.m_C1, ilist21 );
    m_C2.intersect( B.m_C2, ilist22 );

    real_type L0  = m_C0.length();
    real_type L1  = L0+m_C1.length();
    real_type LB0 = B.m_C0.length();
    real_type LB1 = LB0+B.m_C1.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +
                   ilist02.size() +

                   ilist10.size() +
                   ilist11.size() +
                   ilist12.size() +

                   ilist20.size() +
                   ilist21.size() +
                   ilist22.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist02 ) { it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist12 ) { it.first += L0; it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist20 ) { it.first += L1; ilist.push_back( it ); }
    for ( auto & it : ilist21 ) { it.first += L1; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist22 ) { it.first += L1; it.second += LB1; ilist.push_back( it ); }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::intersect_ISO(
    real_type       offs,
    Dubins const  & B,
    real_type       offs_B,
    IntersectList & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist02,
                  ilist10, ilist11, ilist12,
                  ilist20, ilist21, ilist22;

    m_C0.intersect_ISO( offs, B.m_C0, offs_B, ilist00 );
    m_C0.intersect_ISO( offs, B.m_C1, offs_B, ilist01 );
    m_C0.intersect_ISO( offs, B.m_C2, offs_B, ilist02 );

    m_C1.intersect_ISO( offs, B.m_C0, offs_B, ilist10 );
    m_C1.intersect_ISO( offs, B.m_C1, offs_B, ilist11 );
    m_C1.intersect_ISO( offs, B.m_C2, offs_B, ilist12 );

    m_C2.intersect_ISO( offs, B.m_C0, offs_B, ilist20 );
    m_C2.intersect_ISO( offs, B.m_C1, offs_B, ilist21 );
    m_C2.intersect_ISO( offs, B.m_C2, offs_B, ilist22 );

    real_type L0  = m_C0.length();
    real_type L1  = L0+m_C1.length();
    real_type LB0 = B.m_C0.length();
    real_type LB1 = LB0+B.m_C1.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +
                   ilist02.size() +

                   ilist10.size() +
                   ilist11.size() +
                   ilist12.size() +

                   ilist20.size() +
                   ilist21.size() +
                   ilist22.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist02 ) { it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist12 ) { it.first += L0; it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist20 ) { it.first += L1; ilist.push_back( it ); }
    for ( auto & it : ilist21 ) { it.first += L1; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist22 ) { it.first += L1; it.second += LB1; ilist.push_back( it ); }

 }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    integer res  = m_C0.closest_point_ISO( qx, qy, x,  y,  s,  t,  dst  );
    integer res1 = m_C1.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length()+m_C2.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  void
  Dubins::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      this->intersect( C, ilist );
    } else {
      G2lib::intersect( this, pC, ilist );
    }
  }

  void
  Dubins::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    integer res  = m_C0.closest_point_ISO( qx, qy, offs, x,  y,  s,  t,  dst  );
    integer res1 = m_C1.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length_ISO(offs);
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length_ISO(offs)+m_C2.length_ISO(offs);
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Dubins::info() const
  { return fmt::format( "Dubins\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, Dubins const & bi ) {
    stream
      << "C0\n" << bi.m_C0
      << "C1\n" << bi.m_C1
      << "C2\n" << bi.m_C2
      << "\n";
    return stream;
  }

}

///
/// eof: Dubins.cc
///
