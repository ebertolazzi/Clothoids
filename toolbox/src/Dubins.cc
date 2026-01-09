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

///
/// file: Dubins.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "PolynomialRoots.hh"

namespace G2lib
{

  using PolynomialRoots::Quadratic;
  using PolynomialRoots::Quartic;

  /*\
   |   ____        _     _
   |  |  _ \ _   _| |__ (_)_ __  ___
   |  | | | | | | | '_ \| | '_ \/ __|
   |  | |_| | |_| | |_) | | | | \__ \
   |  |____/ \__,_|_.__/|_|_| |_|___/
  \*/

  static void into_0_2pi( real_type & a )
  {
    while ( a < 0 ) a += Utils::m_2pi;
    while ( a > Utils::m_2pi ) a -= Utils::m_2pi;
  }

  static void minus_pi_pi( real_type & a )
  {
    while ( a < -Utils::m_pi ) a += Utils::m_2pi;
    while ( a > Utils::m_pi ) a -= Utils::m_2pi;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer to_integer( DubinsType const d )
  {
    switch ( d )
    {
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

  Dubins::Dubins(
    real_type         x0,
    real_type         y0,
    real_type         theta0,
    real_type         x1,
    real_type         y1,
    real_type         theta1,
    real_type         k_max,
    string_view const name )
    : BaseCurve( name )
  {
    bool const ok{ this->build( x0, y0, theta0, x1, y1, theta1, k_max ) };
    UTILS_ASSERT(
      ok,
      "Dubins::Dubins(\n"
      "  x0     = {},\n"
      "  y0     = {},\n"
      "  theta0 = {},\n"
      "  x1     = {},\n"
      "  y1     = {},\n"
      "  theta1 = {},\n"
      "  k_max  = {},\n"
      "  name   = {}\n"
      ") failed\n",
      x0,
      y0,
      theta0,
      x1,
      y1,
      theta1,
      k_max,
      name );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::setup( GenericContainer const & gc )
  {
    string const    where{ fmt::format( "Dubins[{}]::setup( gc ):", this->name() ) };
    real_type const x0{ gc.get_map_number( "x0", where ) };
    real_type const y0{ gc.get_map_number( "y0", where ) };
    real_type const theta0{ gc.get_map_number( "theta0", where ) };
    real_type const x1{ gc.get_map_number( "x1", where ) };
    real_type const y1{ gc.get_map_number( "y1", where ) };
    real_type const theta1{ gc.get_map_number( "theta1", where ) };
    real_type const kmax{ gc.get_map_number( "kmax", where ) };
    bool const      ok{ this->build( x0, y0, theta0, x1, y1, theta1, kmax ) };
    UTILS_ASSERT( ok, "Dubins[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::build( LineSegment const & )
  {
    UTILS_ERROR( "cannot convert from LineSegment to Dubins\n" );
  }
  void Dubins::build( CircleArc const & )
  {
    UTILS_ERROR( "cannot convert from CircleArc to Dubins\n" );
  }
  void Dubins::build( Biarc const & )
  {
    UTILS_ERROR( "cannot convert from Biarc to Dubins\n" );
  }
  void Dubins::build( ClothoidCurve const & )
  {
    UTILS_ERROR( "cannot convert from ClothoidCurve to Dubins\n" );
  }
  void Dubins::build( PolyLine const & )
  {
    UTILS_ERROR( "cannot convert from PolyLine to Dubins\n" );
  }
  void Dubins::build( BiarcList const & )
  {
    UTILS_ERROR( "cannot convert from BiarcList to Dubins\n" );
  }
  void Dubins::build( ClothoidList const & )
  {
    UTILS_ERROR( "cannot convert from ClothoidList to Dubins\n" );
  }
  void Dubins::build( Dubins3p const & )
  {
    UTILS_ERROR( "cannot convert from Dubins3p to Dubins\n" );
  }

  void Dubins::build( Dubins const & DB )
  {
    m_C0            = DB.m_C0;
    m_C1            = DB.m_C1;
    m_C2            = DB.m_C2;
    m_solution_type = DB.m_solution_type;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //! Given two points
  //!
  //! \f[ (x_0,y_0),\qquad (x_3,y_3) \f]
  //!
  //! two angles
  //!
  //! \f[ \theta_0,\qquad \theta_1 \f]
  //!
  //! and a maximum of curvature \f$ \kappa_{\max} \f$
  //!
  //! compute the solution of the Dubins problem
  //!
  //! \param[in]  x0     coordinate \f$ x_0 \f$
  //! \param[in]  y0     coordinate \f$ y_0 \f$
  //! \param[in]  theta0 angle \f$ \theta_0 \f$
  //! \param[in]  x3     coordinate \f$ x_3 \f$
  //! \param[in]  y3     coordinate \f$ y_3 \f$
  //! \param[in]  theta3 angle \f$ \theta_3 \f$
  //! \param[in]  k_max  max curvature \f$ \kappa_{\max} \f$
  //! \param[out] type   solution type
  //! \param[out] L1     length of first arc
  //! \param[out] L2     length of second arc
  //! \param[out] L3     length of third arc
  //! \param[out] grad   gradient of the solution respect to initial and final angle
  //! \return `true` if a solution is found
  //!
  bool Dubins_build(
    real_type const x0,
    real_type const y0,
    real_type const theta0,
    real_type const x3,
    real_type const y3,
    real_type const theta3,
    real_type const k_max,
    DubinsType &    type,
    real_type &     L1,
    real_type &     L2,
    real_type &     L3,
    real_type       grad[2] )
  {
    using std::abs;

    // setup in standard form
    real_type const dx{ x3 - x0 };
    real_type const dy{ y3 - y0 };
    real_type const th{ atan2( dy, dx ) };
    real_type       alpha{ theta0 - th };
    real_type       beta{ theta3 - th };
    real_type const d{ k_max * hypot( dx, dy ) };

    // metto angolo in -Pi, Pi
    minus_pi_pi( alpha );
    minus_pi_pi( beta );

    // common constants
    real_type const ca{ cos( alpha ) };
    real_type const sa{ sin( alpha ) };
    real_type const cb{ cos( beta ) };
    real_type const sb{ sin( beta ) };
    real_type const ab{ alpha - beta };
    real_type const cab{ cos( ab ) };
    real_type const sab{ sin( ab ) };
    real_type const sasb{ sa + sb };
    real_type const cacb{ ca + cb };
    real_type const dx2{ 2 * d };
    real_type const dd2{ d * d / 2 };
    real_type const dca{ d * ca };
    real_type const dsa{ d * sa };
    real_type const dcb{ d * cb };
    real_type const dsb{ d * sb };

    real_type const a2{ 2 - cacb };
    real_type const b2{ 2 * sasb };
    real_type const c2{ 2 + cacb };

    real_type min_len{ Utils::Inf<real_type>() };

    type = DubinsType::DUBINS_ERROR;

#define CHECK_DUBINS( TYPE )    \
  real_type ll{ l1 + l2 + l3 }; \
  if ( ll < min_len )           \
  {                             \
    min_len = ll;               \
    L1      = l1;               \
    L2      = l2;               \
    L3      = l3;               \
    type    = DubinsType::TYPE; \
  }

    auto LSL = [&]() -> void
    {
      real_type const A{ d + sa - sb };
      real_type const B{ cb - ca };
      real_type       thM{ atan2( B, A ) };
      real_type       l2{ A * cos( thM ) + B * sin( thM ) };
      if ( l2 < 0 )
      {
        thM += Utils::m_pi;
        l2 = -l2;
      }
      real_type l1{ thM - alpha };
      into_0_2pi( l1 );
      real_type l3{ beta - thM };
      into_0_2pi( l3 );
      CHECK_DUBINS( LSL );
    };

    auto RSR = [&]() -> void
    {
      real_type const A{ d - sa + sb };
      real_type const B{ ca - cb };
      real_type       thM{ atan2( B, A ) };
      real_type       l2{ A * cos( thM ) + B * sin( thM ) };
      if ( l2 < 0 )
      {
        thM += Utils::m_pi;
        l2 = -l2;
      }
      real_type l1{ alpha - thM };
      into_0_2pi( l1 );
      real_type l3{ thM - beta };
      into_0_2pi( l3 );
      CHECK_DUBINS( RSR );
    };

    auto LSR = [&]() -> void
    {
      real_type  X[2];
      bool const good{ abs( a2 ) > abs( c2 ) };

      Quadratic const Q( good ? a2 : c2, b2 + dx2, good ? c2 : a2 );
      integer const   nr{ Q.get_real_roots( X ) };

      for ( integer ir{ 0 }; ir < nr; ++ir )
      {
        real_type const Xsol{ X[ir] };
        real_type const Xsol2{ Xsol * Xsol };
        real_type const th{ atan2( -2 * Xsol, good ? Xsol2 - 1 : 1 - Xsol2 ) };
        real_type       l2{ ( d + sasb ) * cos( th ) - cacb * sin( th ) };
        if ( l2 >= 0 )
        {
          real_type l1{ th - alpha };
          into_0_2pi( l1 );
          real_type l3{ th - beta };
          into_0_2pi( l3 );
          CHECK_DUBINS( LSR );
        }
      }
    };

    auto RSL = [&]() -> void
    {
      real_type  X[2];
      bool const good{ abs( a2 ) > abs( c2 ) };

      Quadratic const Q( good ? a2 : c2, b2 - dx2, good ? c2 : a2 );
      integer const   nr{ Q.get_real_roots( X ) };

      for ( integer ir{ 0 }; ir < nr; ++ir )
      {
        real_type const Xsol{ X[ir] };
        real_type const Xsol2{ Xsol * Xsol };
        real_type const th{ atan2( -2 * Xsol, good ? Xsol2 - 1 : 1 - Xsol2 ) };
        real_type const l2{ ( d - sasb ) * cos( th ) + cacb * sin( th ) };
        if ( l2 >= 0 )
        {
          real_type l1{ alpha - th };
          into_0_2pi( l1 );
          real_type l3{ beta - th };
          into_0_2pi( l3 );
          CHECK_DUBINS( RSL );
        }
      }
    };

    auto LRL = [&]() -> void
    {
      Quadratic const Q( cab - 1 + dd2 + d * sasb, 4 * ( dcb + sab ), 3 + dd2 + dsa - 3 * ( dsb + cab ) );
      real_type       Z[2];
      integer const   nr{ Q.get_real_roots( Z ) };

      for ( integer ir{ 0 }; ir < nr; ++ir )
      {
        real_type const Zsol{ Z[ir] };
        real_type const l3{ Utils::m_pi + 2 * atan( Zsol ) };
        real_type const th{ l3 + ab };
        real_type const t{ dca + sab - 2 * sin( th ) };
        real_type const b{ dsa + 1 - cab + 2 * cos( th ) };
        real_type       l1{ atan2( t, b ) };
        into_0_2pi( l1 );
        real_type l2{ l1 + th };
        into_0_2pi( l2 );
        CHECK_DUBINS( LRL );
      }
    };

    auto RLR = [&]() -> void
    {
      Quadratic const Q( cab - 1 + dd2 - d * sasb, 4 * ( dcb - sab ), 3 + dd2 - dsa + 3 * ( dsb - cab ) );
      real_type       Z[2];
      integer const   nr{ Q.get_real_roots( Z ) };

      for ( integer ir = 0; ir < nr; ++ir )
      {
        real_type const Zsol{ Z[ir] };
        real_type const l3{ Utils::m_pi + 2 * atan( Zsol ) };
        real_type const th{ ab - l3 };
        real_type const t{ dca - sab + 2 * sin( th ) };
        real_type const b{ 1 - dsa - cab + 2 * cos( th ) };
        real_type       l1{ atan2( t, b ) };
        into_0_2pi( l1 );
        real_type l2{ l1 - th };
        into_0_2pi( l2 );
        CHECK_DUBINS( RLR );
      }
    };

    LSL();
    RSR();
    LSR();
    RSL();
    LRL();
    RLR();

    switch ( type )
    {
      case DubinsType::LSL:
      case DubinsType::RSR:
      case DubinsType::LSR:
      case DubinsType::RSL:
      {
        real_type const C1{ cos( L1 ) };
        real_type const C3{ cos( L3 ) };
        switch ( type )
        {
          case DubinsType::LSL:
            grad[0] = C1 - 1;
            grad[1] = 1 - C3;
            break;
          case DubinsType::RSR:
            grad[0] = 1 - C1;
            grad[1] = C3 - 1;
            break;
          case DubinsType::LSR:
            grad[0] = C1 - 1;
            grad[1] = C3 - 1;
            break;
          case DubinsType::RSL:
            grad[0] = 1 - C1;
            grad[1] = 1 - C3;
            break;
          default: break;
        }
      }
      break;
      case DubinsType::LRL:
      {
        real_type const t2{ 2 * sin( ab + ( L1 + L3 ) ) };
        grad[0] = ( dca + sab ) / t2 - 1;
        grad[1] = 1 - ( dcb + sab ) / t2;
      }
      break;
      case DubinsType::RLR:
      {
        real_type const t2{ 2 * sin( ab - ( L1 + L3 ) ) };
        grad[0] = ( dca - sab ) / t2 + 1;
        grad[1] = ( sab - dcb ) / t2 - 1;
      }
      break;
      default: grad[0] = grad[1] = 0; break;
    }

    return type != DubinsType::DUBINS_ERROR;
  }

  bool Dubins::build(
    real_type const x0,
    real_type const y0,
    real_type const theta0,

    real_type const x3,
    real_type const y3,
    real_type const theta3,

    real_type const k_max )
  {
    real_type  L1, L2, L3, grad[2];
    bool const ok{ Dubins_build( x0, y0, theta0, x3, y3, theta3, k_max, m_solution_type, L1, L2, L3, grad ) };

    if ( ok )
    {
      real_type s1{ 1 }, s2{ 0 }, s3{ 1 };
      switch ( m_solution_type )
      {
        case DubinsType::LSL: break;
        case DubinsType::RSR: s1 = s3 = -1; break;
        case DubinsType::LSR: s3 = -1; break;
        case DubinsType::RSL: s1 = -1; break;
        case DubinsType::LRL: s2 = -1; break;
        case DubinsType::RLR:
          s1 = s3 = -1;
          s2      = 1;
          break;
        default: break;
      }

      // scale derivative
      m_length        = ( L1 + L2 + L3 ) / k_max;
      m_length_Dalpha = grad[0] / k_max;
      m_length_Dbeta  = grad[1] / k_max;

      // salva e riscala
      m_C0.build( x0, y0, theta0, s1 * k_max, L1 / k_max );
      m_C1.build( m_C0.x_end(), m_C0.y_end(), m_C0.theta_end(), s2 * k_max, L2 / k_max );
      m_C2.build( m_C1.x_end(), m_C1.y_end(), m_C1.theta_end(), s3 * k_max, L3 / k_max );
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type Dubins::length_ISO( real_type const offs ) const
  {
    return m_C0.length_ISO( offs ) + m_C1.length_ISO( offs ) + m_C2.length_ISO( offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT( FUN )                     \
  if ( s < m_C0.length() ) return m_C0.FUN( s ); \
  s -= m_C0.length();                            \
  if ( s < m_C1.length() ) return m_C1.FUN( s ); \
  s -= m_C1.length();                            \
  return m_C2.FUN( s )

  real_type Dubins::theta( real_type s ) const
  {
    DUBINS_SELECT( theta );
  }
  real_type Dubins::theta_D( real_type s ) const
  {
    DUBINS_SELECT( theta_D );
  }
  real_type Dubins::X( real_type s ) const
  {
    DUBINS_SELECT( X );
  }
  real_type Dubins::X_D( real_type s ) const
  {
    DUBINS_SELECT( X_D );
  }
  real_type Dubins::X_DD( real_type s ) const
  {
    DUBINS_SELECT( X_DD );
  }
  real_type Dubins::X_DDD( real_type s ) const
  {
    DUBINS_SELECT( X_DDD );
  }
  real_type Dubins::Y( real_type s ) const
  {
    DUBINS_SELECT( Y );
  }
  real_type Dubins::Y_D( real_type s ) const
  {
    DUBINS_SELECT( Y_D );
  }
  real_type Dubins::Y_DD( real_type s ) const
  {
    DUBINS_SELECT( Y_DD );
  }
  real_type Dubins::Y_DDD( real_type s ) const
  {
    DUBINS_SELECT( Y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT_EVAL( FUN, ... )                       \
  if ( s < m_C0.length() ) { m_C0.FUN( s, __VA_ARGS__ ); }   \
  else                                                       \
  {                                                          \
    s -= m_C0.length();                                      \
    if ( s < m_C1.length() ) { m_C1.FUN( s, __VA_ARGS__ ); } \
    else                                                     \
    {                                                        \
      s -= m_C1.length();                                    \
      m_C2.FUN( s, __VA_ARGS__ );                            \
    }                                                        \
  }

  void Dubins::eval( real_type s, real_type & theta, real_type & kappa, real_type & x, real_type & y ) const
  {
    DUBINS_SELECT_EVAL( evaluate, theta, kappa, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval( real_type s, real_type & x, real_type & y ) const
  {
    DUBINS_SELECT_EVAL( eval, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_D( real_type s, real_type & x_D, real_type & y_D ) const
  {
    DUBINS_SELECT_EVAL( eval_D, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_DD( real_type s, real_type & x_DD, real_type & y_DD ) const
  {
    DUBINS_SELECT_EVAL( eval_DD, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_DDD( real_type s, real_type & x_DDD, real_type & y_DDD ) const
  {
    DUBINS_SELECT_EVAL( eval_DDD, x_DDD, y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void Dubins::eval_ISO( real_type s, real_type const offs, real_type & x, real_type & y ) const
  {
    DUBINS_SELECT_EVAL( eval_ISO, offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_ISO_D( real_type s, real_type const offs, real_type & x_D, real_type & y_D ) const
  {
    DUBINS_SELECT_EVAL( eval_ISO_D, offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_ISO_DD( real_type s, real_type const offs, real_type & x_DD, real_type & y_DD ) const
  {
    DUBINS_SELECT_EVAL( eval_ISO_DD, offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::eval_ISO_DDD( real_type s, real_type const offs, real_type & x_DDD, real_type & y_DDD ) const
  {
    DUBINS_SELECT_EVAL( eval_ISO_DDD, offs, x_DDD, y_DDD );
  }

  void Dubins::reverse()
  {
    CircleArc const TMP( m_C0 );
    m_C0 = m_C2;
    m_C2 = TMP;
    m_C0.reverse();
    m_C1.reverse();
    m_C2.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::scale( real_type const scl )
  {
    m_C0.scale( scl );
    m_C1.scale( scl );
    m_C1.change_origin( m_C0.x_end(), m_C0.y_end() );
  }

  void Dubins::change_origin( real_type const newx0, real_type const newy0 )
  {
    m_C0.change_origin( newx0, newy0 );
    m_C1.change_origin( m_C0.x_end(), m_C0.y_end() );
    m_C2.change_origin( m_C1.x_end(), m_C1.y_end() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::trim( real_type, real_type )
  {
    UTILS_ERROR0( "Dubins::trim not defined, convert to ClothoidList to trim the curve!" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::bbox( real_type & xmin, real_type & ymin, real_type & xmax, real_type & ymax ) const
  {
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

  void Dubins::bbox_ISO( real_type const offs, real_type & xmin, real_type & ymin, real_type & xmax, real_type & ymax )
    const
  {
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

  bool Dubins::collision( Dubins const & B ) const
  {
    return m_C0.collision( B.m_C0 ) || m_C0.collision( B.m_C1 ) || m_C0.collision( B.m_C2 ) ||

           m_C1.collision( B.m_C0 ) || m_C1.collision( B.m_C1 ) || m_C1.collision( B.m_C2 ) ||

           m_C2.collision( B.m_C0 ) || m_C2.collision( B.m_C1 ) || m_C2.collision( B.m_C2 );
  }

  bool Dubins::collision_ISO( real_type const offs, Dubins const & B, real_type const offs_B ) const
  {
    return m_C0.collision_ISO( offs, B.m_C0, offs_B ) || m_C0.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C0.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C1.collision_ISO( offs, B.m_C0, offs_B ) || m_C1.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C2.collision_ISO( offs, B.m_C0, offs_B ) || m_C2.collision_ISO( offs, B.m_C1, offs_B ) ||
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

  bool Dubins::collision( BaseCurve const * pC ) const
  {
    if ( pC->type() == CurveType::DUBINS )
    {
      Dubins const & C = *dynamic_cast<Dubins const *>( pC );
      return this->collision( C );
    }
    return G2lib::collision( this, pC );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool Dubins::collision_ISO( real_type const offs, BaseCurve const * pC, real_type const offs_C ) const
  {
    if ( pC->type() == CurveType::DUBINS )
    {
      Dubins const & C = *dynamic_cast<Dubins const *>( pC );
      return this->collision_ISO( offs, C, offs_C );
    }
    return G2lib::collision_ISO( this, offs, pC, offs_C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::intersect( Dubins const & B, IntersectList & ilist ) const
  {
    IntersectList ilist00, ilist01, ilist02, ilist10, ilist11, ilist12, ilist20, ilist21, ilist22;

    m_C0.intersect( B.m_C0, ilist00 );
    m_C0.intersect( B.m_C1, ilist01 );
    m_C0.intersect( B.m_C2, ilist02 );

    m_C1.intersect( B.m_C0, ilist10 );
    m_C1.intersect( B.m_C1, ilist11 );
    m_C1.intersect( B.m_C2, ilist12 );

    m_C2.intersect( B.m_C0, ilist20 );
    m_C2.intersect( B.m_C1, ilist21 );
    m_C2.intersect( B.m_C2, ilist22 );

    real_type const L0  = m_C0.length();
    real_type const L1  = L0 + m_C1.length();
    real_type const LB0 = B.m_C0.length();
    real_type const LB1 = LB0 + B.m_C1.length();

    ilist.reserve(
      ilist.size() +

      ilist00.size() + ilist01.size() + ilist02.size() +

      ilist10.size() + ilist11.size() + ilist12.size() +

      ilist20.size() + ilist21.size() + ilist22.size() );

    for ( auto const & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 )
    {
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist02 )
    {
      it.second += LB1;
      ilist.push_back( it );
    }

    for ( auto & it : ilist10 )
    {
      it.first += L0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist11 )
    {
      it.first += L0;
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist12 )
    {
      it.first += L0;
      it.second += LB1;
      ilist.push_back( it );
    }

    for ( auto & it : ilist20 )
    {
      it.first += L1;
      ilist.push_back( it );
    }
    for ( auto & it : ilist21 )
    {
      it.first += L1;
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist22 )
    {
      it.first += L1;
      it.second += LB1;
      ilist.push_back( it );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::intersect_ISO( real_type const offs, Dubins const & B, real_type const offs_B, IntersectList & ilist )
    const
  {
    IntersectList ilist00, ilist01, ilist02, ilist10, ilist11, ilist12, ilist20, ilist21, ilist22;

    m_C0.intersect_ISO( offs, B.m_C0, offs_B, ilist00 );
    m_C0.intersect_ISO( offs, B.m_C1, offs_B, ilist01 );
    m_C0.intersect_ISO( offs, B.m_C2, offs_B, ilist02 );

    m_C1.intersect_ISO( offs, B.m_C0, offs_B, ilist10 );
    m_C1.intersect_ISO( offs, B.m_C1, offs_B, ilist11 );
    m_C1.intersect_ISO( offs, B.m_C2, offs_B, ilist12 );

    m_C2.intersect_ISO( offs, B.m_C0, offs_B, ilist20 );
    m_C2.intersect_ISO( offs, B.m_C1, offs_B, ilist21 );
    m_C2.intersect_ISO( offs, B.m_C2, offs_B, ilist22 );

    real_type const L0  = m_C0.length();
    real_type const L1  = L0 + m_C1.length();
    real_type const LB0 = B.m_C0.length();
    real_type const LB1 = LB0 + B.m_C1.length();

    ilist.reserve(
      ilist.size() +

      ilist00.size() + ilist01.size() + ilist02.size() +

      ilist10.size() + ilist11.size() + ilist12.size() +

      ilist20.size() + ilist21.size() + ilist22.size() );

    for ( auto const & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 )
    {
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist02 )
    {
      it.second += LB1;
      ilist.push_back( it );
    }

    for ( auto & it : ilist10 )
    {
      it.first += L0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist11 )
    {
      it.first += L0;
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist12 )
    {
      it.first += L0;
      it.second += LB1;
      ilist.push_back( it );
    }

    for ( auto & it : ilist20 )
    {
      it.first += L1;
      ilist.push_back( it );
    }
    for ( auto & it : ilist21 )
    {
      it.first += L1;
      it.second += LB0;
      ilist.push_back( it );
    }
    for ( auto & it : ilist22 )
    {
      it.first += L1;
      it.second += LB1;
      ilist.push_back( it );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer Dubins::closest_point_ISO(
    real_type const qx,
    real_type const qy,
    real_type &     x,
    real_type &     y,
    real_type &     s,
    real_type &     t,
    real_type &     dst ) const
  {
    real_type x1, y1, s1, t1, dst1;
    integer   res  = m_C0.closest_point_ISO( qx, qy, x, y, s, t, dst );
    integer   res1 = m_C1.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst )
    {
      x   = x1;
      y   = y1;
      s   = s1 + m_C0.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst )
    {
      x   = x1;
      y   = y1;
      s   = s1 + m_C0.length() + m_C2.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  void Dubins::intersect( BaseCurve const * pC, IntersectList & ilist ) const
  {
    if ( pC->type() == CurveType::DUBINS )
    {
      Dubins const & C = *dynamic_cast<Dubins const *>( pC );
      this->intersect( C, ilist );
    }
    else
    {
      G2lib::intersect( this, pC, ilist );
    }
  }

  void Dubins::intersect_ISO(
    real_type const   offs,
    BaseCurve const * pC,
    real_type const   offs_C,
    IntersectList &   ilist ) const
  {
    if ( pC->type() == CurveType::DUBINS )
    {
      Dubins const & C = *dynamic_cast<Dubins const *>( pC );
      this->intersect_ISO( offs, C, offs_C, ilist );
    }
    else
    {
      G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer Dubins::closest_point_ISO(
    real_type const qx,
    real_type const qy,
    real_type const offs,
    real_type &     x,
    real_type &     y,
    real_type &     s,
    real_type &     t,
    real_type &     dst ) const
  {
    real_type x1, y1, s1, t1, dst1;
    integer   res{ m_C0.closest_point_ISO( qx, qy, offs, x, y, s, t, dst ) };
    integer   res1{ m_C1.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 ) };
    if ( dst1 < dst )
    {
      x   = x1;
      y   = y1;
      s   = s1 + m_C0.length_ISO( offs );
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst )
    {
      x   = x1;
      y   = y1;
      s   = s1 + m_C0.length_ISO( offs ) + m_C2.length_ISO( offs );
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string Dubins::info() const
  {
    return fmt::format( "Dubins\n{}\n", *this );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string_view Dubins::solution_type_string() const
  {
    return to_string( m_solution_type );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string Dubins::solution_type_string_short() const
  {
    string const tmp{ to_string( m_solution_type ) };
    string       res;
    if ( m_C0.length() > 0 ) res += tmp[0];
    if ( m_C1.length() > 0 ) res += tmp[1];
    if ( m_C2.length() > 0 ) res += tmp[2];
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define OFFSET_ROOT   \
  real_type p, dp;    \
  Q.eval( x, p, dp ); \
  x += 1e-4 * ( dp > 0 ? 1 : -1 )

  integer Dubins::get_range_angles_begin(
    real_type const x0,
    real_type const y0,
    real_type const x1,
    real_type const y1,
    real_type const theta1,
    real_type const k_max,
    real_type       angles[] ) const
  {
    // setup in standard form
    real_type const dx{ x1 - x0 };
    real_type const dy{ y1 - y0 };
    real_type const th{ atan2( dy, dx ) };
    real_type       beta{ theta1 - th };
    real_type const d{ k_max * hypot( dx, dy ) };

    // metto angolo in -Pi, Pi
    minus_pi_pi( beta );

    integer         npts{ 0 };
    real_type const sb{ sin( beta ) };
    real_type const cb{ cos( beta ) };

    {  // case CCC
      real_type const b2{ 2 * beta };
      real_type const t3{ sin( b2 ) };
      real_type const t2{ cos( b2 ) };
      real_type const t7{ d * d };
      real_type const t8{ t7 * sb };
      real_type const t9{ cb * d };
      real_type const t10{ 24 };
      real_type const t11{ ( 10 - t7 ) * t7 };
      real_type const t12{ 2 * t2 * ( 1 - t7 ) };
      real_type const t14{ 8 };
      real_type const t15{ 16 * cb };
      real_type const t16{ ( t10 * t7 - 48 ) * sb };

      real_type const s_x_d[2]{ d, -d };
      for ( integer i{ 0 }; i < 2; ++i )
      {
        real_type const t5{ s_x_d[i] };
        real_type const t6{ t5 * sb };
        real_type const t13{ t5 * ( t2 - t7 ) };

        real_type       A{ -t10 * ( t6 - cb ) + 4 * ( ( t8 + t3 ) * t5 - d * t9 ) + 26 + t11 - t12 };
        real_type       B{ t14 * ( t13 + t3 ) + t16 + t5 * ( 40 - t15 ) };
        real_type const C{ ( t14 * t7 - 16 ) * sb * t5 - 2 * ( d * d ) * t7 + 12 * t2 + 4 * t7 * ( t2 + 1 ) + 52 };
        real_type       D{ t14 * ( t13 - t3 ) + t16 + t5 * ( t15 + 40 ) };
        real_type       E{ -t10 * ( t6 + cb ) - 4 * ( ( t3 - t8 ) * t5 - d * t9 ) + 26 + t11 - t12 };

        bool const reciprocal{ std::abs( A ) >= std::abs( E ) };
        if ( reciprocal )
        {
          std::swap( A, E );
          std::swap( B, D );
        }
        Quartic const Q( A, B, C, D, E );
        real_type     X[4];
        integer const nr{ Q.get_real_roots( X ) };
        for ( integer ir{ 0 }; ir < nr; ++ir )
        {
          // convert to angle
          real_type y{ X[ir] };
          real_type x{ 1 };
          if ( reciprocal ) std::swap( y, x );
          real_type theta{ 2 * atan2( y, x ) + th };
          minus_pi_pi( theta );
          angles[npts++] = theta;
        }
      }
    }
    {  // case CSC+
      real_type const t{ d * d / 2 - 1 };
      real_type const s_x_d[2]{ d, -d };
      for ( integer i{ 0 }; i < 2; ++i )
      {
        real_type const tmp{ s_x_d[i] * sb + t };
        real_type       A{ tmp - cb };
        real_type const B{ 2 * ( s_x_d[i] + sb ) };
        real_type       C{ tmp + cb };

        bool const reciprocal{ std::abs( A ) >= std::abs( C ) };
        if ( reciprocal ) std::swap( A, C );

        Quadratic     Q( A, B, C );
        real_type     X[2];
        integer const nr{ Q.get_real_roots( X ) };
        for ( integer ir{ 0 }; ir < nr; ++ir )
        {
          // convert to angle
          real_type y{ X[ir] };
          real_type x{ 1 };
          if ( reciprocal ) std::swap( y, x );
          real_type theta{ 2 * atan2( y, x ) + th };
          minus_pi_pi( theta );
          angles[npts++] = theta;
        }
      }
    }
    return npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer Dubins::get_range_angles_end(
    real_type const x0,
    real_type const y0,
    real_type const theta0,
    real_type const x1,
    real_type const y1,
    real_type const k_max,
    real_type       angles[] ) const
  {
    // setup in standard form
    real_type const dx{ x1 - x0 };
    real_type const dy{ y1 - y0 };
    real_type const th{ atan2( dy, dx ) };
    real_type       alpha{ theta0 - th };
    real_type const d{ k_max * hypot( dx, dy ) };

    // metto angolo in -Pi, Pi
    minus_pi_pi( alpha );

    integer         npts{ 0 };
    real_type const sa{ sin( alpha ) };
    real_type const ca{ cos( alpha ) };

    {  // case CCC
      real_type const a2{ 2 * alpha };
      real_type const t3{ sin( a2 ) };
      real_type const t2{ cos( a2 ) };
      real_type const t5{ d * d };
      real_type const t6{ t5 * sa };
      real_type const t7{ ca * d };
      real_type const t10{ 24 };
      real_type const t11{ ( 10 - t5 ) * t5 };
      real_type const t12{ 2 * t2 * ( 1 - t5 ) };
      real_type const t14{ -8 };
      real_type const t15{ 16 * ca };
      real_type const t16{ ( t10 * t5 - 48 ) * sa };

      real_type const s_x_d[2]{ d, -d };
      for ( integer i{ 0 }; i < 2; ++i )
      {
        real_type const t8{ s_x_d[i] };
        real_type const t9{ t8 * sa };
        real_type const t13{ t8 * ( t2 - t5 ) };

        real_type       A{ t10 * ( t9 + ca ) + t11 - t12 - 4 * ( ( t6 + t3 ) * t8 + d * t7 ) + 26 };
        real_type       B{ t14 * ( t13 - t3 ) + t16 + t8 * ( t15 - 40 ) };
        real_type const C{ -2 * ( d * d ) * t5 + ( t14 * t5 + 16 ) * sa * t8 + 12 * t2 + 4 * t5 * ( t2 + 1 ) + 52 };
        real_type       D{ t14 * ( t13 + t3 ) + t16 + t8 * ( -t15 - 40 ) };
        real_type       E{ t10 * ( t9 - ca ) + t11 - t12 - 4 * ( ( t6 - t3 ) * t8 - d * t7 ) + 26 };

        bool const reciprocal{ std::abs( A ) >= std::abs( E ) };
        if ( reciprocal )
        {
          std::swap( A, E );
          std::swap( B, D );
        }

        Quartic const Q( A, B, C, D, E );
        real_type     X[4];
        integer const nr{ Q.get_real_roots( X ) };
        for ( integer ir{ 0 }; ir < nr; ++ir )
        {
          // convert to angle
          real_type y{ X[ir] };
          real_type x{ 1 };
          if ( reciprocal ) std::swap( y, x );
          real_type theta{ 2 * atan2( y, x ) + th };
          minus_pi_pi( theta );
          angles[npts++] = theta;
        }
      }
    }
    {  // case CSC+
      real_type const t{ d * d / 2 - 1 };
      real_type const s_x_d[2]{ d, -d };
      for ( integer i{ 0 }; i < 2; ++i )
      {
        real_type const tmp{ s_x_d[i] * sa + t };
        real_type       A{ tmp - ca };
        real_type const B{ 2 * ( s_x_d[i] + sa ) };
        real_type       C{ tmp + ca };

        bool const reciprocal{ std::abs( A ) >= std::abs( C ) };
        if ( reciprocal ) std::swap( A, C );

        Quadratic     Q( A, B, C );
        real_type     X[2];
        integer const nr{ Q.get_real_roots( X ) };
        for ( integer ir{ 0 }; ir < nr; ++ir )
        {
          // convert to angle
          real_type y{ X[ir] };
          real_type x{ 1 };
          if ( reciprocal ) std::swap( y, x );
          real_type theta{ 2 * atan2( y, x ) + th };
          minus_pi_pi( theta );
          angles[npts++] = theta;
        }
      }
    }
    return npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `Dubins` object
  //!
  //!  \param stream the output stream
  //!  \param bi     an instance of `Dubins` object
  //!  \return the output stream
  //!
  ostream_type & operator<<( ostream_type & stream, Dubins const & bi )
  {
    stream << "C0\n" << bi.m_C0 << "C1\n" << bi.m_C1 << "C2\n" << bi.m_C2 << "\n";
    return stream;
  }

}  // namespace G2lib

///
/// eof: Dubins.cc
///
