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
/// file: Dubins3p.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#include "PolynomialRoots.hh"

namespace G2lib {

  Dubins3pBuildType
  string_to_Dubins3pBuildType( string const & str ) {
    map<string,Dubins3pBuildType> str_to_type {
      {"sample",Dubins3pBuildType::SAMPLE_ONE_DEGREE},
      {"pattern",Dubins3pBuildType::PATTERN_SEARCH},
      {"pattern748",Dubins3pBuildType::PATTERN_SEARCH_WITH_ALGO748},
      {"pattern_search",Dubins3pBuildType::PATTERN_SEARCH},
      {"pattern_search748",Dubins3pBuildType::PATTERN_SEARCH_WITH_ALGO748},
      {"trichotomy",Dubins3pBuildType::PATTERN_TRICHOTOMY},
      {"trichotomy748",Dubins3pBuildType::PATTERN_TRICHOTOMY_WITH_ALGO748},
      {"pattern_trichotomy",Dubins3pBuildType::PATTERN_TRICHOTOMY},
      {"pattern_trichotomy748",Dubins3pBuildType::PATTERN_TRICHOTOMY_WITH_ALGO748},
      {"poly",Dubins3pBuildType::POLYNOMIAL_SYSTEM},
      {"polynomial",Dubins3pBuildType::POLYNOMIAL_SYSTEM},
      {"ellipse",Dubins3pBuildType::ELLIPSE}
    };
    return str_to_type.at( str );
  }

  /*
  //   ____        _     _           _____
  //  |  _ \ _   _| |__ (_)_ __  ___|___ / _ __
  //  | | | | | | | '_ \| | '_ \/ __| |_ \| '_ \
  //  | |_| | |_| | |_) | | | | \__ \___) | |_) |
  //  |____/ \__,_|_.__/|_|_| |_|___/____/| .__/
  //                                      |_|
  */

  bool
  Dubins3p::build(
    real_type         xi,
    real_type         yi,
    real_type         thetai,
    real_type         xm,
    real_type         ym,
    real_type         xf,
    real_type         yf,
    real_type         thetaf,
    real_type         k_max,
    Dubins3pBuildType method
  ) {
    switch ( method ) {
    case Dubins3pBuildType::SAMPLE_ONE_DEGREE:
      return build_sample( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max );
    case Dubins3pBuildType::PATTERN_SEARCH:
      return build_pattern_search( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, m_tolerance, false, false );
    case Dubins3pBuildType::PATTERN_TRICHOTOMY:
      return build_pattern_search( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, m_tolerance, true, false );
    case Dubins3pBuildType::PATTERN_SEARCH_WITH_ALGO748:
      return build_pattern_search( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, m_tolerance, false, true );
    case Dubins3pBuildType::PATTERN_TRICHOTOMY_WITH_ALGO748:
      return build_pattern_search( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, m_tolerance, true, true );
    case Dubins3pBuildType::ELLIPSE:
      return build_ellipse( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max );
    case Dubins3pBuildType::POLYNOMIAL_SYSTEM:
      break;
    }
    return false;
  }

  bool
  Dubins3p::build_sample(
    real_type xi,
    real_type yi,
    real_type thetai,
    real_type xm,
    real_type ym,
    real_type xf,
    real_type yf,
    real_type thetaf,
    real_type k_max
  ) {
    m_evaluation = 0;
    real_type thetam{0};
    m_Dubins0.build( xi, yi, thetai, xm, ym, thetam, k_max );
    m_Dubins1.build( xm, ym, thetam, xf, yf, thetaf, k_max );
    real_type len{m_Dubins0.length()+m_Dubins1.length()};

    Dubins D0{"temporary Dubins A"};
    Dubins D1{"temporary Dubins B"};
    real_type dangle{Utils::m_2pi/m_sample_points};
    for ( thetam = dangle; thetam < Utils::m_2pi; thetam += dangle ) {
      D0.build( xi, yi, thetai, xm, ym, thetam, k_max );
      D1.build( xm, ym, thetam, xf, yf, thetaf, k_max );
      ++m_evaluation;
      real_type len1{D0.length()+D1.length()};
      if ( len1 < len ) { len = len1; m_Dubins0.copy(D0); m_Dubins1.copy(D1); }
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::setup( GenericContainer const & gc ) {
    string cwhere{ fmt::format("Dubins[{}]::setup( gc ):", this->name() ) };
    char const * where{ cwhere.c_str() };
    real_type x0     { gc.get_map_number("x0",     where ) };
    real_type y0     { gc.get_map_number("y0",     where ) };
    real_type theta0 { gc.get_map_number("theta0", where ) };
    real_type xm     { gc.get_map_number("xm",     where ) };
    real_type ym     { gc.get_map_number("ym",     where ) };
    real_type x1     { gc.get_map_number("x1",     where ) };
    real_type y1     { gc.get_map_number("y1",     where ) };
    real_type theta1 { gc.get_map_number("theta1", where ) };
    real_type kmax   { gc.get_map_number("kmax",   where ) };

    string method_str{ gc.get_map_string("method", where ) };

    Dubins3pBuildType method{ string_to_Dubins3pBuildType(method_str) };
    bool ok = this->build( x0, y0, theta0, xm, ym, x1, y1, theta1, kmax, method );
    UTILS_ASSERT( ok, "Dubins[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::set_tolerance( real_type tol ) {
    UTILS_ASSERT(
      tol > 0 && tol < 1,
      "Dubins3p::set_tolerance( tol={} ) tol must be > 0 and less than 1\n",
      tol
    );
    m_tolerance = tol;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::set_sample_angle( real_type ang ) {
    UTILS_ASSERT(
      180*ang > Utils::m_pi && 3*ang <= Utils::m_2pi,
      "Dubins3p::set_sample_angle( ang={} ) ang must be > pi/180 and less than (2/3)*pi\n",
      ang
    );
    m_sample_angle = ang;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::set_sample_points( integer npts ) {
    UTILS_ASSERT(
      npts >= 4 && npts <= 36000000,
      "Dubins3p::set_sample_points( npts={} ) npts must be >= 4 and less 36000000\n",
      npts
    );
    m_sample_points = npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::set_max_evaluation( integer max_eval ) {
    UTILS_ASSERT(
      max_eval > 0 && max_eval < 1000000,
      "Dubins3p::set_max_evaluation( max_eval={} ) max_eval must be > 0 and less than 1000000\n",
      max_eval
    );
    m_max_evaluation = max_eval;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins3p::build( LineSegment const & )   { UTILS_ERROR("can convert from LineSegment to Dubins3p\n"); }
  void Dubins3p::build( CircleArc const & )     { UTILS_ERROR("can convert from CircleArc to Dubins3p\n"); }
  void Dubins3p::build( Biarc const & )         { UTILS_ERROR("can convert from Biarc to Dubins3p\n"); }
  void Dubins3p::build( ClothoidCurve const & ) { UTILS_ERROR("can convert from ClothoidCurve to Dubins3p\n"); }
  void Dubins3p::build( PolyLine const & )      { UTILS_ERROR("can convert from PolyLine to Dubins3p\n"); }
  void Dubins3p::build( BiarcList const & )     { UTILS_ERROR("can convert from BiarcList to Dubins3p\n"); }
  void Dubins3p::build( ClothoidList const & )  { UTILS_ERROR("can convert from ClothoidList to Dubins3p\n"); }
  void Dubins3p::build( Dubins const & )        { UTILS_ERROR("can convert from Dubins to Dubins3p\n"); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins3p::length() const {
    return m_Dubins0.length() + m_Dubins1.length();
  }

  real_type
  Dubins3p::length_ISO( real_type offs ) const {
    return m_Dubins0.length_ISO( offs ) + m_Dubins1.length_ISO( offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT(FUN)                              \
if ( s < m_Dubins0.length() ) return m_Dubins0.FUN(s);  \
s -= m_Dubins0.length();                                \
return m_Dubins1.FUN(s)

  real_type Dubins3p::theta   ( real_type s ) const { DUBINS_SELECT( theta ); }
  real_type Dubins3p::theta_D ( real_type s ) const { DUBINS_SELECT( theta_D ); }
  real_type Dubins3p::X       ( real_type s ) const { DUBINS_SELECT( X ); }
  real_type Dubins3p::X_D     ( real_type s ) const { DUBINS_SELECT( X_D ); }
  real_type Dubins3p::X_DD    ( real_type s ) const { DUBINS_SELECT( X_DD ); }
  real_type Dubins3p::X_DDD   ( real_type s ) const { DUBINS_SELECT( X_DDD ); }
  real_type Dubins3p::Y       ( real_type s ) const { DUBINS_SELECT( Y ); }
  real_type Dubins3p::Y_D     ( real_type s ) const { DUBINS_SELECT( Y_D ); }
  real_type Dubins3p::Y_DD    ( real_type s ) const { DUBINS_SELECT( Y_DD ); }
  real_type Dubins3p::Y_DDD   ( real_type s ) const { DUBINS_SELECT( Y_DDD ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DUBINS_SELECT_EVAL(FUN,...)    \
  if ( s < m_Dubins0.length() ) {      \
    m_Dubins0.FUN( s, __VA_ARGS__ );   \
  } else {                             \
    s -= m_Dubins0.length();           \
    m_Dubins1.FUN( s, __VA_ARGS__ );   \
  }

  void
  Dubins3p::eval(
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
  Dubins3p::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    DUBINS_SELECT_EVAL( eval, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    DUBINS_SELECT_EVAL( eval_D, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    DUBINS_SELECT_EVAL( eval_DD, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    DUBINS_SELECT_EVAL( eval_DDD, x_DDD, y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void
  Dubins3p::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO, offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_D, offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_DD, offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    DUBINS_SELECT_EVAL( eval_ISO_DDD, offs, x_DDD, y_DDD );
  }

  void
  Dubins3p::reverse() {
    Dubins TMP{m_Dubins0}; m_Dubins0.copy(m_Dubins1); m_Dubins1.copy(TMP);
    m_Dubins0.reverse();
    m_Dubins1.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::scale( real_type scl ) {
    m_Dubins0.scale( scl );
    m_Dubins1.scale( scl );
    m_Dubins1.change_origin( m_Dubins0.x_end(), m_Dubins0.y_end() );
  }

  void
  Dubins3p::change_origin( real_type newx0, real_type newy0 ) {
    m_Dubins0.change_origin(newx0,newy0);
    m_Dubins1.change_origin( m_Dubins0.x_end(), m_Dubins0.y_end() );
  }

  void
  Dubins3p::trim( real_type, real_type ) {
    UTILS_ERROR0( "Dubins::trim not defined, convert to ClothoidList to trim the curve!");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_Dubins0.bbox( xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_Dubins1.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::bbox_ISO(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_Dubins0.bbox_ISO( offs, xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_Dubins1.bbox_ISO( offs, xmi1, ymi1, xma1, yma1 );
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
  Dubins3p::collision( Dubins3p const & B ) const {
    return m_Dubins0.collision( B.m_Dubins0 ) ||
           m_Dubins0.collision( B.m_Dubins1 ) ||
           m_Dubins1.collision( B.m_Dubins0 ) ||
           m_Dubins1.collision( B.m_Dubins1 );
  }

  bool
  Dubins3p::collision_ISO(
    real_type        offs,
    Dubins3p const & B,
    real_type        offs_B
  ) const {
    return m_Dubins0.collision_ISO( offs, B.m_Dubins0, offs_B ) ||
           m_Dubins0.collision_ISO( offs, B.m_Dubins1, offs_B ) ||

           m_Dubins1.collision_ISO( offs, B.m_Dubins0, offs_B ) ||
           m_Dubins1.collision_ISO( offs, B.m_Dubins1, offs_B );
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
  Dubins3p::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::DUBINS3P ) {
      Dubins3p const & C = *static_cast<Dubins3p const *>(pC);
      return this->collision( C );
    } else {
      return G2lib::collision( this, pC );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Dubins3p::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins3p const & C = *static_cast<Dubins3p const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      return G2lib::collision_ISO( this, offs, pC, offs_C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::intersect(
    Dubins3p const & B,
    IntersectList  & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist10, ilist11;

    m_Dubins0.intersect( B.m_Dubins0, ilist00 );
    m_Dubins0.intersect( B.m_Dubins1, ilist01 );

    m_Dubins1.intersect( B.m_Dubins0, ilist10 );
    m_Dubins1.intersect( B.m_Dubins1, ilist11 );

    real_type L0  = m_Dubins0.length();
    real_type LB0 = B.m_Dubins0.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +

                   ilist10.size() +
                   ilist11.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins3p::intersect_ISO(
    real_type         offs,
    Dubins3p const  & B,
    real_type         offs_B,
    IntersectList   & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist10, ilist11;

    m_Dubins0.intersect_ISO( offs, B.m_Dubins0, offs_B, ilist00 );
    m_Dubins0.intersect_ISO( offs, B.m_Dubins1, offs_B, ilist01 );

    m_Dubins1.intersect_ISO( offs, B.m_Dubins0, offs_B, ilist10 );
    m_Dubins1.intersect_ISO( offs, B.m_Dubins1, offs_B, ilist11 );

    real_type L0  = m_Dubins0.length();
    real_type LB0 = B.m_Dubins0.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +

                   ilist10.size() +
                   ilist11.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }
 }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins3p::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    integer res  = m_Dubins0.closest_point_ISO( qx, qy, x,  y,  s,  t,  dst  );
    integer res1 = m_Dubins1.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_Dubins0.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  void
  Dubins3p::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS3P ) {
      Dubins3p const & C = *static_cast<Dubins3p const *>(pC);
      this->intersect( C, ilist );
    } else {
      G2lib::intersect( this, pC, ilist );
    }
  }

  void
  Dubins3p::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS3P ) {
      Dubins3p const & C = *static_cast<Dubins3p const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins3p::closest_point_ISO(
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
    integer res  = m_Dubins0.closest_point_ISO( qx, qy, offs, x,  y,  s,  t,  dst  );
    integer res1 = m_Dubins1.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_Dubins0.length_ISO(offs);
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Dubins3p::info() const
  { return fmt::format( "Dubins3p\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Dubins3p::solution_type_string() const {
    return m_Dubins0.solution_type_string()+
           m_Dubins1.solution_type_string();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  Dubins3p::solution_type_string_short() const {
    return m_Dubins0.solution_type_string_short()+
           m_Dubins1.solution_type_string_short();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins3p::get_range_angles(
    real_type xi,
    real_type yi,
    real_type thetai,
    real_type xm,
    real_type ym,
    real_type xf,
    real_type yf,
    real_type thetaf,
    real_type k_max,
    real_type angles[]
  ) const {
    integer npts{0};
    npts += m_Dubins0.get_range_angles_end   ( xi, yi, thetai, xm, ym,         k_max, angles        );
    npts += m_Dubins1.get_range_angles_begin ( xm, ym,         xf, yf, thetaf, k_max, angles + npts );
    // put in [0,2pi]
    for ( integer i{0}; i < npts; ++i ) {
      real_type & a{angles[i]};
      while ( a < 0            ) a += Utils::m_2pi;
      while ( a > Utils::m_2pi ) a -= Utils::m_2pi;
    }
    std::sort( angles, angles + npts );
    return npts;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `Dubins3p` object
  //!
  //!  \param stream the output stream
  //!  \param bi     an instance of `Dubins3p` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, Dubins3p const & bi ) {
    stream
      << "tolerance      = " << bi.m_tolerance      << '\n'
      << "max_evaluation = " << bi.m_max_evaluation << '\n'
      << "evaluation     = " << bi.m_evaluation     << '\n'
      << "length         = " << bi.length()         << '\n'
      << "theta_M        = " << bi.theta3_begin()   << "\n\n"
      << "Dubins0\n"         << bi.m_Dubins0
      << "Dubins1\n"         << bi.m_Dubins1;
    return stream;
  }

}

///
/// eof: Dubins3p.cc
///
