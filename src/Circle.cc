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

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

// Workaround for Visual Studio
#ifdef min
  #undef min
#endif

#ifdef max
  #undef max
#endif

#include <cmath>
#include <algorithm>

namespace G2lib {

  using std::min;
  using std::max;
  using std::abs;
  using std::tan;
  using std::abs;
  using std::ceil;
  using std::floor;
  using std::swap;
  using std::vector;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::setup( GenericContainer const & gc ) {
    string const where{ fmt::format("CircleArc[{}]::setup( gc ):", this->name() ) };
    real_type const x0     = gc.get_map_number("x0",     where );
    real_type const y0     = gc.get_map_number("y0",     where );
    real_type const theta0 = gc.get_map_number("theta0", where );
    real_type const x1     = gc.get_map_number("x1",     where );
    real_type const y1     = gc.get_map_number("y1",     where );
    bool const ok = this->build_G1( x0, y0, theta0, x1, y1 );
    UTILS_ASSERT( ok, "CircleArc[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::build( LineSegment const & LS ) {
    m_x0     = LS.x_begin();
    m_y0     = LS.y_begin();
    m_theta0 = LS.m_theta0;
    m_c0     = LS.m_c0;
    m_s0     = LS.m_s0;
    m_k      = 0;
    m_L      = LS.length();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void CircleArc::build( Biarc         const & ) { UTILS_ERROR("cannot convert from Biarc to CircleArc\n"); }
  void CircleArc::build( ClothoidCurve const & ) { UTILS_ERROR("cannot convert from ClothoidCurve to CircleArc\n"); }
  void CircleArc::build( PolyLine      const & ) { UTILS_ERROR("cannot convert from PolyLine to CircleArc\n"); }
  void CircleArc::build( BiarcList     const & ) { UTILS_ERROR("cannot convert from BiarcList to CircleArc\n"); }
  void CircleArc::build( ClothoidList  const & ) { UTILS_ERROR("cannot convert from ClothoidList to CircleArc\n"); }
  void CircleArc::build( Dubins        const & ) { UTILS_ERROR("cannot convert from Dubins to CircleArc\n"); }
  void CircleArc::build( Dubins3p      const & ) { UTILS_ERROR("cannot convert from Dubins3p to CircleArc\n"); }

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  CircleArc::CircleArc( BaseCurve const * pC ) : CircleArc( pC->name() ) {

    G2LIB_DEBUG_MESSAGE( "CircleArc convert: {}\n", pC->type_name() );

    switch ( pC->type() ) {
    case CurveType::LINE:
      G2LIB_DEBUG_MESSAGE( "LineSegment -> CircleArc\n" );
      this->build( *dynamic_cast<LineSegment const *>(pC) );
      break;
    case CurveType::CIRCLE:
      G2LIB_DEBUG_MESSAGE( "to -> CircleArc\n" );
      *this = *dynamic_cast<CircleArc const *>(pC);
      break;
    default:
      UTILS_ERROR( "CircleArc constructor cannot convert from: {}\n", pC->type_name() );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::build_G1(
    real_type const x0,
    real_type const y0,
    real_type const theta0,
    real_type const x1,
    real_type const y1
  ) {

    real_type const dx { x1 - x0 };
    real_type const dy { y1 - y0 };
    real_type const d  { hypot( dx, dy ) };

    if ( d > 0 ) {
      real_type const th = atan2( dy, dx ) - theta0;
      m_x0     = x0;
      m_y0     = y0;
      m_theta0 = theta0;
      m_k      = 2*sin(th)/d;
      m_L      = d/Sinc(th);
      return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::build_3P(
    real_type const x0,
    real_type const y0,
    real_type const x1,
    real_type const y1,
    real_type const x2,
    real_type const y2
  ) {
    real_type const dxa { x1 - x0 };
    real_type const dya { y1 - y0 };
    real_type const dxb { x2 - x1 };
    real_type const dyb { y2 - y1 };
    real_type const La  { hypot(dya,dxa) };
    real_type const Lb  { hypot(dyb,dxb) };
    real_type cosom{ (dxa*dxb + dya*dyb)/(La*Lb) };
    if      ( cosom >  1 ) cosom = 1;
    else if ( cosom < -1 ) cosom = -1;
    real_type const omega{ acos(cosom) };

    real_type       alpha { omega - atan2(Lb*sin(omega),La+Lb*cos(omega)) };
    real_type const dxc   { x2 - x0 };
    real_type const dyc   { y2 - y0 };
    real_type const Lc    { hypot(dyc,dxc) };
    real_type cosal{ (dxa*dxc + dya*dyc)/(La*Lc) };
    if ( cosal > 1 ) cosal = 1;
    else if ( cosal < -1 ) cosal = -1;
    alpha += acos(cosal);

    if ( dxa*dyb > dya*dxb ) alpha = -alpha;
    real_type const _theta0{ atan2( dyc, dxc ) + alpha };
    return build_G1( x0, y0, _theta0, x2, y2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::theta_min_max( real_type & thMin, real_type & thMax ) const  {
    thMin = m_theta0;
    thMax = m_theta0 + m_L * m_k;
    if ( thMax < thMin ) swap( thMin, thMax );
    return thMax-thMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::X( real_type const s ) const {
    real_type const sk{ (s*m_k)/2 };
    return m_x0+s*Sinc(sk)*cos(m_theta0+sk);
  }

  real_type
  CircleArc::X_D( real_type const s ) const {
    return cos(m_theta0+s*m_k);
  }

  real_type
  CircleArc::X_DD( real_type const s ) const {
    return -m_k*sin(m_theta0+s*m_k);
  }

  real_type
  CircleArc::X_DDD( real_type const s ) const {
    return -(m_k*m_k)*cos(m_theta0+s*m_k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y( real_type const s ) const {
    real_type const sk{ (s*m_k)/2 };
    return m_y0+s*Sinc(sk)*sin(m_theta0+sk);
  }

  real_type
  CircleArc::Y_D( real_type const s ) const {
    return sin(m_theta0+s*m_k);
  }

  real_type
  CircleArc::Y_DD( real_type const s ) const {
    return m_k*cos(m_theta0+s*m_k);
  }

  real_type
  CircleArc::Y_DDD( real_type const s ) const {
    return -(m_k*m_k)*sin(m_theta0+s*m_k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |  _____                   _   _   _
   | |_   _|   __ _ _ __   __| | | \ | |
   |   | |    / _` | '_ \ / _` | |  \| |
   |   | |   | (_| | | | | (_| | | |\  |
   |   |_|    \__,_|_| |_|\__,_| |_| \_|
  \*/

  void
  CircleArc::tg(
    real_type const s,
    real_type &     tx,
    real_type &     ty
  ) const {
    real_type const th{ theta(s) };
    tx = cos(th);
    ty = sin(th);
  }

  void
  CircleArc::tg_D(
    real_type const s,
    real_type &     tx_D,
    real_type &     ty_D
  ) const {
    real_type const th{ theta(s) };
    tx_D = -sin(th)*m_k;
    ty_D = cos(th)*m_k;
  }

  void
  CircleArc::tg_DD(
    real_type const s,
    real_type &     tx_DD,
    real_type &     ty_DD
  ) const {
    real_type const th { theta(s) };
    real_type const k2 { m_k*m_k };
    tx_DD = -cos(th)*k2;
    ty_DD = -sin(th)*k2;
  }

  void
  CircleArc::tg_DDD(
    real_type const s,
    real_type &     tx_DDD,
    real_type &     ty_DDD
  ) const {
    real_type const th { theta(s) };
    real_type const k3 { m_k*m_k*m_k };
    tx_DDD = sin(th)*k3;
    ty_DDD = -cos(th)*k3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval(
    real_type const s,
    real_type &     x,
    real_type &     y
  ) const {
    real_type const sk  { (s*m_k)/2 };
    real_type const LS  { s*Sinc(sk) };
    real_type const arg { m_theta0+sk };
    x = m_x0+LS*cos(arg);
    y = m_y0+LS*sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_D(
    real_type const s,
    real_type &     x_D,
    real_type &     y_D
  ) const {
    real_type const arg{ m_theta0+s*m_k };
    x_D = cos(arg);
    y_D = sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DD(
    real_type const s,
    real_type &     x_DD,
    real_type &     y_DD
  ) const {
    real_type const arg{ m_theta0+s*m_k };
    x_DD = -m_k*sin(arg);
    y_DD = m_k*cos(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DDD(
    real_type const s,
    real_type &     x_DDD,
    real_type &     y_DDD
  ) const {
    real_type const arg { m_theta0+s*m_k };
    real_type const k2  { m_k*m_k };
    x_DDD = -k2*cos(arg);
    y_DDD = -k2*sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::trim( real_type const s_begin, real_type const s_end ) {
    UTILS_ASSERT(
      s_end > s_begin,
      "CircleArc::trim( begin={}, s_end={} ) s_end must be > s_begin\n",
      s_begin, s_end
    );
    real_type x, y;
    eval( s_begin, x, y );
    m_theta0 += s_begin * m_k;
    m_L  = s_end - s_begin;
    m_x0 = x;
    m_y0 = y;
    m_c0 = cos(m_theta0);
    m_s0 = sin(m_theta0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::rotate( real_type const angle, real_type const cx, real_type const cy ) {
    real_type const dx  { m_x0 - cx };
    real_type const dy  { m_y0 - cy };
    real_type const C   { cos(angle) };
    real_type const S   { sin(angle) };
    real_type const ndx { C*dx - S*dy };
    real_type const ndy { C*dy + S*dx };
    m_x0 = cx + ndx;
    m_y0 = cy + ndy;
    m_theta0 += angle;
    m_c0      = cos(m_theta0);
    m_s0      = sin(m_theta0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::scale( real_type const s ) {
    m_k /= s;
    m_L *= s;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::reverse() {
    real_type xx, yy;
    eval( m_L, xx, yy );
    m_theta0 += m_L*m_k+Utils::m_pi;
    while ( m_theta0 >  Utils::m_pi ) m_theta0 -= Utils::m_2pi;
    while ( m_theta0 < -Utils::m_pi ) m_theta0 += Utils::m_2pi;
    m_x0 = xx;
    m_y0 = yy;
    m_c0 = cos(m_theta0);
    m_s0 = sin(m_theta0);
    m_k  = -m_k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::center( real_type & cx, real_type & cy ) const {
    real_type const nx { -sin(m_theta0) };
    real_type const ny { cos(m_theta0)  };
    cx = m_x0 + nx/m_k;
    cy = m_y0 + ny/m_k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::change_curvilinear_origin( real_type const s0, real_type const newL ) {
    real_type new_x0, new_y0;
    eval( s0, new_x0, new_y0 );
    m_x0      = new_x0;
    m_y0      = new_y0;
    m_theta0 += m_k*s0;
    m_L       = newL;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::bbTriangle(
    real_type & x0, real_type & y0,
    real_type & x1, real_type & y1,
    real_type & x2, real_type & y2
  ) const {
    real_type const dtheta { m_L * m_k };
    bool const ok{ abs(dtheta) <= Utils::m_pi/3 };
    if ( ok ) {
      x0 = m_x0;
      y0 = m_y0;
      eval( m_L, x2, y2 );
      x1 = (x0+x2)/2;
      y1 = (y0+y2)/2;
      real_type const nx { y0-y2 };
      real_type const ny { x2-x0 };
      real_type const tg { tan(dtheta/2)/2 };
      x1 -= nx * tg;
      y1 -= ny * tg;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::bbTriangle_ISO(
    real_type const offs,
    real_type & x0, real_type & y0,
    real_type & x1, real_type & y1,
    real_type & x2, real_type & y2
  ) const {
    real_type const dtheta{ m_L * m_k };
    bool const ok{ abs(dtheta) <= Utils::m_pi/3 };
    if ( ok ) {
      eval_ISO( 0,   offs, x0, y0 );
      eval_ISO( m_L, offs, x2, y2 );
      x1 = (x0+x2)/2;
      y1 = (y0+y2)/2;
      real_type const nx { y0-y2 };
      real_type const ny { x2-x0 };
      real_type const tg { tan(dtheta/2)/2 };
      x1 -= nx * tg;
      y1 -= ny * tg;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bb_triangles(
    vector<Triangle2D> & tvec,
    real_type const      max_angle,
    real_type const      max_size,
    integer   const      icurve
  ) const {
    real_type dtheta { abs( min(m_L,max_size) * m_k) };
    integer   n      { 1 };
    if ( dtheta > max_angle ) {
      n       = static_cast<integer>(ceil(dtheta / max_angle));
      dtheta /= n;
    }
    real_type tg { tan(dtheta/2)/2 };
    if ( m_k < 0 ) tg = -tg;
    tvec.reserve( static_cast<size_t>(n) );
    real_type       xx0 { m_x0 };
    real_type       yy0 { m_y0 };
    real_type const ds  { m_L/n };
    real_type       ss  { ds };
    for ( integer iter = 0; iter < n; ++iter, ss += ds ) {
      real_type xx2, yy2;
      eval( ss, xx2, yy2 );
      real_type       xx1 { (xx0+xx2)/2 };
      real_type       yy1 { (yy0+yy2)/2 };
      real_type const nx  { yy0-yy2 };
      real_type const ny  { xx2-xx0 };
      xx1 -= nx * tg;
      yy1 -= ny * tg;
      tvec.emplace_back( xx0, yy0, xx1, yy1, xx2, yy2, 0, 0, icurve );
      xx0 = xx2;
      yy0 = yy2;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bb_triangles_ISO(
    real_type const      offs,
    vector<Triangle2D> & tvec,
    real_type const      max_angle,
    real_type const      max_size,
    integer   const      icurve
  ) const {
    real_type const scale  { 1+m_k*offs };
    real_type       dtheta { abs( min(m_L,max_size/scale) * m_k ) };
    integer         n      { 1 };
    if ( dtheta > max_angle ) {
      n       = static_cast<integer>(ceil(dtheta / max_angle));
      dtheta /= n;
    }
    tvec.reserve( static_cast<size_t>(n) );
    real_type const ds { m_L/n };
    real_type       ss { ds };
    real_type       tg { scale * tan(dtheta/2)/2 };
    if ( m_k < 0 ) tg = -tg;
    real_type xx0, yy0;
    eval_ISO( 0, offs, xx0, yy0 );
    for ( integer iter = 0; iter < n; ++iter, ss += ds ) {
      real_type xx2, yy2;
      eval_ISO( ss, offs, xx2, yy2 );
      real_type       xx1 { (xx0+xx2)/2 };
      real_type       yy1 { (yy0+yy2)/2 };
      real_type const nx  { yy0-yy2 };
      real_type const ny  { xx2-xx0 };
      xx1 -= nx * tg;
      yy1 -= ny * tg;
      tvec.emplace_back( xx0, yy0, xx1, yy1, xx2, yy2, 0, 0, icurve );
      xx0 = xx2;
      yy0 = yy2;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    vector<Triangle2D> tvec;
    this->bb_triangles( tvec, Utils::m_pi/4, 1e100, 0 );
    tvec[0].bbox( xmin, ymin, xmax, ymax );
    for ( integer iter{1}; iter < static_cast<integer>(tvec.size()); ++iter ) {
      real_type xmin1, ymin1, xmax1, ymax1;
      tvec[iter].bbox( xmin1, ymin1, xmax1, ymax1 );
      if ( xmin1 < xmin ) xmin = xmin1;
      if ( ymin1 < ymin ) ymin = ymin1;
      if ( xmax1 > xmax ) xmax = xmax1;
      if ( ymax1 > ymax ) ymax = ymax1;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbox_ISO(
    real_type const offs,
    real_type &     xmin,
    real_type &     ymin,
    real_type &     xmax,
    real_type &     ymax
  ) const {
    vector<Triangle2D> tvec;
    this->bb_triangles_ISO( offs, tvec, Utils::m_pi/4, 1e100, 0);
    tvec[0].bbox( xmin, ymin, xmax, ymax );
    for ( integer iter{1}; iter < static_cast<integer>(tvec.size()); ++iter ) {
      real_type xmin1, ymin1, xmax1, ymax1;
      tvec[iter].bbox( xmin1, ymin1, xmax1, ymax1 );
      if ( xmin1 < xmin ) xmin = xmin1;
      if ( ymin1 < ymin ) ymin = ymin1;
      if ( xmax1 > xmax ) xmax = xmax1;
      if ( ymax1 > ymax ) ymax = ymax1;
    }
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
  CircleArc::collision( CircleArc const & C ) const {
    real_type s1[2], s2[2];
    integer const ni { intersectCircleCircle(
      m_x0, m_y0, m_theta0, m_k,
      C.m_x0, C.m_y0, C.m_theta0, C.m_k, s1, s2
    ) };
    real_type const eps1 { machepsi100*m_L };
    real_type const eps2 { machepsi100*C.m_L };
    for ( integer i{0}; i < ni; ++i ) {
      if ( s1[i] >= -eps1 && s1[i] <= m_L+eps1 &&
           s2[i] >= -eps2 && s2[i] <= m_L+eps2 )
        return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::CIRCLE ) {
      CircleArc const & C{ *dynamic_cast<CircleArc const *>(pC) };
      return this->collision( C );
    }
    CurveType const CT{ curve_promote( this->type(), pC->type() ) };
    if ( CT == CurveType::CIRCLE ) {
      CircleArc const C(pC);
      return this->collision( C );
    }
    return G2lib::collision( this, pC );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::collision_ISO(
    real_type const   offs,
    CircleArc const & C,
    real_type const   offs_C
  ) const {
    real_type s1[2], s2[2];
    real_type const sc1 = 1+m_k*offs;
    real_type const sc2 = 1+C.m_k*offs_C;
    integer const ni = intersectCircleCircle(
      this->X_ISO(0,offs),
      this->Y_ISO(0,offs),
      m_theta0,
      m_k/sc2,
      C.X_ISO(0,offs_C),
      C.Y_ISO(0,offs_C),
      C.m_theta0,
      C.m_k/sc2,
      s1, s2
    );
    real_type const eps1 { machepsi100*m_L };
    real_type const eps2 { machepsi100*C.m_L };
    for ( integer i{0}; i < ni; ++i ) {
      real_type const ss1{ s1[i]/sc1 };
      real_type const ss2{ s2[i]/sc2 };
      if ( ss1 >= -eps1 && ss1 <= m_L+eps1 &&
           ss2 >= -eps2 && ss2 <= C.m_L+eps2 )
        return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::collision_ISO(
    real_type const   offs,
    BaseCurve const * pC,
    real_type const   offs_C
  ) const {
    if ( pC->type() == CurveType::CIRCLE ) {
      CircleArc const & C{ *dynamic_cast<CircleArc const *>(pC) };
      return this->collision_ISO( offs, C, offs_C );
    }
    CurveType const CT{ curve_promote( this->type(), pC->type() ) };
    if ( CT == CurveType::CIRCLE ) {
      CircleArc const C(pC);
      return this->collision_ISO( offs, C, offs_C );
    }
    return G2lib::collision_ISO( this, offs, pC, offs_C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::intersect(
    CircleArc const & C,
    IntersectList   & ilist
  ) const {
    real_type s1[2], s2[2];
    integer const ni { intersectCircleCircle(
      m_x0, m_y0, m_theta0, m_k,
      C.m_x0, C.m_y0, C.m_theta0, C.m_k, s1, s2
    ) };
    real_type const eps1 { machepsi100*m_L   };
    real_type const eps2 { machepsi100*C.m_L };
    for ( integer i{0}; i < ni; ++i ) {
      real_type ss1{ s1[i] };
      real_type ss2{ s2[i] };
      if ( ss1 >= -eps1 && ss1 <= m_L+eps1 &&
           ss2 >= -eps2 && ss2 <= C.m_L+eps2 ) {
        ilist.emplace_back( ss1, ss2 );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::intersect_ISO(
    real_type const   offs,
    CircleArc const & C,
    real_type const   offs_C,
    IntersectList   & ilist
  ) const {
    real_type s1[2], s2[2];
    real_type const sc1 { 1+m_k*offs };
    real_type const sc2 { 1+C.m_k*offs_C };
    integer   const ni  { intersectCircleCircle(
      this->X_ISO(0,offs),
      this->Y_ISO(0,offs),
      m_theta0,
      m_k/sc2,
      C.X_ISO(0,offs_C),
      C.Y_ISO(0,offs_C),
      C.m_theta0,
      C.m_k/sc2,
      s1, s2
    ) };
    real_type const eps1 { machepsi100*m_L   };
    real_type const eps2 { machepsi100*C.m_L };
    for ( integer i{0}; i < ni; ++i ) {
      real_type ss1{ s1[i]/sc1 };
      real_type ss2{ s2[i]/sc2 };
      if ( ss1 >= -eps1 && ss1 <= m_L+eps1 &&
           ss2 >= -eps2 && ss2 <= C.m_L+eps2 ) {
        ilist.emplace_back( ss1, ss2 );
      }
    }
  }

  void
  CircleArc::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CIRCLE ) {
      CircleArc const & C = *dynamic_cast<CircleArc const *>(pC);
      this->intersect( C, ilist );
    } else {
      CurveType const CT{ curve_promote( this->type(), pC->type() ) };
      if ( CT == CurveType::CIRCLE ) {
        CircleArc const C(pC);
        this->intersect( C, ilist );
      } else {
        G2lib::intersect( this, pC, ilist );
      }
    }
  }

  void
  CircleArc::intersect_ISO(
    real_type const   offs,
    BaseCurve const * pC,
    real_type const   offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CIRCLE ) {
      CircleArc const & C = *dynamic_cast<CircleArc const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      CurveType const CT{ curve_promote( this->type(), pC->type() ) };
      if ( CT == CurveType::CIRCLE ) {
        CircleArc const C(pC);
        this->intersect_ISO( offs, C, offs_C, ilist );
      } else {
        G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
      }
    }
  }

  /*\
   |        _                     _   ____       _       _
   |    ___| | ___  ___  ___  ___| |_|  _ \ ___ (_)_ __ | |_
   |   / __| |/ _ \/ __|/ _ \/ __| __| |_) / _ \| | '_ \| __|
   |  | (__| | (_) \__ \  __/\__ \ |_|  __/ (_) | | | | | |_
   |   \___|_|\___/|___/\___||___/\__|_|   \___/|_|_| |_|\__|
  \*/

  integer
  CircleArc::closest_point_ISO(
    real_type const qx,
    real_type const qy,
    real_type &     x,
    real_type &     y,
    real_type &     s,
    real_type &     t,
    real_type &     dst
  ) const {
    real_type const cc0 = cos(m_theta0);
    real_type const ss0 = sin(m_theta0);
    s = projectPointOnCircleArc( m_x0, m_y0, cc0, ss0, m_k, m_L, qx, qy );
    integer res = 1;
    if ( s < 0 || s > m_L ) {
      s = m_L;
      t = 0;
      eval( s, x, y );
      // costruisco piano
      real_type const nx = x-m_x0;
      real_type const ny { y-m_y0 };
      real_type const dx { 2*qx-(m_x0+x) };
      real_type const dy { 2*qy-(m_y0+y) };
      if ( nx*dx + ny*dy <= 0 ) {
        s = 0;
        x = m_x0;
        y = m_y0;
      }
      res = -1;
    } else {
      eval( s, x, y );
    }
    real_type nx, ny;
    nor_ISO( s, nx, ny );
    real_type const dx = qx-x;
    real_type const dy = qy-y;
    t   = dx * nx + dy * ny;
    dst = hypot( dx, dy );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  CircleArc::closest_point_ISO(
    real_type const qx,
    real_type const qy,
    real_type const offs,
    real_type &     x,
    real_type &     y,
    real_type &     s,
    real_type &     t,
    real_type &     dst
  ) const  {
    real_type const cc0 = cos(m_theta0);
    real_type const ss0 = sin(m_theta0);
    real_type const xx0 = m_x0+offs*nx_begin_ISO();
    real_type const yy0 = m_y0+offs*ny_begin_ISO();
    real_type const ff  = 1+m_k*offs;
    real_type const LL  = m_L*ff;
    s = projectPointOnCircleArc( xx0, yy0, cc0, ss0, m_k/ff, LL, qx, qy );
    integer res = 1;
    if ( s < 0 || s > LL ) {
      s = m_L;
      eval_ISO( s, offs, x, y );
      // costruisco piano
      real_type const nx = x-xx0;
      real_type const ny { y-yy0 };
      real_type const dx { 2*qx-(xx0+x) };
      real_type const dy { 2*qy-(yy0+y) };
      if ( nx*dx + ny*dy <= 0 ) {
        s = 0;
        x = xx0;
        y = yy0;
      }
      res = -1;
    } else {
      eval_ISO( s, offs, x, y );
    }
    real_type nx, ny;
    nor_ISO( s, nx, ny );
    real_type const dx = qx-x;
    real_type const dy = qy-y;
    t   = dx * nx + dy * ny + offs;
    dst = hypot( dx, dy );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::paramNURBS(
    integer & n_knots,
    integer & n_pnts
  ) const {
    real_type const dtheta = m_L*m_k;
    auto ns = static_cast<integer>(floor(3 * abs(dtheta) / Utils::m_pi));
    if ( ns < 1 ) ns = 1;
    n_pnts  = 1+2*ns;
    n_knots = n_pnts+3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::toNURBS(
    real_type knots[],
    real_type Poly[][3]
  ) const {

    real_type const dtheta = m_L*m_k;
    auto ns = static_cast<integer>(floor(3 * abs(dtheta) / Utils::m_pi));
    if ( ns < 1 ) ns = 1;

    real_type const th = dtheta/(2*ns);
    real_type const w  = cos(th);
    real_type const tg = tan(th)/2;

    real_type p0[2], p2[2];
    p0[0] = m_x0;
    p0[1] = m_y0;

    knots[0] = knots[1] = knots[2] = 0;
    Poly[0][0] = p0[0];
    Poly[0][1] = p0[1];
    Poly[0][2] = 1;

    real_type       s{0};
    real_type const ds{m_L/ns};
    integer         kk{0};
    for ( integer i{0}; i < ns; ++i ) {
      s += ds;
      eval( s, p2[0], p2[1] );

      real_type const nx = p0[1]-p2[1];
      real_type const ny = p2[0]-p0[0];
      real_type const xm = (p0[0]+p2[0])/2;
      real_type const ym = (p0[1]+p2[1])/2;

      ++kk;
      Poly[kk][0] = w*(xm - nx * tg);
      Poly[kk][1] = w*(ym - ny * tg);
      Poly[kk][2] = w;

      ++kk;
      Poly[kk][0] = p2[0];
      Poly[kk][1] = p2[1];
      Poly[kk][2] = 1;

      knots[kk+1] = i+1;
      knots[kk+2] = i+1;

      p0[0] = p2[0];
      p0[1] = p2[1];

    }
    knots[kk+3] = ns;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::len_tolerance( real_type const tol ) const {
    real_type const absk { abs(m_k) };
    real_type const tmp  { absk*tol };
    if ( tmp > 0 ) {
      real_type const dtheta { 2*(Utils::m_pi-acos(tmp-1)) };
      return dtheta/absk;
    }
    return m_L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  CircleArc::info() const
  { return fmt::format( "CircleArc\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `CircleArc` object
  //!
  //!  \param stream the output stream
  //!  \param bi      an instance of `CircleArc` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, CircleArc const & bi ) {
    fmt::print( stream,
      "x     = {} : {}\n"
      "y     = {} : {}\n"
      "theta = {} : {}\n"
      "k     = {}\n"
      "L     = {}\n",
      bi.m_x0,     bi.x_end(),
      bi.m_y0,     bi.y_end(),
      bi.m_theta0, bi.theta_end(),
      bi.m_k, bi.m_L
    );
    return stream;
  }

}

// EOF: Circle.cc
