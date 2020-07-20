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

#include "ClothoidList.hh"
#include "Biarc.hh"
#include "BiarcList.hh"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <limits>
#include <algorithm>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace G2lib {

  using std::numeric_limits;
  using std::lower_bound;
  using std::vector;
  using std::swap;
  using std::abs;

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/

  ClothoidList::ClothoidList( LineSegment const & LS )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( LS );
  }

  ClothoidList::ClothoidList( CircleArc const & C )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( C );
  }

  ClothoidList::ClothoidList( Biarc const & C )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( C.getC0() );
    this->push_back( C.getC1() );
  }

  ClothoidList::ClothoidList( BiarcList const & c )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( c );
  }

  ClothoidList::ClothoidList( ClothoidCurve const & c )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( c );
  }

  ClothoidList::ClothoidList( PolyLine const & pl )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    this->push_back( pl );
  }

  ClothoidList::ClothoidList( BaseCurve const & C )
  : BaseCurve(G2LIB_CLOTHOID_LIST)
  , curve_is_closed(false)
  , aabb_done(false)
  {
    this->resetLastInterval();
    this->init();
    switch ( C.type() ) {
    case G2LIB_LINE:
      push_back( *static_cast<LineSegment const *>(&C) );
      break;
    case G2LIB_CIRCLE:
      push_back( *static_cast<CircleArc const *>(&C) );
      break;
    case G2LIB_CLOTHOID:
      push_back( *static_cast<ClothoidCurve const *>(&C) );
      break;
    case G2LIB_BIARC:
      push_back( *static_cast<Biarc const *>(&C) );
      break;
    case G2LIB_BIARC_LIST:
      push_back( *static_cast<BiarcList const *>(&C) );
      break;
    case G2LIB_CLOTHOID_LIST:
      copy( *static_cast<ClothoidList const *>(&C) );
      break;
    case G2LIB_POLYLINE:
      push_back( *static_cast<PolyLine const *>(&C) );
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::wrap_in_range( real_type & s ) const {
    real_type a = s0.front();
    real_type b = s0.back();
    real_type L = b-a;
    s -= a;
    s  = fmod( s, L );
    if ( s < 0 ) s += L;
    // while ( s < 0 ) s += L;
    // while ( s > L ) s -= L;
    s += a;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::init() {
    this->s0.clear();
    this->clotoidList.clear();
    this->resetLastInterval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::copy( ClothoidList const & L ) {
    clotoidList.clear();
    clotoidList.reserve(L.clotoidList.size());
    std::copy( L.clotoidList.begin(),
               L.clotoidList.end(),
               back_inserter(clotoidList) );
    s0.clear();
    s0.reserve(L.s0.size());
    std::copy( L.s0.begin(),
               L.s0.end(),
               back_inserter(s0) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reserve( int_type n ) {
    s0.reserve(size_t(n+1));
    clotoidList.reserve(size_t(n));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( LineSegment const & LS ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0);
      s0.push_back(LS.length());
    } else {
      s0.push_back(s0.back()+LS.length());
    }
    clotoidList.push_back(ClothoidCurve(LS));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( CircleArc const & C ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0);
      s0.push_back(C.length());
    } else {
      s0.push_back(s0.back()+C.length());
    }
    clotoidList.push_back(ClothoidCurve(C));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( Biarc const & c ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0);
      s0.push_back(c.length());
    } else {
      s0.push_back(s0.back()+c.getC0().length());
      s0.push_back(s0.back()+c.getC1().length());
    }
    clotoidList.push_back(ClothoidCurve(c.getC0()));
    clotoidList.push_back(ClothoidCurve(c.getC1()));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( ClothoidCurve const & c ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0);
      s0.push_back(c.length());
    } else {
      s0.push_back(s0.back()+c.length());
    }
    clotoidList.push_back(c);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( BiarcList const & c ) {
    s0.reserve( s0.size() + c.biarcList.size() + 1 );
    clotoidList.reserve( clotoidList.size() + 2*c.biarcList.size() );

    if ( s0.empty() ) s0.push_back(0);

    vector<Biarc>::const_iterator ip = c.biarcList.begin();
    for (; ip != c.biarcList.end(); ++ip ) {
      s0.push_back(s0.back()+ip->length());
      Biarc const & b = *ip;
      clotoidList.push_back(ClothoidCurve(b.getC0()));
      clotoidList.push_back(ClothoidCurve(b.getC1()));
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( PolyLine const & c ) {
    s0.reserve( s0.size() + c.polylineList.size() + 1 );
    clotoidList.reserve( clotoidList.size() + c.polylineList.size() );

    if ( s0.empty() ) s0.push_back(0);

    vector<LineSegment>::const_iterator ip = c.polylineList.begin();
    for (; ip != c.polylineList.end(); ++ip ) {
      s0.push_back(s0.back()+ip->length());
      clotoidList.push_back(ClothoidCurve(*ip));
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type kappa0,
    real_type dkappa,
    real_type L
  ) {
    G2LIB_ASSERT(
      !clotoidList.empty(),
      "ClothoidList::push_back_G1(...) empty list!"
    )
    ClothoidCurve c;
    real_type x0     = clotoidList.back().xEnd();
    real_type y0     = clotoidList.back().yEnd();
    real_type theta0 = clotoidList.back().thetaEnd();
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type kappa0,
    real_type dkappa,
    real_type L
  ) {
    ClothoidCurve c;
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    push_back( c );
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    G2LIB_ASSERT(
      !clotoidList.empty(),
      "ClothoidList::push_back_G1(...) empty list!"
    )
    ClothoidCurve c;
    real_type x0     = clotoidList.back().xEnd();
    real_type y0     = clotoidList.back().yEnd();
    real_type theta0 = clotoidList.back().thetaEnd();
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    ClothoidCurve c;
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1(
    int_type        n,
    real_type const x[],
    real_type const y[]
  ) {
    init();
    reserve( n-1 );
    ClothoidCurve c;

    G2LIB_ASSERT(
      n > 1,
      "ClothoidList::build_G1, at least 2 points are necessary"
    )

    if ( n == 2 ) {

      real_type theta = atan2( y[1] - y[0], x[1] - x[0] );
      c.build_G1( x[0], y[0], theta, x[1], y[1], theta );
      push_back(c);

    } else {

      Biarc b;
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10;
      real_type thetaC(0);
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" )
        thetaC = b.thetaMiddle();
      }
      ok = b.build_3P( x[0], y[0], x[1], y[1], x[2], y[2] );
      G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" )
      real_type theta0 = ciclic ? thetaC : b.thetaBegin();
      real_type theta1 = b.thetaMiddle();
      c.build_G1( x[0], y[0], theta0, x[1], y[1], theta1 );
      push_back(c);
      for ( int_type k = 2; k < n-1; ++k ) {
        theta0 = theta1;
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" )
        theta1 = b.thetaMiddle();
        c.build_G1( x[k-1], y[k-1], theta0, x[k], y[k], theta1 );
        push_back(c);
      }
      theta0 = theta1;
      theta1 = ciclic ? thetaC : b.thetaEnd();
      c.build_G1( x[n-2], y[n-2], theta0, x[n-1], y[n-1], theta1 );
      push_back(c);

    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1(
    int_type        n,
    real_type const x[],
    real_type const y[],
    real_type const theta[]
  ) {

    G2LIB_ASSERT(
      n > 1,
      "ClothoidList::build_G1, at least 2 points are necessary"
    )

    init();
    reserve( n-1 );
    ClothoidCurve c;
    for ( int_type k = 1; k < n; ++k ) {
      c.build_G1( x[k-1], y[k-1], theta[k-1], x[k], y[k], theta[k] );
      push_back(c);
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build(
    real_type       x0,
    real_type       y0,
    real_type       theta0,
    int_type        n,
    real_type const s[],
    real_type const kappa[]
  ) {
    if ( n < 2 ) return false;
    init();
    real_type k  = kappa[0];
    real_type L  = s[1]-s[0];
    real_type dk = (kappa[1]-k)/L;

    push_back( x0, y0, theta0, k, dk, L );
    for ( int_type i = 2; i < n; ++i ) {
      k  = kappa[i-1];
      L  = s[i]-s[i-1];
      dk = (kappa[i]-k)/L;
      push_back( k, dk, L );
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::get( int_type idx ) const {
    G2LIB_ASSERT(
      !clotoidList.empty(),
      "ClothoidList::get( " << idx << " ) empty list"
    )
    G2LIB_ASSERT(
      idx >= 0 && idx < int_type(clotoidList.size()),
      "ClothoidList::get( " << idx <<
      " ) bad index, must be in [0," << clotoidList.size()-1 << "]"
    )
    return clotoidList[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::getAtS( real_type s ) const
  { return get(findAtS(s)); }

  /*\
   |   _                  _   _
   |  | | ___ _ __   __ _| |_| |__
   |  | |/ _ \ '_ \ / _` | __| '_ \
   |  | |  __/ | | | (_| | |_| | | |
   |  |_|\___|_| |_|\__, |\__|_| |_|
   |                |___/
  \*/

  real_type
  ClothoidList::length() const {
    return s0.back() - s0.front();
  }

  real_type
  ClothoidList::length_ISO( real_type offs ) const {
    real_type L = 0;
    vector<ClothoidCurve>::const_iterator is = clotoidList.begin();
    for (; is != clotoidList.end(); ++is ) L += is->length_ISO( offs );
    return L;
  }

  real_type
  ClothoidList::segment_length( int_type nseg ) const {
    ClothoidCurve const & c = get( nseg );
    return c.length();
  }

  real_type
  ClothoidList::segment_length_ISO( int_type nseg, real_type offs ) const {
    ClothoidCurve const & c = get( nseg );
    return c.length_ISO( offs );
  }

  /*\
   |  _    _   _____    _                _
   | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
   | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
   | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
   |                               |___/
  \*/

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidList::bbTriangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size
  ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    for ( int_type ipos = 0; ic != clotoidList.end(); ++ic, ++ipos )
      ic->bbTriangles_ISO( offs, tvec, max_angle, max_size, ipos );
  }

  /*\
   |   _     _
   |  | |__ | |__   _____  __
   |  | '_ \| '_ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |_.__/|_.__/ \___/_/\_\
  \*/

  void
  ClothoidList::bbox_ISO(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    vector<Triangle2D> tvec;
    bbTriangles_ISO( offs, tvec, m_pi/18, 1e100 );
    xmin = ymin = numeric_limits<real_type>::infinity();
    xmax = ymax = -xmin;
    vector<Triangle2D>::const_iterator it;
    for ( it = tvec.begin(); it != tvec.end(); ++it ) {
      // - - - - - - - - - - - - - - - - - - - -
      if      ( it->x1() < xmin ) xmin = it->x1();
      else if ( it->x1() > xmax ) xmax = it->x1();
      if      ( it->x2() < xmin ) xmin = it->x2();
      else if ( it->x2() > xmax ) xmax = it->x2();
      if      ( it->x3() < xmin ) xmin = it->x3();
      else if ( it->x3() > xmax ) xmax = it->x3();
      // - - - - - - - - - - - - - - - - - - - -
      if      ( it->y1() < ymin ) ymin = it->y1();
      else if ( it->y1() > ymax ) ymax = it->y1();
      if      ( it->y2() < ymin ) ymin = it->y2();
      else if ( it->y2() > ymax ) ymax = it->y2();
      if      ( it->y3() < ymin ) ymin = it->y3();
      else if ( it->y3() > ymax ) ymax = it->y3();
    }
  }

  /*\
   |  _   _          _
   | | |_| |__   ___| |_ __ _
   | | __| '_ \ / _ \ __/ _` |
   | | |_| | | |  __/ || (_| |
   |  \__|_| |_|\___|\__\__,_|
  \*/

  real_type
  ClothoidList::theta( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.theta( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_D( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_D( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_DD( real_type s ) const  {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_DD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_DDD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_DDD( s - s0[idx] );
  }

  /*\
   |  _____                   _   _   _
   | |_   _|   __ _ _ __   __| | | \ | |
   |   | |    / _` | '_ \ / _` | |  \| |
   |   | |   | (_| | | | | (_| | | |\  |
   |   |_|    \__,_|_| |_|\__,_| |_| \_|
  \*/

  real_type
  ClothoidList::tx( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tx( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.ty( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_D( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_D( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_D( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_D( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_DD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_DD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_DD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_DD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_DDD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_DDD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_DDD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_DDD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg(
    real_type   s,
    real_type & tg_x,
    real_type & tg_y
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tg( s - s0[idx], tg_x, tg_y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_D(
    real_type   s,
    real_type & tg_x_D,
    real_type & tg_y_D
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_D( s - s0[idx], tg_x_D, tg_y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_DD(
    real_type   s,
    real_type & tg_x_DD,
    real_type & tg_y_DD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_DD( s - s0[idx], tg_x_DD, tg_y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_DDD(
    real_type   s,
    real_type & tg_x_DDD,
    real_type & tg_y_DDD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_DDD( s - s0[idx], tg_x_DDD, tg_y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::evaluate(
    real_type   s,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    c.evaluate( s - s0[idx], th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::evaluate_ISO(
    real_type   s,
    real_type   offs,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    c.evaluate_ISO( s - s0[idx], offs, th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y( real_type s ) const  {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_D( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_D( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_D( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_D( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_DD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_DD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_DD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_DD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_DDD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_DDD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_DDD( real_type s ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_DDD( s - s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval( s - s0[idx], x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_D( s - s0[idx], x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_DD( s - s0[idx], x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_DDD( s - s0[idx], x_DDD, y_DDD );
  }

  /*\
   |         __  __          _
   |   ___  / _|/ _|___  ___| |_
   |  / _ \| |_| |_/ __|/ _ \ __|
   | | (_) |  _|  _\__ \  __/ |_
   |  \___/|_| |_| |___/\___|\__|
  \*/

  real_type
  ClothoidList::X_ISO( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_D( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_D( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_D( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_D( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_DD( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_DD( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_DD( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_DD( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_DDD( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_DDD( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_DDD( real_type s, real_type offs ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_DDD( s - s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO( s - s0[idx], offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO_D( s - s0[idx], offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO_DD( s - s0[idx], offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( this->curve_is_closed ) this->wrap_in_range( s );
    int_type idx = findAtS( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO_DDD( s - s0[idx], offs, x_DDD, y_DDD );
  }

  /*\
   |  _                        __
   | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
   | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
   | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
   |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
  \*/

  void
  ClothoidList::translate( real_type tx, real_type ty ) {
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) ic->translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::rotate( real_type angle, real_type cx, real_type cy ) {
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) ic->rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::scale( real_type sfactor ) {
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    real_type newx0 = ic->xBegin();
    real_type newy0 = ic->yBegin();
    s0[0] = 0;
    for ( size_t k=0; ic != clotoidList.end(); ++ic, ++k ) {
      ic->scale( sfactor );
      ic->changeOrigin( newx0, newy0 );
      newx0 = ic->xEnd();
      newy0 = ic->yEnd();
      s0[k+1] = s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reverse() {
    std::reverse( clotoidList.begin(), clotoidList.end() );
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    ic->reverse();
    real_type newx0 = ic->xEnd();
    real_type newy0 = ic->yEnd();
    s0[0] = 0;
    s0[1] = ic->length();
    size_t k = 1;
    for ( ++ic; ic != clotoidList.end(); ++ic, ++k ) {
      ic->reverse();
      ic->changeOrigin( newx0, newy0 );
      newx0   = ic->xEnd();
      newy0   = ic->yEnd();
      s0[k+1] = s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::changeOrigin( real_type newx0, real_type newy0 ) {
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) {
      ic->changeOrigin( newx0, newy0 );
      newx0 = ic->xEnd();
      newy0 = ic->yEnd();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::trim( real_type s_begin, real_type s_end ) {
    G2LIB_ASSERT(
      s_begin >= s0.front() && s_end <= s0.back() && s_end > s_begin,
      "ClothoidList::trim( s_begin=" << s_begin << ", s_end=" << s_end <<
      ") bad range, must be in [ " << s0.front() << ", " << s0.back() << " ]"
    )

    size_t i_begin = size_t(findAtS( s_begin ));
    size_t i_end   = size_t(findAtS( s_end ));
    if ( i_begin == i_end ) {
      clotoidList[i_begin].trim( s_begin-s0[i_begin], s_end-s0[i_begin] );
    } else {
      clotoidList[i_begin].trim( s_begin-s0[i_begin], s0[i_begin+1]-s0[i_begin] );
      clotoidList[i_end].trim( 0, s_end-s0[i_end] );
    }
    clotoidList.erase( clotoidList.begin()+i_end+1, clotoidList.end() );
    clotoidList.erase( clotoidList.begin(), clotoidList.begin()+i_begin );
    if ( clotoidList.back().L <= machepsi100 ) clotoidList.pop_back();
    vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    s0.resize( clotoidList.size() + 1 );
    s0[0] = 0;
    size_t k = 0;
    for (; ic != clotoidList.end(); ++ic, ++k )
      s0[k+1] = s0[k] + ic->length();
    this->resetLastInterval();
  }

  /*\
   |     _        _    ____  ____  _
   |    / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |   / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |  / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   | /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  void
  ClothoidList::build_AABBtree_ISO(
    real_type offs,
    real_type max_angle,
    real_type max_size
  ) const {

    if ( aabb_done &&
         isZero( offs-aabb_offs ) &&
         isZero( max_angle-aabb_max_angle ) &&
         isZero( max_size-aabb_max_size ) ) return;

    #ifdef G2LIB_USE_CXX11
    vector<shared_ptr<BBox const> > bboxes;
    #else
    vector<BBox const *> bboxes;
    #endif

    bbTriangles_ISO( offs, aabb_tri, max_angle, max_size );
    bboxes.reserve(aabb_tri.size());
    vector<Triangle2D>::const_iterator it;
    int_type ipos = 0;
    for ( it = aabb_tri.begin(); it != aabb_tri.end(); ++it, ++ipos ) {
      real_type xmin, ymin, xmax, ymax;
      it->bbox( xmin, ymin, xmax, ymax );
      #ifdef G2LIB_USE_CXX11
      bboxes.push_back( make_shared<BBox const>(
        xmin, ymin, xmax, ymax, G2LIB_CLOTHOID, ipos
      ) );
      #else
      bboxes.push_back(
        new BBox( xmin, ymin, xmax, ymax, G2LIB_CLOTHOID, ipos )
      );
      #endif
    }
    aabb_tree.build(bboxes);
    aabb_done      = true;
    aabb_offs      = offs;
    aabb_max_angle = max_angle;
    aabb_max_size  = max_size;
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
  ClothoidList::collision( ClothoidList const & C ) const {
    this->build_AABBtree_ISO( 0 );
    C.build_AABBtree_ISO( 0 );
    T2D_collision_list_ISO fun( this, 0, &C, 0 );
    return aabb_tree.collision( C.aabb_tree, fun, false );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::collision_ISO(
    real_type            offs,
    ClothoidList const & C,
    real_type            offs_C
  ) const {
    this->build_AABBtree_ISO( offs );
    C.build_AABBtree_ISO( offs_C );
    T2D_collision_list_ISO fun( this, offs, &C, offs_C );
    return aabb_tree.collision( C.aabb_tree, fun, false );
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
   |
  \*/

  void
  ClothoidList::intersect_ISO(
    real_type            offs,
    ClothoidList const & CL,
    real_type            offs_CL,
    IntersectList      & ilist,
    bool                 swap_s_vals
  ) const {
    if ( intersect_with_AABBtree ) {
      this->build_AABBtree_ISO( offs );
      CL.build_AABBtree_ISO( offs_CL );
      AABBtree::VecPairPtrBBox iList;
      aabb_tree.intersect( CL.aabb_tree, iList );

      AABBtree::VecPairPtrBBox::const_iterator ip;
      for ( ip = iList.begin(); ip != iList.end(); ++ip ) {
        size_t ipos1 = size_t(ip->first->Ipos());
        size_t ipos2 = size_t(ip->second->Ipos());

        Triangle2D const & T1 = aabb_tri[ipos1];
        Triangle2D const & T2 = CL.aabb_tri[ipos2];

        ClothoidCurve const & C1 = clotoidList[T1.Icurve()];
        ClothoidCurve const & C2 = CL.clotoidList[T2.Icurve()];

        real_type ss1, ss2;
        bool converged = C1.aabb_intersect_ISO( T1, offs, &C2, T2, offs_CL, ss1, ss2 );

        if ( converged ) {
          ss1 += s0[T1.Icurve()];
          ss2 += CL.s0[T2.Icurve()];
          if ( swap_s_vals ) swap( ss1, ss2 );
          ilist.push_back( Ipair( ss1, ss2 ) );
        }
      }
    } else {
      bbTriangles_ISO( offs, aabb_tri, m_pi/18, 1e100 );
      CL.bbTriangles_ISO( offs_CL, CL.aabb_tri, m_pi/18, 1e100 );
      for ( vector<Triangle2D>::const_iterator i1 = aabb_tri.begin();
            i1 != aabb_tri.end(); ++i1 ) {
        for ( vector<Triangle2D>::const_iterator i2 = CL.aabb_tri.begin();
              i2 != CL.aabb_tri.end(); ++i2 ) {
          Triangle2D const & T1 = *i1;
          Triangle2D const & T2 = *i2;

          ClothoidCurve const & C1 = clotoidList[T1.Icurve()];
          ClothoidCurve const & C2 = CL.clotoidList[T2.Icurve()];

          real_type ss1, ss2;
          bool converged = C1.aabb_intersect_ISO( T1, offs, &C2, T2, offs_CL, ss1, ss2 );

          if ( converged ) {
            ss1 += s0[T1.Icurve()];
            ss2 += CL.s0[T2.Icurve()];
            if ( swap_s_vals ) swap( ss1, ss2 );
            ilist.push_back( Ipair( ss1, ss2 ) );
          }
        }
      }

    }
  }

  /*\
   |      _ _     _
   |   __| (_)___| |_ __ _ _ __   ___ ___
   |  / _` | / __| __/ _` | '_ \ / __/ _ \
   | | (_| | \__ \ || (_| | | | | (_|  __/
   |  \__,_|_|___/\__\__,_|_| |_|\___\___|
  \*/

  int_type
  ClothoidList::closestPoint_ISO(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & DST
  ) const {

    this->build_AABBtree_ISO( offs );

    AABBtree::VecPtrBBox candidateList;
    aabb_tree.min_distance( qx, qy, candidateList );
    AABBtree::VecPtrBBox::const_iterator ic;
    G2LIB_ASSERT(
      candidateList.size() > 0,
      "ClothoidList::closestPoint no candidate"
    )
    int_type icurve = 0;
    DST = numeric_limits<real_type>::infinity();
    for ( ic = candidateList.begin(); ic != candidateList.end(); ++ic ) {
      size_t ipos = size_t((*ic)->Ipos());
      Triangle2D const & T = aabb_tri[ipos];
      real_type dst = T.distMin( qx, qy );
      if ( dst < DST ) {
        // refine distance
        real_type xx, yy, ss;
        clotoidList[T.Icurve()].closestPoint_internal_ISO(
          T.S0(), T.S1(), qx, qy, offs, xx, yy, ss, dst
        );
        if ( dst < DST ) {
          DST    = dst;
          s      = ss + s0[T.Icurve()];
          x      = xx;
          y      = yy;
          icurve = T.Icurve();
        }
      }
    }

    real_type nx, ny;
    clotoidList[icurve].nor_ISO( s - s0[icurve], nx, ny );
    t = (qx-x) * nx + (qy-y) * ny - offs;
    real_type err = abs( abs(t) - DST );
    if ( err > DST*machepsi1000 ) return -1;
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  ClothoidList::closestPoint_ISO(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    return closestPoint_ISO( qx, qy, 0, x, y, s, t, dst );
  }

  /*\
   |      _ _     _
   |   __| (_)___| |_ __ _ _ __   ___ ___
   |  / _` | / __| __/ _` | '_ \ / __/ _ \
   | | (_| | \__ \ || (_| | | | | (_|  __/
   |  \__,_|_|___/\__\__,_|_| |_|\___\___|
  \*/

  int_type
  ClothoidList::closestSegment( real_type qx, real_type qy ) const {
    this->build_AABBtree_ISO( 0 );

    AABBtree::VecPtrBBox candidateList;
    aabb_tree.min_distance( qx, qy, candidateList );
    AABBtree::VecPtrBBox::const_iterator ic;
    G2LIB_ASSERT(
      candidateList.size() > 0,
      "ClothoidList::closestSegment no candidate"
    )
    int_type icurve = 0;
    real_type DST = numeric_limits<real_type>::infinity();
    for ( ic = candidateList.begin(); ic != candidateList.end(); ++ic ) {
      size_t ipos = size_t((*ic)->Ipos());
      Triangle2D const & T = aabb_tri[ipos];
      real_type dst = T.distMin( qx, qy );
      if ( dst < DST ) {
        // refine distance
        real_type xx, yy, ss;
        clotoidList[T.Icurve()].closestPoint_internal_ISO(
          T.S0(), T.S1(), qx, qy, 0, xx, yy, ss, dst
        );
        if ( dst < DST ) {
          DST    = dst;
          icurve = T.Icurve();
        }
      }
    }
    return icurve;
  }

  int_type
  ClothoidList::closestPointInRange_ISO(
    real_type   qx,
    real_type   qy,
    int_type    icurve_begin,
    int_type    icurve_end,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst,
    int_type  & icurve
  ) const {
    G2LIB_ASSERT( !this->clotoidList.empty(), "ClothoidList::closestPointInRange_ISO, empty list" )
    int_type nsegs = this->numSegment();
    if ( nsegs == 1 ) {
      icurve = 0;
      int_type res = this->clotoidList.front().closestPoint_ISO( qx, qy, x, y, s, t, dst );
      s += s0[0];
      return res;
    }

    int_type ib = icurve_begin % nsegs; // to avoid infinite loop in case of bad input
    int_type ie = icurve_end   % nsegs; // to avoid infinite loop in case of bad input
    if ( ib < 0 ) ib += nsegs;
    if ( ie < 0 ) ie += nsegs;
    G2LIB_ASSERT(
      ib >= 0 && ie >= 0,
      "ClothoidList::closestPointInRange_ISO, ib = " << ib  << " ie = " << ie
    )

    icurve = ib;
    int_type res = this->clotoidList[icurve].closestPoint_ISO( qx, qy, x, y, s, t, dst );
    s += this->s0[icurve];

    if ( ib == ie ) return res; // only one segment to check

    int_type iseg = ib;
    do {
      if ( ++iseg >= nsegs ) iseg -= nsegs; // next segment
      real_type C_x, C_y, C_s, C_t, C_dst;
      int_type C_res = this->clotoidList[iseg].closestPoint_ISO( qx, qy, C_x, C_y, C_s, C_t, C_dst );
      if ( C_dst < dst ) {
        dst    = C_dst;
        x      = C_x;
        y      = C_y;
        s      = C_s+s0[iseg];
        t      = C_t;
        icurve = iseg;
        res    = C_res;
      }
    } while ( iseg != ie );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getSK( real_type s[], real_type kappa[] ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type  k  = 0;
    real_type ss = 0;
    while ( ic != clotoidList.end() ) {
      s[k]     = ss;
      kappa[k] = ic->kappaBegin();
      ss      += ic->length();
      ++k;
      ++ic;
    }
    --ic; // last element
    s[k]     = ss;
    kappa[k] = ic->kappaEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getSTK(
    real_type s[],
    real_type theta[],
    real_type kappa[]
  ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type  k  = 0;
    real_type ss = 0;
    while ( ic != clotoidList.end() ) {
      s[k]     = ss;
      theta[k] = ic->thetaBegin();
      kappa[k] = ic->kappaBegin();
      ss      += ic->length();
      ++k;
      ++ic;
    }
    --ic; // last element
    s[k]     = ss;
    theta[k] = ic->thetaEnd();
    kappa[k] = ic->kappaEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getXY( real_type x[], real_type y[] ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k  = 0;
    while ( ic != clotoidList.end() ) {
      x[k] = ic->xBegin();
      y[k] = ic->yBegin();
      ++k; ++ic;
    }
    --ic;
    x[k] = ic->xEnd();
    y[k] = ic->yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaTheta( real_type deltaTheta[] ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k = 0;
    for ( ++ic; ic != clotoidList.end(); ++ic, ++k ) {
      real_type tmp = ic->thetaBegin()-ic[-1].thetaEnd();
      if      ( tmp >  m_pi ) tmp -= m_2pi;
      else if ( tmp < -m_pi ) tmp += m_2pi;
      deltaTheta[k] = tmp;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaKappa( real_type deltaKappa[] ) const {
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k = 0;
    for ( ++ic; ic != clotoidList.end(); ++ic, ++k  )
      deltaKappa[k] = ic->kappaBegin()-ic[-1].kappaEnd();

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  ClothoidList::findST1(
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::findST, empty list" )
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    vector<real_type>::const_iterator     is = s0.begin();

    s = t = 0;
    int_type  ipos = 0;
    int_type  iseg = 0;
    real_type S, T;
    bool ok = ic->findST_ISO( x, y, S, T );
    if ( ok ) {
      s    = *is + S;
      t    = T;
      iseg = 0;
    }

    for ( ++ic, ++is, ++ipos;
          ic != clotoidList.end();
          ++ic, ++is, ++ipos ) {
      bool ok1 = ic->findST_ISO( x, y, S, T );
      if ( ok && ok1 ) ok1 = abs(T) < abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = *is + S;
        t    = T;
        iseg = ipos;
      }
    }

    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  ClothoidList::findST1(
    int_type    ibegin,
    int_type    iend,
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    G2LIB_ASSERT(
      !clotoidList.empty(),
      "ClothoidList::findST, empty list"
    )
    G2LIB_ASSERT(
      ibegin >= 0 && ibegin <= iend && iend < int_type(clotoidList.size()),
      "ClothoidList::findST( ibegin=" << ibegin << ", iend = " << iend <<
      " , x, y, s, t ) bad range not in [0," << clotoidList.size()-1 << "]"
    )
    s = t = 0;
    int_type iseg = 0;
    bool ok = false;
    for ( int_type k = ibegin; k <= iend; ++k ) {
      ClothoidCurve const & ck = clotoidList[k];
      real_type S, T;
      bool ok1 = ck.findST_ISO( x, y, S, T );
      if ( ok && ok1 ) ok1 = abs(T) < abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = s0[k] + S;
        t    = T;
        iseg = k;
      }
    }
    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_table( ostream_type & stream ) const {
    stream << "x\ty\ttheta0\tkappa0\tdkappa\tL\n";
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic )
      stream
        << ic->xBegin()     << '\t'
        << ic->yBegin()     << '\t'
        << ic->thetaBegin() << '\t'
        << ic->kappaBegin() << '\t'
        << ic->dkappa()     << '\t'
        << ic->length()     << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_ruby( ostream_type & stream ) const {
    stream << "data = {\n";
    vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic )
      stream
        << ic->xBegin()     << '\t'
        << ic->yBegin()     << '\t'
        << ic->thetaBegin() << '\t'
        << ic->kappaBegin() << '\t'
        << ic->dkappa()     << '\t'
        << ic->length()     << '\n';
    stream << "}\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, ClothoidList const & CL ) {
    vector<ClothoidCurve>::const_iterator ic = CL.clotoidList.begin();
    for (; ic != CL.clotoidList.end(); ++ic )
      stream << *ic << '\n';
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

// EOF: ClothoidList.cc
