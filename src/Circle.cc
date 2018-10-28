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

#include "Circle.hh"
#include "CubicRootsFlocke.hh"

#include <cmath>

namespace G2lib {

  using namespace std;

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  bool
  CircleArc::build_G1( real_type _x0,
                       real_type _y0,
                       real_type _theta0,
                       real_type _x1,
                       real_type _y1 ) {

    real_type dx = _x1 - _x0;
    real_type dy = _y1 - _y0;
    real_type d  = hypot( dx, dy );

    if ( d > 0 ) {
      real_type th = atan2( dy, dx ) - _theta0;
      x0     = _x0;
      y0     = _y0;
      theta0 = _theta0;
      k      = 2*sin(th)/d;
      L      = d/Sinc(th);
      return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::build_3P( real_type _x0,
                       real_type _y0,
                       real_type _x1,
                       real_type _y1,
                       real_type _x2,
                       real_type _y2 ) {

    real_type dxa   = _x1 - _x0;
    real_type dya   = _y1 - _y0;
    real_type dxb   = _x2 - _x1;
    real_type dyb   = _y2 - _y1;
    real_type La    = hypot(dya,dxa);
    real_type Lb    = hypot(dyb,dxb);
    real_type cosom = (dxa*dxb + dya*dyb)/(La*Lb);
    if      ( cosom >  1 ) cosom = 1;
    else if ( cosom < -1 ) cosom = -1;
    real_type omega = acos(cosom);

    real_type alpha = omega - atan2(Lb*sin(omega),La+Lb*cos(omega));
    real_type dxc   = _x2 - _x0;
    real_type dyc   = _y2 - _y0;
    real_type Lc    = hypot(dyc,dxc);
    real_type cosal = (dxa*dxc + dya*dyc)/(La*Lc);
    if      ( cosal >  1 ) cosal = 1;
    else if ( cosal < -1 ) cosal = -1;
    alpha += acos(cosal);

    if ( dxa*dyb > dya*dxb ) alpha = -alpha;
    real_type _theta0 = atan2( dyc, dxc ) + alpha;
    return build_G1( _x0, _y0, _theta0, _x2, _y2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::thetaMinMax( real_type & thMin, real_type & thMax ) const  {
    thMin = theta0;
    thMax = theta0 + L * k;
    if ( thMax < thMin ) std::swap( thMin, thMax );
    return thMax-thMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::X( real_type s ) const {
    real_type sk = (s*k)/2;
    return x0+s*Sinc(sk)*cos(theta0+sk);
  }

  real_type
  CircleArc::X_D( real_type s ) const {
    return cos(theta0+s*k);
  }

  real_type
  CircleArc::X_DD( real_type s ) const {
    return -k*sin(theta0+s*k);
  }

  real_type
  CircleArc::X_DDD( real_type s ) const {
    return -(k*k)*cos(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y( real_type s ) const {
    real_type sk = (s*k)/2;
    return x0+s*Sinc(sk)*sin(theta0+sk);
  }

  real_type
  CircleArc::Y_D( real_type s ) const {
    return sin(theta0+s*k);
  }

  real_type
  CircleArc::Y_DD( real_type s ) const {
    return k*cos(theta0+s*k);
  }

  real_type
  CircleArc::Y_DDD( real_type s ) const {
    return -(k*k)*sin(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::tg( real_type s, real_type t[2] ) const {
    real_type th = theta(s);
    t[0] = cos(th);
    t[1] = sin(th);
  }

  void
  CircleArc::tg_D( real_type s, real_type t_D[2] ) const {
    real_type th = theta(s);
    t_D[0] = -sin(th)*k;
    t_D[1] = cos(th)*k;
  }

  void
  CircleArc::tg_DD( real_type s, real_type t_DD[2] ) const {
    real_type th = theta(s);
    real_type k2 = k*k;
    t_DD[0] = -cos(th)*k2;
    t_DD[1] = -sin(th)*k2;
  }

  void
  CircleArc::tg_DDD( real_type s, real_type t_DDD[2] ) const {
    real_type th = theta(s);
    real_type k3 = k*k*k;
    t_DDD[0] = sin(th)*k3;
    t_DDD[1] = -cos(th)*k3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::nor( real_type s, real_type n[2] ) const {
    real_type th = theta(s);
    n[0] = -sin(th);
    n[1] = cos(th);
  }

  void
  CircleArc::nor_D( real_type s, real_type n_D[2] ) const {
    real_type th = theta(s);
    n_D[0] = -cos(th)*k;
    n_D[1] = -sin(th)*k;
  }

  void
  CircleArc::nor_DD( real_type s, real_type n_DD[2] ) const {
    real_type th = theta(s);
    real_type k2 = k*k;
    n_DD[0] = -cos(th)*k2;
    n_DD[1] = -sin(th)*k2;
  }

  void
  CircleArc::nor_DDD( real_type s, real_type n_DDD[2] ) const {
    real_type th = theta(s);
    real_type k3 = k*k*k;
    n_DDD[0] = sin(th)*k3;
    n_DDD[1] = -cos(th)*k3;
  }

  /*\
   |  _____                   _   _   _
   | |_   _|   __ _ _ __   __| | | \ | |
   |   | |    / _` | '_ \ / _` | |  \| |
   |   | |   | (_| | | | | (_| | | |\  |
   |   |_|    \__,_|_| |_|\__,_| |_| \_|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval( real_type   s,
                   real_type & x,
                   real_type & y ) const {
    real_type sk  = (s*k)/2;
    real_type LS  = s*Sinc(sk);
    real_type arg = theta0+sk;
    x = x0+LS*cos(arg);
    y = y0+LS*sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_D( real_type   s,
                     real_type & x_D,
                     real_type & y_D ) const {
    real_type arg = theta0+s*k;
    x_D = cos(arg);
    y_D = sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DD( real_type   s,
                      real_type & x_DD,
                      real_type & y_DD ) const {
    real_type arg = theta0+s*k;
    x_DD = -k*sin(arg);
    y_DD = k*cos(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DDD( real_type   s,
                       real_type & x_DDD,
                       real_type & y_DDD ) const {
    real_type arg = theta0+s*k;
    real_type k2  = k*k;
    x_DDD = -k2*cos(arg);
    y_DDD = -k2*sin(arg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::trim( real_type s_begin, real_type s_end ) {
    real_type x, y;
    eval( s_begin, x, y );
    theta0 += s_begin * k;
    L  = s_end - s_begin;
    x0 = x;
    y0 = y;
    c0 = cos(theta0);
    s0 = sin(theta0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::rotate( real_type angle, real_type cx, real_type cy ) {
    real_type dx  = x0 - cx;
    real_type dy  = y0 - cy;
    real_type C   = cos(angle);
    real_type S   = sin(angle);
    real_type ndx = C*dx - S*dy;
    real_type ndy = C*dy + S*dx;
    x0      = cx + ndx;
    y0      = cy + ndy;
    theta0 += angle;
    c0      = cos(theta0);
    s0      = sin(theta0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::scale( real_type s ) {
    k /= s;
    L *= s;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::reverse() {
    theta0 = theta0 + m_pi;
    if ( theta0 > m_pi ) theta0 -= 2*m_pi;
    c0 = cos(theta0);
    s0 = sin(theta0);
    k  = -k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::center( real_type & cx, real_type & cy ) const {
    real_type nx = -sin(theta0);
    real_type ny = cos(theta0);
    cx = x0 + nx/k;
    cy = y0 + ny/k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::findST( real_type   x,
                     real_type   y,
                     real_type & s,
                     real_type & t ) const {
    real_type X, Y, n[2];
    s = projectPointOnCircle( x0, y0, cos(theta0), sin(theta0), k, L, x, y );
    eval( s, X, Y );
    nor( s, n );
    t = n[0]*(x-X) + n[1]*(y-Y);
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::changeCurvilinearOrigin( real_type new_s0, real_type newL ) {
    real_type new_x0, new_y0;
    eval( new_s0,  new_x0, new_y0 );
    x0      = new_x0;
    y0      = new_y0;
    theta0 += k*new_s0;
    L       = newL;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! get the bounding box triangle (if angle variation less that pi/3)
  bool
  CircleArc::bbTriangle( real_type & xx0, real_type & yy0,
                         real_type & xx1, real_type & yy1,
                         real_type & xx2, real_type & yy2 ) const {
    real_type dtheta = L * k;
    bool ok = std::abs(dtheta) <= m_pi/3;
    if ( ok ) {
      xx0 = x0; yy0 = y0;
      eval( L, xx2, yy2 );
      xx1 = (xx0+xx2)/2;
      yy1 = (yy0+yy2)/2;
      real_type nx = yy0-yy2;
      real_type ny = xx2-xx0;
      real_type tg = std::tan(dtheta/2)/2;
      xx1 -= nx * tg;
      yy1 -= ny * tg;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! get the bounding box triangle (if angle variation less that pi/3)
  bool
  CircleArc::bbTriangle( real_type   offs,
                         real_type & xx0, real_type & yy0,
                         real_type & xx1, real_type & yy1,
                         real_type & xx2, real_type & yy2 ) const {
    real_type dtheta = L * k;
    bool ok = std::abs(dtheta) <= m_pi/3;
    if ( ok ) {
      eval( 0, offs, xx0, yy0 );
      eval( L, offs, xx2, yy2 );
      xx1 = (xx0+xx2)/2;
      yy1 = (yy0+yy2)/2;
      real_type nx = yy0-yy2;
      real_type ny = xx2-xx0;
      real_type tg = std::tan(dtheta/2)/2;
      xx1 -= nx * tg;
      yy1 -= ny * tg;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbTriangles( std::vector<Triangle2D> & tvec,
                          real_type max_angle ) const {
    real_type dtheta = std::abs(L * k);
    int_type  n      = 1;
    if ( dtheta > max_angle ) {
      n       = int_type(std::ceil( dtheta/max_angle ));
      dtheta /= n ;
    }
    tvec.reserve( n );
    real_type xx0 = x0;
    real_type yy0 = y0;
    real_type ds  = L/n;
    real_type ss  = ds;
    for ( int_type k = 0; k < n; ++k, ss += ds ) {
      real_type xx2, yy2;
      eval( ss, xx2, yy2 );
      real_type xx1 = (xx0+xx2)/2;
      real_type yy1 = (yy0+yy2)/2;
      real_type nx = yy0-yy2;
      real_type ny = xx2-xx0;
      real_type tg = std::tan(dtheta/2)/2;
      xx1 -= nx * tg;
      yy1 -= ny * tg;
      tvec.push_back( Triangle2D( xx0, yy0, xx1, yy1, xx2, yy2) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbTriangle( real_type offs,
                         std::vector<Triangle2D> & tvec,
                         real_type max_angle ) const {
    real_type dtheta = std::abs(L * k);
    int_type  n      = 1;
    if ( dtheta > max_angle ) {
      n       = int_type(std::ceil( dtheta/max_angle ));
      dtheta /= n ;
    }
    tvec.reserve( n );
    real_type ds    = L/n;
    real_type ss    = ds;
    real_type scale = 1-k*offs;
    real_type xx0, yy0;
    eval( 0, offs, xx0, yy0 );
    for ( int_type k = 0; k < n; ++k, ss += ds ) {
      real_type xx2, yy2;
      eval( ss, offs, xx2, yy2 );
      real_type xx1 = (xx0+xx2)/2;
      real_type yy1 = (yy0+yy2)/2;
      real_type nx = yy0-yy2;
      real_type ny = xx2-xx0;
      real_type tg = scale * std::tan(dtheta/2)/2;
      xx1 -= nx * tg;
      yy1 -= ny * tg;
      tvec.push_back( Triangle2D( xx0, yy0, xx1, yy1, xx2, yy2) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbox( real_type & xmin,
                   real_type & ymin,
                   real_type & xmax,
                   real_type & ymax ) const {
    real_type xx0, yy0, xx1, yy1, xx2, yy2;
    bool ok = bbTriangle( xx0, yy0, xx1, yy1, xx2, yy2 );
    if ( ok ) {
      minmax3( xx0, xx1, xx2, xmin, xmax );
      minmax3( yy0, yy1, yy2, ymin, ymax );
    } else {
      real_type cx, cy, delta;
      center( cx, cy );
      delta = 1/std::abs(k);
      xmin = cx-delta;
      xmax = cx+delta;
      ymin = cy-delta;
      ymax = cy+delta;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::bbox( real_type   offs,
                   real_type & xmin,
                   real_type & ymin,
                   real_type & xmax,
                   real_type & ymax ) const {
    real_type xx0, yy0, xx1, yy1, xx2, yy2;
    bool ok = bbTriangle( offs, xx0, yy0, xx1, yy1, xx2, yy2 );
    if ( ok ) {
      minmax3( xx0, xx1, xx2, xmin, xmax );
      minmax3( yy0, yy1, yy2, ymin, ymax );
    } else {
      real_type cx, cy, delta;
      center( cx, cy );
      delta = std::abs(1/k-offs);
      xmin = cx-delta;
      xmax = cx+delta;
      ymin = cy-delta;
      ymax = cy+delta;
    }
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  bool
  CircleArc::collision( BaseCurve const & obj ) const {
    bool ok = false;
    switch ( obj.type() ) {
    case G2LIB_LINE:
      { // promote
        CircleArc C(*static_cast<LineSegment const*>(&obj)) ;
        ok = this->collision( C );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C(*static_cast<CircleArc const*>(&obj)) ;
        ok = this->collision( C );
      }
      break;
    case G2LIB_POLYLINE:
    case G2LIB_BIARC:
    case G2LIB_CLOTHOID:
    case G2LIB_CLOTHOID_LIST:
      break;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::collision( real_type         offs,
                        BaseCurve const & obj,
                        real_type         offs_obj ) const {
    bool ok = false;
    switch ( obj.type() ) {
    case G2LIB_LINE:
      { // promote to arc
        CircleArc C(*static_cast<LineSegment const*>(&obj)) ;
        ok = this->collision( offs, C, offs_obj );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C(*static_cast<CircleArc const*>(&obj)) ;
        ok = this->collision( offs, C, offs_obj );
      }
      break;
    case G2LIB_POLYLINE:
    case G2LIB_BIARC:
    case G2LIB_CLOTHOID:
    case G2LIB_CLOTHOID_LIST:
      break;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::intersect( BaseCurve const & obj,
                        IntersectList   & ilist ) const {
    switch ( obj.type() ) {
    case G2LIB_LINE:
      { // promote to arc
        CircleArc C(*static_cast<LineSegment const*>(&obj)) ;
        this->intersect( C, ilist );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C(*static_cast<CircleArc const*>(&obj)) ;
        this->intersect( C, ilist );
      }
      break;
    case G2LIB_POLYLINE:
    case G2LIB_BIARC:
    case G2LIB_CLOTHOID:
    case G2LIB_CLOTHOID_LIST:
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::intersect( real_type         offs,
                        BaseCurve const & obj,
                        real_type         offs_obj,
                        IntersectList   & ilist ) const {
    switch ( obj.type() ) {
    case G2LIB_LINE:
      { // promote to arc
        CircleArc C(*static_cast<LineSegment const*>(&obj)) ;
        this->intersect( offs, C, offs_obj, ilist );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C(*static_cast<CircleArc const*>(&obj)) ;
        this->intersect( offs, C, offs_obj, ilist );
      }
      break;
    case G2LIB_POLYLINE:
    case G2LIB_BIARC:
    case G2LIB_CLOTHOID:
    case G2LIB_CLOTHOID_LIST:
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  CircleArc::collision( CircleArc const & C ) const {
    real_type s1[2], s2[2];
    int_type ni = intersectCircleCircle( this->x0,
                                         this->y0,
                                         this->theta0,
                                         this->k,
                                         C.x0,
                                         C.y0,
                                         C.theta0,
                                         C.k,
                                         s1, s2 );
    real_type eps1 = machepsi100*L;
    real_type eps2 = machepsi100*C.L;
    for ( int_type i = 0; i < ni; ++i ) {
      if ( s1[i] >= -eps1 && s1[i] <= L+eps1 &&
           s2[i] >= -eps2 && s2[i] <= L+eps2 )
        return true;
    }
    return false;
  }

  bool
  CircleArc::collision( real_type         offs,
                        CircleArc const & C,
                        real_type         offs_C ) const {
    real_type s1[2], s2[2];
    real_type sc1 = 1-k*offs;
    real_type sc2 = 1-C.k*offs_C;
    int_type ni = intersectCircleCircle( this->X(0,offs),
                                         this->Y(0,offs),
                                         this->theta0,
                                         this->k/sc2,
                                         C.X(0,offs_C),
                                         C.Y(0,offs_C),
                                         C.theta0,
                                         C.k/sc2,
                                         s1, s2 );
    real_type eps1 = machepsi100*L;
    real_type eps2 = machepsi100*C.L;
    for ( int_type i = 0; i < ni; ++i ) {
      real_type ss1 = s1[i]/sc1;
      real_type ss2 = s2[i]/sc2;
      if ( ss1 >= -eps1 && ss1 <= L+eps1 &&
           ss2 >= -eps2 && ss2 <= L+eps2 )
        return true;
    }
    return false;
  }

  void
  CircleArc::intersect( CircleArc const & C,
                        IntersectList   & ilist ) const {
    real_type s1[2], s2[2];
    int_type ni = intersectCircleCircle( this->x0,
                                         this->y0,
                                         this->theta0,
                                         this->k,
                                         C.x0,
                                         C.y0,
                                         C.theta0,
                                         C.k,
                                         s1, s2 );
    real_type eps1 = machepsi100*L;
    real_type eps2 = machepsi100*C.L;
    for ( int_type i = 0; i < ni; ++i ) {
      if ( s1[i] >= -eps1 && s1[i] <= L+eps1 &&
           s2[i] >= -eps2 && s2[i] <= L+eps2 )
        ilist.push_back( Ipair(s1[i],s2[i]) );
    }
  }

  void
  CircleArc::intersect( real_type         offs,
                        CircleArc const & C,
                        real_type         offs_C,
                        IntersectList   & ilist ) const {
    real_type s1[2], s2[2];
    real_type sc1 = 1-k*offs;
    real_type sc2 = 1-C.k*offs_C;
    int_type ni = intersectCircleCircle( this->X(0,offs),
                                         this->Y(0,offs),
                                         this->theta0,
                                         this->k/sc2,
                                         C.X(0,offs_C),
                                         C.Y(0,offs_C),
                                         C.theta0,
                                         C.k/sc2,
                                         s1, s2 );
    real_type eps1 = machepsi100*L;
    real_type eps2 = machepsi100*C.L;
    for ( int_type i = 0; i < ni; ++i ) {
      real_type ss1 = s1[i]/sc1;
      real_type ss2 = s2[i]/sc2;
      if ( ss1 >= -eps1 && ss1 <= L+eps1 &&
           ss2 >= -eps2 && ss2 <= L+eps2 )
        ilist.push_back( Ipair(ss1,ss2) );
    }
  }

  /*\
   |                  _           _   _
   |  _ __  _ __ ___ (_) ___  ___| |_(_) ___  _ __
   | | '_ \| '__/ _ \| |/ _ \/ __| __| |/ _ \| '_ \
   | | |_) | | | (_) | |  __/ (__| |_| | (_) | | | |
   | | .__/|_|  \___// |\___|\___|\__|_|\___/|_| |_|
   | |_|           |__/
  \*/

  int_type
  CircleArc::projection( real_type   qx,
                         real_type   qy,
                         real_type & x,
                         real_type & y,
                         real_type & s ) const {
    real_type c0 = cos(theta0);
    real_type s0 = sin(theta0);
    s = projectPointOnCircle( x0, y0, c0, s0, k, L, qx, qy );
    if ( s < 0 || s > L ) {
      s = L;
      eval( s, x, y );
      // costruisco piano
      real_type nx = x-x0;
      real_type ny = y-y0;
      real_type dx = 2*qx-(x0+x);
      real_type dy = 2*qy-(y0+y);
      if ( nx*dx + ny*dy <= 0 ) {
        s = 0;
        x = x0;
        y = y0;
      }
      return 0;
    }
    eval( s, x, y );
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  CircleArc::projection( real_type   qx,
                         real_type   qy,
                         real_type   offs,
                         real_type & x,
                         real_type & y,
                         real_type & s ) const  {
    real_type c0  = cos(theta0);
    real_type s0  = sin(theta0);
    real_type xx0 = x0-offs*s0;
    real_type yy0 = y0+offs*c0;
    real_type ff = 1-k*offs;
    real_type LL = L*ff;
    s = projectPointOnCircle( xx0, yy0, c0, s0, k/ff, LL, qx, qy );
    if ( s < 0 || s > LL ) {
      s = L;
      eval( s, offs, x, y );
      // costruisco piano
      real_type nx = x-xx0;
      real_type ny = y-yy0;
      real_type dx = 2*qx-(xx0+x);
      real_type dy = 2*qy-(yy0+y);
      if ( nx*dx + ny*dy <= 0 ) {
        s = 0;
        x = xx0;
        y = yy0;
      }
      return 0;
    }
    eval( s, offs, x, y );
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::paramNURBS( int_type & n_knots,
                         int_type & n_pnts ) const {
    real_type dtheta = L*k;
    int_type  ns     = int_type(std::floor(3*std::abs(dtheta)/m_pi));
    if ( ns < 1 ) ns = 1;
    n_pnts  = 1+2*ns;
    n_knots = n_pnts+3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::toNURBS( real_type knots[],
                      real_type Poly[][3] ) const {

    real_type dtheta = L*k;
    int_type  ns     = int_type(std::floor(3*std::abs(dtheta)/m_pi));
    if ( ns < 1 ) ns = 1;

    real_type th = dtheta/(2*ns);
    real_type w  = std::cos(th);
    real_type tg = std::tan(th)/2;

    real_type p0[2], p2[2];
    p0[0] = x0; p0[1] = y0;

    knots[0] = knots[1] = knots[2] = 0;
    Poly[0][0] = p0[0];
    Poly[0][1] = p0[1];
    Poly[0][2] = 1;

    real_type s  = 0;
    real_type ds = L/ns;
    int_type  kk = 0;
    for ( int_type i = 0; i < ns; ++i ) {
      s += ds;
      eval( s, p2[0], p2[1] );

      real_type nx = p0[1]-p2[1];
      real_type ny = p2[0]-p0[0];
      real_type xm = (p0[0]+p2[0])/2;
      real_type ym = (p0[1]+p2[1])/2;

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
  CircleArc::lenTolerance( real_type tol ) const {
    real_type absk = std::abs(k);
    real_type tmp  = absk*tol;
    if ( tmp > 0 ) {
      real_type dtheta = 2*(m_pi-acos(tmp-1));
      return dtheta/absk;
    } else {
      return L;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, CircleArc const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nk      = " << c.k
           << "\nL      = " << c.L
           << "\n";
    return stream;
  }

}

// EOF: Circle.cc
