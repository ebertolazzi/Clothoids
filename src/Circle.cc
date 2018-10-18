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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y( real_type s ) const {
    real_type sk = (s*k)/2;
    return x0+s*Sinc(sk)*sin(theta0+sk);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::X_D( real_type s ) const {
    return cos(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y_D( real_type s ) const {
    return sin(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::X_DD( real_type s ) const {
    return -k*sin(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y_DD( real_type s ) const {
    return k*cos(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::X_DDD( real_type s ) const {
    return -(k*k)*cos(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::Y_DDD( real_type s ) const {
    return -(k*k)*sin(theta0+s*k);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::XY( real_type   s,
                 real_type & x,
                 real_type & y ) const {
    real_type sk  = (s*k)/2;
    real_type tmp = s*Sinc(sk);
    x = x0+tmp*cos(theta0+sk);
    y = y0+tmp*sin(theta0+sk);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::XY( real_type   s,
                 real_type   t,
                 real_type & x,
                 real_type & y ) const {
    real_type sk  = (s*k)/2;
    real_type tmp = s*Sinc(sk);
    real_type th  = theta(s);
    real_type th0 = theta0+sk;
    x = x0+tmp*cos(th0)-t*sin(th);
    y = y0+tmp*sin(th0)+t*cos(th);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::TG( real_type s, real_type & tx, real_type & ty ) const {
    real_type th = theta(s);
    tx = cos(th);
    ty = sin(th);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::NOR( real_type s, real_type & nx, real_type & ny ) const {
    real_type th = theta(s);
    nx = -sin(th);
    ny = cos(th);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::NOR_D( real_type s, real_type & nx_D, real_type & ny_D ) const {
    real_type th = theta(s);
    nx_D = -cos(th)*k;
    ny_D = -sin(th)*k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::NOR_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const {
    real_type th = theta(s);
    real_type k2 = k*k;
    real_type S  = sin(th);
    real_type C  = cos(th);
    nx_DD =  S*k2;
    ny_DD = -C*k2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::NOR_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const {
    real_type th = theta(s);
    real_type S  = sin(th);
    real_type C  = cos(th);
    real_type k3 = k*k*k;
    nx_DDD = C*k3;
    ny_DDD = S*k3;
  }

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
  CircleArc::eval( real_type   s,
                   real_type   t,
                   real_type & x,
                   real_type & y ) const {
    real_type sk  = (s*k)/2;
    real_type LS  = s*Sinc(sk);
    real_type arg = theta0+sk;
    real_type nx, ny;
    NOR( s, nx, ny );
    x = x0 + LS*cos(arg) + nx * t;
    y = y0 + LS*sin(arg) + ny * t;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_D( real_type   s,
                     real_type   t,
                     real_type & x_D,
                     real_type & y_D ) const {
    real_type arg = theta0+s*k;
    real_type nx_D, ny_D;
    NOR_D( s, nx_D, ny_D );
    x_D = cos(arg)+t*nx_D;
    y_D = sin(arg)+t*ny_D;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DD( real_type   s,
                      real_type   t,
                      real_type & x_DD,
                      real_type & y_DD ) const {
    real_type arg = theta0+s*k;
    real_type nx_DD, ny_DD;
    NOR_DD( s, nx_DD, ny_DD );
    x_DD = -k*sin(arg)+t*nx_DD;
    y_DD =  k*cos(arg)+t*ny_DD;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::eval_DDD( real_type   s,
                       real_type   t,
                       real_type & x_DDD,
                       real_type & y_DDD ) const {
    real_type arg = theta0+s*k;
    real_type k2  = k*k;
    real_type nx_DDD, ny_DDD;
    NOR_DDD( s, nx_DDD, ny_DDD );
    x_DDD = -k2*cos(arg)+t*nx_DDD;
    y_DDD = -k2*sin(arg)+t*ny_DDD;
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
    k = -k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  CircleArc::closestPoint( real_type   qx,
                           real_type   qy,
                           real_type & X,
                           real_type & Y,
                           real_type & S ) const {

    S = projectPointOnCircle( x0, y0, cos(theta0), sin(theta0), k, L, qx, qy );

    if ( S < 0 || S > L ) { // minimum distance at the border
      eval( L, X, Y );
      // costruisco piano
      real_type nx = X-x0;
      real_type ny = Y-y0;
      real_type dx = 2*qx-(x0+X);
      real_type dy = 2*qy-(y0+Y);
      if ( nx*dx + ny*dy > 0 ) {
        S = L;
      } else {
        S = 0;
        X = x0;
        Y = y0;
      }
    } else {
      eval( S, X, Y );
    }

    return hypot(qx-X,qy-Y);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  CircleArc::findST( real_type   x,
                     real_type   y,
                     real_type & s,
                     real_type & t ) const {
    real_type X, Y, nx, ny;
    s = projectPointOnCircle( x0, y0, cos(theta0), sin(theta0), k, L, x, y );
    eval( s, X, Y );
    NOR( s, nx, ny );
    t = nx*(x-X) + ny*(y-Y);
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
  CircleArc::bbTriangle( real_type p0[2],
                         real_type p1[2],
                         real_type p2[2] ) const {
    real_type dtheta = L * k;
    bool ok = std::abs(dtheta) <= m_pi/3;
    if ( ok ) {
      p0[0] = x0; p0[1] = y0;
      eval( L, p2[0], p2[1] );
      p1[0] = (p0[0]+p2[0])/2;
      p1[1] = (p0[1]+p2[1])/2;
      real_type nx = p0[1]-p2[1];
      real_type ny = p2[0]-p0[0];
      real_type tg = tan(dtheta/2)/2;
      p1[0] -= nx * tg;
      p1[1] -= ny * tg;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  CircleArc::toNURBS( real_type knots[],
                      real_type Poly[],
                      bool      get_size ) const {

    real_type dtheta = L*k;
    int_type  ns     = int_type(std::floor(3*std::abs(dtheta)/m_pi));
    if ( ns < 1 ) ns = 1;
    if ( get_size ) return 1+2*ns;

    real_type th = dtheta/(2*ns);
    real_type w  = cos(th);
    real_type tg = tan(th)/2;

    real_type p0[2], p2[2];
    p0[0] = x0; p0[1] = y0;

    knots[0] = knots[1] = knots[2] = 0;
    Poly[0] = p0[0];
    Poly[1] = p0[1];
    Poly[2] = 1;

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
      Poly[kk*3+0] = w*(xm - nx * tg);
      Poly[kk*3+1] = w*(ym - ny * tg);
      Poly[kk*3+2] = w;

      ++kk;
      Poly[kk*3+0] = p2[0];
      Poly[kk*3+1] = p2[1];
      Poly[kk*3+2] = 1;

      knots[kk+1] = i+1;
      knots[kk+2] = i+1;

      p0[0] = p2[0];
      p0[1] = p2[1];

    }
    knots[kk+3] = ns;
    return 1+2*ns;
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
