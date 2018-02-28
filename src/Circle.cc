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

  using namespace std ;

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  void
  CircleArc::build_G1( valueType _x0,
                       valueType _y0,
                       valueType _theta0,
                       valueType _x1,
                       valueType _y1 ) {

    valueType dx = _x1 - _x0 ;
    valueType dy = _y1 - _y0 ;
    valueType d  = hypot( dx, dy );
    valueType th = atan2( dy, dx ) - _theta0 ;

    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = _theta0 ;
    c0     = cos(_theta0);
    s0     = sin(_theta0);
    k      = 2*sin(th)/d ;
    L      = d/Sinc(th);
  }

  void
  CircleArc::build_3P( valueType _x0,
                       valueType _y0,
                       valueType _x1,
                       valueType _y1,
                       valueType _x2,
                       valueType _y2 ) {

    valueType _dx = _x1 - _x0 ;
    valueType _dy = _y1 - _y0 ;
    valueType _theta0 = atan2( _dy, _dx ) + m_pi_2 ;
    build_G1( _x0, _y0, _theta0, _x2, _y2 );
  }

  valueType
  CircleArc::thetaMinMax( valueType & thMin, valueType & thMax ) const  {
    thMin = theta0 ;
    thMax = theta0 + L * k ;
    if ( thMax < thMin ) std::swap( thMin, thMax ) ;
    return thMax-thMin ;
  }

  valueType
  CircleArc::X( valueType s ) const {
    valueType sk = s*k;
    valueType S  = Sinc(sk);
    valueType C  = Cosc(sk);
    return x0+s*(c0*S-s0*C);
  }

  valueType
  CircleArc::Y( valueType s ) const {
    valueType sk = s*k;
    valueType S  = Sinc(sk);
    valueType C  = Cosc(sk);
    return y0+s*(c0*C+s0*S);
  }

  void
  CircleArc::eval( valueType s, valueType & x, valueType & y ) const {
    valueType sk = s*k;
    valueType S  = Sinc(sk);
    valueType C  = Cosc(sk);
    x = x0+s*(c0*S-s0*C);
    y = y0+s*(c0*C+s0*S);
  }

  void
  CircleArc::eval_D( valueType s, valueType & x_D, valueType & y_D ) const {
    valueType sk  = s*k;
    valueType S   = Sinc(sk);
    valueType C   = Cosc(sk);
    valueType S_D = Sinc_D(sk);
    valueType C_D = Cosc_D(sk);
    x_D = (c0*S-s0*C)+sk*(c0*S_D-s0*C_D);
    y_D = (c0*C+s0*S)+sk*(c0*C_D+s0*S_D);
  }

  void
  CircleArc::eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const {
    valueType sk   = s*k;
    valueType S_D  = Sinc_D(sk);
    valueType C_D  = Cosc_D(sk);
    valueType S_DD = Sinc_DD(sk);
    valueType C_DD = Cosc_DD(sk);
    x_DD = k*(2*(c0*S_D-s0*C_D)+sk*(c0*S_DD-s0*C_DD));
    y_DD = k*(2*(c0*C_D+s0*S_D)+sk*(c0*C_DD+s0*S_DD));
  }

  void
  CircleArc::eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const {
    valueType sk    = s*k;
    valueType k2    = k*k;
    valueType S_DD  = Sinc_DD(sk);
    valueType C_DD  = Cosc_DD(sk);
    valueType S_DDD = Sinc_DDD(sk);
    valueType C_DDD = Cosc_DDD(sk);
    x_DDD = k2*(3*(c0*S_DD-s0*C_DD)+sk*(c0*S_DDD-s0*C_DDD));
    y_DDD = k2*(3*(c0*C_DD+s0*S_DD)+sk*(c0*C_DDD+s0*S_DDD));
  }

  void
  CircleArc::trim( valueType s_begin, valueType s_end ) {
    valueType x, y ;
    eval( s_begin, x, y ) ;
    theta0 += s_begin * k ;
    s0 = sin(theta0);
    c0 = cos(theta0);
    L  = s_end - s_begin ;
    x0 = x ;
    y0 = y ;
  }

  void
  CircleArc::rotate( valueType angle, valueType cx, valueType cy ) {
    valueType dx  = x0 - cx ;
    valueType dy  = y0 - cy ;
    valueType C   = cos(angle) ;
    valueType S   = sin(angle) ;
    valueType ndx = C*dx - S*dy ;
    valueType ndy = C*dy + S*dx ;
    x0      = cx + ndx ;
    y0      = cy + ndy ;
    theta0 += angle ;
    c0      = cos(theta0) ;
    s0      = sin(theta0) ;
  }

  void
  CircleArc::scale( valueType s ) {
    k /= s ;
    L *= s ;
  }

  void
  CircleArc::reverse() {
    theta0 = theta0 + m_pi ;
    if ( theta0 > m_pi ) theta0 -= 2*m_pi ;
    c0 = cos(theta0) ;
    s0 = sin(theta0) ;
    k  = -k ;
  }

  valueType
  CircleArc::distance( valueType   x,
                       valueType   y,
                       valueType & s ) const {

    valueType dx  = x0 - x ;
    valueType dy  = y0 - y ;
    valueType a0  = c0 * dy - s0 * dx ;
    valueType b0  = s0 * dy + c0 * dx ;
    valueType tmp = a0*k ;

    if ( 1+2*tmp > 0 ) {

      tmp = b0/(1+tmp) ;
      s   = -tmp*Atanc(tmp*k) ;

    } else {

      valueType om = atan2( b0, a0+1/k ) ;
      if ( k < 0 ) om += m_pi ;
      s = -om/k ;
      valueType circ = 2*m_pi/abs(k);
      while ( s < 0 ) s += circ;
      while ( s > L ) s -= circ;

    }

    valueType xx(0), yy(0) ;
    if ( s < 0 || s > L ) { // distanza sul bordo
      valueType d0 = hypot( x0-x, y0-y ) ;
      eval( L, xx, yy ); valueType d1 = hypot( x-xx,y-yy ) ;
      if ( d0 < d1 ) { s = 0 ; return d0 ; }
      else           { s = L ; return d1 ; }
    }

    eval( s, xx, yy );
    return hypot(x-xx,y-yy) ;
  }

  void
  CircleArc::changeCurvilinearOrigin( valueType s0, valueType newL ) {
    valueType new_x0, new_y0 ;
    eval( s0,  new_x0, new_y0 ) ;
    x0      = new_x0 ;
    y0      = new_y0 ;
    theta0 += k*s0 ;
    c0      = cos(theta0) ;
    s0      = sin(theta0) ;
    L       = newL ;
  }

  //! get the bounding box triangle (if angle variation less that pi/3)
  bool
  CircleArc::bbTriangle( valueType p0[2],
                         valueType p1[2],
                         valueType p2[2] ) const {
    valueType dtheta = L * k ;
    bool ok = std::abs(dtheta) <= m_pi/3 ;
    if ( ok ) {
      p0[0] = x0 ; p0[1] = y0 ;
      eval( L, p2[0], p2[1] );
      p1[0] = (p0[0]+p2[0])/2 ;
      p1[1] = (p0[1]+p2[1])/2 ;
      valueType nx = p0[1]-p2[1] ;
      valueType ny = p2[0]-p0[0] ;
      valueType tg = tan(dtheta/2)/2;
      p1[0] -= nx * tg ;
      p1[1] -= ny * tg ;
    }
    return ok ;
  }

  indexType
  CircleArc::toNURBS(  valueType knots[12], valueType Poly[9][3] ) const {

    valueType dtheta = L*k ;
    indexType ns     = indexType(std::floor(3*std::abs(dtheta)/m_pi)) ;
    if      ( ns < 1 ) ns = 1 ;
    else if ( ns > 4 ) ns = 4 ;

    valueType th = dtheta/(2*ns) ;
    valueType w  = cos(th) ;
    valueType tg = tan(th)/2;

    valueType p0[2], p2[2] ;
    p0[0] = x0 ; p0[1] = y0 ;

    knots[0] = knots[1] = knots[2] = 0 ;
    Poly[0][0] = p0[0] ;
    Poly[0][1] = p0[1] ;
    Poly[0][2] = 1  ;

    valueType s  = 0 ;
    valueType ds = L/ns ;
    indexType kk = 0 ;
    for ( indexType i = 0 ; i < ns ; ++i ) {
      s += ds ;
      eval( s, p2[0], p2[1] );

      valueType nx = p0[1]-p2[1] ;
      valueType ny = p2[0]-p0[0] ;
      valueType xm = (p0[0]+p2[0])/2 ;
      valueType ym = (p0[1]+p2[1])/2 ;

      ++kk;
      Poly[kk][0] = w*(xm - nx * tg) ;
      Poly[kk][1] = w*(ym - ny * tg) ;
      Poly[kk][2] = w ;

      ++kk;
      Poly[kk][0] = p2[0] ;
      Poly[kk][1] = p2[1] ;
      Poly[kk][2] = 1 ;

      knots[kk+1] = i+1 ;
      knots[kk+2] = i+1 ;

      p0[0] = p2[0] ;
      p0[1] = p2[1] ;

    }
    knots[kk+3] = ns ;
    return 1+2*ns;
  }

  std::ostream &
  operator << ( std::ostream & stream, CircleArc const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nk      = " << c.k
           << "\nL      = " << c.L
           << "\n" ;
    return stream ;
  }

}

// EOF: Circle.cc
