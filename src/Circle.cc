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

namespace Circle {
  static const valueType m_pi = 3.14159265358979323846264338328  ; // pi

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
    k      = tan(th) / d ;
    s_min  = 0;
    s_max  = 2*d*cos(th)/Sinc(th);
  }

  valueType
  CircleArc::thetaMinMax( valueType & thMin, valueType & thMax ) const  {
    thMin = theta0 + s_min * k ;
    thMax = theta0 + s_max * k ;
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

  //! set the origin of the clothoid to the curvilinear abscissa s0
  void
  CircleArc::change_origin( valueType s0 ) {
    valueType sL = s0*k - theta0 ;
    x0     += s0 * Sinc( sL ) ;
    y0     += s0 * Cosc( sL ) ;
    theta0 += s0 * k ;
    s_min -= s0 ;
    s_max -= s0 ;
    this->s0 = sin(theta0) ;
    this->c0 = cos(theta0) ;
  }

  //! get the bounding box triangle (if angle variation less that pi/3)
  bool
  CircleArc::bbTriangle( valueType p0[2],
                         valueType p1[2],
                         valueType p2[2] ) const {
    valueType dtheta = (s_max - s_min) * k ;
    bool ok = std::abs(dtheta) <= m_pi/3 ;
    if ( ok ) {
      eval( s_min, p0[0], p0[1] );
      eval( s_max, p2[0], p2[1] );
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

    valueType dtheta = (s_max-s_min)*k ;
    indexType ns     = std::floor(3*std::abs(dtheta)/m_pi) ;
    if      ( ns < 1 ) ns = 1 ;
    else if ( ns > 4 ) ns = 4 ;

    valueType th = dtheta/(2*ns) ;
    valueType w  = cos(th) ;
    valueType tg = tan(th)/2;

    valueType p0[2], p2[2] ;
    eval( s_min, p0[0], p0[1] );

    knots[0] = knots[1] = knots[2] = 0 ;
    Poly[0][0] = p0[0] ;
    Poly[0][1] = p0[1] ;
    Poly[0][2] = 1  ;

    valueType s  = s_min ;
    valueType ds = (s_max-s_min)/ns ;
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
           << "\nL      = " << c.s_max-c.s_min
           << "\ns_min  = " << c.s_min
           << "\ns_max  = " << c.s_max
           << "\n" ;
    return stream ;
  }

}

// EOF: Circle.cc
