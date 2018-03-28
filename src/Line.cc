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

#include "Line.hh"

namespace G2lib {

  /*\
   |   _     _
   |  | |   (_)_ __   ___
   |  | |   | | '_ \ / _ \
   |  | |___| | | | |  __/
   |  |_____|_|_| |_|\___| 
  \*/

  void
  LineSegment::build_2P( valueType _x0,
                         valueType _y0,
                         valueType _x1,
                         valueType _y1 ) {
    valueType dx = _x1-_x0 ;
    valueType dy = _y1-_y0 ;
    valueType d  = hypot( dx, dy ) ;
    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = atan2(dy, dx);
    c0     = dx / d ;
    s0     = dy / d ;
    L      = 0 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LineSegment::rotate( valueType angle, valueType cx, valueType cy ) {
    valueType dx  = x0 - cx ;
    valueType dy  = y0 - cy ;
    valueType C   = cos(angle) ;
    valueType S   = sin(angle) ;
    valueType ndx = C*dx - S*dy ;
    valueType ndy = C*dy + S*dx ;
    x0      = cx + ndx ;
    y0      = cy + ndy ;
    theta0 += angle ;
    c0      = cos(theta0);
    s0      = sin(theta0);

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LineSegment::reverse() {
    x0     += c0 * L ;
    y0     += s0 * L ;
    c0      = -c0 ;
    s0      = -s0 ;
    theta0 += m_pi ;
    if ( theta0 > m_pi ) theta0 -= 2*m_pi ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  LineSegment::closestPoint( valueType   x,
                             valueType   y,
                             valueType & X,
                             valueType & Y,
                             valueType & S ) const {

    S = projectPointOnLine( x0, y0, c0, s0, x, y ) ;

    if ( S <= 0 ) { // distanza sul bordo 0
      S = 0 ;
      X = x0 ;
      Y = y0 ;
    } else {
      if ( S >= L ) S = L ;
      eval( S, X, Y );
    }
    return hypot( x-X, y-Y ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int
  LineSegment::toNURBS( valueType knots[4], valueType Poly[2][3] ) const {
    knots[0] = knots[1] = 0 ;
    knots[2] = knots[3] = 1 ;
    Poly[0][0] = x0 ;
    Poly[0][1] = y0 ;
    Poly[0][2] = 1  ;
    Poly[1][0] = x0+L*c0 ;
    Poly[1][1] = y0+L*s0 ;
    Poly[1][2] = 1  ;
    return 2 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::ostream &
  operator << ( std::ostream & stream, LineSegment const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nL      = " << c.L
           << "\n" ;
    return stream ;
  }

}

// EOF: Line.cc
