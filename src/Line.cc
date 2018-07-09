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
  LineSegment::build_2P( real_type _x0,
                         real_type _y0,
                         real_type _x1,
                         real_type _y1 ) {
    real_type dx = _x1-_x0 ;
    real_type dy = _y1-_y0 ;
    L      = hypot( dx, dy ) ;
    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = atan2(dy, dx);
    if ( L > 0 ) {
      c0 = dx / L ;
      s0 = dy / L ;
    } else {
      c0 = s0 = 0 ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LineSegment::rotate( real_type angle, real_type cx, real_type cy ) {
    real_type dx  = x0 - cx ;
    real_type dy  = y0 - cy ;
    real_type C   = cos(angle) ;
    real_type S   = sin(angle) ;
    real_type ndx = C*dx - S*dy ;
    real_type ndy = C*dy + S*dx ;
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

  real_type
  LineSegment::closestPoint( real_type   x,
                             real_type   y,
                             real_type & X,
                             real_type & Y,
                             real_type & S ) const {

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
  LineSegment::toNURBS( real_type knots[4], real_type Poly[2][3] ) const {
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

  int
  LineSegment::toBS( real_type knots[4], real_type Poly[2][2] ) const {
    knots[0] = knots[1] = 0 ;
    knots[2] = knots[3] = 1 ;
    Poly[0][0] = x0 ;
    Poly[0][1] = y0 ;
    Poly[1][0] = x0+L*c0 ;
    Poly[1][1] = y0+L*s0 ;
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
