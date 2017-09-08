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

namespace Line {

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
    s_min  = 0 ;
    s_max  = d ;
  }

  int
  LineSegment::toNURBS( valueType knots[5], valueType Poly[2][3] ) const {
    valueType L = s_max-s_min ;
    knots[0] = knots[1] = knots[2] = 0 ;
    knots[3] = knots[4] = knots[5] = 1;
    Poly[0][0] = x0 ;
    Poly[0][1] = y0 ;
    Poly[0][2] = 1  ;
    Poly[1][0] = x0+(L/2)*c0 ;
    Poly[1][1] = y0+(L/2)*s0 ;
    Poly[1][2] = 1  ;
    Poly[2][0] = x0+L*c0 ;
    Poly[2][1] = y0+L*s0 ;
    Poly[2][2] = 1  ;
    return 3 ;
  }

  std::ostream &
  operator << ( std::ostream & stream, LineSegment const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nL      = " << c.s_max-c.s_min
           << "\ns_min  = " << c.s_min
           << "\ns_max  = " << c.s_max
           << "\n" ;
    return stream ;
  }

}

// EOF: Line.cc
