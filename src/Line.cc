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
#include <algorithm>

// Microsoft visual studio Workaround
#ifdef max
  #undef max
#endif

#ifdef min
  #undef min
#endif

namespace G2lib {

  using std::max;

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
    real_type dx = _x1-_x0;
    real_type dy = _y1-_y0;
    L      = hypot( dx, dy );
    x0     = _x0;
    y0     = _y0;
    theta0 = atan2(dy, dx);
    if ( L > 0 ) {
      c0 = dx / L;
      s0 = dy / L;
    } else {
      c0 = s0 = 0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  LineSegment::rotate( real_type angle, real_type cx, real_type cy ) {
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
  LineSegment::reverse() {
    x0     += c0 * L;
    y0     += s0 * L;
    c0      = -c0;
    s0      = -s0;
    theta0 += m_pi;
    if ( theta0 > m_pi ) theta0 -= 2*m_pi;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  LineSegment::closestPoint( real_type   x,
                             real_type   y,
                             real_type & X,
                             real_type & Y,
                             real_type & S ) const {

    S = projectPointOnLine( x0, y0, c0, s0, x, y );

    if ( S <= 0 ) { // distanza sul bordo 0
      S = 0;
      X = x0;
      Y = y0;
    } else {
      if ( S >= L ) S = L;
      eval( S, X, Y );
    }
    return hypot( x-X, y-Y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int
  LineSegment::toNURBS( real_type knots[4], real_type Poly[2][3] ) const {
    knots[0] = knots[1] = 0;
    knots[2] = knots[3] = 1;
    Poly[0][0] = x0;
    Poly[0][1] = y0;
    Poly[0][2] = 1;
    Poly[1][0] = x0+L*c0;
    Poly[1][1] = y0+L*s0;
    Poly[1][2] = 1;
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int
  LineSegment::toBS( real_type knots[4], real_type Poly[2][2] ) const {
    knots[0] = knots[1] = 0;
    knots[2] = knots[3] = 1;
    Poly[0][0] = x0;
    Poly[0][1] = y0;
    Poly[1][0] = x0+L*c0;
    Poly[1][1] = y0+L*s0;
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Given three colinear points p, q, r, the function checks if
  // point q lies on line segment 'pr'
  static
  bool
  onSegment( real_type const p[2],
             real_type const q[2],
             real_type const r[2],
             real_type const epsi ) {

    real_type mi_x, ma_x;
    if ( p[0] > r[0] ) { mi_x = r[0]; ma_x = p[0]; }
    else               { mi_x = p[0]; ma_x = r[0]; }

    bool ok = q[0] <= ma_x+epsi && q[0] >= mi_x-epsi;
    if ( ok ) {
      real_type mi_y, ma_y;
      if ( p[1] > r[1] ) { mi_y = r[1]; ma_y = p[1]; }
      else               { mi_y = p[1]; ma_y = r[1]; }
      ok = q[1] <= ma_y+epsi && q[1] >= mi_y-epsi;
    }
    return ok;
  }

  // To find orientation of ordered triplet (p, q, r).
  // The function returns following values
  // 0 --> p, q and r are collinear
  // 1 --> Clockwise
  // 2 --> Counterclockwise
  static
  int_type
  orientation( real_type const p[2],
               real_type const q[2],
               real_type const r[2],
               real_type const epsi ) {
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    real_type qp_x = q[0] - p[0];
    real_type qp_y = q[1] - p[1];
    real_type rq_x = r[0] - q[0];
    real_type rq_y = r[1] - q[1];

    real_type det = qp_y * rq_x - qp_x * rq_y;
    if ( abs(det) < epsi ) return 0;  // collinear
    return (det > 0)? 1: 2; // clock or counterclock wise
  }

  bool
  LineSegment::intersect( LineSegment const & S,
                          real_type         & s1,
                          real_type         & s2 ) const {
    // The main function that returns true if line segment 'p1q1'
    // and 'p2q2' intersect.
    real_type const p1[2] = { xBegin(),   yBegin()   };
    real_type const q1[2] = { xEnd(),     yEnd()     };
    real_type const p2[2] = { S.xBegin(), S.yBegin() };
    real_type const q2[2] = { S.xEnd(),   S.yEnd()   };
    real_type const epsi  = max(L,S.L)*machepsi100;

    // Find the four orientations needed for general and special cases
    int_type o1 = orientation( p1, q1, p2, epsi );
    int_type o2 = orientation( p1, q1, q2, epsi );
    int_type o3 = orientation( p2, q2, p1, epsi );
    int_type o4 = orientation( p2, q2, q1, epsi );

    // General case
    if ( o1 != o2 && o3 != o4 ) {
      real_type det = c0 * S.s0 - s0 * S.c0;
      real_type px  = p2[0]-p1[0];
      real_type py  = p2[1]-p1[1];
      s1 = (px * S.s0 - py * S.c0)/ det;
      s2 = (px * s0 - py * c0)/ det;
      return true;
    }

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if ( o1 == 0 && onSegment( p1, p2, q1, epsi ) ) {
      s1 = hypot( p2[0]-p1[0], p2[1]-p1[1] );
      s2 = 0;
      return true;
    }

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if ( o2 == 0 && onSegment( p1, q2, q1, epsi ) ) {
      s1 = hypot( q2[0]-p1[0], q2[1]-p1[1] );
      s2 = S.L;
      return true;
    }

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if ( o3 == 0 && onSegment( p2, p1, q2, epsi ) ) {
      s1 = 0;
      s2 = hypot( p1[0]-p2[0], p1[1]-p2[1] );
      return true;
    }

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if ( o4 == 0 && onSegment( p2, q1, q2, epsi ) ) {
      s1 = L;
      s2 = hypot( q1[0]-p2[0], q1[1]-p2[1] );
      return true;
    }

    s1 = s2 = 0;
    return false; // Doesn't fall in any of the above cases
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  LineSegment::intersect( LineSegment const & S ) const {
    // The main function that returns true if line segment 'p1q1'
    // and 'p2q2' intersect.
    real_type const p1[2] = { xBegin(),   yBegin()   };
    real_type const q1[2] = { xEnd(),     yEnd()     };
    real_type const p2[2] = { S.xBegin(), S.yBegin() };
    real_type const q2[2] = { S.xEnd(),   S.yEnd()   };
    real_type const epsi  = std::max(L,S.L)*machepsi100;

    // Find the four orientations needed for general and special cases
    int_type o1 = orientation( p1, q1, p2, epsi );
    int_type o2 = orientation( p1, q1, q2, epsi );
    int_type o3 = orientation( p2, q2, p1, epsi );
    int_type o4 = orientation( p2, q2, q1, epsi );

    // General case
    if ( o1 != o2 && o3 != o4 ) return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if ( o1 == 0 && onSegment( p1, p2, q1, epsi ) ) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if ( o2 == 0 && onSegment( p1, q2, q1, epsi ) ) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if ( o3 == 0 && onSegment( p2, p1, q2, epsi ) ) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if ( o4 == 0 && onSegment( p2, q1, q2, epsi ) ) return true;

    return false; // Doesn't fall in any of the above cases
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, LineSegment const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nL      = " << c.L
           << "\n";
    return stream;
  }

}

// EOF: Line.cc
