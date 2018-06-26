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

///
/// file: Line.hh
///

#ifndef LINE_HH
#define LINE_HH

#include "G2lib.hh"

namespace G2lib {

  /*\
   |   _     _
   |  | |   (_)_ __   ___
   |  | |   | | '_ \ / _ \
   |  | |___| | | | |  __/
   |  |_____|_|_| |_|\___| 
  \*/

  class LineSegment {

    real_type x0,      //!< initial x coordinate of the line
              y0,      //!< initial y coordinate of the line
              theta0 ; //!< angle of the line

    real_type c0,     //!< `cos(theta0)`
              s0,     //!< `sin(theta0)`
              L ;     //!< length of the segment

  public:

    LineSegment()
    : x0(0)
    , y0(0)
    , theta0(0)
    , c0(1)
    , s0(0)
    , L(0)
    {}

    //! construct a circle curve with the standard parameters
    LineSegment( real_type _x0,
                 real_type _y0,
                 real_type _theta0,
                 real_type _L )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , c0(cos(_theta0))
    , s0(sin(_theta0))
    , L(_L)
    {}

    void
    copy( LineSegment const & c ) {
      x0     = c.x0 ;
      y0     = c.y0 ;
      theta0 = c.theta0 ;
      c0     = c.c0 ;
      s0     = c.s0 ;
      L      = c.L ;
    }

    LineSegment( LineSegment const & s ) { copy(s) ; }

    LineSegment const & operator = ( LineSegment const & s )
    { copy(s) ; return *this ; }

    real_type xBegin()   const { return x0 ; }
    real_type yBegin()   const { return y0 ; }
    real_type theta()    const { return theta0 ; }
    real_type sinTheta() const { return s0 ; }
    real_type cosTheta() const { return c0 ; }
    real_type length()   const { return L ; }

    void
    build( real_type _x0,
           real_type _y0,
           real_type _theta0,
           real_type _L ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      c0     = cos(_theta0);
      s0     = sin(_theta0);
      L      = _L ;
    }

    //! construct a clothoid with the standard parameters
    void
    build_2P( real_type _x0,
              real_type _y0,
              real_type _x1,
              real_type _y1 ) ;

    real_type X( real_type s ) const { return x0 + c0 * s ; }
    real_type Y( real_type s ) const { return y0 + s0 * s ; }

    void
    eval( real_type s, real_type & x, real_type & y ) const {
      x = x0 + c0 * s ;
      y = y0 + s0 * s ;
    }

    void
    eval_D( real_type, real_type & x_D, real_type & y_D ) const {
      x_D = c0 ;
      y_D = s0 ;
    }

    void
    eval_DD( real_type, real_type & x_DD, real_type & y_DD ) const {
      x_DD = y_DD = 0 ;
    }

    void
    eval_DDD( real_type, real_type & x_DDD, real_type & y_DDD ) const {
      x_DDD = y_DDD = 0 ;
    }

    void
    trim( real_type s_begin, real_type s_end ) {
      x0 += c0 * s_begin ;
      y0 += s0 * s_begin ;
      L   = s_end - s_begin ;
    }

    void
    translate( real_type tx, real_type ty )
    { x0 += tx ; y0 += ty ; }

    void
    rotate( real_type angle, real_type cx, real_type cy ) ;

    void
    reverse() ;

    void
    changeOrigin( real_type newx0, real_type newy0 )
    { x0 = newx0 ; y0 = newy0 ; }

    /*!
     * \brief compute the point at minimum distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param X x-coordinate of the closest point
     * \param Y y-coordinate of the closest point
     * \param S param of the closest point
     * \return the distance point-segment
    \*/
    real_type
    closestPoint( real_type   x,
                  real_type   y,
                  real_type & X,
                  real_type & Y,
                  real_type & S ) const ;

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param S param at minimum distance
     * \return the distance point-segment
    \*/
    real_type
    distance( real_type   x,
              real_type   y,
              real_type & S ) const {
      real_type X, Y ;
      return closestPoint( x, y, X, Y, S );
    }

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \return the distance point-segment
    \*/
    real_type
    distance( real_type x, real_type y ) const {
      real_type ss ;
      return distance( x, y, ss );
    }

    /*! \brief Compute rational B-spline coefficients for a line segment
     *
     * \param knots  knots of the B-spline
     * \param Poly   polygon of the B-spline
     * \return       3 the number of polygon points
     */

    int
    toNURBS( real_type knots[5], real_type Poly[2][3] ) const ;

    friend
    std::ostream &
    operator << ( std::ostream & stream, LineSegment const & c ) ;

  } ;

}

#endif

///
/// eof: Line.hh
///
