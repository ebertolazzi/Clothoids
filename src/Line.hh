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

    valueType x0,      //!< initial x coordinate of the line
              y0,      //!< initial y coordinate of the line
              theta0 ; //!< angle of the line

    valueType c0,     //!< `cos(theta0)`
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
    LineSegment( valueType _x0,
                 valueType _y0,
                 valueType _theta0,
                 valueType _L )
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

    valueType getX0()        const { return x0 ; }
    valueType getY0()        const { return y0 ; }
    valueType getTheta0()    const { return theta0 ; }
    valueType getSinTheta0() const { return s0 ; }
    valueType getCosTheta0() const { return c0 ; }
    valueType getL()         const { return L ; }

    void
    build( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _L ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      c0     = cos(_theta0);
      s0     = sin(_theta0);
      L      = _L ;
    }

    //! construct a clothoid with the standard parameters
    void
    build_2P( valueType _x0,
              valueType _y0,
              valueType _x1,
              valueType _y1 ) ;

    valueType X( valueType s ) const { return x0 + c0 * s ; }
    valueType Y( valueType s ) const { return y0 + s0 * s ; }

    void
    eval( valueType s, valueType & x, valueType & y ) const {
      x = x0 + c0 * s ;
      y = y0 + s0 * s ;
    }

    void
    eval_D( valueType, valueType & x_D, valueType & y_D ) const {
      x_D = c0 ;
      y_D = s0 ;
    }

    void
    eval_DD( valueType, valueType & x_DD, valueType & y_DD ) const {
      x_DD = y_DD = 0 ;
    }

    void
    eval_DDD( valueType, valueType & x_DDD, valueType & y_DDD ) const {
      x_DDD = y_DDD = 0 ;
    }

    void
    trim( valueType s_begin, valueType s_end ) {
      x0 += c0 * s_begin ;
      y0 += s0 * s_begin ;
      L   = s_end - s_begin ;
    }

    void
    translate( valueType tx, valueType ty )
    { x0 += tx ; y0 += ty ; }

    void
    rotate( valueType angle, valueType cx, valueType cy ) ;

    void
    reverse() ;

    void
    changeOrigin( valueType newx0, valueType newy0 )
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
    valueType
    closestPoint( valueType   x,
                  valueType   y,
                  valueType & X,
                  valueType & Y,
                  valueType & S ) const ;

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param S param at minimum distance
     * \return the distance point-segment
    \*/
    valueType
    distance( valueType   x,
              valueType   y,
              valueType & S ) const {
      valueType X, Y ;
      return closestPoint( x, y, X, Y, S );
    }

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \return the distance point-segment
    \*/
    valueType
    distance( valueType x, valueType y ) const {
      valueType ss ;
      return distance( x, y, ss );
    }

    /*! \brief Compute rational B-spline coefficients for a line segment
     *
     * \param knots  knots of the B-spline
     * \param Poly   polygon of the B-spline
     * \return       3 the number of polygon points
     */

    int
    toNURBS( valueType knots[5], valueType Poly[2][3] ) const ;

    friend
    std::ostream &
    operator << ( std::ostream & stream, LineSegment const & c ) ;

  } ;

}

#endif

///
/// eof: Line.hh
///
