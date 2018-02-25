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
              s_min,  //!< initial curvilinear coordinate of the line
              s_max ; //!< final curvilinear coordinate of the line

  public:

    LineSegment()
    : x0(0)
    , y0(0)
    , theta0(0)
    , c0(1)
    , s0(0)
    , s_min(0)
    , s_max(0)
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
    , s_min(0)
    , s_max(_L)
    {}

    LineSegment( valueType _x0,
                 valueType _y0,
                 valueType _theta0,
                 valueType _smin ,
                 valueType _smax )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , c0(cos(_theta0))
    , s0(sin(_theta0))
    , s_min(_smin)
    , s_max(_smax)
    {}

    void
    copy( LineSegment const & c ) {
      x0     = c.x0 ;
      y0     = c.y0 ;
      theta0 = c.theta0 ;
      c0     = c.c0 ;
      s0     = c.s0 ;
      s_min  = c.s_min ;
      s_max  = c.s_max ;
    }

    LineSegment( LineSegment const & s ) { copy(s) ; }

    LineSegment const & operator = ( LineSegment const & s )
    { copy(s) ; return *this ; }

    valueType getX0()        const { return x0 ; }
    valueType getY0()        const { return y0 ; }
    valueType getTheta0()    const { return theta0 ; }
    valueType getSinTheta0() const { return s0 ; }
    valueType getCosTheta0() const { return c0 ; }
    valueType getSmin()      const { return s_min ; }
    valueType getSmax()      const { return s_max ; }
    valueType getL()         const { return s_max-s_min ; }

    void
    build( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _smin,
           valueType _smax ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      c0     = cos(_theta0);
      s0     = sin(_theta0);
      s_min  = _smin ;
      s_max  = _smax ;
    }

    //! construct a clothoid with the standard parameters
    void
    build_2P( valueType _x0,
              valueType _y0,
              valueType _x1,
              valueType _y1 ) ;

    valueType
    totalLength() const { return s_max-s_min ; }

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
      s_min = s_begin ;
      s_max = s_end ;
    }

    //! set the origin of the clothoid to the curvilinear abscissa s0
    void
    changeCurvilinearOrigin( valueType s0_new ) {
      x0 += c0 * s0_new ;
      y0 += s0 * s0_new ;
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
