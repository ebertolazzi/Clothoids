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
/// file: Circle.hh
///

#ifndef CIRCLE_HH
#define CIRCLE_HH

#include "G2lib.hh"

//! Clothoid computations routine
namespace G2lib {

  /*\
   |    ____ _          _
   |   / ___(_)_ __ ___| | ___  ___
   |  | |   | | '__/ __| |/ _ \/ __|
   |  | |___| | | | (__| |  __/\__ \
   |   \____|_|_|  \___|_|\___||___/
  \*/

  void
  CircleTangentPoints( valueType PA[2],
                       valueType rA,
                       valueType PB[2],
                       valueType rB,
                       bool &    external_tangents,
                       valueType PTE0[2][2],
                       valueType PTE1[2][2],
                       bool &    internal_tangents,
                       valueType PTI0[2][2],
                       valueType PTI1[2][2] ) ;

  bool
  CircleLineTransition( valueType C[2],
                        valueType r,
                        valueType P[2],
                        valueType theta,
                        valueType C0[2],
                        valueType C1[2] ) ;

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  //! \brief Class to manage Clothoid Curve
  class CircleArc {

    valueType x0,     //!< initial x coordinate of the clothoid
              y0,     //!< initial y coordinate of the clothoid
              theta0, //!< initial angle of the clothoid
              k ;     //!< curvature

    valueType c0,     //!< cos(theta0)
              s0,     //!< sin(theta0)
              L ;     //!< length of the circle segment

  public:

    CircleArc()
    : x0(0)
    , y0(0)
    , theta0(0)
    , k(0)
    , c0(1)
    , s0(0)
    , L(0)
    {}

    //! construct a circle curve with the standard parameters
    CircleArc( valueType _x0,
               valueType _y0,
               valueType _theta0,
               valueType _k,
               valueType _L )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , k(_k)
    , c0(cos(_theta0))
    , s0(sin(_theta0))
    , L(_L)
    {}

    void
    copy( CircleArc const & c ) {
      x0     = c.x0 ;
      y0     = c.y0 ;
      theta0 = c.theta0 ;
      c0     = c.c0 ;
      s0     = c.s0 ;
      k      = c.k ;
      L      = c.L ;
    }

    CircleArc( CircleArc const & s ) { copy(s) ; }

    CircleArc const & operator = ( CircleArc const & s )
    { copy(s) ; return *this ; }

    valueType sinTheta0() const { return s0 ; }
    valueType cosTheta0() const { return c0 ; }
    valueType kappa()     const { return k ; }
    valueType length()    const { return L ; }

    valueType xBegin()     const { return x0 ; }
    valueType yBegin()     const { return y0 ; }
    valueType thetaBegin() const { return theta0 ; }

    valueType xEnd()     const { return X(L) ; }
    valueType yEnd()     const { return Y(L) ; }
    valueType thetaEnd() const { return theta(L) ; }

    //! construct a circle with the standard parameters
    void
    build( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _k,
           valueType _L ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      c0     = cos(_theta0);
      s0     = sin(_theta0);
      k      = _k ;
      L      = _L ;
    }

    //! build a circle by solving the hermite G1 problem
    bool
    build_G1( valueType _x0,
              valueType _y0,
              valueType _theta0,
              valueType _x1,
              valueType _y1 );

    //! build a circle passing by 3 points
    bool
    build_3P( valueType _x0,
              valueType _y0,
              valueType _x1,
              valueType _y1,
              valueType _x2,
              valueType _y2 );

    valueType
    delta_theta() const { return L*k ; }

    valueType
    theta( valueType s ) const { return theta0 + s*k ; }

    valueType
    theta_D( valueType ) const { return k ; }

    valueType
    theta_DD( valueType ) const { return 0 ; }

    valueType
    theta_DDD( valueType ) const { return 0 ; }

    valueType
    totalLength() const { return L ; }

    valueType
    thetaTotalVariation() const
    { return std::abs(L*k) ; }

    valueType
    thetaMinMax( valueType & thMin, valueType & thMax ) const ;

    valueType
    deltaTheta() const
    { valueType thMin, thMax ; return thetaMinMax( thMin, thMax ) ; }

    valueType X( valueType s ) const ;
    valueType Y( valueType s ) const ;

    void eval( valueType s, valueType & x, valueType & y ) const ;
    void eval_D( valueType s, valueType & x_D, valueType & y_D ) const ;
    void eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const ;
    void eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const ;

    void
    trim( valueType s_begin, valueType s_end ) ;

    void
    changeCurvilinearOrigin( valueType s0, valueType newL ) ;

    void
    changeOrigin( valueType newx0, valueType newy0 )
    { x0 = newx0 ; y0 = newy0 ; }

    void
    rotate( valueType angle, valueType cx, valueType cy );

    void
    scale( valueType s );

    void
    reverse();

    //! get the bounding box triangle (if angle variation less that pi/3)
    bool
    bbTriangle( valueType p0[2],
                valueType p1[2],
                valueType p2[2] ) const ;

    void
    translate( valueType tx, valueType ty )
    { x0 += tx ; y0 += ty ; }

    /*!
     * \brief compute the point at minimum distance from a point `[x,y]` and the circle arc
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param X x-coordinate of the closest point
     * \param Y y-coordinate of the closest point
     * \param S param of the closest point
     * \return the distance point-circle
    \*/
    valueType
    closestPoint( valueType   x,
                  valueType   y,
                  valueType & X,
                  valueType & Y,
                  valueType & S ) const ;

    /*!
     * \brief compute the distance from a point `[x,y]` and the circle arc
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param S param at minimum distance
     * \return the distance point-circle
    \*/
    valueType
    distance( valueType   x,
              valueType   y,
              valueType & S ) const {
      valueType X, Y ;
      return closestPoint( x, y, X, Y, S ) ;
    }

    /*!
     * \brief compute the distance from a point `[x,y]` and the circle arc
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \return the distance point-circle
    \*/
    valueType
    distance( valueType x, valueType y ) const {
      valueType ss ;
      return distance( x, y, ss );
    }

    /*! \brief Compute rational B-spline coefficients for a circle arc
     *
     * \param knots  knots of the B-spline
     * \param Poly   polygon of the B-spline
     * \return       3 up to 9 the number of polygon points
    \*/

    indexType
    toNURBS( valueType knots[], valueType Poly[], bool get_size ) const ; // Poly 3 x n matrix

    friend
    std::ostream &
    operator << ( std::ostream & stream, CircleArc const & c ) ;

  } ;

}

#endif

///
/// eof: Circle.hh
///
