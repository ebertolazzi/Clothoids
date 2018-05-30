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
/// file: Biarc.hh
///

#ifndef BIARC_HH
#define BIARC_HH

#include "G2lib.hh"
#include "Circle.hh"

//! Clothoid computations routine
namespace G2lib {

  /*\
   |   ____  _
   |  | __ )(_) __ _ _ __ ___
   |  |  _ \| |/ _` | '__/ __|
   |  | |_) | | (_| | | | (__
   |  |____/|_|\__,_|_|  \___|
  \*/

  /*! \brief Compute biarc fitting by Hemite data
   *
   */

  class Biarc {
    CircleArc C0, C1 ;
    valueType xs, ys, thetas, cs, ss, alpha ;
  public:
  
    Biarc()
    : xs(0)
    , ys(0)
    , thetas(0)
    , cs(0)
    , ss(0)
    , alpha(0)
    {}

    //! construct a clothoid with the standard parameters
    Biarc( valueType x0,
           valueType y0,
           valueType theta0,
           valueType x1,
           valueType y1,
           valueType theta1 ) {
      bool ok = build( x0, y0, theta0, x1, y1, theta1 );
      G2LIB_ASSERT( ok, "Biarc( x0 = " << x0     <<
                        ", y0 = "      << y0     <<
                        ", theta0 = "  << theta0 <<
                        ", x1 = "      << x1     <<
                        ", y1 = "      << y1     <<
                        ", theta1 = "  << theta1 <<
                        ") cannot be computed" ) ;
    }

    void
    copy( Biarc const & c ) {
      C0.copy(c.C0);
      C1.copy(c.C1);
      xs     = c.xs ;
      ys     = c.ys ;
      thetas = c.thetas ;
      cs     = c.cs ;
      ss     = c.ss ;
      alpha  = c.alpha ;
    }

    Biarc( Biarc const & ba ) { copy(ba) ; }

    Biarc const & operator = ( Biarc const & ba )
    { copy(ba) ; return *this ; }

    CircleArc const & getC0() const { return C0 ; }
    CircleArc const & getC1() const { return C1 ; }

    //! construct a clothoid with the standard parameters
    bool
    build( valueType x0,
           valueType y0,
           valueType theta0,
           valueType x1,
           valueType y1,
           valueType theta1 ) ;

    valueType X( valueType s ) const ;
    valueType Y( valueType s ) const ;
    valueType theta( valueType s ) const ;
    valueType kappa( valueType s ) const ;

    valueType Xstar()     const { return xs ; }
    valueType Ystar()     const { return ys ; }
    valueType thetaStar() const { return thetas ; }

    void eval( valueType s, valueType & th, valueType & k, valueType & x, valueType & y ) const ;
    void eval( valueType s, valueType & x, valueType & y ) const ;
    void eval_D( valueType s, valueType & x_D, valueType & y_D ) const ;
    void eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const ;
    void eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const ;

    valueType xBegin0()     const { return C0.xBegin(); }
    valueType xEnd0()       const { return C0.xEnd(); }
    valueType yBegin0()     const { return C0.yBegin(); }
    valueType yEnd0()       const { return C0.yEnd(); }
    valueType thetaBegin0() const { return C0.thetaBegin(); }
    valueType thetaEnd0()   const { return C0.thetaEnd(); }
    valueType kappa0()      const { return C0.kappa(); }
    valueType length0()     const { return C0.length(); }

    valueType xBegin1()     const { return C1.xBegin(); }
    valueType xEnd1()       const { return C1.xEnd(); }
    valueType yBegin1()     const { return C1.yBegin(); }
    valueType yEnd1()       const { return C1.yEnd(); }
    valueType thetaBegin1() const { return C1.thetaBegin(); }
    valueType thetaEnd1()   const { return C1.thetaEnd(); }
    valueType kappa1()      const { return C1.kappa(); }
    valueType length1()     const { return C1.length(); }

    valueType length()      const { return C0.length() + C1.length() ; }
    valueType delta_theta() const { return C0.delta_theta() + C1.delta_theta() ; }

    void changeOrigin( valueType newx0, valueType newy0 );
    void translate( valueType tx, valueType ty );
    void rotate( valueType angle, valueType cx, valueType cy );
    void reverse();
    void scale( valueType s );

    /*!
     * \brief compute the point at minimum distance from a point `[x,y]` and the biarc
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
     * \brief compute the distance from a point `[x,y]` and the biarc
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
     * \brief compute the distance from a point `[x,y]` and the biarc
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

    friend
    std::ostream &
    operator << ( std::ostream & stream, Biarc const & bi ) ;

  } ;

}

#endif

///
/// eof: Biarc.hh
///
