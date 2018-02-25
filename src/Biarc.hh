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

    valueType Xstar()     const { return xs ; }
    valueType Ystar()     const { return ys ; }
    valueType thetaStar() const { return thetas ; }

    void eval( valueType s, valueType & x, valueType & y ) const ;
    void eval_D( valueType s, valueType & x_D, valueType & y_D ) const ;
    void eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const ;
    void eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const ;

    valueType getL() const { return C0.getL() + C1.getL() ; }
    valueType delta_theta() const { return C0.delta_theta() + C1.delta_theta() ; }

    friend
    std::ostream &
    operator << ( std::ostream & stream, Biarc const & bi ) ;

  } ;

}

#endif

///
/// eof: Biarc.hh
///
