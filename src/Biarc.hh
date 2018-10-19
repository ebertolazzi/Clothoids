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

  /*!
   * \brief Compute biarc fitting by Hemite data
   *
  \*/

  class Biarc {
    CircleArc C0, C1;
    real_type omega;

    void
    gfun( real_type alpha, real_type g[3] ) const {
      real_type so  = sin(alpha);
      real_type co  = cos(alpha);
      real_type oco = alpha*co;
      g[0] = so + oco;
      g[1] = 2*co - alpha*so;
      g[2] = -3*so - oco;
    }

  public:

    Biarc()
    : omega(0)
    {}

    //! construct a clothoid with the standard parameters
    Biarc( real_type x0,
           real_type y0,
           real_type theta0,
           real_type x1,
           real_type y1,
           real_type theta1 ) {
      bool ok = build( x0, y0, theta0, x1, y1, theta1 );
      G2LIB_ASSERT( ok, "Biarc( x0 = " << x0     <<
                        ", y0 = "      << y0     <<
                        ", theta0 = "  << theta0 <<
                        ", x1 = "      << x1     <<
                        ", y1 = "      << y1     <<
                        ", theta1 = "  << theta1 <<
                        ") cannot be computed" );
    }

    void
    copy( Biarc const & c ) {
      C0.copy(c.C0);
      C1.copy(c.C1);
      omega = c.omega;
    }

    Biarc( Biarc const & ba ) { copy(ba); }

    Biarc const & operator = ( Biarc const & ba )
    { copy(ba); return *this; }

    CircleArc const & getC0() const { return C0; }
    CircleArc const & getC1() const { return C1; }

    //! construct a biarc with the standard parameters
    bool
    build( real_type x0,
           real_type y0,
           real_type theta0,
           real_type x1,
           real_type y1,
           real_type theta1 );

    /*!
    //  \brief
    //  construct a biarc by 3 point at "minimum energy"
    //  - Planar point set fairing and fitting by arc splines
    //  - Xunnian Yang and Guozhao Wang
    //  - Computer-Aided Design, vol 33, 2001
    */
    bool
    build_3P( real_type x0,
              real_type y0,
              real_type x1,
              real_type y1,
              real_type x2,
              real_type y2 );

    real_type kappa    ( real_type s ) const;
    real_type kappa_D  ( real_type s ) const;
    real_type kappa_DD ( real_type s ) const;
    real_type kappa_DDD( real_type s ) const;

    real_type theta    ( real_type s ) const;
    real_type theta_D  ( real_type s ) const;
    real_type theta_DD ( real_type s ) const;
    real_type theta_DDD( real_type s ) const;

    real_type X    ( real_type s ) const;
    real_type X_D  ( real_type s ) const;
    real_type X_DD ( real_type s ) const;
    real_type X_DDD( real_type s ) const;

    real_type Y    ( real_type s ) const;
    real_type Y_D  ( real_type s ) const;
    real_type Y_DD ( real_type s ) const;
    real_type Y_DDD( real_type s ) const;

    real_type tg_x( real_type s ) const { return cos(theta(s)); }
    real_type tg_y( real_type s ) const { return sin(theta(s)); }

    real_type nor_x( real_type s ) const { return -sin(theta(s)); }
    real_type nor_y( real_type s ) const { return cos(theta(s)); }

    void XY( real_type s, real_type & x, real_type & y ) const;
    void XY( real_type s, real_type t, real_type & x, real_type & y ) const;
    void TG( real_type s, real_type & tx, real_type & ty ) const;

    void NOR    ( real_type s, real_type & nx,     real_type & ny ) const;
    void NOR_D  ( real_type s, real_type & nx_D,   real_type & ny_D ) const;
    void NOR_DD ( real_type s, real_type & nx_DD,  real_type & ny_DD ) const;
    void NOR_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const;

    real_type Xstar()     const { return C1.xBegin(); }
    real_type Ystar()     const { return C1.yBegin(); }
    real_type thetaStar() const { return C1.thetaBegin(); }

    void
    eval( real_type   s,
          real_type & th,
          real_type & k,
          real_type & x,
          real_type & y ) const;

    void eval    ( real_type s, real_type & x,     real_type & y ) const;
    void eval_D  ( real_type s, real_type & x_D,   real_type & y_D ) const;
    void eval_DD ( real_type s, real_type & x_DD,  real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type & x_DDD, real_type & y_DDD ) const;

    void eval    ( real_type s, real_type t, real_type & x,     real_type & y ) const;
    void eval_D  ( real_type s, real_type t, real_type & x_D,   real_type & y_D ) const;
    void eval_DD ( real_type s, real_type t, real_type & x_DD,  real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type t, real_type & x_DDD, real_type & y_DDD ) const;

    real_type xBegin0()     const { return C0.xBegin(); }
    real_type xEnd0()       const { return C0.xEnd(); }
    real_type yBegin0()     const { return C0.yBegin(); }
    real_type yEnd0()       const { return C0.yEnd(); }
    real_type thetaBegin0() const { return C0.thetaBegin(); }
    real_type thetaEnd0()   const { return C0.thetaEnd(); }
    real_type kappa0()      const { return C0.kappa(); }
    real_type length0()     const { return C0.length(); }

    real_type xBegin1()     const { return C1.xBegin(); }
    real_type xEnd1()       const { return C1.xEnd(); }
    real_type yBegin1()     const { return C1.yBegin(); }
    real_type yEnd1()       const { return C1.yEnd(); }
    real_type thetaBegin1() const { return C1.thetaBegin(); }
    real_type thetaEnd1()   const { return C1.thetaEnd(); }
    real_type kappa1()      const { return C1.kappa(); }
    real_type length1()     const { return C1.length(); }

    real_type length()      const { return C0.length() + C1.length(); }
    real_type delta_theta() const { return C0.delta_theta() + C1.delta_theta(); }

    void changeOrigin( real_type newx0, real_type newy0 );
    void translate( real_type tx, real_type ty );
    void rotate( real_type angle, real_type cx, real_type cy );
    void reverse();
    void scale( real_type s );

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
    real_type
    closestPoint( real_type   x,
                  real_type   y,
                  real_type & X,
                  real_type & Y,
                  real_type & S ) const;

    /*!
     * \brief compute the distance from a point `[x,y]` and the biarc
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param S param at minimum distance
     * \return the distance point-circle
    \*/
    real_type
    distance( real_type   x,
              real_type   y,
              real_type & S ) const {
      real_type X, Y;
      return closestPoint( x, y, X, Y, S );
    }

    /*!
     * \brief compute the distance from a point `[x,y]` and the biarc
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \return the distance point-circle
    \*/
    real_type
    distance( real_type x, real_type y ) const {
      real_type ss;
      return distance( x, y, ss );
    }

    /*!
     * \brief Find parametric coordinate.
     *
     * \param x x-coordinate point
     * \param y y-coordinate point
     * \param s value \f$ s \f$
     * \param t value \f$ t \f$
    \*/
    void
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const;

    void
    info( ostream_type & stream ) const
    { stream << "BiArc\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, Biarc const & bi );

  };

}

#endif

///
/// eof: Biarc.hh
///
