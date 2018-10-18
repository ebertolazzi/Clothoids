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
  CircleTangentPoints( real_type PA[2],
                       real_type rA,
                       real_type PB[2],
                       real_type rB,
                       bool &    external_tangents,
                       real_type PTE0[2][2],
                       real_type PTE1[2][2],
                       bool &    internal_tangents,
                       real_type PTI0[2][2],
                       real_type PTI1[2][2] );

  bool
  CircleLineTransition( real_type C[2],
                        real_type r,
                        real_type P[2],
                        real_type theta,
                        real_type C0[2],
                        real_type C1[2] );

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  //! \brief Class to manage Clothoid Curve
  class CircleArc {

    real_type x0,     //!< initial x coordinate of the clothoid
              y0,     //!< initial y coordinate of the clothoid
              theta0, //!< initial angle of the clothoid
              k;      //!< curvature

    real_type L;      //!< length of the circle segment

  public:

    CircleArc()
    : x0(0)
    , y0(0)
    , theta0(0)
    , k(0)
    , L(0)
    {}

    //! construct a circle curve with the standard parameters
    CircleArc( real_type _x0,
               real_type _y0,
               real_type _theta0,
               real_type _k,
               real_type _L )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , k(_k)
    , L(_L)
    {}

    void
    copy( CircleArc const & c ) {
      x0     = c.x0;
      y0     = c.y0;
      theta0 = c.theta0;
      k      = c.k;
      L      = c.L;
    }

    CircleArc( CircleArc const & s ) { copy(s); }

    CircleArc const &
    operator = ( CircleArc const & s )
    { copy(s); return *this; }

    real_type sinTheta0() const { return sin(theta0); }
    real_type cosTheta0() const { return cos(theta0); }
    real_type kappa()     const { return k; }
    real_type length()    const { return L; }

    real_type xBegin()     const { return x0; }
    real_type yBegin()     const { return y0; }
    real_type thetaBegin() const { return theta0; }

    real_type xEnd()     const { return X(L); }
    real_type yEnd()     const { return Y(L); }
    real_type thetaEnd() const { return theta(L); }

    // return the length of the arc that can approximated
    // by a line segment
    real_type lenTolerance( real_type tol ) const;

    //! construct a circle with the standard parameters
    void
    build( real_type _x0,
           real_type _y0,
           real_type _theta0,
           real_type _k,
           real_type _L ) {
      x0     = _x0;
      y0     = _y0;
      theta0 = _theta0;
      k      = _k;
      L      = _L;
    }

    //! build a circle by solving the hermite G1 problem
    bool
    build_G1( real_type _x0,
              real_type _y0,
              real_type _theta0,
              real_type _x1,
              real_type _y1 );

    //! build a circle passing by 3 points
    bool
    build_3P( real_type _x0,
              real_type _y0,
              real_type _x1,
              real_type _y1,
              real_type _x2,
              real_type _y2 );

    real_type
    delta_theta() const { return L*k; }

    real_type
    theta( real_type s ) const { return theta0 + s*k; }

    real_type
    theta_D( real_type ) const { return k; }

    real_type
    theta_DD( real_type ) const { return 0; }

    real_type
    theta_DDD( real_type ) const { return 0; }

    real_type
    totalLength() const { return L; }

    real_type
    thetaTotalVariation() const
    { return std::abs(L*k); }

    real_type
    thetaMinMax( real_type & thMin, real_type & thMax ) const;

    real_type
    deltaTheta() const
    { real_type thMin, thMax; return thetaMinMax( thMin, thMax ); }

    real_type X( real_type s ) const;
    real_type X_D( real_type s ) const;
    real_type X_DD( real_type s ) const;
    real_type X_DDD( real_type s ) const;

    real_type Y( real_type s ) const;
    real_type Y_D( real_type s ) const;
    real_type Y_DD( real_type s ) const;
    real_type Y_DDD( real_type s ) const;

    real_type tg_x( real_type s ) const { return cos(theta(s)); }
    real_type tg_y( real_type s ) const { return sin(theta(s)); }

    real_type nor_x( real_type s ) const { return -sin(theta(s)); }
    real_type nor_y( real_type s ) const { return cos(theta(s)); }

    void XY( real_type s, real_type & x, real_type & y ) const;
    void XY( real_type s, real_type t, real_type & x, real_type & y ) const;
    void TG( real_type s, real_type & tx, real_type & ty ) const;
    void NOR( real_type s, real_type & nx, real_type & ny ) const;
    void NOR_D( real_type s, real_type & nx_D, real_type & ny_D ) const;
    void NOR_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const;
    void NOR_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const;

    void eval    ( real_type s, real_type & x,     real_type & y ) const;
    void eval_D  ( real_type s, real_type & x_D,   real_type & y_D ) const;
    void eval_DD ( real_type s, real_type & x_DD,  real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type & x_DDD, real_type & y_DDD ) const;

    void eval    ( real_type s, real_type t, real_type & x,     real_type & y ) const;
    void eval_D  ( real_type s, real_type t, real_type & x_D,   real_type & y_D ) const;
    void eval_DD ( real_type s, real_type t, real_type & x_DD,  real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type t, real_type & x_DDD, real_type & y_DDD ) const;

    void
    trim( real_type s_begin, real_type s_end );

    void
    changeCurvilinearOrigin( real_type s0, real_type newL );

    void
    changeOrigin( real_type newx0, real_type newy0 )
    { x0 = newx0; y0 = newy0; }

    void
    rotate( real_type angle, real_type cx, real_type cy );

    void
    scale( real_type s );

    void
    reverse();

    //! get the bounding box triangle (if angle variation less that pi/3)
    bool
    bbTriangle( real_type p0[2],
                real_type p1[2],
                real_type p2[2] ) const;

    void
    translate( real_type tx, real_type ty )
    { x0 += tx; y0 += ty; }

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
    real_type
    closestPoint( real_type   x,
                  real_type   y,
                  real_type & X,
                  real_type & Y,
                  real_type & S ) const;

    /*!
     * \brief compute the distance from a point `[x,y]` and the circle arc
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
     * \brief compute the distance from a point `[x,y]` and the circle arc
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

    /*! \brief Find parametric coordinate.
     *
     * \param x x-coordinate point
     * \param y y-coordinate point
     * \param s value \f$ s \f$
     * \param t value \f$ t \f$
     */
    void
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const;

    /*! \brief Compute rational B-spline coefficients for a circle arc
     *
     * \param knots  knots of the B-spline0
     * \param Poly   polygon of the B-spline
     * \return       3 up to 9 the number of polygon points
    \*/
    int_type
    toNURBS( real_type knots[], real_type Poly[], bool get_size ) const;
    // Poly 3 x n matrix

    void
    info( ostream_type & stream ) const
    { stream << "CircleArc\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, CircleArc const & c );

    friend class ClothoidCurve;

  };

}

#endif

///
/// eof: Circle.hh
///
