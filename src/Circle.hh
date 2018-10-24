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
#include "Triangle2D.hh"

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
  class CircleArc { // : public BaseCurve {

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


    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#if 0
    virtual
    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const G2LIB_OVERRIDE;

    virtual
    void
    bbox( real_type   offs,
          real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const G2LIB_OVERRIDE;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    length() const G2LIB_OVERRIDE
    { return L; }

    virtual
    real_type
    length( real_type offs ) const G2LIB_OVERRIDE

    virtual
    real_type
    xBegin() const G2LIB_OVERRIDE
    { return x0; }

    virtual
    real_type
    yBegin() const G2LIB_OVERRIDE
    { return y0; }

    virtual
    real_type
    xEnd() const G2LIB_OVERRIDE
    { return X(L); }

    virtual
    real_type
    yEnd() const G2LIB_OVERRIDE
    { return Y(L); }

    virtual
    real_type
    xBegin( real_type offs ) const G2LIB_OVERRIDE
    { return X(0,offs); }

    virtual
    real_type
    yBegin( real_type offs ) const G2LIB_OVERRIDE
    { return Y(0,offs); }

    virtual
    real_type
    xEnd( real_type offs ) const G2LIB_OVERRIDE
    { return X(L,offs); }

    virtual
    real_type
    yEnd( real_type offs ) const G2LIB_OVERRIDE
    { return Y(L,offs); }

    virtual
    real_type
    tx_Begin() const G2LIB_OVERRIDE
    { return tx(0); }

    virtual
    real_type
    ty_Begin() const G2LIB_OVERRIDE
    { return ty(0); }

    virtual
    real_type
    tx_End() const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    ty_End()   const G2LIB_OVERRIDE
    { return s0; }

    virtual
    real_type
    nx_Begin() const G2LIB_OVERRIDE
    { return -s0; }

    virtual
    real_type
    ny_Begin() const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    nx_End() const G2LIB_OVERRIDE
    { return -s0; }

    virtual
    real_type
    ny_End()   const G2LIB_OVERRIDE
    { return c0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    X( real_type s ) const G2LIB_OVERRIDE
    { return x0+s*c0; }

    virtual
    real_type
    Y( real_type s ) const G2LIB_OVERRIDE
    { return y0+s*s0; }

    virtual
    real_type
    X_D( real_type ) const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    Y_D( real_type ) const G2LIB_OVERRIDE
    { return s0; }

    virtual
    real_type
    X_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    X_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const G2LIB_OVERRIDE {
      x = x0+s*c0;
      y = y0+s*s0;
    }

    virtual
    void
    eval_D( real_type,
            real_type & x_D,
            real_type & y_D ) const G2LIB_OVERRIDE {
      x_D = c0;
      y_D = s0;
    }

    virtual
    void
    eval_DD( real_type,
             real_type & x_DD,
             real_type & y_DD ) const G2LIB_OVERRIDE {
      x_DD = 0;
      y_DD = 0;
    }

    virtual
    void
    eval_DDD( real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const G2LIB_OVERRIDE {
      x_DDD = 0;
      y_DDD = 0;
    }

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    virtual
    real_type
    nx( real_type ) const G2LIB_OVERRIDE
    { return -s0; }

    virtual
    real_type
    ny( real_type ) const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    nx_D( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ny_D( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    nx_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ny_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    nx_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ny_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    real_type
    tx( real_type ) const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    ty( real_type ) const G2LIB_OVERRIDE
    { return s0; }

    virtual
    real_type
    tx_D( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ty_D( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    tx_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ty_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    tx_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    ty_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    nor( real_type, real_type n[2] ) const G2LIB_OVERRIDE
    { n[0] = -s0; n[1] = c0; }

    virtual
    void
    nor_D( real_type, real_type n_D[2] ) const G2LIB_OVERRIDE
    { n_D[0] = n_D[1] = 0; }

    virtual
    void
    nor_DD( real_type, real_type n_DD[2] ) const G2LIB_OVERRIDE
    { n_DD[0] = n_DD[1] = 0; }

    virtual
    void
    nor_DDD( real_type, real_type n_DDD[2] ) const G2LIB_OVERRIDE
    { n_DDD[0] = n_DDD[1] = 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    tan( real_type, real_type t[2] ) const G2LIB_OVERRIDE
    { t[0] = c0; t[1] = s0; }

    virtual
    void
    tan_D( real_type, real_type t_D[2] ) const G2LIB_OVERRIDE
    { t_D[0] = t_D[1] = 0; }

    virtual
    void
    tan_DD( real_type, real_type t_DD[2] ) const G2LIB_OVERRIDE
    { t_DD[0] = t_DD[1] = 0; }

    virtual
    void
    tan_DDD( real_type, real_type t_DDD[2] ) const G2LIB_OVERRIDE
    { t_DDD[0] = t_DDD[1] = 0; }

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    virtual
    real_type
    X( real_type s, real_type offs ) const G2LIB_OVERRIDE
    { return x0 + s*c0 - offs*s0; }

    virtual
    real_type
    Y( real_type s, real_type offs ) const G2LIB_OVERRIDE
    { return y0 + s*s0 + offs*c0; }

    virtual
    real_type
    X_D( real_type, real_type ) const G2LIB_OVERRIDE
    { return c0; }

    virtual
    real_type
    Y_D( real_type, real_type ) const G2LIB_OVERRIDE
    { return s0; }

    virtual
    real_type
    X_DD( real_type, real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y_DD( real_type, real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    X_DDD( real_type, real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y_DDD( real_type, real_type ) const G2LIB_OVERRIDE
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    virtual
    void
    eval( real_type   s,
          real_type   offs,
          real_type & x,
          real_type & y ) const G2LIB_OVERRIDE {
      x = x0 + s*c0 - offs*s0;
      y = y0 + s*s0 + offs*c0;
    }

    virtual
    void
    eval_D( real_type,
            real_type,
            real_type & x_D,
            real_type & y_D ) const G2LIB_OVERRIDE {
      x_D = c0;
      y_D = s0;
    }

    virtual
    void
    eval_DD( real_type,
             real_type,
             real_type & x_DD,
             real_type & y_DD ) const G2LIB_OVERRIDE {
      x_DD = y_DD = 0;
    }

    virtual
    void
    eval_DDD( real_type,
              real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const G2LIB_OVERRIDE {
      x_DDD = y_DDD = 0;
    }

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    virtual
    bool
    collision( BaseCurve const & ) const G2LIB_OVERRIDE;

    virtual
    bool
    collision( real_type         offs,
               BaseCurve const & obj,
               real_type         offs_obj ) const G2LIB_OVERRIDE;

    virtual
    void
    intersect( BaseCurve const & obj,
               IntersectList   & ilist ) const G2LIB_OVERRIDE;

    virtual
    void
    intersect( real_type         offs,
               BaseCurve const & obj,
               real_type         offs_obj,
               IntersectList   & ilist ) const G2LIB_OVERRIDE;
    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    /*!
     * \brief compute the point at minimum distance from a point `[x,y]` and the line segment
     *
     * \param qx x-coordinate
     * \param qy y-coordinate
     * \param x  x-coordinate of the closest point
     * \param y  y-coordinate of the closest point
     * \param s  param of the closest point
     * \return the distance point-segment
    \*/

    virtual
    real_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type & x,
                  real_type & y,
                  real_type & s ) const G2LIB_OVERRIDE;

    virtual
    real_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type   offs,
                  real_type & x,
                  real_type & y,
                  real_type & s ) const G2LIB_OVERRIDE;

    /*!
     | \param  qx  x-coordinate of the point
     | \param  qy  y-coordinate of the point
     | \param  x   x-coordinate of the projected point on the curve
     | \param  y   y-coordinate of the projected point on the curve
     | \param  s   parameter on the curve of the projection
     | \return 1  = unique orthogonal projection
     |         0  = more than one projection (first returned)
     |         -1 = projection line not othogonal to curve
     |         -2 = projection line not othogonal andnot unique
    \*/
    virtual
    int_type
    projection( real_type   qx,
                real_type   qy,
                real_type & x,
                real_type & y,
                real_type & s ) const G2LIB_OVERRIDE;

    /*!
     | \param  qx   x-coordinate of the point
     | \param  qy   y-coordinate of the point
     | \param  offs offset of the curve
     | \param  x    x-coordinate of the projected point on the curve
     | \param  y    y-coordinate of the projected point on the curve
     | \param  s    parameter on teh curve of the projection
     | \return 1  = unique orthogonal projection
     |         0  = more than one projection (first returned)
     |         -1 = projection line not othogonal to curve
     |         -2 = projection line not othogonal andnot unique
    \*/
    virtual
    int_type // true if projection is unique and orthogonal
    projection( real_type   qx,
                real_type   qy,
                real_type   offs,
                real_type & x,
                real_type & y,
                real_type & s ) const G2LIB_OVERRIDE;

    /*\
     |    __ _           _ ____ _____
     |   / _(_)_ __   __| / ___|_   _|
     |  | |_| | '_ \ / _` \___ \ | |
     |  |  _| | | | | (_| |___) || |
     |  |_| |_|_| |_|\__,_|____/ |_|
    \*/

    /*! \brief Find parametric coordinate.
     *
     * We consider the line passing to the point \f$ P \f$
     * with tangent \f$ T \f$ and a point \f$ Q \f$
     * compute the coordinte \f$ s \f$ and \f$ t \f$ such that
     * \f$ Q = P + T s + N t \f$
     * where \f$ P + T s \f$ is the point on the line at coordinate
     * \f$ s \f$ and \f$ N \f$ is the normal to the line obtained by
     * rotating by `90` degree counterclockwise the tangent \f$ T \f$.
     *
     * \param x x-coordinate point
     * \param y y-coordinate point
     * \param s value \f$ s \f$
     * \param t value \f$ t \f$
     */

    virtual
    bool
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const G2LIB_OVERRIDE {
      real_type dx = x - x0;
      real_type dy = y - y0;
      s = c0 * dx + s0 * dy;
      t = c0 * dy - s0 * dx;
      return true;
    }

    /*\
     |   _   _ _   _ ____  ____ ____
     |  | \ | | | | |  _ \| __ ) ___|
     |  |  \| | | | | |_) |  _ \___ \
     |  | |\  | |_| |  _ <| |_) |__) |
     |  |_| \_|\___/|_| \_\____/____/
    \*/

    virtual
    void
    paramNURBS( int_type & n_knots,
                int_type & n_pnts ) const G2LIB_OVERRIDE;

    virtual
    void
    toNURBS( real_type knots[],
             real_type Poly[][3] ) const G2LIB_OVERRIDE;

    virtual
    void
    toBS( real_type knots[],
          real_type Poly[][2] ) const G2LIB_OVERRIDE;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

#endif










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
    delta_theta() const
    { return L*k; }

    real_type
    theta( real_type s ) const
    { return theta0 + s*k; }

    real_type
    theta_D( real_type ) const
    { return k; }

    real_type
    theta_DD( real_type ) const
    { return 0; }

    real_type
    theta_DDD( real_type ) const
    { return 0; }

    real_type
    totalLength() const
    { return L; }

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

    void
    XY( real_type s, real_type & x, real_type & y ) const;

    void
    XY( real_type s, real_type t, real_type & x, real_type & y ) const;

    void
    TG( real_type s, real_type & tx, real_type & ty ) const;

    void
    NOR( real_type s, real_type & nx, real_type & ny ) const;

    void
    NOR_D( real_type s, real_type & nx_D, real_type & ny_D ) const;

    void
    NOR_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const;

    void
    NOR_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const;

    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const;

    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const;

    void
    eval_DD( real_type   s,
             real_type & x_DD,
             real_type & y_DD ) const;

    void
    eval_DDD( real_type   s,
              real_type & x_DDD,
              real_type & y_DDD ) const;

    void
    eval( real_type   s,
          real_type   t,
          real_type & x,
          real_type & y ) const;

    void
    eval_D( real_type   s,
            real_type   t,
            real_type & x_D,
            real_type & y_D ) const;

    void
    eval_DD( real_type   s,
             real_type   t,
             real_type & x_DD,
             real_type & y_DD ) const;

    void
    eval_DDD( real_type   s,
              real_type   t,
              real_type & x_DDD,
              real_type & y_DDD ) const;

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

    void
    center( real_type & cx, real_type & cy ) const;

    real_type ray() const { return 1/std::abs(k); }

    //! get the bounding box triangle (if angle variation less that pi/3)
    bool
    bbTriangle( real_type & x0, real_type & y0,
                real_type & x1, real_type & y1,
                real_type & x2, real_type & y2 ) const;

    bool
    bbTriangle( real_type p0[2],
                real_type p1[2],
                real_type p2[2] ) const {
      return bbTriangle( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );
    }

    bool
    bbTriangle( Triangle2D & t ) const {
      real_type p0[2], p1[2], p2[2];
      bool ok = bbTriangle( p0, p1, p2 );
      if ( ok ) t.build( p0, p1, p2 );
      return ok;
    }

    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const;

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

    void
    paramNURBS( int_type & n_knots,
                int_type & n_pnts ) const ;

    /*!
     | \brief Compute rational B-spline coefficients for a circle arc
     |
     | \param knots  knots of the B-spline0
     | \param Poly   polygon of the B-spline
    \*/
    void
    toNURBS( real_type knots[], real_type Poly[] ) const;
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
