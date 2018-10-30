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
/// file: ClothoidList.hh
///

#ifndef CLOTHOID_LIST_HH
#define CLOTHOID_LIST_HH

#include "Clothoid.hh"

//! Clothoid computations routine
namespace G2lib {

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/

  // Clothoid-clothoid
  class G2solve2arc {

    real_type tolerance;
    int_type  maxIter;

    real_type x0;
    real_type y0;
    real_type theta0;
    real_type kappa0;

    real_type x1;
    real_type y1;
    real_type theta1;
    real_type kappa1;

    // standard problem
    real_type lambda, phi, xbar, ybar;
    real_type th0, th1;
    real_type k0, k1;
    real_type DeltaK;
    real_type DeltaTheta;

    ClothoidCurve S0, S1;

    void
    evalA( real_type   alpha,
           real_type   L,
           real_type & A ) const;

    void
    evalA( real_type   alpha,
           real_type   L,
           real_type & A,
           real_type & A_1,
           real_type & A_2 ) const;

    void
    evalG( real_type alpha,
           real_type L,
           real_type th,
           real_type k,
           real_type G[2] ) const;

    void
    evalG( real_type alpha,
           real_type L,
           real_type th,
           real_type k,
           real_type G[2],
           real_type G_1[2],
           real_type G_2[2] ) const;

    void
    evalF( real_type const vars[2], real_type F[2] ) const;

    void
    evalFJ( real_type const vars[2],
            real_type       F[2],
            real_type       J[2][2] ) const;

    void
    buildSolution( real_type alpha, real_type L );

  public:

    G2solve2arc()
    : tolerance(1e-10)
    , maxIter(20)
    , x0(0)
    , y0(0)
    , theta0(0)
    , kappa0(0)
    , x1(0)
    , y1(0)
    , theta1(0)
    , kappa1(0)
    , lambda(0)
    , phi(0)
    , xbar(0)
    , ybar(0)
    , th0(0)
    , th1(0)
    , k0(0)
    , k1(0)
    {}

    ~G2solve2arc() {}

    int
    build( real_type x0, real_type y0, real_type theta0, real_type kappa0,
           real_type x1, real_type y1, real_type theta1, real_type kappa1 );

    void
    setTolerance( real_type tol );

    void
    setMaxIter( int tol );

    int
    solve();

    ClothoidCurve const & getS0() const { return S0; }
    ClothoidCurve const & getS1() const { return S1; }

  };

  /*\
   |    ____ ____            _            ____ _     ____
   |   / ___|___ \ ___  ___ | |_   _____ / ___| |   / ___|
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |   | |  | |
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/ |___| |__| |___
   |   \____|_____|___/\___/|_| \_/ \___|\____|_____\____|
  \*/

  // Clothoid-line-clothoid
  class G2solveCLC {

    real_type tolerance;
    int       maxIter;

    real_type x0;
    real_type y0;
    real_type theta0;
    real_type kappa0;
    real_type x1;
    real_type y1;
    real_type theta1;
    real_type kappa1;

    // standard problem
    real_type lambda, phi, xbar, ybar;
    real_type th0, th1;
    real_type k0, k1;

    ClothoidCurve S0, SM, S1;

    bool
    buildSolution( real_type sM, real_type thM );

  public:

    G2solveCLC()
    : tolerance(1e-10)
    , maxIter(20)
    , x0(0)
    , y0(0)
    , theta0(0)
    , kappa0(0)
    , x1(0)
    , y1(0)
    , theta1(0)
    , kappa1(0)
    , lambda(0)
    , phi(0)
    , xbar(0)
    , ybar(0)
    , th0(0)
    , th1(0)
    , k0(0)
    , k1(0)
    {}

    ~G2solveCLC() {}

    int
    build( real_type x0, real_type y0, real_type theta0, real_type kappa0,
           real_type x1, real_type y1, real_type theta1, real_type kappa1 );

    void
    setTolerance( real_type tol );

    void
    setMaxIter( int tol );

    int
    solve();

    ClothoidCurve const & getS0() const { return S0; }
    ClothoidCurve const & getSM() const { return SM; }
    ClothoidCurve const & getS1() const { return S1; }

  };

  /*\
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/
  // Clothoid-clothoid-clothoid with G2 continuity
  class G2solve3arc {

    ClothoidCurve S0, SM, S1;

    real_type tolerance;
    int       maxIter;

    // G2 interpolation data
    real_type x0;
    real_type y0;
    real_type theta0;
    real_type kappa0;
    real_type x1;
    real_type y1;
    real_type theta1;
    real_type kappa1;

    // standard scaled problem
    real_type phi, Lscale;
    real_type th0, th1;
    real_type s0, s1;

    // precomputed values
    real_type K0, K1, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14;

    void
    evalFJ( real_type const vars[2],
            real_type       F[2],
            real_type       J[2][2] ) const;

    void
    evalF( real_type const vars[2], real_type F[2] ) const;

    void
    buildSolution( real_type sM, real_type thM );

    int
    solve( real_type sM_guess, real_type thM_guess );

  public:

    G2solve3arc()
    : tolerance(1e-10)
    , maxIter(100)
    {}

    ~G2solve3arc() {}

    void setTolerance( real_type tol );
    void setMaxIter( int miter );

    /*!
     | Compute the 3 arc clothoid spline that fit the data
     |
     | \param[in] x0      initial `x` position
     | \param[in] y0      initial `y` position
     | \param[in] theta0  initial angle
     | \param[in] kappa0  initial curvature
     | \param[in] x1      final `x` position
     | \param[in] y1      final `y` position
     | \param[in] theta1  final angle
     | \param[in] kappa1  final curvature
     | \param[in] Dmax    rough desidered maximum angle variation, if 0 computed automatically
     | \param[in] dmax    rough desidered maximum angle divergence from guess, if 0 computed automatically
     | \return number of iteration, -1 if fails
     |
    \*/
    int
    build( real_type x0,
           real_type y0,
           real_type theta0,
           real_type kappa0,
           real_type x1,
           real_type y1,
           real_type theta1,
           real_type kappa1,
           real_type Dmax = 0,
           real_type dmax = 0 );

    /*!
     | Compute the 3 arc clothoid spline that fit the data
     |
     | \param[in] s0      length of the first segment
     | \param[in] x0      initial `x` position
     | \param[in] y0      initial `y` position
     | \param[in] theta0  initial angle
     | \param[in] kappa0  initial curvature
     | \param[in] s1      length of the last segment
     | \param[in] x1      final `x` position
     | \param[in] y1      final `y` position
     | \param[in] theta1  final angle
     | \param[in] kappa1  final curvature
     | \return number of iteration, -1 if fails
     |
    \*/
    int
    build_fixed_length( real_type s0,
                        real_type x0,
                        real_type y0,
                        real_type theta0,
                        real_type kappa0,
                        real_type s1,
                        real_type x1,
                        real_type y1,
                        real_type theta1,
                        real_type kappa1 );

    /*!
     | \return get the first clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getS0() const { return S0; }
    /*!
     | \return get the last clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getS1() const { return S1; }
    /*!
     | \return get the middle clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getSM() const { return SM; }

    /*!
     | \return get the length of the 3 arc G2 fitting
    \*/
    real_type
    totalLength() const {
      return S0.length() + S1.length() + SM.length();
    }

    /*!
     | \return get the total angle variation of the 3 arc G2 fitting
    \*/
    real_type
    thetaTotalVariation() const {
      return S0.thetaTotalVariation() +
             S1.thetaTotalVariation() +
             SM.thetaTotalVariation();
    }

    /*!
     | \return get the total curvature variation of the 3 arc G2 fitting
    \*/
    real_type
    curvatureTotalVariation() const {
      return S0.curvatureTotalVariation() +
             S1.curvatureTotalVariation() +
             SM.curvatureTotalVariation();
    }

    /*!
     | \return get the integral of the curvature squared of the 3 arc G2 fitting
    \*/
    real_type
    integralCurvature2() const {
      return S0.integralCurvature2() +
             S1.integralCurvature2() +
             SM.integralCurvature2();
    }

    /*!
     | \return get the integral of the jerk squared of the 3 arc G2 fitting
    \*/
    real_type
    integralJerk2() const {
      return S0.integralJerk2() +
             S1.integralJerk2() +
             SM.integralJerk2();
    }

    /*!
     | \return get the integral of the snap squared of the 3 arc G2 fitting
    \*/
    real_type
    integralSnap2() const {
      return S0.integralSnap2() +
             S1.integralSnap2() +
             SM.integralSnap2();
    }

    /*!
     | \param[out] thMin minimum angle in the 3 arc G2 fitting curve
     | \param[out] thMax maximum angle in the 3 arc G2 fitting curve
     | \return the difference of `thMax` and `thMin`
    \*/
    real_type
    thetaMinMax( real_type & thMin, real_type & thMax ) const;

    /*!
     | \return the difference of maximum-minimum angle in the 3 arc G2 fitting curve
    \*/
    real_type
    deltaTheta() const
    { real_type thMin, thMax; return thetaMinMax( thMin, thMax ); }

    /*!
     | \param[out] kMin minimum curvature in the 3 arc G2 fitting curve
     | \param[out] kMax maximum curvature in the 3 arc G2 fitting curve
     | \return the difference of `kMax` and `kMin`
    \*/
    real_type
    curvatureMinMax( real_type & kMin, real_type & kMax ) const;

    /*!
     | \return angle as a function of curvilinear coordinate
    \*/
    real_type theta( real_type s ) const;

    /*!
     | \return angle derivative (curvature) as a function of curvilinear coordinate
    \*/
    real_type theta_D( real_type s ) const;

    /*!
     | \return angle second derivative (curvature derivative) as a function of curvilinear coordinate
    \*/
    real_type theta_DD( real_type s ) const;

    /*!
     | \return angle third derivative as a function of curvilinear coordinate
    \*/
    real_type theta_DDD( real_type s ) const;

    /*!
     | \return x coordinate of the3 arc clothoid as a function of curvilinear coordinate
    \*/
    real_type X( real_type s ) const;

    /*!
     | \return y coordinate of the3 arc clothoid as a function of curvilinear coordinate
    \*/
    real_type Y( real_type s ) const;

    /*!
     | \return initial x coordinate of the 3 arc clothoid
    \*/
    real_type xBegin() const { return S0.xBegin(); }

    /*!
     | \return initial y coordinate of the 3 arc clothoid
    \*/
    real_type yBegin() const { return S0.yBegin(); }

    /*!
     | \return initial curvature of the 3 arc clothoid
    \*/
    real_type kappaBegin() const { return S0.kappaBegin(); }

    /*!
     | \return initial angle of the 3 arc clothoid
    \*/
    real_type thetaBegin() const { return S0.thetaBegin(); }

    /*!
     | \return final x coordinate of the 3 arc clothoid
    \*/
    real_type xEnd()const { return S1.xEnd(); }

    /*!
     | \return final y coordinate of the 3 arc clothoid
    \*/
    real_type yEnd() const { return S1.yEnd(); }

    /*!
     | \return final curvature of the 3 arc clothoid
    \*/
    real_type kappaEnd() const { return S1.kappaEnd(); }

    /*!
     | \return final angle of the 3 arc clothoid
    \*/
    real_type thetaEnd() const { return S1.thetaEnd(); }

    /*!
     | Compute parameters of 3 arc clothoid at curvilinear coordinate `s`
     |
     | \param[in]  s     curvilinear coordinate of where curve is computed
     | \param[out] theta the curve angle
     | \param[out] kappa the curve curvature
     | \param[out] x     the curve x-coordinate
     | \param[out] y     the curve y-coordinate
    \*/
    void
    eval( real_type   s,
          real_type & theta,
          real_type & kappa,
          real_type & x,
          real_type & y ) const;

    void eval( real_type s, real_type & x, real_type & y ) const;
    void eval_D( real_type s, real_type & x_D, real_type & y_D ) const;
    void eval_DD( real_type s, real_type & x_DD, real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type & x_DDD, real_type & y_DDD ) const;

    // offset curve
    void eval( real_type s, real_type offs, real_type & x, real_type & y ) const;
    void eval_D( real_type s, real_type offs, real_type & x_D, real_type & y_D ) const;
    void eval_DD( real_type s, real_type offs, real_type & x_DD, real_type & y_DD ) const;
    void eval_DDD( real_type s, real_type offs, real_type & x_DDD, real_type & y_DDD ) const;

    void
    rotate( real_type angle, real_type cx, real_type cy ) {
      S0.rotate( angle, cx, cy );
      S1.rotate( angle, cx, cy );
      SM.rotate( angle, cx, cy );
    }

    void
    translate( real_type tx, real_type ty ){
      S0.translate( tx, ty );
      S1.translate( tx, ty );
      SM.translate( tx, ty );
    }

    void
    reverse() {
      std::swap( S0, S1 );
      S0.reverse();
      S1.reverse();
      SM.reverse();
    }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidCurve const & c );

  };

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/
  //! \brief Class to manage a list Clothoid Curve (not necessarily G2 or G1 connected)
  class ClothoidList {

    std::vector<real_type>     s0;
    std::vector<ClothoidCurve> clotoidList;
    mutable int_type           last_idx;

  public:

    ClothoidList() : last_idx(0) {}
    ~ClothoidList();

    ClothoidList( ClothoidList const & s ) { copy(s); }

    ClothoidList const & operator = ( ClothoidList const & s )
    { copy(s); return *this; }

    void
    init() {
      s0.clear();
      clotoidList.clear();
      last_idx = 0;
    }

    void reserve( int_type n );
    void copy( ClothoidList const & L );

    void push_back( ClothoidCurve const & c );
    void push_back( real_type kappa0, real_type dkappa, real_type L );
    void push_back( real_type x0,     real_type y0,     real_type theta0,
                    real_type kappa0, real_type dkappa, real_type L );

    void push_back_G1( real_type x1, real_type y1, real_type theta1 );
    void push_back_G1( real_type x0, real_type y0, real_type theta0,
                       real_type x1, real_type y1, real_type theta1 );

    bool
    build_G1( int_type        n,
              real_type const x[],
              real_type const y[] );

    bool
    build_G1( int_type        n,
              real_type const x[],
              real_type const y[],
              real_type const theta[] );

    bool
    build_theta( int_type        n,
                 real_type const x[],
                 real_type const y[],
                 real_type       theta[] ) const;

    ClothoidCurve const & get( int_type idx ) const;
    ClothoidCurve const & getAtS( real_type s ) const;

    int_type numSegment() const { return int_type(clotoidList.size()); }

    bool findAtS( real_type s ) const;

    real_type theta( real_type s ) const;
    real_type theta_D( real_type s ) const;
    real_type theta_DD( real_type s ) const;
    real_type theta_DDD( real_type, int_type & ) const { return 0; }

    real_type totalLength() const {
      if ( s0.empty() ) return 0;
      return s0.back() - s0.front();
    }

    real_type X( real_type s ) const;
    real_type Y( real_type s ) const;

    real_type sBegin() const { return s0.front(); }
    real_type sEnd()   const { return s0.back(); }

    real_type xBegin() const { return clotoidList.front().xBegin(); }
    real_type xEnd()   const { return clotoidList.back().xEnd(); }

    real_type yBegin() const { return clotoidList.front().yBegin(); }
    real_type yEnd()   const { return clotoidList.back().yEnd(); }

    real_type thetaBegin() const { return clotoidList.front().thetaBegin(); }
    real_type thetaEnd()   const { return clotoidList.back().thetaEnd(); }

    real_type kappaBegin() const { return clotoidList.front().kappaBegin(); }
    real_type kappaEnd()   const { return clotoidList.back().kappaEnd(); }

    real_type length( int_type idx ) const { return s0[unsigned(idx+1)] - s0[unsigned(idx)]; }

    real_type sBegin( int_type idx ) const { return s0[unsigned(idx)]; }
    real_type sEnd  ( int_type idx ) const { return s0[unsigned(idx+1)]; }

    real_type xBegin( int_type idx ) const { return clotoidList[unsigned(idx)].xBegin(); }
    real_type xEnd  ( int_type idx ) const { return clotoidList[unsigned(idx)].xEnd(); }

    real_type yBegin( int_type idx ) const { return clotoidList[unsigned(idx)].yBegin(); }
    real_type yEnd  ( int_type idx ) const { return clotoidList[unsigned(idx)].yEnd(); }

    real_type thetaBegin( int_type idx ) const { return clotoidList[unsigned(idx)].thetaBegin(); }
    real_type thetaEnd  ( int_type idx ) const { return clotoidList[unsigned(idx)].thetaEnd(); }

    real_type kappaBegin( int_type idx ) const { return clotoidList[unsigned(idx)].kappaBegin(); }
    real_type kappaEnd  ( int_type idx ) const { return clotoidList[unsigned(idx)].kappaEnd(); }

    void
    eval( real_type   s,
          real_type & theta,
          real_type & kappa,
          real_type & x,
          real_type & y ) const;

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

    // offset curve
    void
    eval( real_type   s,
          real_type   offs,
          real_type & x,
          real_type & y ) const;
    void
    eval_D( real_type   s,
            real_type   offs,
            real_type & x_D,
            real_type & y_D ) const;
    void
    eval_DD( real_type   s,
             real_type   offs,
             real_type & x_DD,
             real_type & y_DD ) const;
    void
    eval_DDD( real_type   s,
              real_type   offs,
              real_type & x_DDD,
              real_type & y_DDD ) const;

    void
    getSTK( real_type s[],
            real_type theta[],
            real_type kappa[] ) const;

    void
    getXY( real_type x[], real_type y[] ) const;

    void
    getDeltaTheta( real_type deltaTheta[] ) const;

    void
    getDeltaKappa( real_type deltaKappa[] ) const;

    real_type
    closestPoint( real_type   qx,
                  real_type   qy,
                  real_type & X,
                  real_type & Y,
                  real_type & S ) const;

    real_type
    distance( real_type qx, real_type qy, real_type & S ) const
    { real_type X, Y; return closestPoint( qx, qy, X, Y, S ); }

    real_type
    distance( real_type qx, real_type qy ) const
    { real_type S; return distance( qx, qy, S ); }

    /*!
     | \brief Find parametric coordinate.
     |
     | \param  x    x-coordinate point
     | \param  y    y-coordinate point
     | \param  s    value \f$ s \f$
     | \param  t    value \f$ t \f$
     | \return idx  the segment with point at minimal distance, otherwise
     |              -(idx+1) if (x,y) cannot be projected orthogonally on the segment
     |
    \*/
    int_type
    findST( real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const;

    /*!
     | \brief Find parametric coordinate.
     |
     | \param  ibegin initial segment to compute the distance
     | \param  iend   final segment to compute the distance
     | \param  x      x-coordinate point
     | \param  y      y-coordinate point
     | \param  s      value \f$ s \f$
     | \param  t      value \f$ t \f$
     | \return idx    the segment with point at minimal distance, otherwise
     |                -(idx+1) if (x,y) cannot be projected orthogonally on the segment
    \*/
    int_type
    findST( int_type    ibegin,
            int_type    iend,
            real_type   x,
            real_type   y,
            real_type & s,
            real_type & t ) const;

    void rotate( real_type angle, real_type cx, real_type cy );
    void translate( real_type tx, real_type ty );
    void changeOrigin( real_type newx0, real_type newy0 );
    void scale( real_type sfactor );
    void reverse();

    /*! \brief split clothois in smaller segments
     |
     | \param split_angle maximum angle variation
     | \param split_size  maximum height of the triangle
     | \param split_offs  curve offset
     | \param bb          splitting data structures vector
    \*/
    void
    bbSplit( real_type split_angle,
             real_type split_size,
             real_type split_offs,
             std::vector<ClothoidCurve::bbData> & bb ) const {
      bb.clear();
      std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
      for (; ic != clotoidList.end(); ++ic )
        ic->bbSplit( split_angle, split_size, split_offs, bb, false );
    }
    /*\
     |  _     _                      _
     | (_)_ _| |_ ___ _ _ ___ ___ __| |_
     | | | ' \  _/ -_) '_(_-</ -_) _|  _|
     | |_|_||_\__\___|_| /__/\___\__|\__|
     |
    \*/
    /*! \brief intersect two clothoid arcs
     *
     * \param offs      offset of the first arc
     * \param c         the second clothoid arc
     * \param c_offs    offset of the second arc
     * \param s1        intersection parameters of the first arc
     * \param s2        intersection parameters of the second arc
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( real_type                offs,
               ClothoidCurve const &    c,
               real_type                c_offs,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const;

    /*! \brief intersect two clothoid curve
     *
     * \param c         the clothoid arc
     * \param s1        intersection parameters of the arc
     * \param s2        intersection parameters of the clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( ClothoidCurve const    & c,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      intersect( 0, c, 0, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect a clothoid curve to a circle arc
     *
     * \param c_in      the circle arc
     * \param s1        intersection parameters of the arc
     * \param s2        intersection parameters of the clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( CircleArc const        & c_in,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      ClothoidCurve c(c_in);
      intersect( 0, c, 0, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect a clothoid curve to a circle arc
     *
     * \param offs      offset of the clothoid
     * \param c_in      the circle arc
     * \param c_offs    offset of the circle
     * \param s1        intersection parameters of the arc
     * \param s2        intersection parameters of the clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( real_type                offs,
               CircleArc const        & c_in,
               real_type                c_offs,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      ClothoidCurve c(c_in);
      intersect( offs, c, c_offs, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect a clothoid curve to a line segment
     *
     * \param c_in      the line segment
     * \param s1        intersection parameters of the arc
     * \param s2        intersection parameters of the clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( LineSegment const      & c_in,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      ClothoidCurve c(c_in);
      intersect( 0, c, 0, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect a clothoid curve to a line segment
     *
     * \param offs      offset of the clothoid
     * \param c_in      the line segment
     * \param c_offs    offset of the line segment
     * \param s1        intersection parameters of the arc
     * \param s2        intersection parameters of the clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( real_type                offs,
               LineSegment const      & c_in,
               real_type                c_offs,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      ClothoidCurve c(c_in);
      intersect( offs, c, c_offs, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect two clothoid list
     *
     * \param CL        the second clothoid list
     * \param s1        intersection parameters of the first clothoid list
     * \param s2        intersection parameters of the second clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( ClothoidList const     & CL,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const {
      intersect( 0, CL, 0, s1, s2, max_iter, tolerance );
    }

    /*! \brief intersect two clothoid list
     *
     * \param offs      offset of the first clothoid list
     * \param CL        the second clothoid list
     * \param c_offs    offset of the second clothoid list
     * \param s1        intersection parameters of the first clothoid list
     * \param s2        intersection parameters of the second clothoid list
     * \param max_iter  max allowed iteration
     * \param tolerance admitted tolerance
     */
    void
    intersect( real_type                offs,
               ClothoidList const &     CL,
               real_type                c_offs,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2,
               int_type                 max_iter,
               real_type                tolerance ) const;

    /*! \brief Save Clothoid list to a stream
     *
     * \param stream stream to save
     */
    void
    export_table( ostream_type & stream ) const;

    /*! \brief Save Clothoid list to a stream
     *
     * \param stream streamstream to save
     */
    void
    export_ruby( ostream_type & stream ) const;

  };

  /*\
   |
   |    ___ _     _   _        _    _ ___      _ _           ___ ___
   |   / __| |___| |_| |_  ___(_)__| / __|_ __| (_)_ _  ___ / __|_  )
   |  | (__| / _ \  _| ' \/ _ \ / _` \__ \ '_ \ | | ' \/ -_) (_ |/ /
   |   \___|_\___/\__|_||_\___/_\__,_|___/ .__/_|_|_||_\___|\___/___|
   |                                     |_|
  \*/

  class ClothoidSplineG2 {
  public:
    typedef enum { P1 = 1, P2, P3, P4, P5, P6, P7, P8, P9 } TargetType;

  private:

    std::vector<real_type> x;
    std::vector<real_type> y;
    TargetType             tt;
    real_type              theta_I;
    real_type              theta_F;
    int_type               npts;

    // work vector
    mutable std::vector<real_type> k, dk, L, kL, L_1, L_2, k_1, k_2, dk_1, dk_2;

    real_type
    diff2pi( real_type in ) const {
      return in-m_2pi*round(in/m_2pi);
    }

  public:

    ClothoidSplineG2() : tt(P1) {}
    ~ClothoidSplineG2() {}

    void
    setP1( real_type theta0, real_type thetaN )
    { tt = P1; theta_I = theta0; theta_F = thetaN; }

    void setP2() { tt = P2; }
    void setP3() { tt = P3; }
    void setP4() { tt = P4; }
    void setP5() { tt = P5; }
    void setP6() { tt = P6; }
    void setP7() { tt = P7; }
    void setP8() { tt = P8; }
    void setP9() { tt = P9; }

    void
    build( real_type const xvec[],
           real_type const yvec[],
           int_type        npts );

    int_type numPnts() const { return npts; }
    int_type numTheta() const;
    int_type numConstraints() const;

    void
    guess( real_type theta_guess[],
           real_type theta_min[],
           real_type theta_max[] ) const;

    bool
    objective( real_type const theta[], real_type & f ) const;

    bool
    gradient( real_type const theta[], real_type g[] ) const;

    bool
    constraints( real_type const theta[], real_type c[] ) const;

    int_type
    jacobian_nnz() const;

    bool
    jacobian_pattern( int_type i[], int_type j[] ) const;

    bool
    jacobian_pattern_matlab( real_type i[], real_type j[] ) const;

    bool
    jacobian( real_type const theta[], real_type vals[] ) const;

    void
    info( ostream_type & stream ) const
    { stream << "ClothoidSplineG2\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidSplineG2 const & c );

  };

}

#endif

///
/// eof: ClothoidList.hh
///
