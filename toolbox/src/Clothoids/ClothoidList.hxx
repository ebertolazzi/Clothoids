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
/// file: ClothoidList.hxx
///

namespace G2lib {

  using std::vector;

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/


  //!
  //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
  //! two clothoids arc that solve the G2 problem
  //!
  //! \f[
  //! \begin{array}{ll}
  //! \textrm{endpoints:}\quad &
  //! \begin{cases}
  //!   G(0) = \mathbf{p}_0 & \\[0.5em]
  //!   G(L) = \mathbf{p}_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{angles:}\quad &
  //! \begin{cases}
  //!   \theta(0) = \theta_0 & \\[0.5em]
  //!   \theta(L) = \theta_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{curvature:}\quad &
  //! \begin{cases}
  //!   \kappa(0) = \kappa_0 & \\[0.5em]
  //!   \kappa(L) = \kappa_1 &
  //! \end{cases}
  //! \end{array}
  //! \f]
  //!
  //! **note**
  //!
  //! The solution do not exist for all the combination of points/angle/curvature
  //!
  //!
  class G2solve2arc {

    real_type tolerance{real_type(1e-10)};
    integer   maxIter{20};

    real_type x0{real_type(0)};
    real_type y0{real_type(0)};
    real_type theta0{real_type(0)};
    real_type kappa0{real_type(0)};

    real_type x1{real_type(0)};
    real_type y1{real_type(0)};
    real_type theta1{real_type(0)};
    real_type kappa1{real_type(0)};

    // standard problem
    real_type lambda{real_type(0)};
    real_type phi{real_type(0)};
    real_type xbar{real_type(0)};
    real_type ybar{real_type(0)};
    real_type th0{real_type(0)};
    real_type th1{real_type(0)};
    real_type k0{real_type(0)};
    real_type k1{real_type(0)};
    real_type DeltaK{real_type(0)};
    real_type DeltaTheta{real_type(0)};

    ClothoidCurve S0, S1;

    void
    evalA(
      real_type   alpha,
      real_type   L,
      real_type & A
    ) const;

    void
    evalA(
      real_type   alpha,
      real_type   L,
      real_type & A,
      real_type & A_1,
      real_type & A_2
    ) const;

    void
    evalG(
      real_type alpha,
      real_type L,
      real_type th,
      real_type k,
      real_type G[2]
    ) const;

    void
    evalG(
      real_type alpha,
      real_type L,
      real_type th,
      real_type k,
      real_type G[2],
      real_type G_1[2],
      real_type G_2[2]
    ) const;

    void
    evalF( real_type const vars[2], real_type F[2] ) const;

    void
    evalFJ(
      real_type const vars[2],
      real_type       F[2],
      real_type       J[2][2]
    ) const;

    void
    buildSolution( real_type alpha, real_type L );

  public:

    //!
    //! Build an empty clothoid list
    //!
    G2solve2arc() = default;
    ~G2solve2arc() = default;

    //!
    //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
    //! two clothoids arc that solve the G2 problem, with data
    //!
    //! \f[ \mathbf{p}_0 = (x_0,y_0)^T, \qquad \mathbf{p}_1 = (x_1,y_1)^T \f]
    //! \f[ \theta_0, \qquad \theta_1, \qquad \kappa_0, \qquad \kappa_1 \f]
    //!
    //! \param[in] x0     \f$ x_0      \f$
    //! \param[in] y0     \f$ y_0      \f$
    //! \param[in] theta0 \f$ \theta_0 \f$
    //! \param[in] kappa0 \f$ \kappa_0 \f$
    //! \param[in] x1     \f$ x_1      \f$
    //! \param[in] y1     \f$ y_1      \f$
    //! \param[in] theta1 \f$ \theta_1 \f$
    //! \param[in] kappa1 \f$ \kappa_1 \f$
    //! \return number of iterations of -1 if failed
    //!
    int
    build(
      real_type x0, real_type y0, real_type theta0, real_type kappa0,
      real_type x1, real_type y1, real_type theta1, real_type kappa1
    );

    //!
    //! Fix tolerance for the G2 problem
    //!
    void set_tolerance( real_type tol );
    void setTolerance( real_type tol ) { set_tolerance( tol ); }

    //!
    //! Fix maximum number of iteration for the G2 problem
    //!
    void set_max_iter( integer miter );
    void setMaxIter( integer miter ) { set_max_iter( miter ); }

    //!
    //! Solve the G2 problem
    //!
    //! \return number of iterations of -1 if failed
    //!
    int solve();

    //!
    //! Return the first clothoid of the G2 clothoid list
    //!
    ClothoidCurve const & getS0() const { return S0; }

    //!
    //! Return the second clothoid of the G2 clothoid list
    //!
    ClothoidCurve const & getS1() const { return S1; }

  };

  /*\
   |    ____ ____            _            ____ _     ____
   |   / ___|___ \ ___  ___ | |_   _____ / ___| |   / ___|
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |   | |  | |
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/ |___| |__| |___
   |   \____|_____|___/\___/|_| \_/ \___|\____|_____\____|
  \*/

  //!
  //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
  //! 2 clothoid and one line segment that solve the G2 problem
  //!
  //! \f[
  //! \begin{array}{ll}
  //! \textrm{endpoints:}\quad &
  //! \begin{cases}
  //!   G(0) = \mathbf{p}_0 & \\[0.5em]
  //!   G(L) = \mathbf{p}_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{angles:}\quad &
  //! \begin{cases}
  //!   \theta(0) = \theta_0 & \\[0.5em]
  //!   \theta(L) = \theta_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{curvature:}\quad &
  //! \begin{cases}
  //!   \kappa(0) = \kappa_0 & \\[0.5em]
  //!   \kappa(L)  = \kappa_1 &
  //! \end{cases}
  //! \end{array}
  //! \f]
  //!
  //! **note**
  //!
  //! The solution do not exist for all the combination of points/angle/curvature
  //!
  class G2solveCLC {

    real_type tolerance{real_type(1e-10)};
    int       maxIter{20};

    real_type x0{real_type(0)};
    real_type y0{real_type(0)};
    real_type theta0{real_type(0)};
    real_type kappa0{real_type(0)};
    real_type x1{real_type(0)};
    real_type y1{real_type(0)};
    real_type theta1{real_type(0)};
    real_type kappa1{real_type(0)};

    // standard problem
    real_type lambda{real_type(0)};
    real_type phi{real_type(0)};
    real_type xbar{real_type(0)};
    real_type ybar{real_type(0)};
    real_type th0{real_type(0)};
    real_type th1{real_type(0)};
    real_type k0{real_type(0)};
    real_type k1{real_type(0)};

    ClothoidCurve S0, SM, S1;

    bool
    buildSolution( real_type sM, real_type thM );

  public:

    //!
    //! Build an empty clothoid list
    //!
    G2solveCLC() = default;
    ~G2solveCLC() = default;

    //!
    //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
    //! two clothoids and one line segment that solve the G2 problem, with data
    //!
    //! \f[ \mathbf{p}_0 = (x_0,y_0)^T, \qquad \mathbf{p}_1 = (x_1,y_1)^T \f]
    //! \f[ \theta_0, \qquad \theta_1, \qquad \kappa_0, \qquad \kappa_1 \f]
    //!
    //! \param[in] x0     \f$ x_0      \f$
    //! \param[in] y0     \f$ y_0      \f$
    //! \param[in] theta0 \f$ \theta_0 \f$
    //! \param[in] kappa0 \f$ \kappa_0 \f$
    //! \param[in] x1     \f$ x_1      \f$
    //! \param[in] y1     \f$ y_1      \f$
    //! \param[in] theta1 \f$ \theta_1 \f$
    //! \param[in] kappa1 \f$ \kappa_1 \f$
    //! \return number of iterations of -1 if failed
    //!
    int
    build(
      real_type x0, real_type y0, real_type theta0, real_type kappa0,
      real_type x1, real_type y1, real_type theta1, real_type kappa1
    );

    //!
    //! Fix tolerance for the G2 problem
    //!
    void set_tolerance( real_type tol );
    void setTolerance( real_type tol ) { set_tolerance( tol ); }

    //!
    //! Fix maximum number of iteration for the G2 problem
    //!
    void set_max_iter( integer miter );
    void setMaxIter( integer miter ) { set_max_iter( miter ); }

    //!
    //! Solve the G2 problem
    //!
    //! \return number of iterations of -1 if failed
    //!
    int solve();

    //!
    //! Return the first clothoid of the G2 clothoid list
    //!
    ClothoidCurve const & getS0() const { return S0; }

    //!
    //! Return the second segment (the line) as a clothoid
    //!
    ClothoidCurve const & getSM() const { return SM; }

    //!
    //! Return the third clothoid of the G2 clothoid list
    //!
    ClothoidCurve const & getS1() const { return S1; }

    void save( ostream_type & stream ) const;
  };

  /*\
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/

  //!
  //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
  //! 3 clothoid and one line segment that solve the G2 problem
  //!
  //! **match**
  //!
  //! \f[
  //! \begin{array}{ll}
  //! \textrm{endpoints:}\quad&
  //! \begin{cases}
  //!   G(0) = \mathbf{p}_0 & \\[0.5em]
  //!   G(L) = \mathbf{p}_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{angles:}\quad&
  //! \begin{cases}
  //!   \theta(0) = \theta_0 & \\[0.5em]
  //!   \theta(L) = \theta_1 &
  //! \end{cases}
  //! \\[1em]
  //! \textrm{curvature:}\quad&
  //! \begin{cases}
  //!   \kappa(0) = \kappa_0 & \\[0.5em]
  //!   \kappa(L)  = \kappa_1 &
  //! \end{cases}
  //! \end{array}
  //! \f]
  //!
  //! **Reference**
  //!
  //! The solution algorithm is described in
  //!
  //! - **E.Bertolazzi, M.Frego**, On the G2 Hermite Interpolation Problem with clothoids
  //!   Journal of Computational and Applied Mathematics, vol 341, pp. 99-116, 2018
  //!
  //! \rst
  //!
  //!   .. image:: ../../images/G2problem3arc.jpg
  //!      :width: 80%
  //!      :align: center
  //!
  //! \endrst
  //!
  class G2solve3arc {

    ClothoidCurve S0, SM, S1;

    real_type tolerance{real_type(1e-10)};
    int       maxIter{100};

    // G2 interpolation data
    real_type x0{real_type(0)};
    real_type y0{real_type(0)};
    real_type theta0{real_type(0)};
    real_type kappa0{real_type(0)};
    real_type x1{real_type(0)};
    real_type y1{real_type(0)};
    real_type theta1{real_type(0)};
    real_type kappa1{real_type(0)};

    // standard scaled problem
    real_type phi{real_type(0)};
    real_type Lscale{real_type(0)};
    real_type th0{real_type(0)};
    real_type th1{real_type(0)};
    real_type s0{real_type(0)};
    real_type s1{real_type(0)};

    // precomputed values
    real_type K0{real_type(0)},
              K1{real_type(0)},
              c0{real_type(0)},
              c1{real_type(0)},
              c2{real_type(0)},
              c3{real_type(0)},
              c4{real_type(0)},
              c5{real_type(0)},
              c6{real_type(0)},
              c7{real_type(0)},
              c8{real_type(0)},
              c9{real_type(0)},
              c10{real_type(0)},
              c11{real_type(0)},
              c12{real_type(0)},
              c13{real_type(0)},
              c14{real_type(0)};

    void
    evalFJ(
      real_type const vars[2],
      real_type       F[2],
      real_type       J[2][2]
    ) const;

    void
    evalF( real_type const vars[2], real_type F[2] ) const;

    void
    buildSolution( real_type sM, real_type thM );

    int
    solve( real_type sM_guess, real_type thM_guess );

  public:

    G2solve3arc() = default;
    ~G2solve3arc() = default;

    //!
    //! Fix tolerance for the G2 problem
    //!
    void set_tolerance( real_type tol );
    void setTolerance( real_type tol ) { set_tolerance( tol ); }

    //!
    //! Fix maximum number of iteration for the G2 problem
    //!
    void set_max_iter( integer miter );
    void setMaxIter( integer miter ) { set_max_iter( miter ); }

    //!
    //! Compute the 3 arc clothoid spline that fit the data
    //!
    //! \param[in] x0      initial `x` position
    //! \param[in] y0      initial `y` position
    //! \param[in] theta0  initial angle
    //! \param[in] kappa0  initial curvature
    //! \param[in] x1      final `x` position
    //! \param[in] y1      final `y` position
    //! \param[in] theta1  final angle
    //! \param[in] kappa1  final curvature
    //! \param[in] Dmax    rough desidered maximum angle variation, if 0 computed automatically
    //! \param[in] dmax    rough desidered maximum angle divergence from guess, if 0 computed automatically
    //! \return number of iteration, -1 if fails
    //!
    //!
    int
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type kappa0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type kappa1,
      real_type Dmax = 0,
      real_type dmax = 0
    );

    //!
    //! Compute the 3 arc clothoid spline that fit the data
    //!
    //! \param[in] s0      length of the first segment
    //! \param[in] x0      initial `x` position
    //! \param[in] y0      initial `y` position
    //! \param[in] theta0  initial angle
    //! \param[in] kappa0  initial curvature
    //! \param[in] s1      length of the last segment
    //! \param[in] x1      final `x` position
    //! \param[in] y1      final `y` position
    //! \param[in] theta1  final angle
    //! \param[in] kappa1  final curvature
    //! \return number of iteration, -1 if fails
    //!
    //!
    int
    build_fixed_length(
      real_type s0,
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type kappa0,
      real_type s1,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type kappa1
    );

    //!
    //! \return get the first clothoid for the 3 arc G2 fitting
    //!
    ClothoidCurve const & getS0() const { return S0; }

    //!
    //! \return get the last clothoid for the 3 arc G2 fitting
    //!
    ClothoidCurve const & getS1() const { return S1; }

    //!
    //! \return get the middle clothoid for the 3 arc G2 fitting
    //!
    ClothoidCurve const & getSM() const { return SM; }

    //!
    //! \return get the length of the 3 arc G2 fitting
    //!
    real_type
    total_length() const {
      return S0.length() + S1.length() + SM.length();
    }

    real_type totalLength() const { return total_length(); }

    //!
    //! \return get the total angle variation of the 3 arc G2 fitting
    //!
    real_type
    theta_total_variation() const {
      return S0.theta_total_variation() +
             S1.theta_total_variation() +
             SM.theta_total_variation();
    }

    //!
    //! \return get the total curvature variation of the 3 arc G2 fitting
    //!
    real_type
    curvature_total_variation() const {
      return S0.curvature_total_variation() +
             S1.curvature_total_variation() +
             SM.curvature_total_variation();
    }

    //!
    //! \return get the integral of the curvature squared of the 3 arc G2 fitting
    //!
    real_type
    integral_curvature2() const {
      return S0.integral_curvature2() +
             S1.integral_curvature2() +
             SM.integral_curvature2();
    }

    real_type
    integralCurvature2() const
    { return integral_curvature2(); }

    //!
    //! \return get the integral of the jerk squared of the 3 arc G2 fitting
    //!
    real_type
    integralJerk2() const {
      return S0.integralJerk2() +
             S1.integralJerk2() +
             SM.integralJerk2();
    }

    //!
    //! \return get the integral of the snap squared of the 3 arc G2 fitting
    //!
    real_type
    integral_snap2() const {
      return S0.integral_snap2() +
             S1.integral_snap2() +
             SM.integral_snap2();
    }

    real_type integralSnap2() const { return integral_snap2(); }

    //!
    //! \param[out] thMin minimum angle in the 3 arc G2 fitting curve
    //! \param[out] thMax maximum angle in the 3 arc G2 fitting curve
    //! \return the difference of `thMax` and `thMin`
    //!
    real_type
    theta_min_max( real_type & thMin, real_type & thMax ) const;

    //!
    //! Return the difference of maximum-minimum angle in the 3 arc G2 fitting curve
    //!
    real_type
    delta_theta() const
    { real_type thMin, thMax; return theta_min_max( thMin, thMax ); }

    real_type deltaTheta() const { return delta_theta(); }

    //!
    //! \param[out] kMin minimum curvature in the 3 arc G2 fitting curve
    //! \param[out] kMax maximum curvature in the 3 arc G2 fitting curve
    //! \return the difference of `kMax` and `kMin`
    //!
    real_type
    curvature_min_max( real_type & kMin, real_type & kMax ) const;

    //!
    //! Return angle as a function of curvilinear coordinate
    //!
    real_type theta( real_type s ) const;

    //!
    //! Return angle derivative (curvature) as a function of curvilinear coordinate
    //!
    real_type theta_D( real_type s ) const;

    //!
    //! Return angle second derivative (curvature derivative) as a function of curvilinear coordinate
    //!
    real_type theta_DD( real_type s ) const;

    //!
    //! Return angle third derivative as a function of curvilinear coordinate
    //!
    real_type theta_DDD( real_type s ) const;

    //!
    //! Return x coordinate of the3 arc clothoid as a function of curvilinear coordinate
    //!
    real_type X( real_type s ) const;

    //!
    //! Return y coordinate of the3 arc clothoid as a function of curvilinear coordinate
    //!
    real_type Y( real_type s ) const;

    //!
    //! Return initial x coordinate of the 3 arc clothoid
    //!
    real_type x_begin() const { return S0.x_begin(); }

    //!
    //! Return initial y coordinate of the 3 arc clothoid
    //!
    real_type y_begin() const { return S0.y_begin(); }

    //!
    //! Return initial curvature of the 3 arc clothoid
    //!
    real_type kappa_begin() const { return S0.kappa_begin(); }

    //!
    //! Return initial angle of the 3 arc clothoid
    //!
    real_type theta_begin() const { return S0.theta_begin(); }

    //!
    //! Return final x coordinate of the 3 arc clothoid
    //!
    real_type x_end()const { return S1.x_end(); }

    //!
    //! Return final y coordinate of the 3 arc clothoid
    //!
    real_type y_end() const { return S1.y_end(); }

    //!
    //! Return final curvature of the 3 arc clothoid
    //!
    real_type kappa_end() const { return S1.kappa_end(); }

    //!
    //! Return final angle of the 3 arc clothoid
    //!
    real_type theta_end() const { return S1.theta_end(); }

    //!
    //! Compute parameters of 3 arc clothoid at curvilinear coordinate `s`
    //!
    //! \param[in]  s     curvilinear coordinate of where curve is computed
    //! \param[out] theta the curve angle
    //! \param[out] kappa the curve curvature
    //! \param[out] x     the curve x-coordinate
    //! \param[out] y     the curve y-coordinate
    //!
    void
    eval(
      real_type   s,
      real_type & theta,
      real_type & kappa,
      real_type & x,
      real_type & y
    ) const;

    //!
    //! x and y-coordinate at curvilinear coordinate `s`
    //!
    void eval( real_type s, real_type & x, real_type & y ) const;

    //!
    //! x and y-coordinate derivative at curvilinear coordinate `s`
    //!
    void eval_D( real_type s, real_type & x_D, real_type & y_D ) const;

    //!
    //! x and y-coordinate second derivative at curvilinear coordinate `s`
    //!
    void eval_DD( real_type s, real_type & x_DD, real_type & y_DD ) const;

    //!
    //! x and y-coordinate third derivative at curvilinear coordinate `s`
    //!
    void eval_DDD( real_type s, real_type & x_DDD, real_type & y_DDD ) const;

    //!
    //! x and y-coordinate at curvilinear coordinate `s` with offset
    //!
    void eval_ISO( real_type s, real_type offs, real_type & x, real_type & y ) const;

    //!
    //! x and y-coordinate derivative at curvilinear coordinate `s` with offset
    //!
    void eval_ISO_D( real_type s, real_type offs, real_type & x_D, real_type & y_D ) const;

    //!
    //! x and y-coordinate second derivative at curvilinear coordinate `s` with offset
    //!
    void eval_ISO_DD( real_type s, real_type offs, real_type & x_DD, real_type & y_DD ) const;

    //!
    //! x and y-coordinate third derivative at curvilinear coordinate `s` with offset
    //!
    void eval_ISO_DDD( real_type s, real_type offs, real_type & x_DDD, real_type & y_DDD ) const;

    //!
    //! Rotate curve by angle \f$ theta \f$ centered at point  \f$ (c_x,c_y)\f$
    //!
    //! \param[in] angle angle \f$ theta \f$
    //! \param[in] cx    \f$ c_x\f$
    //! \param[in] cy    \f$ c_y\f$
    //!
    void
    rotate( real_type angle, real_type cx, real_type cy ) {
      S0.rotate( angle, cx, cy );
      S1.rotate( angle, cx, cy );
      SM.rotate( angle, cx, cy );
    }

    //!
    //! Translate curve by \f$ (t_x,t_y) \f$
    //!
    void
    translate( real_type tx, real_type ty ){
      S0.translate( tx, ty );
      S1.translate( tx, ty );
      SM.translate( tx, ty );
    }

    //!
    //! Reverse curve parameterization
    //!
    void
    reverse() {
      ClothoidCurve tmp(S0); S1 = S0; S0 = tmp;
      S0.reverse();
      S1.reverse();
      SM.reverse();
    }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidCurve const & c );

    //! save clothoid list of a file stream
    void save( ostream_type & stream ) const;

    // BACK COMPATIBILITY
    real_type
    thetaTotalVariation() const
    { return theta_total_variation(); }

    real_type
    thetaMinMax( real_type & thMin, real_type & thMax ) const
    { return theta_min_max(thMin,thMax); }
  };

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/
  //!
  //! Manage a piecewise clothoids \f$ G(s) \f$ composed by
  //! n clothoids (not necessarily G2 or G1 connected)
  //!
  //! \rst
  //!
  //!   .. image:: ../../images/G2problem3arc.jpg
  //!      :width: 80%
  //!      :align: center
  //!
  //! \endrst
  //!
  class ClothoidList : public BaseCurve {

    bool                  m_curve_is_closed{false};
    vector<real_type>     m_s0;
    vector<ClothoidCurve> m_clotoid_list;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable Utils::BinarySearch<integer> m_lastInterval;
    #else
    mutable integer m_lastInterval{0};
    #endif

    mutable bool               m_aabb_done{false};
    mutable AABB_TREE          m_aabb_tree;
    mutable real_type          m_aabb_offs{real_type(0)};
    mutable real_type          m_aabb_max_angle{real_type(0)};
    mutable real_type          m_aabb_max_size{real_type(0)};
    mutable vector<Triangle2D> m_aabb_triangles;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable std::mutex m_aabb_mutex;
    #endif

    void
    resetLastInterval() {
      #ifdef CLOTHOIDS_USE_THREADS
      bool ok;
      integer & lastInterval = *m_lastInterval.search( std::this_thread::get_id(), ok );
      #else
      integer & lastInterval = m_lastInterval;
      #endif
      lastInterval = 0;
    }

    integer
    closest_point_internal(
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & DST
    ) const;

  public:

    #include "BaseCurve_using.hxx"

    //!
    //! Build an empty clothoid list
    //!
    ClothoidList()
    { this->resetLastInterval(); }

    ~ClothoidList() override {
      m_s0.clear();
      m_clotoid_list.clear();
      m_aabb_triangles.clear();
    }

    //!
    //! Build a copy of an existing clothoid list
    //!
    ClothoidList( ClothoidList const & s )
    { this->resetLastInterval(); this->copy(s); }

    CurveType type() const override { return CurveType::CLOTHOID_LIST; }

    //!
    //! Initialize the clothoid list
    //!
    void init();

    //!
    //! Reserve memory for `n` clothoid
    //!
    void reserve( integer n );

    //!
    //! Build a clothoid list copying an existing one
    //!
    void copy( ClothoidList const & L );

    //!
    //! Copy an existing clothoid list
    //!
    ClothoidList const & operator = ( ClothoidList const & s )
    { this->copy(s); return *this; }

    //!
    //! Build a clothoid from a line segment
    //!
    explicit ClothoidList( LineSegment const & LS );

    //!
    //! Build a clothoid from a circle arc
    //!
    explicit ClothoidList( CircleArc const & C );

    //!
    //! Build a clothoid from a biarc
    //!
    explicit ClothoidList( Biarc const & B );

    //!
    //! Build a clothoid from a list of biarc
    //!
    explicit ClothoidList( BiarcList const & BL );

    //!
    //! Build a clothoid from a clothoid curve
    //!
    explicit ClothoidList( ClothoidCurve const & CL );

    //!
    //! Build a clothoid from a list line segment
    //!
    explicit ClothoidList( PolyLine const & PL );

    //!
    //! Build a clothoid from G2solve2arc
    //!
    explicit ClothoidList( G2solve2arc const & C );

    //!
    //! Build a clothoid from G2solve3arc
    //!
    explicit ClothoidList( G2solve3arc const & C );

    //!
    //! Build a clothoid from G2solveCLC
    //!
    explicit ClothoidList( G2solveCLC const & C );

    //!
    //! Build a clothoid from a curve
    //!
    explicit ClothoidList( BaseCurve const * pC );

    void build( LineSegment const & );
    void build( CircleArc const & );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( G2solve2arc const & );
    void build( G2solve3arc const & );
    void build( G2solveCLC const & );

    //!
    //! Add a line segment to the tail of clothoid list
    //!
    void push_back( LineSegment const & c );

    //!
    //! Add a circle arc to the tail of clothoid list
    //!
    void push_back( CircleArc const & c );

    //!
    //! Add a biarc to the tail of clothoid list
    //!
    void push_back( Biarc const & c );

    //!
    //! Add a biarc list to the tail of clothoid list
    //!
    void push_back( BiarcList const & c );

    //!
    //! Add a clothoid curve to the tail of clothoid list
    //!
    void push_back( ClothoidCurve const & c );

    //!
    //! Add a clothoid list to the tail of clothoid list
    //!
    void push_back( ClothoidList const & c );

    //!
    //! Add a G2solve2arc to the tail of clothoid list
    //!
    void push_back( G2solve2arc const & c );

    //!
    //! Add a G2solve3arc list to the tail of clothoid list
    //!
    void push_back( G2solve3arc const & c );

    //!
    //! Add a clothoid list to the tail of clothoid list
    //!
    void push_back( G2solveCLC const & c );

    //!
    //! Add a list of line segment to the tail of clothoid list
    //!
    void push_back( PolyLine const & c );

    //!
    //! Add a clothoid to the tail of the clothoid list.
    //!
    //! \param kappa0 initial curvature
    //! \param dkappa derivative of the curvature
    //! \param L      length of the segment
    //!
    void
    push_back(
      real_type kappa0,
      real_type dkappa,
      real_type L
    );

    //!
    //! Add a clothoid to the tail of the clothoid list.
    //! The builded clothoid is translated to the tail of the clothioid list.
    //!
    //! \param x0     initial x
    //! \param y0     initial y
    //! \param theta0 initial angle
    //! \param kappa0 initial curvature
    //! \param dkappa derivative of the curvature
    //! \param L      length of the segment
    //!
    void
    push_back(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type kappa0,
      real_type dkappa,
      real_type L
    );

    //!
    //! Add a clothoid to the tail of the clothoid list solving the G1 problem.
    //! The initial point and angle are taken from the tail of the clothoid list.
    //!
    //! \param x1     final x
    //! \param y1     final y
    //! \param theta1 final angle
    //!
    void
    push_back_G1(
      real_type x1,
      real_type y1,
      real_type theta1
    );

    //!
    //! Add a clothoid to the tail of the clothoid list solving the G1 problem.
    //! The initial point and angle are taken from the tail of the clothoid list.
    //! The builded clothoid is translated to the tail of the clothioid list.
    //!
    //! \param x0     initial x
    //! \param y0     initial y
    //! \param theta0 initial angle
    //! \param x1     final x
    //! \param y1     final y
    //! \param theta1 final angle
    //!
    void
    push_back_G1(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1
    );

    //!
    //! True if curve is closed
    //!
    bool is_closed() const { return m_curve_is_closed; }

    //!
    //! Set clousure flag to true
    //!
    void make_closed() { m_curve_is_closed = true; }

    //!
    //! Set clousure flag to false
    //!
    void make_open() { m_curve_is_closed = false; }

    //!
    //! Difference initial final point x component
    //!
    real_type closure_gap_x() const { return this->x_end() - this->x_begin(); }

    //!
    //! Difference initial final point y component
    //!
    real_type closure_gap_y() const { return this->y_end() - this->y_begin(); }

    //!
    //! Difference initial final tangent x component
    //!
    real_type closure_gap_tx() const { return this->tx_end() - this->tx_begin(); }

    //!
    //! Difference initial final tangent y component
    //!
    real_type closure_gap_ty() const { return this->ty_end() - this->ty_begin(); }

    //!
    //! check if clothoid list is closed
    //!
    //! \param[in] tol_xy position tolerance
    //! \param[in] tol_tg angle (tangent) tolerance
    //!
    //! \return true if curve is closed
    //!
    bool
    closure_check( real_type tol_xy = 1e-6, real_type tol_tg = 1e-6 ) const {
      return std::abs(closure_gap_x())  < tol_xy &&
             std::abs(closure_gap_y())  < tol_xy &&
             std::abs(closure_gap_tx()) < tol_tg &&
             std::abs(closure_gap_ty()) < tol_tg;
    }

    //!
    //! Build clothoid list passing to a list of points
    //! solving a series of G1 fitting problems.
    //! The angle at points are estimated using the routine `xy_to_guess_angle`
    //!
    //! \param[in] n number of points
    //! \param[in] x x-coordinates
    //! \param[in] y y-coordinates
    //!
    //! \return false if routine fails
    //!
    bool
    build_G1(
      integer           n,
      real_type const * x,
      real_type const * y
    );

    //!
    //! Build clothoid list passing to a list of points
    //! solving a series of G1 fitting problems.
    //!
    //! \param[in] n     number of points
    //! \param[in] x     x-coordinates
    //! \param[in] y     y-coordinates
    //! \param[in] theta angles at the points
    //!
    //! \return false if routine fails
    //!
    bool
    build_G1(
      integer           n,
      real_type const * x,
      real_type const * y,
      real_type const * theta
    );

    //!
    //! Build clothoid list with G2 continuity.
    //! The vector `s` contains the breakpoints of the curve.
    //! Between two breakpoint the curvature change linearly (is a clothoid)
    //!
    //! \param[in] x0     initial x
    //! \param[in] y0     initial y
    //! \param[in] theta0 initial angle
    //! \param[in] n      number of segments
    //! \param[in] s      break point of the piecewise curve
    //! \param[in] kappa  curvature at the break point
    //!
    //! \return true if curve is closed
    //!
    bool
    build(
      real_type         x0,
      real_type         y0,
      real_type         theta0,
      integer           n,
      real_type const * s,
      real_type const * kappa
    );

    //!
    //! Build clothoid list with G2 continuity.
    //! The vector `s` contains the breakpoints of the curve.
    //! Between two breakpoint the curvature change linearly (is a clothoid)
    //!
    //! \param[in] x0     initial x
    //! \param[in] y0     initial y
    //! \param[in] theta0 initial angle
    //! \param[in] s      break point of the piecewise curve
    //! \param[in] kappa  curvature at the break point
    //!
    //! \return true if curve is closed
    //!
    bool
    build(
      real_type                 x0,
      real_type                 y0,
      real_type                 theta0,
      vector<real_type> const & s,
      vector<real_type> const & kappa
    ) {
      if ( s.size() != kappa.size() ) return false;
      return build(
        x0, y0, theta0,
        integer(s.size()),
        &s.front(), &kappa.front()
      );
    }

    //!
    //! Build clothoid listy using raw data.
    //!
    //! \param[in] n        number of points
    //! \param[in] x        x-coordinates
    //! \param[in] y        y-coordinates
    //! \param[in] abscissa break point of the piecewise curve
    //! \param[in] theta    angles at breakpoints
    //! \param[in] kappa    curvature at the break point
    //!
    //! \return false if fails
    //!
    bool
    build_raw(
      integer           n,
      real_type const * x,
      real_type const * y,
      real_type const * abscissa,
      real_type const * theta,
      real_type const * kappa
    );

    //!
    //! Build clothoid listy using raw data.
    //!
    //! \param[in] x        x-coordinates
    //! \param[in] y        y-coordinates
    //! \param[in] abscissa break point of the piecewise curve
    //! \param[in] theta    angles at breakpoints
    //! \param[in] kappa    curvature at the break point
    //!
    //! \return false if fails
    //!
    bool
    build_raw(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & abscissa,
      vector<real_type> const & theta,
      vector<real_type> const & kappa
    ) {
      integer n = integer(x.size());
      if ( n != integer(y.size())        ||
           n != integer(abscissa.size()) ||
           n != integer(theta.size())    ||
           n != integer(kappa.size()) ) return false;
      return build_raw(
        n, &x.front(), &y.front(),
        &abscissa.front(), &theta.front(), &kappa.front()
      );
    }

    //!
    //! Get the `idx`-th clothoid of the list
    //!
    ClothoidCurve const & get( integer idx ) const;

    //!
    //! Get the `idx`-th clothoid of the list where `idx` is the clothoid at parameter `s`
    //!
    ClothoidCurve const & get_at_s( real_type s ) const;

    //!
    //! Return the numbber of clothoid of the list
    //!
    integer num_segments() const { return integer(m_clotoid_list.size()); }

    //!
    //! The list of clothoid has total length \f$ L \f$
    //! the parameter \f$ s \f$ us recomputed as \f$ s+kL\f$ in such a way
    //! \f$ s+kL\in[0,L)\f$ with \f$ k\in\mathbb{Z} \f$.
    //!
    void wrap_in_range( real_type & s ) const;

    //!
    //! Find the clothoid segment whose definiton range contains `s`
    //!
    integer find_at_s( real_type & s ) const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type length() const override;
    real_type length_ISO( real_type offs ) const override;

    //!
    //! Return the length of the `nseg`-th clothoid of the list
    //!
    real_type
    segment_length( integer nseg ) const;

    //!
    //! Return the length of the `nseg`-th clothoid of the list with offset
    //!
    real_type
    segment_length_ISO( integer nseg, real_type offs ) const;

    //!
    //! Return the length of the `nseg`-th clothoid of the list with offset
    //!
    real_type
    segment_length_SAE( integer nseg, real_type offs ) const
    { return segment_length_ISO( nseg, -offs ); }

    /*\
     |  _    _   _____    _                _
     | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
     | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
     | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
     |                               |___/
    \*/

    void
    bbTriangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override;

    void
    bbTriangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override;

    void
    bbTriangles_SAE(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      this->bbTriangles_ISO( -offs, tvec, max_angle, max_size, icurve );
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    void
    build_AABBtree_ISO(
      real_type offs,
      real_type max_angle = Utils::m_pi/6, // 30 degree
      real_type max_size  = 1e100
    ) const;
    #endif

    /*\
     |   _     _
     |  | |__ | |__   _____  __
     |  | '_ \| '_ \ / _ \ \/ /
     |  | |_) | |_) | (_) >  <
     |  |_.__/|_.__/ \___/_/\_\
    \*/

    void
    bbox(
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
    ) const override {
      bbox_ISO( 0, xmin, ymin, xmax, ymax );
    }

    void
    bbox_ISO(
      real_type   offs,
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
    ) const override;

    /*\
     |   ____             _          _______           _
     |  | __ )  ___  __ _(_)_ __    / / ____|_ __   __| |
     |  |  _ \ / _ \/ _` | | '_ \  / /|  _| | '_ \ / _` |
     |  | |_) |  __/ (_| | | | | |/ / | |___| | | | (_| |
     |  |____/ \___|\__, |_|_| |_/_/  |_____|_| |_|\__,_|
     |              |___/
    \*/

    real_type
    theta_begin() const override
    { return m_clotoid_list.front().theta_begin(); }

    real_type
    theta_end() const override
    { return m_clotoid_list.back().theta_end(); }

    real_type
    x_begin() const override
    { return m_clotoid_list.front().x_begin(); }

    real_type
    y_begin() const override
    { return m_clotoid_list.front().y_begin(); }

    real_type
    x_end() const override
    { return m_clotoid_list.back().x_end(); }

    real_type
    y_end() const override
    { return m_clotoid_list.back().y_end(); }

    real_type
    x_begin_ISO( real_type offs ) const override
    { return m_clotoid_list.front().x_begin_ISO( offs ); }

    real_type
    y_begin_ISO( real_type offs ) const override
    { return m_clotoid_list.front().y_begin_ISO( offs ); }

    real_type
    x_end_ISO( real_type offs ) const override
    { return m_clotoid_list.back().x_end_ISO( offs ); }

    real_type
    y_end_ISO( real_type offs ) const override
    { return m_clotoid_list.back().y_end_ISO( offs ); }

    real_type
    tx_begin() const override
    { return m_clotoid_list.front().tx_begin(); }

    real_type
    ty_begin() const override
    { return m_clotoid_list.front().ty_begin(); }

    real_type
    tx_end() const override
    { return m_clotoid_list.back().tx_end(); }

    real_type
    ty_end() const override
    { return m_clotoid_list.back().ty_end(); }

    real_type
    nx_begin_ISO() const override
    { return m_clotoid_list.front().nx_begin_ISO(); }

    real_type
    ny_begin_ISO() const override
    { return m_clotoid_list.front().ny_begin_ISO(); }

    real_type
    nx_end_ISO() const override
    { return m_clotoid_list.back().nx_end_ISO(); }

    real_type
    ny_end_ISO() const override
    { return m_clotoid_list.back().ny_end_ISO(); }

    /*\
     |  _   _          _
     | | |_| |__   ___| |_ __ _
     | | __| '_ \ / _ \ __/ _` |
     | | |_| | | |  __/ || (_| |
     |  \__|_| |_|\___|\__\__,_|
    \*/

    real_type theta    ( real_type s ) const override;
    real_type theta_D  ( real_type s ) const override;
    real_type theta_DD ( real_type s ) const override;
    real_type theta_DDD( real_type s ) const override;

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type s ) const override;
    real_type ty    ( real_type s ) const override;
    real_type tx_D  ( real_type s ) const override;
    real_type ty_D  ( real_type s ) const override;
    real_type tx_DD ( real_type s ) const override;
    real_type ty_DD ( real_type s ) const override;
    real_type tx_DDD( real_type s ) const override;
    real_type ty_DDD( real_type s ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg(
      real_type   s,
      real_type & tg_x,
      real_type & tg_y
    ) const override;

    void
    tg_D(
      real_type   s,
      real_type & tg_x_D,
      real_type & tg_y_D
    ) const override;

    void
    tg_DD(
      real_type   s,
      real_type & tg_x_DD,
      real_type & tg_y_DD
    ) const override;

    void
    tg_DDD(
      real_type   s,
      real_type & tg_x_DDD,
      real_type & tg_y_DDD
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate(
      real_type   s,
      real_type & th,
      real_type & k,
      real_type & x,
      real_type & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate_ISO(
      real_type   s,
      real_type   offs,
      real_type & th,
      real_type & k,
      real_type & x,
      real_type & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X    ( real_type s ) const override;
    real_type Y    ( real_type s ) const override;
    real_type X_D  ( real_type s ) const override;
    real_type Y_D  ( real_type s ) const override;
    real_type X_DD ( real_type s ) const override;
    real_type Y_DD ( real_type s ) const override;
    real_type X_DDD( real_type s ) const override;
    real_type Y_DDD( real_type s ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const override;

    void
    eval_D(
      real_type   s,
      real_type & x_D,
      real_type & y_D
    ) const override;

    void
    eval_DD(
      real_type   s,
      real_type & x_DD,
      real_type & y_DD
    ) const override;

    void
    eval_DDD(
      real_type   s,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override;

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    real_type X_ISO    ( real_type s, real_type offs ) const override;
    real_type Y_ISO    ( real_type s, real_type offs ) const override;
    real_type X_ISO_D  ( real_type s, real_type offs ) const override;
    real_type Y_ISO_D  ( real_type s, real_type offs ) const override;
    real_type X_ISO_DD ( real_type s, real_type offs ) const override;
    real_type Y_ISO_DD ( real_type s, real_type offs ) const override;
    real_type X_ISO_DDD( real_type s, real_type offs ) const override;
    real_type Y_ISO_DDD( real_type s, real_type offs ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const override;

    void
    eval_ISO_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const override;

    void
    eval_ISO_DD(
      real_type   s,
      real_type   offs,
      real_type & x_DD,
      real_type & y_DD
    ) const override;

    void
    eval_ISO_DDD(
      real_type   s,
      real_type   offs,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override;

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void translate( real_type tx, real_type ty ) override;
    void rotate( real_type angle, real_type cx, real_type cy ) override;
    void scale( real_type sc ) override;
    void reverse() override;
    void change_origin( real_type newx0, real_type newy0 ) override;
    void trim( real_type s_begin, real_type s_end ) override;
    void trim( real_type s_begin, real_type s_end, ClothoidList & newCL ) const;

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    //!
    //! \param  qx     x-coordinate of the point
    //! \param  qy     y-coordinate of the point
    //! \param  x      x-coordinate of the projected point on the curve
    //! \param  y      y-coordinate of the projected point on the curve
    //! \param  s      parameter on the curve of the projection
    //! \param  t      curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst    distance point projected point
    //! \return n >= 0 point is projected orthogonal, n is the number of the segment at minimum distance<br>
    //!        -(n+1)  minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_ISO(
      real_type   qx,
      real_type   qy,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst
    ) const override;

    //!
    //! \param  qx    x-coordinate of the point
    //! \param  qy    y-coordinate of the point
    //! \param  offs  offset of the curve
    //! \param  x     x-coordinate of the projected point on the curve
    //! \param  y     y-coordinate of the projected point on the curve
    //! \param  s     parameter on the curve of the projection
    //! \param  t     curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst   distance point projected point
    //! \return n > 0 point is projected orthogonal, n-1 is the number of the segment at minimum distance<br>
    //!        -(n+1) minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_ISO(
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst
    ) const override;

    //!
    //! Compute the point on clothoid at minimal distance from a given point
    //! using the optimized algorithm described in the publication:
    //!
    //! \param  ds sampling step
    //! \param  qx x-coordinate of the given point
    //! \param  qy y-coordinate of the given point
    //! \param  X  x-coordinate of the point on clothoid at minimal distance
    //! \param  Y  y-coordinate of the point on clothoid at minimal distance
    //! \param  S  curvilinear coordinate of the point (X,Y) on the clothoid
    //! \return the distance of the point from the clothoid
    //!
    real_type
    closest_point_by_sample(
      real_type   ds,
      real_type   qx,
      real_type   qy,
      real_type & X,
      real_type & Y,
      real_type & S
    ) const;

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    //!
    //! \param  qx  x-coordinate of the point
    //! \param  qy  y-coordinate of the point
    //! \return the segment at minimal distance from point (qx,qy)
    //!
    integer
    closest_segment( real_type qx, real_type qy ) const;

    //!
    //! \param  qx           x-coordinate of the point
    //! \param  qy           y-coordinate of the point
    //! \param  icurve_begin index of the initial segment
    //! \param  icurve_end   index of the past to the last segment
    //! \param  x            x-coordinate of the projected point on the curve
    //! \param  y            y-coordinate of the projected point on the curve
    //! \param  s            parameter on the curve of the projection
    //! \param  t            curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst          distance point projected point
    //! \param  icurve       number of the segment with the projected point
    //! \return 1            point is projected orthogonal<br>
    //!         0 =          more than one projection (first returned)<br>
    //!        -1 =          minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_in_range_ISO(
      real_type   qx,
      real_type   qy,
      integer     icurve_begin,
      integer     icurve_end,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst,
      integer   & icurve
    ) const;

    //!
    //! \param  qx           x-coordinate of the point
    //! \param  qy           y-coordinate of the point
    //! \param  icurve_begin index of the initial segment
    //! \param  icurve_end   index of the past to the last segment
    //! \param  x            x-coordinate of the projected point on the curve
    //! \param  y            y-coordinate of the projected point on the curve
    //! \param  s            parameter on the curve of the projection
    //! \param  t            curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst          distance point projected point
    //! \param  icurve       number of the segment with the projected point
    //! \return 1            point is projected orthogonal<br>
    //!         0            = more than one projection (first returned)<br>
    //!        -1            = minimum point is not othogonal projection to curve<br>
    //!
    integer
    closest_point_in_range_SAE(
      real_type   qx,
      real_type   qy,
      integer     icurve_begin,
      integer     icurve_end,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst,
      integer   & icurve
    ) const {
      integer res = this->closest_point_in_range_ISO(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
      t = -t;
      return res;
    }

    //!
    //! \param  qx      x-coordinate of the point
    //! \param  qy      y-coordinate of the point
    //! \param  s_begin initial curvilinear coordinate of the search range
    //! \param  s_end   final curvilinear coordinate of the search range
    //! \param  x       x-coordinate of the projected point on the curve
    //! \param  y       y-coordinate of the projected point on the curve
    //! \param  s       parameter on the curve of the projection
    //! \param  t       curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst     distance point projected point
    //! \param  icurve  number of the segment with the projected point
    //!
    //! \return 1 ok -1 projection failed
    //!
    integer
    closest_point_in_s_range_ISO(
      real_type   qx,
      real_type   qy,
      real_type   s_begin,
      real_type   s_end,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst,
      integer   & icurve
    ) const;

    //!
    //! \param  qx      x-coordinate of the point
    //! \param  qy      y-coordinate of the point
    //! \param  s_begin initial curvilinear coordinate of the search range
    //! \param  s_end   final curvilinear coordinate of the search range
    //! \param  x       x-coordinate of the projected point on the curve
    //! \param  y       y-coordinate of the projected point on the curve
    //! \param  s       parameter on the curve of the projection
    //! \param  t       curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst     distance point projected point
    //! \param  icurve  number of the segment with the projected point
    //!
    //! \return 1 ok -1 projection failed
    //!
    integer
    closest_point_in_s_range_SAE(
      real_type   qx,
      real_type   qy,
      integer     s_begin,
      integer     s_end,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst,
      integer   & icurve
    ) const {
      integer res = this->closest_point_in_s_range_ISO(
        qx, qy, s_begin, s_end, x, y, s, t, dst, icurve
      );
      t = -t;
      return res;
    }

    void
    info( ostream_type & stream ) const override
    { stream << "ClothoidList\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidList const & CL );

    //!
    //! Return the clothoid list as a list of nodes and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] kappa curvature
    //!
    void
    get_SK( real_type * s, real_type * kappa ) const;

    //!
    //! Return the clothoid list as a list of nodes and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] kappa curvature
    //!
    void
    get_SK(
      vector<real_type> & s,
      vector<real_type> & kappa
    ) const {
      s.resize( m_clotoid_list.size()+1 );
      kappa.resize( m_clotoid_list.size()+1 );
      get_SK( &s.front(), &kappa.front() );
    }

    //!
    //! Return the clothoid list as a list of nodes angles and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] theta angles
    //! \param[out] kappa curvature
    //!
    void
    get_STK(
      real_type * s,
      real_type * theta,
      real_type * kappa
    ) const;

    //!
    //! Return the clothoid list as a list of nodes angles and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] theta angles
    //! \param[out] kappa curvature
    //!
    void
    get_STK(
      vector<real_type> & s,
      vector<real_type> & theta,
      vector<real_type> & kappa
    ) const {
      s.resize( m_clotoid_list.size()+1 );
      theta.resize( m_clotoid_list.size()+1 );
      kappa.resize( m_clotoid_list.size()+1 );
      get_STK( &s.front(), &theta.front(), &kappa.front() );
    }

    //!
    //! Return the points of the clothoid list at breakpoints
    //!
    //! \param[out] x x-coordinates
    //! \param[out] y y-coordinates
    //!
    void
    get_XY( real_type * x, real_type * y ) const;

    void
    get_delta_theta( real_type * delta_theta ) const;

    void
    get_delta_kappa( real_type * deltaKappa ) const;

    //!
    //! Find parametric coordinate.
    //!
    //! \param  x    x-coordinate point
    //! \param  y    y-coordinate point
    //! \param  s    value \f$ s \f$
    //! \param  t    value \f$ t \f$
    //! \return idx  the segment with point at minimal distance, otherwise
    //!              -(idx+1) if (x,y) cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      real_type   x,
      real_type   y,
      real_type & s,
      real_type & t
    ) const;

    //!
    //! Find parametric coordinate.
    //!
    //! \param  ibegin initial segment to compute the distance
    //! \param  iend   final segment to compute the distance
    //! \param  x      x-coordinate point
    //! \param  y      y-coordinate point
    //! \param  s      value \f$ s \f$
    //! \param  t      value \f$ t \f$
    //! \return idx    the segment with point at minimal distance, otherwise
    //!                -(idx+1) if (x,y) cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      integer     ibegin,
      integer     iend,
      real_type   x,
      real_type   y,
      real_type & s,
      real_type & t
    ) const;

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    //!
    //! Detect a collision with another clothoid list
    //!
    bool
    collision( ClothoidList const & CL ) const {
      return collision_ISO( 0, CL, 0 );
    }

    //!
    //! Detect a collision with another clothoid list with offset
    //!
    //! \param[in] offs   offset of first clothoid list
    //! \param[in] CL     second clothoid list
    //! \param[in] offs_C offset of second clothoid list
    //!
    bool
    collision_ISO(
      real_type            offs,
      ClothoidList const & CL,
      real_type            offs_C
    ) const;

    bool
    collision( BaseCurve const * pC ) const override;

    bool
    collision_ISO(
      real_type         offs,
      BaseCurve const * pC,
      real_type         offs_C
    ) const override;

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    //!
    //! Intersect a clothoid list with another clothoid list
    //!
    //! \param[in]  CL    second clothoid list
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      ClothoidList const & CL,
      IntersectList      & ilist
    ) const {
      this->intersect_ISO( 0, CL, 0, ilist );
    }

    //!
    //! Intersect a clothoid list with another clothoid list with offset (ISO)
    //!
    //! \param[in]  offs     offset of first clothoid list
    //! \param[in]  CL       second clothoid list
    //! \param[in]  offs_obj offset of second clothoid list
    //! \param[out] ilist    list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type            offs,
      ClothoidList const & CL,
      real_type            offs_obj,
      IntersectList      & ilist
    ) const;

    void
    intersect(
      BaseCurve const * pC,
      IntersectList   & ilist
    ) const override;

    void
    intersect_ISO(
      real_type         offs,
      BaseCurve const * pC,
      real_type         offs_LS,
      IntersectList   & ilist
    ) const override;

    //!
    //! Save Clothoid list to a stream
    //!
    //!\param stream stream to save
    //!
    void
    export_table( ostream_type & stream ) const;

    //!
    //! Save Clothoid list to a stream
    //!
    //!\param stream streamstream to save
    //!
    void
    export_ruby( ostream_type & stream ) const;

    //!
    //! Save the clothoid list on a stream.
    //! The data is saved as follows
    //!
    //!       # x y theta kappa
    //!       x0 y0 theta0 kappa0
    //!       x1 y1 theta1 kappa1
    //!       ...
    //!       xn yn thetan kappan
    //!
    void save( ostream_type & stream ) const;

    //!
    //! Read the clothoid list from a stream.
    //! The data is assumed to be saved as follows
    //!
    //!       # x y theta kappa
    //!       x0 y0 theta0 kappa0
    //!       x1 y1 theta1 kappa1
    //!       ...
    //!       xn yn thetan kappan
    //!
    void load( istream_type & stream, real_type epsi = 1e-8 );

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "ClothoidList_compatibility.hxx"
#endif

  };

  /*\
   |
   |    ___ _     _   _        _    _ ___      _ _           ___ ___
   |   / __| |___| |_| |_  ___(_)__| / __|_ __| (_)_ _  ___ / __|_  )
   |  | (__| / _ \  _| ' \/ _ \ / _` \__ \ '_ \ | | ' \/ -_) (_ |/ /
   |   \___|_\___/\__|_||_\___/_\__,_|___/ .__/_|_|_||_\___|\___/___|
   |                                     |_|
  \*/

  //!
  //! Class for the computation of G2 spljne of clothoids
  //!
  class ClothoidSplineG2 {
  public:
    using TargetType = enum class TargetType : integer
    { P1, P2, P3, P4, P5, P6, P7, P8, P9 };

    static
    inline
    string
    to_string( TargetType n ) {
      string res = "";
      switch ( n ) {
      case TargetType::P1: res = "P1"; break;
      case TargetType::P2: res = "P2"; break;
      case TargetType::P3: res = "P3"; break;
      case TargetType::P4: res = "P4"; break;
      case TargetType::P5: res = "P5"; break;
      case TargetType::P6: res = "P6"; break;
      case TargetType::P7: res = "P7"; break;
      case TargetType::P8: res = "P8"; break;
      case TargetType::P9: res = "P9"; break;
      }
      return res;
    };

  private:

    Utils::Malloc<real_type> real_values{"ClothoidSplineG2"};

    real_type * m_x{nullptr};
    real_type * m_y{nullptr};
    TargetType  m_tt{TargetType::P1};
    real_type   m_theta_I{real_type(0)};
    real_type   m_theta_F{real_type(0)};
    integer     m_npts{0};

    // work vector
    mutable real_type * m_k{nullptr};
    mutable real_type * m_dk{nullptr};
    mutable real_type * m_L{nullptr};
    mutable real_type * m_kL{nullptr};
    mutable real_type * m_L_1{nullptr};
    mutable real_type * m_L_2{nullptr};
    mutable real_type * m_k_1{nullptr};
    mutable real_type * m_k_2{nullptr};
    mutable real_type * m_dk_1{nullptr};
    mutable real_type * m_dk_2{nullptr};

    real_type
    diff2pi( real_type in ) const {
      return in-Utils::m_2pi*round(in/Utils::m_2pi);
    }

  public:

    ClothoidSplineG2() = default;
    ~ClothoidSplineG2() = default;

    void
    setP1( real_type theta0, real_type thetaN ) {
      m_tt      = TargetType::P1;
      m_theta_I = theta0;
      m_theta_F = thetaN;
    }

    void setP2() { m_tt = TargetType::P2; }
    void setP3() { m_tt = TargetType::P3; }
    void setP4() { m_tt = TargetType::P4; }
    void setP5() { m_tt = TargetType::P5; }
    void setP6() { m_tt = TargetType::P6; }
    void setP7() { m_tt = TargetType::P7; }
    void setP8() { m_tt = TargetType::P8; }
    void setP9() { m_tt = TargetType::P9; }

    void
    build(
      real_type const * xvec,
      real_type const * yvec,
      integer           npts
    );

    integer numPnts() const { return m_npts; }
    integer numTheta() const;
    integer numConstraints() const;

    void
    guess(
      real_type * theta_guess,
      real_type * theta_min,
      real_type * theta_max
    ) const;

    bool
    objective( real_type const * theta, real_type & f ) const;

    bool
    gradient( real_type const * theta, real_type * g ) const;

    bool
    constraints( real_type const * theta, real_type * c ) const;

    integer
    jacobian_nnz() const;

    bool
    jacobian_pattern( integer * i, integer * j ) const;

    bool
    jacobian_pattern_matlab( real_type * i, real_type * j ) const;

    bool
    jacobian( real_type const * theta, real_type * vals ) const;

    void
    info( ostream_type & stream ) const
    { stream << "ClothoidSplineG2\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidSplineG2 const & c );

  };

}

///
/// eof: ClothoidList.hxx
///
