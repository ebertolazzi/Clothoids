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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: ClothoidList_G2solve3arc.hxx
///

namespace G2lib {

  /*\
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/

  //!
  //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
  //! 3 clothoid and one line segment that solve the \f$ G^2 \f$ problem
  //!
  //! **match**
  //!
  //! \f[
  //! \begin{array}{ll}
  //! \textrm{endpoints:}\quad&
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   G(0) &=& \mathbf{p}_0 \\[0.5em]
  //!   G(L) &=& \mathbf{p}_1
  //! \end{array}
  //! \right.
  //! \\[1em]
  //! \textrm{angles:}\quad&
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   \theta(0) &=& \theta_0 \\[0.5em]
  //!   \theta(L) &=& \theta_1
  //! \end{array}
  //! \right.
  //! \\[1em]
  //! \textrm{curvature:}\quad&
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   \kappa(0) &=& \kappa_0 \\[0.5em]
  //!   \kappa(L) &=& \kappa_1
  //! \end{array}
  //! \right.
  //! \end{array}
  //! \f]
  //!
  //! **Reference**
  //!
  //! The solution algorithm is described in
  //!
  //! - **E.Bertolazzi, M.Frego**, On the \f$ G^2 \f$ Hermite Interpolation Problem with clothoids
  //!   Journal of Computational and Applied Mathematics, vol 341, pp. 99-116, 2018
  //!
  //! @html_image{G2problem3arc.png,width=60%}
  //!
  class G2solve3arc {

    ClothoidCurve m_S0{"G2solve3arc_S0"};
    ClothoidCurve m_SM{"G2solve3arc_SM"};
    ClothoidCurve m_S1{"G2solve3arc_S1"};

    real_type m_tolerance{real_type(1e-10)};
    int       m_max_iter{100};

    // \f$ G^2 \f$ interpolation data
    real_type m_x0{real_type(0)};
    real_type m_y0{real_type(0)};
    real_type m_theta0{real_type(0)};
    real_type m_kappa0{real_type(0)};
    real_type m_x1{real_type(0)};
    real_type m_y1{real_type(0)};
    real_type m_theta1{real_type(0)};
    real_type m_kappa1{real_type(0)};

    // standard scaled problem
    real_type m_phi{real_type(0)};
    real_type m_Lscale{real_type(0)};
    real_type m_th0{real_type(0)};
    real_type m_th1{real_type(0)};
    real_type m_s0{real_type(0)};
    real_type m_s1{real_type(0)};

    // precomputed values
    real_type m_K0{real_type(0)},
              m_K1{real_type(0)},
              m_c0{real_type(0)},
              m_c1{real_type(0)},
              m_c2{real_type(0)},
              m_c3{real_type(0)},
              m_c4{real_type(0)},
              m_c5{real_type(0)},
              m_c6{real_type(0)},
              m_c7{real_type(0)},
              m_c8{real_type(0)},
              m_c9{real_type(0)},
              m_c10{real_type(0)},
              m_c11{real_type(0)},
              m_c12{real_type(0)},
              m_c13{real_type(0)},
              m_c14{real_type(0)};

    void
    evalFJ(
      real_type const vars[2],
      real_type       F[2],
      real_type       J[2][2]
    ) const;

    void evalF( real_type const vars[2], real_type F[2] ) const;

    void build_solution( real_type const sM, real_type const thM );

    int solve( real_type const sM_guess, real_type const thM_guess );

  public:

    G2solve3arc() = default;

    ~G2solve3arc() = default;

    //void setup( GenericContainer const & gc );

    //!
    //! Fix tolerance for the \f$ G^2 \f$ problem
    //!
    void set_tolerance( real_type const tol );

    //!
    //! Fix maximum number of iteration for the \f$ G^2 \f$ problem
    //!
    void set_max_iter( integer const miter );

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
      real_type const s0,
      real_type const x0,
      real_type const y0,
      real_type const theta0,
      real_type const kappa0,
      real_type const s1,
      real_type const x1,
      real_type const y1,
      real_type const theta1,
      real_type const kappa1
    );

    //!
    //! \return get the first clothoid for the 3 arc \f$ G^2 \f$ fitting
    //!
    ClothoidCurve const & S0() const { return m_S0; }

    //!
    //! \return get the last clothoid for the 3 arc \f$ G^2 \f$ fitting
    //!
    ClothoidCurve const & S1() const { return m_S1; }

    //!
    //! \return get the middle clothoid for the 3 arc \f$ G^2 \f$ fitting
    //!
    ClothoidCurve const & SM() const { return m_SM; }

    //!
    //! \return get the length of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    total_length() const {
      return m_S0.length() +
             m_S1.length() +
             m_SM.length();
    }


    //!
    //! \return get the total angle variation of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    theta_total_variation() const {
      return m_S0.theta_total_variation() +
             m_S1.theta_total_variation() +
             m_SM.theta_total_variation();
    }

    //!
    //! \return get the total curvature variation of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    curvature_total_variation() const {
      return m_S0.curvature_total_variation() +
             m_S1.curvature_total_variation() +
             m_SM.curvature_total_variation();
    }

    //!
    //! \return get the integral of the curvature squared of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    integral_curvature2() const {
      return m_S0.integral_curvature2() +
             m_S1.integral_curvature2() +
             m_SM.integral_curvature2();
    }

    //!
    //! \return get the integral of the jerk squared of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    integral_jerk2() const {
      return m_S0.integral_jerk2() +
             m_S1.integral_jerk2() +
             m_SM.integral_jerk2();
    }

    //!
    //! \return get the integral of the snap squared of the 3 arc \f$ G^2 \f$ fitting
    //!
    real_type
    integral_snap2() const {
      return m_S0.integral_snap2() +
             m_S1.integral_snap2() +
             m_SM.integral_snap2();
    }

    //!
    //! \param[out] thMin minimum angle in the 3 arc \f$ G^2 \f$ fitting curve
    //! \param[out] thMax maximum angle in the 3 arc \f$ G^2 \f$ fitting curve
    //! \return the difference of `thMax` and `thMin`
    //!
    real_type
    theta_min_max( real_type & thMin, real_type & thMax ) const;

    //!
    //! Return the difference of maximum-minimum angle in the 3 arc \f$ G^2 \f$ fitting curve
    //!
    real_type
    delta_theta() const
    { real_type thMin, thMax; return theta_min_max( thMin, thMax ); }

    //!
    //! \param[out] kMin minimum curvature in the 3 arc \f$ G^2 \f$ fitting curve
    //! \param[out] kMax maximum curvature in the 3 arc \f$ G^2 \f$ fitting curve
    //! \return the difference of `kMax` and `kMin`
    //!
    real_type curvature_min_max( real_type & kMin, real_type & kMax ) const;

    //!
    //! Return angle as a function of curvilinear coordinate
    //!
    real_type theta( real_type const s ) const;

    //!
    //! Return angle derivative (curvature) as a function of curvilinear coordinate
    //!
    real_type theta_D( real_type const s ) const;

    //!
    //! Return angle second derivative (curvature derivative) as a function of curvilinear coordinate
    //!
    real_type theta_DD( real_type const s ) const;

    //!
    //! Return angle third derivative as a function of curvilinear coordinate
    //!
    real_type theta_DDD( real_type const s ) const;

    //!
    //! Return \f$x\f$-coordinate of the arc clothoid as a function of curvilinear coordinate
    //!
    real_type X( real_type const s ) const;

    //!
    //! Return \f$y\f$-coordinate of the arc clothoid as a function of curvilinear coordinate
    //!
    real_type Y( real_type const s ) const;

    //!
    //! Return initial \f$x\f$-coordinate of the 3 arc clothoid
    //!
    real_type x_begin() const { return m_S0.x_begin(); }

    //!
    //! Return initial \f$y\f$-coordinate of the 3 arc clothoid
    //!
    real_type y_begin() const { return m_S0.y_begin(); }

    //!
    //! Return initial curvature of the 3 arc clothoid
    //!
    real_type kappa_begin() const { return m_S0.kappa_begin(); }

    //!
    //! Return initial angle of the 3 arc clothoid
    //!
    real_type theta_begin() const { return m_S0.theta_begin(); }

    //!
    //! Return final \f$x\f$-coordinate of the 3 arc clothoid
    //!
    real_type x_end()const { return m_S1.x_end(); }

    //!
    //! Return final \f$y\f$-coordinate of the 3 arc clothoid
    //!
    real_type y_end() const { return m_S1.y_end(); }

    //!
    //! Return final curvature of the 3 arc clothoid
    //!
    real_type kappa_end() const { return m_S1.kappa_end(); }

    //!
    //! Return final angle of the 3 arc clothoid
    //!
    real_type theta_end() const { return m_S1.theta_end(); }

    //!
    //! Compute parameters of 3 arc clothoid at curvilinear coordinate \f$s\f$
    //!
    //! \param[in]  s     curvilinear coordinate of where curve is computed
    //! \param[out] theta the curve angle
    //! \param[out] kappa the curve curvature
    //! \param[out] x     the curve \f$x\f$-coordinate
    //! \param[out] y     the curve \f$y\f$-coordinate
    //!
    void
    eval(
      real_type const s,
      real_type     & theta,
      real_type     & kappa,
      real_type     & x,
      real_type     & y
    ) const;

    //!
    //! x and \f$y\f$-coordinate at curvilinear coordinate \f$s\f$
    //!
    void eval( real_type const s, real_type & x, real_type & y ) const;

    //!
    //! x and \f$y\f$-coordinate derivative at curvilinear coordinate \f$s\f$
    //!
    void eval_D( real_type const s, real_type & x_D, real_type & y_D ) const;

    //!
    //! x and \f$y\f$-coordinate second derivative at curvilinear coordinate \f$s\f$
    //!
    void eval_DD( real_type const s, real_type & x_DD, real_type & y_DD ) const;

    //!
    //! x and \f$y\f$-coordinate third derivative at curvilinear coordinate \f$s\f$
    //!
    void eval_DDD( real_type const s, real_type & x_DDD, real_type & y_DDD ) const;

    //!
    //! x and \f$y\f$-coordinate at curvilinear coordinate \f$s\f$ with offset
    //!
    void eval_ISO( real_type const s, real_type offs, real_type & x, real_type & y ) const;

    //!
    //! x and \f$y\f$-coordinate derivative at curvilinear coordinate \f$s\f$ with offset
    //!
    void eval_ISO_D( real_type const s, real_type offs, real_type & x_D, real_type & y_D ) const;

    //!
    //! x and \f$y\f$-coordinate second derivative at curvilinear coordinate \f$s\f$ with offset
    //!
    void eval_ISO_DD( real_type const s, real_type offs, real_type & x_DD, real_type & y_DD ) const;

    //!
    //! x and \f$y\f$-coordinate third derivative at curvilinear coordinate \f$s\f$ with offset
    //!
    void eval_ISO_DDD( real_type const s, real_type offs, real_type & x_DDD, real_type & y_DDD ) const;

    //!
    //! Rotate curve by angle \f$ \theta \f$ centered at point \f$ (c_x,c_y) \f$
    //!
    //! \param[in] angle angle \f$ \theta \f$
    //! \param[in] cx    \f$ c_x \f$
    //! \param[in] cy    \f$ c_y \f$
    //!
    void
    rotate( real_type const angle, real_type const cx, real_type const cy ) {
      m_S0.rotate( angle, cx, cy );
      m_S1.rotate( angle, cx, cy );
      m_SM.rotate( angle, cx, cy );
    }

    //!
    //! Translate curve by \f$ (t_x,t_y) \f$
    //!
    void
    translate( real_type const tx, real_type const ty ){
      m_S0.translate( tx, ty );
      m_S1.translate( tx, ty );
      m_SM.translate( tx, ty );
    }

    //!
    //! Reverse curve parameterization
    //!
    void
    reverse() {
      ClothoidCurve tmp(m_S0); m_S1 = m_S0; m_S0 = tmp;
      m_S0.reverse();
      m_S1.reverse();
      m_SM.reverse();
    }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidCurve const & c );

    //! save clothoid list of a file stream
    void save( ostream_type & stream ) const;

    // BACK COMPATIBILITY
    #ifdef CLOTHOIDS_BACK_COMPATIBILITY
    void setTolerance( real_type tol ) { set_tolerance( tol ); }
    void setMaxIter( integer miter ) { set_max_iter( miter ); }
    real_type thetaTotalVariation() const { return theta_total_variation(); }
    real_type thetaMinMax( real_type & thMin, real_type & thMax ) const { return theta_min_max(thMin,thMax); }
    real_type totalLength() const { return total_length(); }
    real_type integralCurvature2() const { return integral_curvature2(); }
    real_type integralJerk2() const { return integral_jerk2(); }
    real_type integralSnap2() const { return integral_snap2(); }
    real_type deltaTheta() const { return delta_theta(); }
    #endif

  };
}

///
/// eof: ClothoidList_G2solve3arc.hxx
///
