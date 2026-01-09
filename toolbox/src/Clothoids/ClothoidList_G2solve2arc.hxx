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
/// file: ClothoidList_G2solve2arc.hxx
///

namespace G2lib
{

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/


  //!
  //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
  //! two clothoids arc that solve the \f$ G^2 \f$ problem
  //!
  //! \f[
  //! \begin{array}{ll}
  //! \textrm{endpoints:}\quad &
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   G(0) &=& \mathbf{p}_0 \\[0.5em]
  //!   G(L) &=& \mathbf{p}_1
  //! \end{array}
  //! \right.
  //! \\[1em]
  //! \textrm{angles:}\quad &
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   \theta(0) &=& \theta_0 \\[0.5em]
  //!   \theta(L) &=& \theta_1
  //! \end{array}
  //! \right.
  //! \\[1em]
  //! \textrm{curvature:}\quad &
  //! \left\{
  //! \begin{array}{r@{~}c@{~}l}
  //!   \kappa(0) &=& \kappa_0 \\[0.5em]
  //!   \kappa(L) &=& \kappa_1
  //! \end{array}
  //! \right.
  //! \end{array}
  //! \f]
  //!
  //! **note**
  //!
  //! The solution do not exist for all the combination of points/angle/curvature
  //!
  //!
  class G2solve2arc
  {
    real_type m_tolerance{ real_type( 1e-10 ) };
    integer   m_max_iter{ 20 };

    real_type m_x0{ real_type( 0 ) };
    real_type m_y0{ real_type( 0 ) };
    real_type m_theta0{ real_type( 0 ) };
    real_type m_kappa0{ real_type( 0 ) };

    real_type m_x1{ real_type( 0 ) };
    real_type m_y1{ real_type( 0 ) };
    real_type m_theta1{ real_type( 0 ) };
    real_type m_kappa1{ real_type( 0 ) };

    // standard problem
    real_type m_lambda{ real_type( 0 ) };
    real_type m_phi{ real_type( 0 ) };
    real_type m_xbar{ real_type( 0 ) };
    real_type m_ybar{ real_type( 0 ) };
    real_type m_th0{ real_type( 0 ) };
    real_type m_th1{ real_type( 0 ) };
    real_type m_k0{ real_type( 0 ) };
    real_type m_k1{ real_type( 0 ) };
    real_type m_DeltaK{ real_type( 0 ) };
    real_type m_DeltaTheta{ real_type( 0 ) };

    ClothoidCurve m_S0{ "G2solve2arc_S0" };
    ClothoidCurve m_S1{ "G2solve2arc_S1" };

    void evalA( real_type const alpha, real_type const L, real_type & A ) const;

    void evalA( real_type const alpha, real_type const L, real_type & A, real_type & A_1, real_type & A_2 ) const;

    void evalG( real_type const alpha, real_type const L, real_type const th, real_type const k, real_type G[2] ) const;

    void evalG(
      real_type const alpha,
      real_type const L,
      real_type const th,
      real_type const k,
      real_type       G[2],
      real_type       G_1[2],
      real_type       G_2[2] ) const;

    void evalF( real_type const vars[2], real_type F[2] ) const;

    void evalFJ( real_type const vars[2], real_type F[2], real_type J[2][2] ) const;

    void build_solution( real_type alpha, real_type L );

  public:
    //!
    //! Build an empty clothoid list
    //!
    G2solve2arc() = default;

    ~G2solve2arc() = default;

    //!
    //! Construct a piecewise clothoids \f$ G(s) \f$ composed by
    //! two clothoids arc that solve the \f$ G^2 \f$ problem, with data
    //!
    //! \f[
    //!   \mathbf{p}_0 = (x_0,y_0)^T, \qquad \mathbf{p}_1 = (x_1,y_1)^T
    //! \f]
    //! \f[
    //!   \theta_0, \qquad \theta_1, \qquad \kappa_0, \qquad \kappa_1
    //! \f]
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
    int build(
      real_type const x0,
      real_type const y0,
      real_type const theta0,
      real_type const kappa0,
      real_type const x1,
      real_type const y1,
      real_type const theta1,
      real_type const kappa1 );

    //!
    //! Fix tolerance for the \f$ G^2 \f$ problem
    //!
    void set_tolerance( real_type const tol );

    //!
    //! Fix maximum number of iteration for the \f$ G^2 \f$ problem
    //!
    void set_max_iter( integer const miter );

    //!
    //! Solve the \f$ G^2 \f$ problem
    //!
    //! \return number of iterations of -1 if failed
    //!
    int solve();

    //!
    //! Return the first clothoid of the \f$ G^2 \f$ clothoid list
    //!
    ClothoidCurve const & S0() const { return m_S0; }

    //!
    //! Return the second clothoid of the \f$ G^2 \f$ clothoid list
    //!
    ClothoidCurve const & S1() const { return m_S1; }

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
    void setTolerance( real_type const tol ) { set_tolerance( tol ); }
    void setMaxIter( integer const miter ) { set_max_iter( miter ); }
#endif
  };

}  // namespace G2lib

///
/// eof: ClothoidList_G2solve2arc.hxx
///
