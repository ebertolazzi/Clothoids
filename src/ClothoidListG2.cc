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

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wswitch-enum"
#endif

namespace G2lib {

  constexpr bool      VERBOSE{true};
  constexpr real_type SOLVER_TOLERANCE{1.0e-9};
  constexpr integer   MAX_ITERATIONS{100};

  // vecchio per testare

  void
  ClothoidSplineG2::allocate( integer const n ) {
    m_npts = n;
    m_x.resize( n );
    m_y.resize( n );
    m_G1_vec.resize( n-1 );
  }

  void
  ClothoidSplineG2::build(
    real_type const xvec[],
    real_type const yvec[],
    integer   const n
  ) {
    allocate( n );
    std::copy_n( xvec, n, m_x.data() );
    std::copy_n( yvec, n, m_y.data() );
  }

  void
  ClothoidSplineG2::evaluate_for_NLP( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP temporary cc"};
    integer const ne{ m_npts - 1 };
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1( m_x[j],   m_y[j],   theta[j],
                   m_x[j+1], m_y[j+1], theta[j+1] );
      auto & G{ m_G1_vec[j] };
      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end();
    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_D( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_D temporary cc"};
    integer const ne{ m_npts - 1 };
    real_type L_D[2], k_D[2], dk_D[2];
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1_D( m_x[j],   m_y[j],   theta[j],
                     m_x[j+1], m_y[j+1], theta[j+1],
                     L_D, k_D, dk_D );

      auto & G{ m_G1_vec[j] };

      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end();

      G.L__L = L_D[0];
      G.L__R = L_D[1];

      G.k__L = k_D[0];
      G.k__R = k_D[1];

      G.dk__L = dk_D[0];
      G.dk__R = dk_D[1];

    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_DD( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_DD temporary cc"};
    integer const ne{ m_npts - 1 };
    real_type L_D[2], k_D[2], dk_D[2];
    real_type L_DD[3], k_DD[3], dk_DD[3];
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1_DD( m_x[j],   m_y[j],   theta[j],
                      m_x[j+1], m_y[j+1], theta[j+1],
                      L_D, k_D, dk_D, L_DD, k_DD, dk_DD );

      auto & G{ m_G1_vec[j] };

      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end();

      G.L__L  = L_D[0];
      G.L__R  = L_D[1];
      G.L__LL = L_DD[0];
      G.L__LR = L_DD[1];
      G.L__RR = L_DD[2];

      G.k__L  = k_D[0];
      G.k__R  = k_D[1];
      G.k__LL = k_DD[0];
      G.k__LR = k_DD[1];
      G.k__RR = k_DD[2];

      G.dk__L  = dk_D[0];
      G.dk__R  = dk_D[1];
      G.dk__LL = dk_DD[0];
      G.dk__LR = dk_DD[1];
      G.dk__RR = dk_DD[2];

    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_BC( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_BC temporary cc"};
    integer const N[]{ 0, m_npts-2 };
    for ( integer j : N ) {
      cc.build_G1( m_x[j],   m_y[j],   theta[j],
                   m_x[j+1], m_y[j+1], theta[j+1] );

      auto & G{ m_G1_vec[j] };

      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end(); // m_k0[j]+m_dk[j]*m_L[j];
    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_D_BC( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_D_BC temporary cc"};
    real_type L_D[2], k_D[2], dk_D[2];
    integer const N[]{ 0, m_npts-2 };
    for ( integer j : N ) {
      cc.build_G1_D( m_x[j],   m_y[j],   theta[j],
                     m_x[j+1], m_y[j+1], theta[j+1],
                     L_D, k_D, dk_D );

      auto & G{ m_G1_vec[j] };

      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end();

      G.L__L = L_D[0];
      G.L__R = L_D[1];

      G.k__L = k_D[0];
      G.k__R = k_D[1];

      G.dk__L = dk_D[0];
      G.dk__R = dk_D[1];
    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_DD_BC( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_DD_BC temporary cc"};
    real_type L_D[2], k_D[2], dk_D[2];
    real_type L_DD[3], k_DD[3], dk_DD[3];
    integer const N[]{ 0, m_npts-2 };
    for ( integer j : N ) {
      cc.build_G1_DD( m_x[j],   m_y[j],   theta[j],
                      m_x[j+1], m_y[j+1], theta[j+1],
                      L_D, k_D, dk_D, L_DD, k_DD, dk_DD );

      auto & G{ m_G1_vec[j] };

      G.k0 = cc.kappa_begin();
      G.dk = cc.dkappa();
      G.L  = cc.length();
      G.k1 = cc.kappa_end();

      G.L__L  = L_D[0];
      G.L__R  = L_D[1];
      G.L__LL = L_DD[0];
      G.L__LR = L_DD[1];
      G.L__RR = L_DD[2];

      G.k__L  = k_D[0];
      G.k__R  = k_D[1];
      G.k__LL = k_DD[0];
      G.k__LR = k_DD[1];
      G.k__RR = k_DD[2];

      G.dk__L  = dk_D[0];
      G.dk__R  = dk_D[1];
      G.dk__LL = dk_DD[0];
      G.dk__LR = dk_DD[1];
      G.dk__RR = dk_DD[2];

    }
  }
  /*\
   |
   |    ___ _     _   _        _    _ ___      _ _           ___ ___
   |   / __| |___| |_| |_  ___(_)__| / __|_ __| (_)_ _  ___ / __|_  )
   |  | (__| / _ \  _| ' \/ _ \ / _` \__ \ '_ \ | | ' \/ -_) (_ |/ /
   |   \___|_\___/\__|_||_\___/_\__,_|___/ .__/_|_|_||_\___|\___/___|
   |                                     |_|
  \*/

  void
  ClothoidSplineG2::guess(
    real_type theta_guess[],
    real_type theta_min[],
    real_type theta_max[]
  ) const {
    Eigen::Vector<real_type,Eigen::Dynamic> omega, len;
    omega.resize( m_npts );
    len.resize( m_npts );
    G2lib::xy_to_guess_angle( m_npts, m_x.data(), m_y.data(), theta_guess, theta_min, theta_max, omega.data(), len.data() );
  }

  void
  ClothoidSplineG2::build(
    integer   const n,
    real_type const xvec[],
    real_type const yvec[],
    real_type       theta[]
  ) {
  
    using Vector       = Pipal::Vector<real_type>;
    using SparseMatrix = Pipal::SparseMatrix<real_type>;
    
    for ( integer i{0}; i < n; ++i )
      std::cout << fmt::format( "x[{}] = {}, y[{}] = {}\n", i, xvec[i], i, yvec[i] );
  
    allocate( n );

    std::copy_n( xvec, n, m_x.data() );
    std::copy_n( yvec, n, m_y.data() );
    
    Vector theta_guess( m_npts ),
           theta_min( m_npts ),
           theta_max( m_npts ),
           theta_sol( m_npts );

    this->guess( theta_guess.data(), theta_min.data(), theta_max.data() );

    auto solver = std::make_unique<Pipal::Solver<real_type>>("ClothoidSplineG2::build",
      // Objective function
      [this] ( Vector const & theta, real_type & f ) -> bool
      { return this->objective( theta.data(), f ); },
      
      // Gradient of the objective function
      [this] ( Vector const & theta, Vector & out ) -> bool {
        out.resize( this->numTheta() );
        return this->gradient( theta.data(), out.data() );
      },
      
      // Constraints function
      [this] ( Vector const & theta, Vector & out ) -> bool {
        out.resize( this->numConstraints() );
        return this->constraints( theta.data(), out.data() );
      },
      
      // Jacobian of the constraints function
      // jacobian_nnz
      // jacobian_pattern
      // jacobian( theta, vals )
      [this] ( Vector const & theta, SparseMatrix & out ) -> bool {
        return this->jacobian( theta.data(), out );
      },
      
      // Hessian of the Lagrangian
      [this] ( Vector const & theta, Vector const & lambda, SparseMatrix & out ) -> bool {
        return this->lagrangian_hessian( theta.data(), lambda.data(), out );
      },

      // Lower bounds on the primal variables
      [this,&theta_min] ( Vector & out ) -> bool
      { out.resize( m_npts ); out.noalias() = theta_min; return true; },

      // Upper bounds on the primal variables
      [this,&theta_max] ( Vector & out ) -> bool
      { out.resize( m_npts ); out.noalias() = theta_max; return true; },

      // Lower bounds on the constraints
      [this] ( Vector & out ) -> bool
      { out.resize( this->numConstraints() ); out.setZero(); return true; },

      // Upper bounds on the constraints
      [this] ( Vector & out ) -> bool
      { out.resize( this->numConstraints() ); out.setZero(); return true; }
    );
    solver->algorithm(Pipal::Algorithm::CONSERVATIVE);
    solver->verbose_mode(VERBOSE);
    solver->tolerance(SOLVER_TOLERANCE);
    solver->max_iterations(MAX_ITERATIONS);
    bool ok{ solver->optimize( theta_guess, theta_sol ) };
    UTILS_ASSERT( ok, "ClothoidSplineG2::build( n={}, x, y ) failed PIPAL solver", m_npts );
    
    std::copy_n( theta_sol.data(), n, theta );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidSplineG2::numTheta() const { return m_npts; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidSplineG2::numConstraints() const {
    switch (m_tt) {
      case TargetType::P1:
      case TargetType::P2: return m_npts;
      default: break;
    }
    return m_npts-2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::objective( real_type const theta[], real_type & f ) const {
    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };
    switch (m_tt) {
    case TargetType::P1:
    case TargetType::P2:
      f = 0;
      break;
    case TargetType::P3:
      // forward target
      break;
    case TargetType::P4: {
        evaluate_for_NLP_BC(theta);
        auto const & G0{ m_G1_vec[0] };
        auto const & GE{ m_G1_vec[ne1] };
        f = G0.dk*G0.dk+GE.dk*GE.dk;
      }
      break;
    case TargetType::P5: {
        evaluate_for_NLP_BC(theta);
        auto const & G0{ m_G1_vec[0] };
        auto const & GE{ m_G1_vec[ne1] };
        f = G0.L+GE.L;
      }
      break;
    case TargetType::P6:
      f = 0;
      evaluate_for_NLP(theta);
      for ( integer j{0}; j < ne; ++j ) f += m_G1_vec[j].L;
      break;
    case TargetType::P7:
      f = 0;
      evaluate_for_NLP(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        f = f + G.L * ( G.L* ( G.dk*( (G.dk*G.L)/3 + G.k0) ) + G.k0*G.k0 );
      }
      break;
    case TargetType::P8:
      f = 0;
      evaluate_for_NLP(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        f += G.L*G.dk*G.dk;
      }
      break;
    case TargetType::P9:
      f = 0;
      evaluate_for_NLP(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        real_type const L    { G.L };
        real_type const kur  { G.k0 };
        real_type const dk   { G.dk };
        real_type const k2   { kur*kur };
        real_type const k3   { k2*kur };
        real_type const k4   { k2*k2 };
        real_type const dk2  { dk*dk };
        real_type const dk3  { dk*dk2 };
        f += (k4+dk2+(2*k3*dk+(2*k2*dk2+(dk3*(kur+dk*L/5))*L)*L)*L)*L;
      }
      break;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::gradient(
    real_type const theta[],
    real_type       g[]
  ) const {
    std::fill_n( g, m_npts, 0 );
    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };
    switch (m_tt) {
    case TargetType::P1:
    case TargetType::P2:
      break;
    case TargetType::P3:
      break;
    case TargetType::P4:
      evaluate_for_NLP_D_BC(theta);
      {
        auto const & G0{ m_G1_vec[0] };
        auto const & GE{ m_G1_vec[ne1] };
        g[0]   = 2*G0.dk*G0.dk__L;
        g[1]   = 2*G0.dk*G0.dk__R;
        g[ne1] = 2*GE.dk*GE.dk__L;
        g[ne]  = 2*GE.dk*GE.dk__R;
      }
      break;
    case TargetType::P5: {
        evaluate_for_NLP_D_BC(theta);
        auto const & G0{ m_G1_vec[0] };
        auto const & GE{ m_G1_vec[ne1] };
        g[0]   = G0.L__L;
        g[1]   = G0.L__R;
        g[ne1] = GE.L__L;
        g[ne]  = GE.L__R;
      }
      break;
    case TargetType::P6:
      evaluate_for_NLP_D(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        g[j]   += G.L__L;
        g[j+1] += G.L__R;
      }
      break;
    case TargetType::P7:
      evaluate_for_NLP_D(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        real_type const L     { G.L };
        real_type const L__L  { G.L__L };
        real_type const L__R  { G.L__R };
        real_type const L2    { L*L };
        real_type const L3    { L*L2 };
        real_type const kur   { G.k0 };
        real_type const k__L  { G.k__L };
        real_type const k__R  { G.k__R };
        real_type const k2    { kur*kur };
        real_type const dk    { G.dk };
        real_type const dk__L { G.dk__L };
        real_type const dk__R { G.dk__R };
        real_type const dk2   { dk*dk };
        g[j]   += 2*(dk*dk__L*L3)/3
                  + (dk2*L2*L__L)
                  + dk__L*L2*kur
                  + 2*dk*L*L__L*kur
                  + dk*L2*k__L
                  + L__L*k2
                  + 2*L*kur*k__L;
        g[j+1] += 2*(dk*dk__R*L3)/3
                  + (dk2*L2*L__R)
                  + dk__R*L2*kur
                  + 2*dk*L*L__R*kur
                  + dk*L2*k__R
                  + L__R*k2
                  + 2*L*kur*k__R;
      }
      break;
    case TargetType::P8:
      evaluate_for_NLP_D(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        g[j]   += (2*G.L*G.dk__L + G.L__L*G.dk)*G.dk;
        g[j+1] += (2*G.L*G.dk__R + G.L__R*G.dk)*G.dk;
      }
      break;
    case TargetType::P9:
      evaluate_for_NLP_D(theta);
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G{ m_G1_vec[j] };
        real_type const L     { G.L };
        real_type const L__L  { G.L__L };
        real_type const L__R  { G.L__R };
        real_type const kur   { G.k0 };
        real_type const k__L  { G.k__L };
        real_type const k__R  { G.k__R };
        real_type const k2    { kur*kur };
        real_type const k3    { kur*k2 };
        real_type const dk    { G.dk };
        real_type const dk__L { G.dk__L };
        real_type const dk__R { G.dk__R };
        real_type const dk2   { dk*dk };
        real_type const dkL   { dk*L };
        real_type const A     { ( ( (dkL+4*kur)*dkL + 6*k2)*dkL + 4*k3) * dkL + dk2 + k2*k2 };
        real_type const B     { ( ( ( ( 3*kur + 0.8*dkL ) * dkL + 4*k2 ) * dkL +2*k3 ) * L + 2*dk ) * L };
        real_type const C     { ( ( ( dkL + 4*kur ) * dkL + 6*k2 ) * dkL + 4*k3 ) * L };
        g[j]   += A*L__L + B*dk__L + C*k__L;
        g[j+1] += A*L__R + B*dk__R + C*k__R;
      }
      break;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::constraints(
    real_type const theta[],
    real_type       c[]
  ) const {

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    evaluate_for_NLP( theta );

    for ( integer j{0}; j < ne1; ++j ) c[j] = m_G1_vec[j].k1-m_G1_vec[j+1].k0;

    switch (m_tt) {
    case TargetType::P1:
      c[ne1] = diff2pi( theta[0]  - m_theta_I );
      c[ne]  = diff2pi( theta[ne] - m_theta_F );
      break;
    case TargetType::P2:
      {
        auto const & G0{ m_G1_vec[0]   };
        auto const & GE{ m_G1_vec[ne1] };
        c[ne1] = GE.k1 - G0.k0;
        c[ne]  = diff2pi( theta[0] - theta[ne] );
      }
      break;
    default:
      break;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidSplineG2::jacobian_nnz() const {
    integer nnz{ 3*(m_npts-2) };
    switch (m_tt) {
    case TargetType::P1: nnz += 2; break;
    case TargetType::P2: nnz += 6; break;
    default:                       break;
    }
    return nnz;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::jacobian_pattern(
    integer ii[],
    integer jj[]
  ) const {

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    integer kk{0};
    for ( integer j{0}; j < ne1; ++j ) {
      ii[kk] = j; jj[kk] = j;   ++kk;
      ii[kk] = j; jj[kk] = j+1; ++kk;
      ii[kk] = j; jj[kk] = j+2; ++kk;
    }

    switch (m_tt) {
    case TargetType::P1:
      ii[kk] = ne1; jj[kk] = 0; ++kk;
      ii[kk] = ne;  jj[kk] = ne; // ++kk;
      break;
    case TargetType::P2:
      ii[kk] = ne1; jj[kk] = 0;   ++kk;
      ii[kk] = ne1; jj[kk] = 1;   ++kk;
      ii[kk] = ne1; jj[kk] = ne1; ++kk;
      ii[kk] = ne1; jj[kk] = ne;  ++kk;
      ii[kk] = ne;  jj[kk] = 0;   ++kk;
      ii[kk] = ne;  jj[kk] = ne; // ++kk;
      break;
    default:
      break;
    }

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::jacobian_pattern_matlab(
    real_type ii[],
    real_type jj[]
  ) const {

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    integer kk{0};
    for ( integer j{1}; j <= ne1; ++j ) {
      ii[kk] = j; jj[kk] = j;   ++kk;
      ii[kk] = j; jj[kk] = j+1; ++kk;
      ii[kk] = j; jj[kk] = j+2; ++kk;
    }

    switch (m_tt) {
    case TargetType::P1:
      ii[kk] = ne;     jj[kk] = 1;      ++kk;
      ii[kk] = m_npts; jj[kk] = m_npts; ++kk;
      break;
    case TargetType::P2:
      ii[kk] = ne;     jj[kk] = 1;      ++kk;
      ii[kk] = ne;     jj[kk] = 2;      ++kk;
      ii[kk] = ne;     jj[kk] = ne;     ++kk;
      ii[kk] = ne;     jj[kk] = m_npts; ++kk;
      ii[kk] = m_npts; jj[kk] = 1;      ++kk;
      ii[kk] = m_npts; jj[kk] = m_npts; ++kk;
      break;
    default:
      break;
    }

    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::jacobian(
    real_type const theta[],
    real_type       vals[]
  ) const {
    integer const ne1 { m_npts - 2 };

    evaluate_for_NLP_D( theta );

    integer kk{0};
    for ( integer j{0}; j < ne1; ++j ) {
      auto const & G  { m_G1_vec[j]   };
      auto const & G1 { m_G1_vec[j+1] };
      vals[kk++] =  G.k__L + G.dk__L * G.L + G.dk * G.L__L;
      vals[kk++] =  G.k__R + G.dk__R * G.L + G.dk * G.L__R - G1.k__L;
      vals[kk++] = -G1.k__R;
    }

    switch (m_tt) {
    case TargetType::P1:
      vals[kk++] = 1;
      vals[kk++] = 1;
      break;
    case TargetType::P2: {
        auto const & G0 { m_G1_vec[0]   };
        auto const & GE { m_G1_vec[ne1] };
        vals[kk++] = -G0.k__L;
        vals[kk++] = -G0.k__R;
        vals[kk++] = GE.k__L + GE.L__L * GE.dk + GE.L * GE.dk__L;
        vals[kk++] = GE.k__R + GE.L__R * GE.dk + GE.L * GE.dk__R;
        vals[kk++] = 1;
        vals[kk++] = -1;
      }
      break;
    default:
      break;
    }
    return true;
  }



  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidSplineG2::jacobian( real_type const theta[], Pipal::SparseMatrix<real_type> & J ) const {

    evaluate_for_NLP_D(theta);
    J.setZero();
    J.resize( this->numConstraints(), m_npts );
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(10*m_npts);

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    evaluate_for_NLP_D( theta );

    for ( integer j{0}; j < ne1; ++j ) {
      auto const & G  { m_G1_vec[j] };
      auto const & G1 { m_G1_vec[j+1] };
      real_type const A{ G.k__L + G.dk__L * G.L + G.dk * G.L__L };
      real_type const B{ G.k__R + G.dk__R * G.L + G.dk * G.L__R - G1.k__L };
      real_type const C{                                        - G1.k__R };
      triplets.emplace_back(j, j,   A );
      triplets.emplace_back(j, j+1, B );
      triplets.emplace_back(j, j+2, C );
    }

    switch (m_tt) {
    case TargetType::P1:
      triplets.emplace_back( ne1, 0,  1 );
      triplets.emplace_back( ne,  ne, 1 );
      break;
    case TargetType::P2:
      {
        auto const & G0 { m_G1_vec[0]   };
        auto const & GE { m_G1_vec[ne1] };
        triplets.emplace_back( ne1, 0,   -G0.k__L );
        triplets.emplace_back( ne1, 1,   -G0.k__R );
        triplets.emplace_back( ne1, ne1, GE.k__L + GE.L__L * GE.dk + GE.L * GE.dk__L );
        triplets.emplace_back( ne1, ne,  GE.k__R + GE.L__R * GE.dk + GE.L * GE.dk__R );
        triplets.emplace_back( ne,  0,   1 );
        triplets.emplace_back( ne,  ne, -1 );
      }
      break;
    default:
      break;
    }

    J.setFromTriplets(triplets.begin(), triplets.end());
    J.makeCompressed();
    return true;
  }


  bool
  ClothoidSplineG2::lagrangian_hessian(
    real_type const                  theta[],
    real_type const                  lambda[],
    Pipal::SparseMatrix<real_type> & H
  ) const {
    evaluate_for_NLP_DD(theta);
    H.setZero();
    H.resize( m_npts, m_npts );
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(10*m_npts);

    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };
    switch (m_tt) {
    case TargetType::P1:
    case TargetType::P2:
    case TargetType::P3:
      break;
    case TargetType::P4:
      {
        auto const & G { m_G1_vec[0]   };
        real_type const tmp{ 2*(G.dk__L*G.dk__R+G.dk*G.dk__LR) };
        triplets.emplace_back(0, 0, 2*(G.dk__L*G.dk__L+G.dk*G.dk__LL) );
        triplets.emplace_back(0, 1, tmp );
        triplets.emplace_back(1, 0, tmp );
        triplets.emplace_back(1, 1, 2*(G.dk__R*G.dk__R+G.dk*G.dk__RR) );
      }
      {
        auto const & G { m_G1_vec[ne1] };
        real_type const tmp { 2*(G.dk__L*G.dk__R+G.dk*G.dk__LR) };
        triplets.emplace_back(ne1, ne1, 2*(G.dk__L*G.dk__L+G.dk*G.dk__LL) );
        triplets.emplace_back(ne1, ne,  tmp );
        triplets.emplace_back(ne,  ne1, tmp );
        triplets.emplace_back(ne,  ne,  2*(G.dk__R*G.dk__R+G.dk*G.dk__RR) );
      }
      break;
    case TargetType::P5:
      {
        auto const & G { m_G1_vec[0]   };
        triplets.emplace_back(0, 0, G.L__LL );
        triplets.emplace_back(0, 1, G.L__LR );
        triplets.emplace_back(1, 0, G.L__LR );
        triplets.emplace_back(1, 1, G.L__RR );
      }
      {
        auto const & G { m_G1_vec[ne1] };
        triplets.emplace_back(ne1, ne1, G.L__LL );
        triplets.emplace_back(ne1, ne,  G.L__LR );
        triplets.emplace_back(ne,  ne1, G.L__LR );
        triplets.emplace_back(ne,  ne,  G.L__RR );
      }
      break;
    case TargetType::P6:
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G { m_G1_vec[j] };
        triplets.emplace_back( j,   j,   G.L__LL );
        triplets.emplace_back( j,   j+1, G.L__LR );
        triplets.emplace_back( j+1, j,   G.L__LR );
        triplets.emplace_back( j+1, j+1, G.L__RR );
      }
      break;
    case TargetType::P7:
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G { m_G1_vec[j] };
        real_type const L      { G.L };
        real_type const L__L   { G.L__L };
        real_type const L__R   { G.L__R };
        real_type const L__LL  { G.L__LL };
        real_type const L__LR  { G.L__LR };
        real_type const L__RR  { G.L__RR };
        real_type const k0     { G.k0 };
        real_type const k__L   { G.k__L };
        real_type const k__R   { G.k__R };
        real_type const k__LL  { G.k__LL };
        real_type const k__LR  { G.k__LR };
        real_type const k__RR  { G.k__RR };
        real_type const dk     { G.dk };
        real_type const dk__L  { G.dk__L };
        real_type const dk__R  { G.dk__R };
        real_type const dk__LL { G.dk__LL };
        real_type const dk__LR { G.dk__LR };
        real_type const dk__RR { G.dk__RR };
        real_type const t2  = dk__L*dk__L;
        real_type const t4  = L*L;
        real_type const t5  = t4*L;
        real_type const t11 = dk*dk;
        real_type const t15 = L__L*L__L;
        real_type const t18 = k__L*k__L;
        real_type const t22 = k0*k0;
        real_type const t55 = dk__R*dk__R;
        real_type const t65 = L__R*L__R;
        real_type const t68 = k__R*k__R;
        real_type const LL = 2.0/3.0*t5*(dk*dk__LL+t2)+t4*(4.0*L__L*dk*dk__L+t11*L__LL)+2.0*L*(k0*k__LL+t11*t15+t18)+t22*L__LL+4.0*L__L*k0*k__L;
        real_type const LR = 2.0/3.0*t5*(dk*dk__LR+dk__L*dk__R)+2.0*t4*dk*(L__L*dk__R+L__LR*dk/2.0+L__R*dk__L)
                            +2.0*L*(L__L*t11*L__R+k0*k__LR+k__R*k__L)+k0*(2.0*L__L*k__R+L__LR*k0+2.0*L__R*k__L);
        real_type const RR = 2.0/3.0*t5*(dk*dk__RR+t55)+t4*(4.0*L__R*dk*dk__R+t11*L__RR)+2.0*L*(k0*k__RR+t11*t65+t68)+t22*L__RR+4.0*L__R*k0*k__R;
        triplets.emplace_back( j,   j,   LL );
        triplets.emplace_back( j,   j+1, LR );
        triplets.emplace_back( j+1, j,   LR );
        triplets.emplace_back( j+1, j+1, RR );
      }
      break;
    case TargetType::P8:
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G { m_G1_vec[j] };
        real_type const L      { G.L };
        real_type const L__L   { G.L__L };
        real_type const L__R   { G.L__R };
        real_type const L__LL  { G.L__LL };
        real_type const L__LR  { G.L__LR };
        real_type const L__RR  { G.L__RR };
        real_type const dk     { G.dk };
        real_type const dk__L  { G.dk__L };
        real_type const dk__R  { G.dk__R };
        real_type const dk__LL { G.dk__LL };
        real_type const dk__LR { G.dk__LR };
        real_type const dk__RR { G.dk__RR };
        real_type const t1 = dk*L;
        real_type const t4 = dk__L*dk__L;
        real_type const t10 = dk*dk;
        real_type const t24 = dk__R*dk__R;
        real_type const LL = 4.0*L__L*dk*dk__L+2.0*t4*L+t10*L__LL+2.0*dk__LL*t1;
        real_type const LR = t10*L__LR+(2.0*L*dk__LR+2.0*L__L*dk__R+2.0*L__R*dk__L)*dk+2.0*L*dk__R*dk__L;
        real_type const RR = 4.0*L__R*dk*dk__R+2.0*t24*L+t10*L__RR+2.0*dk__RR*t1;

        triplets.emplace_back( j,   j,   LL );
        triplets.emplace_back( j,   j+1, LR );
        triplets.emplace_back( j+1, j,   LR );
        triplets.emplace_back( j+1, j+1, RR );
      }
      break;
    case TargetType::P9:
      for ( integer j{0}; j < ne; ++j ) {
        auto const & G { m_G1_vec[j] };
        real_type const L      { G.L };
        real_type const L__L   { G.L__L };
        real_type const L__R   { G.L__R };
        real_type const L__LL  { G.L__LL };
        real_type const L__LR  { G.L__LR };
        real_type const L__RR  { G.L__RR };
        real_type const k0     { G.k0 };
        real_type const k__L   { G.k__L };
        real_type const k__R   { G.k__R };
        real_type const k__LL  { G.k__LL };
        real_type const k__LR  { G.k__LR };
        real_type const k__RR  { G.k__RR };
        real_type const dk     { G.dk };
        real_type const dk__L  { G.dk__L };
        real_type const dk__R  { G.dk__R };
        real_type const dk__LL { G.dk__LL };
        real_type const dk__LR { G.dk__LR };
        real_type const dk__RR { G.dk__RR };
        real_type const t1 = dk*dk;
        real_type const t2 = t1*dk;
        real_type const t4 = dk__L*dk__L;
        real_type const t8 = L*L;
        real_type const t9 = t8*t8;
        real_type const t10 = t9*L;
        real_type const t14 = L__L*dk__L;
        real_type const t18 = dk__L*k__L;
        real_type const t29 = t1*t1;
        real_type const t30 = L__L*L__L;
        real_type const t33 = 2.0*L__L*k__L;
        real_type const t34 = k0*L__LL;
        real_type const t40 = k__L*k__L;
        real_type const t43 = k0*t18;
        real_type const t45 = k0*k0;
        real_type const t46 = t45*dk__LL;
        real_type const t51 = t8*L;
        real_type const t58 = 2.0*t14;
        real_type const t75 = t45*k0;
        real_type const t89 = t75*L__L;
        real_type const t94 = t45*t45;
        real_type const LL = 4.0/5.0*t10*(dk__LL*t2+3.0*t4*t1)+t9*dk*(L__LL*t2+t1*(8.0*t14+k__LL)
+3.0*dk*(dk__LL*k0+2.0*t18)+6.0*t4*k0)+4.0*t51*(t30*t29+t2*(t33+t34)+t1*(k0*(
6.0*t14+k__LL)+t40)+dk*(4.0*t43+t46)+t4*t45)+12.0*t8*k0*(t30*t2+t1*(t34/2.0+t33
)+dk*(k0*(t58+k__LL/2.0)+t40)+t46/6.0+t43)+2.0*L*(6.0*t1*t45*t30+dk*(12.0*k__L*
t45*L__L+2.0*t75*L__LL+dk__LL)+2.0*t75*(t58+k__LL)+6.0*t40*t45+t4)+t1*L__LL+4.0
*dk*(t89+dk__L)*L__L+t94*L__LL+8.0*k__L*t89;
        real_type const t107 = L__L*dk__R;
        real_type const t108 = L__R*dk__L;
        real_type const t112 = dk__L*k__R;
        real_type const t113 = dk__LR*k0;
        real_type const t114 = dk__R*k__L;
        real_type const t127 = L__L*k__R;
        real_type const t128 = L__LR*k0;
        real_type const t129 = L__R*k__L;
        real_type const t136 = k__R*k__L;
        real_type const t167 = L__R*L__L;
        real_type const LR = 4.0/5.0*t10*(3.0*dk__R*dk__L*t1+dk__LR*t2)+4.0*t9*dk*(L__LR*t2/4.0+
t1*(t107+t108+k__LR/4.0)+3.0/4.0*dk*(t112+t113+t114)+3.0/2.0*dk__R*k0*dk__L)+
4.0*t51*(L__L*L__R*t29+t2*(t127+t128+t129)+t1*(k0*(3.0*t107+3.0*t108+k__LR)+
t136)+2.0*dk*(t114+t113/2.0+t112)*k0+dk__L*dk__R*t45)+12.0*t8*(L__R*L__L*t2+t1*
(t128/2.0+t127+t129)+dk*(k0*(t107+t108+k__LR/2.0)+t136)+(t114+t113/3.0+t112)*k0
/2.0)*k0+2.0*L*(6.0*t1*t45*t167+dk*(2.0*L__LR*t75+6.0*t45*(t127+t129)+dk__LR)+
2.0*t75*(t107+t108+k__LR)+6.0*k__R*k__L*t45+dk__L*dk__R)+t1*L__LR+2.0*dk*(2.0*
t75*t167+t107+t108)+(4.0*t127+t128+4.0*t129)*t75;
        real_type const t199 = dk__R*dk__R;
        real_type const t206 = L__R*dk__R;
        real_type const t210 = dk__R*k__R;
        real_type const t221 = L__R*L__R;
        real_type const t224 = 2.0*L__R*k__R;
        real_type const t225 = k0*L__RR;
        real_type const t231 = k__R*k__R;
        real_type const t234 = k0*t210;
        real_type const t236 = t45*dk__RR;
        real_type const t247 = 2.0*t206;
        real_type const t277 = t75*L__R;
        real_type const RR = 4.0/5.0*t10*(dk__RR*t2+3.0*t199*t1)+t9*dk*(L__RR*t2+t1*(8.0*t206+
k__RR)+3.0*dk*(dk__RR*k0+2.0*t210)+6.0*t199*k0)+4.0*t51*(t221*t29+t2*(t224+t225
)+t1*(k0*(6.0*t206+k__RR)+t231)+dk*(4.0*t234+t236)+t199*t45)+12.0*t8*k0*(t221*
t2+t1*(t225/2.0+t224)+dk*(k0*(t247+k__RR/2.0)+t231)+t236/6.0+t234)+2.0*L*(6.0*
t1*t45*t221+dk*(12.0*k__R*t45*L__R+2.0*t75*L__RR+dk__RR)+2.0*t75*(t247+k__RR)+
6.0*t231*t45+t199)+t1*L__RR+4.0*dk*(t277+dk__R)*L__R+t94*L__RR+8.0*k__R*t277;

        triplets.emplace_back( j,   j,   LL );
        triplets.emplace_back( j,   j+1, LR );
        triplets.emplace_back( j+1, j,   LR );
        triplets.emplace_back( j+1, j+1, RR );
      }
      break;
    }
    for ( integer j{0}; j < ne1; ++j ) {
      auto const & G  { m_G1_vec[j] };
      auto const & G1 { m_G1_vec[j+1] };
      real_type const L      { G.L };
      real_type const L__L   { G.L__L };
      real_type const L__R   { G.L__R };
      real_type const L__LL  { G.L__LL };
      real_type const L__LR  { G.L__LR };
      real_type const L__RR  { G.L__RR };
      real_type const k__LL  { G.k__LL };
      real_type const k__LR  { G.k__LR };
      real_type const k__RR  { G.k__RR };
      real_type const k1__LL { G1.k__LL };
      real_type const k1__LR { G1.k__LR };
      real_type const k1__RR { G1.k__RR };
      real_type const dk     { G.dk };
      real_type const dk__L  { G.dk__L };
      real_type const dk__R  { G.dk__R };
      real_type const dk__LL { G.dk__LL };
      real_type const dk__LR { G.dk__LR };
      real_type const dk__RR { G.dk__RR };

      real_type m00 = dk__LL*L+2.0*L__L*dk__L+L__LL*dk+k__LL;
      real_type m01 = L*dk__LR+L__L*dk__R+L__LR*dk+L__R*dk__L+k__LR;
      real_type m02 = 0.0;
      real_type m11 = dk__RR*L+2.0*L__R*dk__R+L__RR*dk-k1__LL+k__RR;
      real_type m12 = -k1__LR;
      real_type m22 = -k1__RR;

      m00 *= -lambda[j];
      m01 *= -lambda[j];
      m02 *= -lambda[j];
      m11 *= -lambda[j];
      m12 *= -lambda[j];
      m22 *= -lambda[j];

      triplets.emplace_back( j,   j,   m00 );
      triplets.emplace_back( j,   j+1, m01 );
      triplets.emplace_back( j,   j+2, m02 );
      triplets.emplace_back( j+1, j,   m01 );
      triplets.emplace_back( j+1, j+1, m11 );
      triplets.emplace_back( j+1, j+2, m12 );
      triplets.emplace_back( j+2, j,   m02 );
      triplets.emplace_back( j+2, j+1, m12 );
      triplets.emplace_back( j+2, j+2, m22 );
    }

    switch (m_tt) {
    case TargetType::P1:
      //c[ne1] = diff2pi( theta[0]  - m_theta_I );
      //c[ne]  = diff2pi( theta[ne] - m_theta_F );
      break;
    case TargetType::P2:
      //c[ne1] = m_k1[ne1] - m_k0[0];
      //c[ne]  = diff2pi( theta[0] - theta[ne] );
      break;
    default:
      break;
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    H.makeCompressed();
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `ClothoidSplineG2` object
  //!
  //!  \param stream the output stream
  //!  \param c     an instance of `ClothoidSplineG2` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, ClothoidSplineG2 const & c ) {
    fmt::print( stream,
      "npts   = {}\n"
      "target = {}\n",
      c.m_npts, ClothoidSplineG2::to_string(c.m_tt)
    );
    return stream;
  }

}

// EOF: ClothoidListG2.cc
