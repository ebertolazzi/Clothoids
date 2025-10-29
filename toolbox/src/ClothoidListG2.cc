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

#include "Utils/3rd/Eigen/SparseCholesky"

#define PIPAL_EIGEN_EXTERNAL
#include "Pipal.hh"

#include <cfloat>

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
  ClothoidSplineG2::build(
    real_type const xvec[],
    real_type const yvec[],
    integer   const n
  ) {
    m_npts = n;

    m_x . resize( n );
    m_y . resize( n );

    m_L  . resize( n-1 );
    m_k0 . resize( n-1 );
    m_k1 . resize( n-1 );
    m_dk . resize( n-1 );

    m_L__L  . resize( n-1 );
    m_L__R  . resize( n-1 );
    m_L__LL . resize( n-1 );
    m_L__LR . resize( n-1 );
    m_L__RR . resize( n-1 );

    m_k__L  . resize( n-1 );
    m_k__R  . resize( n-1 );
    m_k__LL . resize( n-1 );
    m_k__LR . resize( n-1 );
    m_k__RR . resize( n-1 );

    m_dk__L . resize( n-1 );
    m_dk__R . resize( n-1 );
    m_dk__LL . resize( n-1 );
    m_dk__LR . resize( n-1 );
    m_dk__RR . resize( n-1 );

    std::copy_n( xvec, n, m_x.data() );
    std::copy_n( yvec, n, m_y.data() );
  }


  void
  ClothoidSplineG2::evaluate_for_NLP( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP temporary cc"};
    integer const ne{ m_npts - 1 };
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1( m_x(j),   m_y(j),   theta[j],
                   m_x(j+1), m_y(j+1), theta[j+1] );
      m_k0[j] = cc.kappa_begin();
      m_dk[j] = cc.dkappa();
      m_L[j]  = cc.length();
      m_k1[j] = cc.kappa_end(); // m_k0[j]+m_dk[j]*m_L[j];
    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_D( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_D temporary cc"};
    integer const ne{ m_npts - 1 };
    real_type L_D[2], k_D[2], dk_D[2];
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1_D( m_x(j),   m_y(j),   theta[j],
                     m_x(j+1), m_y(j+1), theta[j+1],
                     L_D, k_D, dk_D );
      m_k0[j] = cc.kappa_begin();
      m_dk[j] = cc.dkappa();
      m_L[j]  = cc.length();
      m_k1[j] = cc.kappa_end();

      m_L__L[j] = L_D[0];
      m_L__R[j] = L_D[1];

      m_k__L[j] = k_D[0];
      m_k__R[j] = k_D[1];

      m_dk__L[j] = dk_D[0];
      m_dk__R[j] = dk_D[1];

    }
  }

  void
  ClothoidSplineG2::evaluate_for_NLP_DD( real_type const theta[] ) const {
    ClothoidCurve cc{"ClothoidSplineG2::evaluate_for_NLP_DD temporary cc"};
    integer const ne{ m_npts - 1 };
    real_type L_D[2], k_D[2], dk_D[2];
    real_type L_DD[3], k_DD[3], dk_DD[3];
    for ( integer j{0}; j < ne; ++j ) {
      cc.build_G1_DD( m_x(j),   m_y(j),   theta[j],
                      m_x(j+1), m_y(j+1), theta[j+1],
                      L_D, k_D, dk_D, L_DD, k_DD, dk_DD );
      m_k0[j] = cc.kappa_begin();
      m_dk[j] = cc.dkappa();
      m_L[j]  = cc.length();
      m_k1[j] = cc.kappa_end();

      m_L__L[j]  = L_D[0];
      m_L__R[j]  = L_D[1];
      m_L__LL[j] = L_DD[0];
      m_L__LR[j] = L_DD[1];
      m_L__RR[j] = L_DD[2];

      m_k__L[j]  = k_D[0];
      m_k__R[j]  = k_D[1];
      m_k__LL[j] = k_DD[0];
      m_k__LR[j] = k_DD[1];
      m_k__RR[j] = k_DD[2];

      m_dk__L[j]  = dk_D[0];
      m_dk__R[j]  = dk_D[1];
      m_dk__LL[j] = dk_DD[0];
      m_dk__LR[j] = dk_DD[1];
      m_dk__RR[j] = dk_DD[2];

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
    Eigen::Vector<real_type,Eigen::Dynamic> omega( m_npts );
    Eigen::Vector<real_type,Eigen::Dynamic> len( m_npts );
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
    using Indices      = Pipal::Indices;
  
    build( xvec, yvec, n );
    
    Vector theta_guess( m_npts ), theta_min( m_npts ), theta_max( m_npts ), theta_sol( m_npts );
    this->guess( theta_guess.data(), theta_min.data(), theta_max.data() );
    
    integer nnz{ this->jacobian_nnz() };
    Vector  V( nnz );
    Indices I( nnz ), J( nnz );
    
    this->jacobian_pattern( I.data(), J.data() );

    Pipal::Solver<real_type> solver("ClothoidSplineG2::build",
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
      [this,&V,&I,&J,&nnz] ( Vector const & theta, SparseMatrix & out ) -> bool {
        bool ok{ this->jacobian( theta.data(), V.data() ) };
        out.resize( this->numConstraints(), this->numTheta() );
        out.setZero();
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(nnz);
        for ( integer i{0}; i < nnz; ++i ) triplets.emplace_back(I(i),J(i),V(i));  // solo elementi diagonali a 0
        out.setFromTriplets(triplets.begin(), triplets.end());
        return ok;
      },
      
      // Hessian of the Lagrangian
      [this] ( Vector const &, Vector const &, SparseMatrix & out ) -> bool {
        out.resize(m_npts,m_npts);
        out.setZero();
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(m_npts);
        for ( integer i{0}; i < m_npts; ++i ) triplets.emplace_back(i,i,0);  // solo elementi diagonali a 0
        out.setFromTriplets(triplets.begin(), triplets.end());
        return true;
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
    solver.algorithm(Pipal::Algorithm::CONSERVATIVE);
    solver.verbose_mode(VERBOSE);
    solver.tolerance(SOLVER_TOLERANCE);
    solver.max_iterations(MAX_ITERATIONS);
    bool ok{ solver.optimize( theta_guess, theta_sol ) };
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
  ClothoidSplineG2::objective(
    real_type const theta[],
    real_type     & f
  ) const {
    ClothoidCurve cL{"ClothoidSplineG2::objective temporary cL"};
    ClothoidCurve cR{"ClothoidSplineG2::objective temporary cR"};
    ClothoidCurve c{"ClothoidSplineG2::objective temporary c"};
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
    case TargetType::P4:
      cL.build_G1( m_x[0],   m_y[0],   theta[0],   m_x[1],  m_y[1],  theta[1] );
      cR.build_G1( m_x[ne1], m_y[ne1], theta[ne1], m_x[ne], m_y[ne], theta[ne] );
      { real_type const dk_L{ cL.dkappa() };
        real_type const dk_R{ cR.dkappa() };
        f = dk_L*dk_L+dk_R*dk_R;
      }
      break;
    case TargetType::P5:
      cL.build_G1( m_x[0],   m_y[0],   theta[0],   m_x[1],  m_y[1],  theta[1] );
      cR.build_G1( m_x[ne1], m_y[ne1], theta[ne1], m_x[ne], m_y[ne], theta[ne] );
      f = cL.length()+cR.length();
      break;
    case TargetType::P6:
      f = 0;
      for ( integer j{0}; j < ne; ++j ) {
        c.build_G1( m_x[j], m_y[j], theta[j], m_x[j+1], m_y[j+1], theta[j+1] );
        f += c.length();
      }
      break;
    case TargetType::P7:
      f = 0;
      for ( integer j{0}; j < ne; ++j ) {
        c.build_G1( m_x[j], m_y[j], theta[j], m_x[j+1], m_y[j+1], theta[j+1] );
        real_type const Len  { c.length() };
        real_type const kur  { c.kappa_begin() };
        real_type const dkur { c.dkappa() };
        f = f + Len * ( Len * ( dkur*( (dkur*Len)/3 + kur) ) + kur*kur );
      }
      break;
    case TargetType::P8:
      f = 0;
      for ( integer j{0}; j < ne; ++j ) {
        c.build_G1( m_x[j], m_y[j], theta[j], m_x[j+1], m_y[j+1], theta[j+1] );
        real_type const Len  { c.length() };
        real_type const dkur { c.dkappa() };
        f += Len*dkur*dkur;
      }
      break;
    case TargetType::P9:
      f = 0;
      for ( integer j{0}; j < ne; ++j ) {
        c.build_G1( m_x[j], m_y[j], theta[j], m_x[j+1], m_y[j+1], theta[j+1] );
        real_type const Len  { c.length() };
        real_type const kur  { c.kappa_begin() };
        real_type const k2   { kur*kur };
        real_type const k3   { k2*kur };
        real_type const k4   { k2*k2 };
        real_type const dkur { c.dkappa() };
        real_type const dk2  { dkur*dkur };
        real_type const dk3  { dkur*dk2 };
        f += (k4+dk2+(2*k3*dkur+(2*k2*dk2+(dk3*(kur+dkur*Len/5))*Len)*Len)*Len)*Len;
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
    ClothoidCurve cL{"ClothoidSplineG2::objective temporary cL"};
    ClothoidCurve cR{"ClothoidSplineG2::objective temporary cR"};
    ClothoidCurve c{"ClothoidSplineG2::objective temporary c"};
    real_type     LL_D[2], kL_D[2], dkL_D[2];
    real_type     LR_D[2], kR_D[2], dkR_D[2];
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
      cL.build_G1_D(
        m_x[0], m_y[0], theta[0],
        m_x[1], m_y[1], theta[1],
        LL_D, kL_D, dkL_D
      );
      cR.build_G1_D(
        m_x[ne1], m_y[ne1], theta[ne1],
        m_x[ne],  m_y[ne],  theta[ne],
        LR_D, kR_D, dkR_D
      );
      {
        real_type const dkL { cL.dkappa() };
        real_type const dkR { cR.dkappa() };
        g[0]   = 2*dkL*dkL_D[0];
        g[1]   = 2*dkL*dkL_D[1];
        g[ne1] = 2*dkR*dkR_D[0];
        g[ne]  = 2*dkR*dkR_D[1];
      }
      break;
    case TargetType::P5:
      cL.build_G1_D(
        m_x[0], m_y[0], theta[0],
        m_x[1], m_y[1], theta[1],
        LL_D, kL_D, dkL_D
      );
      cR.build_G1_D(
        m_x[ne1], m_y[ne1], theta[ne1],
        m_x[ne],  m_y[ne],  theta[ne],
        LR_D, kR_D, dkR_D
      );
      g[0]   = LL_D[0];
      g[1]   = LL_D[1];
      g[ne1] = LR_D[0];
      g[ne]  = LR_D[1];
      break;
    case TargetType::P6:
      for ( integer j{0}; j < ne; ++j ) {
        real_type L_D[2], k_D[2], dk_D[2];
        c.build_G1_D(
          m_x[j],   m_y[j],   theta[j],
          m_x[j+1], m_y[j+1], theta[j+1],
          L_D, k_D, dk_D
        );
        g[j]   += L_D[0];
        g[j+1] += L_D[1];
      }
      break;
    case TargetType::P7:
      for ( integer j{0}; j < ne; ++j ) {
        real_type L_D[2], k_D[2], dk_D[2];
        c.build_G1_D(
          m_x[j],   m_y[j],   theta[j],
          m_x[j+1], m_y[j+1], theta[j+1],
          L_D, k_D, dk_D
        );
        real_type const Len  { c.length() };
        real_type const L2   { Len*Len };
        real_type const L3   { Len*L2 };
        real_type const kur  { c.kappa_begin() };
        real_type const k2   { kur*kur };
        real_type const dkur { c.dkappa() };
        real_type const dk2  { dkur*dkur };
        g[j]   += 2*(dkur*dk_D[0]*L3)/3
                  + (dk2*L2*L_D[0])
                  + dk_D[0]*L2*kur
                  + 2*dkur*Len*L_D[0]*kur
                  + dkur*L2*k_D[0]
                  + L_D[0]*k2
                  + 2*Len*kur*k_D[0];
        g[j+1] += 2*(dkur*dk_D[1]*L3)/3
                  + (dk2*L2*L_D[1])
                  + dk_D[1]*L2*kur
                  + 2*dkur*Len*L_D[1]*kur
                  + dkur*L2*k_D[1]
                  + L_D[1]*k2
                  + 2*Len*kur*k_D[1];
      }
      break;
    case TargetType::P8:
      for ( integer j{0}; j < ne; ++j ) {
        real_type L_D[2], k_D[2], dk_D[2];
        c.build_G1_D(
          m_x[j],   m_y[j],   theta[j],
          m_x[j+1], m_y[j+1], theta[j+1],
          L_D, k_D, dk_D
        );
        real_type const Len  { c.length() };
        real_type const dkur { c.dkappa() };
        g[j]   += (2*Len*dk_D[0] + L_D[0]*dkur)*dkur;
        g[j+1] += (2*Len*dk_D[1] + L_D[1]*dkur)*dkur;
      }
      break;
    case TargetType::P9:
      for ( integer j{0}; j < ne; ++j ) {
        real_type L_D[2], k_D[2], dk_D[2];
        c.build_G1_D(
          m_x[j],   m_y[j],   theta[j],
          m_x[j+1], m_y[j+1], theta[j+1],
          L_D, k_D, dk_D
        );
        real_type const Len  { c.length() };
        real_type const kur  { c.kappa_begin() };
        real_type const k2   { kur*kur };
        real_type const k3   { kur*k2 };
        real_type const dkur { c.dkappa() };
        real_type const dk2  { dkur*dkur };
        real_type const dkL  { dkur*Len };
        real_type const A    { ( ( (dkL+4*kur)*dkL + 6*k2)*dkL + 4*k3) * dkL + dk2 + k2*k2 };
        real_type const B    { ( ( ( ( 3*kur + 0.8*dkL ) * dkL + 4*k2 ) * dkL +2*k3 ) * Len + 2*dkur ) * Len };
        real_type const C    { ( ( ( dkL + 4*kur ) * dkL + 6*k2 ) * dkL + 4*k3 ) * Len };
        g[j]   += A*L_D[0] + B*dk_D[0] + C*k_D[0];
        g[j+1] += A*L_D[1] + B*dk_D[1] + C*k_D[1];
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

    for ( integer j{0}; j < ne1; ++j ) c[j] = m_k1[j]-m_k0[j+1];

    switch (m_tt) {
    case TargetType::P1:
      c[ne1] = diff2pi( theta[0]  - m_theta_I );
      c[ne]  = diff2pi( theta[ne] - m_theta_F );
      break;
    case TargetType::P2:
      c[ne1] = m_k1[ne1] - m_k0[0];
      c[ne]  = diff2pi( theta[0] - theta[ne] );
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
    ClothoidCurve cc{"ClothoidSplineG2::jacobian_pattern temporary cc"};
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
    ClothoidCurve cc{"ClothoidSplineG2::jacobian_pattern_matlab temporary cc"};
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
    integer const ne  { m_npts - 1 };
    integer const ne1 { m_npts - 2 };

    evaluate_for_NLP_D( theta );

    integer kk{0};
    for ( integer j{0}; j < ne1; ++j ) {
      vals[kk++] =  m_k__L[j] + m_dk__L[j]*m_L[j] + m_dk[j]*m_L__L[j];
      vals[kk++] =  m_k__R[j] + m_dk__R[j]*m_L[j] + m_dk[j]*m_L__R[j] - m_k__L[j+1];
      vals[kk++] = -m_k__R[j+1];
    }

    switch (m_tt) {
    case TargetType::P1:
      vals[kk++] = 1;
      vals[kk++] = 1;
      break;
    case TargetType::P2:
      vals[kk++] = -m_k__L[0];
      vals[kk++] = -m_k__R[0];
      vals[kk++] = m_k__L[ne1]+m_L__L[ne1]*m_dk[ne1]+m_L[ne1]*m_dk__L[ne1];
      vals[kk++] = m_k__R[ne1]+m_L__R[ne1]*m_dk[ne1]+m_L[ne1]*m_dk__R[ne1];
      vals[kk++] = 1;
      vals[kk++] = -1;
      break;
    default:
      break;
    }
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
