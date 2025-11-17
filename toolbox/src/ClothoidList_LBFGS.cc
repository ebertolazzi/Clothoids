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
#include "Utils_LBFGS.hh"


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

  static real_type const h_fraction{ 1e-3 };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G2_with_target(
    integer   const n,
    real_type const x[],
    real_type const y[],
    real_type const theta[],
    real_type const w_min[],
    real_type const w_max[],
    real_type const theta_init,
    real_type const theta_end,
    std::function<real_type( ClothoidList const & lst )> const & target
  ) {
    using Vector = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    Utils::LBFGS_minimizer<real_type>::Options opts;
    opts.max_iter       = 100;
    opts.m              = 20;
    opts.verbose        = true; // Disabilita verbose per output più pulito
    opts.use_projection = true;
    opts.g_tol          = 1e-4;
    opts.f_tol          = 1e-8;
    opts.x_tol          = 1e-6;

    Utils::LBFGS_minimizer<real_type> minimizer(opts);
    minimizer.set_bounds( n, w_min, w_max );
    
    Vector nx( n ), ny( n );
    Eigen::Map<Vector const> THETA( theta, n );
    Eigen::Map<Vector const> X0( x, n );
    Eigen::Map<Vector const> Y0( y, n );
    Eigen::Map<Vector const> WMIN( w_min, n );
    Eigen::Map<Vector const> WMAX( w_max, n );
    nx = -THETA.array().sin();
    ny =  THETA.array().cos();
    real_type w_h{ (WMAX-WMIN).maxCoeff()*h_fraction };

    auto TARGET = [ &n, &X0, &Y0, &nx, &ny, &theta_init, &theta_end, &target, &w_h](
      Vector const & offs, Vector * grad
    ) -> real_type {
      // costruisco punti offsettati
      Vector X( n ), Y( n );
      X = X0.array() + offs.array() * nx.array();
      Y = Y0.array() + offs.array() * ny.array();
      // costruisco clotoide G2
      ClothoidList tmp("temporary");
      bool ok{ tmp.build_G2( n, X.data(), Y.data(), theta_init, theta_end ) };

      real_type value{ target(tmp) };
      if ( grad != nullptr ) {
        for ( integer i{0}; i < n; ++i ) {
          real_type x_save{ X.coeff(i) };
          real_type y_save{ Y.coeff(i) };
          real_type nx_save{ nx.coeff(i) };
          real_type ny_save{ ny.coeff(i) };
          X.coeffRef(i) = x_save + w_h * nx_save;
          Y.coeffRef(i) = y_save + w_h * ny_save;
          ok = tmp.build_G2( n, X.data(), Y.data(), theta_init, theta_end );
          real_type vp{ target(tmp) };
          X.coeffRef(i) = x_save - w_h * nx_save;
          Y.coeffRef(i) = y_save - w_h * ny_save;
          ok = tmp.build_G2( n, X.data(), Y.data(), theta_init, theta_end );
          real_type vm{ target(tmp) };
          grad->coeffRef(i) = (vp-vm)/(2*w_h);
        }
      }

      return value;
    };

    Vector offs0(n);
    offs0.setZero();

    // Inizializza le line search
    //Utils::StrongWolfeLineSearch<real_type> line_search;
    //Utils::HagerZhangLineSearch<real_type>  line_search;
    Utils::MoreThuenteLineSearch<real_type> line_search;
    
    auto iter_data = minimizer.minimize( offs0, TARGET, line_search );

    bool ok{ iter_data.status == Utils::LBFGS_minimizer<real_type>::Status::CONVERGED ||
             iter_data.status == Utils::LBFGS_minimizer<real_type>::Status::GRADIENT_TOO_SMALL };
    
    //if ( ok ) {
      Vector const & offs{ minimizer.solution() };
      // costruisco punti offsettati
      Vector X( n ), Y( n );
      X = X0.array() + offs.array() * nx.array();
      Y = Y0.array() + offs.array() * ny.array();
      // costruisco clotoide G2
      ok = this->build_G2( n, X.data(), Y.data(), theta_init, theta_end );
    //}
    return ok;
  }

  bool
  ClothoidList::build_G2_cyclic_with_target(
    integer   const n,
    real_type const x[],
    real_type const y[],
    real_type const theta[],
    real_type const w_min[],
    real_type const w_max[],
    std::function<real_type(ClothoidList const & lst)> const & target
  ) {
  
    using Vector = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    Utils::LBFGS_minimizer<real_type>::Options opts;
    opts.max_iter       = 100;
    opts.m              = 20;
    opts.verbose        = true; // Disabilita verbose per output più pulito
    opts.use_projection = true;
    opts.g_tol          = 1e-4;
    opts.f_tol          = 1e-8;
    opts.x_tol          = 1e-6;

    Utils::LBFGS_minimizer<real_type> minimizer(opts);
    minimizer.set_bounds( n, w_min, w_max );
    
    Vector nx( n ), ny( n );
    Eigen::Map<Vector const> THETA( theta, n );
    Eigen::Map<Vector const> X0( x, n );
    Eigen::Map<Vector const> Y0( y, n );
    Eigen::Map<Vector const> WMIN( w_min, n );
    Eigen::Map<Vector const> WMAX( w_max, n );
    nx = -THETA.array().sin();
    ny = THETA.array().cos();
    real_type w_h{ (WMAX-WMIN).maxCoeff()*h_fraction };

    auto TARGET = [ &n, &X0, &Y0, &nx, &ny, &target, &w_h]( Vector const & offs, Vector * grad ) -> real_type {
      // costruisco punti offsettati
      Vector X( n ), Y( n );
      X = X0.array() + offs.array() * nx.array();
      Y = Y0.array() + offs.array() * ny.array();
      // costruisco clotoide G2
      ClothoidList tmp("temporary");
      bool ok{ tmp.build_G2_cyclic( n, X.data(), Y.data() ) };

      real_type value{ target(tmp) };
      if ( grad != nullptr ) {
        for ( integer i{0}; i < n; ++i ) {
          real_type x_save{ X.coeff(i) };
          real_type y_save{ Y.coeff(i) };
          real_type nx_save{ nx.coeff(i) };
          real_type ny_save{ ny.coeff(i) };
          X.coeffRef(i) = x_save + w_h * nx_save;
          Y.coeffRef(i) = y_save + w_h * ny_save;
          ok = tmp.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vp{ target(tmp) };
          X.coeffRef(i) = x_save - w_h * nx_save;
          Y.coeffRef(i) = y_save - w_h * ny_save;
          ok = tmp.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vm{ target(tmp) };
          grad->coeffRef(i) = (vp-vm)/(2*w_h);
        }
      }

      return value;
    };

    Vector offs0(n);
    offs0.setZero();

    // Inizializza le line search
    //Utils::StrongWolfeLineSearch<real_type> line_search;
    //Utils::HagerZhangLineSearch<real_type>  line_search;
    Utils::MoreThuenteLineSearch<real_type> line_search;
    
    auto iter_data = minimizer.minimize( offs0, TARGET, line_search );

    bool ok{ iter_data.status == Utils::LBFGS_minimizer<real_type>::Status::CONVERGED ||
             iter_data.status == Utils::LBFGS_minimizer<real_type>::Status::GRADIENT_TOO_SMALL };
    
    //if ( ok ) {
      Vector const & offs{ minimizer.solution() };
      // costruisco punti offsettati
      Vector X( n ), Y( n );
      X = X0.array() + offs.array() * nx.array();
      Y = Y0.array() + offs.array() * ny.array();
      // costruisco clotoide G2
      ok = this->build_G2_cyclic( n, X.data(), Y.data() );
    //}
    return ok;
  }

}

// EOF: ClothoidList_LBFGS.cc
