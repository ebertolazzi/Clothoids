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
#include "Utils_SPSA.hh"


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
  static auto minmod = [] ( real_type a,  real_type b ) -> real_type {
    if (a * b <= 0) return 0;
    return (a > 0) ? std::min(a, b) : std::max(a, b);
  };

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
    opts.max_iter        = 100;
    opts.m               = 20;
    opts.iter_reset      = 50;
    opts.verbose         = true; // Disabilita verbose per output più pulito
    opts.use_projection  = true;
    opts.g_tol           = 1e-6;
    opts.f_tol           = 1e-12;
    opts.x_tol           = 1e-6;
    opts.step_max        = 10;
    opts.sty_min_factor  = 1e-12;
    opts.very_small_step = 1e-8;

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
      ClothoidList tmp1("temporary");
      ClothoidList tmp2("temporary");
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
          ok = tmp1.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vp{ target(tmp1) };
          X.coeffRef(i) = x_save - w_h * nx_save;
          Y.coeffRef(i) = y_save - w_h * ny_save;
          ok = tmp2.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vm{ target(tmp2) };
          grad->coeffRef(i) = (vp-vm)/(2*w_h); // minmod( vp - value, value - vm)/w_h;
        }
      }

      return value;
    };

    Utils::SPSAGradientEstimator<real_type>::Options opts2;
    opts2.repeats = 40;        // consigliato 10–40
    opts2.c_base  = w_h;       // spesso si usa la stessa scala del passo finito
    opts2.use_rademacher = true;
    
    Utils::SPSAGradientEstimator<real_type> spsa(opts2);

    auto TARGET_new = [
       &n, &X0, &Y0, &nx, &ny,
       &theta_init, &theta_end,
       &target,
       &spsa   // <--- aggiungi qui il tuo SPSAGradientEstimator
    ]( Vector const & offs, Vector * grad ) -> real_type {
      // funzione che costruisce la curva e ritorna il valore del tuo target
      auto eval_f = [&](Vector const & x)->real_type { Vector X(n), Y(n);
        X = X0.array() + x.array() * nx.array();
        Y = Y0.array() + x.array() * ny.array();
        ClothoidList tmp("temporary");
        bool ok = tmp.build_G2(n, X.data(), Y.data(), theta_init, theta_end);
        return target(tmp);
      };

      // se NON serve il gradiente → calcolo diretto
      if (grad == nullptr) return eval_f(offs);

      // Se serve il gradiente → usiamo SPSA
      // Definiamo un callback che rispetti l’interfaccia richiesta da SPSA:
      typename Utils::SPSAGradientEstimator<real_type>::Callback cb = [&](
        Vector const & x, Vector * g_unused
      )->real_type { return eval_f(x); };

      // stimiamo gradiente tramite SPSA
      grad->resize(n);
      spsa.estimate(offs, cb, *grad);

      // ritorniamo anche il valore in offs
      return eval_f(offs);
    };

    Vector offs0(n);
    offs0.setZero();

    // Inizializza le line search
    //Utils::StrongWolfeLineSearch<real_type> line_search;
    //Utils::HagerZhangLineSearch<real_type>  line_search;
    //Utils::MoreThuenteLineSearch<real_type> line_search;
    //Utils::WeakWolfeLineSearch<real_type> line_search;
    Utils::ArmijoLineSearch<real_type> line_search;
    
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
    opts.max_iter        = 100;
    opts.m               = 20;
    opts.iter_reset      = 50;
    opts.verbose         = true; // Disabilita verbose per output più pulito
    opts.use_projection  = true;
    opts.g_tol           = 1e-6;
    opts.f_tol           = 1e-12;
    opts.x_tol           = 1e-6;
    opts.step_max        = 10;
    opts.sty_min_factor  = 1e-12;
    opts.very_small_step = 1e-8;

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
      ClothoidList tmp1("temporary");
      ClothoidList tmp2("temporary");
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
          ok = tmp1.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vp{ target(tmp1) };
          X.coeffRef(i) = x_save - w_h * nx_save;
          Y.coeffRef(i) = y_save - w_h * ny_save;
          ok = tmp2.build_G2_cyclic( n, X.data(), Y.data() );
          real_type vm{ target(tmp2) };
          grad->coeffRef(i) = (vp-vm)/(2*w_h); // minmod( vp - value, value - vm)/w_h;
        }
      }

      return value;
    };

    Vector offs0(n);
    offs0.setZero();

    // Inizializza le line search
    //Utils::StrongWolfeLineSearch<real_type> line_search;
    //Utils::HagerZhangLineSearch<real_type>  line_search;
    //Utils::MoreThuenteLineSearch<real_type> line_search;
    //Utils::WeakWolfeLineSearch<real_type> line_search;
    Utils::ArmijoLineSearch<real_type> line_search;
    
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
