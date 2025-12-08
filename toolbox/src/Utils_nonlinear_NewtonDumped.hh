// File: Utils_NewtonDumped.hh

/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022-2025                                                 |
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

#pragma once

#ifndef UTILS_NEWTON_DUMPED_dot_HH
#define UTILS_NEWTON_DUMPED_dot_HH

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"
#include "Utils_nonlinear_system.hh"
#include "Utils_pseudoinverse.hh"

namespace Utils
{

  /**
   * @brief Newton's method with multiple damping strategies for solving
   * nonlinear systems
   *
   * This class implements various damping strategies for Newton's method to
   * solve nonlinear systems of equations F(x) = 0. The available strategies
   * include:
   *
   * - **DEUFLHARD**: Multiplicative damping with Armijo-type condition
   * - **L2_CLASSIC**: Classical Levenberg-Marquardt style L2 regularization
   * - **L2_ADAPTIVE**: Adaptive L2 regularization with trust region management
   * - **L2_HYBRID**: Hybrid strategy switching between L2 and Deuflhard
   * - **DOGLEG**: Trust region dogleg method combining steepest descent and
   * Gauss-Newton
   * - **WOLFE_LINE_SEARCH**: Line search satisfying strong Wolfe conditions
   * - **CUBIC_REGULARIZATION**: Adaptive Regularization by Cubics (ARC) method
   * - **BACKTRACKING_QUADRATIC**: Quadratic interpolation backtracking line
   * search
   * - **BANK_ROSE**: Bank and Rose (1981) affine-invariant damping strategy
   * - **GRIEWANK**: Griewank (1980) adaptive damping with directional
   * derivative
   * - **FILTER_METHODS**: Filter-based acceptance mechanism for constrained
   * problems
   * - **CUBIC_TRUST_REGION**: Trust region with cubic models
   *
   * @section usage_sec Usage Example
   * @code{.cpp}
   * NewtonDumped solver;
   * solver.set_tolerance(1e-8);
   * solver.set_max_iterations(100);
   * solver.set_damping_strategy_wolfe();
   * solver.set_verbose(true);
   *
   * Vector x_initial(n);
   * // ... initialize x_initial ...
   *
   * if (solver.solve(my_system, x_initial)) {
   *   std::cout << "Converged in " << solver.get_num_iterations() << "
   * iterations\n";
   * }
   * @endcode
   *
   * @see NonlinearSystem for the system interface requirements
   */
  class NewtonDumped
  {
  public:
    /**
     * @brief Enumeration of available damping strategies
     */
    enum DampingStrategy
    {
      DEUFLHARD,               ///< Multiplicative Deuflhard damping
      L2_CLASSIC,              ///< Classical L2 norm damping (Levenberg-Marquardt)
      L2_ADAPTIVE,             ///< Adaptive L2 with trust region
      L2_HYBRID,               ///< Hybrid Deuflhard-L2 strategy
      DOGLEG,                  ///< Trust region dogleg method
      WOLFE_LINE_SEARCH,       ///< Line search with Wolfe conditions
      CUBIC_REGULARIZATION,    ///< Adaptive cubic regularization (ARC)
      BACKTRACKING_QUADRATIC,  ///< Quadratic/cubic backtracking
      BANK_ROSE,               ///< Bank and Rose (1981) strategy
      GRIEWANK,                ///< Griewank (1980) strategy
      FILTER_METHODS,          ///< Filter methods for nonlinear equations
      CUBIC_TRUST_REGION       ///< Trust region with cubic models
    };

    using integer      = Eigen::Index;                                 ///< Integer type for indexing
    using real_type    = double;                                       ///< Real number type
    using Vector       = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;  ///< Vector type
    using SparseMatrix = Eigen::SparseMatrix<real_type>;               ///< Sparse matrix type

  private:
    // =========================================================================
    // SOLVER CONFIGURATION
    // =========================================================================

    DampingStrategy m_damping_strategy = DEUFLHARD;  ///< Selected damping strategy

    real_type m_tolerance              = 1e-8;  ///< Stopping tolerance on ||f||
    integer   m_max_iterations         = 100;   ///< Maximum Newton iterations
    integer   m_max_damping_iterations = 20;    ///< Maximum damping tries per iteration

    // =========================================================================
    // DEUFLHARD PARAMETERS
    // =========================================================================

    real_type m_min_lambda          = 1e-8;  ///< Minimum allowed damping factor
    real_type m_lambda_factor       = 0.5;   ///< Reduction factor for lambda
    real_type m_sufficient_decrease = 0.01;  ///< Armijo-like constant for sufficient decrease

    // =========================================================================
    // L2 DAMPING PARAMETERS
    // =========================================================================

    real_type m_mu_init               = 0.01;   ///< Initial L2 damping parameter
    real_type m_mu_min                = 1e-8;   ///< Minimum mu
    real_type m_mu_max                = 1e4;    ///< Maximum mu
    real_type m_mu_increase_factor    = 4.0;    ///< Factor to increase mu on failure
    real_type m_mu_decrease_factor    = 0.1;    ///< Factor to decrease mu on success
    real_type m_trust_region_radius   = 1.0;    ///< Initial trust region radius
    real_type m_trust_region_min      = 1e-6;   ///< Minimum trust region radius
    real_type m_trust_region_max      = 100.0;  ///< Maximum trust region radius
    real_type m_acceptance_ratio_good = 0.8;    ///< Threshold for good steps
    real_type m_acceptance_ratio_bad  = 0.2;    ///< Threshold for rejecting steps

    // =========================================================================
    // BANK & ROSE PARAMETERS
    // =========================================================================

    real_type m_bank_rose_alpha     = 0.5;   ///< Damping reduction factor
    real_type m_bank_rose_beta      = 0.1;   ///< Sufficient decrease constant
    real_type m_bank_rose_gamma     = 0.9;   ///< Contraction factor
    real_type m_bank_rose_theta_min = 1e-4;  ///< Minimum damping factor
    real_type m_bank_rose_theta_max = 1.0;   ///< Maximum damping factor

    // =========================================================================
    // GRIEWANK PARAMETERS
    // =========================================================================

    real_type m_griewank_eta   = 0.1;   ///< Acceptance threshold
    real_type m_griewank_omega = 0.5;   ///< Reduction factor
    real_type m_griewank_tau   = 1e-4;  ///< Minimum step size
    real_type m_griewank_zeta  = 0.9;   ///< Contraction factor

    // =========================================================================
    // FILTER METHOD PARAMETERS
    // =========================================================================

    real_type m_filter_theta_min   = 1e-6;  ///< Minimum constraint violation
    real_type m_filter_gamma_theta = 0.01;  ///< Filter parameter for theta
    real_type m_filter_gamma_f     = 0.01;  ///< Filter parameter for f
    real_type m_filter_alpha       = 0.5;   ///< Backtracking factor
    real_type m_filter_beta        = 0.8;   ///< Margin for filter acceptance

    // =========================================================================
    // CUBIC TRUST REGION PARAMETERS
    // =========================================================================

    real_type m_ctr_delta     = 1.0;    ///< Trust region radius
    real_type m_ctr_delta_min = 1e-6;   ///< Minimum trust region radius
    real_type m_ctr_delta_max = 100.0;  ///< Maximum trust region radius
    real_type m_ctr_eta1      = 0.1;    ///< Unsuccessful step threshold
    real_type m_ctr_eta2      = 0.75;   ///< Very successful step threshold
    real_type m_ctr_gamma1    = 0.5;    ///< Shrink factor
    real_type m_ctr_gamma2    = 2.0;    ///< Expand factor
    real_type m_ctr_sigma     = 0.1;    ///< Cubic regularization parameter

    // =========================================================================
    // DOGLEG PARAMETERS
    // =========================================================================

    real_type m_dogleg_delta     = 1.0;    ///< Trust region radius
    real_type m_dogleg_delta_min = 1e-6;   ///< Minimum trust region radius
    real_type m_dogleg_delta_max = 100.0;  ///< Maximum trust region radius
    real_type m_dogleg_eta1      = 0.1;    ///< Unsuccessful step threshold
    real_type m_dogleg_eta2      = 0.75;   ///< Very successful step threshold
    real_type m_dogleg_gamma1    = 0.5;    ///< Shrink factor
    real_type m_dogleg_gamma2    = 2.0;    ///< Expand factor

    // =========================================================================
    // WOLFE LINE SEARCH PARAMETERS
    // =========================================================================

    real_type m_wolfe_c1         = 1e-4;  ///< Armijo condition constant (sufficient decrease)
    real_type m_wolfe_c2         = 0.9;   ///< Curvature condition constant
    real_type m_wolfe_alpha_init = 1.0;   ///< Initial step size
    real_type m_wolfe_alpha_min  = 1e-8;  ///< Minimum step size
    real_type m_wolfe_alpha_max  = 10.0;  ///< Maximum step size
    real_type m_wolfe_rho        = 0.5;   ///< Backtracking factor

    // =========================================================================
    // CUBIC REGULARIZATION (ARC) PARAMETERS
    // =========================================================================

    real_type m_cubic_sigma          = 0.1;   ///< Cubic regularization parameter
    real_type m_cubic_sigma_min      = 1e-8;  ///< Minimum sigma
    real_type m_cubic_sigma_max      = 1e8;   ///< Maximum sigma
    real_type m_cubic_gamma_decrease = 0.5;   ///< Decrease factor for sigma
    real_type m_cubic_gamma_increase = 2.0;   ///< Increase factor for sigma
    real_type m_cubic_eta1           = 0.1;   ///< Unsuccessful step threshold
    real_type m_cubic_eta2           = 0.75;  ///< Very successful step threshold

    // =========================================================================
    // QUADRATIC BACKTRACKING PARAMETERS
    // =========================================================================

    real_type m_quad_alpha_init = 1.0;   ///< Initial step size
    real_type m_quad_alpha_min  = 1e-8;  ///< Minimum step size
    real_type m_quad_rho        = 0.5;   ///< Reduction factor
    real_type m_quad_c1         = 1e-4;  ///< Sufficient decrease constant
    real_type m_quad_c2         = 0.9;   ///< Curvature condition constant (compatibility)
    integer   m_quad_max_interp = 10;    ///< Maximum interpolation attempts

    // =========================================================================
    // STATISTICS
    // =========================================================================

    integer   m_num_iterations     = 0;      ///< Number of Newton iterations performed
    integer   m_num_function_evals = 0;      ///< Number of function evaluations
    integer   m_num_jacobian_evals = 0;      ///< Number of Jacobian evaluations
    bool      m_converged          = false;  ///< Convergence flag
    real_type m_final_residual     = 0;      ///< Final residual norm

    bool m_verbose = false;  ///< Verbose output flag

  public:
    // =========================================================================
    // CONSTRUCTOR
    // =========================================================================

    /**
     * @brief Default constructor
     */
    NewtonDumped() = default;

    // =========================================================================
    // GENERAL PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set the convergence tolerance
     * @param tol Tolerance on ||f|| for convergence
     */
    void
    set_tolerance( real_type tol )
    {
      m_tolerance = tol;
    }

    /**
     * @brief Set maximum number of Newton iterations
     * @param max_iter Maximum iterations
     */
    void
    set_max_iterations( integer max_iter )
    {
      m_max_iterations = max_iter;
    }

    /**
     * @brief Set maximum damping iterations per Newton step
     * @param max_damp Maximum damping iterations
     */
    void
    set_max_damping_iterations( integer max_damp )
    {
      m_max_damping_iterations = max_damp;
    }

    /**
     * @brief Enable or disable verbose output
     * @param val True to enable verbose output
     */
    void
    set_verbose( bool val )
    {
      m_verbose = val;
    }

    // =========================================================================
    // BANK & ROSE PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set Bank-Rose alpha parameter
     * @param alpha Damping reduction factor
     */
    void
    set_bank_rose_alpha( real_type alpha )
    {
      m_bank_rose_alpha = alpha;
    }

    /**
     * @brief Set Bank-Rose beta parameter
     * @param beta Sufficient decrease constant
     */
    void
    set_bank_rose_beta( real_type beta )
    {
      m_bank_rose_beta = beta;
    }

    /**
     * @brief Set Bank-Rose gamma parameter
     * @param gamma Contraction factor
     */
    void
    set_bank_rose_gamma( real_type gamma )
    {
      m_bank_rose_gamma = gamma;
    }

    /**
     * @brief Set Bank-Rose minimum theta
     * @param theta Minimum damping factor
     */
    void
    set_bank_rose_theta_min( real_type theta )
    {
      m_bank_rose_theta_min = theta;
    }

    /**
     * @brief Set Bank-Rose maximum theta
     * @param theta Maximum damping factor
     */
    void
    set_bank_rose_theta_max( real_type theta )
    {
      m_bank_rose_theta_max = theta;
    }

    // =========================================================================
    // GRIEWANK PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set Griewank eta parameter
     * @param eta Acceptance threshold
     */
    void
    set_griewank_eta( real_type eta )
    {
      m_griewank_eta = eta;
    }

    /**
     * @brief Set Griewank omega parameter
     * @param omega Reduction factor
     */
    void
    set_griewank_omega( real_type omega )
    {
      m_griewank_omega = omega;
    }

    /**
     * @brief Set Griewank tau parameter
     * @param tau Minimum step size
     */
    void
    set_griewank_tau( real_type tau )
    {
      m_griewank_tau = tau;
    }

    /**
     * @brief Set Griewank zeta parameter
     * @param zeta Contraction factor
     */
    void
    set_griewank_zeta( real_type zeta )
    {
      m_griewank_zeta = zeta;
    }

    // =========================================================================
    // FILTER METHOD PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set filter minimum theta
     * @param theta Minimum constraint violation
     */
    void
    set_filter_theta_min( real_type theta )
    {
      m_filter_theta_min = theta;
    }

    /**
     * @brief Set filter gamma_theta parameter
     * @param gamma Filter parameter for theta
     */
    void
    set_filter_gamma_theta( real_type gamma )
    {
      m_filter_gamma_theta = gamma;
    }

    /**
     * @brief Set filter gamma_f parameter
     * @param gamma Filter parameter for f
     */
    void
    set_filter_gamma_f( real_type gamma )
    {
      m_filter_gamma_f = gamma;
    }

    /**
     * @brief Set filter alpha parameter
     * @param alpha Backtracking factor
     */
    void
    set_filter_alpha( real_type alpha )
    {
      m_filter_alpha = alpha;
    }

    /**
     * @brief Set filter beta parameter
     * @param beta Margin for filter acceptance
     */
    void
    set_filter_beta( real_type beta )
    {
      m_filter_beta = beta;
    }

    // =========================================================================
    // CUBIC TRUST REGION PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set cubic trust region delta
     * @param delta Trust region radius
     */
    void
    set_ctr_delta( real_type delta )
    {
      m_ctr_delta = delta;
    }

    /**
     * @brief Set cubic trust region minimum delta
     * @param delta Minimum trust region radius
     */
    void
    set_ctr_delta_min( real_type delta )
    {
      m_ctr_delta_min = delta;
    }

    /**
     * @brief Set cubic trust region maximum delta
     * @param delta Maximum trust region radius
     */
    void
    set_ctr_delta_max( real_type delta )
    {
      m_ctr_delta_max = delta;
    }

    /**
     * @brief Set cubic trust region eta1
     * @param eta Unsuccessful step threshold
     */
    void
    set_ctr_eta1( real_type eta )
    {
      m_ctr_eta1 = eta;
    }

    /**
     * @brief Set cubic trust region eta2
     * @param eta Very successful step threshold
     */
    void
    set_ctr_eta2( real_type eta )
    {
      m_ctr_eta2 = eta;
    }

    /**
     * @brief Set cubic trust region gamma1
     * @param gamma Shrink factor
     */
    void
    set_ctr_gamma1( real_type gamma )
    {
      m_ctr_gamma1 = gamma;
    }

    /**
     * @brief Set cubic trust region gamma2
     * @param gamma Expand factor
     */
    void
    set_ctr_gamma2( real_type gamma )
    {
      m_ctr_gamma2 = gamma;
    }

    /**
     * @brief Set cubic trust region sigma
     * @param sigma Cubic regularization parameter
     */
    void
    set_ctr_sigma( real_type sigma )
    {
      m_ctr_sigma = sigma;
    }

    // =========================================================================
    // DOGLEG PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set dogleg delta
     * @param delta Trust region radius
     */
    void
    set_dogleg_delta( real_type delta )
    {
      m_dogleg_delta = delta;
    }

    /**
     * @brief Set dogleg minimum delta
     * @param delta Minimum trust region radius
     */
    void
    set_dogleg_delta_min( real_type delta )
    {
      m_dogleg_delta_min = delta;
    }

    /**
     * @brief Set dogleg maximum delta
     * @param delta Maximum trust region radius
     */
    void
    set_dogleg_delta_max( real_type delta )
    {
      m_dogleg_delta_max = delta;
    }

    /**
     * @brief Set dogleg eta1
     * @param eta Unsuccessful step threshold
     */
    void
    set_dogleg_eta1( real_type eta )
    {
      m_dogleg_eta1 = eta;
    }

    /**
     * @brief Set dogleg eta2
     * @param eta Very successful step threshold
     */
    void
    set_dogleg_eta2( real_type eta )
    {
      m_dogleg_eta2 = eta;
    }

    /**
     * @brief Set dogleg gamma1
     * @param gamma Shrink factor
     */
    void
    set_dogleg_gamma1( real_type gamma )
    {
      m_dogleg_gamma1 = gamma;
    }

    /**
     * @brief Set dogleg gamma2
     * @param gamma Expand factor
     */
    void
    set_dogleg_gamma2( real_type gamma )
    {
      m_dogleg_gamma2 = gamma;
    }

    // =========================================================================
    // WOLFE LINE SEARCH PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set Wolfe c1 parameter
     * @param c1 Armijo condition constant
     */
    void
    set_wolfe_c1( real_type c1 )
    {
      m_wolfe_c1 = c1;
    }

    /**
     * @brief Set Wolfe c2 parameter
     * @param c2 Curvature condition constant
     */
    void
    set_wolfe_c2( real_type c2 )
    {
      m_wolfe_c2 = c2;
    }

    /**
     * @brief Set Wolfe initial alpha
     * @param alpha Initial step size
     */
    void
    set_wolfe_alpha_init( real_type alpha )
    {
      m_wolfe_alpha_init = alpha;
    }

    /**
     * @brief Set Wolfe minimum alpha
     * @param alpha Minimum step size
     */
    void
    set_wolfe_alpha_min( real_type alpha )
    {
      m_wolfe_alpha_min = alpha;
    }

    /**
     * @brief Set Wolfe maximum alpha
     * @param alpha Maximum step size
     */
    void
    set_wolfe_alpha_max( real_type alpha )
    {
      m_wolfe_alpha_max = alpha;
    }

    /**
     * @brief Set Wolfe rho parameter
     * @param rho Backtracking factor
     */
    void
    set_wolfe_rho( real_type rho )
    {
      m_wolfe_rho = rho;
    }

    // =========================================================================
    // CUBIC REGULARIZATION (ARC) PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set cubic sigma
     * @param sigma Cubic regularization parameter
     */
    void
    set_cubic_sigma( real_type sigma )
    {
      m_cubic_sigma = sigma;
    }

    /**
     * @brief Set cubic minimum sigma
     * @param sigma Minimum sigma
     */
    void
    set_cubic_sigma_min( real_type sigma )
    {
      m_cubic_sigma_min = sigma;
    }

    /**
     * @brief Set cubic maximum sigma
     * @param sigma Maximum sigma
     */
    void
    set_cubic_sigma_max( real_type sigma )
    {
      m_cubic_sigma_max = sigma;
    }

    /**
     * @brief Set cubic gamma decrease
     * @param gamma Decrease factor for sigma
     */
    void
    set_cubic_gamma_decrease( real_type gamma )
    {
      m_cubic_gamma_decrease = gamma;
    }

    /**
     * @brief Set cubic gamma increase
     * @param gamma Increase factor for sigma
     */
    void
    set_cubic_gamma_increase( real_type gamma )
    {
      m_cubic_gamma_increase = gamma;
    }

    /**
     * @brief Set cubic eta1
     * @param eta Unsuccessful step threshold
     */
    void
    set_cubic_eta1( real_type eta )
    {
      m_cubic_eta1 = eta;
    }

    /**
     * @brief Set cubic eta2
     * @param eta Very successful step threshold
     */
    void
    set_cubic_eta2( real_type eta )
    {
      m_cubic_eta2 = eta;
    }

    // =========================================================================
    // QUADRATIC BACKTRACKING PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set quadratic backtracking initial alpha
     * @param alpha Initial step size
     */
    void
    set_quadratic_backtracking_alpha_init( real_type alpha )
    {
      m_quad_alpha_init = alpha;
    }

    /**
     * @brief Set quadratic backtracking minimum alpha
     * @param alpha Minimum step size
     */
    void
    set_quadratic_backtracking_alpha_min( real_type alpha )
    {
      m_quad_alpha_min = alpha;
    }

    /**
     * @brief Set quadratic backtracking rho
     * @param rho Reduction factor
     */
    void
    set_quadratic_backtracking_rho( real_type rho )
    {
      m_quad_rho = rho;
    }

    /**
     * @brief Set quadratic backtracking c1
     * @param c1 Sufficient decrease constant
     */
    void
    set_quadratic_backtracking_c1( real_type c1 )
    {
      m_quad_c1 = c1;
    }

    /**
     * @brief Set quadratic c1 parameter
     * @param c1 Sufficient decrease constant
     */
    void
    set_quadratic_c1( real_type c1 )
    {
      m_quad_c1 = c1;
    }

    /**
     * @brief Set quadratic c2 parameter
     * @param c2 Curvature condition constant
     */
    void
    set_quadratic_c2( real_type c2 )
    {
      m_quad_c2 = c2;
    }

    /**
     * @brief Set quadratic maximum interpolation attempts
     * @param max_interp Maximum interpolation attempts
     */
    void
    set_quadratic_max_interp( integer max_interp )
    {
      m_quad_max_interp = max_interp;
    }

    // =========================================================================
    // DAMPING STRATEGY SETTERS
    // =========================================================================

    /**
     * @brief Set damping strategy by enum
     * @param strategy The damping strategy to use
     */
    void
    set_damping_strategy( DampingStrategy strategy )
    {
      m_damping_strategy = strategy;
    }

    /**
     * @brief Set strategy to DEUFLHARD
     */
    void
    set_damping_strategy_deuflhard()
    {
      m_damping_strategy = DEUFLHARD;
    }

    /**
     * @brief Set strategy to L2_CLASSIC
     */
    void
    set_damping_strategy_l2_classic()
    {
      m_damping_strategy = L2_CLASSIC;
    }

    /**
     * @brief Set strategy to L2_ADAPTIVE
     */
    void
    set_damping_strategy_l2_adaptive()
    {
      m_damping_strategy = L2_ADAPTIVE;
    }

    /**
     * @brief Set strategy to L2_HYBRID
     */
    void
    set_damping_strategy_l2_hybrid()
    {
      m_damping_strategy = L2_HYBRID;
    }

    /**
     * @brief Set strategy to DOGLEG
     */
    void
    set_damping_strategy_dogleg()
    {
      m_damping_strategy = DOGLEG;
    }

    /**
     * @brief Set strategy to WOLFE_LINE_SEARCH
     */
    void
    set_damping_strategy_wolfe()
    {
      m_damping_strategy = WOLFE_LINE_SEARCH;
    }

    /**
     * @brief Set strategy to CUBIC_REGULARIZATION
     */
    void
    set_damping_strategy_cubic()
    {
      m_damping_strategy = CUBIC_REGULARIZATION;
    }

    /**
     * @brief Set strategy to BACKTRACKING_QUADRATIC
     */
    void
    set_damping_strategy_quadratic_backtracking()
    {
      m_damping_strategy = BACKTRACKING_QUADRATIC;
    }

    /**
     * @brief Set strategy to BANK_ROSE
     */
    void
    set_damping_strategy_bank_rose()
    {
      m_damping_strategy = BANK_ROSE;
    }

    /**
     * @brief Set strategy to GRIEWANK
     */
    void
    set_damping_strategy_griewank()
    {
      m_damping_strategy = GRIEWANK;
    }

    /**
     * @brief Set strategy to FILTER_METHODS
     */
    void
    set_damping_strategy_filter()
    {
      m_damping_strategy = FILTER_METHODS;
    }

    /**
     * @brief Set strategy to CUBIC_TRUST_REGION
     */
    void
    set_damping_strategy_cubic_trust_region()
    {
      m_damping_strategy = CUBIC_TRUST_REGION;
    }

    // =========================================================================
    // DEUFLHARD PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set minimum lambda for Deuflhard damping
     * @param lambda Minimum damping factor
     */
    void
    set_min_lambda( real_type lambda )
    {
      m_min_lambda = lambda;
    }

    // =========================================================================
    // L2 PARAMETER SETTERS
    // =========================================================================

    /**
     * @brief Set initial mu for L2 damping
     * @param mu Initial damping parameter
     */
    void
    set_mu_init( real_type mu )
    {
      m_mu_init = mu;
    }

    /**
     * @brief Set minimum mu for L2 damping
     * @param mu Minimum mu
     */
    void
    set_mu_min( real_type mu )
    {
      m_mu_min = mu;
    }

    /**
     * @brief Set maximum mu for L2 damping
     * @param mu Maximum mu
     */
    void
    set_mu_max( real_type mu )
    {
      m_mu_max = mu;
    }

    /**
     * @brief Set mu increase factor for L2 damping
     * @param factor Increase factor
     */
    void
    set_mu_increase_factor( real_type factor )
    {
      m_mu_increase_factor = factor;
    }

    /**
     * @brief Set mu decrease factor for L2 damping
     * @param factor Decrease factor
     */
    void
    set_mu_decrease_factor( real_type factor )
    {
      m_mu_decrease_factor = factor;
    }

    /**
     * @brief Set trust region radius for L2 adaptive
     * @param radius Initial trust region radius
     */
    void
    set_trust_region_radius( real_type radius )
    {
      m_trust_region_radius = radius;
    }

    /**
     * @brief Set minimum trust region radius
     * @param radius Minimum trust region radius
     */
    void
    set_trust_region_min( real_type radius )
    {
      m_trust_region_min = radius;
    }

    /**
     * @brief Set maximum trust region radius
     * @param radius Maximum trust region radius
     */
    void
    set_trust_region_max( real_type radius )
    {
      m_trust_region_max = radius;
    }

    // =========================================================================
    // STATISTICS GETTERS
    // =========================================================================

    /**
     * @brief Get number of iterations performed
     * @return Number of Newton iterations
     */
    integer
    get_num_iterations() const
    {
      return m_num_iterations;
    }

    /**
     * @brief Get number of function evaluations
     * @return Total function evaluations
     */
    integer
    get_num_function_evals() const
    {
      return m_num_function_evals;
    }

    /**
     * @brief Get number of Jacobian evaluations
     * @return Total Jacobian evaluations
     */
    integer
    get_num_jacobian_evals() const
    {
      return m_num_jacobian_evals;
    }

    /**
     * @brief Check if solver converged
     * @return True if converged within tolerance
     */
    bool
    has_converged() const
    {
      return m_converged;
    }

    /**
     * @brief Get final residual norm
     * @return Final ||f|| value
     */
    real_type
    get_final_residual() const
    {
      return m_final_residual;
    }

    /**
     * @brief Get current damping strategy
     * @return The active damping strategy
     */
    DampingStrategy
    get_damping_strategy() const
    {
      return m_damping_strategy;
    }

  private:
    // -------------------------------------------------------------------------
    // Deuflhard damping strategy (original implementation)
    // -------------------------------------------------------------------------

    /**
     * @brief Deuflhard damping strategy implementation
     *
     * Uses multiplicative damping λ with Armijo-type sufficient decrease
     * condition:
     * ||f(x + λdx)|| ≤ ||f(x)|| - c₁ λ ||f(x)||
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged or acceptable step found
     */
    bool
    solve_deuflhard( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Deuflhard" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        // Solve linear system J dx = -f
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Damping loop (Deuflhard)
        real_type lambda        = 1.0;
        bool      step_accepted = false;
        integer   damping_iter  = 0;

        while ( !step_accepted && damping_iter < m_max_damping_iterations )
        {
          // Candidate update
          x_new.noalias() = x + lambda * dx;

          try
          {
            system.check_if_admissible( x_new );

            // Compute trial residual
            system.evaluate( x_new, f_new );
            ++m_num_function_evals;

            real_type norm_f_new         = f_new.norm();
            real_type actual_reduction   = norm_f - norm_f_new;
            real_type expected_reduction = m_sufficient_decrease * lambda * norm_f;

            // Print damping info
            step_accepted = actual_reduction >= expected_reduction;
            if ( m_verbose ) print_damping_info( lambda, norm_f_new, step_accepted );

            // Accept?
            if ( step_accepted )
            {
              // Accept step
              x      = x_new;
              f      = f_new;
              norm_f = norm_f_new;
            }
            else
            {
              // Reduce lambda
              lambda *= m_lambda_factor;
              ++damping_iter;

              if ( lambda < m_min_lambda )
              {
                if ( m_verbose ) print_damping_failure();
                m_final_residual = norm_f;
                return false;
              }
            }
          }
          catch ( ... )
          {
            // x_new invalid → reduce lambda
            lambda *= m_lambda_factor;
            ++damping_iter;

            if ( m_verbose ) print_invalid_step( lambda );

            if ( lambda < m_min_lambda )
            {
              m_final_residual = norm_f;
              return false;
            }
          }
        }

        if ( !step_accepted )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Classical L2 norm damping (Levenberg-Marquardt style)
    // -------------------------------------------------------------------------

    /**
     * @brief Classical L2 norm damping (Levenberg-Marquardt style)
     *
     * Solves (J^T J + μI) dx = -J^T f with adaptive μ adjustment.
     * Increases μ on failure, decreases on success.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_l2_classic( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type mu = m_mu_init;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "L2-Classic" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        SP_TikhonovSolver2 TS( J, mu );
        dx.noalias() = -TS.solve( f );

        if ( m_verbose ) print_jacobian_ok();

        // Try the step
        x_new.noalias() = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          ++m_num_function_evals;

          real_type norm_f_new = f_new.norm();

          bool accept = norm_f_new < norm_f;

          if ( m_verbose ) print_l2_damping_info( mu, norm_f_new, accept );

          if ( accept )
          {
            // Step accepted - decrease mu
            x.noalias() = x_new;
            f.noalias() = f_new;
            norm_f      = norm_f_new;
            mu *= m_mu_decrease_factor;
            if ( mu < m_mu_min ) mu = m_mu_min;
          }
          else
          {
            // Step rejected - increase mu
            mu *= m_mu_increase_factor;
            if ( mu >= m_mu_max )
            {
              if ( m_verbose ) print_l2_damping_failure( mu );
              m_final_residual = norm_f;
              return false;
            }
          }
        }
        catch ( ... )
        {
          // Invalid step - increase mu
          mu *= m_mu_increase_factor;

          if ( m_verbose ) print_invalid_step_l2( mu );

          if ( mu >= m_mu_max )
          {
            m_final_residual = norm_f;
            return false;
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Adaptive L2 damping with trust region
    // -------------------------------------------------------------------------

    /**
     * @brief Adaptive L2 damping with trust region management
     *
     * Combines L2 regularization with trust region approach for better
     * robustness. Adjusts both μ parameter and trust region radius based on
     * step quality.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_l2_adaptive( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type mu           = m_mu_init;
      real_type trust_radius = m_trust_region_radius;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "L2-Adaptive" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        SP_TikhonovSolver2 TS( J, mu );
        dx.noalias() = -TS.solve( f );

        // Check if step is within trust region
        real_type step_norm = dx.norm();
        if ( step_norm > trust_radius )
        {
          // Scale step to trust region
          dx *= trust_radius / step_norm;
          step_norm = trust_radius;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Try the step
        x_new.noalias() = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new      = f_new.norm();
          real_type reduction_ratio = ( norm_f - norm_f_new ) / ( 0.5 * dx.dot( mu * dx - J.transpose() * f ) );

          if ( m_verbose ) print_trust_region_info( mu, trust_radius, step_norm, norm_f_new, reduction_ratio );

          if ( reduction_ratio > m_acceptance_ratio_bad )
          {
            // Step accepted
            x.noalias() = x_new;
            f           = f_new;
            norm_f      = norm_f_new;

            // Adjust trust region
            if ( reduction_ratio > m_acceptance_ratio_good )
            {
              trust_radius *= 2;
              if ( trust_radius > m_trust_region_max ) trust_radius = m_trust_region_max;
            }

            // Adjust mu
            if ( reduction_ratio > m_acceptance_ratio_good )
            {
              mu *= m_mu_decrease_factor;
              if ( mu < m_mu_min ) mu = m_mu_min;
            }
            else if ( reduction_ratio < m_acceptance_ratio_bad )
            {
              mu *= m_mu_increase_factor;
              if ( mu > m_mu_max ) mu = m_mu_max;
            }
          }
          else
          {
            // Step rejected - reduce trust region
            trust_radius *= 0.5;
            if ( trust_radius < m_trust_region_min ) trust_radius = m_trust_region_min;

            if ( trust_radius <= m_trust_region_min && mu >= m_mu_max )
            {
              if ( m_verbose ) print_trust_region_failure();
              m_final_residual = norm_f;
              return false;
            }
          }
        }
        catch ( ... )
        {
          // Invalid step
          trust_radius *= 0.5;
          if ( trust_radius < m_trust_region_min ) trust_radius = m_trust_region_min;

          mu *= m_mu_increase_factor;
          if ( mu > m_mu_max ) mu = m_mu_max;

          if ( m_verbose ) print_invalid_step_trust( mu, trust_radius );

          if ( trust_radius <= m_trust_region_min && mu >= m_mu_max )
          {
            m_final_residual = norm_f;
            return false;
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Hybrid Deuflhard-L2 strategy
    // -------------------------------------------------------------------------

    /**
     * @brief Hybrid Deuflhard-L2 strategy
     *
     * Dynamically switches between Deuflhard damping and L2 regularization
     * based on performance. Uses L2 when far from solution, switches to
     * Deuflhard when close to convergence.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_l2_hybrid( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type mu     = m_mu_init;
      real_type lambda = 1.0;
      bool      use_l2 = true;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, use_l2 ? "Hybrid-L2" : "Hybrid-Deuflhard" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        if ( use_l2 )
        {
          // L2 step
          SP_TikhonovSolver2 TS( J, mu );
          dx.noalias() = -TS.solve( f );
        }
        else
        {
          // Deuflhard step (plain Newton)
          Eigen::SparseLU<SparseMatrix> solver;
          solver.compute( J );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_jacobian_failure();
            m_final_residual = norm_f;
            return false;
          }

          dx = -solver.solve( f );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_linear_solve_failure();
            m_final_residual = norm_f;
            return false;
          }

          // Scale by lambda
          dx *= lambda;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Try the step
        x_new = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new = f_new.norm();

          if ( use_l2 )
          {
            if ( m_verbose ) print_l2_damping_info( mu, norm_f_new, norm_f_new < norm_f );

            if ( norm_f_new < norm_f )
            {
              // L2 step accepted
              x      = x_new;
              f      = f_new;
              norm_f = norm_f_new;
              mu     = std::max( m_mu_min, mu * m_mu_decrease_factor );

              // Switch to Deuflhard if mu is small (close to Newton)
              if ( mu < 10.0 * m_mu_min )
              {
                use_l2 = false;
                lambda = 1.0;
                if ( m_verbose ) print_strategy_switch( "L2 → Deuflhard" );
              }
            }
            else
            {
              // L2 step rejected
              mu = std::min( m_mu_max, mu * m_mu_increase_factor );

              if ( mu >= m_mu_max )
              {
                // Switch to Deuflhard if L2 fails
                use_l2 = false;
                lambda = 1.0;
                mu     = m_mu_init;
                if ( m_verbose ) print_strategy_switch( "L2 failed → Deuflhard" );
              }
            }
          }
          else
          {
            // Deuflhard step
            if ( m_verbose ) print_damping_info( lambda, norm_f_new, norm_f_new < norm_f );

            if ( norm_f_new < norm_f )
            {
              // Deuflhard step accepted
              x      = x_new;
              f      = f_new;
              norm_f = norm_f_new;
              lambda = 1.0;  // Reset lambda for next iteration
            }
            else
            {
              // Deuflhard step rejected - reduce lambda
              lambda *= m_lambda_factor;

              if ( lambda < m_min_lambda )
              {
                // Switch to L2 if Deuflhard fails
                use_l2 = true;
                lambda = 1.0;
                if ( m_verbose ) print_strategy_switch( "Deuflhard failed → L2" );
              }
            }
          }
        }
        catch ( ... )
        {
          // Invalid step
          if ( use_l2 )
          {
            mu = std::min( m_mu_max, mu * m_mu_increase_factor );
            if ( m_verbose ) print_invalid_step_l2( mu );
          }
          else
          {
            lambda *= m_lambda_factor;
            if ( m_verbose ) print_invalid_step( lambda );
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Bank & Rose (1981) strategy
    // -------------------------------------------------------------------------

    /**
     * @brief Bank and Rose (1981) affine-invariant damping strategy
     *
     * Implements the affine-invariant damping strategy from Bank and Rose
     * (1981). Uses a damping parameter θ that is adjusted based on step
     * quality.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_bank_rose( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type theta = 1.0;  // damping parameter

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Bank-Rose" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        // Solve linear system
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Bank & Rose damping loop
        bool    step_accepted = false;
        integer damping_iter  = 0;

        while ( !step_accepted && damping_iter < m_max_damping_iterations )
        {
          // Compute trial step
          x_new = x + theta * dx;

          try
          {
            system.check_if_admissible( x_new );

            // Evaluate trial point
            system.evaluate( x_new, f_new );
            m_num_function_evals++;

            real_type norm_f_new       = f_new.norm();
            real_type actual_reduction = norm_f - norm_f_new;

            // Bank & Rose condition: sufficient decrease
            real_type expected_reduction = m_bank_rose_beta * theta * norm_f;

            if ( m_verbose )
            {
              fmt::color  col  = ( actual_reduction >= expected_reduction ) ? fmt::color::green : fmt::color::yellow;
              std::string mark = ( actual_reduction >= expected_reduction ) ? "✓" : "…";
              fmt::print( fg( col ), "  θ={:.2e} → {:.2e} Δ={:.2e}/{:.2e} {}\n", theta, norm_f_new, actual_reduction,
                          expected_reduction, mark );
            }

            if ( actual_reduction >= expected_reduction )
            {
              // Accept step
              x             = x_new;
              f             = f_new;
              norm_f        = norm_f_new;
              step_accepted = true;

              // Increase theta for next iteration (more aggressive)
              theta = std::min( m_bank_rose_theta_max, theta / m_bank_rose_gamma );
            }
            else
            {
              // Reduce theta
              theta *= m_bank_rose_alpha;
              damping_iter++;

              if ( theta < m_bank_rose_theta_min )
              {
                if ( m_verbose ) print_damping_failure();
                m_final_residual = norm_f;
                return false;
              }
            }
          }
          catch ( ... )
          {
            // Invalid step
            theta *= m_bank_rose_alpha;
            damping_iter++;

            if ( m_verbose ) print_invalid_step( theta );

            if ( theta < m_bank_rose_theta_min )
            {
              m_final_residual = norm_f;
              return false;
            }
          }
        }

        if ( !step_accepted )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Griewank (1980) strategy
    // -------------------------------------------------------------------------

    /**
     * @brief Griewank (1980) adaptive damping strategy
     *
     * Implements the adaptive damping strategy from Griewank (1980).
     * Uses directional derivative information for step acceptance.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_griewank( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      // Precompute F(x) = 0.5 * ||f||^2
      real_type F = 0.5 * norm_f * norm_f;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Griewank" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian and gradient
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        Vector g = J.transpose() * f;  // Gradient of F

        // Solve for Newton direction
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = -solver.solve( f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Griewank adaptive damping
        real_type lambda        = 1.0;
        bool      step_accepted = false;
        integer   damping_iter  = 0;

        while ( !step_accepted && damping_iter < m_max_damping_iterations )
        {
          x_new.noalias() = x + lambda * dx;

          try
          {
            system.check_if_admissible( x_new );

            system.evaluate( x_new, f_new );
            ++m_num_function_evals;

            real_type norm_f_new = f_new.norm();
            real_type F_new      = 0.5 * norm_f_new * norm_f_new;

            // Griewank condition: check relative improvement
            real_type actual_reduction       = F - F_new;
            real_type directional_derivative = g.dot( dx );
            real_type expected_reduction     = -m_griewank_eta * directional_derivative;

            if ( m_verbose )
            {
              fmt::color  col  = ( actual_reduction >= expected_reduction ) ? fmt::color::green : fmt::color::yellow;
              std::string mark = ( actual_reduction >= expected_reduction ) ? "✓" : "…";
              fmt::print( fg( col ), "  λ={:.2e} → {:.2e} ΔF={:.2e}/{:.2e} {}\n", lambda, norm_f_new, actual_reduction,
                          expected_reduction, mark );
            }

            if ( actual_reduction >= expected_reduction )
            {
              // Accept step
              x             = x_new;
              f             = f_new;
              norm_f        = norm_f_new;
              F             = F_new;
              step_accepted = true;

              // Adjust lambda for next iteration
              if ( actual_reduction > 2.0 * expected_reduction ) { lambda = std::min( 1.0, lambda / m_griewank_zeta ); }
            }
            else
            {
              // Reduce lambda
              lambda *= m_griewank_omega;
              ++damping_iter;

              if ( lambda < m_griewank_tau )
              {
                if ( m_verbose ) print_damping_failure();
                m_final_residual = norm_f;
                return false;
              }
            }
          }
          catch ( ... )
          {
            // Invalid step
            lambda *= m_griewank_omega;
            damping_iter++;

            if ( m_verbose ) print_invalid_step( lambda );

            if ( lambda < m_griewank_tau )
            {
              m_final_residual = norm_f;
              return false;
            }
          }
        }

        if ( !step_accepted )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Filter Methods for nonlinear equations
    // -------------------------------------------------------------------------

    /**
     * @brief Filter methods for nonlinear equations
     *
     * Implements a filter-based acceptance mechanism that maintains
     * a Pareto front of (θ, F) pairs, where θ measures constraint
     * violation and F = 0.5||f||² is the objective.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_filter( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      // Initialize filter
      std::vector<std::pair<real_type, real_type>> filter;
      // Add initial point to filter
      real_type theta = norm_f;                 // constraint violation measure
      real_type F     = 0.5 * norm_f * norm_f;  // objective function
      filter.emplace_back( theta, F );

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Filter" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        // Solve for Newton direction
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Filter line search
        real_type alpha         = 1.0;
        bool      step_accepted = false;
        integer   filter_iter   = 0;

        while ( !step_accepted && filter_iter < m_max_damping_iterations )
        {
          x_new = x + alpha * dx;

          try
          {
            system.check_if_admissible( x_new );

            system.evaluate( x_new, f_new );
            m_num_function_evals++;

            real_type norm_f_new = f_new.norm();
            real_type theta_new  = norm_f_new;
            real_type F_new      = 0.5 * norm_f_new * norm_f_new;

            // Check if (theta_new, F_new) is acceptable to the filter
            bool acceptable_to_filter = true;
            for ( const auto & [theta_f, F_f] : filter )
            {
              if ( theta_new >= theta_f - m_filter_gamma_theta * theta_f && F_new >= F_f - m_filter_gamma_f * F_f )
              {
                acceptable_to_filter = false;
                break;
              }
            }

            // Armijo condition for sufficient decrease
            Vector    g                      = J.transpose() * f;
            real_type directional_derivative = g.dot( dx );
            bool      armijo                 = ( F_new <= F + m_filter_alpha * alpha * directional_derivative );

            if ( m_verbose )
            {
              fmt::color  col  = ( acceptable_to_filter && armijo ) ? fmt::color::green : fmt::color::yellow;
              std::string mark = ( acceptable_to_filter && armijo ) ? "✓" : "…";
              fmt::print( fg( col ), "  α={:.2e} → θ={:.2e}, F={:.2e} Filter:{}, Armijo:{} {}\n", alpha, theta_new,
                          F_new, acceptable_to_filter ? "Y" : "N", armijo ? "Y" : "N", mark );
            }

            if ( acceptable_to_filter && armijo )
            {
              // Accept step
              x             = x_new;
              f             = f_new;
              norm_f        = norm_f_new;
              step_accepted = true;

              // Add to filter if not dominated
              filter.emplace_back( theta_new, F_new );

              // Optionally clean filter (remove dominated entries)
              if ( filter.size() > 10 )
              {  // Keep filter size manageable
                std::sort( filter.begin(), filter.end() );
                filter.erase( std::unique( filter.begin(), filter.end() ), filter.end() );
              }
            }
            else
            {
              // Backtrack
              alpha *= m_filter_beta;
              filter_iter++;

              if ( alpha < m_min_lambda )
              {
                if ( m_verbose ) print_filter_failure();
                m_final_residual = norm_f;
                return false;
              }
            }
          }
          catch ( ... )
          {
            // Invalid step
            alpha *= m_filter_beta;
            filter_iter++;

            if ( m_verbose ) print_invalid_step_filter( alpha );

            if ( alpha < m_min_lambda )
            {
              m_final_residual = norm_f;
              return false;
            }
          }
        }

        if ( !step_accepted )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Cubic Trust Region (CTR) method
    // -------------------------------------------------------------------------

    /**
     * @brief Cubic Trust Region method
     *
     * Implements a trust region method with cubic regularization.
     * Solves a cubic subproblem within a trust region to compute steps.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_cubic_trust_region( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type delta = m_ctr_delta;
      real_type sigma = m_ctr_sigma;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Cubic-TR" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        // Compute gradient g = J^T * f
        Vector g = J.transpose() * f;

        // Solve cubic subproblem within trust region
        // min m(s) = 0.5||f + J s||^2 + (sigma/3)||s||^3 subject to ||s|| <=
        // delta
        Vector s = Vector::Zero( n );  // initial guess

        // Use iterative solver for cubic subproblem
        integer   subproblem_iter      = 0;
        bool      subproblem_converged = false;
        real_type s_norm_old           = 0.0;

        while ( subproblem_iter < 10 && !subproblem_converged )
        {
          real_type s_norm = s.norm();

          // Solve (J^T J + sigma * s_norm * I) ds = -J^T (f + J s) - sigma *
          // s_norm * s
          SparseMatrix A = J.transpose() * J;
          for ( integer i = 0; i < n; ++i ) { A.coeffRef( i, i ) += sigma * s_norm; }

          Vector rhs = -J.transpose() * ( f + J * s ) - sigma * s_norm * s;

          Eigen::SparseLU<SparseMatrix> solver;
          solver.compute( A );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_jacobian_failure();
            m_final_residual = norm_f;
            return false;
          }

          Vector ds = solver.solve( rhs );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_linear_solve_failure();
            m_final_residual = norm_f;
            return false;
          }

          s += ds;

          // Project onto trust region if needed
          real_type new_s_norm = s.norm();
          if ( new_s_norm > delta ) { s = ( delta / new_s_norm ) * s; }

          // Check convergence of subproblem
          real_type rel_change = ds.norm() / ( s_norm + 1.0 );
          subproblem_converged = ( rel_change < 1e-4 || std::abs( new_s_norm - s_norm_old ) < 1e-6 );
          s_norm_old           = new_s_norm;
          subproblem_iter++;
        }

        if ( m_verbose ) print_jacobian_ok();

        // Trial step
        dx    = s;
        x_new = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new       = f_new.norm();
          real_type F                = 0.5 * norm_f * norm_f;
          real_type F_new            = 0.5 * norm_f_new * norm_f_new;
          real_type actual_reduction = F - F_new;

          // Compute model reduction
          real_type s_norm          = dx.norm();
          Vector    Js              = J * dx;
          real_type model_reduction = -g.dot( dx ) - 0.5 * Js.squaredNorm() -
                                      ( sigma / 3.0 ) * s_norm * s_norm * s_norm;

          real_type rho = actual_reduction / model_reduction;

          if ( m_verbose )
          {
            fmt::color  col  = ( rho > m_ctr_eta1 ) ? fmt::color::green
                               : ( rho > 0 )        ? fmt::color::yellow
                                                    : fmt::color::red;
            std::string mark = ( rho > m_ctr_eta2 ) ? "✓✓" : ( rho > m_ctr_eta1 ) ? "✓" : "✗";
            fmt::print( fg( col ), "  Δ={:.2e} σ={:.2e} ρ={:.3f} {} ‖s‖={:.2e}\n", delta, sigma, rho, mark, s_norm );
          }

          // Update trust region and cubic parameter
          if ( rho > m_ctr_eta1 )
          {
            // Successful step
            x      = x_new;
            f      = f_new;
            norm_f = norm_f_new;

            if ( rho > m_ctr_eta2 )
            {
              // Very successful step - expand trust region
              delta = std::min( m_ctr_delta_max, m_ctr_gamma2 * delta );
              // Reduce cubic regularization
              sigma = std::max( 1e-8, sigma * 0.5 );
            }
          }
          else
          {
            // Unsuccessful step - shrink trust region
            delta = std::max( m_ctr_delta_min, m_ctr_gamma1 * delta );
            // Increase cubic regularization
            sigma = std::min( 1e8, sigma * 2.0 );
          }

          if ( delta <= m_ctr_delta_min && sigma >= 1e6 )
          {
            if ( m_verbose ) print_ctr_failure();
            m_final_residual = norm_f;
            return false;
          }
        }
        catch ( ... )
        {
          // Invalid step
          delta = std::max( m_ctr_delta_min, m_ctr_gamma1 * delta );
          sigma = std::min( 1e8, sigma * 2.0 );

          if ( m_verbose ) print_invalid_step_ctr( delta, sigma );

          if ( delta <= m_ctr_delta_min && sigma >= 1e6 )
          {
            m_final_residual = norm_f;
            return false;
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Dogleg method
    // -------------------------------------------------------------------------

    /**
     * @brief Dogleg trust region method
     *
     * Implements the dogleg method that combines steepest descent
     * and Gauss-Newton directions within a trust region.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_dogleg( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), dx_sd( n ), dx_gn( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type delta = m_dogleg_delta;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Dogleg" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        // Compute gradient g = J^T * f
        Vector g = J.transpose() * f;

        // Steepest descent direction
        dx_sd = -g;

        // Gauss-Newton direction
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx_gn = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        // Dogleg path
        real_type alpha     = g.squaredNorm() / ( J * g ).squaredNorm();
        Vector    dx_cauchy = alpha * dx_sd;

        // Choose point on dogleg path
        if ( dx_gn.norm() <= delta )
        {
          dx = dx_gn;  // Full Newton step inside trust region
        }
        else if ( dx_cauchy.norm() >= delta )
        {
          dx = ( delta / dx_cauchy.norm() ) * dx_cauchy;  // Steepest descent
        }
        else
        {
          // Dogleg: linear combination
          Vector    p = dx_gn - dx_cauchy;
          real_type a = p.squaredNorm();
          real_type b = 2.0 * dx_cauchy.dot( p );
          real_type c = dx_cauchy.squaredNorm() - delta * delta;

          real_type tau = ( -b + std::sqrt( b * b - 4.0 * a * c ) ) / ( 2.0 * a );
          dx            = dx_cauchy + tau * p;
        }

        // Trial step
        x_new = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new       = f_new.norm();
          real_type F                = 0.5 * norm_f * norm_f;
          real_type F_new            = 0.5 * norm_f_new * norm_f_new;
          real_type actual_reduction = F - F_new;

          // Predicted reduction from quadratic model
          Vector    Jdx                 = J * dx;
          real_type predicted_reduction = -g.dot( dx ) - 0.5 * Jdx.squaredNorm();

          real_type rho = actual_reduction / predicted_reduction;

          if ( m_verbose )
          {
            fmt::color  col  = ( rho > m_dogleg_eta1 ) ? fmt::color::green
                               : ( rho > 0 )           ? fmt::color::yellow
                                                       : fmt::color::red;
            std::string mark = ( rho > m_dogleg_eta2 ) ? "✓✓" : ( rho > m_dogleg_eta1 ) ? "✓" : "✗";
            fmt::print( fg( col ), "  Δ={:.2e} ρ={:.3f} {} ‖dx‖={:.2e}\n", delta, rho, mark, dx.norm() );
          }

          // Update trust region
          if ( rho > m_dogleg_eta1 )
          {
            // Successful step
            x      = x_new;
            f      = f_new;
            norm_f = norm_f_new;

            if ( rho > m_dogleg_eta2 )
            {
              // Very successful step - expand trust region
              delta = std::min( m_dogleg_delta_max, m_dogleg_gamma2 * delta );
            }
          }
          else
          {
            // Unsuccessful step - shrink trust region
            delta = std::max( m_dogleg_delta_min, m_dogleg_gamma1 * delta );
          }

          if ( delta <= m_dogleg_delta_min )
          {
            if ( m_verbose ) print_dogleg_failure();
            m_final_residual = norm_f;
            return false;
          }
        }
        catch ( ... )
        {
          // Invalid step
          delta = std::max( m_dogleg_delta_min, m_dogleg_gamma1 * delta );

          if ( m_verbose ) print_invalid_step_dogleg( delta );

          if ( delta <= m_dogleg_delta_min )
          {
            m_final_residual = norm_f;
            return false;
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Wolfe line search
    // -------------------------------------------------------------------------

    /**
     * @brief Wolfe line search method
     *
     * Implements line search with strong Wolfe conditions:
     * 1. Armijo condition (sufficient decrease)
     * 2. Curvature condition
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_wolfe( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n ), g( n ), g_new( n );
      SparseMatrix J( n, n );

      // Compute initial objective value
      real_type F = 0.5 * norm_f * norm_f;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Wolfe" );

        // Convergence test
        m_converged = norm_f < m_tolerance;
        if ( m_converged )
        {
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian and gradient
        system.jacobian( x, J );
        ++m_num_jacobian_evals;

        g = J.transpose() * f;

        // Compute Newton direction
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        // Line search with Wolfe conditions
        real_type alpha                  = m_wolfe_alpha_init;
        real_type directional_derivative = g.dot( dx );
        bool      wolfe1 = false, wolfe2 = false;
        integer   ls_iter = 0;

        while ( ls_iter < m_max_damping_iterations )
        {
          x_new = x + alpha * dx;

          try
          {
            system.check_if_admissible( x_new );

            system.evaluate( x_new, f_new );
            m_num_function_evals++;

            real_type norm_f_new = f_new.norm();
            real_type F_new      = 0.5 * norm_f_new * norm_f_new;

            // Armijo condition (Wolfe 1)
            wolfe1 = ( F_new <= F + m_wolfe_c1 * alpha * directional_derivative );

            if ( wolfe1 )
            {
              // Compute gradient at new point for curvature condition
              system.jacobian( x_new, J );
              m_num_jacobian_evals++;
              g_new = J.transpose() * f_new;

              // Curvature condition (Wolfe 2)
              real_type new_directional_derivative = g_new.dot( dx );
              wolfe2 = ( std::abs( new_directional_derivative ) <= m_wolfe_c2 * std::abs( directional_derivative ) );
            }

            if ( m_verbose )
            {
              fmt::color col = ( wolfe1 && wolfe2 ) ? fmt::color::green : wolfe1 ? fmt::color::yellow : fmt::color::red;
              std::string mark = ( wolfe1 && wolfe2 ) ? "✓✓" : wolfe1 ? "✓" : "✗";
              fmt::print( fg( col ), "  α={:.2e} → {:.2e} W1:{}, W2:{} {}\n", alpha, norm_f_new, wolfe1 ? "Y" : "N",
                          wolfe2 ? "Y" : "N", mark );
            }

            if ( wolfe1 && wolfe2 )
            {
              // Wolfe conditions satisfied
              x      = x_new;
              f      = f_new;
              norm_f = norm_f_new;
              F      = F_new;
              break;
            }
            else if ( wolfe1 && !wolfe2 )
            {
              // Curvature condition not satisfied - increase step
              alpha = std::min( m_wolfe_alpha_max, 2.0 * alpha );
            }
            else
            {
              // Armijo condition not satisfied - decrease step
              alpha *= m_wolfe_rho;
            }

            if ( alpha < m_wolfe_alpha_min )
            {
              if ( m_verbose ) print_wolfe_failure();
              m_final_residual = norm_f;
              return false;
            }
          }
          catch ( ... )
          {
            // Invalid step
            alpha *= m_wolfe_rho;

            if ( m_verbose ) print_invalid_step_wolfe( alpha );

            if ( alpha < m_wolfe_alpha_min )
            {
              m_final_residual = norm_f;
              return false;
            }
          }

          ls_iter++;
        }

        if ( ls_iter >= m_max_damping_iterations )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Cubic regularization (ARC)
    // -------------------------------------------------------------------------

    /**
     * @brief Adaptive Regularization by Cubics (ARC) method
     *
     * Implements cubic regularization of the objective function
     * to ensure global convergence while maintaining fast local rate.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_cubic( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      real_type sigma = m_cubic_sigma;

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Cubic-ARC" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian and gradient
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        Vector g = J.transpose() * f;

        // Solve cubic subproblem: min m(s) = f + g^T s + 0.5 s^T J^T J s +
        // (sigma/3) ||s||^3 We use an iterative solver for the cubic subproblem
        dx                      = Vector::Zero( n );
        real_type lambda        = sigma;
        integer   sub_iter      = 0;
        bool      sub_converged = false;

        while ( sub_iter < 10 && !sub_converged )
        {
          // Solve (J^T J + lambda I) dx = -g
          SparseMatrix A = J.transpose() * J;
          for ( integer i = 0; i < n; ++i ) { A.coeffRef( i, i ) += lambda; }

          Vector b = -g;

          Eigen::SparseLU<SparseMatrix> solver;
          solver.compute( A );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_jacobian_failure();
            m_final_residual = norm_f;
            return false;
          }

          Vector dx_new = solver.solve( b );

          if ( solver.info() != Eigen::Success )
          {
            if ( m_verbose ) print_linear_solve_failure();
            m_final_residual = norm_f;
            return false;
          }

          // Update lambda using cubic regularization
          real_type dx_norm = dx_new.norm();
          lambda            = sigma * dx_norm;

          real_type rel_change = ( dx_new - dx ).norm() / ( dx_norm + 1.0 );
          dx                   = dx_new;
          sub_converged        = ( rel_change < 1e-4 );
          sub_iter++;
        }

        // Trial step
        x_new = x + dx;

        try
        {
          system.check_if_admissible( x_new );

          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new       = f_new.norm();
          real_type F                = 0.5 * norm_f * norm_f;
          real_type F_new            = 0.5 * norm_f_new * norm_f_new;
          real_type actual_reduction = F - F_new;

          // Compute model reduction
          real_type dx_norm         = dx.norm();
          Vector    Jdx             = J * dx;
          real_type model_reduction = -g.dot( dx ) - 0.5 * Jdx.squaredNorm() -
                                      ( sigma / 3.0 ) * dx_norm * dx_norm * dx_norm;

          real_type rho = actual_reduction / model_reduction;

          if ( m_verbose )
          {
            fmt::color  col  = ( rho > m_cubic_eta1 ) ? fmt::color::green
                               : ( rho > 0 )          ? fmt::color::yellow
                                                      : fmt::color::red;
            std::string mark = ( rho > m_cubic_eta2 ) ? "✓✓" : ( rho > m_cubic_eta1 ) ? "✓" : "✗";
            fmt::print( fg( col ), "  σ={:.2e} ρ={:.3f} {} ‖dx‖={:.2e}\n", sigma, rho, mark, dx_norm );
          }

          // Update cubic regularization parameter
          if ( rho > m_cubic_eta1 )
          {
            // Successful step
            x      = x_new;
            f      = f_new;
            norm_f = norm_f_new;

            if ( rho > m_cubic_eta2 )
            {
              // Very successful step - decrease sigma
              sigma = std::max( m_cubic_sigma_min, m_cubic_gamma_decrease * sigma );
            }
          }
          else
          {
            // Unsuccessful step - increase sigma
            sigma = std::min( m_cubic_sigma_max, m_cubic_gamma_increase * sigma );
          }

          if ( sigma >= m_cubic_sigma_max )
          {
            if ( m_verbose ) print_cubic_failure();
            m_final_residual = norm_f;
            return false;
          }
        }
        catch ( ... )
        {
          // Invalid step
          sigma = std::min( m_cubic_sigma_max, m_cubic_gamma_increase * sigma );

          if ( m_verbose ) print_invalid_step_cubic( sigma );

          if ( sigma >= m_cubic_sigma_max )
          {
            m_final_residual = norm_f;
            return false;
          }
        }
      }

      return false;
    }

    // -------------------------------------------------------------------------
    // Quadratic backtracking
    // -------------------------------------------------------------------------

    /**
     * @brief Quadratic backtracking line search
     *
     * Implements backtracking line search with quadratic interpolation
     * for step size selection.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Current solution (updated in-place)
     * @param[in,out] f Current residual (updated in-place)
     * @param[in,out] norm_f Current residual norm (updated in-place)
     * @return True if converged
     */
    bool
    solve_quadratic_backtracking( NonlinearSystem & system, Vector & x, Vector & f, real_type & norm_f )
    {
      integer      n = system.num_equations();
      Vector       dx( n ), x_new( n ), f_new( n );
      SparseMatrix J( n, n );

      for ( m_num_iterations = 1; m_num_iterations <= m_max_iterations; ++m_num_iterations )
      {
        if ( m_verbose ) print_iteration_info( m_num_iterations, norm_f, "Quad-Back" );

        // Convergence test
        if ( norm_f < m_tolerance )
        {
          m_converged      = true;
          m_final_residual = norm_f;
          if ( m_verbose ) print_convergence_success();
          return true;
        }

        // Compute Jacobian
        system.jacobian( x, J );
        m_num_jacobian_evals++;

        // Compute Newton direction
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute( J );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_jacobian_failure();
          m_final_residual = norm_f;
          return false;
        }

        dx = solver.solve( -f );

        if ( solver.info() != Eigen::Success )
        {
          if ( m_verbose ) print_linear_solve_failure();
          m_final_residual = norm_f;
          return false;
        }

        // Quadratic backtracking line search
        real_type alpha                  = m_quad_alpha_init;
        real_type F                      = 0.5 * norm_f * norm_f;
        real_type directional_derivative = ( J.transpose() * f ).dot( dx );
        bool      step_accepted          = false;
        integer   bt_iter                = 0;

        while ( !step_accepted && bt_iter < m_max_damping_iterations )
        {
          x_new = x + alpha * dx;

          try
          {
            system.check_if_admissible( x_new );

            system.evaluate( x_new, f_new );
            m_num_function_evals++;

            real_type norm_f_new = f_new.norm();
            real_type F_new      = 0.5 * norm_f_new * norm_f_new;

            // Armijo condition
            real_type expected_reduction = m_quad_c1 * alpha * directional_derivative;
            real_type actual_reduction   = F - F_new;

            if ( m_verbose )
            {
              fmt::color  col  = ( actual_reduction >= expected_reduction ) ? fmt::color::green : fmt::color::yellow;
              std::string mark = ( actual_reduction >= expected_reduction ) ? "✓" : "…";
              fmt::print( fg( col ), "  α={:.2e} → {:.2e} Δ={:.2e}/{:.2e} {}\n", alpha, norm_f_new, actual_reduction,
                          expected_reduction, mark );
            }

            if ( actual_reduction >= expected_reduction )
            {
              // Accept step
              x             = x_new;
              f             = f_new;
              norm_f        = norm_f_new;
              step_accepted = true;
            }
            else
            {
              // Quadratic interpolation for new alpha
              real_type alpha_new = -0.5 * directional_derivative * alpha * alpha /
                                    ( F_new - F - directional_derivative * alpha );

              // Safeguard
              alpha_new = std::max( 0.1 * alpha, std::min( 0.9 * alpha, alpha_new ) );
              alpha     = alpha_new;

              if ( alpha < m_quad_alpha_min )
              {
                if ( m_verbose ) print_quadratic_backtracking_failure();
                m_final_residual = norm_f;
                return false;
              }
            }
          }
          catch ( ... )
          {
            // Invalid step
            alpha *= m_quad_rho;

            if ( m_verbose ) print_invalid_step_quad( alpha );

            if ( alpha < m_quad_alpha_min )
            {
              m_final_residual = norm_f;
              return false;
            }
          }

          bt_iter++;
        }

        if ( !step_accepted )
        {
          if ( m_verbose ) print_no_acceptable_step();
          m_final_residual = norm_f;
          return false;
        }
      }

      return false;
    }


    // -------------------------------------------------------------------------
    // Helper methods for verbose output
    // -------------------------------------------------------------------------

    /**
     * @brief Print iteration information
     * @param iter Current iteration number
     * @param norm_f Current residual norm
     * @param strategy Current strategy name
     */
    void
    print_iteration_info( integer iter, real_type norm_f, const std::string & strategy )
    {
      fmt::print( fmt::fg( fmt::color::light_blue ), "[{:3}][{}] ‖f‖ = {:.2e}", iter, strategy, norm_f );
    }

    /**
     * @brief Print convergence success message
     */
    void
    print_convergence_success()
    {
      fmt::print( fmt::fg( fmt::color::green ), "  ✓ Converged below tolerance ({:.2e})\n", m_tolerance );
    }

    /**
     * @brief Print Jacobian factorization failure
     */
    void
    print_jacobian_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "  ✗ Jacobian factorization failed\n" );
    }

    /**
     * @brief Print linear solve failure
     */
    void
    print_linear_solve_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "  ✗ Linear solve failed\n" );
    }

    /**
     * @brief Print Jacobian OK message
     */
    void
    print_jacobian_ok()
    {
      fmt::print( fmt::fg( fmt::color::green ), "  J " );
    }

    /**
     * @brief Print damping information
     * @param lambda Current damping factor
     * @param norm_f_new New residual norm
     * @param accepted Whether step was accepted
     */
    void
    print_damping_info( real_type lambda, real_type norm_f_new, bool accepted )
    {
      fmt::color  col  = accepted ? fmt::color::green : ( lambda < 1.0 ? fmt::color::yellow : fmt::color::red );
      std::string mark = accepted ? "✓" : ( lambda < 1.0 ? "…" : "✗" );
      fmt::print( fg( col ), "  λ={:.2e} → {:.2e} {}\n", lambda, norm_f_new, mark );
    }

    /**
     * @brief Print L2 damping information
     * @param mu Current L2 damping parameter
     * @param norm_f_new New residual norm
     * @param accepted Whether step was accepted
     */
    void
    print_l2_damping_info( real_type mu, real_type norm_f_new, bool accepted )
    {
      fmt::color  col  = accepted ? fmt::color::green : fmt::color::red;
      std::string mark = accepted ? "✓" : "✗";
      fmt::print( fg( col ), "  μ={:.2e} → {:.2e} {}\n", mu, norm_f_new, mark );
    }

    /**
     * @brief Print trust region information
     * @param mu Current L2 damping parameter
     * @param radius Trust region radius
     * @param step_norm Step norm
     * @param norm_f_new New residual norm
     * @param ratio Reduction ratio
     */
    void
    print_trust_region_info( real_type mu,
                             real_type radius,
                             real_type step_norm,
                             real_type norm_f_new,
                             real_type ratio )
    {
      fmt::color  col  = ( ratio > m_acceptance_ratio_bad ) ? fmt::color::green : fmt::color::yellow;
      std::string mark = ( ratio > m_acceptance_ratio_bad ) ? "✓" : "…";
      fmt::print( fg( col ), "  μ={:.2e} Δ={:.2e} ‖dx‖={:.2e} → {:.2e} ρ={:.2f} {}\n", mu, radius, step_norm,
                  norm_f_new, ratio, mark );
    }

    /**
     * @brief Print damping failure message
     */
    void
    print_damping_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Damping failed (lambda too small)\n" );
    }

    /**
     * @brief Print L2 damping failure message
     * @param mu Current mu value
     */
    void
    print_l2_damping_failure( real_type mu )
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ L2 damping failed (mu={:.2e} too large)\n", mu );
    }

    /**
     * @brief Print trust region failure message
     */
    void
    print_trust_region_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Trust region failed (radius too small, mu too large)\n" );
    }

    /**
     * @brief Print no acceptable step message
     */
    void
    print_no_acceptable_step()
    {
      fmt::print( fmt::fg( fmt::color::red ), "  ✗ No acceptable damping step\n" );
    }

    /**
     * @brief Print invalid step message for lambda
     * @param lambda Current lambda value
     */
    void
    print_invalid_step( real_type lambda )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  λ invalid → reduce to {:.2e}\n", lambda );
    }

    /**
     * @brief Print invalid step message for L2
     * @param mu Current mu value
     */
    void
    print_invalid_step_l2( real_type mu )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → increase μ to {:.2e}\n", mu );
    }

    /**
     * @brief Print invalid step message for trust region
     * @param mu Current mu value
     * @param radius Current trust region radius
     */
    void
    print_invalid_step_trust( real_type mu, real_type radius )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → μ={:.2e}, Δ={:.2e}\n", mu, radius );
    }

    /**
     * @brief Print strategy switch message
     * @param message Switch message
     */
    void
    print_strategy_switch( const std::string & message )
    {
      fmt::print( fmt::fg( fmt::color::cyan ), "  ↳ {} \n", message );
    }

    /**
     * @brief Print filter method failure
     */
    void
    print_filter_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Filter method failed (no acceptable step)\n" );
    }

    /**
     * @brief Print invalid step message for filter
     * @param alpha Current step size
     */
    void
    print_invalid_step_filter( real_type alpha )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → reduce α to {:.2e}\n", alpha );
    }

    /**
     * @brief Print cubic trust region failure
     */
    void
    print_ctr_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ),
                  "    ✗ Cubic trust region failed (trust region too small, "
                  "sigma too large)\n" );
    }

    /**
     * @brief Print invalid step message for cubic trust region
     * @param delta Current trust region radius
     * @param sigma Current sigma value
     */
    void
    print_invalid_step_ctr( real_type delta, real_type sigma )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → Δ={:.2e}, σ={:.2e}\n", delta, sigma );
    }

    /**
     * @brief Print dogleg failure message
     */
    void
    print_dogleg_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Dogleg method failed (trust region too small)\n" );
    }

    /**
     * @brief Print invalid step message for dogleg
     * @param delta Current trust region radius
     */
    void
    print_invalid_step_dogleg( real_type delta )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → reduce Δ to {:.2e}\n", delta );
    }

    /**
     * @brief Print Wolfe line search failure
     */
    void
    print_wolfe_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Wolfe line search failed (step size too small)\n" );
    }

    /**
     * @brief Print invalid step message for Wolfe
     * @param alpha Current step size
     */
    void
    print_invalid_step_wolfe( real_type alpha )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → reduce α to {:.2e}\n", alpha );
    }

    /**
     * @brief Print cubic regularization failure
     */
    void
    print_cubic_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Cubic regularization failed (sigma too large)\n" );
    }

    /**
     * @brief Print invalid step message for cubic
     * @param sigma Current sigma value
     */
    void
    print_invalid_step_cubic( real_type sigma )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → increase σ to {:.2e}\n", sigma );
    }

    /**
     * @brief Print quadratic backtracking failure
     */
    void
    print_quadratic_backtracking_failure()
    {
      fmt::print( fmt::fg( fmt::color::red ), "    ✗ Quadratic backtracking failed (step size too small)\n" );
    }

    /**
     * @brief Print invalid step message for quadratic backtracking
     * @param alpha Current step size
     */
    void
    print_invalid_step_quad( real_type alpha )
    {
      fmt::print( fmt::fg( fmt::color::yellow ), "  step invalid → reduce α to {:.2e}\n", alpha );
    }

  public:
    // -------------------------------------------------------------------------
    // Main solve method - dispatches to appropriate strategy
    // -------------------------------------------------------------------------

    /**
     * @brief Main solve method
     *
     * Solves the nonlinear system F(x) = 0 using the selected damping strategy.
     * Resets statistics, evaluates initial residual, and dispatches to the
     * appropriate solver implementation.
     *
     * @param[in] system The nonlinear system to solve
     * @param[in,out] x Initial guess (input), solution (output)
     * @return True if solver converged successfully
     */
    bool
    solve( NonlinearSystem & system, Vector & x )
    {
      // Reset statistics
      m_num_iterations     = 0;
      m_num_function_evals = 0;
      m_num_jacobian_evals = 0;
      m_converged          = false;
      m_final_residual     = 0.0;

      // Evaluate initial residual
      integer n = system.num_equations();
      Vector  f( n );
      system.evaluate( x, f );
      m_num_function_evals++;
      real_type norm_f = f.norm();

      if ( m_verbose )
      {
        fmt::print( fmt::fg( fmt::color::light_blue ), "[START] ‖f‖ = {:.2e}, Strategy: ", norm_f );
        switch ( m_damping_strategy )
        {
          case DEUFLHARD:
            fmt::print( "Deuflhard\n" );
            break;
          case L2_CLASSIC:
            fmt::print( "L2-Classic\n" );
            break;
          case L2_ADAPTIVE:
            fmt::print( "L2-Adaptive\n" );
            break;
          case L2_HYBRID:
            fmt::print( "L2-Hybrid\n" );
            break;
          case DOGLEG:
            fmt::print( "Dogleg\n" );
            break;
          case WOLFE_LINE_SEARCH:
            fmt::print( "Wolfe Line Search\n" );
            break;
          case CUBIC_REGULARIZATION:
            fmt::print( "Cubic Regularization\n" );
            break;
          case BACKTRACKING_QUADRATIC:
            fmt::print( "Quadratic Backtracking\n" );
            break;
          case BANK_ROSE:
            fmt::print( "Bank & Rose\n" );
            break;
          case GRIEWANK:
            fmt::print( "Griewank\n" );
            break;
          case FILTER_METHODS:
            fmt::print( "Filter Methods\n" );
            break;
          case CUBIC_TRUST_REGION:
            fmt::print( "Cubic Trust Region\n" );
            break;
        }
      }

      // Dispatch to appropriate solver
      bool success = false;
      switch ( m_damping_strategy )
      {
        case DEUFLHARD:
          success = solve_deuflhard( system, x, f, norm_f );
          break;
        case L2_CLASSIC:
          success = solve_l2_classic( system, x, f, norm_f );
          break;
        case L2_ADAPTIVE:
          success = solve_l2_adaptive( system, x, f, norm_f );
          break;
        case L2_HYBRID:
          success = solve_l2_hybrid( system, x, f, norm_f );
          break;
        case DOGLEG:
          success = solve_dogleg( system, x, f, norm_f );
          break;
        case WOLFE_LINE_SEARCH:
          success = solve_wolfe( system, x, f, norm_f );
          break;
        case CUBIC_REGULARIZATION:
          success = solve_cubic( system, x, f, norm_f );
          break;
        case BACKTRACKING_QUADRATIC:
          success = solve_quadratic_backtracking( system, x, f, norm_f );
          break;
        case BANK_ROSE:
          success = solve_bank_rose( system, x, f, norm_f );
          break;
        case GRIEWANK:
          success = solve_griewank( system, x, f, norm_f );
          break;
        case FILTER_METHODS:
          success = solve_filter( system, x, f, norm_f );
          break;
        case CUBIC_TRUST_REGION:
          success = solve_cubic_trust_region( system, x, f, norm_f );
          break;
      }

      // Final output
      if ( !success && !m_converged && m_verbose )
      {
        fmt::print( fmt::fg( fmt::color::yellow ), "  ⚠ Max iterations ({}) reached\n", m_max_iterations );
      }

      return success;
    }
  };

}  // namespace Utils

#endif
