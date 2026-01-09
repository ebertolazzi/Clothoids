/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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

/*--------------------------------------------------------------------------*\
 |  Sequential Quadratic Programming (SQP) Optimizer                        |
 |  Improved and extensively documented version                             |
 |                                                                          |
 |  References:                                                             |
 |  [1] Nocedal & Wright, "Numerical Optimization", 2nd ed., Springer 2006  |
 |  [2] Gill et al., "Practical Optimization", Academic Press 1981          |
 |  [3] Powell, "A fast algorithm for nonlinearly constrained               |
 |      optimization calculations", 1978                                    |
 |  [4] Byrd et al., "A limited memory algorithm for bound constrained      |
 |      optimization", SIAM J. Sci. Comput., 1995                           |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/**
 * @file Utils_minimize_SQP.hh
 * @brief Advanced Sequential Quadratic Programming optimizer with comprehensive
 *        box constraint handling and robust convergence guarantees
 */

#pragma once

#ifndef UTILS_MINIMIZE_SQP_HH
#define UTILS_MINIMIZE_SQP_HH

#include "Utils_minimize.hh"
#include "Utils_ssolver.hh"

namespace Utils
{

  /**
   * @class SQP_minimizer
   * @brief Production-grade Sequential Quadratic Programming optimizer
   *
   * Solves box-constrained nonlinear optimization:
   * \f[
   * \begin{aligned}
   *   \min_{x \in \mathbb{R}^n} \quad & f(x) \\
   *   \text{subject to} \quad & \ell \leq x \leq u
   * \end{aligned}
   * \f]
   *
   * @tparam Scalar Floating point type (float, double, long double)
   */
  template <typename Scalar = double> class SQP_minimizer
  {
  public:
    using integer      = Eigen::Index;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using VectorI      = Eigen::Matrix<integer, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using BOX          = BoxConstraintHandler<Scalar>;

    /**
     * @brief User-provided objective function callback
     *
     * @param x Current point
     * @param g [out] Gradient ∇f(x) (nullptr to skip)
     * @param H [out] Hessian ∇²f(x) (nullptr to skip)
     * @return f(x)
     */
    using Callback = std::function<Scalar( Vector const &, Vector *, SparseMatrix * )>;

    /// @brief Optimization termination status
    enum class Status
    {
      CONVERGED          = 0,  ///< KKT conditions satisfied
      MAX_ITERATIONS     = 1,  ///< Maximum iterations reached
      QP_FAILED          = 2,  ///< QP subproblem failed
      GRADIENT_TOO_SMALL = 3,  ///< Gradient norm below threshold
      LINE_SEARCH_FAILED = 4,  ///< Line search failed
      FAILED             = 5   ///< Generic failure
    };

    static std::string to_string( Status s )
    {
      switch ( s )
      {
        case Status::CONVERGED: return "CONVERGED";
        case Status::MAX_ITERATIONS: return "MAX_ITERATIONS";
        case Status::QP_FAILED: return "QP_FAILED";
        case Status::GRADIENT_TOO_SMALL: return "GRADIENT_TOO_SMALL";
        case Status::LINE_SEARCH_FAILED: return "LINE_SEARCH_FAILED";
        case Status::FAILED: return "FAILED";
        default: return "UNKNOWN";
      }
    }

    /**
     * @struct Options
     * @brief Algorithm configuration parameters
     */
    struct Options
    {
      // Convergence tolerances
      integer max_iter  = 500;
      Scalar  g_tol     = 1e-8;
      Scalar  f_tol     = 1e-10;
      Scalar  x_tol     = 1e-18;
      Scalar  kkt_tol   = 1e-6;
      Scalar  max_angle = 89.0;

      // Line search parameters
      Scalar  alpha_init      = 1.0;
      Scalar  alpha_min       = 1e-10;
      Scalar  rho             = 0.5;
      Scalar  c1              = 1e-4;
      Scalar  c2              = 0.9;
      integer max_line_search = 30;

      // Regularization parameters
      Scalar lambda_init            = 1e-6;
      Scalar lambda_min             = 1e-10;
      Scalar lambda_max             = 1e8;
      Scalar lambda_increase_factor = 10.0;
      Scalar lambda_decrease_factor = 0.5;

      // Active set parameters
      Scalar active_tol = 1e-8;

      // Verbosity level
      integer verbosity = 1;
    };

  private:
    Options m_opts;
    BOX     m_box_handler;
    Scalar  m_epsi;

    // Optimization state
    Scalar m_lambda;
    Scalar m_best_f{ std::numeric_limits<Scalar>::max() };
    Vector m_best_x;

    // Results storage
    Status      m_status{ Status::FAILED };
    integer     m_iterations{ 0 };
    integer     m_function_evals{ 0 };
    integer     m_hessian_evals{ 0 };
    integer     m_qp_solves{ 0 };
    integer     m_line_searches{ 0 };
    Scalar      m_final_f{ 0 };
    Scalar      m_initial_f{ 0 };
    Scalar      m_final_grad_norm{ 0 };
    Scalar      m_final_kkt_norm{ 0 };
    Vector      m_solution;
    Vector      m_multipliers;
    integer     m_active_constraints{ 0 };
    std::string m_method_used{ "SQP" };

    // Active set
    VectorI m_active;

    // Pre-allocated working buffers (performance optimization)
    mutable Vector m_work_x;
    mutable Vector m_work_g;

    /**
     * @brief Initialize working buffers
     * @param n Problem dimension
     */
    void initialize_buffers( integer n )
    {
      m_work_x.resize( n );
      m_work_g.resize( n );
    }

    /**
     * @struct QPSolution
     * @brief QP subproblem result
     */
    struct QPSolution
    {
      bool   success{ false };
      Vector p;
      Vector multipliers;
    };

    /**
     * @brief Solve QP subproblem with box constraints
     *
     * Solves: \f$ \min_p \frac{1}{2} p^T (H + \lambda I) p + g^T p \f$
     * subject to: \f$ \ell - x \leq p \leq u - x \f$
     *
     * Delegated to BoxConstraintHandler::solve_qp_subproblem
     */
    bool solve_qp_subproblem( Vector const & x, Vector const & g, SparseMatrix const & H, Vector & p )
    {
      return m_box_handler.solve_qp_subproblem( x, g, H, m_lambda, p, m_multipliers, m_opts.verbosity );
    }

    /**
     * @struct LineSearchResult
     * @brief Line search result
     */
    struct LineSearchResult
    {
      bool    success{ false };
      Scalar  alpha{ 0 };
      Scalar  f_new{ 0 };
      integer evals{ 0 };
    };

    /**
     * @brief Backtracking line search with Armijo-Wolfe conditions
     *
     * Finds α satisfying:
     * - Armijo: f(x + αp) ≤ f(x) + c₁α∇f^Tp
     * - Wolfe: ∇f(x + αp)^Tp ≥ c₂∇f^Tp
     *
     * Utilizes BoxConstraintHandler::project for maintaining feasibility
     */
    LineSearchResult line_search(
      Vector const &   x,
      Vector const &   g,
      Vector const &   p,
      Scalar           f,
      Callback const & callback )
    {
      LineSearchResult result;
      result.evals = 0;

      Scalar directional_deriv = g.dot( p );

      // Check descent direction
      if ( directional_deriv >= 0 )
      {
        if ( m_opts.verbosity >= 3 ) { fmt::print( "  LS: Not descent: {:.2e}\n", directional_deriv ); }
        result.f_new = f;
        return result;
      }

      // Backtracking line search
      integer strict_iters = 0;
      Scalar  alpha        = m_opts.alpha_init;

      for ( integer i = 0; i < m_opts.max_line_search; ++i )
      {
        m_work_x.noalias() = x + alpha * p;
        m_box_handler.project( m_work_x );

        Scalar f_new = callback( m_work_x, nullptr, nullptr );
        result.evals++;

        // Armijo condition
        bool armijo = ( f_new <= f + m_opts.c1 * alpha * directional_deriv );

        if ( armijo )
        {
          // Check Wolfe curvature
          callback( m_work_x, &m_work_g, nullptr );
          result.evals++;

          Scalar new_directional_deriv = m_work_g.dot( p );
          bool   wolfe                 = ( new_directional_deriv >= m_opts.c2 * directional_deriv );

          // Accept if both conditions met or first iteration
          if ( wolfe || i == 0 )
          {
            result.success = true;
            result.alpha   = alpha;
            result.f_new   = f_new;

            if ( m_opts.verbosity >= 3 ) { fmt::print( "  LS: α={:.2e}, Δf={:.2e}\n", alpha, f_new - f ); }
            return result;
          }

          strict_iters++;
        }

        // Relax to Armijo-only after several strict attempts
        if ( strict_iters > 5 && armijo )
        {
          result.success = true;
          result.alpha   = alpha;
          result.f_new   = f_new;

          if ( m_opts.verbosity >= 3 ) { fmt::print( "  LS: Armijo-only α={:.2e}\n", alpha ); }
          return result;
        }

        alpha *= m_opts.rho;
        if ( alpha < m_opts.alpha_min ) break;
      }

      // Last resort: accept any improvement
      alpha = m_opts.alpha_min * 10;
      for ( integer i = 0; i < 3; ++i )
      {
        m_work_x.noalias() = x + alpha * p;
        m_box_handler.project( m_work_x );

        Scalar f_new = callback( m_work_x, nullptr, nullptr );
        result.evals++;

        if ( f_new < f )
        {
          result.success = true;
          result.alpha   = alpha;
          result.f_new   = f_new;

          if ( m_opts.verbosity >= 3 ) { fmt::print( "  LS: Fallback α={:.2e}\n", alpha ); }
          return result;
        }
        alpha *= 0.5;
      }

      if ( m_opts.verbosity >= 3 ) { fmt::print( "  LS: Failed\n" ); }

      result.f_new = f;
      return result;
    }

    // ========== Printing Methods ==========

    void print_header( integer n, Scalar f0 ) const
    {
      if ( m_opts.verbosity < 1 ) return;
      fmt::print(
        "╔═══════════════════════════════════════════════════════╗\n"
        "║      Sequential Quadratic Programming (SQP)           ║\n"
        "╠═══════════════════════════════════════════════════════╣\n"
        "║ Dimension: {:<4d}   Max Iter: {:<6d}                    ║\n"
        "║ Initial F = {:<8.2g} GTol: {:<8.2g} KKT Tol: {:<8.2g} ║\n"
        "╚═══════════════════════════════════════════════════════╝\n",
        n,
        m_opts.max_iter,
        f0,
        m_opts.g_tol,
        m_opts.kkt_tol );
    }

    void print_iteration( Scalar f, Scalar gnorm, Scalar kkt_norm, Scalar alpha, bool ok ) const
    {
      if ( m_opts.verbosity < 2 ) return;

      auto color = ok ? PrintColors::SUCCESS : PrintColors::WARNING;
      fmt::print(
        color,
        "[{:4d}] {} f={:12.6e} ‖∇f‖={:8.2e} KKT={:8.2e} α={:6.4f} Act={:3d}\n",
        m_iterations,
        ok ? "✓" : "✗",
        f,
        gnorm,
        kkt_norm,
        alpha,
        m_active_constraints );
    }

    void print_convergence( Scalar f, Scalar gnorm, Scalar kkt ) const
    {
      if ( m_opts.verbosity < 1 ) return;

      auto converged = m_status == Status::CONVERGED || m_status == Status::GRADIENT_TOO_SMALL;
      auto color     = converged ? PrintColors::SUCCESS : PrintColors::WARNING;

      fmt::print(
        color,
        "╔═══════════════════════════════════════════════════════╗\n"
        "║ Status: {:20s} Iterations: {:<4d}         ║\n"
        "║ f = {:<10.5g}   ‖∇f‖ = {:<10.5g}   KKT = {:<10.5g} ║\n"
        "╚═══════════════════════════════════════════════════════╝\n",
        to_string( m_status ),
        m_iterations,
        f,
        gnorm,
        kkt );
    }

  public:
    /**
     * @brief Constructor
     * @param opts Algorithm options
     */
    SQP_minimizer( Options opts = Options() )
      : m_opts( opts ), m_epsi( std::numeric_limits<Scalar>::epsilon() ), m_lambda( opts.lambda_init )
    {
      // Inizializza i vettori di lavoro con dimensione minima
      m_work_x.resize( 0 );
      m_work_g.resize( 0 );
    }

    /**
     * @brief Set box constraints
     * @param lower Lower bounds ℓ
     * @param upper Upper bounds u
     */
    void set_bounds( Vector const & lower, Vector const & upper ) { m_box_handler.set_bounds( lower, upper ); }

    /**
     * @brief Minimize objective function
     *
     * Main SQP algorithm with all improvements applied. Uses BoxConstraintHandler for:
     * - Projection operations to maintain feasibility
     * - Active set identification
     * - QP subproblem solving with box constraints
     * - KKT condition checking
     *
     * @param x0 Initial point
     * @param callback Objective function
     */
    void minimize( Vector const & x0, Callback const & callback )
    {
      // Initialize
      m_lambda = m_opts.lambda_init;
      m_best_f = std::numeric_limits<Scalar>::max();
      m_active.resize( x0.size() );
      m_active.setZero();

      integer n = x0.size();
      initialize_buffers( n );

      Vector x = x0, g( n ), g_new( n );

      m_multipliers.resize( 2 * n );
      m_multipliers.setZero();
      SparseMatrix H( n, n );

      // Project to feasible region using BoxConstraintHandler
      m_box_handler.project( x );

      // Initial evaluation
      Scalar  f       = callback( x, &g, &H );
      integer f_evals = 1;
      m_hessian_evals = 1;
      m_initial_f     = f;

      m_best_f = f;
      m_best_x = x;

      print_header( n, m_initial_f );

      // Main loop
      m_status               = Status::MAX_ITERATIONS;
      m_qp_solves            = 0;
      m_line_searches        = 0;
      integer total_ls_evals = 0;

      for ( m_iterations = 0; m_iterations < m_opts.max_iter; ++m_iterations )
      {
        // Step 1: Identify active constraints using BoxConstraintHandler
        m_active_constraints = m_box_handler.identify_active_set( x, g, m_active );

        // Step 2: Update multipliers and compute convergence metrics
        m_box_handler.update_multipliers( x, g );
        Scalar proj_gnorm = m_box_handler.projected_gradient_norm_inf( x, g );

        // Compute KKT norm using BoxConstraintHandler's total KKT error
        Scalar kkt_norm = m_box_handler.total_kkt_error( x, g );

        // Update stored multipliers
        integer n_vars               = x.size();
        m_multipliers.head( n_vars ) = m_box_handler.lambda_lower();
        m_multipliers.tail( n_vars ) = m_box_handler.lambda_upper();

        // Step 3: Check convergence
        if ( proj_gnorm <= m_opts.g_tol )
        {
          if ( kkt_norm <= m_opts.kkt_tol ) { m_status = Status::CONVERGED; }
          else
          {
            m_status = Status::GRADIENT_TOO_SMALL;
          }
          print_convergence( f, proj_gnorm, kkt_norm );
          break;
        }

        callback( x, nullptr, &H );
        ++m_hessian_evals;

        // Step 5: Solve QP subproblem using BoxConstraintHandler
        Vector qp_p;
        bool   qp_sol_success = solve_qp_subproblem( x, g, H, qp_p );
        ++m_qp_solves;

        if ( !qp_sol_success )
        {
          m_lambda = std::min( m_opts.lambda_max, m_lambda * m_opts.lambda_increase_factor );

          if ( m_lambda >= m_opts.lambda_max )
          {
            m_status = Status::QP_FAILED;
            print_convergence( f, proj_gnorm, kkt_norm );
            break;
          }
          continue;
        }

        // Step 6: Check descent direction
        Scalar directional_deriv = g.dot( qp_p );

        Scalar cos_angle = -directional_deriv / ( g.norm() * qp_p.norm() + Scalar( 1e-100 ) );
        Scalar angle     = Scalar( 180 ) * std::acos( std::clamp( cos_angle, Scalar( -1 ), Scalar( 1 ) ) ) / m_pi;

        if ( m_opts.verbosity >= 3 ) { fmt::print( "  angle = {:.2f}°, λ = {:.2e}\n", angle, m_lambda ); }

        bool   use_projected_gradient = false;
        Vector pg;

        // If the QP direction is not a descent direction or the angle is too large, try projected gradient
        if ( angle > m_opts.max_angle || directional_deriv >= 0 )
        {
          pg             = m_box_handler.projected_gradient( x, g );
          Scalar pg_norm = pg.template lpNorm<Eigen::Infinity>();

          if ( pg_norm > m_opts.g_tol )
          {
            qp_p                   = -pg;
            directional_deriv      = g.dot( qp_p );
            use_projected_gradient = true;

            // Aggiorna i moltiplicatori quando usi il gradiente proiettato
            m_box_handler.update_multipliers( x, g );
            integer n               = x.size();
            m_multipliers.head( n ) = m_box_handler.lambda_lower();
            m_multipliers.tail( n ) = m_box_handler.lambda_upper();

            if ( m_opts.verbosity >= 3 ) { fmt::print( "  Using projected gradient, ‖pg‖_∞ = {:.2e}\n", pg_norm ); }
          }
          else
          {
            // Projected gradient is too small, consider convergence
            // Assicurati che i moltiplicatori siano aggiornati
            m_box_handler.update_multipliers( x, g );
            integer n               = x.size();
            m_multipliers.head( n ) = m_box_handler.lambda_lower();
            m_multipliers.tail( n ) = m_box_handler.lambda_upper();

            // Ricalcola KKT con i moltiplicatori aggiornati
            kkt_norm = m_box_handler.total_kkt_error( x, g );

            if ( kkt_norm <= m_opts.kkt_tol ) { m_status = Status::CONVERGED; }
            else
            {
              m_status = Status::GRADIENT_TOO_SMALL;
            }
            print_convergence( f, proj_gnorm, kkt_norm );
            break;
          }
        }

        // If we used projected gradient, skip the regularization increase and angle check
        if ( !use_projected_gradient )
        {
          if ( angle > m_opts.max_angle )
          {
            m_lambda = std::min( m_opts.lambda_max, m_lambda * m_opts.lambda_increase_factor );

            if ( m_lambda >= m_opts.lambda_max )
            {
              m_status = Status::QP_FAILED;
              print_convergence( f, proj_gnorm, kkt_norm );
              break;
            }
            continue;
          }
        }

        // Step 7: Check step size
        Scalar step_norm = qp_p.template lpNorm<Eigen::Infinity>();

        if ( step_norm < m_opts.x_tol )
        {
          // If step is too small, try projected gradient if not already used
          if ( !use_projected_gradient )
          {
            pg             = m_box_handler.projected_gradient( x, g );
            Scalar pg_norm = pg.template lpNorm<Eigen::Infinity>();

            if ( pg_norm > m_opts.g_tol )
            {
              qp_p                   = -pg;
              step_norm              = qp_p.template lpNorm<Eigen::Infinity>();
              directional_deriv      = g.dot( qp_p );
              use_projected_gradient = true;

              if ( m_opts.verbosity >= 3 )
              {
                fmt::print( "  Step too small, using projected gradient, ‖pg‖_∞ = {:.2e}\n", pg_norm );
              }
            }
          }

          // Check again after possibly switching to projected gradient
          if ( step_norm < m_opts.x_tol )
          {
            if ( proj_gnorm <= m_opts.g_tol )
            {
              m_status = Status::CONVERGED;
              print_convergence( f, proj_gnorm, kkt_norm );
              break;
            }

            // If step is still too small and gradient is not small, increase regularization
            if ( !use_projected_gradient )
            {
              m_lambda = std::min( m_opts.lambda_max, m_lambda * m_opts.lambda_increase_factor );

              if ( m_lambda >= m_opts.lambda_max )
              {
                m_status = Status::FAILED;
                print_convergence( f, proj_gnorm, kkt_norm );
                break;
              }
              continue;
            }
            // If using projected gradient and step is too small, but gradient is not small, we still proceed to line
            // search The line search will handle very small steps by reducing alpha.
          }
        }
        // Step 8: Line search (uses BoxConstraintHandler::project internally)
        LineSearchResult ls_result = line_search( x, g, qp_p, f, callback );
        ++m_line_searches;
        total_ls_evals += ls_result.evals;

        if ( !ls_result.success || ls_result.alpha < m_opts.alpha_min )
        {
          m_status = Status::LINE_SEARCH_FAILED;
          print_convergence( f, proj_gnorm, kkt_norm );
          break;
        }

        // Step 9: Update solution
        Vector x_new = x + ls_result.alpha * qp_p;

        m_box_handler.project( x_new );

        Scalar f_new = callback( x_new, &g_new, nullptr );
        ++f_evals;

        // Step 10: Accept or reject step
        bool accept_step = ( f_new < f );  // Strict decrease requirement
        print_iteration( f_new, proj_gnorm, kkt_norm, ls_result.alpha, accept_step );

        if ( accept_step )
        {
          // Update state (no BFGS updates anymore)
          x = x_new;
          g = g_new;
          f = f_new;

          // Track best
          if ( f < m_best_f )
          {
            m_best_f = f;
            m_best_x = x;
          }

          // Decrease regularization on success
          m_lambda = std::max( m_opts.lambda_min, m_lambda * m_opts.lambda_decrease_factor );
        }
        else
        {
          // Increase regularization on failure
          m_lambda = std::min( m_opts.lambda_max, m_lambda * m_opts.lambda_increase_factor );
        }

        // Step 11: Check function convergence
        if ( accept_step )
        {
          Scalar f_change = std::abs( f - f_new );
          if ( f_change <= m_opts.f_tol && proj_gnorm <= m_opts.g_tol )
          {
            m_status = Status::CONVERGED;
            ++m_iterations;
            print_convergence( f, proj_gnorm, kkt_norm );
            break;
          }
        }
      }

      // Final: Restore best solution if needed
      if ( m_best_f < f )
      {
        x = m_best_x;
        f = m_best_f;
        callback( x, &g, nullptr );
        f_evals++;
      }

      // Assicurati che i moltiplicatori siano aggiornati per la soluzione finale
      m_box_handler.update_multipliers( x, g );
      integer n_vars               = x.size();
      m_multipliers.head( n_vars ) = m_box_handler.lambda_lower();
      m_multipliers.tail( n_vars ) = m_box_handler.lambda_upper();

      // Store results
      m_function_evals  = f_evals + total_ls_evals;
      m_final_f         = f;
      m_final_grad_norm = m_box_handler.projected_gradient_norm_inf( x, g );
      m_final_kkt_norm  = m_box_handler.total_kkt_error( x, g );
      m_solution        = x;
    }

    /// @brief Get optimization status
    Status status() const { return m_status; }

    /// @brief Get number of iterations performed
    integer iterations() const { return m_iterations; }

    /// @brief Get number of function evaluations
    integer function_evals() const { return m_function_evals; }

    /// @brief Get number of Hessian evaluations
    integer hessian_evals() const { return m_hessian_evals; }

    /// @brief Get number of QP subproblems solved
    integer qp_solves() const { return m_qp_solves; }

    /// @brief Get number of line searches performed
    integer line_searches() const { return m_line_searches; }

    /// @brief Get final objective function value
    Scalar final_f() const { return m_final_f; }

    /// @brief Get initial objective function value
    Scalar initial_f() const { return m_initial_f; }

    /// @brief Get final gradient norm
    Scalar final_grad_norm() const { return m_final_grad_norm; }

    /// @brief Get final KKT norm
    Scalar final_kkt_norm() const { return m_final_kkt_norm; }

    /// @brief Get final solution
    Vector const & solution() const { return m_solution; }

    /// @brief Get Lagrange multipliers
    Vector const & multipliers() const { return m_multipliers; }

    /// @brief Get number of active constraints
    integer active_constraints() const { return m_active_constraints; }

    /// @brief Get method used
    std::string const & method_used() const { return m_method_used; }

    /// @brief Get best solution found
    Vector const & best_solution() const { return m_best_x; }

    /// @brief Get current options
    Options const & options() const { return m_opts; }

    /// @brief Set new options
    void set_options( Options const & opts ) { m_opts = opts; }
  };

}  // namespace Utils

#endif  // UTILS_SQP_IMPROVED_HH

//
// eof: Utils_minimize_SQP.hh
//
