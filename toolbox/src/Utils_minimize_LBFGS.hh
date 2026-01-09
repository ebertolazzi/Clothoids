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
 |  CORRECTED VERSION - Fixed gradient residual issues                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef UTILS_MINIMIZE_LBFGS_HH
#define UTILS_MINIMIZE_LBFGS_HH

#include "Utils_minimize.hh"
#include "Utils_LBFGS.hh"
#include "Utils_Linesearch.hh"

/**
 * @file Utils_minimize_LBFGS.hh
 * @brief Complete header-only implementation of L-BFGS optimization algorithms
 *
 * This file provides a comprehensive implementation of the Limited-memory
 * Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) quasi-Newton optimization method,
 * including multiple line search strategies and support for box constraints.
 *
 * ## Main Components
 *
 * - **LBFGS**: Core L-BFGS storage and two-loop recursion algorithm
 * - **Line Search Policies**: Multiple strategies (Armijo, Wolfe, Strong Wolfe,
 *   Goldstein, Hager-Zhang, More-Thuente)
 * - **LBFGS_minimizer**: High-level minimizer with automatic line search and
 *   projected gradient methods for box-constrained optimization
 *
 * ## Theoretical Background
 *
 * The L-BFGS method approximates the inverse Hessian matrix using a limited
 * number of recent gradient differences. This makes it suitable for large-scale
 * optimization where storing the full Hessian is impractical.
 *
 * ### Algorithm Overview
 *
 * At iteration k, L-BFGS computes a search direction \f$p_k\f$ by:
 * \f[
 *   p_k = -H_k \nabla f(x_k)
 * \f]
 * where \f$H_k\f$ is an approximation to \f$[\nabla^2 f(x_k)]^{-1}\f$
 * constructed from the m most recent correction pairs \f$(s_i, y_i)\f$:
 * \f[
 *   s_i = x_{i+1} - x_i, \quad y_i = \nabla f(x_{i+1}) - \nabla f(x_i)
 * \f]
 *
 * The key advantage is that \f$H_k\f$ is never formed explicitly; instead,
 * the product \f$H_k g\f$ is computed efficiently via two-loop recursion in
 * O(mn) time.
 *
 * ## Corrections in This Version
 *
 * 1. **Proper gradient norm computation**: Uses L2 norm instead of L∞
 * 2. **Stricter convergence criteria**: Primary check on gradient, secondary
 *    on stagnation with tighter tolerances
 * 3. **Improved fallback steps**: Require both function and gradient improvement
 * 4. **Less aggressive memory reset**: Preserves more curvature information
 * 5. **Stricter curvature condition**: Better LBFGS update acceptance
 * 6. **Stagnation tracking**: Counts consecutive weak steps before accepting
 *
 * ## References
 *
 * -# J. Nocedal (1980). "Updating Quasi-Newton Matrices with Limited Storage".
 *    Mathematics of Computation, 35(151), 773-782.
 * -# D.C. Liu and J. Nocedal (1989). "On the Limited Memory BFGS Method for
 *    Large Scale Optimization". Mathematical Programming, 45(1-3), 503-528.
 * -# R.H. Byrd, P. Lu, J. Nocedal, and C. Zhu (1995). "A Limited Memory
 *    Algorithm for Bound Constrained Optimization". SIAM Journal on Scientific
 *    Computing, 16(5), 1190-1208.
 * -# J. Nocedal and S.J. Wright (2006). "Numerical Optimization", 2nd Edition.
 * -# J.J. Moré and D.J. Thuente (1994). "Line Search Algorithms with Guaranteed
 *    Sufficient Decrease". ACM TOMS, 20(3), 286-307.
 * -# W.W. Hager and H. Zhang (2006). "A New Conjugate Gradient Method with
 *    Guaranteed Descent and an Efficient Line Search". SIOPT, 16(1), 170-192.
 *
 * @author Enrico Bertolazzi
 * @date 2025
 */

namespace Utils
{

  // ---------------------------------------------------------------------------
  // LBFGSMinimizer: high-level optimizer with line-search and box constraints
  // ---------------------------------------------------------------------------

  template <typename Scalar = double> class LBFGS_minimizer
  {
  public:
    using integer  = Eigen::Index;
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar( Vector const &, Vector * )>;
    using LBFGS    = Utils::LBFGS<Scalar>;
    using BOX      = BoxConstraintHandler<Scalar>;

    enum class Status
    {
      CONVERGED          = 0,
      MAX_ITERATIONS     = 1,
      LINE_SEARCH_FAILED = 2,
      FAILED             = 3
    };

    static string to_string( Status status )
    {
      switch ( status )
      {
        case Status::CONVERGED: return "CONVERGED";
        case Status::MAX_ITERATIONS: return "MAX_ITER";
        case Status::LINE_SEARCH_FAILED: return "LS_FAILED";
        case Status::FAILED: return "FAILED";
        default: return "UNKNOWN";
      }
    }

    struct Options
    {
      integer max_iter   = 500;
      integer iter_reset = 200;  // FIXED: Less aggressive reset (was 50)
      integer m          = 20;

      Scalar g_tol      = 1e-8;
      Scalar f_tol      = 1e-12;  // FIXED: Tighter (was 1e-10)
      Scalar x_tol      = 1e-20;  // FIXED: Tighter (was 1e-18)
      Scalar g_tol_weak = 1e-6;   // FIXED: Much tighter (was 1e-4)

      Scalar step_max        = 10;
      Scalar very_small_step = 1e-10;  // FIXED: Smaller (was 1e-8)

      integer verbosity_level = 1;  // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
    };

  private:
    Scalar  m_epsi{ std::numeric_limits<Scalar>::epsilon() };
    Options m_options;
    LBFGS   m_LBFGS;
    BOX     m_box_handler;

    // Cache vectors to avoid allocations
    mutable Vector  m_x, m_g, m_p, m_x_new, m_g_new;
    mutable integer m_iter_since_reset{ 0 };
    mutable integer m_function_evaluations{ 0 };
    mutable integer m_consecutive_weak_steps{ 0 };  // NEW: Track stagnation

    // Optimization results
    Status  m_status{ Status::FAILED };
    integer m_total_iterations{ 0 };
    Scalar  m_final_gradient_norm{ 0 };
    Scalar  m_final_function_value{ 0 };
    Scalar  m_initial_function_value{ 0 };

    // ===========================================================================
    // HELPER METHODS
    // ===========================================================================

    /**
     * @brief FIXED: Compute proper gradient norm for convergence checking
     *
     * Uses L2 norm instead of L∞. For box-constrained problems, computes
     * the norm of the projected gradient.
     */
    Scalar compute_convergence_norm() const
    {
      if ( m_box_handler.is_active() )
      {
        // Use L2 norm of projected gradient for bounded problems
        Vector pg = m_box_handler.projected_gradient( m_x, m_g );
        return pg.norm();
      }
      else
      {
        // Use L2 norm for unconstrained problems
        return m_g.norm();
      }
    }

    /**
     * @brief Compute search direction with BoxConstraintHandler integration
     */
    bool compute_search_direction()
    {
      // Compute L-BFGS direction
      Scalar h0 = m_LBFGS.compute_initial_h0( 1 );
      m_p       = -m_LBFGS.two_loop_recursion( m_g, h0 );

      // Handle box constraints
      if ( m_box_handler.is_active() )
      {
        m_box_handler.project_direction( m_x, m_p );

        // Use projected gradient if direction becomes zero
        if ( m_p.isZero( m_epsi * m_p.size() ) ) { m_p = -m_box_handler.projected_gradient( m_x, m_g ); }
      }

      // Check descent condition
      Scalar pg = m_p.dot( m_g );
      if ( !std::isfinite( pg ) ) return false;

      // FIXED: More robust descent check using proper norms
      Scalar gnorm = compute_convergence_norm();
      Scalar pnorm = m_p.norm();

      if ( pnorm < m_epsi || pg > -m_epsi * gnorm * pnorm )
      {
        // Fallback to (projected) negative gradient
        m_p = -m_g;
        if ( m_box_handler.is_active() ) { m_box_handler.project_direction( m_x, m_p ); }
      }

      return true;
    }

    /**
     * @brief FIXED: Improved fallback steps with gradient checking
     *
     * Now requires BOTH function improvement AND gradient not increasing too much.
     * Tries more steps with slower decrement.
     */
    bool try_fallback_steps( Callback const & callback, Scalar f_current, Scalar gnorm_current, Scalar & f_new )
    {
      Scalar base_step = std::min( m_options.very_small_step, Scalar( 1.0 ) / ( m_p.norm() + m_epsi ) );

      Scalar best_f     = f_current;
      Scalar best_gnorm = gnorm_current;
      bool   improved   = false;
      Vector best_x_new = m_x_new;
      Vector best_g_new = m_g_new;

      for ( int i = 0; i < 8; ++i )  // Try more steps (was 5)
      {
        Scalar step       = base_step * std::pow( Scalar( 0.5 ), i );  // Slower decrease (was 0.1)
        m_x_new.noalias() = m_x + step * m_p;
        m_box_handler.project( m_x_new );

        Scalar f_test = callback( m_x_new, &m_g_new );
        m_function_evaluations++;

        // Temporarily update m_x to compute proper gradient norm
        Vector x_backup   = m_x;
        m_x               = m_x_new;
        Scalar gnorm_test = compute_convergence_norm();
        m_x               = x_backup;

        // FIXED: Accept only if BOTH function AND gradient conditions met
        if (
          f_test < best_f - m_epsi * std::abs( f_current ) &&
          gnorm_test < gnorm_current * Scalar( 1.2 ) )  // Allow 20% gradient increase
        {
          best_f     = f_test;
          best_gnorm = gnorm_test;
          best_x_new = m_x_new;
          best_g_new = m_g_new;
          improved   = true;

          // If gradient also decreased significantly, accept immediately
          if ( gnorm_test < gnorm_current * Scalar( 0.9 ) ) break;
        }
      }

      if ( improved )
      {
        f_new   = best_f;
        m_x_new = best_x_new;
        m_g_new = best_g_new;
        return true;
      }

      return false;
    }

    /**
     * @brief FIXED: Stricter convergence logic
     *
     * Primary criterion: gradient norm below g_tol
     * Secondary criterion: multiple consecutive stagnant steps with gradient below g_tol_weak
     */
    Status check_convergence( integer iteration, Scalar f_old, Scalar f_new, Vector const & s, Scalar gnorm )
    {
      // Check iteration limit
      if ( iteration >= m_options.max_iter - 1 ) { return Status::MAX_ITERATIONS; }

      // PRIMARY: Gradient convergence (strict)
      if ( gnorm <= m_options.g_tol ) { return Status::CONVERGED; }

      // SECONDARY: Check for stagnation
      Scalar f_change          = std::abs( f_new - f_old );
      Scalar f_scale           = std::max( std::abs( f_old ), std::abs( f_new ) );
      Scalar f_change_relative = ( f_scale > m_epsi ) ? f_change / f_scale : f_change;

      Scalar x_change          = s.norm();  // FIXED: Use L2 norm (was L∞)
      Scalar x_scale           = m_x.norm();
      Scalar x_change_relative = ( x_scale > m_epsi ) ? x_change / x_scale : x_change;

      // Track consecutive weak steps
      if ( f_change_relative < m_options.f_tol && x_change_relative < m_options.x_tol ) { m_consecutive_weak_steps++; }
      else
      {
        m_consecutive_weak_steps = 0;
      }

      // FIXED: Only accept weak convergence if:
      // 1. Gradient is reasonably small (1e-6, not 1e-4)
      // 2. At least 5 consecutive stagnant steps
      // 3. Both function and variable changes are tiny
      if (
        m_consecutive_weak_steps >= 5 && gnorm <= m_options.g_tol_weak && f_change_relative <= m_options.f_tol &&
        x_change_relative <= m_options.x_tol )
      {
        return Status::CONVERGED;
      }

      return Status::FAILED;  // Not converged yet
    }

    // ===========================================================================
    // PRINTING METHODS
    // ===========================================================================

    void print_optimization_header( integer n, Scalar f0 ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
        "╔══════════════════════════════════════════════════════════════╗\n"
        "║                       L-BFGS Optimization                      ║\n"
        "╠══════════════════════════════════════════════════════════════╣\n"
        "║ {:62} ║\n"
        "║ {:62} ║\n"
        "║ {:62} ║\n"
        "║ {:62} ║\n"
        "║ {:62} ║\n"
        "╚══════════════════════════════════════════════════════════════╝\n",
        fmt::format( "Dimension: {:d}", n ),
        fmt::format( "Max Iterations: {:d}", m_options.max_iter ),
        fmt::format( "Memory (m): {:d}", m_options.m ),
        fmt::format( "Gradient Tolerance: {:.2e}", m_options.g_tol ),
        fmt::format( "Bounds: {}", ( m_box_handler.is_active() ? "Active" : "None" ) ) );
      fmt::print( "Initial F = {:.6e}\n", f0 );
    }

    void print_iteration_summary( integer iter, Scalar f, Scalar gnorm, Scalar step, Scalar pg, bool improved ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      auto   color = improved ? PrintColors::SUCCESS : PrintColors::WARNING;
      string icon  = improved ? "↗" : "→";

      if ( m_options.verbosity_level >= 3 )
      {
        fmt::print(
          color,
          "[{:4d}] {} F = {:<12.6e} | ‖pg‖ = {:<12.6e} | Step = {:<12.6e} | p·g = {:<12.6e}\n",
          iter,
          icon,
          f,
          gnorm,
          step,
          pg );
      }
      else
      {
        fmt::print(
          color,
          "[{:4d}] {} F = {:<12.6e} | ‖pg‖ = {:<12.6e} | Step = {:<12.6e}\n",
          iter,
          icon,
          f,
          gnorm,
          step );
      }
    }

    void print_line_search_result( Scalar step, integer evals, bool success ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      auto   color  = success ? PrintColors::SUCCESS : PrintColors::ERROR;
      string status = success ? "success" : "failure";

      fmt::print( color, "Line search: step = {:<12.6e}, evals = {}, {}\n", step, evals, status );
    }

    void print_lbfgs_update( Scalar sty, bool accepted ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      auto   color  = accepted ? PrintColors::SUCCESS : PrintColors::WARNING;
      string status = accepted ? "accepted" : "rejected";

      fmt::print( color, "L-BFGS update: sᵀy = {:<12.6e}, {}\n", sty, status );
    }

    void print_optimization_statistics() const
    {
      if ( m_options.verbosity_level < 1 ) return;

      string status_str   = to_string( m_status );
      auto   status_color = PrintColors::INFO;

      switch ( m_status )
      {
        case Status::CONVERGED: status_color = PrintColors::SUCCESS; break;
        case Status::MAX_ITERATIONS: status_color = PrintColors::WARNING; break;
        default: status_color = PrintColors::ERROR;
      }

      fmt::print(
        PrintColors::HEADER,
        "╔══════════════════════════════════════════════════════════════╗\n"
        "║                    Optimization Finished                       ║\n"
        "╠══════════════════════════════════════════════════════════════╣\n"
        "║  Final Status       : {:<39}  ║\n"
        "║  Final Value        : {:<39.6e}  ║\n"
        "║  Initial Value      : {:<39.6e}  ║\n"
        "║  Total Iterations   : {:<39}  ║\n"
        "║  Total Evals        : {:<39}  ║\n"
        "║  Final ‖pg‖         : {:<39.6e}  ║\n"
        "║  L-BFGS Memory      : {:<39}  ║\n"
        "╚══════════════════════════════════════════════════════════════╝\n",
        status_str,
        m_final_function_value,
        m_initial_function_value,
        m_total_iterations,
        m_function_evaluations,
        m_final_gradient_norm,
        m_LBFGS.size() );
    }

  public:
    LBFGS_minimizer( Options opts = Options() )
      : m_options( opts )
      , m_LBFGS( opts.m )
      , m_x( 0 )  // Initialize as empty
      , m_g( 0 )
      , m_p( 0 )
      , m_x_new( 0 )
      , m_g_new( 0 )
    {
    }

    // Accessor methods
    Vector const & solution() const { return m_x; }
    Status         status() const { return m_status; }
    integer        total_iterations() const { return m_total_iterations; }
    integer        total_evaluations() const { return m_function_evaluations; }
    Scalar         final_gradient_norm() const { return m_final_gradient_norm; }
    Scalar         final_function_value() const { return m_final_function_value; }
    Scalar         initial_function_value() const { return m_initial_function_value; }

    void set_bounds( integer n, Scalar const lower[], Scalar const upper[] )
    {
      Vector lower_vec = Eigen::Map<const Vector>( lower, n );
      Vector upper_vec = Eigen::Map<const Vector>( upper, n );
      m_box_handler.set_bounds( lower_vec, upper_vec );
    }

    void set_bounds( Vector const & lower, Vector const & upper ) { m_box_handler.set_bounds( lower, upper ); }

    void reset_memory()
    {
      if ( m_options.verbosity_level >= 2 ) fmt::print( "[LBFGS] Periodic memory reset\n" );
      m_LBFGS.clear();
      m_iter_since_reset = 0;
    }

    template <typename Linesearch> void minimize(
      Vector const &     x0,
      Callback const &   callback,
      Linesearch const & linesearch = MoreThuenteLineSearch<Scalar>() )
    {
      // Initialize
      m_function_evaluations   = 0;
      m_consecutive_weak_steps = 0;  // NEW: Reset stagnation counter
      m_status                 = Status::FAILED;
      integer const n          = x0.size();

      // Allocate memory (only if size changed)
      if ( m_x.size() != n )
      {
        m_x.resize( n );
        m_g.resize( n );
        m_x_new.resize( n );
        m_g_new.resize( n );
        m_p.resize( n );
      }

      // Initialize and project
      m_x = x0;
      m_box_handler.project( m_x );

      // Initial evaluation
      Scalar f = callback( m_x, &m_g );
      m_function_evaluations++;
      m_initial_function_value = f;

      print_optimization_header( n, f );

      // Main optimization loop
      Scalar gnorm = compute_convergence_norm();  // FIXED: Use proper norm from start

      for ( integer iter = 0; iter < m_options.max_iter; ++iter )
      {
        ++m_iter_since_reset;

        // Periodic memory reset
        if ( m_iter_since_reset >= m_options.iter_reset ) reset_memory();

        // Compute search direction
        if ( !compute_search_direction() )
        {
          m_status           = Status::FAILED;
          m_total_iterations = iter;
          break;
        }

        Scalar pg    = m_p.dot( m_g );
        Scalar f_old = f;

        // Line search
        auto step_opt = linesearch( f, pg, m_x, m_p, callback, m_options.step_max );

        if ( !step_opt.has_value() )
        {
          if ( m_options.verbosity_level >= 1 )
            fmt::print( PrintColors::WARNING, "[LBFGS] Line search failed, trying fallback\n" );

          m_LBFGS.clear();

          // FIXED: Pass current gradient norm to fallback
          if ( !try_fallback_steps( callback, f, gnorm, f ) )
          {
            m_status           = Status::LINE_SEARCH_FAILED;
            m_total_iterations = iter;
            break;
          }

          // Update with fallback
          Vector s = m_x_new - m_x;
          Vector y = m_g_new - m_g;

          // Check if we can update L-BFGS
          Scalar sty     = s.dot( y );
          bool   updated = false;
          // FIXED: Stricter curvature condition (was m_epsi)
          if ( sty > Scalar( 1e-8 ) * s.norm() * y.norm() ) { updated = m_LBFGS.add_correction( s, y ); }
          else if ( m_options.verbosity_level >= 2 )
          {
            fmt::print( PrintColors::WARNING, "[LBFGS] Skipping update: sty={:.2e} too small\n", sty );
          }

          print_lbfgs_update( sty, updated );

          m_x.swap( m_x_new );
          m_g.swap( m_g_new );

          // FIXED: Recalculate properly after swap
          gnorm = compute_convergence_norm();

          // Check convergence after fallback
          Status conv_status = check_convergence( iter, f_old, f, s, gnorm );
          print_iteration_summary( iter, f, gnorm, Scalar( 0 ), pg, f < f_old );

          if ( conv_status != Status::FAILED )
          {
            m_status           = conv_status;
            m_total_iterations = iter + 1;
            break;
          }

          continue;
        }

        // Successful line search
        auto [step, n_evals] = *step_opt;
        m_function_evaluations += n_evals;

        print_line_search_result( step, n_evals, true );

        // Update point
        m_x_new.noalias() = m_x + step * m_p;
        m_box_handler.project( m_x_new );
        f = callback( m_x_new, &m_g_new );
        m_function_evaluations++;

        // L-BFGS update
        Vector s       = m_x_new - m_x;
        Vector y       = m_g_new - m_g;
        Scalar sty     = s.dot( y );
        bool   updated = false;

        // FIXED: Stricter curvature condition (was m_epsi)
        if ( sty > Scalar( 1e-8 ) * s.norm() * y.norm() ) { updated = m_LBFGS.add_correction( s, y ); }
        else if ( m_options.verbosity_level >= 2 )
        {
          fmt::print( PrintColors::WARNING, "[LBFGS] Skipping update: sty={:.2e} too small\n", sty );
        }

        print_lbfgs_update( sty, updated );

        // Move to new point for gradient norm calculation
        m_x.swap( m_x_new );
        m_g.swap( m_g_new );

        // FIXED: Compute proper gradient norm after swap
        gnorm = compute_convergence_norm();

        // Check convergence
        Status conv_status = check_convergence( iter, f_old, f, s, gnorm );

        print_iteration_summary( iter, f, gnorm, step, pg, f < f_old );

        if ( conv_status != Status::FAILED )
        {
          m_status           = conv_status;
          m_total_iterations = iter + 1;
          break;
        }
      }

      // If we exited because of iteration limit
      if ( m_status == Status::FAILED && m_total_iterations >= m_options.max_iter )
      {
        m_status           = Status::MAX_ITERATIONS;
        m_total_iterations = m_options.max_iter;
      }

      // Store final results
      m_final_gradient_norm  = gnorm;
      m_final_function_value = f;

      print_optimization_statistics();
    }
  };

}  // namespace Utils

#endif

//
// eof: Utils_minimize_LBFGS.hh
//
