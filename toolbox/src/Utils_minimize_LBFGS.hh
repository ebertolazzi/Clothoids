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

//
// file: Utils_LBFGS.hh
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_LBFGS_MINIMIZE_dot_HH
#define UTILS_LBFGS_MINIMIZE_dot_HH

#include <set>

#include "Utils_LBFGS.hh"
#include "Utils_fmt.hh"
#include "Utils_nonlinear_linesearch.hh"

/**
 * @file Utils_LBFGS.hh
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
 * - **LBFGS_BlockCoordinate**: Enhanced block coordinate L-BFGS with
 * overlapping random consecutive blocks
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
 * ## References
 *
 * -# J. Nocedal (1980). "Updating Quasi-Newton Matrices with Limited Storage".
 *    Mathematics of Computation, 35(151), 773-782.
 *    DOI: 10.1090/S0025-5718-1980-0572855-7
 *
 * -# D.C. Liu and J. Nocedal (1989). "On the Limited Memory BFGS Method for
 *    Large Scale Optimization". Mathematical Programming, 45(1-3), 503-528.
 *    DOI: 10.1007/BF01589116
 *
 * -# R.H. Byrd, P. Lu, J. Nocedal, and C. Zhu (1995). "A Limited Memory
 * Algorithm for Bound Constrained Optimization". SIAM Journal on Scientific
 * Computing, 16(5), 1190-1208. DOI: 10.1137/0916069
 *
 * -# J. Nocedal and S.J. Wright (2006). "Numerical Optimization", 2nd Edition,
 *    Springer. ISBN: 978-0-387-30303-1
 *
 * -# J.J. Moré and D.J. Thuente (1994). "Line Search Algorithms with Guaranteed
 *    Sufficient Decrease". ACM Transactions on Mathematical Software, 20(3),
 * 286-307. DOI: 10.1145/192115.192132
 *
 * -# W.W. Hager and H. Zhang (2006). "A New Conjugate Gradient Method with
 *    Guaranteed Descent and an Efficient Line Search". SIAM Journal on
 * Optimization, 16(1), 170-192. DOI: 10.1137/030601880
 *
 * ## Implementation Features
 *
 * - **Robustness**: Includes Powell damping, curvature checks, and fallback
 * strategies
 * - **Flexibility**: Template-based design supports different scalar types
 * (float/double)
 * - **Efficiency**: Uses Eigen for vectorized operations, minimal memory
 * allocations
 * - **Box Constraints**: Implements projected gradient methods for
 * bound-constrained problems
 * - **Multiple Line Searches**: Choose the best strategy for your problem
 * characteristics
 * - **Block Coordinate**: Enhanced block coordinate descent with overlapping
 * blocks
 *
 * ## Usage Example
 *
 * @code{.cpp}
 * using namespace Utils;
 *
 * // Define objective function
 * auto rosenbrock = [](Vector const& x, Vector* g) -> double {
 *   double f = 100*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0)) + (1-x(0))*(1-x(0));
 *   if (g) {
 *     (*g)(0) = -400*(x(1)-x(0)*x(0))*x(0) - 2*(1-x(0));
 *     (*g)(1) = 200*(x(1)-x(0)*x(0));
 *   }
 *   return f;
 * };
 *
 * // Setup minimizer
 * LBFGS_minimizer<double>::Options opts;
 * opts.max_iter = 1000;
 * opts.g_tol = 1e-6;
 * opts.verbosity_level = 2;  // Verbosity like NelderMead
 *
 * LBFGS_minimizer<double> minimizer(opts);
 *
 * // Initial point
 * Vector x0(2);
 * x0 << -1.2, 1.0;
 *
 * // Minimize
 * auto [status, x_opt, f_opt, data] = minimizer.minimize(
 *   x0, rosenbrock, StrongWolfeLineSearch<double>()
 * );
 *
 * std::cout << "Solution: " << x_opt.transpose() << std::endl;
 * std::cout << "Minimum: " << f_opt << std::endl;
 * @endcode
 *
 * ## Block Coordinate L-BFGS Example
 *
 * @code{.cpp}
 * using namespace Utils;
 *
 * // Define objective function
 * auto rosenbrock_high_dim = [](Vector const& x, Vector* g) -> double {
 *   double f = 0.0;
 *   for (int i = 0; i < x.size()-1; ++i) {
 *     f += 100*(x(i+1)-x(i)*x(i))*(x(i+1)-x(i)*x(i)) + (1-x(i))*(1-x(i));
 *   }
 *   if (g) {
 *     g->setZero(x.size());
 *     for (int i = 0; i < x.size()-1; ++i) {
 *       (*g)(i) += -400*(x(i+1)-x(i)*x(i))*x(i) - 2*(1-x(i));
 *       (*g)(i+1) += 200*(x(i+1)-x(i)*x(i));
 *     }
 *   }
 *   return f;
 * };
 *
 * // Setup block coordinate minimizer
 * LBFGS_BlockCoordinate<double>::Options opts;
 * opts.block_size = 20;
 * opts.max_outer_iterations = 50;
 * opts.block_selection = "sliding_window";
 * opts.overlap_ratio = 0.3;
 * opts.verbosity_level = 2;
 *
 * LBFGS_BlockCoordinate<double> block_minimizer(opts);
 *
 * // Initial point
 * Vector x0(30);
 * x0.setConstant(-1.2);
 *
 * // Minimize
 * auto result = block_minimizer.minimize(x0, rosenbrock_high_dim);
 *
 * std::cout << "Solution norm: " << result.solution.norm() << std::endl;
 * std::cout << "Minimum: " << result.final_function_value << std::endl;
 * @endcode
 *
 * @author Enrico Bertolazzi
 * @date 2025
 */

namespace Utils
{

  // ---------------------------------------------------------------------------
  // LBFGSMinimizer: high-level optimizer with line-search and optional box
  // bounds
  // ---------------------------------------------------------------------------

  template <typename Scalar = double>
  class LBFGS_minimizer
  {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar( Vector const &, Vector * )>;

    enum class Status
    {
      CONVERGED          = 0,  // Convergenza raggiunta
      MAX_ITERATIONS     = 1,  // Massimo numero di iterazioni raggiunto
      LINE_SEARCH_FAILED = 2,  // Line search fallita
      GRADIENT_TOO_SMALL = 3,  // Gradiente troppo piccolo
      FAILED             = 4   // Fallimento generico
    };

    static string
    to_string( Status status )
    {
      switch ( status )
      {
        case Status::CONVERGED:
          return "CONVERGED";
        case Status::MAX_ITERATIONS:
          return "MAX_ITER";
        case Status::LINE_SEARCH_FAILED:
          return "LINE_SEARCH_FAILED";
        case Status::GRADIENT_TOO_SMALL:
          return "GRAD_SMALL";
        case Status::FAILED:
          return "FAILED";
        default:
          return "UNKNOWN";
      }
    }

    struct Result
    {
      Status status{ Status::FAILED };
      size_t total_iterations{ 0 };
      size_t total_evaluations{ 0 };
      Scalar final_gradient_norm{ 0 };
      Scalar final_function_value{ 0 };
      Scalar initial_function_value{ 0 };
      size_t line_search_evaluations{ 0 };  // Aggiunto per tracciare eval line search
    };

    struct Options
    {
      size_t max_iter{ 200 };
      size_t iter_reset{ 50 };
      size_t m{ 20 };

      Scalar g_tol{ 1e-8 };
      Scalar g_tol_weak{ 1e-4 };
      Scalar f_tol{ 1e-12 };
      Scalar x_tol{ 1e-10 };

      Scalar step_max{ 10 };
      Scalar sty_min_factor{ 1e-12 };
      Scalar very_small_step{ 1e-8 };

      bool use_projection{ false };

      // Verbosity system like NelderMead
      size_t verbosity_level{ 1 };         // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      bool   use_unicode_borders{ true };  // Usare bordi Unicode come NelderMead
    };

  private:
    Scalar               m_epsi{ std::numeric_limits<Scalar>::epsilon() };
    Options              m_options;
    Utils::LBFGS<Scalar> m_LBFGS;
    Vector               m_lower;
    Vector               m_upper;
    Vector               m_tol_lower;
    Vector               m_tol_upper;

    string m_indent{ "" };  // Per indentazione consistente

    // Helper functions for bound checks with tolerances
    bool
    is_on_lower_bound( Scalar x_i, Scalar lower_i ) const
    {
      return ( x_i <= lower_i + m_epsi * ( 1 + std::abs( lower_i ) ) );
    }

    bool
    is_on_upper_bound( Scalar x_i, Scalar upper_i ) const
    {
      return ( x_i >= upper_i - m_epsi * ( 1 + std::abs( upper_i ) ) );
    }

    void
    check_bounds_consistency() const
    {
      if ( m_options.use_projection && ( m_lower.array() > m_upper.array() ).any() )
      {
        throw std::invalid_argument( "Lower bounds must be <= upper bounds" );
      }
    }

    // Project x into bounds
    void
    project_inplace( Vector & x ) const
    {
      x = x.cwiseMax( m_lower ).cwiseMin( m_upper );
    }

    Vector
    projected_gradient( Vector const & x, Vector const & g ) const
    {
      Vector result = g;
      projected_gradient_inplace( x, result );
      return result;
    }
    void
    projected_gradient_inplace( const Vector & x, Vector & g ) const
    {
      const Scalar eps = Scalar( 10 ) * m_epsi;  // tolleranza più robusta

      for ( int i = 0; i < x.size(); ++i )
      {
        bool on_lower = ( x[i] <= m_lower[i] + m_tol_lower[i] );
        bool on_upper = ( x[i] >= m_upper[i] - m_tol_upper[i] );

        // Gradiente numericamente significativo?
        bool grad_neg = ( g[i] < -eps );  // spinge verso il basso
        bool grad_pos = ( g[i] > eps );   // spinge verso l'alto

        // Stato attivo robusto: solo se gradiente spinge fuori in modo
        // significativo
        if ( ( on_lower && grad_neg ) || ( on_upper && grad_pos ) ) { g[i] = Scalar( 0 ); }
      }
    }
    void
    projected_direction_inplace( const Vector & x, Vector & d ) const
    {
      const Scalar eps = Scalar( 10 ) * m_epsi;

      for ( int i = 0; i < x.size(); ++i )
      {
        bool on_lower = ( x[i] <= m_lower[i] + m_tol_lower[i] );
        bool on_upper = ( x[i] >= m_upper[i] - m_tol_upper[i] );

        bool dir_neg = ( d[i] < -eps );
        bool dir_pos = ( d[i] > eps );

        // Annulla solo le componenti che violerebbero seriamente il bound
        if ( ( on_lower && dir_neg ) || ( on_upper && dir_pos ) ) { d[i] = Scalar( 0 ); }
      }
    }

    Scalar
    projected_gradient_norm( Vector const & x, Vector const & g ) const
    {
      return projected_gradient( x, g ).template lpNorm<Eigen::Infinity>();
    }

    // Helper function for descent direction validation
    bool
    is_valid_descent_direction( Scalar pg ) const
    {
      if ( !std::isfinite( pg ) ) return false;
      if ( pg >= -m_epsi * ( 1 + std::abs( pg ) ) ) return false;
      return true;
    }

    // ===========================================================================
    // PRINTING METHODS (LBFGS style - similar to NelderMead)
    // ===========================================================================

    void
    print_optimization_header( size_t n, Scalar f0 ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "{}"
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "{}║                       L-BFGS Optimization                      "
                  "║\n"
                  "{}"
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}"
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_indent, m_indent, m_indent, m_indent, fmt::format( "Dimension: {:d}", n ), m_indent,
                  fmt::format( "Max Iterations: {:d}", m_options.max_iter ), m_indent,
                  fmt::format( "Memory (m): {:d}", m_options.m ), m_indent,
                  fmt::format( "Gradient Tolerance: {:.2e}", m_options.g_tol ), m_indent,
                  fmt::format( "Bounds: {}", ( m_options.use_projection ? "Active" : "None" ) ), m_indent );
      fmt::print( "{}Initial F = {:.6e}\n", m_indent, f0 );
    }

    void
    print_iteration_summary( size_t iter, Scalar f, Scalar gnorm, Scalar step, Scalar pg, bool improved ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      bool show_detailed = m_options.verbosity_level >= 3;

      auto   color = improved ? LBFGS_utils::PrintColors::SUCCESS : LBFGS_utils::PrintColors::WARNING;
      string icon  = improved ? "↗" : "→";

      if ( show_detailed )
      {
        // Versione dettagliata per livello 3+
        fmt::print( color,
                    "{}[{:4d}] {} F = {:<12.6e} | ‖pg‖ = {:<12.6e} | Step = "
                    "{:<12.6e} | pg = {:<12.6e}\n",
                    m_indent, iter, icon, f, gnorm, step, pg );
      }
      else
      {
        // Versione compatta per livello 2
        fmt::print( color,
                    "{}[{:4d}] {} F = {:<12.6e} | ‖pg‖ = {:<12.6e} | Step = "
                    "{:<12.6e}\n",
                    m_indent, iter, icon, f, gnorm, step );
      }
    }

    void
    print_line_search_result( size_t iter, Scalar step, size_t evals, bool success ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      auto   color  = success ? LBFGS_utils::PrintColors::SUCCESS : LBFGS_utils::PrintColors::ERROR;
      string status = success ? "success" : "failure";

      fmt::print( color, "{}Line search: step = {:<12.6e}, evals = {}, {}\n", m_indent, step, evals, status );
    }

    void
    print_lbfgs_update( size_t iter, Scalar sty, bool accepted ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      auto   color  = accepted ? LBFGS_utils::PrintColors::SUCCESS : LBFGS_utils::PrintColors::WARNING;
      string status = accepted ? "accepted" : "rejected";

      fmt::print( color, "{}L-BFGS update: sᵀy = {:<12.6e}, {}\n", m_indent, sty, status );
    }

    void
    print_convergence_info( Status status, Scalar gnorm, Scalar f_change, Scalar x_change ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      auto color = LBFGS_utils::PrintColors::INFO;
      fmt::print( color,
                  "{}Convergence: status = {}, ‖pg‖ = {:.2e}, Δf = {:.2e}, Δx "
                  "= {:.2e}\n",
                  m_indent, static_cast<int>( status ), gnorm, f_change, x_change );
    }

    void
    print_optimization_statistics( Result const & data ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      string status_str;
      auto   status_color = LBFGS_utils::PrintColors::INFO;

      status_str = to_string( data.status );
      switch ( data.status )
      {
        case Status::CONVERGED:
          status_color = LBFGS_utils::PrintColors::SUCCESS;
          break;
        case Status::GRADIENT_TOO_SMALL:
          status_color = LBFGS_utils::PrintColors::SUCCESS;
          break;
        case Status::MAX_ITERATIONS:
          status_color = LBFGS_utils::PrintColors::WARNING;
          break;
        case Status::LINE_SEARCH_FAILED:
          status_color = LBFGS_utils::PrintColors::ERROR;
          break;
        default:
          status_color = LBFGS_utils::PrintColors::ERROR;
      }

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "{}"
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "{}║                    Optimization Finished                       "
                  "║\n"
                  "{}"
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "{}║  Final Status       : {:<39}  ║\n"
                  "{}║  Final Value        : {:<39.6e}  ║\n"
                  "{}║  Initial Value      : {:<39.6e}  ║\n"
                  "{}║  Total Iterations   : {:<39}  ║\n"
                  "{}║  Total Evals        : {:<39}  ║\n"
                  "{}║  Final ‖pg‖         : {:<39.6e}  ║\n"
                  "{}║  L-BFGS Memory      : {:<39}  ║\n"
                  "{}"
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_indent, m_indent, m_indent, m_indent, status_str, m_indent, data.final_function_value, m_indent,
                  data.initial_function_value, m_indent, data.total_iterations, m_indent, data.total_evaluations,
                  m_indent, data.final_gradient_norm, m_indent, m_LBFGS.size(), m_indent );
    }

    // Mutable state variables
    mutable Vector m_x, m_g, m_p, m_x_new, m_g_new, m_s, m_y;
    mutable size_t m_iter_since_reset{ 0 };
    mutable size_t m_function_evaluations{ 0 };
    mutable size_t m_line_search_evaluations{ 0 };  // Aggiunto per tracciare eval line search

  public:
    LBFGS_minimizer( Options opts = Options() ) : m_options( opts ), m_LBFGS( opts.m ) {}

    void
    set_indent( string const & indent )
    {
      m_indent = indent;
    }

    Vector const &
    solution() const
    {
      return m_x;
    }

    void
    set_bounds( size_t n, Scalar const lower[], Scalar const upper[] )
    {
      m_lower.resize( n );
      m_upper.resize( n );
      std::copy_n( lower, n, m_lower.data() );
      std::copy_n( upper, n, m_upper.data() );

      m_tol_lower.resize( n );
      m_tol_upper.resize( n );

      // Calcola tolleranze vettoriali
      m_tol_lower = m_epsi * ( Vector::Ones( n ).array() + m_lower.array().abs() );
      m_tol_upper = m_epsi * ( Vector::Ones( n ).array() + m_upper.array().abs() );

      m_options.use_projection = true;
      check_bounds_consistency();
    }

    void
    set_bounds( Vector const & lower, Vector const & upper )
    {
      assert( lower.size() == upper.size() );
      set_bounds( lower.size(), lower.data(), upper.data() );
    }

    void
    reset_memory()
    {
      if ( m_options.verbosity_level >= 2 ) fmt::print( "{}[LBFGS] Periodic memory reset\n", m_indent );
      m_LBFGS.clear();
      m_iter_since_reset = 0;
    }

    template <typename Linesearch>
    Result
    minimize( Vector const &     x0,
              Callback const &   callback,
              Linesearch const & linesearch = MoreThuenteLineSearch<Scalar>() )
    {
      Status status{ Status::MAX_ITERATIONS };
      Scalar gnorm{ 0 };

      // Reset counters
      m_function_evaluations    = 0;
      m_line_search_evaluations = 0;

      auto const n{ x0.size() };
      m_x.resize( n );
      m_g.resize( n );
      m_x_new.resize( n );
      m_g_new.resize( n );
      m_p.resize( n );
      m_s.resize( n );
      m_y.resize( n );

      if ( m_options.use_projection )
      {
        assert( m_lower.size() == n );
        assert( m_upper.size() == n );
        check_bounds_consistency();
      }

      // Initialize and project if needed
      m_x.noalias() = x0;
      if ( m_options.use_projection ) project_inplace( m_x );

      // Initial evaluation
      Scalar f = callback( m_x, &m_g );
      m_function_evaluations++;

      Scalar f_prev{ f };
      Scalar f_initial{ f };
      size_t iteration{ 0 };

      print_optimization_header( n, f_initial );

      // Main optimization loop
      for ( ; iteration < m_options.max_iter; ++iteration )
      {
        ++m_iter_since_reset;

        // Check gradient norm
        gnorm = projected_gradient_norm( m_x, m_g );

        if ( gnorm <= m_options.g_tol )
        {
          status = Status::GRADIENT_TOO_SMALL;
          goto exit_position;
        }

        // Periodic memory reset
        if ( m_iter_since_reset >= m_options.iter_reset ) reset_memory();

        // Compute search direction using L-BFGS
        Scalar h0 = m_LBFGS.compute_initial_h0( 1 );
        m_p       = -m_LBFGS.two_loop_recursion( m_g, h0 );

        // Project direction if bounds are active
        if ( m_options.use_projection )
        {
          projected_direction_inplace( m_x, m_p );
          if ( m_p.isZero( m_epsi ) ) m_p = -projected_gradient( m_x, m_g );
        }

        // Robust descent direction check with fallback strategies
        Scalar pg = m_p.dot( m_g );
        if ( !is_valid_descent_direction( pg ) )
        {
          // Try gradient direction
          m_p = -m_g;
          if ( m_options.use_projection ) projected_direction_inplace( m_x, m_p );
          pg = m_p.dot( m_g );
          if ( !is_valid_descent_direction( pg ) )
          {
            // Reset L-BFGS memory and try again
            m_LBFGS.clear();
            h0  = m_LBFGS.compute_initial_h0( 1 );
            m_p = -m_LBFGS.two_loop_recursion( m_g, h0 );
            if ( m_options.use_projection ) projected_direction_inplace( m_x, m_p );
            pg = m_p.dot( m_g );
            if ( m_options.verbosity_level >= 1 )
              fmt::print( LBFGS_utils::PrintColors::ERROR, "{}[LBFGS] Cannot find descent direction, stopping\n",
                          m_indent );
            gnorm  = projected_gradient_norm( m_x, m_g );
            status = Status::FAILED;
            goto exit_position;
          }
        }

        // Line search
        auto step_opt = linesearch( f, pg, m_x, m_p, callback, m_options.step_max );

        if ( !step_opt.has_value() )
        {
          if ( m_options.verbosity_level >= 1 )
            fmt::print( LBFGS_utils::PrintColors::WARNING, "{}[LBFGS] line search failed, trying fallback steps\n",
                        m_indent );

          m_LBFGS.clear();

          // Try progressively smaller fixed steps
          bool                        fallback_success = false;
          const std::array<Scalar, 3> fallback_steps   = { m_options.very_small_step,
                                                           m_options.very_small_step * Scalar( 0.1 ),
                                                           m_options.very_small_step * Scalar( 0.01 ) };

          for ( Scalar fallback_step : fallback_steps )
          {
            m_x_new.noalias() = m_x + fallback_step * m_p;
            if ( m_options.use_projection ) project_inplace( m_x_new );
            Scalar f_test = callback( m_x_new, &m_g_new );
            m_function_evaluations++;

            if ( f_test < f )
            {
              m_s.noalias() = m_x_new - m_x;
              m_y.noalias() = m_g_new - m_g;
              m_x.swap( m_x_new );
              m_g.swap( m_g_new );
              f                = f_test;
              fallback_success = true;
              print_line_search_result( iteration, fallback_step, 1, true );
              break;
            }
          }

          if ( !fallback_success )
          {
            gnorm  = projected_gradient_norm( m_x, m_g );
            status = Status::LINE_SEARCH_FAILED;
            goto exit_position;
          }
          continue;
        }

        auto [step, n_evals] = *step_opt;
        m_function_evaluations += n_evals;
        m_line_search_evaluations += n_evals;

        print_line_search_result( iteration, step, n_evals, true );

        // Evaluate new point
        m_x_new.noalias() = m_x + step * m_p;
        if ( m_options.use_projection ) project_inplace( m_x_new );
        Scalar f_new = callback( m_x_new, &m_g_new );
        m_function_evaluations++;

        // Compute differences for L-BFGS update
        m_s.noalias() = m_x_new - m_x;
        m_y.noalias() = m_g_new - m_g;

        // Robust curvature check for L-BFGS update
        Scalar sty           = m_s.dot( m_y );
        Scalar sty_tolerance = std::max( m_options.sty_min_factor, m_epsi * m_s.squaredNorm() * m_y.squaredNorm() );

        bool update_accepted = false;
        if ( sty > sty_tolerance ) { update_accepted = m_LBFGS.add_correction( m_s, m_y ); }

        print_lbfgs_update( iteration, sty, update_accepted );

        if ( !update_accepted && m_options.verbosity_level >= 2 && m_LBFGS.size() > 0 )
        {
          fmt::print( LBFGS_utils::PrintColors::WARNING,
                      "{}[LBFGS] curvature condition failed (s^T y = {:.2e}), "
                      "skipping update\n",
                      m_indent, sty );
        }

        // Move to new point
        m_x.swap( m_x_new );
        m_g.swap( m_g_new );
        Scalar f_old = f;
        f            = f_new;

        // Check function value change
        Scalar f_change = std::abs( f - f_old );
        if ( f_change <= m_options.f_tol )
        {
          gnorm = projected_gradient_norm( m_x, m_g );
          if ( gnorm <= m_options.g_tol_weak )
          {
            if ( m_options.verbosity_level >= 1 )
              fmt::print( LBFGS_utils::PrintColors::SUCCESS,
                          "{}[LBFGS] Converged by function change: {:.2e} < {:.2e}\n", m_indent, f_change,
                          m_options.f_tol );
            status = Status::CONVERGED;
            goto exit_position;
          }
        }

        // Check for stagnation in variables
        {
          Scalar x_change = m_s.template lpNorm<Eigen::Infinity>();
          if ( x_change < m_options.x_tol )
          {
            gnorm = projected_gradient_norm( m_x, m_g );
            if ( gnorm <= m_options.g_tol_weak )
            {
              if ( m_options.verbosity_level >= 1 )
                fmt::print( LBFGS_utils::PrintColors::SUCCESS, "{}[LBFGS] Converged by x change: {:.2e} < {:.2e}\n",
                            m_indent, x_change, m_options.x_tol );
              status = Status::CONVERGED;
              goto exit_position;
            }
          }
        }

        // Print iteration summary
        print_iteration_summary( iteration, f, gnorm, step, pg, f < f_old );
        f_prev = f;
      }

      // Final gradient norm if not converged
      gnorm = projected_gradient_norm( m_x, m_g );

    exit_position:
      Result result{ status, iteration, m_function_evaluations, gnorm, f, f_initial };

      // Aggiungi informazioni line search al risultato
      result.line_search_evaluations = m_line_search_evaluations;

      print_optimization_statistics( result );

      return result;
    }
  };

  // ===========================================================================
  // CLASS: LBFGS_BlockCoordinate (Enhanced Block Coordinate L-BFGS)
  // ===========================================================================

  /**
   * @class LBFGS_BlockCoordinate
   * @brief Enhanced block coordinate L-BFGS with consecutive interval blocks
   *
   * This version implements a single block selection strategy based on
   * consecutive intervals that cover the entire coordinate space with
   * configurable overlap between blocks.
   *
   * ### Block Selection Strategy
   *
   * Generates intervals [start, end] such that:
   * - The union of all intervals covers [0, n-1] (all coordinates)
   * - Interval sizes range from min_block_size to max_block_size
   * - Statistical overlap between consecutive blocks is approximately
   * overlap_ratio
   * - Blocks are generated as consecutive integer intervals
   *
   * @tparam Scalar Floating-point type (typically double)
   */
  template <typename Scalar = double>
  class LBFGS_BlockCoordinate
  {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar( Vector const &, Vector * )>;

    enum class Status
    {
      CONVERGED,
      MAX_OUTER_ITERATIONS,
      MAX_INNER_ITERATIONS,
      GRADIENT_TOO_SMALL,
      LINE_SEARCH_FAILED,
      FAILED
    };

    static string
    to_string( Status status )
    {
      switch ( status )
      {
        case Status::CONVERGED:
          return "CONVERGED";
        case Status::MAX_OUTER_ITERATIONS:
          return "MAX_OUTER_ITERATIONS";
        case Status::MAX_INNER_ITERATIONS:
          return "MAX_INNER_ITERATIONS";
        case Status::LINE_SEARCH_FAILED:
          return "LINE_SEARCH_FAILED";
        case Status::FAILED:
          return "FAILED";
        default:
          return "UNKNOWN";
      }
    }

    struct Result
    {
      Scalar final_function_value{ 0 };
      Scalar initial_function_value{ 0 };
      Status status{ Status::FAILED };

      // Statistics
      size_t outer_iterations{ 0 };
      size_t inner_iterations{ 0 };
      size_t total_iterations{ 0 };

      size_t outer_evaluations{ 0 };
      size_t inner_evaluations{ 0 };
      size_t total_evaluations{ 0 };

      // Enhanced statistics
      size_t blocks_processed{ 0 };
      Scalar max_block_improvement{ 0 };
    };

    struct Options
    {
      // Block configuration
      size_t block_size{ 20 };      ///< Base block size
      double overlap_ratio{ 0.3 };  ///< Fraction of block overlap (0-1)
      size_t min_block_size{ 5 };   ///< Minimum block size
      size_t max_block_size{ 15 };  ///< Maximum block size

      // Outer loop control
      size_t max_outer_iterations{ 50 };
      size_t max_inner_iterations{ 200 };
      Scalar outer_tolerance{ 1e-6 };

      // Inner L-BFGS options
      size_t lbfgs_m{ 10 };
      Scalar lbfgs_g_tol{ 1e-6 };
      Scalar lbfgs_f_tol{ 1e-12 };  // Much tighter for inner convergence

      // Adaptive parameters
      bool   adaptive_block_size{ true };
      Scalar progress_threshold{ 1e-6 };            // Relative threshold for block acceptance
      Scalar absolute_progress_threshold{ 1e-12 };  // Absolute threshold for tiny improvements

      // Convergence control
      size_t max_stagnation_count{ 3 };
      Scalar relative_stagnation_tol{ 1e-8 };  // Tighter relative tolerance

      // Verbosity (enhanced)
      size_t verbosity_level{ 1 };  // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      bool   use_unicode_borders{ true };
    };

  private:
    Options m_options;
    Vector  m_lower;
    Vector  m_upper;
    Vector  m_x;
    bool    m_use_bounds{ false };

    size_t m_outer_iteration_count{ 0 };

    // Random number generation
    std::random_device m_rd;
    std::mt19937       m_gen;

    // Progress tracking
    Scalar m_previous_best{ std::numeric_limits<Scalar>::max() };
    size_t m_stagnation_count{ 0 };

    string m_indent{ "" };  // Per indentazione consistente

    // ===========================================================================
    // PRINTING METHODS for Block Coordinate
    // ===========================================================================

    void
    print_outer_iteration_header( size_t                 outer_iter,
                                  size_t                 total_cycles,
                                  vector<size_t> const & block_indices,
                                  size_t                 block_size ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "\n"
                  "    "
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "    ║ Outer Iteration {:3d} - Block {:2d}/{:2d}                     "
                  "  "
                  "       ║\n"
                  "    "
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "    ║ Size: {:<5} Indices: {:<41} ║\n"
                  "    "
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  outer_iter, ( outer_iter % total_cycles ) + 1, total_cycles, block_size,
                  LBFGS_utils::format_index_vector_compact<size_t>( block_indices, 7 ) );
    }

    void
    print_outer_statistics( size_t total_outer_iters,
                            size_t total_inner_iters,
                            size_t total_evals,
                            Scalar best_value,
                            bool   converged ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::INFO,
                  "\n"
                  "  "
                  "╔═══════════════════════════════════════════════════════════════╗\n"
                  "  ║                      Outer Iteration Summary                  "
                  "║\n"
                  "  "
                  "╠═══════════════════════════════════════════════════════════════╣\n"
                  "  ║ Completed:        {:<43} ║\n"
                  "  ║ Outer Iterations: {:<43} ║\n"
                  "  ║ Inner Iterations: {:<43} ║\n"
                  "  ║ Total Evals:      {:<43} ║\n"
                  "  ║ Best Value:       {:<43.6e} ║\n"
                  "  ║ Status:           {:<43} ║\n"
                  "  "
                  "╚═══════════════════════════════════════════════════════════════╝\n"
                  "\n",
                  fmt::format( "{}/{}", m_outer_iteration_count, total_outer_iters ), m_outer_iteration_count,
                  total_inner_iters, total_evals, best_value, ( converged ? "CONVERGED" : "RUNNING" ) );
    }

    void
    print_block_coordinate_header( size_t n, Scalar f0, Scalar g0 ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "║                   Block Coordinate L-BFGS                      ║\n"
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "║ {:62} ║\n"
                  "║ {:62} ║\n"
                  "║ {:62} ║\n"
                  "║ {:62} ║\n"
                  "║ {:62} ║\n"
                  "║ {:62} ║\n"
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  fmt::format( "Dimension:            {}", n ),
                  fmt::format( "Block Size:           {}", m_options.block_size ),
                  fmt::format( "Max Outer Iterations: {}", m_options.max_outer_iterations ),
                  fmt::format( "Overlap Ratio:        {:.4e}", m_options.overlap_ratio ),
                  fmt::format( "Initial f:            {:.6e}", f0 ),
                  fmt::format( "‖g‖∞                  {:.6e}", g0 ) );
    }

    // ===========================================================================
    // ENHANCED PRINTING METHODS
    // ===========================================================================

    void
    print_outer_iteration_start( Scalar f, Scalar gnorm ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "\n"
                  "    "
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "    ║ OUTER ITERATION {:3d} START F = {:<12.6e} ‖pg‖ = {:<12.6e} ║\n"
                  "    "
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_outer_iteration_count, f, gnorm );
    }

    void
    print_outer_iteration_end( Scalar f, Scalar gnorm, Scalar improvement ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( LBFGS_utils::PrintColors::HEADER,
                  "\n"
                  "    "
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "    ║ OUTER ITERATION {:3d} END                                     "
                  "  "
                  " ║\n"
                  "    ║ F = {:<12.6e}    ‖pg‖ = {:<12.6e}    ΔF = {:<12.6e}   ║\n"
                  "    "
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_outer_iteration_count, f, gnorm, improvement );
    }

    void
    print_blocks_summary( std::vector<std::vector<size_t>> const & blocks ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      fmt::print( LBFGS_utils::PrintColors::INFO,
                  "    "
                  "┌────────────────────────────────────────────────────────────────┐\n"
                  "    │                      BLOCKS FOR ITERATION {:3d}               "
                  "  "
                  " │\n"
                  "    "
                  "├────────────────────────────────────────────────────────────────┤"
                  "\n",
                  m_outer_iteration_count );

      for ( size_t i = 0; i < blocks.size(); ++i )
      {
        fmt::print( LBFGS_utils::PrintColors::INFO, "    │ Block {:3d}: Size={:3d}  Indices: {:<32} │\n", i + 1,
                    blocks[i].size(), LBFGS_utils::format_index_vector_compact<size_t>( blocks[i], 5 ) );
      }

      fmt::print( LBFGS_utils::PrintColors::INFO,
                  "    "
                  "└───────────────────────────────────────────────────────────"
                  "─────┘\n" );
    }

    void
    print_block_result_detail( size_t block_idx,
                               size_t total_blocks,
                               Scalar f_before,
                               Scalar f_after,
                               size_t inner_iters,
                               size_t inner_evals,
                               bool   success ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      Scalar improvement = f_before - f_after;
      auto   color       = success && improvement > 0 ? LBFGS_utils::PrintColors::SUCCESS
                           : success                  ? LBFGS_utils::PrintColors::WARNING
                                                      : LBFGS_utils::PrintColors::ERROR;

      string status_icon = success ? ( improvement > 0 ? "✅" : "⚡" ) : "❌";
      string status_text = success ? ( improvement > 0 ? "IMPROVED" : "NO CHANGE" ) : "FAILED";

      fmt::print( color,
                  "      ┌────────────────────── BLOCK {:>2d}/{:<2d} RESULT "
                  "──────────────────────┐\n"
                  "      │ {} Status: {:<10}   Inner Iters: {:<6}   Evals: {:<6}    │\n"
                  "      │ F: {:>16.8e} → {:<16.8e}  ΔF: {:<16.8e}   │\n"
                  "      "
                  "└────────────────────────────────────────────────────────────────┘"
                  "\n",
                  block_idx + 1, total_blocks, status_icon, status_text, inner_iters, inner_evals, f_before, f_after,
                  improvement );
    }

    void
    print_block_header( size_t                      block_idx,
                        size_t                      total_blocks,
                        std::vector<size_t> const & block_indices,
                        size_t                      block_size ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      fmt::print( LBFGS_utils::PrintColors::INFO,
                  "\n"
                  "      "
                  "┌────────────────────────────────────────────────────────────────┐\n"
                  "      │ Block: {:>4d}/{:<4d}    Iteration: {:<4d}                   "
                  "  "
                  "       │\n"
                  "      │ Size: {:<4} Indices: {:<42} │\n"
                  "      "
                  "└────────────────────────────────────────────────────────────────┘"
                  "\n",
                  block_idx + 1, total_blocks, m_outer_iteration_count, block_size,
                  LBFGS_utils::format_index_vector_compact<size_t>( block_indices, 5 ) );
    }

    // ===========================================================================
    // CONSECUTIVE INTERVAL BLOCK GENERATION
    // ===========================================================================

    /**
     * @brief Generate consecutive interval blocks covering [0, n-1] with
     * guaranteed overlap
     *
     * Generates blocks as consecutive integer intervals [start, end] such that:
     * - The union of all intervals covers [0, n-1]
     * - Interval sizes are exactly block_size (except possibly the last one)
     * - Consecutive blocks overlap by approximately overlap_ratio * block_size
     *
     * @param n_total Total number of coordinates
     * @return Vector of blocks, each block is a vector of consecutive indices
     */
    std::vector<std::vector<size_t>>
    generate_consecutive_interval_blocks( size_t n_total )
    {
      std::vector<std::vector<size_t>> blocks;

      if ( n_total == 0 ) return blocks;

      size_t block_size = std::clamp( m_options.block_size, m_options.min_block_size,
                                      std::min( m_options.max_block_size, n_total ) );

      // Calculate step size based on overlap ratio
      size_t step = static_cast<size_t>( block_size * ( 1.0 - m_options.overlap_ratio ) );
      step        = std::max<size_t>( step, 1 );

      // If step >= block_size, we get no overlap - prevent this
      if ( step >= block_size )
      {
        step = block_size / 2;
        if ( step < 1 ) step = 1;
      }

      // Generate blocks with overlap
      size_t start = 0;
      while ( start < n_total )
      {
        size_t              end = std::min( start + block_size - 1, n_total - 1 );
        std::vector<size_t> block;
        for ( size_t i = start; i <= end; ++i ) { block.push_back( i ); }
        blocks.push_back( block );

        // If we have reached the end, break
        if ( end == n_total - 1 ) break;

        start += step;
      }

      // Now, if the last block is too small, merge it with the previous one
      if ( blocks.size() > 1 )
      {
        auto & last_block = blocks.back();
        if ( last_block.size() < m_options.min_block_size )
        {
          auto & prev_block = blocks[blocks.size() - 2];
          // Merge the last block into the previous one
          prev_block.insert( prev_block.end(), last_block.begin(), last_block.end() );
          // Remove duplicates? The blocks are consecutive, so no duplicates?
          // Actually, they might overlap. Sort and remove duplicates?
          std::sort( prev_block.begin(), prev_block.end() );
          prev_block.erase( std::unique( prev_block.begin(), prev_block.end() ), prev_block.end() );
          // Now remove the last block
          blocks.pop_back();
        }
      }

      return blocks;
    }

    /**
     * @brief Adjust block size based on progress and coverage
     */
    void
    adapt_block_size( Scalar improvement, size_t n_total )
    {
      if ( !m_options.adaptive_block_size ) return;

      Scalar relative_improvement = std::abs( improvement ) / ( 1.0 + std::abs( m_previous_best ) );

      // Adjust based on progress
      if ( relative_improvement > 0.1 )
      {
        // Good progress - increase block size
        m_options.block_size = std::min( m_options.max_block_size, m_options.block_size + 5 );
      }
      else if ( relative_improvement < 0.001 )
      {
        // Poor progress - decrease block size
        m_options.block_size = std::max( m_options.min_block_size, m_options.block_size - 2 );
      }

      // Clamp to valid range
      m_options.block_size = std::clamp( m_options.block_size, m_options.min_block_size,
                                         std::min( m_options.max_block_size, n_total ) );

      m_previous_best = improvement;

      if ( m_options.verbosity_level >= 2 )
      {
        fmt::print( "    [BlockLBFGS] Adjusted block_size to {}\n", m_options.block_size );
      }
    }

    Vector
    extract_subvector( Vector const & full, std::vector<size_t> const & indices )
    {
      Vector sub( indices.size() );
      for ( size_t i = 0; i < indices.size(); ++i ) { sub( i ) = full( indices[i] ); }
      return sub;
    }

    void
    update_full_vector( Vector & full, Vector const & sub, std::vector<size_t> const & indices )
    {
      for ( size_t i = 0; i < indices.size(); ++i ) { full( indices[i] ) = sub( i ); }
    }

    void
    extract_subbounds( std::vector<size_t> const & indices, Vector & sub_lower, Vector & sub_upper )
    {
      sub_lower.resize( indices.size() );
      sub_upper.resize( indices.size() );

      for ( size_t i = 0; i < indices.size(); ++i )
      {
        sub_lower( i ) = m_lower( indices[i] );
        sub_upper( i ) = m_upper( indices[i] );
      }
    }

    void
    project_point( Vector & x ) const
    {
      if ( m_use_bounds ) { x = x.cwiseMax( m_lower ).cwiseMin( m_upper ); }
    }

    Scalar
    projected_gradient_norm( Vector const & x, Vector const & g ) const
    {
      if ( !m_use_bounds ) return g.template lpNorm<Eigen::Infinity>();

      Vector pg = g;

      // Zero gradient components for variables at bounds pointing outward
      for ( Eigen::Index i = 0; i < x.size(); ++i )
      {
        if ( ( x( i ) <= m_lower( i ) + 1e-12 && g( i ) > 0 ) || ( x( i ) >= m_upper( i ) - 1e-12 && g( i ) < 0 ) )
        {
          pg( i ) = 0;
        }
      }

      return pg.template lpNorm<Eigen::Infinity>();
    }

  public:
    explicit LBFGS_BlockCoordinate( Options const & opts = Options() ) : m_options( opts ), m_gen( m_rd() ) {}

    void
    set_bounds( Vector const & lower, Vector const & upper )
    {
      UTILS_ASSERT( lower.size() == upper.size(),
                    "BlockLBFGS::set_bounds: lower and upper bounds must have "
                    "same dimension" );

      UTILS_ASSERT( ( lower.array() <= upper.array() ).all(),
                    "BlockLBFGS::set_bounds: lower bounds must be <= upper bounds "
                    "for all coordinates" );

      m_lower      = lower;
      m_upper      = upper;
      m_use_bounds = true;
    }

    void
    clear_bounds()
    {
      m_use_bounds = false;
    }

    Vector const &
    solution() const
    {
      return m_x;
    }

    /**
     * @brief Enhanced optimization routine with consecutive interval blocks
     */
    template <typename Linesearch>
    Result
    minimize( Vector const &     x0,
              Callback const &   global_callback,
              Linesearch const & linesearch = MoreThuenteLineSearch<Scalar>() )
    {
      Result result;
      size_t n = x0.size();

      // If dimension <= block_size, use standard LBFGS
      if ( n <= m_options.block_size )
      {
        if ( m_options.verbosity_level >= 1 )
        {
          fmt::print(
              "    [BlockLBFGS] Problem dimension {} <= block_size {}, "
              "using standard LBFGS\n",
              n, m_options.block_size );
        }

        // Use standard LBFGS instead of block coordinate
        typename LBFGS_minimizer<Scalar>::Options opts;
        opts.max_iter        = m_options.max_outer_iterations * m_options.max_inner_iterations;
        opts.g_tol           = m_options.outer_tolerance;
        opts.m               = m_options.lbfgs_m;
        opts.verbosity_level = m_options.verbosity_level;

        LBFGS_minimizer<Scalar> minimizer( opts );
        if ( m_use_bounds ) { minimizer.set_bounds( m_lower, m_upper ); }

        auto standard_result = minimizer.minimize( x0, global_callback, linesearch );

        m_x = minimizer.solution();

        // Convert result
        result.final_function_value   = standard_result.final_function_value;
        result.initial_function_value = standard_result.initial_function_value;
        result.status                 = ( standard_result.status == LBFGS_minimizer<Scalar>::Status::CONVERGED ||
                          standard_result.status == LBFGS_minimizer<Scalar>::Status::GRADIENT_TOO_SMALL )
                                            ? Status::CONVERGED
                                            : Status::FAILED;
        result.outer_iterations       = 1;
        result.inner_iterations       = standard_result.total_iterations;
        result.total_iterations       = standard_result.total_iterations;
        result.total_evaluations      = standard_result.total_evaluations;
        result.blocks_processed       = 1;
        result.max_block_improvement  = result.initial_function_value - result.final_function_value;

        return result;
      }

      // Initialize solution
      m_x.resize( x0.size() );
      m_x = x0;
      if ( m_use_bounds ) project_point( m_x );

      // Initial evaluation
      Vector g( n );
      Scalar f                      = global_callback( m_x, &g );
      result.initial_function_value = f;

      size_t total_evals           = 1;
      size_t total_inner_iters     = 0;
      size_t blocks_processed      = 0;
      Scalar max_block_improvement = 0;

      print_block_coordinate_header( n, f, g.template lpNorm<Eigen::Infinity>() );

      // Reset stagnation counter
      m_stagnation_count = 0;
      m_previous_best    = f;

      // Main outer loop
      m_outer_iteration_count = 0;
      while ( m_outer_iteration_count < m_options.max_outer_iterations )
      {
        ++m_outer_iteration_count;
        Scalar f_start_cycle = f;

        // STAMPA INIZIO OUTER ITERATION
        print_outer_iteration_start( f, projected_gradient_norm( m_x, g ) );

        // Generate consecutive interval blocks for this iteration
        auto all_blocks = generate_consecutive_interval_blocks( n );
        print_blocks_summary( all_blocks );

        // Track if any block made progress
        bool   made_progress_in_cycle    = false;
        Scalar best_improvement_in_cycle = 0;

        // Process each block
        for ( size_t block_idx = 0; block_idx < all_blocks.size(); ++block_idx )
        {
          auto block_indices = all_blocks[block_idx];

          if ( block_indices.empty() ) continue;

          // STAMPA INIZIO BLOCCO
          print_block_header( block_idx, all_blocks.size(), block_indices, block_indices.size() );

          // Extract block subproblem
          Vector x_block = extract_subvector( m_x, block_indices );
          Vector g_block = extract_subvector( g, block_indices );

          // Check if block gradient is significant
          Scalar block_gnorm = g_block.template lpNorm<Eigen::Infinity>();

          // Skip if block gradient is too small relative to global gradient
          if ( block_gnorm < m_options.lbfgs_g_tol * 0.1 )
          {
            if ( m_options.verbosity_level >= 2 )
            {
              fmt::print(
                  "    [BlockLBFGS] Block gradient too small ({:.2e} < "
                  "{:.2e}), skipping\n",
                  block_gnorm, m_options.lbfgs_g_tol * 0.1 );
            }
            continue;
          }

          // Create block callback
          auto block_callback = [&]( Vector const & x_sub, Vector * g_sub ) -> Scalar
          {
            Vector x_full = m_x;
            update_full_vector( x_full, x_sub, block_indices );

            if ( m_use_bounds ) project_point( x_full );

            Vector g_full( n );
            Scalar f_val = global_callback( x_full, &g_full );
            total_evals++;

            if ( g_sub ) { *g_sub = extract_subvector( g_full, block_indices ); }

            return f_val;
          };

          // Scala le tolleranze in base al valore corrente
          Scalar f_scale        = std::max( std::abs( f ), Scalar( 1.0 ) );
          Scalar adaptive_f_tol = std::max( m_options.lbfgs_f_tol, m_options.absolute_progress_threshold / f_scale );

          Scalar global_gnorm = g.template lpNorm<Eigen::Infinity>();

          // Set up inner L-BFGS
          typename LBFGS_minimizer<Scalar>::Options inner_opts;
          inner_opts.max_iter        = m_options.max_inner_iterations;
          inner_opts.g_tol           = std::max( m_options.lbfgs_g_tol,
                                                 global_gnorm * Scalar( 0.1 ) );  // m_options.lbfgs_g_tol;
          inner_opts.f_tol           = m_options.lbfgs_f_tol;
          inner_opts.m               = m_options.lbfgs_m;
          inner_opts.verbosity_level = max( m_options.verbosity_level, static_cast<size_t>( 1 ) );

          // Improved inner L-BFGS options
          inner_opts.step_max        = 1.0;
          inner_opts.very_small_step = 1e-10;
          inner_opts.iter_reset      = std::min( m_options.lbfgs_m, size_t( 10 ) );
          inner_opts.x_tol           = 1e-12;
          inner_opts.sty_min_factor  = 1e-14;

          // Debug print before inner optimization
          if ( m_options.verbosity_level >= 3 )
          {
            fmt::print(
                "    [BlockLBFGS] Starting inner L-BFGS on block [{}] "
                "(size: {}, ‖g‖∞: {:.2e})\n",
                LBFGS_utils::format_index_vector_compact<size_t>( block_indices, 10 ), block_indices.size(),
                g_block.template lpNorm<Eigen::Infinity>() );
          }

          LBFGS_minimizer<Scalar> inner_minimizer( inner_opts );
          inner_minimizer.set_indent( "      " );  // 6 spaces indentation

          // Set bounds for block if needed
          if ( m_use_bounds )
          {
            Vector block_lower, block_upper;
            extract_subbounds( block_indices, block_lower, block_upper );
            inner_minimizer.set_bounds( block_lower, block_upper );
          }

          // Save current state for comparison
          Vector x_before_block = m_x;
          Scalar f_before_block = f;

          // Optimize block
          auto inner_result = inner_minimizer.minimize( x_block, block_callback, linesearch );

          total_inner_iters += inner_result.total_iterations;
          blocks_processed++;

          // Update solution
          update_full_vector( m_x, inner_minimizer.solution(), block_indices );
          if ( m_use_bounds ) project_point( m_x );

          // Re-evaluate at new point
          f = global_callback( m_x, &g );
          total_evals++;

          // Track block improvement
          Scalar block_improvement = f_before_block - f;
          max_block_improvement    = std::max( max_block_improvement, block_improvement );

          // Enhanced acceptance criteria: both relative AND absolute
          Scalar relative_improvement = block_improvement / ( std::abs( f_before_block ) + Scalar( 1e-16 ) );
          Scalar absolute_improvement = block_improvement;

          bool block_made_progress = ( relative_improvement > m_options.progress_threshold ) ||
                                     ( absolute_improvement > m_options.absolute_progress_threshold );

          // Check if this block made significant progress
          // bool block_made_progress = (block_improvement >
          // m_options.progress_threshold * std::abs(f_before_block));
          if ( block_made_progress )
          {
            made_progress_in_cycle    = true;
            best_improvement_in_cycle = std::max( best_improvement_in_cycle, block_improvement );
          }

          // STAMPA RISULTATO DETTAGLIATO BLOCCO
          print_block_result_detail( block_idx, all_blocks.size(), f_before_block, f, inner_result.total_iterations,
                                     inner_result.total_evaluations, block_made_progress );
        }

        // STAMPA FINE OUTER ITERATION
        Scalar total_improvement = f_start_cycle - f;
        print_outer_iteration_end( f, projected_gradient_norm( m_x, g ), total_improvement );

        // Adaptive block size adjustment
        Scalar f_improvement = f_start_cycle - f;

        // More conservative adaptive block sizing
        if ( made_progress_in_cycle )
        {
          // Good progress - consider increasing block size if we're making
          // progress
          if ( best_improvement_in_cycle > 0.01 * std::abs( f_start_cycle ) )
          {
            m_options.block_size = std::min( m_options.max_block_size, m_options.block_size + 2 );
          }
        }
        else
        {
          // Poor progress - decrease block size more cautiously
          m_options.block_size = std::max( m_options.min_block_size, m_options.block_size - 1 );
        }

        // Clamp to valid range
        m_options.block_size = std::clamp( m_options.block_size, m_options.min_block_size,
                                           std::min( m_options.max_block_size, n ) );

        if ( m_options.verbosity_level >= 2 && m_options.adaptive_block_size )
        {
          fmt::print( "    [BlockLBFGS] Adjusted block_size to {}\n", m_options.block_size );
        }

        // Convergence checking
        Scalar gnorm = projected_gradient_norm( m_x, g );

        // Calculate current coverage for display
        // Print outer statistics periodically
        if ( m_options.verbosity_level >= 1 &&
             ( m_outer_iteration_count % 5 == 0 || m_outer_iteration_count == m_options.max_outer_iterations - 1 ) )
        {
          print_outer_statistics( m_options.max_outer_iterations, total_inner_iters, total_evals, f, false );
        }

        // Improved convergence checking
        bool converged = false;

        // Check gradient norm convergence
        if ( gnorm < m_options.outer_tolerance )
        {
          if ( m_options.verbosity_level >= 1 )
          {
            fmt::print(
                "    [BlockLBFGS] Converged by gradient norm: {:.2e} < "
                "{:.2e}\n",
                gnorm, m_options.outer_tolerance );
          }
          result.status = Status::GRADIENT_TOO_SMALL;
          converged     = true;
        }

        // Check function value stagnation with improved criteria
        Scalar cycle_improvement   = f_start_cycle - f;
        Scalar relative_stagnation = cycle_improvement / ( std::abs( f_start_cycle ) + Scalar( 1e-16 ) );

        // Scalar relative_stagnation_tol{1e-8};  // Tighter relative tolerance

        // Check function value stagnation
        Scalar relative_improvement = f_improvement / ( 1.0 + std::abs( f_start_cycle ) );
        if ( !converged && relative_improvement < m_options.progress_threshold )
        {
          ++m_stagnation_count;
          if ( m_stagnation_count >= m_options.max_stagnation_count )
          {
            // Require more consecutive stagnations
            if ( m_options.verbosity_level >= 1 )
            {
              fmt::print(
                  "    [BlockLBFGS] Converged by stagnation ({} "
                  "iterations with Δf < {:.2e})\n",
                  m_stagnation_count, m_options.progress_threshold );
            }
            result.status = Status::CONVERGED;
            converged     = true;
          }
        }
        else
        {
          m_stagnation_count = 0;
        }

        if ( converged ) { break; }

        m_previous_best = f;
      }

      // Final results
      result.final_function_value  = f;
      result.outer_iterations      = m_outer_iteration_count;
      result.inner_iterations      = total_inner_iters;
      result.total_iterations      = result.outer_iterations + total_inner_iters;
      result.outer_evaluations     = total_evals;
      result.inner_evaluations     = total_evals;
      result.total_evaluations     = total_evals;
      result.blocks_processed      = blocks_processed;
      result.max_block_improvement = max_block_improvement;

      if ( result.status == Status::FAILED && result.outer_iterations >= m_options.max_outer_iterations )
      {
        result.status = Status::MAX_OUTER_ITERATIONS;
      }

      // Print final outer statistics
      if ( m_options.verbosity_level >= 1 )
      {
        print_outer_statistics( m_options.max_outer_iterations, total_inner_iters, total_evals, f, true );

        fmt::print(
            "    [BlockLBFGS] Finished: status={}, f_final={:.6e}, "
            "total_evals={}\n",
            static_cast<int>( result.status ), result.final_function_value, total_evals );
      }

      return result;
    }
  };

}  // namespace Utils

#endif

#endif

//
// eof: Utils_LBFGS.hh
//
