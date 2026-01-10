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
// file: Utils_minimize_LBFGS_BlockCoordinate.hh
//

#pragma once

#ifndef UTILS_MINIMIZE_LBFGS_BLOCK_COORDINATE_dot_HH
#define UTILS_MINIMIZE_LBFGS_BLOCK_COORDINATE_dot_HH

#include "Utils_minimize_LBFGS.hh"

/**
 * @file Utils_minimize_LBFGS_BlockCoordinate.hh
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
 * minimizer.minimize(x0, rosenbrock, StrongWolfeLineSearch<double>());
 *
 * std::cout << "Solution: " << minimizer.solution().transpose() << std::endl;
 * std::cout << "Minimum: " << minimizer.final_function_value() << std::endl;
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
 * block_minimizer.minimize(x0, rosenbrock_high_dim);
 *
 * std::cout << "Solution norm: " << block_minimizer.solution().norm() << std::endl;
 * std::cout << "Minimum: " << block_minimizer.final_function_value() << std::endl;
 * @endcode
 *
 * @author Enrico Bertolazzi
 * @date 2025
 */

namespace Utils
{

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
  template <typename Scalar = double> class LBFGS_BlockCoordinate
  {
  public:
    using integer  = Eigen::Index;
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

    static string to_string( Status status )
    {
      switch ( status )
      {
        case Status::CONVERGED: return "CONVERGED";
        case Status::MAX_OUTER_ITERATIONS: return "MAX_OUTER_ITERATIONS";
        case Status::MAX_INNER_ITERATIONS: return "MAX_INNER_ITERATIONS";
        case Status::LINE_SEARCH_FAILED: return "LS_FAILED";
        case Status::FAILED: return "FAILED";
        default: return "UNKNOWN";
      }
    }

    struct Options
    {
      // Block configuration
      integer block_size{ 20 };      ///< Base block size
      double  overlap_ratio{ 0.3 };  ///< Fraction of block overlap (0-1)
      integer min_block_size{ 5 };   ///< Minimum block size
      integer max_block_size{ 15 };  ///< Maximum block size

      // Outer loop control
      integer max_outer_iterations{ 50 };
      integer max_inner_iterations{ 200 };
      Scalar  outer_tolerance{ 1e-6 };

      // Inner L-BFGS options
      integer lbfgs_m{ 10 };
      Scalar  lbfgs_g_tol{ 1e-6 };
      Scalar  lbfgs_f_tol{ 1e-12 };  // Much tighter for inner convergence

      // Adaptive parameters
      bool   adaptive_block_size{ true };
      Scalar progress_threshold{ 1e-6 };            // Relative threshold for block acceptance
      Scalar absolute_progress_threshold{ 1e-12 };  // Absolute threshold for tiny improvements

      // Convergence control
      integer max_stagnation_count{ 3 };
      Scalar  relative_stagnation_tol{ 1e-8 };  // Tighter relative tolerance

      // Verbosity (enhanced)
      integer verbosity_level{ 1 };  // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      bool    use_unicode_borders{ true };
    };

  private:
    Options                      m_options;
    BoxConstraintHandler<Scalar> m_box_handler;
    Vector                       m_x;

    // Optimization results storage
    Status  m_status{ Status::FAILED };
    Scalar  m_final_function_value{ 0 };
    Scalar  m_initial_function_value{ 0 };
    integer m_outer_iterations{ 0 };
    integer m_inner_iterations{ 0 };
    integer m_total_iterations{ 0 };
    integer m_outer_evaluations{ 0 };
    integer m_inner_evaluations{ 0 };
    integer m_total_evaluations{ 0 };
    integer m_blocks_processed{ 0 };
    Scalar  m_max_block_improvement{ 0 };

    integer m_outer_iteration_count{ 0 };

    // Random number generation
    std::random_device m_rd;
    std::mt19937       m_gen;

    // Progress tracking
    Scalar  m_previous_best{ std::numeric_limits<Scalar>::max() };
    integer m_stagnation_count{ 0 };

    string m_indent{ "" };  // For consistent indentation

    // ===========================================================================
    // PRINTING METHODS for Block Coordinate
    // ===========================================================================

    void print_outer_iteration_header(
      integer                 outer_iter,
      integer                 total_cycles,
      vector<integer> const & block_indices,
      integer                 block_size ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
        "\n"
        "    ╔════════════════════════════════════════════════════════════════╗\n"
        "    ║ Outer Iteration {:3d} - Block {:2d}/{:2d}                              ║\n"
        "    ╠════════════════════════════════════════════════════════════════╣\n"
        "    ║ Size: {:<5} Indices: {:<41} ║\n"
        "    ╚════════════════════════════════════════════════════════════════╝"
        "\n",
        outer_iter,
        ( outer_iter % total_cycles ) + 1,
        total_cycles,
        block_size,
        format_index_vector_compact<integer>( block_indices, 7 ) );
    }

    void print_outer_statistics(
      integer total_outer_iters,
      integer total_inner_iters,
      integer total_evals,
      Scalar  best_value,
      bool    converged ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::INFO,
        "\n"
        "  ╔═══════════════════════════════════════════════════════════════╗\n"
        "  ║                      Outer Iteration Summary                  ║\n"
        "  ╠═══════════════════════════════════════════════════════════════╣\n"
        "  ║ Completed:        {:<43} ║\n"
        "  ║ Outer Iterations: {:<43} ║\n"
        "  ║ Inner Iterations: {:<43} ║\n"
        "  ║ Total Evals:      {:<43} ║\n"
        "  ║ Best Value:       {:<43.6e} ║\n"
        "  ║ Status:           {:<43} ║\n"
        "  ╚═══════════════════════════════════════════════════════════════╝\n"
        "\n",
        fmt::format( "{}/{}", m_outer_iteration_count, total_outer_iters ),
        m_outer_iteration_count,
        total_inner_iters,
        total_evals,
        best_value,
        ( converged ? "CONVERGED" : "RUNNING" ) );
    }

    void print_block_coordinate_header( integer n, Scalar f0, Scalar g0 ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
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

    void print_outer_iteration_start( Scalar f, Scalar gnorm ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
        "\n"
        "    ╔════════════════════════════════════════════════════════════════╗\n"
        "    ║ OUTER ITERATION {:3d} START F = {:<12.6e} ‖pg‖ = {:<12.6e} ║\n"
        "    ╚════════════════════════════════════════════════════════════════╝"
        "\n",
        m_outer_iteration_count,
        f,
        gnorm );
    }

    void print_outer_iteration_end( Scalar f, Scalar gnorm, Scalar improvement ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
        "\n"
        "    ╔════════════════════════════════════════════════════════════════╗\n"
        "    ║ OUTER ITERATION {:3d} END                                        ║\n"
        "    ║ F = {:<12.6e}    ‖pg‖ = {:<12.6e}    ΔF = {:<12.6e}   ║\n"
        "    ╚════════════════════════════════════════════════════════════════╝"
        "\n",
        m_outer_iteration_count,
        f,
        gnorm,
        improvement );
    }

    void print_blocks_summary( std::vector<std::vector<integer>> const & blocks ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      fmt::print(
        PrintColors::INFO,
        "    ┌────────────────────────────────────────────────────────────────┐\n"
        "    │                      BLOCKS FOR ITERATION {:3d}                  │\n"
        "    ├────────────────────────────────────────────────────────────────┤\n",
        m_outer_iteration_count );

      for ( integer i = 0; i < blocks.size(); ++i )
      {
        fmt::print(
          PrintColors::INFO,
          "    │ Block {:3d}: Size={:3d}  Indices: {:<32} │\n",
          i + 1,
          blocks[i].size(),
          format_index_vector_compact<integer>( blocks[i], 5 ) );
      }

      fmt::print( PrintColors::INFO, "    └────────────────────────────────────────────────────────────────┘\n" );
    }

    void print_block_result_detail(
      integer block_idx,
      integer total_blocks,
      Scalar  f_before,
      Scalar  f_after,
      integer inner_iters,
      integer inner_evals,
      bool    success ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      Scalar improvement = f_before - f_after;
      auto   color       = success && improvement > 0 ? PrintColors::SUCCESS
                           : success                  ? PrintColors::WARNING
                                                      : PrintColors::ERROR;

      string status_icon = success ? ( improvement > 0 ? "✅" : "⚡" ) : "❌";
      string status_text = success ? ( improvement > 0 ? "IMPROVED" : "NO CHANGE" ) : "FAILED";

      fmt::print(
        color,
        "      ┌────────────────────── BLOCK {:>2d}/{:<2d} RESULT ──────────────────────┐\n"
        "      │ {} Status: {:<10}   Inner Iters: {:<6}   Evals: {:<6}    │\n"
        "      │ F: {:>16.8e} → {:<16.8e}  ΔF: {:<16.8e}   │\n"
        "      └────────────────────────────────────────────────────────────────┘\n",
        block_idx + 1,
        total_blocks,
        status_icon,
        status_text,
        inner_iters,
        inner_evals,
        f_before,
        f_after,
        improvement );
    }

    void print_block_header(
      integer                      block_idx,
      integer                      total_blocks,
      std::vector<integer> const & block_indices,
      integer                      block_size ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      fmt::print(
        PrintColors::INFO,
        "\n"
        "      ┌────────────────────────────────────────────────────────────────┐\n"
        "      │ Block: {:>4d}/{:<4d}    Iteration: {:<4d}                            │\n"
        "      │ Size: {:<4} Indices: {:<42} │\n"
        "      └────────────────────────────────────────────────────────────────┘\n",
        block_idx + 1,
        total_blocks,
        m_outer_iteration_count,
        block_size,
        format_index_vector_compact<integer>( block_indices, 5 ) );
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
    std::vector<std::vector<integer>> generate_consecutive_interval_blocks( integer n_total )
    {
      std::vector<std::vector<integer>> blocks;

      if ( n_total == 0 ) return blocks;

      integer block_size =
        std::clamp( m_options.block_size, m_options.min_block_size, std::min( m_options.max_block_size, n_total ) );

      // Calculate step size based on overlap ratio
      integer step = static_cast<integer>( block_size * ( 1.0 - m_options.overlap_ratio ) );
      step         = std::max<integer>( step, 1 );

      // If step >= block_size, we get no overlap - prevent this
      if ( step >= block_size )
      {
        step = block_size / 2;
        if ( step < 1 ) step = 1;
      }

      // Generate blocks with overlap
      integer start = 0;
      while ( start < n_total )
      {
        integer              end = std::min( start + block_size - 1, n_total - 1 );
        std::vector<integer> block;
        for ( integer i = start; i <= end; ++i ) { block.push_back( i ); }
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
    void adapt_block_size( Scalar improvement, integer n_total )
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
      m_options.block_size =
        std::clamp( m_options.block_size, m_options.min_block_size, std::min( m_options.max_block_size, n_total ) );

      m_previous_best = improvement;

      if ( m_options.verbosity_level >= 2 )
      {
        fmt::print( "    [BlockLBFGS] Adjusted block_size to {}\n", m_options.block_size );
      }
    }

    Vector extract_subvector( Vector const & full, std::vector<integer> const & indices )
    {
      Vector sub( indices.size() );
      for ( integer i = 0; i < indices.size(); ++i ) { sub( i ) = full( indices[i] ); }
      return sub;
    }

    void update_full_vector( Vector & full, Vector const & sub, std::vector<integer> const & indices )
    {
      for ( integer i = 0; i < indices.size(); ++i ) { full( indices[i] ) = sub( i ); }
    }

    void extract_subbounds( std::vector<integer> const & indices, Vector & sub_lower, Vector & sub_upper )
    {
      if ( !m_box_handler.is_active() )
      {
        sub_lower.resize( indices.size() );
        sub_upper.resize( indices.size() );
        sub_lower.setConstant( -std::numeric_limits<Scalar>::infinity() );
        sub_upper.setConstant( std::numeric_limits<Scalar>::infinity() );
        return;
      }

      // Extract global bounds
      Vector const & global_lower = m_box_handler.get_lower();
      Vector const & global_upper = m_box_handler.get_upper();

      sub_lower.resize( indices.size() );
      sub_upper.resize( indices.size() );

      for ( integer i = 0; i < indices.size(); ++i )
      {
        sub_lower( i ) = global_lower( indices[i] );
        sub_upper( i ) = global_upper( indices[i] );
      }
    }

    void project_point( Vector & x ) const
    {
      if ( m_box_handler.is_active() ) { m_box_handler.project( x ); }
    }

    Scalar projected_gradient_norm( Vector const & x, Vector const & g ) const
    {
      if ( !m_box_handler.is_active() ) return g.template lpNorm<Eigen::Infinity>();
      return m_box_handler.projected_gradient_norm_inf( x, g );
    }

  public:
    explicit LBFGS_BlockCoordinate( Options const & opts = Options() ) : m_options( opts ), m_gen( m_rd() ) {}

    void set_bounds( Vector const & lower, Vector const & upper )
    {
      UTILS_ASSERT(
        lower.size() == upper.size(),
        "BlockLBFGS::set_bounds: lower and upper bounds must have "
        "same dimension" );

      UTILS_ASSERT(
        ( lower.array() <= upper.array() ).all(),
        "BlockLBFGS::set_bounds: lower bounds must be <= upper bounds "
        "for all coordinates" );

      m_box_handler.set_bounds( lower, upper );
    }

    void clear_bounds() { m_box_handler.clear(); }

    Vector const & solution() const { return m_x; }

    // Accessor methods for optimization results
    Status  status() const { return m_status; }
    Scalar  final_function_value() const { return m_final_function_value; }
    Scalar  initial_function_value() const { return m_initial_function_value; }
    integer outer_iterations() const { return m_outer_iterations; }
    integer inner_iterations() const { return m_inner_iterations; }
    integer total_iterations() const { return m_total_iterations; }
    integer outer_evaluations() const { return m_outer_evaluations; }
    integer inner_evaluations() const { return m_inner_evaluations; }
    integer total_evaluations() const { return m_total_evaluations; }
    integer blocks_processed() const { return m_blocks_processed; }
    Scalar  max_block_improvement() const { return m_max_block_improvement; }

    /**
     * @brief Enhanced optimization routine with consecutive interval blocks
     */
    template <typename Linesearch> void minimize(
      Vector const &     x0,
      Callback const &   global_callback,
      Linesearch const & linesearch = MoreThuenteLineSearch<Scalar>() )
    {
      integer n = x0.size();

      // Reset results
      m_status                 = Status::FAILED;
      m_final_function_value   = 0;
      m_initial_function_value = 0;
      m_outer_iterations       = 0;
      m_inner_iterations       = 0;
      m_total_iterations       = 0;
      m_outer_evaluations      = 0;
      m_inner_evaluations      = 0;
      m_total_evaluations      = 0;
      m_blocks_processed       = 0;
      m_max_block_improvement  = 0;

      // If dimension <= block_size, use standard LBFGS
      if ( n <= m_options.block_size )
      {
        if ( m_options.verbosity_level >= 1 )
        {
          fmt::print(
            "    [BlockLBFGS] Problem dimension {} <= block_size {}, "
            "using standard LBFGS\n",
            n,
            m_options.block_size );
        }

        // Use standard LBFGS instead of block coordinate
        typename LBFGS_minimizer<Scalar>::Options opts;
        opts.max_iter        = m_options.max_outer_iterations * m_options.max_inner_iterations;
        opts.g_tol           = m_options.outer_tolerance;
        opts.m               = m_options.lbfgs_m;
        opts.verbosity_level = m_options.verbosity_level;

        LBFGS_minimizer<Scalar> minimizer( opts );
        if ( m_box_handler.is_active() )
        {
          minimizer.set_bounds( m_box_handler.get_lower(), m_box_handler.get_upper() );
        }

        minimizer.minimize( x0, global_callback, linesearch );

        m_x = minimizer.solution();

        // Store results from standard LBFGS
        m_final_function_value   = minimizer.final_function_value();
        m_initial_function_value = minimizer.initial_function_value();
        m_status                 = ( minimizer.status() == LBFGS_minimizer<Scalar>::Status::CONVERGED ||
                     minimizer.status() == LBFGS_minimizer<Scalar>::Status::GRADIENT_TOO_SMALL )
                                     ? Status::CONVERGED
                                     : Status::FAILED;
        m_outer_iterations       = 1;
        m_inner_iterations       = minimizer.total_iterations();
        m_total_iterations       = minimizer.total_iterations();
        m_total_evaluations      = minimizer.total_evaluations();
        m_outer_evaluations      = minimizer.total_evaluations();
        m_inner_evaluations      = minimizer.total_evaluations();
        m_blocks_processed       = 1;
        m_max_block_improvement  = m_initial_function_value - m_final_function_value;

        return;
      }

      // Initialize solution
      m_x.resize( x0.size() );
      m_x = x0;
      project_point( m_x );

      // Initial evaluation
      Vector g( n );
      Scalar f                 = global_callback( m_x, &g );
      m_initial_function_value = f;

      integer total_evals           = 1;
      integer total_inner_iters     = 0;
      integer blocks_processed      = 0;
      Scalar  max_block_improvement = 0;

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

        // PRINT OUTER ITERATION START
        print_outer_iteration_start( f, projected_gradient_norm( m_x, g ) );

        // Generate consecutive interval blocks for this iteration
        auto all_blocks = generate_consecutive_interval_blocks( n );
        print_blocks_summary( all_blocks );

        // Track if any block made progress
        bool   made_progress_in_cycle    = false;
        Scalar best_improvement_in_cycle = 0;

        // Process each block
        for ( integer block_idx = 0; block_idx < all_blocks.size(); ++block_idx )
        {
          auto block_indices = all_blocks[block_idx];

          if ( block_indices.empty() ) continue;

          // PRINT BLOCK HEADER
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
                block_gnorm,
                m_options.lbfgs_g_tol * 0.1 );
            }
            continue;
          }

          // Create block callback
          auto block_callback = [&]( Vector const & x_sub, Vector * g_sub ) -> Scalar
          {
            Vector x_full = m_x;
            update_full_vector( x_full, x_sub, block_indices );
            project_point( x_full );

            Vector g_full( n );
            Scalar f_val = global_callback( x_full, &g_full );
            total_evals++;

            if ( g_sub ) { *g_sub = extract_subvector( g_full, block_indices ); }

            return f_val;
          };

          // Scale tolerances based on current value
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
          inner_opts.verbosity_level = max( m_options.verbosity_level, static_cast<integer>( 1 ) );

          // Improved inner L-BFGS options
          inner_opts.step_max        = 1.0;
          inner_opts.very_small_step = 1e-10;
          inner_opts.iter_reset      = std::min( m_options.lbfgs_m, integer( 10 ) );
          inner_opts.x_tol           = 1e-12;
          inner_opts.sty_min_factor  = 1e-14;

          // Debug print before inner optimization
          if ( m_options.verbosity_level >= 3 )
          {
            fmt::print(
              "    [BlockLBFGS] Starting inner L-BFGS on block [{}] "
              "(size: {}, ‖g‖∞: {:.2e})\n",
              format_index_vector_compact<integer>( block_indices, 10 ),
              block_indices.size(),
              g_block.template lpNorm<Eigen::Infinity>() );
          }

          LBFGS_minimizer<Scalar> inner_minimizer( inner_opts );
          inner_minimizer.set_indent( "      " );  // 6 spaces indentation

          // Set bounds for block if needed
          if ( m_box_handler.is_active() )
          {
            Vector block_lower, block_upper;
            extract_subbounds( block_indices, block_lower, block_upper );
            inner_minimizer.set_bounds( block_lower, block_upper );
          }

          // Save current state for comparison
          Vector x_before_block = m_x;
          Scalar f_before_block = f;

          // Optimize block
          inner_minimizer.minimize( x_block, block_callback, linesearch );

          total_inner_iters += inner_minimizer.total_iterations();
          blocks_processed++;

          // Update solution
          update_full_vector( m_x, inner_minimizer.solution(), block_indices );
          project_point( m_x );

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

          if ( block_made_progress )
          {
            made_progress_in_cycle    = true;
            best_improvement_in_cycle = std::max( best_improvement_in_cycle, block_improvement );
          }

          // PRINT DETAILED BLOCK RESULT
          print_block_result_detail(
            block_idx,
            all_blocks.size(),
            f_before_block,
            f,
            inner_minimizer.total_iterations(),
            inner_minimizer.total_evaluations(),
            block_made_progress );
        }

        // PRINT OUTER ITERATION END
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
        m_options.block_size =
          std::clamp( m_options.block_size, m_options.min_block_size, std::min( m_options.max_block_size, n ) );

        if ( m_options.verbosity_level >= 2 && m_options.adaptive_block_size )
        {
          fmt::print( "    [BlockLBFGS] Adjusted block_size to {}\n", m_options.block_size );
        }

        // Convergence checking
        Scalar gnorm = projected_gradient_norm( m_x, g );

        // Calculate current coverage for display
        // Print outer statistics periodically
        if (
          m_options.verbosity_level >= 1 &&
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
              gnorm,
              m_options.outer_tolerance );
          }
          m_status  = Status::GRADIENT_TOO_SMALL;
          converged = true;
        }

        // Check function value stagnation with improved criteria
        Scalar cycle_improvement   = f_start_cycle - f;
        Scalar relative_stagnation = cycle_improvement / ( std::abs( f_start_cycle ) + Scalar( 1e-16 ) );

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
                m_stagnation_count,
                m_options.progress_threshold );
            }
            m_status  = Status::CONVERGED;
            converged = true;
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
      m_final_function_value  = f;
      m_outer_iterations      = m_outer_iteration_count;
      m_inner_iterations      = total_inner_iters;
      m_total_iterations      = m_outer_iterations + total_inner_iters;
      m_outer_evaluations     = total_evals;
      m_inner_evaluations     = total_evals;
      m_total_evaluations     = total_evals;
      m_blocks_processed      = blocks_processed;
      m_max_block_improvement = max_block_improvement;

      if ( m_status == Status::FAILED && m_outer_iterations >= m_options.max_outer_iterations )
      {
        m_status = Status::MAX_OUTER_ITERATIONS;
      }

      // Print final outer statistics
      if ( m_options.verbosity_level >= 1 )
      {
        print_outer_statistics( m_options.max_outer_iterations, total_inner_iters, total_evals, f, true );

        fmt::print(
          "    [BlockLBFGS] Finished: status={}, f_final={:.6e}, "
          "total_evals={}\n",
          static_cast<int>( m_status ),
          m_final_function_value,
          total_evals );
      }
    }
  };

}  // namespace Utils

#endif

//
// eof: Utils_minimize_LBFGS_BlockCoordinate.hh
//
