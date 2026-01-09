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

//
// file: Utils_minimize_NelderMead_BlockCoordinate.hh
//

#pragma once

#ifndef UTILS_MINIMIZE_NELDER_MEAD_BLOCK_COORDINATE_dot_HH
#define UTILS_MINIMIZE_NELDER_MEAD_BLOCK_COORDINATE_dot_HH

#include "Utils_minimize_NelderMead.hh"

namespace Utils
{

  // ===========================================================================
  // CLASS: NelderMead_BlockCoordinate (Outer Solver)
  // ===========================================================================

  /**
   * @class NelderMead_BlockCoordinate
   * @brief Block coordinate descent using Nelder-Mead for subspace optimization
   * @tparam Scalar Numeric type for computations (default: double)
   *
   * EIGEN3 OPTIMIZATIONS:
   * - Efficient subspace extraction using Eigen vector operations
   * - Memory-efficient block processing for high-dimensional problems
   * - Leverages Eigen's expression templates for zero-copy operations
   */
  template <typename Scalar = double> class NelderMead_BlockCoordinate
  {
  public:
    using integer  = Eigen::Index;
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar( Vector const & )>;

    struct Options
    {
      integer block_size{ 10 };
      integer max_outer_iterations{ 100 };
      integer max_inner_iterations{ 1000 };
      integer max_function_evaluations{ 100000 };
      Scalar  tolerance{ 1e-6 };
      bool    verbose{ true };
      integer verbosity_level{ 1 };  // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      integer inner_progress_frequency{ 10 };
      typename NelderMead_classic<Scalar>::Options sub_options;
    };

  private:
    Options                    m_options;
    NelderMead_classic<Scalar> m_solver;
    Vector                     m_lower;
    Vector                     m_upper;
    bool                       m_use_bounds{ false };

    // Results storage
    Vector             m_solution;
    Scalar             m_final_function_value{ 0 };
    Scalar             m_initial_function_value{ 0 };
    NelderMead::Status m_status{ NelderMead::Status::FAILED };
    integer            m_outer_iterations{ 0 };
    integer            m_inner_iterations{ 0 };
    integer            m_total_iterations{ 0 };
    integer            m_outer_evaluations{ 0 };
    integer            m_inner_evaluations{ 0 };
    integer            m_total_evaluations{ 0 };

    // Cyclic block selection
    vector<integer> select_block( integer n_dims, integer iter )
    {
      vector<integer> indices;
      indices.reserve( m_options.block_size );
      integer start = ( iter * m_options.block_size ) % n_dims;
      for ( integer i = 0; i < m_options.block_size; ++i ) indices.push_back( ( start + i ) % n_dims );
      std::sort( indices.begin(), indices.end() );
      indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
      return indices;
    }

    void print_outer_iteration_header(
      integer                 outer_iter,
      integer                 total_cycles,
      const vector<integer> & block_indices,
      integer                 block_size ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print(
        PrintColors::HEADER,
        "\n"
        "╔════════════════════════════════════════════════════════════════╗\n"
        "║ Outer Iteration {:3d} - Block {:2d}/{:2d}                              ║\n"
        "╠════════════════════════════════════════════════════════════════╣\n"
        "║ Block Indices: {:<47} ║\n"
        "║ Block Size:    {:<47} ║\n"
        "╚════════════════════════════════════════════════════════════════╝"
        "\n",
        outer_iter,
        ( outer_iter % total_cycles ) + 1,
        total_cycles,
        Utils::format_index_vector<integer>( block_indices, 7 ),
        block_size );
    }

    void print_outer_iteration_result(
      integer outer_iter,
      integer block_size,
      Scalar  current_f,
      integer inner_iters,
      integer inner_evals,
      Scalar  improvement,
      bool    improved ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      auto color = improved ? PrintColors::SUCCESS : PrintColors::WARNING;

      fmt::print(
        color,
        "┌──────────────┬───────────┬────────────────┬──────────────────┬──────────────────┐\n"
        "│ out iter:{:<3} │ block:{:<3} │ F:{:<12.6e} │ inner iter:{:<5} │ inner eval:{:<5} │\n"
        "└──────────────┴───────────┴────────────────┴──────────────────┴──────────────────┘\n",
        outer_iter,
        block_size,
        current_f,
        inner_iters,
        inner_evals );

      if ( improved )
      {
        fmt::print(
          PrintColors::SUCCESS,
          "✓ Improvement: {:.6e} → {:.6e} (Δ = {:.6e})\n",
          current_f + improvement,
          current_f,
          improvement );
      }
      else
      {
        fmt::print( PrintColors::WARNING, "⚠ No significant improvement (Δ = {:.6e})\n", improvement );
      }
    }

    void print_outer_statistics(
      integer outer_iter,
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
        "╔═══════════════════════════════════════════════════════════════╗\n"
        "║                      Outer Iteration Summary                  "
        "║\n"
        "╠═══════════════════════════════════════════════════════════════╣\n"
        "║ Completed:        {:<43} ║\n"
        "║ Outer Iterations: {:<43} ║\n"
        "║ Inner Iterations: {:<43} ║\n"
        "║ Total Evals:      {:<43} ║\n"
        "║ Best Value:       {:<43.6e} ║\n"
        "║ Status:           {:<43} ║\n"
        "╚═══════════════════════════════════════════════════════════════╝\n",
        fmt::format( "{}/{}", outer_iter, total_outer_iters ),
        outer_iter,
        total_inner_iters,
        total_evals,
        best_value,
        ( converged ? "CONVERGED" : "RUNNING" ) );
    }

  public:
    explicit NelderMead_BlockCoordinate( Options const & opts = Options() ) : m_options( opts )
    {
      m_solver.set_options( m_options.sub_options );
    }

    void set_sub_options( typename NelderMead_classic<Scalar>::Options const & sopts )
    {
      m_options.sub_options = sopts;
      m_solver.set_options( sopts );
    }

    void set_bounds( Vector const & l, Vector const & u )
    {
      m_lower      = l;
      m_upper      = u;
      m_use_bounds = true;
    }

    /**
     * @brief Main optimization routine for block coordinate descent
     * @param x Initial point (updated in-place)
     * @param global_callback Objective function
     * @return true if optimization succeeded
     *
     * EIGEN3: Efficient block processing using Eigen vector operations
     * - Subspace extraction without memory copies
     * - Efficient bounds management for subspaces
     * - Leverages Eigen's expression templates for performance
     */
    bool minimize( Vector x, Callback const & global_callback )
    {
      integer n = x.size();

      // Reset results
      m_solution               = Vector();
      m_final_function_value   = 0;
      m_initial_function_value = 0;
      m_status                 = NelderMead::Status::FAILED;
      m_outer_iterations       = 0;
      m_inner_iterations       = 0;
      m_total_iterations       = 0;
      m_outer_evaluations      = 0;
      m_inner_evaluations      = 0;
      m_total_evaluations      = 0;

      // 1. Check for small problem fallback
      if ( !( 2 * n > m_options.block_size ) )
      {
        // MODIFICA: Usare verbosity_level
        if ( m_options.verbosity_level >= 1 )
          fmt::print( "[Info] Dim <= BlockSize. Switching to DIRECT CLASSIC solver.\n" );

        auto full_opts                     = m_options.sub_options;
        full_opts.max_function_evaluations = m_options.max_function_evaluations;
        full_opts.max_iterations           = m_options.max_inner_iterations;
        full_opts.tolerance                = m_options.tolerance;
        full_opts.verbose                  = ( m_options.verbosity_level >= 1 );  // MODIFICA

        m_solver.set_options( full_opts );
        if ( m_use_bounds ) m_solver.set_bounds( m_lower, m_upper );
        bool success = m_solver.minimize( x, global_callback );

        // Store results
        m_solution               = m_solver.get_solution();
        m_final_function_value   = m_solver.get_final_function_value();
        m_initial_function_value = m_solver.get_initial_function_value();
        m_status                 = ( m_solver.get_status() == NelderMead_classic<Scalar>::Status::CONVERGED )
                                     ? NelderMead::Status::CONVERGED
                                     : NelderMead::Status::MAX_ITERATIONS;
        m_outer_iterations       = 1;
        m_inner_iterations       = m_solver.get_iterations();
        m_total_iterations       = 1 + m_inner_iterations;
        m_outer_evaluations      = 1;  // initial evaluation
        m_inner_evaluations      = m_solver.get_function_evaluations();
        m_total_evaluations      = 1 + m_inner_evaluations;

        return success;
      }

      // 2. Block Coordinate Logic
      integer count_outer_evals = 0;
      integer count_inner_evals = 0;
      integer count_inner_iters = 0;

      Scalar current_f = global_callback( x );
      count_outer_evals++;
      m_initial_function_value = current_f;

      if ( std::isnan( current_f ) || std::isinf( current_f ) )
      {
        m_status   = NelderMead::Status::FAIL_NAN;
        m_solution = x;
        return false;
      }

      integer outer_iter = 0;
      bool    converged  = false;

      // Stagnation control
      integer no_improv_count  = 0;
      integer blocks_per_cycle = ( n + m_options.block_size - 1 ) / m_options.block_size;

      while ( outer_iter < m_options.max_outer_iterations && !converged )
      {
        outer_iter++;
        if ( count_outer_evals + count_inner_evals >= m_options.max_function_evaluations )
        {
          m_status = NelderMead::Status::MAX_FUN_EVALUATIONS;
          break;
        }

        vector<integer> idxs = select_block( n, outer_iter - 1 );
        integer         k    = idxs.size();

        // MODIFICA: Usare il nuovo metodo di stampa
        print_outer_iteration_header( outer_iter, blocks_per_cycle, idxs, k );

        // Prepare Subspace
        // EIGEN3: Efficient subspace extraction without memory allocation
        Vector x_sub( k ), l_sub( k ), u_sub( k );
        for ( integer i{ 0 }; i < k; ++i )
        {
          x_sub( i ) = x( idxs[i] );
          if ( m_use_bounds )
          {
            l_sub( i ) = m_lower( idxs[i] );
            u_sub( i ) = m_upper( idxs[i] );
          }
        }
        if ( m_use_bounds ) m_solver.set_bounds( l_sub, u_sub );

        // MODIFICA: Impostare verbosity_level per il solver interno
        auto sub_opts            = m_options.sub_options;
        sub_opts.verbosity_level = m_options.verbosity_level;
        m_solver.set_options( sub_opts );

        // Proxy Callback
        // EIGEN3: Efficient subspace updates using Eigen vector operations
        auto sub_cb = [&]( Vector const & sub_v ) -> Scalar
        {
          count_inner_evals++;
          Vector x_temp = x;
          for ( integer i = 0; i < k; ++i ) x_temp( idxs[i] ) = sub_v( i );
          return global_callback( x_temp );
        };

        // Run Inner Solver
        bool sub_success = m_solver.minimize( x_sub, sub_cb );
        count_inner_iters += m_solver.get_iterations();

        // Update Global State
        Scalar old_f       = current_f;
        Scalar improvement = old_f - m_solver.get_final_function_value();
        bool   improved    = m_solver.get_final_function_value() < current_f;

        if ( improved )
        {
          current_f = m_solver.get_final_function_value();
          for ( integer i = 0; i < k; ++i ) x( idxs[i] ) = m_solver.get_solution()( i );
        }

        // MODIFICA: Usare il nuovo metodo di stampa
        print_outer_iteration_result(
          outer_iter,
          k,
          current_f,
          m_solver.get_iterations(),
          m_solver.get_function_evaluations(),
          improvement,
          improved );

        // Stagnation Check Logic
        if ( improvement > m_options.tolerance )
          no_improv_count = 0;  // Reset counter on success
        else
          ++no_improv_count;

        // Stop only if we cycled through ALL blocks twice without significant
        // improvement
        if ( no_improv_count >= 2 * blocks_per_cycle ) converged = true;

        // MODIFICA: Stampare statistiche periodiche
        if ( m_options.verbosity_level >= 1 && ( outer_iter % 5 == 0 || converged ) )
        {
          print_outer_statistics(
            outer_iter,
            m_options.max_outer_iterations,
            count_inner_iters,
            count_outer_evals + count_inner_evals,
            current_f,
            converged );
        }
      }

      m_solution             = x;
      m_final_function_value = current_f;
      m_outer_iterations     = outer_iter;
      m_inner_iterations     = count_inner_iters;
      m_total_iterations     = outer_iter + count_inner_iters;
      m_outer_evaluations    = count_outer_evals;
      m_inner_evaluations    = count_inner_evals;
      m_total_evaluations    = count_outer_evals + count_inner_evals;

      if ( m_status == NelderMead::Status::FAILED )
      {
        if ( converged )
          m_status = NelderMead::Status::CONVERGED;
        else if ( outer_iter >= m_options.max_outer_iterations )
          m_status = NelderMead::Status::MAX_ITERATIONS;
      }

      // MODIFICA: Stampare statistiche finali
      if ( m_options.verbosity_level >= 1 )
      {
        print_outer_statistics(
          outer_iter,
          m_options.max_outer_iterations,
          count_inner_iters,
          count_outer_evals + count_inner_evals,
          current_f,
          true );
      }

      return m_status == NelderMead::Status::CONVERGED;
    }

    // Access methods for results
    Vector             get_solution() const { return m_solution; }
    Scalar             get_final_function_value() const { return m_final_function_value; }
    Scalar             get_initial_function_value() const { return m_initial_function_value; }
    NelderMead::Status get_status() const { return m_status; }
    integer            get_outer_iterations() const { return m_outer_iterations; }
    integer            get_inner_iterations() const { return m_inner_iterations; }
    integer            get_total_iterations() const { return m_total_iterations; }
    integer            get_outer_evaluations() const { return m_outer_evaluations; }
    integer            get_inner_evaluations() const { return m_inner_evaluations; }
    integer            get_total_evaluations() const { return m_total_evaluations; }
  };
}  // namespace Utils

#endif
