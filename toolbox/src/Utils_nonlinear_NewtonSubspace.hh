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

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
#include <string>
#include <iomanip>
#include <numeric>
#include <unordered_map>
#include <deque>

// Assicurati che questi file siano nel path
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"
#include "Utils_nonlinear_system.hh"
#include "Utils_pseudoinverse.hh"

namespace Utils
{

  /*\
   |   _____       _
   |  / ____|     | |
   | | (___  _   _| |__  ___ _ __   __ _  ___ ___
   |  \___ \| | | | '_ \/ __| '_ \ / _` |/ __/ _ \
   |  ____) | |_| | |_) \__ \ |_) | (_| | (_|  __/
   | |_____/ \__,_|_.__/|___/ .__/ \__,_|\___\___|
   |                        | |
   |                        |_|
  \*/

  class SubspaceInnerSolver
  {
  public:
    using Scalar       = double;
    using Vector       = Eigen::VectorXd;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

  private:
    int    m_max_inner_iter{ 10 };
    Scalar m_inner_tol{ 1e-4 };
    Scalar m_inner_lambda{ 1e-6 };
    Scalar m_lambda_min{ 1e-10 };
    Scalar m_lambda_max{ 1e10 };
    bool   m_use_linesearch{ true };
    int    m_verbose{ 0 };

    // Lazy Jacobian parameters
    bool m_lazy_jacobian{ true };
    int  m_lazy_update_freq{ 3 };

    // Counters
    int m_total_inner_iter{ 0 };
    int m_inner_feval{ 0 };
    int m_inner_jeval{ 0 };

    // Adaptive damping parameters
    bool   m_adaptive_damping{ true };
    Scalar m_damping_increase_factor{ 2.0 };
    Scalar m_damping_decrease_factor{ 0.5 };
    int    m_max_consecutive_rejects{ 3 };
    int    m_consecutive_rejects{ 0 };

  public:
    SubspaceInnerSolver() = default;

    void
    set_max_iterations( int it )
    {
      m_max_inner_iter = it;
    }
    void
    set_tolerance( Scalar t )
    {
      m_inner_tol = t;
    }
    void
    set_damping( Scalar l )
    {
      m_inner_lambda = l;
    }
    void
    enable_line_search( bool v )
    {
      m_use_linesearch = v;
    }
    void
    set_verbose_level( int v )
    {
      m_verbose = v;
    }
    void
    enable_lazy_jacobian( bool v )
    {
      m_lazy_jacobian = v;
    }
    void
    enable_adaptive_damping( bool v )
    {
      m_adaptive_damping = v;
    }

    // Getters for stats
    int
    get_total_inner_iterations() const
    {
      return m_total_inner_iter;
    }
    int
    get_inner_function_evals() const
    {
      return m_inner_feval;
    }
    int
    get_inner_jacobian_evals() const
    {
      return m_inner_jeval;
    }

    void
    reset_counters()
    {
      m_total_inner_iter    = 0;
      m_inner_feval         = 0;
      m_inner_jeval         = 0;
      m_consecutive_rejects = 0;
    }

    // Risolve il problema ridotto mantenendo fisse le variabili non in active_indices
    bool
    solve_in_subspace( NonlinearSystem const &  sys,
                       Vector &                 x_full,
                       std::vector<int> const & active_indices,
                       Vector &                 f,
                       int                      outer_iter_idx )
    {
      int n_full   = x_full.size();
      int n_active = static_cast<int>( active_indices.size() );
      if ( n_active == 0 ) return true;

      // Mapping: active_subset -> full_vector
      Vector x_active( n_active );
      for ( int i = 0; i < n_active; ++i ) { x_active[i] = x_full[active_indices[i]]; }

      // Mappa colonne sparse per estrazione rapida
      std::vector<int> col_map( n_full, -1 );
      for ( int i = 0; i < n_active; ++i ) { col_map[active_indices[i]] = i; }

      Scalar current_norm = f.norm();
      Scalar initial_norm = current_norm;
      Scalar lambda       = m_inner_lambda;
      Scalar nu           = 2.0;

      SparseMatrix J_sub;  // Jacobiano persistente per lazy evaluation
      m_consecutive_rejects = 0;

      int inner_iter = 0;

      for ( ; inner_iter < m_max_inner_iter; ++inner_iter )
      {
        // Modified Newton: Aggiorna J solo alla prima iter o ogni m_lazy_update_freq
        bool force_update = ( inner_iter == 0 ) || ( !m_lazy_jacobian ) || ( inner_iter % m_lazy_update_freq == 0 );

        if ( force_update )
        {
          SparseMatrix J_full;
          try
          {
            sys.jacobian( x_full, J_full );
            ++m_inner_jeval;
          }
          catch ( std::exception const & e )
          {
            if ( m_verbose > 0 ) fmt::print( fg( fmt::color::red ), "    [Jac Error] {}\n", e.what() );
            return false;
          }

          // Estrazione efficiente delle colonne attive
          J_sub.resize( sys.num_equations(), n_active );
          std::vector<Eigen::Triplet<Scalar>> triplets;
          // Stima non zeri per reserve (euristica)
          triplets.reserve( J_full.nonZeros() * n_active / std::max( 1, n_full ) );

          for ( int k = 0; k < J_full.outerSize(); ++k )
          {
            for ( typename SparseMatrix::InnerIterator it( J_full, k ); it; ++it )
            {
              int col_idx = col_map[it.col()];
              if ( col_idx >= 0 ) { triplets.emplace_back( it.row(), col_idx, it.value() ); }
            }
          }
          J_sub.setFromTriplets( triplets.begin(), triplets.end() );
        }

        Vector dx_active( n_active );
        bool   direction_ok = compute_newton_direction_with_tikhonov( J_sub, f, lambda, dx_active );

        if ( !direction_ok )
        {
          m_consecutive_rejects++;
          if ( m_consecutive_rejects > m_max_consecutive_rejects )
          {
            if ( m_verbose > 0 )
            {
              fmt::print( fg( fmt::color::red ), "    Too many consecutive direction failures\n" );
            }
            break;
          }

          // Se la direzione fallisce con J vecchio, forziamo l'aggiornamento
          if ( !force_update )
          {
            inner_iter--;             // Ripeti iterazione con J fresco
            m_lazy_jacobian = false;  // Disabilita lazy temporaneamente
            continue;
          }
          break;
        }
        else
        {
          m_consecutive_rejects = 0;  // Reset su successo
        }

        Vector x_active_new = x_active;
        Scalar alpha        = 1.0;
        bool   step_success = false;

        if ( m_use_linesearch )
        {
          step_success = subspace_line_search( sys, x_full, x_active_new, active_indices, f, dx_active, alpha );
        }
        else
        {
          x_active_new = x_active + dx_active;
          step_success = update_and_evaluate( sys, x_full, x_active_new, active_indices, f );
        }

        Scalar new_norm  = f.norm();
        Scalar reduction = 0.0;
        if ( initial_norm > 1e-16 ) reduction = 100.0 * ( 1.0 - new_norm / initial_norm );

        // --- STAMPA MIGLIORATA ---
        if ( m_verbose > 1 )
        {
          bool        is_good    = step_success && ( new_norm < current_norm );
          auto        color_code = is_good ? fg( fmt::color::green ) : fg( fmt::color::red );
          std::string status_msg = is_good ? "converging" : "diverging/rejected";

          fmt::print( "  [outer:{:3d} | inner:{:3d}] |f|={:12.3e} red={:6.2f}% ", outer_iter_idx, inner_iter, new_norm,
                      reduction );

          fmt::print( color_code, "--> {} (α={:.1e}, λ={:.1e})\n", status_msg, alpha, lambda );
        }
        // ---------------------

        if ( step_success && new_norm < current_norm )
        {
          x_active     = x_active_new;
          current_norm = new_norm;

          // Adaptive damping: reduce lambda when successful
          if ( m_adaptive_damping ) { lambda = std::max( lambda * m_damping_decrease_factor, m_lambda_min ); }
          else
          {
            lambda = std::max( lambda / 3.0, m_lambda_min );
          }
          nu = 2.0;

          if ( new_norm < m_inner_tol * initial_norm || new_norm < 1e-12 ) { break; }
        }
        else
        {
          // Increase damping when step fails
          if ( m_adaptive_damping ) { lambda = std::min( lambda * m_damping_increase_factor, m_lambda_max ); }
          else
          {
            lambda = std::min( lambda * nu, m_lambda_max );
          }
          nu = 2.0 * nu;
          if ( lambda >= m_lambda_max ) break;
        }
      }

      m_total_inner_iter += ( inner_iter + 1 );

      // Assicuriamoci che x_full sia sincronizzato
      for ( int i = 0; i < n_active; ++i ) { x_full[active_indices[i]] = x_active[i]; }

      // Valutazione finale per sicurezza
      try
      {
        sys.evaluate( x_full, f );
        ++m_inner_feval;
      }
      catch ( ... )
      {
        return false;
      }

      return true;
    }

  private:
    // Usa SP_TikhonovSolver2 per calcolare la direzione di Newton con regolarizzazione diagonale
    bool
    compute_newton_direction_with_tikhonov( SparseMatrix const & J_sub, Vector const & f, Scalar lambda, Vector & dx )
    {
      int n = J_sub.cols();

      // Calcola regolarizzazione diagonale basata sulle norme delle colonne
      Vector D = Vector::Ones( n );
      for ( int i = 0; i < n; ++i )
      {
        // Norma al quadrato della colonna i-esima
        Scalar col_norm_sq = 0.0;
        for ( typename SparseMatrix::InnerIterator it( J_sub, i ); it; ++it )
        {
          col_norm_sq += it.value() * it.value();
        }
        // Regolarizzazione proporzionale alla norma della colonna
        D[i] = lambda * lambda * ( col_norm_sq + 1e-8 );
      }

      try
      {
        // Usa SP_TikhonovSolver2 con regolarizzazione diagonale D
        // Questo risolve: min ||J_sub * dx + f||^2 + ||D^{1/2} * dx||^2
        SP_TikhonovSolver2<Scalar> solver( J_sub, D );
        dx = solver.solve( -f );

        return dx.norm() > 1e-14 && dx.allFinite();
      }
      catch ( std::exception const & e )
      {
        // Fallback al metodo originale se SP_TikhonovSolver2 fallisce
        if ( m_verbose > 1 )
        {
          fmt::print( fg( fmt::color::yellow ), "    SP_TikhonovSolver2 failed: {}, using fallback\n", e.what() );
        }

        // Fallback: risolve direttamente (J^T J + diag(D)) dx = -J^T f
        Vector       rhs = -J_sub.transpose() * f;
        SparseMatrix JtJ = J_sub.transpose() * J_sub;

        for ( int i = 0; i < n; ++i ) { JtJ.coeffRef( i, i ) += D[i]; }

        Eigen::SimplicialLDLT<SparseMatrix> solver;
        solver.compute( JtJ );
        if ( solver.info() == Eigen::Success ) { dx = solver.solve( rhs ); }
        else
        {
          // Fallback a decomposizione densa
          Eigen::MatrixXd JtJ_dense = Eigen::MatrixXd( JtJ );
          dx                        = JtJ_dense.ldlt().solve( rhs );
        }

        return dx.norm() > 1e-14 && dx.allFinite();
      }
    }

    bool
    subspace_line_search( NonlinearSystem const &  sys,
                          Vector const &           x_full,
                          Vector &                 x_active,
                          std::vector<int> const & active_indices,
                          Vector &                 f,
                          Vector const &           dx_active,
                          Scalar &                 alpha_out )
    {
      Scalar alpha    = 1.0;
      Scalar c        = 1e-4;
      Scalar tau      = 0.5;
      Scalar fnorm_sq = f.squaredNorm();

      Vector trial_active( x_active.size() );
      Vector x_trial( x_full.size() );
      Vector f_trial( f.size() );

      int max_ls_iter = 10;

      for ( int ls_iter = 0; ls_iter < max_ls_iter; ++ls_iter )
      {
        trial_active = x_active + alpha * dx_active;

        x_trial = x_full;
        for ( size_t i = 0; i < active_indices.size(); ++i ) { x_trial[active_indices[i]] = trial_active[i]; }

        try
        {
          sys.evaluate( x_trial, f_trial );
          ++m_inner_feval;
        }
        catch ( ... )
        {
          alpha *= tau;
          continue;
        }

        Scalar new_fnorm_sq = f_trial.squaredNorm();
        if ( std::isfinite( new_fnorm_sq ) && new_fnorm_sq <= fnorm_sq + c * alpha * 0.001 )
        {
          x_active  = trial_active;
          f         = f_trial;
          alpha_out = alpha;
          return true;
        }
        alpha *= tau;
      }
      return false;
    }

    bool
    update_and_evaluate( NonlinearSystem const &  sys,
                         Vector &                 x_full,
                         Vector const &           x_active_new,
                         std::vector<int> const & active_indices,
                         Vector &                 f )
    {
      for ( size_t i = 0; i < active_indices.size(); ++i ) { x_full[active_indices[i]] = x_active_new[i]; }
      try
      {
        sys.evaluate( x_full, f );
        ++m_inner_feval;
        return std::isfinite( f.norm() );
      }
      catch ( ... )
      {
        return false;
      }
    }
  };

  /*\
   |    ____        _
   |   / __ \      | |
   |  | |  | |_   _| |_ ___ _ __
   |  | |  | | | | | __/ _ \ '__|
   |  | |__| | |_| | ||  __/ |
   |   \____/ \__,_|\__\___|_|
   |
  \*/

  class TwoLevelSubspaceNewton
  {
  public:
    using Scalar       = double;
    using Vector       = Eigen::VectorXd;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

    enum SelectionStrategy
    {
      RANDOM_UNIFORM  = 0,
      CYCLIC          = 1,
      GREEDY          = 2,
      SCALED_GRADIENT = 3,
      ADAPTIVE        = 4
    };

    enum FallbackStrategy
    {
      NO_FALLBACK       = 0,
      INCREASE_BLOCK    = 1,
      FULL_NEWTON_STEP  = 2,
      GRADIENT_FALLBACK = 3,
      RANDOM_RESTART    = 4
    };

  private:
    int               m_max_outer_iter{ 100 };
    Scalar            m_outer_tol{ 1e-8 };
    int               m_block_size{ 1 };
    SelectionStrategy m_strategy{ SCALED_GRADIENT };
    int               m_verbose{ 0 };

    SubspaceInnerSolver m_inner_solver;
    std::mt19937        m_rng{ std::random_device{}() };
    int                 m_cyclic_index{ 0 };
    int                 m_outer_iter{ 0 };

    // Stats tracking for Outer Loop
    int    m_outer_jeval{ 0 };
    Scalar m_final_residual{ 0.0 };

    SparseMatrix m_cached_jacobian;
    bool         m_jacobian_valid{ false };

    // Stagnation detection
    std::vector<Scalar> m_residual_history;
    std::deque<int>     m_selection_history;
    int                 m_stagnation_counter{ 0 };
    int                 m_max_stagnation_before_reset{ 10 };
    Scalar              m_stagnation_tolerance{ 1e-3 };
    bool                m_adaptive_block_size{ true };
    int                 m_min_block_size{ 1 };
    int                 m_max_block_size{ 50 };
    Scalar              m_progress_threshold{ 0.1 };
    FallbackStrategy    m_fallback_strategy{ INCREASE_BLOCK };

    // Adaptive strategy switching
    bool   m_adaptive_strategy{ true };
    int    m_strategy_switch_frequency{ 5 };
    Scalar m_current_progress{ 0.0 };

  public:
    TwoLevelSubspaceNewton()
    {
      m_inner_solver.set_max_iterations( 5 );
      m_inner_solver.set_tolerance( 1e-2 );
      m_inner_solver.enable_adaptive_damping( true );
    }

    // --- Configuration API (API compatibile con i Test) ---
    void
    set_max_iterations( int it )
    {
      m_max_outer_iter = it;
    }
    void
    set_tolerance( Scalar t )
    {
      m_outer_tol = t;
    }
    void
    set_block_size( int b )
    {
      m_block_size = std::max( 1, b );
    }
    void
    set_verbose_level( int v )
    {
      m_verbose = v;
      m_inner_solver.set_verbose_level( v );
    }
    void
    set_strategy( SelectionStrategy s )
    {
      m_strategy = s;
    }

    // Stagnation and adaptive configuration
    void
    enable_adaptive_block_size( bool enable )
    {
      m_adaptive_block_size = enable;
    }
    void
    set_min_block_size( int s )
    {
      m_min_block_size = std::max( 1, s );
    }
    void
    set_max_block_size( int s )
    {
      m_max_block_size = std::min( s, 1000 );
    }
    void
    set_fallback_strategy( FallbackStrategy s )
    {
      m_fallback_strategy = s;
    }
    void
    set_stagnation_tolerance( Scalar tol )
    {
      m_stagnation_tolerance = tol;
    }
    void
    set_max_stagnation_before_reset( int n )
    {
      m_max_stagnation_before_reset = n;
    }
    void
    enable_adaptive_strategy( bool enable )
    {
      m_adaptive_strategy = enable;
    }

    // Proxy methods for Inner Solver configuration
    void
    set_inner_max_iterations( int i )
    {
      m_inner_solver.set_max_iterations( i );
    }
    void
    set_inner_tolerance( Scalar t )
    {
      m_inner_solver.set_tolerance( t );
    }
    void
    set_inner_damping( Scalar d )
    {
      m_inner_solver.set_damping( d );
    }
    void
    enable_inner_line_search( bool b )
    {
      m_inner_solver.enable_line_search( b );
    }
    void
    enable_inner_adaptive_damping( bool b )
    {
      m_inner_solver.enable_adaptive_damping( b );
    }

    // Direct Access if needed
    SubspaceInnerSolver &
    inner_solver()
    {
      return m_inner_solver;
    }

    // --- Getters for Statistics (API compatibile con i Test) ---
    int
    get_outer_iterations() const
    {
      return m_outer_iter;
    }
    int
    get_total_iterations() const
    {
      return m_inner_solver.get_total_inner_iterations();
    }

    int
    get_function_evals() const
    {
      return 1 + m_inner_solver.get_inner_function_evals();
    }

    int
    get_jacobian_evals() const
    {
      return m_outer_jeval + m_inner_solver.get_inner_jacobian_evals();
    }

    Scalar
    final_residual() const
    {
      return m_final_residual;
    }

    // Get current progress information
    Scalar
    get_current_progress() const
    {
      return m_current_progress;
    }
    int
    get_stagnation_counter() const
    {
      return m_stagnation_counter;
    }

    // --- Main Solver ---
    bool
    solve( NonlinearSystem const & sys, Vector & x )
    {
      int n = x.size();
      m_inner_solver.reset_counters();
      m_outer_iter         = 0;
      m_outer_jeval        = 0;
      m_jacobian_valid     = false;
      m_stagnation_counter = 0;
      m_residual_history.clear();
      m_selection_history.clear();
      m_current_progress = 0.0;

      Vector f( sys.num_equations() );
      try
      {
        sys.evaluate( x, f );
      }
      catch ( ... )
      {
        return false;
      }
      m_final_residual = f.norm();
      m_residual_history.push_back( m_final_residual );

      if ( m_verbose > 0 )
      {
        fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold,
                    "==========================================================\n" );
        fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold,
                    "=== Two-Level Subspace Newton Solver                    ===\n" );
        fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "=== Strategy: {:15} Block size: {:3d} {:>17}===\n",
                    strategy_name( m_strategy ), m_block_size, m_adaptive_block_size ? "(adaptive)" : "(fixed)" );
        fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold,
                    "==========================================================\n" );
        fmt::print( "Initial residual: {:12.3e} | Variables: {:4d} | Equations: {:4d}\n", m_final_residual, n,
                    sys.num_equations() );
        fmt::print( "----------------------------------------------------------\n" );
      }

      int               dynamic_block_size = m_block_size;
      SelectionStrategy current_strategy   = m_strategy;

      for ( m_outer_iter = 0; m_outer_iter < m_max_outer_iter; ++m_outer_iter )
      {
        m_final_residual = f.norm();

        // Stagnation detection
        if ( m_outer_iter > 0 && detect_stagnation() )
        {
          if ( !apply_fallback_strategy( sys, x, f, dynamic_block_size, n, current_strategy ) )
          {
            if ( m_verbose > 0 )
            {
              fmt::print( fg( fmt::color::red ), "Fallback strategies exhausted at iteration {:3d}\n", m_outer_iter );
            }
            break;
          }
          continue;
        }

        m_residual_history.push_back( m_final_residual );
        if ( m_residual_history.size() > 20 ) { m_residual_history.erase( m_residual_history.begin() ); }

        if ( m_final_residual < m_outer_tol )
        {
          if ( m_verbose > 0 )
          {
            fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                        "==========================================================\n" );
            fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                        " CONVERGED at iteration {:3d} with residual {:12.3e}\n", m_outer_iter, m_final_residual );
            fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                        "==========================================================\n" );
          }
          return true;
        }

        // Adaptive strategy switching
        if ( m_adaptive_strategy && m_outer_iter > 0 && m_outer_iter % m_strategy_switch_frequency == 0 )
        {
          SelectionStrategy old_strategy = current_strategy;
          current_strategy               = select_adaptive_strategy();
          if ( m_verbose > 0 && old_strategy != current_strategy )
          {
            fmt::print( fg( fmt::color::yellow ), "Switching strategy from {} to {} at iteration {:3d}\n",
                        strategy_name( old_strategy ), strategy_name( current_strategy ), m_outer_iter );
          }
        }

        std::vector<int> active_indices = select_variables_enhanced( sys, x, f, dynamic_block_size, current_strategy );

        // --- STAMPA MIGLIORATA DELL'ITERATA OUTER ---
        if ( m_verbose > 0 )
        {
          fmt::print( fg( fmt::color::cyan ) | fmt::emphasis::bold, "┌─ Outer iteration {:3d} ", m_outer_iter );
          fmt::print( "|f|={:12.3e} | Block size: {:3d} | Strategy: {:12} \n", f.norm(), dynamic_block_size,
                      strategy_name( current_strategy ) );

          if ( m_verbose > 1 && !active_indices.empty() )
          {
            fmt::print( "├─ Selected indices ({:3d}): ", static_cast<int>( active_indices.size() ) );

            // Formatta gli indici in modo compatto
            if ( active_indices.size() <= 10 )
            {
              for ( size_t i = 0; i < active_indices.size(); ++i )
              {
                if ( i > 0 ) fmt::print( ", " );
                fmt::print( "{:3d}", active_indices[i] );
              }
              fmt::print( "\n" );
            }
            else
            {
              // Mostra i primi 5 e gli ultimi 5
              for ( int i = 0; i < 5; ++i )
              {
                if ( i > 0 ) fmt::print( ", " );
                fmt::print( "{:3d}", active_indices[i] );
              }
              fmt::print( ", ..., " );
              for ( size_t i = active_indices.size() - 5; i < active_indices.size(); ++i )
              {
                if ( i > active_indices.size() - 5 ) fmt::print( ", " );
                fmt::print( "{:3d}", active_indices[i] );
              }
              fmt::print( " (total {:3d})\n", static_cast<int>( active_indices.size() ) );
            }
          }

          if ( m_verbose == 1 )
          {
            // Per verbose=1, stampa solo un riepilogo ogni 5 iterazioni
            if ( m_outer_iter % 5 == 0 )
            {
              fmt::print( "├─ Progress: {:6.2f}% | Stagnation: {:2d}/{:2d}\n", m_current_progress * 100,
                          m_stagnation_counter, m_max_stagnation_before_reset );
            }
          }

          if ( m_verbose > 1 ) { fmt::print( "└─ Starting inner solver...\n" ); }
        }
        // --------------------------------------------

        // Update selection history
        for ( int idx : active_indices )
        {
          m_selection_history.push_back( idx );
          if ( m_selection_history.size() > 5 * n ) { m_selection_history.pop_front(); }
        }

        bool inner_ok = m_inner_solver.solve_in_subspace( sys, x, active_indices, f, m_outer_iter );

        m_jacobian_valid = false;

        if ( !inner_ok )
        {
          if ( m_verbose > 0 )
          {
            fmt::print( fg( fmt::color::orange ),
                        "  Inner solver failed at outer iteration {:3d}. Activating fallback.\n", m_outer_iter );
          }
          m_stagnation_counter++;
        }
        else
        {
          // Reset stagnation counter on success
          if ( m_final_residual > f.norm() ) { m_stagnation_counter = 0; }
        }

        // Progress report per verbose=1 (ogni 5 iterazioni)
        if ( m_verbose == 1 && m_outer_iter % 5 == 0 && m_outer_iter > 0 )
        {
          fmt::print( "[{:3d}] |f|={:12.3e} | Block={:3d} | Progress={:6.2f}% | Stagn={:2d}\n", m_outer_iter, f.norm(),
                      dynamic_block_size, m_current_progress * 100, m_stagnation_counter );
        }

        // Progress report per verbose>1 (ogni iterazione, dopo il solver interno)
        if ( m_verbose > 1 )
        {
          fmt::print( fg( fmt::color::dark_gray ), "  ──────────────────────────────────────────────────────\n" );
        }
      }

      m_final_residual = f.norm();
      bool converged   = m_final_residual < m_outer_tol;

      if ( m_verbose > 0 )
      {
        if ( converged )
        {
          fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                      "==========================================================\n" );
          fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                      " CONVERGED after {:3d} iterations | Final residual: {:12.3e}\n", m_outer_iter,
                      m_final_residual );
          fmt::print( fg( fmt::color::lime_green ) | fmt::emphasis::bold,
                      "==========================================================\n" );
        }
        else
        {
          fmt::print( fg( fmt::color::orange ) | fmt::emphasis::bold,
                      "==========================================================\n" );
          fmt::print( fg( fmt::color::orange ) | fmt::emphasis::bold,
                      " STOPPED after {:3d} iterations | Final residual: {:12.3e}\n", m_outer_iter, m_final_residual );
          fmt::print( fg( fmt::color::orange ) | fmt::emphasis::bold,
                      "==========================================================\n" );
        }

        // Stampa statistiche finali
        fmt::print( "\nStatistics:\n" );
        fmt::print( "  Outer iterations:     {:6d}\n", m_outer_iter );
        fmt::print( "  Total inner iters:    {:6d}\n", m_inner_solver.get_total_inner_iterations() );
        fmt::print( "  Function evaluations: {:6d}\n", get_function_evals() );
        fmt::print( "  Jacobian evaluations: {:6d}\n", get_jacobian_evals() );
        fmt::print( "  Final residual:       {:12.3e}\n", m_final_residual );
        fmt::print( "  Tolerance:            {:12.3e}\n", m_outer_tol );
        fmt::print( "\n" );
      }

      return converged;
    }

  private:
    std::vector<int>
    select_variables_enhanced( NonlinearSystem const & sys,
                               Vector const &          x,
                               Vector const &          f,
                               int &                   dynamic_block_size,
                               SelectionStrategy       current_strategy )
    {
      int n = x.size();

      // Adaptive block sizing based on progress
      if ( m_adaptive_block_size && m_outer_iter > 0 )
      {
        dynamic_block_size = adjust_block_size_based_on_progress( dynamic_block_size, n );
      }

      if ( dynamic_block_size >= n )
      {
        std::vector<int> all_indices( n );
        std::iota( all_indices.begin(), all_indices.end(), 0 );
        return all_indices;
      }

      std::vector<int> indices;

      if ( current_strategy == RANDOM_UNIFORM )
      {
        std::vector<int> all( n );
        std::iota( all.begin(), all.end(), 0 );
        std::shuffle( all.begin(), all.end(), m_rng );
        indices.assign( all.begin(), all.begin() + dynamic_block_size );
      }
      else if ( current_strategy == CYCLIC )
      {
        for ( int i = 0; i < dynamic_block_size; ++i ) indices.push_back( ( m_cyclic_index + i ) % n );
        m_cyclic_index = ( m_cyclic_index + dynamic_block_size ) % n;
      }
      else
      {
        if ( !m_jacobian_valid )
        {
          try
          {
            sys.jacobian( x, m_cached_jacobian );
            m_jacobian_valid = true;
            m_outer_jeval++;
          }
          catch ( ... )
          {
            // Fallback a random se Jacobian fallisce
            std::vector<int> all( n );
            std::iota( all.begin(), all.end(), 0 );
            std::shuffle( all.begin(), all.end(), m_rng );
            indices.assign( all.begin(), all.begin() + dynamic_block_size );
            return indices;
          }
        }

        Vector criteria( n );
        Vector grad = m_cached_jacobian.transpose() * f;

        if ( current_strategy == GREEDY ) { criteria = grad.cwiseAbs(); }
        else if ( current_strategy == SCALED_GRADIENT || current_strategy == ADAPTIVE )
        {
          // SCALED: |Gradient_i| / ||Column_i||
          Vector col_norms = Vector::Zero( n );
          for ( int k = 0; k < m_cached_jacobian.outerSize(); ++k )
          {
            for ( typename SparseMatrix::InnerIterator it( m_cached_jacobian, k ); it; ++it )
            {
              col_norms[it.col()] += it.value() * it.value();
            }
          }
          col_norms = col_norms.cwiseSqrt();

          for ( int i = 0; i < n; ++i ) { criteria[i] = std::abs( grad[i] ) / ( col_norms[i] + 1e-8 ); }
        }

        // Penalizza variabili selezionate recentemente
        std::unordered_map<int, int> recent_count;
        for ( int idx : m_selection_history ) { recent_count[idx]++; }

        // Aggiungi penalizzazione basata su selezione recente
        if ( m_selection_history.size() > n )
        {
          for ( int i = 0; i < n; ++i ) { criteria[i] /= ( 1.0 + 0.1 * recent_count[i] ); }
        }

        std::vector<std::pair<Scalar, int>> importance( n );
        for ( int i = 0; i < n; ++i ) importance[i] = { criteria[i], i };

        std::partial_sort( importance.begin(), importance.begin() + dynamic_block_size, importance.end(),
                           std::greater<std::pair<Scalar, int>>() );

        for ( int i = 0; i < dynamic_block_size; ++i ) indices.push_back( importance[i].second );
      }

      return indices;
    }

    int
    adjust_block_size_based_on_progress( int current_size, int n )
    {
      if ( m_residual_history.size() < 3 ) return current_size;

      m_current_progress = compute_recent_progress();

      if ( m_current_progress < m_progress_threshold )
      {
        // Progresso lento: aumenta il block size
        int new_size = std::min( current_size * 2, m_max_block_size );
        new_size     = std::min( new_size, n );

        if ( m_verbose > 1 )
        {
          fmt::print( fg( fmt::color::yellow ),
                      "  Slow progress ({:6.2f}%), increasing block size from {:3d} to {:3d}\n",
                      m_current_progress * 100, current_size, new_size );
        }
        return new_size;
      }
      else if ( m_current_progress > 2.0 * m_progress_threshold && current_size > m_min_block_size )
      {
        // Progresso buono: riduci il block size per efficienza
        int new_size = std::max( current_size / 2, m_min_block_size );
        if ( m_verbose > 1 )
        {
          fmt::print( fg( fmt::color::green ),
                      "  Good progress ({:6.2f}%), decreasing block size from {:3d} to {:3d}\n",
                      m_current_progress * 100, current_size, new_size );
        }
        return new_size;
      }

      return current_size;
    }

    Scalar
    compute_recent_progress() const
    {
      if ( m_residual_history.size() < 3 ) return 0.0;

      size_t last     = m_residual_history.size() - 1;
      Scalar progress = 0.0;
      int    count    = 0;

      // Media del progresso negli ultimi 3 passi
      for ( size_t i = std::max( 0, (int) last - 2 ); i < last; ++i )
      {
        if ( m_residual_history[i] > 1e-16 )
        {
          progress += ( m_residual_history[i] - m_residual_history[i + 1] ) / m_residual_history[i];
          count++;
        }
      }

      return count > 0 ? progress / count : 0.0;
    }

    bool
    detect_stagnation()
    {
      if ( m_residual_history.size() < 5 ) return false;

      size_t last          = m_residual_history.size() - 1;
      Scalar avg_reduction = 0.0;
      int    count         = 0;

      // Calcola riduzione media negli ultimi 5 passi
      for ( size_t i = std::max( 0, (int) last - 4 ); i < last; ++i )
      {
        if ( m_residual_history[i] > 1e-16 )
        {
          Scalar reduction = ( m_residual_history[i] - m_residual_history[i + 1] ) / m_residual_history[i];
          if ( reduction > 0 )
          {
            avg_reduction += reduction;
            count++;
          }
        }
      }

      if ( count > 0 )
      {
        avg_reduction /= count;
        if ( avg_reduction < m_stagnation_tolerance )
        {
          m_stagnation_counter++;
          if ( m_verbose > 1 )
          {
            fmt::print( fg( fmt::color::yellow ),
                        "  Stagnation detected: avg reduction = {:8.3e}, counter = {:2d}/{:2d}\n", avg_reduction,
                        m_stagnation_counter, m_max_stagnation_before_reset );
          }
          return true;
        }
      }

      return false;
    }

    bool
    apply_fallback_strategy( NonlinearSystem const & sys,
                             Vector &                x,
                             Vector &                f,
                             int &                   block_size,
                             int                     n,
                             SelectionStrategy &     current_strategy )
    {
      if ( m_stagnation_counter > m_max_stagnation_before_reset )
      {
        if ( m_verbose > 0 )
        {
          fmt::print( fg( fmt::color::orange ) | fmt::emphasis::bold,
                      "  Maximum stagnation reached ({:2d}/{:2d}). Applying fallback strategy: {}.\n",
                      m_stagnation_counter, m_max_stagnation_before_reset,
                      fallback_strategy_name( m_fallback_strategy ) );
        }

        switch ( m_fallback_strategy )
        {
          case INCREASE_BLOCK:
            block_size = std::min( block_size * 2, n );
            if ( m_verbose > 0 ) { fmt::print( "  Increased block size to {:3d}\n", block_size ); }
            m_stagnation_counter = 0;
            return true;

          case FULL_NEWTON_STEP:
            // Un passo di Newton completo usando SP_TikhonovSolver2
            try
            {
              SparseMatrix J;
              sys.jacobian( x, J );

              // Usa regolarizzazione diagonale basata sulle norme delle colonne
              Vector D = Vector::Ones( J.cols() );
              for ( int i = 0; i < J.cols(); ++i )
              {
                Scalar col_norm_sq = J.col( i ).squaredNorm();
                D[i]               = 1e-6 * ( col_norm_sq + 1e-8 );
              }

              SP_TikhonovSolver2<Scalar> solver( J, D );
              Vector                     dx = solver.solve( -f );

              // Line search per il passo completo
              Scalar alpha = 1.0;
              if ( line_search_full( sys, x, f, dx, alpha ) )
              {
                x += alpha * dx;
                sys.evaluate( x, f );
                m_residual_history.push_back( f.norm() );
                m_stagnation_counter = 0;
                if ( m_verbose > 0 )
                {
                  fmt::print( fg( fmt::color::green ), "  Full Newton step successful, new residual: {:12.3e}\n",
                              f.norm() );
                }
                return true;
              }
              else
              {
                if ( m_verbose > 0 )
                {
                  fmt::print( fg( fmt::color::orange ), "  Full Newton step line search failed\n" );
                }
              }
            }
            catch ( ... )
            {
              // Fallback al caso successivo
              if ( m_verbose > 0 )
              {
                fmt::print( fg( fmt::color::orange ), "  Full Newton step computation failed\n" );
              }
            }
            break;

          case GRADIENT_FALLBACK:
            // Usa direzione del gradiente
            try
            {
              SparseMatrix J;
              sys.jacobian( x, J );
              Vector grad  = J.transpose() * f;
              Scalar gnorm = grad.norm();
              if ( gnorm > 1e-12 )
              {
                Vector dx = -grad / gnorm;  // Direzione normalizzata

                Scalar alpha = 1.0;
                for ( int ls = 0; ls < 10; ++ls )
                {
                  Vector x_trial = x + alpha * dx;
                  Vector f_trial;
                  sys.evaluate( x_trial, f_trial );

                  if ( f_trial.norm() < f.norm() )
                  {
                    x = x_trial;
                    f = f_trial;
                    m_residual_history.push_back( f.norm() );
                    m_stagnation_counter = 0;
                    if ( m_verbose > 0 )
                    {
                      fmt::print( fg( fmt::color::green ), "  Gradient fallback successful, new residual: {:12.3e}\n",
                                  f.norm() );
                    }
                    return true;
                  }
                  alpha *= 0.5;
                }
              }
            }
            catch ( ... )
            {
              // Fallback
            }
            break;

          case RANDOM_RESTART:
            // Piccola perturbazione random
            std::normal_distribution<Scalar> dist( 0.0, 0.01 );
            for ( int i = 0; i < x.size(); ++i ) { x[i] += dist( m_rng ); }
            sys.evaluate( x, f );
            m_residual_history.push_back( f.norm() );
            m_stagnation_counter = 0;

            // Cambia strategia per evitare loop
            current_strategy = static_cast<SelectionStrategy>( ( current_strategy + 1 ) % 5 );
            if ( m_verbose > 0 )
            {
              fmt::print( fg( fmt::color::yellow ), "  Random restart applied, new strategy: {}\n",
                          strategy_name( current_strategy ) );
            }
            return true;
        }

        // Tutte le strategie fallite: aumenta counter e continua
        return false;
      }

      return false;
    }

    bool
    line_search_full( NonlinearSystem const & sys,
                      Vector const &          x,
                      Vector const &          f,
                      Vector const &          dx,
                      Scalar &                alpha )
    {
      Scalar f0      = f.squaredNorm();
      Scalar c       = 1e-4;
      Vector x_trial = x;
      Vector f_trial;

      for ( int i = 0; i < 10; ++i )
      {
        x_trial = x + alpha * dx;
        sys.evaluate( x_trial, f_trial );

        // Armijo condition
        if ( f_trial.squaredNorm() <= f0 + c * alpha * f.dot( dx ) ) { return true; }
        alpha *= 0.5;
      }
      return false;
    }

    SelectionStrategy
    select_adaptive_strategy() const
    {
      if ( m_current_progress < 0.01 )
      {
        // Progresso molto lento: prova strategia più aggressiva
        return SCALED_GRADIENT;
      }
      else if ( m_current_progress > 0.1 )
      {
        // Progresso buono: usa strategia più efficiente
        return GREEDY;
      }
      else
      {
        // Progresso moderato: mantieni strategia corrente o ciclica
        return ( m_outer_iter % 3 == 0 ) ? CYCLIC : m_strategy;
      }
    }

    const char *
    strategy_name( SelectionStrategy s ) const
    {
      switch ( s )
      {
        case RANDOM_UNIFORM:
          return "Random";
        case CYCLIC:
          return "Cyclic";
        case GREEDY:
          return "Greedy";
        case SCALED_GRADIENT:
          return "ScaledGrad";
        case ADAPTIVE:
          return "Adaptive";
        default:
          return "?";
      }
    }

    const char *
    fallback_strategy_name( FallbackStrategy s ) const
    {
      switch ( s )
      {
        case NO_FALLBACK:
          return "NoFallback";
        case INCREASE_BLOCK:
          return "IncreaseBlock";
        case FULL_NEWTON_STEP:
          return "FullNewtonStep";
        case GRADIENT_FALLBACK:
          return "GradientFallback";
        case RANDOM_RESTART:
          return "RandomRestart";
        default:
          return "?";
      }
    }
  };
}  // namespace Utils
