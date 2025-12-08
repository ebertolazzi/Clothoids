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
// file: Utils_NonlinearKaczmarz.hh
//

#pragma once

#ifndef UTILS_NONLINEAR_KACZMARZ_dot_HH
#define UTILS_NONLINEAR_KACZMARZ_dot_HH

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"
#include "Utils_nonlinear_system.hh"

namespace Utils
{

  class NonlinearKaczmarz
  {
  public:
    using integer      = Eigen::Index;
    using real_type    = double;
    using Vector       = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    using Matrix       = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<real_type>;

    enum SelectionStrategy
    {
      CYCLIC,           // Selezione ciclica delle equazioni
      RANDOM_UNIFORM,   // Selezione casuale uniforme
      RANDOM_WEIGHTED,  // Selezione casuale pesata per residuo
      GREEDY            // Selezione dell'equazione con residuo massimo
    };

  private:
    // Parametri del metodo
    real_type         m_tolerance          = 1e-8;
    real_type         m_relative_tolerance = 1e-8;  // Tolleranza relativa
    integer           m_max_iterations     = 10000;
    integer           m_max_function_evals = 100000;
    real_type         m_relaxation         = 1;  // Parametro di rilassamento (0 < λ ≤ 2)
    SelectionStrategy m_strategy           = RANDOM_UNIFORM;
    bool              m_use_line_search    = false;  // Usa line search per ottimizzare λ
    real_type         m_line_search_beta   = 0.5;    // Parametro di riduzione per line search
    real_type         m_line_search_c1     = 1e-4;   // Condizione di Armijo
    integer           m_verbose_level      = 1;      // 0=silent, 1=summary, 2=detailed

    // Statistiche
    integer   m_num_iterations     = 0;
    integer   m_num_function_evals = 0;
    integer   m_num_jacobian_evals = 0;
    integer   m_num_line_searches  = 0;
    bool      m_converged          = false;
    real_type m_final_residual     = 0;

    // Generatore casuale
    mutable std::mt19937                              m_random_engine;
    mutable std::uniform_real_distribution<real_type> m_uniform_dist{ 0.0, 1.0 };

    // Cache (per evitare allocazioni ripetute)
    mutable std::vector<real_type> m_residuals;
    mutable Vector                 m_grad_cache;
    mutable Vector                 m_x_new_cache;
    mutable Vector                 m_f_new_cache;

    // Costanti
    static constexpr real_type EPSILON              = 1e-16;
    static constexpr real_type MIN_RELAXATION       = 1e-8;
    static constexpr real_type MAX_RELAXATION       = 2.0;
    static constexpr integer   MAX_LINE_SEARCH_ITER = 10;

  public:
    NonlinearKaczmarz()
    {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      m_random_engine.seed( seed );
    }

    // ========================================================================
    // SETTERS CON VALIDAZIONE
    // ========================================================================

    void
    set_tolerance( real_type tol )
    {
      UTILS_ASSERT( tol > 0, "Tolerance must be positive" );
      m_tolerance = tol;
    }

    void
    set_relative_tolerance( real_type tol )
    {
      UTILS_ASSERT( tol > 0, "Relative tolerance must be positive" );
      m_relative_tolerance = tol;
    }

    void
    set_max_iterations( integer max_iter )
    {
      UTILS_ASSERT( max_iter > 0, "Max iterations must be positive" );
      m_max_iterations = max_iter;
    }

    void
    set_max_function_evals( integer max_feval )
    {
      UTILS_ASSERT( max_feval > 0, "Max function evaluations must be positive" );
      m_max_function_evals = max_feval;
    }

    void
    set_relaxation( real_type lambda )
    {
      UTILS_ASSERT( lambda > 0 && lambda <= MAX_RELAXATION, "Relaxation parameter must be in (0, 2]" );
      m_relaxation = lambda;
    }

    void
    set_strategy( SelectionStrategy s )
    {
      m_strategy = s;
    }

    void
    set_use_line_search( bool use )
    {
      m_use_line_search = use;
    }

    void
    set_line_search_beta( real_type beta )
    {
      UTILS_ASSERT( beta > 0 && beta < 1, "Beta must be in (0, 1)" );
      m_line_search_beta = beta;
    }

    void
    set_line_search_c1( real_type c1 )
    {
      UTILS_ASSERT( c1 > 0 && c1 < 0.5, "c1 must be in (0, 0.5)" );
      m_line_search_c1 = c1;
    }

    void
    set_random_seed( unsigned seed )
    {
      m_random_engine.seed( seed );
    }

    void
    set_verbose_level( integer level )
    {
      UTILS_ASSERT( level >= 0 && level <= 2, "Verbose level must be 0, 1, or 2" );
      m_verbose_level = level;
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    integer
    get_num_iterations() const
    {
      return m_num_iterations;
    }
    integer
    get_num_function_evals() const
    {
      return m_num_function_evals;
    }
    integer
    get_num_jacobian_evals() const
    {
      return m_num_jacobian_evals;
    }
    integer
    get_num_line_searches() const
    {
      return m_num_line_searches;
    }
    bool
    has_converged() const
    {
      return m_converged;
    }
    real_type
    get_final_residual() const
    {
      return m_final_residual;
    }
    real_type
    get_tolerance() const
    {
      return m_tolerance;
    }
    real_type
    get_relaxation() const
    {
      return m_relaxation;
    }

    // ========================================================================
    // METODI PRIVATI
    // ========================================================================

  private:
    //! Seleziona l'equazione in base alla strategia scelta
    integer
    select_equation( const Vector & f ) const
    {
      integer m = f.size();

      switch ( m_strategy )
      {
        case CYCLIC:
          return m_num_iterations % m;

        case RANDOM_UNIFORM:
        {
          std::uniform_int_distribution<integer> dist( 0, m - 1 );
          return dist( m_random_engine );
        }

        case RANDOM_WEIGHTED:
        {
          // Calcola residui assoluti
          m_residuals.resize( m );
          real_type sum_residuals = 0.0;

          for ( integer i = 0; i < m; ++i )
          {
            m_residuals[i] = std::abs( f( i ) );
            sum_residuals += m_residuals[i];
          }

          // Se tutti i residui sono quasi zero, selezione uniforme
          if ( sum_residuals < EPSILON )
          {
            std::uniform_int_distribution<integer> dist( 0, m - 1 );
            return dist( m_random_engine );
          }

          // Normalizza per creare distribuzione di probabilità
          for ( integer i = 0; i < m; ++i ) { m_residuals[i] /= sum_residuals; }

          std::discrete_distribution<integer> dist( m_residuals.begin(), m_residuals.end() );
          return dist( m_random_engine );
        }

        case GREEDY:
        {
          // Trova l'equazione con residuo massimo
          integer   max_idx = 0;
          real_type max_val = std::abs( f( 0 ) );

          for ( integer i = 1; i < m; ++i )
          {
            real_type val = std::abs( f( i ) );
            if ( val > max_val )
            {
              max_val = val;
              max_idx = i;
            }
          }
          return max_idx;
        }

        default:
          return m_num_iterations % m;
      }
    }

    //! Calcola il gradiente della i-esima equazione (efficiente)
    void
    compute_gradient( NonlinearSystem & system, const Vector & x, integer eq_idx, Vector & grad )
    {
      integer n = system.num_equations();

      // Alloca Jacobiano solo se necessario
      if ( m_grad_cache.size() != n ) { m_grad_cache.resize( n ); }

      SparseMatrix J( n, n );
      system.jacobian( x, J );
      m_num_jacobian_evals++;

      // Estrai la riga corrispondente all'equazione
      grad = J.row( eq_idx ).transpose();
    }

    //! Line search di Armijo per Kaczmarz
    real_type
    armijo_line_search( NonlinearSystem & system,
                        const Vector &    x,
                        const Vector &    f,
                        real_type         norm_f,
                        integer           eq_idx,
                        const Vector &    grad,
                        real_type         lambda_init )
    {
      integer n = x.size();

      // Preallocazione cache
      if ( m_x_new_cache.size() != n )
      {
        m_x_new_cache.resize( n );
        m_f_new_cache.resize( n );
      }

      real_type fi           = f( eq_idx );
      real_type grad_norm_sq = grad.squaredNorm();

      if ( grad_norm_sq < EPSILON ) return 0.0;

      // Derivata direzionale: ∇f_i(x)·d = -|f_i(x)|²/||∇f_i(x)||²
      real_type directional_derivative = -fi * fi / grad_norm_sq;
      real_type step_direction         = fi / grad_norm_sq;

      real_type lambda = lambda_init;

      for ( integer ls_iter = 0; ls_iter < MAX_LINE_SEARCH_ITER; ++ls_iter )
      {
        // Nuovo punto
        m_x_new_cache.noalias() = x - lambda * step_direction * grad;

        // Valuta funzione
        try
        {
          system.check_if_admissible( m_x_new_cache );
          system.evaluate( m_x_new_cache, m_f_new_cache );
          m_num_function_evals++;

          real_type norm_f_new         = m_f_new_cache.norm();
          real_type actual_reduction   = norm_f - norm_f_new;
          real_type required_reduction = m_line_search_c1 * lambda * directional_derivative;

          if ( actual_reduction >= required_reduction )
          {
            m_num_line_searches++;
            return lambda;
          }
        }
        catch ( ... )
        {
          // Punto non ammissibile, continua con passo ridotto
        }

        // Riduci il passo
        lambda *= m_line_search_beta;

        if ( lambda < MIN_RELAXATION ) break;
      }

      return 0.0;  // Line search fallita
    }

    //! Verifica convergenza
    bool
    check_convergence( real_type norm_f, real_type initial_norm ) const
    {
      // Convergenza assoluta
      if ( norm_f < m_tolerance ) return true;

      // Convergenza relativa
      if ( initial_norm > EPSILON && norm_f / initial_norm < m_relative_tolerance ) { return true; }

      return false;
    }

    //! Nome della strategia (per output)
    const char *
    strategy_name( SelectionStrategy s ) const
    {
      switch ( s )
      {
        case CYCLIC:
          return "Cyclic";
        case RANDOM_UNIFORM:
          return "Random Uniform";
        case RANDOM_WEIGHTED:
          return "Random Weighted";
        case GREEDY:
          return "Greedy";
        default:
          return "Unknown";
      }
    }

  public:
    // ========================================================================
    // METODO PRINCIPALE DI RISOLUZIONE
    // ========================================================================

    bool
    solve( NonlinearSystem & system, Vector & x )
    {
      integer n = system.num_equations();

      UTILS_ASSERT( n > 0, "System must have at least one equation" );
      UTILS_ASSERT( x.size() == n, "Initial guess size mismatch" );

      Vector f( n );
      Vector grad( n );

      // Reset statistiche
      m_num_iterations     = 0;
      m_num_function_evals = 0;
      m_num_jacobian_evals = 0;
      m_num_line_searches  = 0;
      m_converged          = false;

      // Valuta funzione iniziale
      try
      {
        system.check_if_admissible( x );
        system.evaluate( x, f );
        m_num_function_evals++;
      }
      catch ( const std::exception & e )
      {
        if ( m_verbose_level > 0 )
        {
          fmt::print( fmt::fg( fmt::color::red ), "✗ Initial point not admissible: {}\n", e.what() );
        }
        return false;
      }

      real_type norm_f       = f.norm();
      real_type initial_norm = norm_f;

      if ( m_verbose_level > 0 )
      {
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
        fmt::print( fmt::fg( fmt::color::light_blue ), "Starting Nonlinear Kaczmarz\n" );
        fmt::print( "Initial residual: {:.6e}\n", initial_norm );
        fmt::print( "Strategy: {}\n", strategy_name( m_strategy ) );
        if ( m_use_line_search ) fmt::print( "Line search: Enabled\n" );
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
      }

      // Iterazioni principali
      for ( m_num_iterations = 0; m_num_iterations < m_max_iterations && m_num_function_evals < m_max_function_evals;
            ++m_num_iterations )
      {
        // Check convergenza
        if ( check_convergence( norm_f, initial_norm ) )
        {
          m_converged      = true;
          m_final_residual = norm_f;

          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::green ), "✓ Converged below tolerance ({:.2e})\n", m_tolerance );
            print_summary( initial_norm, norm_f );
          }
          return true;
        }

        // Stampa progresso (simile a NewtonDumped)
        if ( m_verbose_level == 2 )
        {
          real_type reduction = initial_norm > EPSILON ? initial_norm / norm_f : 1.0;
          fmt::print( fmt::fg( fmt::color::light_blue ), "[{:5}] ‖f‖ = {:.2e}",
                      fmt::format( "{}", m_num_iterations + 1 ), norm_f );
          fmt::print( ", reduction = {:.3e}\n", reduction );
        }

        // Seleziona equazione
        integer eq_idx = select_equation( f );

        // Calcola gradiente dell'equazione selezionata
        compute_gradient( system, x, eq_idx, grad );

        real_type grad_norm_sq = grad.squaredNorm();

        // Skip se gradiente troppo piccolo
        if ( grad_norm_sq < EPSILON )
        {
          if ( m_verbose_level > 1 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "  ⚠ Near-zero gradient at iteration {}, skipping\n",
                        m_num_iterations );
          }
          continue;
        }

        // Calcola passo di Kaczmarz
        real_type fi     = f( eq_idx );
        real_type lambda = m_relaxation;

        // Applica line search se richiesto
        if ( m_use_line_search )
        {
          if ( m_verbose_level == 2 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "  Line search: λ₀={:.3f}\n", lambda );
          }
          lambda = armijo_line_search( system, x, f, norm_f, eq_idx, grad, m_relaxation );
          if ( lambda <= 0.0 )
          {
            if ( m_verbose_level > 1 ) { fmt::print( fmt::fg( fmt::color::red ), "  ✗ Line search failed\n" ); }
            lambda = m_relaxation;  // Fallback
          }
          else if ( m_verbose_level == 2 )
          {
            fmt::print( fmt::fg( fmt::color::green ), "  ✓ Line search: λ={:.3f}\n", lambda );
          }
        }

        // Aggiorna la soluzione: x_{k+1} = x_k - λ * (f_i / ||∇f_i||²) * ∇f_i
        x.noalias() -= lambda * ( fi / grad_norm_sq ) * grad;

        // Rivaluta funzione
        try
        {
          system.check_if_admissible( x );
          system.evaluate( x, f );
          m_num_function_evals++;
          norm_f = f.norm();
        }
        catch ( const std::exception & e )
        {
          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::red ), "✗ Error at iteration {}: {}\n", m_num_iterations, e.what() );
          }
          m_final_residual = norm_f;
          return false;
        }

        // Criterio di arresto per progresso insufficiente
        if ( m_num_iterations > 100 && m_num_iterations % 100 == 0 )
        {
          real_type avg_reduction = std::pow( norm_f / initial_norm, 1.0 / m_num_iterations );
          if ( avg_reduction > 0.9995 )
          {
            if ( m_verbose_level > 0 )
            {
              fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Stopping: insufficient progress (avg reduction: {:.4f})\n",
                          avg_reduction );
            }
            break;
          }
        }
      }

      m_final_residual = norm_f;

      if ( m_verbose_level > 0 )
      {
        if ( m_num_iterations >= m_max_iterations )
        {
          fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Max iterations ({}) reached\n", m_max_iterations );
        }
        else if ( m_num_function_evals >= m_max_function_evals )
        {
          fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Max function evaluations ({}) reached\n",
                      m_max_function_evals );
        }
        print_summary( initial_norm, norm_f );
      }

      return m_converged;
    }

    // ========================================================================
    // VARIANTE BLOCK KACZMARZ
    // ========================================================================

    bool
    solve_block( NonlinearSystem & system, Vector & x, integer block_size = 5 )
    {
      integer n = system.num_equations();

      UTILS_ASSERT( n > 0, "System must have at least one equation" );
      UTILS_ASSERT( x.size() == n, "Initial guess size mismatch" );
      UTILS_ASSERT( block_size > 0 && block_size <= n, "Block size must be in [1, n]" );

      Vector f( n );
      Vector grad( n );

      // Reset statistiche
      m_num_iterations     = 0;
      m_num_function_evals = 0;
      m_num_jacobian_evals = 0;
      m_converged          = false;

      // Valuta funzione iniziale
      try
      {
        system.check_if_admissible( x );
        system.evaluate( x, f );
      }
      catch ( const std::exception & e )
      {
        if ( m_verbose_level > 0 )
        {
          fmt::print( fmt::fg( fmt::color::red ), "✗ Initial point not admissible: {}\n", e.what() );
        }
        return false;
      }

      m_num_function_evals++;
      real_type norm_f       = f.norm();
      real_type initial_norm = norm_f;

      if ( m_verbose_level > 0 )
      {
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
        fmt::print( fmt::fg( fmt::color::light_blue ), "Starting Block Kaczmarz (block_size={})\n", block_size );
        fmt::print( "Initial residual: {:.6e}\n", initial_norm );
        fmt::print( "Strategy: {}\n", strategy_name( m_strategy ) );
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
      }

      for ( m_num_iterations = 0; m_num_iterations < m_max_iterations; ++m_num_iterations )
      {
        if ( check_convergence( norm_f, initial_norm ) )
        {
          m_converged      = true;
          m_final_residual = norm_f;

          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::green ), "✓ Converged below tolerance ({:.2e})\n", m_tolerance );
            print_summary_block( initial_norm, norm_f );
          }
          return true;
        }

        // Stampa progresso
        if ( m_verbose_level == 2 )
        {
          real_type reduction = initial_norm > EPSILON ? initial_norm / norm_f : 1.0;
          fmt::print( fmt::fg( fmt::color::light_blue ), "[{:5}] ‖f‖ = {:.2e}",
                      fmt::format( "{}", m_num_iterations + 1 ), norm_f );
          fmt::print( ", reduction = {:.3e}\n", reduction );
        }

        // Seleziona blocco di equazioni
        std::vector<integer> block_indices = select_block( f, n, block_size );

        // Per ogni equazione nel blocco, applica Kaczmarz
        for ( integer eq_idx : block_indices )
        {
          compute_gradient( system, x, eq_idx, grad );
          real_type grad_norm_sq = grad.squaredNorm();

          if ( grad_norm_sq > EPSILON )
          {
            real_type fi = f( eq_idx );
            x.noalias() -= m_relaxation * ( fi / grad_norm_sq ) * grad;
          }
        }

        // Rivaluta funzione dopo il blocco
        try
        {
          system.check_if_admissible( x );
          system.evaluate( x, f );
        }
        catch ( const std::exception & e )
        {
          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::red ), "✗ Error at iteration {}: {}\n", m_num_iterations, e.what() );
          }
          m_final_residual = norm_f;
          return false;
        }

        m_num_function_evals++;
        norm_f = f.norm();
      }

      m_final_residual = norm_f;

      if ( m_verbose_level > 0 )
      {
        if ( m_num_iterations >= m_max_iterations )
        {
          fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Max iterations ({}) reached\n", m_max_iterations );
        }
        print_summary_block( initial_norm, norm_f );
      }

      return false;
    }

  private:
    //! Seleziona un blocco di equazioni
    std::vector<integer>
    select_block( const Vector & f, integer n, integer block_size ) const
    {
      std::vector<integer> block_indices;
      block_indices.reserve( block_size );

      switch ( m_strategy )
      {
        case RANDOM_UNIFORM:
        {
          std::uniform_int_distribution<integer> dist( 0, n - 1 );
          for ( integer i = 0; i < block_size; ++i ) { block_indices.push_back( dist( m_random_engine ) ); }
          break;
        }

        case CYCLIC:
        {
          integer start = ( m_num_iterations * block_size ) % n;
          for ( integer i = 0; i < block_size; ++i ) { block_indices.push_back( ( start + i ) % n ); }
          break;
        }

        case GREEDY:
        case RANDOM_WEIGHTED:
        {
          // Seleziona equazioni con residuo massimo
          m_residuals.resize( n );
          std::vector<integer> sorted_indices( n );

          for ( integer i = 0; i < n; ++i )
          {
            m_residuals[i]    = std::abs( f( i ) );
            sorted_indices[i] = i;
          }

          // Ordina parzialmente per trovare i primi block_size elementi
          std::partial_sort( sorted_indices.begin(), sorted_indices.begin() + block_size, sorted_indices.end(),
                             [this]( integer i, integer j ) { return m_residuals[i] > m_residuals[j]; } );

          for ( integer i = 0; i < block_size; ++i ) { block_indices.push_back( sorted_indices[i] ); }
          break;
        }
      }

      return block_indices;
    }

    //! Stampa riepilogo per il metodo standard
    void
    print_summary( real_type initial_norm, real_type final_norm ) const
    {
      if ( m_verbose_level == 0 ) return;

      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n" );
      fmt::print( "Nonlinear Kaczmarz - Summary\n" );
      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n" );

      if ( m_converged ) { fmt::print( fmt::fg( fmt::color::green ), "Status:             ✓ CONVERGED\n" ); }
      else
      {
        fmt::print( fmt::fg( fmt::color::red ), "Status:             ✗ NOT CONVERGED\n" );
      }

      fmt::print( "Iterations:         {}\n", m_num_iterations );
      fmt::print( "Function evals:     {}\n", m_num_function_evals );
      fmt::print( "Jacobian evals:     {}\n", m_num_jacobian_evals );
      if ( m_use_line_search ) { fmt::print( "Line searches:      {}\n", m_num_line_searches ); }
      fmt::print( "Strategy:           {}\n", strategy_name( m_strategy ) );
      fmt::print( "Relaxation:         {:.3f}\n", m_relaxation );
      fmt::print( "Initial residual:   {:.6e}\n", initial_norm );
      fmt::print( "Final residual:     {:.6e}\n", final_norm );

      if ( initial_norm > EPSILON )
      {
        real_type reduction = initial_norm / final_norm;
        fmt::print( "Residual reduction: {:.3e} ({:.1f}%)\n", reduction, ( 1.0 - final_norm / initial_norm ) * 100.0 );
      }

      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n\n" );
    }

    //! Stampa riepilogo per il metodo block
    void
    print_summary_block( real_type initial_norm, real_type final_norm ) const
    {
      if ( m_verbose_level == 0 ) return;

      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n" );
      fmt::print( "Block Kaczmarz - Summary\n" );
      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n" );

      if ( m_converged ) { fmt::print( fmt::fg( fmt::color::green ), "Status:             ✓ CONVERGED\n" ); }
      else
      {
        fmt::print( fmt::fg( fmt::color::red ), "Status:             ✗ NOT CONVERGED\n" );
      }

      fmt::print( "Iterations:         {}\n", m_num_iterations );
      fmt::print( "Function evals:     {}\n", m_num_function_evals );
      fmt::print( "Jacobian evals:     {}\n", m_num_jacobian_evals );
      fmt::print( "Strategy:           {}\n", strategy_name( m_strategy ) );
      fmt::print( "Relaxation:         {:.3f}\n", m_relaxation );
      fmt::print( "Initial residual:   {:.6e}\n", initial_norm );
      fmt::print( "Final residual:     {:.6e}\n", final_norm );

      if ( initial_norm > EPSILON )
      {
        real_type reduction = initial_norm / final_norm;
        fmt::print( "Residual reduction: {:.3e} ({:.1f}%)\n", reduction, ( 1.0 - final_norm / initial_norm ) * 100.0 );
      }

      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n\n" );
    }
  };

}  // namespace Utils

#endif
