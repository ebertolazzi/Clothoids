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
// file: Utils_IncrementalLevenbergMarquardt.hh
//

#pragma once

#ifndef UTILS_INCREMENTAL_LEVENBERG_MARQUARDT_dot_HH
#define UTILS_INCREMENTAL_LEVENBERG_MARQUARDT_dot_HH

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
#include "Utils_pseudoinverse.hh"

namespace Utils
{

  class IncrementalLevenbergMarquardt
  {
  public:
    using integer      = Eigen::Index;
    using real_type    = double;
    using Vector       = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    using Matrix       = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<real_type>;

    enum SelectionStrategy
    {
      CYCLIC,           // Selezione ciclica dei blocchi
      RANDOM_UNIFORM,   // Selezione casuale uniforme
      RANDOM_WEIGHTED,  // Selezione casuale pesata per residuo
      GREEDY,           // Selezione del blocco con residuo massimo
      RANDOM_PARTITION  // Partizione casuale ad ogni epoca
    };

  private:
    // Parametri del metodo
    real_type         m_tolerance;
    real_type         m_relative_tolerance;  // Tolleranza relativa
    integer           m_max_iterations;
    integer           m_max_function_evals;
    real_type         m_lambda;          // Parametro di regolarizzazione di Levenberg-Marquardt
    real_type         m_lambda_factor;   // Fattore per adattare lambda (tipicamente 10)
    real_type         m_lambda_min;      // Valore minimo per lambda
    real_type         m_lambda_max;      // Valore massimo per lambda
    real_type         m_good_reduction;  // Soglia per riduzione buona (aumenta trust region)
    real_type         m_bad_reduction;   // Soglia per riduzione cattiva (riduce trust region)
    SelectionStrategy m_strategy;
    integer           m_block_size;  // Dimensione del blocco (numero di equazioni per
                                     // iterazione)
    bool      m_adaptive_lambda;     // Adatta lambda automaticamente
    bool      m_use_line_search;     // Usa line search per ottimizzare il passo
    real_type m_line_search_beta;    // Parametro di riduzione per line search
    real_type m_line_search_c1;      // Condizione di Armijo

    // Controllo output
    integer m_verbose_level;    // 0=silent, 1=summary, 2=detailed, 3=debug
    integer m_print_frequency;  // Stampa ogni N iterazioni

    // Statistiche
    integer   m_num_iterations;
    integer   m_num_function_evals;
    integer   m_num_jacobian_evals;
    integer   m_num_line_searches;
    integer   m_num_lambda_updates;
    bool      m_converged;
    real_type m_final_residual;
    real_type m_initial_residual;

    // Generatore casuale
    mutable std::mt19937                              m_random_engine;
    mutable std::uniform_real_distribution<real_type> m_uniform_dist;

    // Cache
    mutable std::vector<real_type> m_residuals;
    mutable std::vector<integer>   m_permutation;

    // Costanti
    static constexpr real_type EPSILON = 1e-16;

  public:
    IncrementalLevenbergMarquardt()
      : m_tolerance( 1e-8 )
      , m_relative_tolerance( 1e-8 )
      , m_max_iterations( 1000 )
      , m_max_function_evals( 10000 )
      , m_lambda( 0.1 )         // Valore iniziale più piccolo
      , m_lambda_factor( 2.0 )  // Fattore più piccolo per adattamento graduale
      , m_lambda_min( 1e-12 )
      , m_lambda_max( 1e6 )       // Massimo più piccolo per evitare overflow
      , m_good_reduction( 0.25 )  // Soglie più conservative
      , m_bad_reduction( 0.1 )
      , m_strategy( RANDOM_WEIGHTED )  // Strategia pesata per residuo
      , m_block_size( 10 )             // Blocco più grande per stabilità
      , m_adaptive_lambda( true )
      , m_use_line_search( true )  // Line search abilitata
      , m_line_search_beta( 0.5 )
      , m_line_search_c1( 1e-4 )
      , m_verbose_level( 1 )
      , m_print_frequency( 50 )
      , m_num_iterations( 0 )
      , m_num_function_evals( 0 )
      , m_num_jacobian_evals( 0 )
      , m_num_line_searches( 0 )
      , m_num_lambda_updates( 0 )
      , m_converged( false )
      , m_final_residual( 0.0 )
      , m_initial_residual( 0.0 )
      , m_uniform_dist( 0.0, 1.0 )
    {
      // Inizializza generatore random con seed temporale
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
    set_lambda( real_type lambda )
    {
      UTILS_ASSERT( lambda >= 0, "Lambda must be non-negative" );
      m_lambda = lambda;
    }

    void
    set_lambda_factor( real_type factor )
    {
      UTILS_ASSERT( factor > 1, "Lambda factor must be > 1" );
      m_lambda_factor = factor;
    }

    void
    set_lambda_min( real_type min_val )
    {
      UTILS_ASSERT( min_val >= 0, "Lambda min must be non-negative" );
      m_lambda_min = min_val;
    }

    void
    set_lambda_max( real_type max_val )
    {
      UTILS_ASSERT( max_val > 0, "Lambda max must be positive" );
      m_lambda_max = max_val;
    }

    void
    set_good_reduction( real_type good )
    {
      UTILS_ASSERT( good > 0 && good < 1, "Good reduction must be in (0,1)" );
      m_good_reduction = good;
    }

    void
    set_bad_reduction( real_type bad )
    {
      UTILS_ASSERT( bad > 0 && bad < 1, "Bad reduction must be in (0,1)" );
      m_bad_reduction = bad;
    }

    void
    set_strategy( SelectionStrategy s )
    {
      m_strategy = s;
    }

    void
    set_block_size( integer size )
    {
      UTILS_ASSERT( size > 0, "Block size must be positive" );
      m_block_size = size;
    }

    void
    set_adaptive_lambda( bool adaptive )
    {
      m_adaptive_lambda = adaptive;
    }

    void
    set_use_line_search( bool use )
    {
      m_use_line_search = use;
    }

    void
    set_line_search_beta( real_type beta )
    {
      UTILS_ASSERT( beta > 0 && beta < 1, "Beta must be in (0,1)" );
      m_line_search_beta = beta;
    }

    void
    set_line_search_c1( real_type c1 )
    {
      UTILS_ASSERT( c1 > 0 && c1 < 0.5, "c1 must be in (0,0.5)" );
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
      UTILS_ASSERT( level >= 0 && level <= 3, "Verbose level must be 0, 1, 2, or 3" );
      m_verbose_level = level;
    }

    void
    set_print_frequency( integer freq )
    {
      UTILS_ASSERT( freq > 0, "Print frequency must be positive" );
      m_print_frequency = freq;
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
    integer
    get_num_lambda_updates() const
    {
      return m_num_lambda_updates;
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
    get_initial_residual() const
    {
      return m_initial_residual;
    }
    real_type
    get_lambda() const
    {
      return m_lambda;
    }
    real_type
    get_tolerance() const
    {
      return m_tolerance;
    }
    integer
    get_block_size() const
    {
      return m_block_size;
    }
    SelectionStrategy
    get_strategy() const
    {
      return m_strategy;
    }

    // ========================================================================
    // METODI PRIVATI
    // ========================================================================

  private:
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
        case RANDOM_PARTITION:
          return "Random Partition";
        default:
          return "Unknown";
      }
    }

    //! Seleziona un blocco di equazioni in base alla strategia scelta
    std::vector<integer>
    select_block( const Vector & f, integer epoch = -1 )
    {
      integer              m = f.size();
      std::vector<integer> block_indices;

      // Assicurati che block_size non sia più grande di m
      integer actual_block_size = std::min( m_block_size, m );

      switch ( m_strategy )
      {
        case CYCLIC:
        {
          integer start = ( m_num_iterations * actual_block_size ) % m;
          for ( integer i = 0; i < actual_block_size; ++i ) { block_indices.push_back( ( start + i ) % m ); }
          break;
        }

        case RANDOM_UNIFORM:
        {
          std::uniform_int_distribution<integer> dist( 0, m - 1 );
          for ( integer i = 0; i < actual_block_size; ++i ) { block_indices.push_back( dist( m_random_engine ) ); }
          break;
        }

        case RANDOM_WEIGHTED:
        {
          // Calcola residui assoluti
          m_residuals.resize( static_cast<size_t>( m ) );
          real_type sum_residuals_sq = 0.0;
          for ( integer i = 0; i < m; ++i )
          {
            real_type fi                          = std::abs( f( i ) );
            m_residuals[static_cast<size_t>( i )] = fi * fi;
            sum_residuals_sq += m_residuals[static_cast<size_t>( i )];
          }

          if ( sum_residuals_sq < EPSILON )
          {
            std::uniform_int_distribution<integer> dist( 0, m - 1 );
            for ( integer i = 0; i < actual_block_size; ++i ) { block_indices.push_back( dist( m_random_engine ) ); }
          }
          else
          {
            // Normalizza per creare distribuzione di probabilità
            std::vector<real_type> probabilities( static_cast<size_t>( m ) );
            for ( integer i = 0; i < m; ++i )
            {
              probabilities[static_cast<size_t>( i )] = m_residuals[static_cast<size_t>( i )] / sum_residuals_sq;
            }

            std::discrete_distribution<integer> dist( probabilities.begin(), probabilities.end() );
            for ( integer i = 0; i < actual_block_size; ++i ) { block_indices.push_back( dist( m_random_engine ) ); }
          }
          break;
        }

        case GREEDY:
        {
          // Trova le equazioni con residuo massimo
          std::vector<std::pair<real_type, integer>> residuals_with_idx;
          residuals_with_idx.reserve( static_cast<size_t>( m ) );

          for ( integer i = 0; i < m; ++i ) { residuals_with_idx.emplace_back( std::abs( f( i ) ), i ); }

          // Ordina per residuo decrescente
          std::partial_sort( residuals_with_idx.begin(), residuals_with_idx.begin() + actual_block_size,
                             residuals_with_idx.end(),
                             []( const auto & a, const auto & b ) { return a.first > b.first; } );

          for ( integer i = 0; i < actual_block_size; ++i )
          {
            block_indices.push_back( residuals_with_idx[static_cast<size_t>( i )].second );
          }
          break;
        }

        case RANDOM_PARTITION:
        {
          // Crea una permutazione casuale delle equazioni
          if ( epoch >= 0 && m_permutation.size() != static_cast<size_t>( m ) )
          {
            m_permutation.resize( static_cast<size_t>( m ) );
            std::iota( m_permutation.begin(), m_permutation.end(), 0 );
            std::shuffle( m_permutation.begin(), m_permutation.end(), m_random_engine );
          }

          integer start = ( m_num_iterations * actual_block_size ) % m;
          for ( integer i = 0; i < actual_block_size; ++i )
          {
            block_indices.push_back( m_permutation[static_cast<size_t>( ( start + i ) % m )] );
          }
          break;
        }

        default:
          std::uniform_int_distribution<integer> dist( 0, m - 1 );
          for ( integer i = 0; i < actual_block_size; ++i ) { block_indices.push_back( dist( m_random_engine ) ); }
      }

      return block_indices;
    }

    //! Estrae il blocco di Jacobiano e residui
    void
    extract_block( NonlinearSystem &            system,
                   const Vector &               x,
                   const Vector &               f_full,
                   const std::vector<integer> & block_indices,
                   Matrix &                     J_block,
                   Vector &                     f_block )
    {
      integer m          = system.num_equations();
      integer block_size = block_indices.size();

      SparseMatrix J_full( m, m );  // Sistema square
      system.jacobian( x, J_full );
      m_num_jacobian_evals++;

      // Estrai il blocco di Jacobiano
      J_block.resize( block_size, m );
      f_block.resize( block_size );

      for ( integer i = 0; i < block_size; ++i )
      {
        integer row_idx = block_indices[static_cast<size_t>( i )];
        f_block( i )    = f_full( row_idx );

        for ( SparseMatrix::InnerIterator it( J_full, row_idx ); it; ++it ) { J_block( i, it.col() ) = it.value(); }
      }
    }

    //! Calcola il passo LM incrementale usando pseudo-inversa regolarizzata
    Vector
    compute_step( const Matrix & J_block, const Vector & f_block, real_type lambda )
    {
      // Aggiungi regolarizzazione diagonale piccola per stabilità
      Matrix A = J_block.transpose() * J_block;
      A.diagonal().array() += lambda;

      // Risolvi (J^T J + lambda*I) dx = -J^T f
      Eigen::LDLT<Matrix> solver;
      solver.compute( A );

      if ( solver.info() != Eigen::Success )
      {
        // Fallback: usa SVD
        Eigen::JacobiSVD<Matrix> svd( J_block, Eigen::ComputeThinU | Eigen::ComputeThinV );
        Vector                   singular_values = svd.singularValues();
        Vector                   inv_singular = singular_values.array() / ( singular_values.array().square() + lambda );
        return -svd.matrixV() * inv_singular.asDiagonal() * svd.matrixU().transpose() * f_block;
      }

      Vector rhs = -J_block.transpose() * f_block;
      return solver.solve( rhs );
    }

    //! Calcola il rapporto di riduzione (rho) per LM - VERSIONE CORRETTA
    real_type
    compute_reduction_ratio( NonlinearSystem & system,
                             const Vector &    x,
                             real_type         norm_f_old,
                             const Vector &    dx,
                             const Matrix &    J_block,
                             const Vector &    f_block,
                             real_type         lambda )
    {
      // Predicted reduction: phi(0) - phi(dx)
      // phi(dx) = ||f_block + J_block*dx||^2 + lambda*||dx||^2
      Vector    Jdx                 = J_block * dx;
      real_type predicted_f_norm_sq = ( f_block + Jdx ).squaredNorm();
      real_type phi_dx              = predicted_f_norm_sq + lambda * dx.squaredNorm();
      real_type phi_0               = f_block.squaredNorm();
      real_type predicted_reduction = phi_0 - phi_dx;

      // Actual reduction: ||f_old||^2 - ||f_new||^2
      Vector x_new = x + dx;
      Vector f_new( system.num_equations() );
      system.evaluate( x_new, f_new );
      m_num_function_evals++;

      real_type norm_f_new       = f_new.norm();
      real_type actual_reduction = norm_f_old * norm_f_old - norm_f_new * norm_f_new;

      // Evita divisione per zero e valori problematici
      if ( std::abs( predicted_reduction ) < 1e-12 ) { return ( actual_reduction > 0 ) ? 1.0 : 0.0; }

      real_type rho = actual_reduction / predicted_reduction;

      // Limita rho per stabilità numerica
      if ( rho > 2.0 ) rho = 2.0;
      if ( rho < -1.0 ) rho = -1.0;

      return rho;
    }

    //! Line search di Armijo per LM
    real_type
    armijo_line_search( NonlinearSystem & system,
                        const Vector &    x,
                        real_type         norm_f_old,
                        const Vector &    dx,
                        const Matrix &    J_block,
                        const Vector &    f_block )
    {
      integer   n     = x.size();
      real_type alpha = 1.0;
      Vector    x_new( n );
      Vector    f_new( system.num_equations() );

      // Calcola la derivata direzionale: 2 * f_block^T * J_block * dx
      Vector    Jdx                    = J_block * dx;
      real_type directional_derivative = 2.0 * f_block.dot( Jdx );

      // Fino a 20 tentativi di riduzione del passo
      for ( int ls_iter = 0; ls_iter < 20; ++ls_iter )
      {
        // Nuovo punto
        x_new = x + alpha * dx;

        // Valuta funzione
        try
        {
          system.check_if_admissible( x_new );
          system.evaluate( x_new, f_new );
          m_num_function_evals++;

          real_type norm_f_new         = f_new.norm();
          real_type actual_reduction   = norm_f_old * norm_f_old - norm_f_new * norm_f_new;
          real_type required_reduction = m_line_search_c1 * alpha * directional_derivative;

          if ( actual_reduction >= required_reduction )
          {
            m_num_line_searches++;
            return alpha;  // Step accettato
          }
        }
        catch ( ... )
        {
          // Punto non ammissibile
        }

        // Riduci il passo
        alpha *= m_line_search_beta;
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

    //! Stampa informazioni sull'iterazione
    void
    print_iteration_info( integer   iter,
                          real_type norm_f,
                          real_type initial_norm,
                          real_type lambda,
                          real_type rho,
                          real_type dx_norm ) const
    {
      if ( m_verbose_level < 2 ) return;
      if ( iter % m_print_frequency != 0 ) return;

      real_type reduction = initial_norm > EPSILON ? initial_norm / norm_f : 1.0;

      fmt::print( fmt::fg( fmt::color::light_blue ), "[{:5}] ‖f‖ = {:.2e}", fmt::format( "{}", iter + 1 ), norm_f );
      fmt::print( ", λ = {:.2e}, ρ = {:.3f}", lambda, rho );
      fmt::print( ", ‖dx‖ = {:.2e}", dx_norm );
      fmt::print( ", reduction = {:.3e}\n", reduction );
    }

    //! Stampa riepilogo finale
    void
    print_summary( real_type initial_norm, real_type final_norm ) const
    {
      if ( m_verbose_level == 0 ) return;

      fmt::print( "\n" );
      fmt::print( fmt::fg( fmt::color::light_blue ),
                  "════════════════════════════════════════════════════════════"
                  "═══════\n" );
      fmt::print( "Incremental Levenberg-Marquardt - Summary\n" );
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
      fmt::print( "Lambda updates:     {}\n", m_num_lambda_updates );
      fmt::print( "Strategy:           {}\n", strategy_name( m_strategy ) );
      fmt::print( "Block size:         {}\n", m_block_size );
      fmt::print( "Final lambda:       {:.2e}\n", m_lambda );
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

  public:
    // ========================================================================
    // METODO PRINCIPALE DI RISOLUZIONE - VERSIONE MIGLIORATA
    // ========================================================================

    bool
    solve( NonlinearSystem & system, Vector & x )
    {
      integer m = system.num_equations();
      Vector  f( m );

      // Reset statistiche
      m_num_iterations     = 0;
      m_num_function_evals = 0;
      m_num_jacobian_evals = 0;
      m_num_line_searches  = 0;
      m_num_lambda_updates = 0;
      m_converged          = false;

      // Valuta funzione iniziale
      system.evaluate( x, f );
      m_num_function_evals++;
      real_type norm_f      = f.norm();
      m_initial_residual    = norm_f;
      real_type best_norm_f = norm_f;
      Vector    best_x      = x;

      if ( m_verbose_level > 0 )
      {
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
        fmt::print( fmt::fg( fmt::color::light_blue ), "Starting Incremental Levenberg-Marquardt\n" );
        fmt::print( "Initial residual:   {:.2e}\n", norm_f );
        fmt::print( "Block size:         {}\n", m_block_size );
        fmt::print( "Strategy:           {}\n", strategy_name( m_strategy ) );
        fmt::print( "Initial lambda:     {:.2e}\n", m_lambda );
        fmt::print( fmt::fg( fmt::color::light_blue ),
                    "══════════════════════════════════════════════════════════"
                    "═════════\n" );
      }

      // Buffer per permutazione se necessario
      if ( m_strategy == RANDOM_PARTITION )
      {
        m_permutation.resize( static_cast<size_t>( m ) );
        std::iota( m_permutation.begin(), m_permutation.end(), 0 );
        std::shuffle( m_permutation.begin(), m_permutation.end(), m_random_engine );
      }

      integer   epoch             = 0;
      real_type last_good_norm_f  = norm_f;
      integer   no_progress_count = 0;

      // Iterazioni principali
      for ( m_num_iterations = 0; m_num_iterations < m_max_iterations && m_num_function_evals < m_max_function_evals;
            ++m_num_iterations )
      {
        // Check convergenza
        if ( check_convergence( norm_f, m_initial_residual ) )
        {
          m_converged      = true;
          m_final_residual = norm_f;

          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::green ), "✓ Converged below tolerance ({:.2e})\n", m_tolerance );
            print_summary( m_initial_residual, norm_f );
          }
          return true;
        }

        // Nuova epoca per RANDOM_PARTITION
        if ( m_strategy == RANDOM_PARTITION && ( m_num_iterations % ( m / m_block_size + 1 ) ) == 0 )
        {
          std::shuffle( m_permutation.begin(), m_permutation.end(), m_random_engine );
          epoch++;
        }

        // Seleziona blocco di equazioni
        std::vector<integer> block_indices = select_block( f, epoch );

        // Estrai Jacobiano e residui per il blocco
        Matrix J_block;
        Vector f_block;
        extract_block( system, x, f, block_indices, J_block, f_block );

        // Calcola il passo LM incrementale
        Vector dx = compute_step( J_block, f_block, m_lambda );

        // Controlla se il passo è troppo piccolo
        real_type dx_norm = dx.norm();
        if ( dx_norm < 1e-12 )
        {
          // Incrementa lambda e continua
          m_lambda = std::min( m_lambda * m_lambda_factor, m_lambda_max );
          m_num_lambda_updates++;

          if ( m_verbose_level >= 3 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "  ⚠ Small step (‖dx‖={:.2e}), increasing λ to {:.2e}\n",
                        dx_norm, m_lambda );
          }
          continue;
        }

        // Line search
        real_type alpha = 1.0;
        if ( m_use_line_search )
        {
          if ( m_verbose_level >= 3 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "  Line search: α₀={:.3f}\n", alpha );
          }

          alpha = armijo_line_search( system, x, norm_f, dx, J_block, f_block );
          if ( alpha <= 0.0 )
          {
            if ( m_verbose_level >= 2 ) { fmt::print( fmt::fg( fmt::color::red ), "  ✗ Line search failed\n" ); }
            alpha = 0.1;  // Fallback a passo piccolo
          }
          else if ( m_verbose_level >= 3 )
          {
            fmt::print( fmt::fg( fmt::color::green ), "  ✓ Line search: α={:.3f}\n", alpha );
          }
        }

        // Calcola il nuovo punto candidato
        Vector x_candidate = x + alpha * dx;
        Vector f_candidate( m );
        system.evaluate( x_candidate, f_candidate );
        m_num_function_evals++;
        real_type norm_f_candidate = f_candidate.norm();

        // Calcola il rapporto di riduzione
        real_type rho = compute_reduction_ratio( system, x, norm_f, alpha * dx, J_block, f_block, m_lambda );

        // Stampa info dettagliate
        print_iteration_info( m_num_iterations, norm_f, m_initial_residual, m_lambda, rho, dx_norm );

        // Adatta lambda in base al rapporto di riduzione
        bool step_accepted = false;
        if ( m_adaptive_lambda )
        {
          if ( rho > m_good_reduction )
          {
            // Riduzione buona
            step_accepted        = true;
            real_type old_lambda = m_lambda;
            m_lambda             = std::max( m_lambda / m_lambda_factor, m_lambda_min );
            m_num_lambda_updates++;
            no_progress_count = 0;

            if ( m_verbose_level >= 3 )
            {
              fmt::print( fmt::fg( fmt::color::green ),
                          "  ✓ Good reduction (ρ={:.3f}), decreasing λ: {:.2e} "
                          "→ {:.2e}\n",
                          rho, old_lambda, m_lambda );
            }
          }
          else if ( rho > m_bad_reduction )
          {
            // Riduzione moderata
            step_accepted = true;
            // Mantieni lambda
            no_progress_count++;

            if ( m_verbose_level >= 3 )
            {
              fmt::print( fmt::fg( fmt::color::yellow ), "  ~ Moderate reduction (ρ={:.3f}), keeping λ={:.2e}\n", rho,
                          m_lambda );
            }
          }
          else
          {
            // Riduzione cattiva
            step_accepted        = false;
            real_type old_lambda = m_lambda;
            m_lambda             = std::min( m_lambda * m_lambda_factor, m_lambda_max );
            m_num_lambda_updates++;
            no_progress_count++;

            if ( m_verbose_level >= 2 )
            {
              fmt::print( fmt::fg( fmt::color::red ),
                          "  ✗ Bad reduction (ρ={:.3f}), increasing λ: {:.2e} "
                          "→ {:.2e}\n",
                          rho, old_lambda, m_lambda );
            }
          }
        }
        else
        {
          step_accepted = ( norm_f_candidate < norm_f );
        }

        // Accetta o rifiuta il passo
        if ( step_accepted )
        {
          x      = x_candidate;
          f      = f_candidate;
          norm_f = norm_f_candidate;

          // Mantieni la migliore soluzione trovata
          if ( norm_f < best_norm_f )
          {
            best_norm_f = norm_f;
            best_x      = x;
          }
        }
        else
        {
          if ( m_verbose_level >= 2 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "  ⚠ Step rejected, trying with larger λ\n" );
          }
        }

        // Criterio di arresto per progresso insufficiente
        if ( no_progress_count > 20 )
        {
          if ( norm_f / last_good_norm_f > 0.99 )
          {
            if ( m_verbose_level > 0 )
            {
              fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Stopping: insufficient progress for 20 iterations\n" );
            }
            break;
          }
          last_good_norm_f  = norm_f;
          no_progress_count = 0;
        }

        // Criterio di arresto per lambda troppo grande
        if ( m_lambda > m_lambda_max * 0.5 )
        {
          if ( m_verbose_level > 0 )
          {
            fmt::print( fmt::fg( fmt::color::yellow ), "⚠ Lambda too large ({:.2e}), resetting to initial\n",
                        m_lambda );
          }
          m_lambda          = 0.1;
          no_progress_count = 0;
        }
      }

      // Ripristina la migliore soluzione trovata
      if ( best_norm_f < norm_f )
      {
        x      = best_x;
        norm_f = best_norm_f;
        if ( m_verbose_level >= 2 )
        {
          fmt::print( fmt::fg( fmt::color::green ), "  Restored best solution with residual {:.2e}\n", best_norm_f );
        }
      }

      m_final_residual = norm_f;
      m_converged      = ( norm_f < m_tolerance );

      if ( m_verbose_level > 0 )
      {
        if ( m_converged )
        {
          fmt::print( fmt::fg( fmt::color::green ),
                      "✓ Finished: {} iterations, {} function evals, final "
                      "residual: {:.2e}\n",
                      m_num_iterations, m_num_function_evals, norm_f );
        }
        else
        {
          fmt::print( fmt::fg( fmt::color::yellow ),
                      "⚠ Finished: {} iterations, {} function evals, final "
                      "residual: {:.2e}\n",
                      m_num_iterations, m_num_function_evals, norm_f );
        }
        print_summary( m_initial_residual, norm_f );
      }

      return m_converged;
    }
  };

}  // namespace Utils

#endif
