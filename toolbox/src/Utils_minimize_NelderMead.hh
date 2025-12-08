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
// file: Utils_NelderMead.hh
//

#pragma once

#ifndef UTILS_NELDER_MEAD_dot_HH
#define UTILS_NELDER_MEAD_dot_HH

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

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::string;
  using std::vector;

  // ===========================================================================
  // ENUMS & RESULT STRUCT
  // ===========================================================================

  namespace NelderMead
  {

    // Utilizzo di Eigen::Matrix per vettori dinamici - ottimizzato per
    // operazioni vettoriali
    template <typename Scalar>
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Utilizzo di Eigen::Matrix per matrici dinamiche - efficiente per
    // operazioni lineari
    template <typename Scalar>
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    template <typename Scalar>
    using Callback = std::function<Scalar( Vector<Scalar> const & )>;

    // Color definitions for consistent styling
    namespace PrintColors
    {
      constexpr auto HEADER    = fmt::fg( fmt::color::light_blue );
      constexpr auto SUCCESS   = fmt::fg( fmt::color::green );
      constexpr auto WARNING   = fmt::fg( fmt::color::yellow );
      constexpr auto ERROR     = fmt::fg( fmt::color::red );
      constexpr auto INFO      = fmt::fg( fmt::color::cyan );
      constexpr auto ITERATION = fmt::fg( fmt::color::white );
    }  // namespace PrintColors

    enum class Status
    {
      RUNNING,
      CONVERGED,            // Function tolerance satisfied
      CONVERGED_TOL_X,      // Simplex size tolerance satisfied
      MAX_ITERATIONS,       // Max iterations reached
      MAX_FUN_EVALUATIONS,  // Max evaluations reached
      FAIL_NAN,             // NaN detected
      FAILED
    };

    inline string
    status_to_string( Status s )
    {
      switch ( s )
      {
        case Status::RUNNING:
          return "RUNNING";
        case Status::CONVERGED:
          return "CONVERGED";
        case Status::CONVERGED_TOL_X:
          return "CONVERGED_X";
        case Status::MAX_ITERATIONS:
          return "MAX_ITER";
        case Status::MAX_FUN_EVALUATIONS:
          return "MAX_EVAL";
        case Status::FAIL_NAN:
          return "FAIL_NAN";
        case Status::FAILED:
          return "FAILED";
        default:
          return "UNKNOWN";
      }
    }

    template <typename Scalar>
    struct Result
    {
      Vector<Scalar> solution;
      Scalar         final_function_value{ 0 };
      Scalar         initial_function_value{ 0 };
      Status         status{ Status::FAILED };

      // Statistics
      size_t outer_iterations{ 0 };  // Main loop cycles
      size_t inner_iterations{ 0 };  // Total simplex steps (sum of all sub-runs)
      size_t total_iterations{ 0 };  // Combined count

      size_t outer_evaluations{ 0 };
      size_t inner_evaluations{ 0 };
      size_t total_evaluations{ 0 };
    };

    // Helper for vector formatting
    template <typename Scalar>
    inline string
    format_vector( Vector<Scalar> const & v, size_t max_size = 10 )
    {
      string tmp{ "[" };
      size_t v_size = v.size();
      if ( v_size <= max_size )
      {
        for ( size_t i = 0; i < v_size; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
      }
      else
      {
        for ( size_t i{ 0 }; i < max_size - 3; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
        tmp.pop_back();
        tmp += "...";
        for ( size_t i{ v_size - 3 }; i < v_size; ++i ) tmp += fmt::format( "{:.4f}, ", v( i ) );
      }
      tmp.pop_back();
      tmp.pop_back();
      tmp += "]";
      return tmp;
    }

  }  // namespace NelderMead

  // ===========================================================================
  // CLASS: NelderMead_classic (Inner Solver)
  // ===========================================================================

  /**
   * @class NelderMead_classic
   * @brief Implementation of the classic Nelder-Mead simplex optimization
   * algorithm
   *
   * @tparam Scalar Numeric type for computations (default: double)
   *
   * This class implements the Nelder-Mead simplex algorithm for unconstrained
   * and bound-constrained optimization. It includes advanced features like
   * adaptive parameters, restart mechanisms, and robust convergence checking.
   *
   * EIGEN3 OPTIMIZATIONS:
   * - Uses Eigen::Matrix for efficient vector/matrix operations
   * - Leverages Eigen's expression templates for lazy evaluation
   * - Utilizes Eigen's efficient memory management and SIMD operations
   * - Employs Eigen's numerical linear algebra routines for stability
   */
  template <typename Scalar = double>
  class NelderMead_classic
  {
  public:
    // Eigen3 type aliases for efficient linear algebra
    using Vector   = NelderMead::Vector<Scalar>;  // Eigen::Matrix column vector
    using Matrix   = NelderMead::Matrix<Scalar>;  // Eigen::Matrix dynamic matrix
    using Callback = NelderMead::Callback<Scalar>;

    /**
     * @brief Optimization status codes
     */
    enum class Status
    {
      CONVERGED,            ///< Successfully converged to optimum
      MAX_ITERATIONS,       ///< Reached maximum iteration limit
      MAX_FUN_EVALUATIONS,  ///< Reached maximum function evaluation limit
      SIMPLEX_TOO_SMALL,    ///< Simplex became too small to continue
      STAGNATED,            ///< Algorithm stagnated (no improvement)
      FAILED                ///< Optimization failed
    };

    /**
     * @brief Configuration options for Nelder-Mead algorithm
     */
    struct Options
    {
      // Global budget (including restarts)
      size_t max_iterations{ 10000 };            ///< Maximum total iterations
      size_t max_function_evaluations{ 50000 };  ///< Maximum function evaluations
      Scalar tolerance{ 1e-8 };                  ///< Convergence tolerance
      Scalar stagnation_tolerance{ 1e-8 };       ///< Stagnation detection tolerance
      Scalar simplex_tolerance{ 1e-10 };         ///< Minimum simplex size tolerance

      // Standard Nelder-Mead parameters
      Scalar rho{ 1.0 };    ///< Reflection coefficient
      Scalar chi{ 2.0 };    ///< Expansion coefficient
      Scalar gamma{ 0.5 };  ///< Contraction coefficient
      Scalar sigma{ 0.5 };  ///< Shrink coefficient

      Scalar initial_step{ 0.1 };  ///< Initial step size for simplex construction

      bool   adaptive_parameters{ true };  ///< Enable adaptive parameter adjustment
      bool   verbose{ true };              ///< Enable verbose output
      size_t progress_frequency{ 100 };    ///< Progress reporting frequency

      // Restart mechanism
      bool   enable_restart{ true };          ///< Enable restart strategy
      size_t max_restarts{ 5 };               ///< Maximum number of restarts
      size_t stagnation_threshold{ 30 };      ///< Iterations before stagnation detection
      bool   use_relative_tolerance{ true };  ///< Use relative convergence criteria

      Scalar restart_perturbation_ratio{ 0.25 };  ///< Restart perturbation scale
      bool   track_best_point{ true };            ///< Track best point across restarts

      // Enhanced convergence parameters
      Scalar min_step_size{ 1e-10 };          ///< Minimum step size to avoid degeneracy
      bool   use_robust_convergence{ true };  ///< Use robust convergence criteria
      Scalar convergence_relaxation{ 10.0 };  ///< Convergence criterion relaxation factor

      // Restart condition thresholds
      Scalar restart_relative_improvement_threshold{ 0.05 };   ///< 5% min improvement for restart
      Scalar restart_progress_per_eval_threshold{ 1e-6 };      ///< Progress per evaluation threshold
      Scalar restart_simplex_geometry_threshold{ 0.3 };        ///< Improvement threshold for small simplex
      Scalar restart_shrink_count_threshold{ 8 };              ///< Shrink operations before restart
      Scalar restart_after_shrink_improvement{ 0.02 };         ///< Improvement after shrink threshold
      Scalar restart_expected_progress_ratio{ 0.05 };          ///< Expected progress ratio for high dimensions
      Scalar restart_quality_metric_threshold{ 1e-4 };         ///< Quality metric threshold
      Scalar restart_degenerate_improvement_threshold{ 0.1 };  ///< Improvement threshold for degenerate simplex
      Scalar restart_improvement_ratio{ 0.95 };                ///< Improvement ratio to accept restart
      Scalar restart_absolute_improvement_threshold{ 0.1 };    ///< Absolute improvement threshold

      // Geometry factors for restart conditions
      Scalar restart_simplex_diameter_factor1{ 50.0 };  ///< Diameter factor for condition 3
      Scalar restart_simplex_diameter_factor2{ 20.0 };  ///< Diameter factor for condition 7
      Scalar restart_std_dev_factor{ 10.0 };            ///< Standard deviation factor

      size_t verbosity_level{ 1 };            // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      size_t inner_progress_frequency{ 10 };  // Frequency for level 2
    };

    /**
     * @brief Inner result structure for compatibility with BlockCoordinate
     */
    struct InnerResult
    {
      Vector solution;
      Scalar final_function_value{ 0 };
      Scalar initial_function_value{ 0 };
      Status status{ Status::FAILED };

      // Statistics for single run
      size_t iterations{ 0 };            ///< Iterations in this run
      size_t function_evaluations{ 0 };  ///< Function evaluations in this run
      Scalar simplex_volume{ 0 };
      Scalar simplex_diameter{ 0 };
      size_t restarts_performed{ 0 };
      size_t shrink_operations{ 0 };

      // Compatibility fields for BlockCoordinate
      size_t inner_iterations{ 0 };   ///< Same as iterations
      size_t inner_evaluations{ 0 };  ///< Same as function_evaluations
    };

  private:
    /**
     * @brief Statistics for simplex analysis
     */
    struct SimplexStats
    {
      Scalar diameter;           ///< Maximum distance between vertices
      Scalar std_dev;            ///< Standard deviation of function values
      Scalar value_range;        ///< Range of function values (max-min)
      Scalar centroid_distance;  ///< Average distance to centroid
    };

    Options m_options;              ///< Algorithm configuration
    Vector  m_lower;                ///< Lower bounds (if used)
    Vector  m_upper;                ///< Upper bounds (if used)
    bool    m_use_bounds{ false };  ///< Whether bounds are active

    Callback const * m_callback{ nullptr };     ///< Objective function callback
    size_t           m_global_iterations{ 0 };  ///< Global iteration counter
    size_t           m_global_evals{ 0 };       ///< Global function evaluation counter

    // EIGEN3: vector of Eigen vectors for efficient simplex storage
    vector<Vector> m_simplex;      ///< Simplex vertices
    vector<Scalar> m_values;       ///< Function values at vertices
    Vector         m_centroid;     ///< Current centroid (excluding worst)
    Vector         m_trial_point;  ///< Trial point for operations
    size_t         m_dim{ 0 };     ///< Problem dimension

    mutable bool           m_simplex_ordered{ false };  ///< Whether simplex is sorted
    mutable vector<size_t> m_sorted_indices;            ///< Indices sorted by function value

    size_t m_stagnation_count{ 0 };                                ///< Consecutive stagnation iterations
    Scalar m_previous_best{ std::numeric_limits<Scalar>::max() };  ///< Previous best value
    size_t m_shrink_count{ 0 };                                    ///< Shrink operation counter

    Vector m_best_point;                                        ///< Best point found (across restarts)
    Scalar m_best_value{ std::numeric_limits<Scalar>::max() };  ///< Best value found

    Scalar m_current_rho, m_current_chi, m_current_gamma,
        m_current_sigma;  ///< Current adaptive parameters

    string m_indent{ "" };

    struct ConvergenceFlags
    {
      bool value_converged;
      bool geometry_converged;
      bool variance_converged;
    };

    ConvergenceFlags
    compute_convergence_flags( Scalar                  best_value,
                               [[maybe_unused]] Scalar worst_value,
                               const SimplexStats &    stats ) const
    {
      ConvergenceFlags flags;

      // 1. Convergenza per i valori della funzione
      if ( m_options.use_relative_tolerance )
      {
        Scalar relative_range = stats.value_range / ( 1.0 + abs( best_value ) );
        flags.value_converged = relative_range < m_options.tolerance;
      }
      else
      {
        flags.value_converged = stats.value_range < m_options.tolerance;
      }

      // 2. Convergenza per la geometria del simplesso
      flags.geometry_converged = stats.diameter < m_options.simplex_tolerance * std::sqrt( m_dim );

      // 3. Convergenza per la varianza
      flags.variance_converged = stats.std_dev < m_options.tolerance * m_options.convergence_relaxation;

      return flags;
    }

    void
    print_inner_iteration_summary( size_t iter_count, Scalar best_value, const SimplexStats & stats ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      bool show_summary = ( m_options.verbosity_level >= 3 ) ||
                          ( iter_count % m_options.inner_progress_frequency == 0 );

      if ( show_summary )
      {
        fmt::print( NelderMead::PrintColors::ITERATION,
                    "{}┌─ Inner Iteration {:5d} "
                    "────────────────────────────────────────────────────────────┐\n"
                    "{}│ Best F:{:<12.6e} Diameter:{:<12.6e} Std Dev:{:<12.6e} "
                    "Volume:{:<12.6e} │\n"
                    "{}"
                    "└─────────────────────────────────────────────────────────────────"
                    "───────────────────┘\n",
                    m_indent, iter_count, m_indent, best_value, stats.diameter, stats.std_dev, compute_volume(),
                    m_indent );
      }
    }

    void
    print_iteration_summary( size_t iter, Scalar best_value, Scalar diameter ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      bool show_detailed = ( m_options.verbosity_level >= 3 ) ||
                           ( m_options.verbosity_level >= 2 && iter % m_options.inner_progress_frequency == 0 );

      if ( !show_detailed ) return;

      auto color = NelderMead::PrintColors::ITERATION;
      fmt::print( color, "{}[{:4d}] F = {:<12.6e} | Diam = {:<12.6e} | Vol = {:<12.6e}\n", m_indent, iter, best_value,
                  diameter, compute_volume() );
    }

    void
    print_inner_operation( string const & operation, Scalar fval, bool improved ) const
    {
      if ( m_options.verbosity_level < 3 ) return;

      auto   color = improved ? NelderMead::PrintColors::SUCCESS : NelderMead::PrintColors::WARNING;
      string icon  = improved ? "✓" : "✗";

      string coords = m_dim <= 5 ? fmt::format( " | x = {}", NelderMead::format_vector<Scalar>( m_trial_point ) ) : "";

      fmt::print( color, "{}  {} {:>12}: F = {:<12.6e}{}\n", m_indent, icon, operation, fval, coords );
    }

    void
    print_convergence_info( Scalar best_value, Scalar worst_value, const SimplexStats & stats ) const
    {
      if ( m_options.verbosity_level < 2 ) return;

      auto flags = compute_convergence_flags( best_value, worst_value, stats );

      bool show_output = ( m_options.verbosity_level >= 3 ) ||
                         ( m_options.verbosity_level >= 2 &&
                           m_global_iterations % m_options.inner_progress_frequency == 0 );

      if ( show_output )
      {
        // Scegli il colore in base allo stato di convergenza
        auto color = ( flags.value_converged && flags.geometry_converged ) ? NelderMead::PrintColors::SUCCESS
                     : ( flags.value_converged || flags.geometry_converged || flags.variance_converged )
                         ? NelderMead::PrintColors::WARNING
                         : NelderMead::PrintColors::INFO;

        // Crea simboli più descrittivi
        string value_symbol = flags.value_converged ? "✓V" : "✗V";
        string geom_symbol  = flags.geometry_converged ? "✓G" : "✗G";
        string var_symbol   = flags.variance_converged ? "✓S" : "✗S";

        fmt::print( color,
                    "{}[Conv{:4d}] {} {} {} | Range={:<10.3e} Diam={:<10.3e} "
                    "StdDev={:<10.3e}",
                    m_indent, m_global_iterations, value_symbol, geom_symbol, var_symbol, stats.value_range,
                    stats.diameter, stats.std_dev );

        // Aggiungi informazioni aggiuntive per verbosità alta
        if ( m_options.verbosity_level >= 3 )
        {
          fmt::print( color, " | CentDist={:<10.3e} Vol={:<10.3e}", stats.centroid_distance, compute_volume() );
        }

        fmt::print( "\n" );

        // Stampa un warning se siamo vicini alla convergenza ma non
        // completamente
        if ( ( flags.value_converged || flags.geometry_converged ) &&
             !( flags.value_converged && flags.geometry_converged ) )
        {
          if ( flags.value_converged && !flags.geometry_converged )
          {
            fmt::print( NelderMead::PrintColors::WARNING,
                        "{}       ⚠ Convergenza valori raggiunta, ma geometria "
                        "simplex ancora ampia\n",
                        m_indent );
          }
          else if ( !flags.value_converged && flags.geometry_converged )
          {
            fmt::print( NelderMead::PrintColors::WARNING,
                        "{}       ⚠ Geometria simplex convergente, ma valori "
                        "funzione ancora dispersi\n",
                        m_indent );
          }
        }

        // Stampa messaggio di successo se completamente convergente
        if ( flags.value_converged && flags.geometry_converged )
        {
          fmt::print( NelderMead::PrintColors::SUCCESS, "{}       🎯 CONVERGENZA RAGGIUNTA!\n", m_indent );
        }
      }
    }

    /**
     * @brief Safely evaluate objective function with bounds checking
     * @param x Point to evaluate
     * @return Function value or large value if out of bounds/invalid
     *
     * EIGEN3: Uses Eigen's array operations for efficient bounds checking
     */
    Scalar
    safe_evaluate( Vector const & x )
    {
      UTILS_ASSERT( m_callback != nullptr, "NelderMead_classic::safe_evaluate(x) Callback not set!" );

      // EIGEN3: Use array operations for efficient bounds checking
      // .array() enables element-wise operations without temporary copies
      if ( m_use_bounds )
      {
        bool out_of_bound = ( x.array() < m_lower.array() ).any() || ( x.array() > m_upper.array() ).any();
        if ( out_of_bound )
        {
          if ( m_options.verbose )
          {
            fmt::print( "{}Warning: Point outside bounds, returning large value\n", m_indent );
          }
          return std::numeric_limits<Scalar>::max();
        }
      }

      Scalar value{ ( *m_callback )( x ) };
      ++m_global_evals;

      if ( !std::isfinite( value ) )
      {
        if ( m_options.verbose )
        {
          fmt::print( "{}Warning: Non-finite function value at x={}\n", m_indent, x.transpose() );
        }
        return std::numeric_limits<Scalar>::max();
      }
      return value;
    }

    /**
     * @brief Project point to feasible region (respect bounds)
     * @param x Point to project (modified in-place)
     *
     * EIGEN3: Uses cwiseMin/cwiseMax for efficient element-wise bounds
     * projection These operations are optimized and avoid explicit loops
     */
    void
    project_point( Vector & x ) const
    {
      if ( m_use_bounds ) x = x.cwiseMax( m_lower ).cwiseMin( m_upper );
    }

    /**
     * @brief Initialize adaptive parameters based on problem dimension
     *
     * EIGEN3: Parameters are tuned for Eigen's efficient vector operations
     */
    void
    initialize_adaptive_parameters()
    {
      if ( m_options.adaptive_parameters && m_dim > 0 )
      {
        Scalar n = static_cast<Scalar>( m_dim );
        // More conservative parameters for high-dimensional problems
        if ( m_dim > 10 )
        {
          m_current_rho   = 1.0;
          m_current_chi   = 1.0 + 1.0 / n;  // Less aggressive expansion
          m_current_gamma = 0.75 - 0.5 / n;
          m_current_sigma = 0.8 - 0.5 / n;  // Less aggressive shrink
        }
        else
        {
          m_current_rho   = 1.0;
          m_current_chi   = 1.0 + 2.0 / n;
          m_current_gamma = 0.75 - 0.5 / n;
          m_current_sigma = 1.0 - 1.0 / n;
        }
      }
      else
      {
        m_current_rho   = m_options.rho;
        m_current_chi   = m_options.chi;
        m_current_gamma = m_options.gamma;
        m_current_sigma = m_options.sigma;
      }
    }

    /**
     * @brief Get sorted indices of simplex vertices by function value
     * @return Const reference to sorted indices vector
     */
    vector<size_t> const &
    get_sorted_indices() const
    {
      if ( !m_simplex_ordered )
      {
        m_sorted_indices.resize( m_dim + 1 );
        std::iota( m_sorted_indices.begin(), m_sorted_indices.end(), 0 );
        std::sort( m_sorted_indices.begin(), m_sorted_indices.end(),
                   [this]( size_t i, size_t j ) { return m_values[i] < m_values[j]; } );
        m_simplex_ordered = true;
      }
      return m_sorted_indices;
    }

    /**
     * @brief Mark simplex as unordered (needs re-sorting)
     */
    void
    mark_simplex_unordered()
    {
      m_simplex_ordered = false;
    }

    /**
     * @brief Compute smart step size for simplex initialization
     * @param dimension_index Dimension index
     * @param current_val Current coordinate value
     * @return Appropriate step size
     */
    Scalar
    get_smart_step( size_t dimension_index, Scalar current_val ) const
    {
      Scalar step = m_options.initial_step;

      if ( m_use_bounds )
      {
        Scalar lower = m_lower( dimension_index );
        Scalar upper = m_upper( dimension_index );

        if ( std::isfinite( lower ) && std::isfinite( upper ) )
        {
          Scalar range = upper - lower;
          if ( range > std::numeric_limits<Scalar>::epsilon() )
          {
            step = std::clamp( range * 0.05, m_options.min_step_size, Scalar( 10.0 ) );
          }
        }
        else if ( std::isfinite( lower ) && !std::isfinite( upper ) )
        {
          Scalar dist_from_lower = current_val - lower;
          if ( dist_from_lower > 0 ) { step = max( m_options.initial_step, dist_from_lower * 0.1 ); }
        }
        else if ( !std::isfinite( lower ) && std::isfinite( upper ) )
        {
          Scalar dist_from_upper = upper - current_val;
          if ( dist_from_upper > 0 ) { step = max( m_options.initial_step, dist_from_upper * 0.1 ); }
        }
      }

      Scalar abs_val = abs( current_val );
      if ( abs_val > 1.0 ) step = max( m_options.initial_step, abs_val * 0.05 );

      return max( step, m_options.min_step_size );
    }

    /**
     * @brief Initialize simplex around starting point
     * @param x0 Starting point
     *
     * EIGEN3: Uses Eigen vector operations for efficient simplex construction
     * All vector operations are optimized and may use SIMD instructions
     */
    void
    initialize_simplex( Vector const & x0 )
    {
      m_dim = x0.size();

      m_simplex.resize( m_dim + 1 );
      m_values.resize( m_dim + 1 );
      m_centroid.resize( m_dim );
      m_trial_point.resize( m_dim );
      m_sorted_indices.resize( m_dim + 1 );

      // EIGEN3: Pre-allocate all vectors for efficient memory management
      for ( auto & vec : m_simplex ) vec.resize( m_dim );

      m_simplex[0] = x0;
      project_point( m_simplex[0] );
      m_values[0] = safe_evaluate( m_simplex[0] );

      Vector x_base   = m_simplex[0];
      Scalar min_step = m_options.min_step_size;

      for ( size_t i{ 0 }; i < m_dim; ++i )
      {
        Scalar step = get_smart_step( i, x_base( i ) );

        // Ensure minimum step to avoid degenerate simplex
        if ( abs( step ) < min_step ) step = ( step >= 0 ) ? min_step : -min_step;

        Vector x_next = x_base;
        x_next( i ) += step;

        project_point( x_next );

        // EIGEN3: Use .norm() for efficient distance computation
        // Eigen's norm() is optimized and may use vectorized operations
        if ( ( x_next - x_base ).norm() < min_step )
        {
          // If still too close, try opposite direction
          x_next = x_base;
          x_next( i ) -= step;
          project_point( x_next );
        }

        m_simplex[i + 1] = x_next;
        m_values[i + 1]  = safe_evaluate( x_next );
      }

      mark_simplex_unordered();

      // Verify simplex is not degenerate
      if ( compute_volume() < std::numeric_limits<Scalar>::epsilon() )
      {
        if ( m_options.verbose )
        {
          fmt::print(
              "{}[Warning] Initial simplex has zero volume, adding "
              "perturbation\n",
              m_indent );
        }
        // EIGEN3: Use Vector::Random() for efficient random vector generation
        // This is optimized and may use vectorized random number generation
        for ( size_t i{ 1 }; i <= m_dim; ++i )
        {
          Vector perturbation = Vector::Random( m_dim ) * min_step * 10;
          m_simplex[i] += perturbation;
          project_point( m_simplex[i] );
          m_values[i] = safe_evaluate( m_simplex[i] );
        }
      }

      if ( m_options.track_best_point )
      {
        auto indices = get_sorted_indices();
        if ( m_values[indices[0]] < m_best_value )
        {
          m_best_value = m_values[indices[0]];
          m_best_point = m_simplex[indices[0]];
        }
        m_previous_best = m_best_value;
      }
    }

    /**
     * @brief Update centroid excluding worst point
     * @param worst_index Index of worst point to exclude
     *
     * EIGEN3: Uses efficient vector accumulation with Eigen
     * setZero() and vector addition are optimized operations
     */
    void
    update_centroid( size_t worst_index )
    {
      m_centroid.setZero();
      for ( size_t i{ 0 }; i <= m_dim; ++i )
      {
        if ( i != worst_index ) m_centroid += m_simplex[i];
      }
      m_centroid /= static_cast<Scalar>( m_dim );
    }

    /**
     * @brief Compute simplex diameter (maximum vertex distance)
     * @return Simplex diameter
     *
     * EIGEN3: Uses Eigen's norm() function which is highly optimized
     * and may use SIMD instructions for distance computation
     */
    Scalar
    compute_diameter() const
    {
      Scalar max_dist{ 0 };
      for ( size_t i{ 0 }; i <= m_dim; ++i )
      {
        for ( size_t j{ i + 1 }; j <= m_dim; ++j )
        {
          Scalar dist = ( m_simplex[i] - m_simplex[j] ).norm();
          max_dist    = max( max_dist, dist );
        }
      }
      return max_dist;
    }

    /**
     * @brief Compute simplex volume
     * @return Simplex volume (approximated for high dimensions)
     *
     * EIGEN3: Uses Eigen's linear algebra routines for volume computation
     * For small dimensions: QR decomposition with column pivoting for numerical
     * stability For large dimensions: approximation using diameter to avoid
     * precision issues
     */
    Scalar
    compute_volume() const
    {
      if ( m_dim == 0 ) return 0;
      // EIGEN3: For high dimensions, use approximation to avoid numerical
      // instability
      if ( m_dim > 100 ) return std::pow( compute_diameter(), m_dim ) * std::exp( -std::lgamma( m_dim + 1 ) );

      // EIGEN3: Construct basis matrix using Eigen
      Matrix basis( m_dim, m_dim );
      for ( size_t i{ 0 }; i < m_dim; ++i ) { basis.col( i ) = m_simplex[i + 1] - m_simplex[0]; }

      // EIGEN3: Use QR decomposition with column pivoting for numerical
      // stability This handles rank-deficient cases gracefully
      Eigen::ColPivHouseholderQR<Matrix> qr( basis );
      if ( static_cast<size_t>( qr.rank() ) < m_dim ) return 0;

      // EIGEN3: Use Eigen's efficient determinant computation
      // logAbsDeterminant() is more stable for large matrices
      Scalar log_det   = qr.logAbsDeterminant();
      Scalar log_gamma = std::lgamma( m_dim + 1 );

      return std::exp( log_det - log_gamma );
    }

    /**
     * @brief Compute comprehensive simplex statistics
     * @return Simplex statistics structure
     *
     * EIGEN3: Uses Eigen's efficient statistical computations
     * Eigen::Map for zero-copy vector views and array operations for statistics
     */
    SimplexStats
    compute_simplex_stats() const
    {
      SimplexStats stats;

      stats.diameter = compute_diameter();

      // EIGEN3: Use Eigen::Map to view std::vector as Eigen vector without
      // copying This allows using Eigen's efficient statistical functions
      Scalar mean = Eigen::Map<const Eigen::VectorXd>( m_values.data(), m_values.size() ).mean();

      Scalar variance = 0;
      for ( auto const & v : m_values ) { variance += ( v - mean ) * ( v - mean ); }
      stats.std_dev = std::sqrt( variance / m_values.size() );

      auto indices      = get_sorted_indices();
      stats.value_range = m_values[indices.back()] - m_values[indices[0]];

      // EIGEN3: Efficient centroid computation using Eigen vector operations
      Vector centroid = Vector::Zero( m_dim );
      for ( auto const & v : m_simplex ) centroid += v;
      centroid /= m_simplex.size();

      Scalar total_dist = 0;
      for ( auto const & v : m_simplex ) { total_dist += ( v - centroid ).norm(); }
      stats.centroid_distance = total_dist / m_simplex.size();

      return stats;
    }

    /**
     * @brief Check convergence using robust criteria
     * @param best_value Best function value in simplex
     * @param worst_value Worst function value in simplex
     * @return True if converged
     */
    bool
    check_convergence_robust( Scalar best_value, Scalar worst_value ) const
    {
      auto const & tolerance              = m_options.tolerance;
      auto const & simplex_tolerance      = m_options.simplex_tolerance;
      auto const & convergence_relaxation = m_options.convergence_relaxation;

      auto                  stats = compute_simplex_stats();
      [[maybe_unused]] auto flags = compute_convergence_flags( best_value, worst_value, stats );

      if ( !m_options.use_robust_convergence )
      {
        // Usa worst_value qui
        bool value_converged    = ( worst_value - best_value ) < tolerance;
        bool geometry_converged = stats.diameter < simplex_tolerance;
        return value_converged && geometry_converged;
      }

      // 1. Primary convergence: function values
      bool value_converged;
      if ( m_options.use_relative_tolerance )
      {
        Scalar relative_range = stats.value_range / ( 1.0 + abs( best_value ) );
        value_converged       = relative_range < tolerance;
      }
      else
      {
        value_converged = stats.value_range < tolerance;
      }

      // 2. Secondary convergence: simplex geometry
      bool geometry_converged = stats.diameter < simplex_tolerance * std::sqrt( m_dim );

      // 3. Variance-based convergence
      bool variance_converged = stats.std_dev < tolerance * convergence_relaxation;

      // Converge if primary condition OR all secondary conditions are satisfied
      bool converged = value_converged ||
                       ( geometry_converged && stats.value_range < tolerance * convergence_relaxation ) ||
                       ( variance_converged && geometry_converged );

      if ( m_options.verbose && ( converged || m_global_iterations % m_options.progress_frequency == 0 ) )
      {
        fmt::print(
            "{}[Conv Check] V:{} G:{} S:{} | Range={:<12.4e} Diam={:<12.4e} "
            "StdDev={:<12.4e}\n",
            m_indent, ( value_converged ? "✓" : "✗" ), ( geometry_converged ? "✓" : "✗" ),
            ( variance_converged ? "✓" : "✗" ), stats.value_range, stats.diameter, stats.std_dev );
      }
      return converged;
    }

    /**
     * @brief Check for optimization stagnation
     * @param current_best Current best function value
     * @return True if stagnated
     */
    bool
    check_stagnation( Scalar current_best )
    {
      if ( !m_options.enable_restart ) return false;

      Scalar improvement          = abs( current_best - m_previous_best );
      Scalar relative_improvement = improvement / ( 1.0 + abs( m_previous_best ) );

      if ( relative_improvement < m_options.stagnation_tolerance )
      {
        ++m_stagnation_count;
        return m_stagnation_count >= m_options.stagnation_threshold;
      }
      else
      {
        m_stagnation_count = 0;
        m_previous_best    = current_best;
        return false;
      }
    }

    /**
     * @brief Determine if restart is worthwhile based on current results
     * @param current_result Current optimization result
     * @return True if restart should be performed
     */
    bool
    is_restart_worthwhile( InnerResult const & current_result ) const
    {
      // Don't restart if we've reached satisfactory tolerance
      if ( current_result.status == Status::CONVERGED ) return false;

      // Calculate relative improvement from start
      Scalar absolute_improvement = abs( current_result.initial_function_value - current_result.final_function_value );
      Scalar relative_improvement = absolute_improvement / ( 1.0 + abs( current_result.initial_function_value ) );

      // 1. Restart for stagnation with insufficient improvement
      if ( current_result.status == Status::STAGNATED )
      {
        if ( relative_improvement < m_options.restart_relative_improvement_threshold ) { return true; }
      }

      // 2. Restart for insufficient progress relative to resources used
      Scalar progress_per_eval = relative_improvement / ( 1.0 + current_result.function_evaluations );
      if ( progress_per_eval < m_options.restart_progress_per_eval_threshold && current_result.iterations > 100 )
      {
        return true;
      }

      // 3. Restart for problematic simplex geometry
      // EIGEN3: Use .norm() for efficient vector norm computation
      Scalar scale               = 1.0 + m_best_point.norm();
      Scalar normalized_diameter = current_result.simplex_diameter / scale;

      if ( normalized_diameter < m_options.simplex_tolerance * m_options.restart_simplex_diameter_factor1 &&
           relative_improvement < m_options.restart_simplex_geometry_threshold )
      {
        return true;
      }

      // 4. Restart for too many shrink operations without progress
      if ( m_shrink_count > m_options.restart_shrink_count_threshold &&
           relative_improvement < m_options.restart_after_shrink_improvement )
      {
        return true;
      }

      // 5. Restart for high-dimensional problems with slow progress
      if ( m_dim > 10 )
      {
        Scalar expected_progress = 1.0 / std::sqrt( 1.0 + current_result.iterations );
        if ( relative_improvement < expected_progress * m_options.restart_expected_progress_ratio &&
             current_result.iterations > 200 )
        {
          return true;
        }
      }

      // 6. Restart based on solution quality metric
      Scalar quality_metric = relative_improvement / ( 1.0 + std::log1p( current_result.function_evaluations ) );
      if ( quality_metric < m_options.restart_quality_metric_threshold && current_result.function_evaluations > 500 )
      {
        return true;
      }

      // 7. Restart if simplex is degenerate but not converged
      auto   stats              = compute_simplex_stats();
      Scalar normalized_std_dev = stats.std_dev / ( 1.0 + abs( m_best_value ) );
      if ( normalized_diameter < m_options.simplex_tolerance * m_options.restart_simplex_diameter_factor2 &&
           normalized_std_dev < m_options.tolerance * m_options.restart_std_dev_factor &&
           relative_improvement < m_options.restart_degenerate_improvement_threshold )
      {
        return true;
      }

      return false;
    }

    /**
     * @brief Adjust adaptive parameters based on algorithm behavior
     */
    void
    adaptive_parameter_adjustment()
    {
      // More conservative adaptive parameters
      if ( m_shrink_count > 8 )
      {
        m_current_sigma = max( 0.3, m_current_sigma * 0.9 );
        m_current_gamma = max( 0.3, m_current_gamma * 0.95 );
        m_shrink_count  = 0;
      }
    }

    /**
     * @brief Perform reflection operation
     * @param worst_index Index of worst point
     * @return Function value at reflected point
     *
     * EIGEN3: Uses Eigen vector arithmetic for efficient reflection computation
     * Expression templates enable efficient computation without temporaries
     */
    Scalar
    reflect_point( size_t worst_index )
    {
      m_trial_point = m_centroid + m_current_rho * ( m_centroid - m_simplex[worst_index] );
      project_point( m_trial_point );
      return safe_evaluate( m_trial_point );
    }

    /**
     * @brief Perform expansion operation
     * @return Function value at expanded point
     *
     * EIGEN3: Vector operations are optimized through expression templates
     */
    Scalar
    expand_point()
    {
      m_trial_point = m_centroid + m_current_chi * ( m_trial_point - m_centroid );
      project_point( m_trial_point );
      return safe_evaluate( m_trial_point );
    }

    /**
     * @brief Perform contraction operation
     * @param worst_index Index of worst point
     * @param outside Whether to contract outside or inside
     * @return Function value at contracted point
     *
     * EIGEN3: Efficient vector arithmetic using Eigen's expression system
     */
    Scalar
    contract_point( size_t worst_index, bool outside )
    {
      if ( outside ) { m_trial_point = m_centroid + m_current_gamma * ( m_trial_point - m_centroid ); }
      else
      {
        m_trial_point = m_centroid - m_current_gamma * ( m_centroid - m_simplex[worst_index] );
      }
      project_point( m_trial_point );
      return safe_evaluate( m_trial_point );
    }

    /**
     * @brief Perform shrink operation towards best point
     * @param best_index Index of best point
     *
     * EIGEN3: Uses efficient vector scaling and addition operations
     * Eigen's expression templates optimize these operations
     */
    void
    shrink_simplex( size_t best_index )
    {
      Vector best = m_simplex[best_index];
      for ( size_t i{ 0 }; i <= m_dim; ++i )
      {
        if ( i != best_index )
        {
          m_simplex[i] = best + m_current_sigma * ( m_simplex[i] - best );
          project_point( m_simplex[i] );
          m_values[i] = safe_evaluate( m_simplex[i] );
        }
      }
      ++m_shrink_count;
      // mark_simplex_unordered();
    }

    /**
     * @brief Perform one iteration of Nelder-Mead algorithm
     * @return True if converged
     *
     * EIGEN3: All vector operations in this method use Eigen's optimized
     * routines
     */
    bool
    nelder_mead_iteration()
    {
      auto   indices          = get_sorted_indices();
      size_t best_idx         = indices[0];
      size_t second_worst_idx = indices[m_dim - 1];
      size_t worst_idx        = indices[m_dim];

      Scalar best_value  = m_values[best_idx];
      Scalar worst_value = m_values[worst_idx];

      if ( m_options.track_best_point && best_value < m_best_value )
      {
        m_best_value = best_value;
        m_best_point = m_simplex[best_idx];
      }

      if ( check_convergence_robust( best_value, worst_value ) ) { return true; }

      update_centroid( worst_idx );
      Scalar f_reflect       = reflect_point( worst_idx );
      Vector reflected_point = m_trial_point;

      bool improve{ f_reflect < best_value };
      print_inner_operation( "Reflect", f_reflect, improve );
      if ( improve )
      {
        Scalar f_expand = expand_point();
        improve         = f_expand < f_reflect;
        print_inner_operation( "Expand", f_expand, improve );
        if ( f_expand < f_reflect )
        {
          m_simplex[worst_idx] = m_trial_point;
          m_values[worst_idx]  = f_expand;
        }
        else
        {
          m_simplex[worst_idx] = reflected_point;
          m_values[worst_idx]  = f_reflect;
        }
      }
      else if ( f_reflect < m_values[second_worst_idx] )
      {
        m_simplex[worst_idx] = reflected_point;
        m_values[worst_idx]  = f_reflect;
      }
      else
      {
        if ( f_reflect < worst_value )
        {
          Scalar f_contract = contract_point( worst_idx, true );
          improve           = f_contract <= f_reflect;
          print_inner_operation( "Contract(out)", f_contract, improve );
          if ( improve )
          {
            m_simplex[worst_idx] = m_trial_point;
            m_values[worst_idx]  = f_contract;
          }
          else
          {
            shrink_simplex( best_idx );
            adaptive_parameter_adjustment();
          }
        }
        else
        {
          Scalar f_contract = contract_point( worst_idx, false );
          improve           = f_contract < worst_value;
          print_inner_operation( "Contract(in)", f_contract, improve );
          if ( f_contract < worst_value )
          {
            m_simplex[worst_idx] = m_trial_point;
            m_values[worst_idx]  = f_contract;
          }
          else
          {
            shrink_simplex( best_idx );
            adaptive_parameter_adjustment();
          }
        }
      }
      mark_simplex_unordered();

      // Add convergence info printing
      if ( m_options.verbosity_level >= 2 )
      {
        auto stats = compute_simplex_stats();
        print_convergence_info( best_value, worst_value, stats );
      }

      return false;
    }

    void
    print_iteration_header( size_t iter, Scalar best_value, Scalar diameter ) const
    {
      if ( !m_options.verbose ) return;
      fmt::print( NelderMead::PrintColors::INFO,
                  "{}┌─────────────────────────┬───────────────┬───────────────┐\n"
                  "{}│ {:^23} │ {:^13} │ {:^13} │\n"
                  "{}└─────────────────────────┴───────────────┴───────────────┘\n",
                  m_indent, m_indent, fmt::format( "Iteration {}", iter ), fmt::format( "F = {:.4e}", best_value ),
                  fmt::format( "Diam = {:.4e}", diameter ), m_indent );
    }

    void
    print_inner_step( string const & step_name, Scalar fval, bool improved ) const
    {
      if ( !m_options.verbose ) return;

      auto   color = improved ? NelderMead::PrintColors::SUCCESS : NelderMead::PrintColors::ERROR;
      string icon  = improved ? "↗" : "↘";

      // Stampare le coordinate solo per problemi piccoli
      string tmp = m_dim > 5 ? "" : fmt::format( " | x = {}", NelderMead::format_vector<Scalar>( m_trial_point ) );
      fmt::print( color, "{}{:4} {} {:>12}: F = {:<12.6e}{}\n", m_indent, m_global_iterations, icon, step_name, fval,
                  tmp );
    }

    /**
     * @brief Run single Nelder-Mead optimization (without restarts)
     * @param x0 Starting point
     * @return Optimization result
     *
     * EIGEN3: All vector operations in this method leverage Eigen's
     * optimizations
     */
    InnerResult
    minimize_single_run( Vector const & x0 )
    {
      InnerResult result;
      size_t      initial_evals = m_global_evals;
      initialize_adaptive_parameters();

      if ( m_global_iterations >= m_options.max_iterations )
      {
        result.status               = Status::MAX_ITERATIONS;
        result.solution             = x0;
        result.final_function_value = safe_evaluate( x0 );
        result.function_evaluations = m_global_evals - initial_evals;
        // Set compatibility fields
        result.inner_iterations  = 0;
        result.inner_evaluations = result.function_evaluations;
        return result;
      }

      initialize_simplex( x0 );
      auto indices                  = get_sorted_indices();
      result.initial_function_value = m_values[indices[0]];

      // MODIFICA: Usare verbosity_level invece di verbose
      if ( m_options.verbosity_level >= 1 )
      {
        fmt::print( "{}[NM-Run]  Start | Dim={:<10} | F_0={:<12.6e}\n", m_indent, m_dim,
                    result.initial_function_value );
      }

      size_t local_iter = 0;

      while ( true )
      {
        if ( m_global_iterations >= m_options.max_iterations )
        {
          result.status = Status::MAX_ITERATIONS;
          break;
        }

        if ( m_global_evals >= m_options.max_function_evaluations )
        {
          result.status = Status::MAX_FUN_EVALUATIONS;
          break;
        }

        ++m_global_iterations;
        ++local_iter;

        // MODIFICA: Aggiungere stampa riepilogo periodico
        auto stats = compute_simplex_stats();
        indices    = get_sorted_indices();
        print_inner_iteration_summary( local_iter, m_values[indices[0]], stats );

        if ( m_options.enable_restart && local_iter % 50 == 0 )
        {
          indices = get_sorted_indices();
          if ( check_stagnation( m_values[indices[0]] ) )
          {
            result.status = Status::STAGNATED;
            break;
          }
        }

        if ( nelder_mead_iteration() )
        {
          result.status = Status::CONVERGED;
          break;
        }

        result.simplex_diameter = compute_diameter();
        if ( result.simplex_diameter < m_options.simplex_tolerance )
        {
          result.status = Status::SIMPLEX_TOO_SMALL;
          break;
        }

        // MODIFICA: Mantenere stampa progresso per compatibilità
        if ( m_options.verbose && ( m_global_iterations % m_options.progress_frequency ) == 0 )
        {
          indices = get_sorted_indices();
          fmt::print( "{}[NM-Iter] {:>5} | F={:<12.6e} | Diam={:<12.6e}\n", m_indent, m_global_iterations,
                      m_values[indices[0]], result.simplex_diameter );
        }
      }

      if ( result.status == Status::FAILED ) { result.status = Status::MAX_ITERATIONS; }

      indices                     = get_sorted_indices();
      result.solution             = m_simplex[indices[0]];
      result.final_function_value = m_values[indices[0]];
      result.simplex_volume       = compute_volume();
      result.simplex_diameter     = compute_diameter();
      result.iterations           = local_iter;
      result.function_evaluations = m_global_evals - initial_evals;
      result.shrink_operations    = m_shrink_count;

      // Set compatibility fields for BlockCoordinate
      result.inner_iterations  = result.iterations;
      result.inner_evaluations = result.function_evaluations;

      return result;
    }

    /**
     * @brief Print optimization header information
     * @param x0 Starting point
     */
    void
    print_header( Vector const & x0 ) const
    {
      if ( !m_options.verbose ) return;
      fmt::print( NelderMead::PrintColors::HEADER,
                  "{}"
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "{}║                    Nelder-Mead Optimization                    "
                  "║\n"
                  "{}"
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}║ {:62} ║\n"
                  "{}"
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_indent, m_indent, m_indent, m_indent, fmt::format( "Dimension: {:d}", x0.size() ), m_indent,
                  fmt::format( "Max Iterations: {:d}", m_options.max_iterations ), m_indent,
                  fmt::format( "Max Evaluations: {:d}", m_options.max_function_evaluations ), m_indent,
                  fmt::format( "Tolerance: {:.2e}", m_options.tolerance ), m_indent,
                  fmt::format( "Bounds: {}", ( m_use_bounds ? "Active" : "None" ) ), m_indent,
                  fmt::format( "Adaptive Parameters: {}", ( m_options.adaptive_parameters ? "Yes" : "No" ) ),
                  m_indent );
      fmt::print( "{}Initial point: {}\n", m_indent, NelderMead::format_vector<Scalar>( x0 ) );
    }

    /**
     * @brief Print optimization statistics
     * @param res Optimization result to print
     */
    void
    print_statistics( InnerResult const & res ) const
    {
      if ( !m_options.verbose ) return;
      fmt::print(
          "{}"
          "╔════════════════════════════════════════════════════════════════╗\n"
          "{}║                    Optimization Finished                       "
          "║\n"
          "{}"
          "╠════════════════════════════════════════════════════════════════╣\n"
          "{}║  Final Status       : {:<39}  ║\n"
          "{}║  Final Value        : {:<39.6e}  ║\n"
          "{}║  Total Iterations   : {:<39}  ║\n"
          "{}║  Total Evals        : {:<39}  ║\n"
          "{}║  Restarts           : {:<39}  ║\n"
          "{}║  Shrink Operations  : {:<39}  ║\n"
          "{}║  Simplex Diameter   : {:<39.6e}  ║\n"
          "{}"
          "╚════════════════════════════════════════════════════════════════╝\n"
          "\n",
          m_indent, m_indent, m_indent, m_indent, status_to_string( res.status ), m_indent, res.final_function_value,
          m_indent, res.iterations, m_indent, res.function_evaluations, m_indent, res.restarts_performed, m_indent,
          res.shrink_operations, m_indent, res.simplex_diameter, m_indent );
    }

  public:
    /**
     * @brief Construct Nelder-Mead optimizer with given options
     * @param opts Optimization options
     */
    explicit NelderMead_classic( Options const & opts = Options() ) : m_options( opts ) {}

    /**
     * @brief Set optimization options
     * @param opts New options
     */
    void
    set_options( Options const & opts )
    {
      m_options = opts;
    }

    /**
     * @brief Set optimization bounds
     * @param lower Lower bounds
     * @param upper Upper bounds
     *
     * EIGEN3: Uses Eigen's array operations for bounds validation
     */
    void
    set_bounds( Vector const & lower, Vector const & upper )
    {
      UTILS_ASSERT( lower.size() == upper.size(), "Bounds size mismatch" );
      UTILS_ASSERT( ( lower.array() <= upper.array() ).all(), "Lower <= Upper" );
      m_lower      = lower;
      m_upper      = upper;
      m_use_bounds = true;
    }

    /**
     * @brief Clear optimization bounds
     */
    void
    clear_bounds()
    {
      m_use_bounds = false;
    }

    /**
     * @brief Minimize objective function starting from x0
     * @param x0 Starting point
     * @param callback Objective function
     * @return Optimization result
     *
     * EIGEN3: Main optimization routine leveraging all Eigen optimizations
     */
    InnerResult
    minimize( Vector const & x0, Callback const & callback )
    {
      m_callback = &callback;

      m_stagnation_count = 0;
      m_shrink_count     = 0;
      // Initialize with NaN to handle negative values properly
      m_best_value        = std::numeric_limits<Scalar>::max();
      m_previous_best     = std::numeric_limits<Scalar>::max();
      m_global_iterations = 0;
      m_global_evals      = 0;
      m_simplex_ordered   = false;

      size_t restarts = 0;

      print_header( x0 );

      InnerResult best_result = minimize_single_run( x0 );

      // Initialize m_best_value if not yet set
      if ( std::isnan( m_best_value ) )
      {
        m_best_value = best_result.final_function_value;
        m_best_point = best_result.solution;
      }

      while ( m_options.enable_restart && restarts < m_options.max_restarts && is_restart_worthwhile( best_result ) )
      {
        if ( m_global_evals >= m_options.max_function_evaluations ) break;
        if ( m_global_iterations >= m_options.max_iterations ) break;

        if ( m_options.verbose )
        {
          fmt::print( "{}[NM-Restart] #{}/{} | Reason: {:<16} | F={:12.6e}\n", m_indent, ( restarts + 1 ),
                      m_options.max_restarts, status_to_string( best_result.status ),
                      best_result.final_function_value );
        }

        Scalar perturbation_scale = m_options.restart_perturbation_ratio * ( 1.0 + restarts * 0.1 );

        Vector restart_x0;
        Vector scale_vec = Vector::Ones( x0.size() );

        size_t x_size = static_cast<size_t>( x0.size() );
        if ( m_use_bounds )
        {
          for ( size_t i{ 0 }; i < x_size; ++i )
          {
            Scalar r       = m_upper( i ) - m_lower( i );
            scale_vec( i ) = std::isfinite( r ) ? r : max( Scalar( 1.0 ), abs( best_result.solution( i ) ) );
          }
        }
        else
        {
          for ( size_t i{ 0 }; i < x_size; ++i )
          {
            scale_vec( i ) = max( Scalar( 1.0 ), abs( best_result.solution( i ) ) );
          }
        }

        // EIGEN3: Use Vector::Random() and cwiseProduct for efficient random
        // perturbation
        Vector perturbation = Vector::Random( x0.size() ).cwiseProduct( scale_vec ) * perturbation_scale;

        if ( m_options.track_best_point ) { restart_x0 = m_best_point + perturbation; }
        else
        {
          restart_x0 = best_result.solution + perturbation;
        }

        project_point( restart_x0 );

        m_stagnation_count = 0;
        m_shrink_count     = 0;
        m_simplex_ordered  = false;

        InnerResult current_result = minimize_single_run( restart_x0 );
        ++restarts;

        // ROBUST IMPROVEMENT CONDITION FOR POSITIVE AND NEGATIVE VALUES
        bool improvement = false;

        if ( std::isnan( current_result.final_function_value ) ) { improvement = false; }
        else if ( std::isnan( best_result.final_function_value ) ) { improvement = true; }
        else
        {
          // Calculate normalized relative improvement
          Scalar abs_best = abs( best_result.final_function_value );
          // Scalar abs_current = abs(current_result.final_function_value);

          // Case 1: Both positive or zero - standard improvement
          if ( best_result.final_function_value >= 0 && current_result.final_function_value >= 0 )
          {
            improvement = current_result.final_function_value <
                          best_result.final_function_value * m_options.restart_improvement_ratio;
          }
          // Case 2: Both negative - improvement means more negative
          else if ( best_result.final_function_value < 0 && current_result.final_function_value < 0 )
          {
            improvement = current_result.final_function_value < best_result.final_function_value;
          }
          // Case 3: Transition from positive to negative - always improvement
          else if ( best_result.final_function_value >= 0 && current_result.final_function_value < 0 )
          {
            improvement = true;
          }
          // Case 4: Transition from negative to positive - usually worsening
          else
          {
            improvement = false;
          }

          // Additional check: significant absolute improvement
          Scalar absolute_improvement = best_result.final_function_value - current_result.final_function_value;
          Scalar relative_improvement = absolute_improvement / ( 1.0 + abs_best );

          if ( !improvement && relative_improvement > m_options.restart_absolute_improvement_threshold )
          {
            improvement = true;
          }
        }

        if ( improvement )
        {
          best_result                    = current_result;
          best_result.restarts_performed = restarts;
        }
        else if ( m_options.verbose ) { fmt::print( "{}[Restart rejected: no improvement]\n", m_indent ); }
      }

      best_result.iterations           = m_global_iterations;
      best_result.function_evaluations = m_global_evals;
      best_result.restarts_performed   = restarts;
      best_result.shrink_operations    = m_shrink_count;

      // Final best point update - correct comparison for minimization
      if ( m_options.track_best_point )
      {
        if ( std::isnan( best_result.final_function_value ) ||
             ( !std::isnan( m_best_value ) && m_best_value < best_result.final_function_value ) )
        {
          best_result.final_function_value = m_best_value;
          best_result.solution             = m_best_point;
        }
      }

      print_statistics( best_result );

      return best_result;
    }

    /**
     * @brief Convert status enum to string
     * @param status Optimization status
     * @return String representation of status
     */
    static string
    status_to_string( Status status )
    {
      switch ( status )
      {
        case Status::CONVERGED:
          return "CONVERGED";
        case Status::MAX_ITERATIONS:
          return "MAX_ITERATIONS";
        case Status::MAX_FUN_EVALUATIONS:
          return "MAX_FUN_EVALUATIONS";
        case Status::SIMPLEX_TOO_SMALL:
          return "SIMPLEX_TOO_SMALL";
        case Status::STAGNATED:
          return "STAGNATED";
        case Status::FAILED:
          return "FAILED";
        default:
          return "UNKNOWN";
      }
    }

    // Access methods for debugging and monitoring
    size_t
    get_total_evaluations() const
    {
      return m_global_evals;
    }  ///< Get total function evaluations
    size_t
    get_total_iterations() const
    {
      return m_global_iterations;
    }  ///< Get total iterations
    Scalar
    get_best_value() const
    {
      return m_best_value;
    }  ///< Get best function value found
    Vector
    get_best_point() const
    {
      return m_best_point;
    }  ///< Get best point found
  };

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
  template <typename Scalar = double>
  class NelderMead_BlockCoordinate
  {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar( Vector const & )>;
    using Result   = NelderMead::Result<Scalar>;

    struct Options
    {
      size_t block_size{ 10 };
      size_t max_outer_iterations{ 100 };
      size_t max_inner_iterations{ 1000 };
      size_t max_function_evaluations{ 100000 };
      Scalar tolerance{ 1e-6 };
      bool   verbose{ true };
      size_t verbosity_level{ 1 };  // 0: quiet, 1: outer stats, 2: inner progress, 3: detailed
      size_t inner_progress_frequency{ 10 };
      typename NelderMead_classic<Scalar>::Options sub_options;
    };

  private:
    Options                    m_options;
    NelderMead_classic<Scalar> m_solver;
    Vector                     m_lower;
    Vector                     m_upper;
    bool                       m_use_bounds{ false };
    string                     m_indent{ "" };

    // Cyclic block selection
    vector<size_t>
    select_block( size_t n_dims, size_t iter )
    {
      vector<size_t> indices;
      indices.reserve( m_options.block_size );
      size_t start = ( iter * m_options.block_size ) % n_dims;
      for ( size_t i = 0; i < m_options.block_size; ++i ) indices.push_back( ( start + i ) % n_dims );
      std::sort( indices.begin(), indices.end() );
      indices.erase( std::unique( indices.begin(), indices.end() ), indices.end() );
      return indices;
    }

    void
    print_outer_iteration_header( size_t                 outer_iter,
                                  size_t                 total_cycles,
                                  const vector<size_t> & block_indices,
                                  size_t                 block_size ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( NelderMead::PrintColors::HEADER,
                  "\n"
                  "{}"
                  "╔════════════════════════════════════════════════════════════════╗\n"
                  "{}║ Outer Iteration {:3d} - Block {:2d}/{:2d}                       "
                  "       ║\n"
                  "{}"
                  "╠════════════════════════════════════════════════════════════════╣\n"
                  "{}║ Block Indices: {:<47} ║\n"
                  "{}║ Block Size:    {:<47} ║\n"
                  "{}"
                  "╚════════════════════════════════════════════════════════════════╝"
                  "\n",
                  m_indent, m_indent, outer_iter, ( outer_iter % total_cycles ) + 1, total_cycles, m_indent, m_indent,
                  Utils::format_index_vector<size_t>( block_indices, 7 ), m_indent, block_size, m_indent );
    }

    void
    print_outer_iteration_result( size_t outer_iter,
                                  size_t block_size,
                                  Scalar current_f,
                                  size_t inner_iters,
                                  size_t inner_evals,
                                  Scalar improvement,
                                  bool   improved ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      auto color = improved ? NelderMead::PrintColors::SUCCESS : NelderMead::PrintColors::WARNING;

      fmt::print( color,
                  "{}"
                  "┌──────────────┬───────────┬────────────────┬───────────────"
                  "───┬──────────────────┐\n"
                  "{}│ out iter:{:<3} │ block:{:<3} │ F:{:<12.6e} │ inner "
                  "iter:{:<5} │ inner eval:{:<5} │\n"
                  "{}"
                  "└──────────────┴───────────┴────────────────┴───────────────"
                  "───┴──────────────────┘\n",
                  m_indent, m_indent, outer_iter, block_size, current_f, inner_iters, inner_evals, m_indent );

      if ( improved )
      {
        fmt::print( NelderMead::PrintColors::SUCCESS, "{}✓ Improvement: {:.6e} → {:.6e} (Δ = {:.6e})\n", m_indent,
                    current_f + improvement, current_f, improvement );
      }
      else
      {
        fmt::print( NelderMead::PrintColors::WARNING, "{}⚠ No significant improvement (Δ = {:.6e})\n", m_indent,
                    improvement );
      }
    }

    void
    print_outer_statistics( size_t outer_iter,
                            size_t total_outer_iters,
                            size_t total_inner_iters,
                            size_t total_evals,
                            Scalar best_value,
                            bool   converged ) const
    {
      if ( m_options.verbosity_level < 1 ) return;

      fmt::print( NelderMead::PrintColors::INFO,
                  "\n"
                  "{}"
                  "╔═══════════════════════════════════════════════════════════════╗\n"
                  "{}║                      Outer Iteration Summary                  "
                  "║\n"
                  "{}"
                  "╠═══════════════════════════════════════════════════════════════╣\n"
                  "{}║ Completed:        {:<43} ║\n"
                  "{}║ Outer Iterations: {:<43} ║\n"
                  "{}║ Inner Iterations: {:<43} ║\n"
                  "{}║ Total Evals:      {:<43} ║\n"
                  "{}║ Best Value:       {:<43.6e} ║\n"
                  "{}║ Status:           {:<43} ║\n"
                  "{}"
                  "╚═══════════════════════════════════════════════════════════════╝\n",
                  m_indent, m_indent, m_indent, m_indent, fmt::format( "{}/{}", outer_iter, total_outer_iters ),
                  m_indent, outer_iter, m_indent, total_inner_iters, m_indent, total_evals, m_indent, best_value,
                  m_indent, ( converged ? "CONVERGED" : "RUNNING" ), m_indent );
    }

  public:
    explicit NelderMead_BlockCoordinate( Options const & opts = Options() ) : m_options( opts )
    {
      m_solver.set_options( m_options.sub_options );
    }

    void
    set_sub_options( typename NelderMead_classic<Scalar>::Options const & sopts )
    {
      m_options.sub_options = sopts;
      m_solver.set_options( sopts );
    }

    void
    set_bounds( Vector const & l, Vector const & u )
    {
      m_lower      = l;
      m_upper      = u;
      m_use_bounds = true;
    }

    /**
     * @brief Main optimization routine for block coordinate descent
     * @param x Initial point (updated in-place)
     * @param global_callback Objective function
     * @return Optimization result
     *
     * EIGEN3: Efficient block processing using Eigen vector operations
     * - Subspace extraction without memory copies
     * - Efficient bounds management for subspaces
     * - Leverages Eigen's expression templates for performance
     */
    Result
    minimize( Vector x, Callback const & global_callback )
    {
      size_t n = x.size();

      // 1. Check for small problem fallback
      if ( !( 2 * n > m_options.block_size ) )
      {
        // MODIFICA: Usare verbosity_level
        if ( m_options.verbosity_level >= 1 )
          fmt::print(
              "{}[Info] Dim <= BlockSize. Switching to DIRECT CLASSIC "
              "solver.\n",
              m_indent );

        auto full_opts                     = m_options.sub_options;
        full_opts.max_function_evaluations = m_options.max_function_evaluations;
        full_opts.max_iterations           = m_options.max_inner_iterations;
        full_opts.tolerance                = m_options.tolerance;
        full_opts.verbose                  = ( m_options.verbosity_level >= 1 );  // MODIFICA

        m_solver.set_options( full_opts );
        if ( m_use_bounds ) m_solver.set_bounds( m_lower, m_upper );
        auto inner_res = m_solver.minimize( x, global_callback );

        // Convert InnerResult to Result
        Result res;
        res.solution               = inner_res.solution;
        res.final_function_value   = inner_res.final_function_value;
        res.initial_function_value = inner_res.initial_function_value;
        res.outer_iterations       = 1;
        res.inner_iterations       = inner_res.inner_iterations;
        res.total_iterations       = 1 + inner_res.inner_iterations;
        res.outer_evaluations      = 1;  // initial evaluation
        res.inner_evaluations      = inner_res.inner_evaluations;
        res.total_evaluations      = 1 + inner_res.inner_evaluations;
        res.status                 = inner_res.status == NelderMead_classic<Scalar>::Status::CONVERGED
                                         ? NelderMead::Status::CONVERGED
                                         : NelderMead::Status::MAX_ITERATIONS;
        return res;
      }

      // 2. Block Coordinate Logic
      Result res;
      size_t count_outer_evals = 0;
      size_t count_inner_evals = 0;
      size_t count_inner_iters = 0;

      Scalar current_f = global_callback( x );
      count_outer_evals++;
      res.initial_function_value = current_f;

      if ( std::isnan( current_f ) || std::isinf( current_f ) )
      {
        res.status   = NelderMead::Status::FAIL_NAN;
        res.solution = x;
        return res;
      }

      size_t outer_iter = 0;
      bool   converged  = false;

      // Stagnation control
      size_t no_improv_count  = 0;
      size_t blocks_per_cycle = ( n + m_options.block_size - 1 ) / m_options.block_size;

      while ( outer_iter < m_options.max_outer_iterations && !converged )
      {
        outer_iter++;
        if ( count_outer_evals + count_inner_evals >= m_options.max_function_evaluations )
        {
          res.status = NelderMead::Status::MAX_FUN_EVALUATIONS;
          break;
        }

        vector<size_t> idxs = select_block( n, outer_iter - 1 );
        size_t         k    = idxs.size();

        // MODIFICA: Usare il nuovo metodo di stampa
        print_outer_iteration_header( outer_iter, blocks_per_cycle, idxs, k );

        // Prepare Subspace
        // EIGEN3: Efficient subspace extraction without memory allocation
        Vector x_sub( k ), l_sub( k ), u_sub( k );
        for ( size_t i{ 0 }; i < k; ++i )
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
          for ( size_t i = 0; i < k; ++i ) x_temp( idxs[i] ) = sub_v( i );
          return global_callback( x_temp );
        };

        // Run Inner Solver
        auto sub_res = m_solver.minimize( x_sub, sub_cb );
        count_inner_iters += sub_res.inner_iterations;

        // Update Global State
        Scalar old_f       = current_f;
        Scalar improvement = old_f - sub_res.final_function_value;
        bool   improved    = sub_res.final_function_value < current_f;

        if ( improved )
        {
          current_f = sub_res.final_function_value;
          for ( size_t i = 0; i < k; ++i ) x( idxs[i] ) = sub_res.solution( i );
        }

        // MODIFICA: Usare il nuovo metodo di stampa
        print_outer_iteration_result( outer_iter, k, current_f, sub_res.inner_iterations, sub_res.inner_evaluations,
                                      improvement, improved );

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
          print_outer_statistics( outer_iter, m_options.max_outer_iterations, count_inner_iters,
                                  count_outer_evals + count_inner_evals, current_f, converged );
        }
      }

      res.solution             = x;
      res.final_function_value = current_f;
      res.outer_iterations     = outer_iter;
      res.inner_iterations     = count_inner_iters;
      res.total_iterations     = outer_iter + count_inner_iters;
      res.outer_evaluations    = count_outer_evals;
      res.inner_evaluations    = count_inner_evals;
      res.total_evaluations    = count_outer_evals + count_inner_evals;

      if ( res.status == NelderMead::Status::FAILED )
      {
        if ( converged )
          res.status = NelderMead::Status::CONVERGED;
        else if ( outer_iter >= m_options.max_outer_iterations )
          res.status = NelderMead::Status::MAX_ITERATIONS;
      }

      // MODIFICA: Stampare statistiche finali
      if ( m_options.verbosity_level >= 1 )
      {
        print_outer_statistics( outer_iter, m_options.max_outer_iterations, count_inner_iters,
                                count_outer_evals + count_inner_evals, current_f, true );
      }

      return res;
    }
  };
}  // namespace Utils

#endif
