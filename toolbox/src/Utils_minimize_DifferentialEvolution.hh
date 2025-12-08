// File: Utils_DifferentialEvolution.hh

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

#ifndef UTILS_DIFFERENTIAL_EVOLUTION_dot_HH
#define UTILS_DIFFERENTIAL_EVOLUTION_dot_HH

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"
#include "Utils_nonlinear_system.hh"

namespace Utils
{

  // =============================================================================
  // DifferentialEvolution
  // Header-only implementation of Differential Evolution algorithm
  // Features:
  // - Memory layout optimized: each individual is stored as a COLUMN (Dim x NP)
  // - Multiple mutation strategies (DE/rand/1, DE/best/1, etc.)
  // - Advanced strategies: current-to-pbest, JADE-style adaptation
  // - Boundary constraint handling using Eigen vector operations
  // - Configurable population size, crossover rate, and differential weight
  // - Parameter adaptation (JADE/SHADE style)
  // - Population size reduction (L-SHADE style)
  // - Archive of inferior solutions for increased diversity
  // - Parallel fitness evaluation (OpenMP)
  // - Convergence monitoring and verbose output
  // - Easy integration with existing codebases
  // =============================================================================

  template <typename RealType = double>
  class DifferentialEvolution
  {
  public:
    using integer   = Eigen::Index;
    using real_type = RealType;
    using Vector    = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
    using Matrix    = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

    // Mutation strategies
    enum Strategy
    {
      RAND_1                        = 1,  // DE/rand/1
      BEST_1                        = 2,  // DE/best/1
      RAND_TO_BEST_1                = 3,  // DE/rand-to-best/1
      BEST_2                        = 4,  // DE/best/2
      RAND_2                        = 5,  // DE/rand/2
      RAND_1_DITHER                 = 6,  // DE/rand/1 with dither
      EITHER_OR                     = 7,  // Either-or algorithm
      CURRENT_TO_PBEST              = 8,  // DE/current-to-pbest/1 (JADE)
      CURRENT_TO_PBEST_WITH_ARCHIVE = 9   // JADE with archive
    };

    // Constraint handling methods
    enum ConstraintMethod
    {
      BOUNCE_BACK      = 1,
      CLAMPING         = 2,
      RANDOM_REINIT    = 3,
      COSINE_TRANSFORM = 4,
      HYBRID           = 5
    };

  private:
    // Algorithm parameters
    integer   m_population_size = 50;      // NP
    integer   m_dimension       = 1;       // Dim
    integer   m_max_iterations  = 1000;    // maxIter
    Strategy  m_strategy        = RAND_1;  // mutation strategy
    real_type m_weight          = 0.8;     // F (differential weight)
    real_type m_crossover_rate  = 0.9;     // CR
    real_type m_dither_min      = 0.5;     // minimum dither factor
    real_type m_dither_max      = 1.0;     // maximum dither factor

    // Boundary constraints
    Vector           m_lower_bounds;
    Vector           m_upper_bounds;
    bool             m_use_bounds        = false;
    ConstraintMethod m_constraint_method = BOUNCE_BACK;

    // Advanced features flags
    bool m_use_parameter_adaptation = false;
    bool m_use_population_reduction = false;
    bool m_use_archive              = false;
    bool m_enable_parallel          = false;
    bool m_use_local_search         = false;

    // Algorithm state - MEMORY LAYOUT: Dim x NP (each column is an individual)
    Matrix               m_population;     // Dim x NP
    Vector               m_fitness;        // NP x 1
    Vector               m_best_solution;  // Dim x 1
    real_type            m_best_fitness      = std::numeric_limits<real_type>::max();
    integer              m_current_iteration = 0;
    std::atomic<integer> m_function_evaluations{ 0 };

    // Random number generation
    std::mt19937                              m_random_engine;
    std::uniform_real_distribution<real_type> m_uniform_dist{ 0.0, 1.0 };

    // Convergence criteria
    real_type m_tolerance           = 1e-8;
    integer   m_stagnation_window   = 50;
    real_type m_value_to_reach      = -std::numeric_limits<real_type>::infinity();
    real_type m_diversity_threshold = 1e-12;

    // Output control
    bool    m_verbose        = false;
    integer m_print_interval = 100;

    // Statistics
    std::vector<real_type> m_best_history;
    bool                   m_converged = false;
    std::string            m_convergence_reason;

    // =========================================================================
    // Advanced features data structures
    // =========================================================================

    // For parameter adaptation (JADE/SHADE style)
    std::vector<real_type> m_successful_F;
    std::vector<real_type> m_successful_CR;
    real_type              m_memory_F    = 0.5;
    real_type              m_memory_CR   = 0.5;
    integer                m_memory_size = 100;

    // For population reduction (L-SHADE style)
    integer m_initial_population_size = 100;
    integer m_min_population_size     = 4;

    // Archive of inferior solutions
    Matrix  m_archive;
    integer m_archive_capacity     = 100;
    integer m_current_archive_size = 0;

    // Cache for fitness evaluations (optional)
    std::unordered_map<size_t, real_type> m_fitness_cache;
    bool                                  m_use_cache = false;
    mutable std::mutex                    m_cache_mutex;

    // Local search parameters
    real_type m_local_search_threshold  = 1e-6;
    integer   m_local_search_iterations = 100;

    // Callbacks
    std::function<void( integer, real_type )>        m_iteration_callback;
    std::function<void( const Vector &, real_type )> m_new_best_callback;

  private:
    // =========================================================================
    // Helper methods
    // =========================================================================

    // Hash function for caching
    size_t
    hash_vector( const Vector & v ) const
    {
      size_t hash = 0;
      for ( integer i = 0; i < v.size(); ++i )
      {
        hash ^= std::hash<real_type>{}( v[i] ) + 0x9e3779b9 + ( hash << 6 ) + ( hash >> 2 );
      }
      return hash;
    }

    // Generate random permutation
    void
    generate_permutation( std::vector<integer> & perm )
    {
      std::iota( perm.begin(), perm.end(), 0 );
      std::shuffle( perm.begin(), perm.end(), m_random_engine );
    }

    // Apply boundary constraints with multiple methods
    void
    apply_boundary_constraints( Vector & individual, const Vector & original )
    {
      if ( !m_use_bounds ) return;

      switch ( m_constraint_method )
      {
        case BOUNCE_BACK:
          for ( integer j = 0; j < m_dimension; ++j )
          {
            real_type & value = individual[j];
            if ( value < m_lower_bounds[j] )
            {
              value = m_lower_bounds[j] + m_uniform_dist( m_random_engine ) * ( original[j] - m_lower_bounds[j] );
            }
            else if ( value > m_upper_bounds[j] )
            {
              value = m_upper_bounds[j] + m_uniform_dist( m_random_engine ) * ( original[j] - m_upper_bounds[j] );
            }
          }
          break;

        case CLAMPING:
          for ( integer j = 0; j < m_dimension; ++j )
          {
            individual[j] = std::max( m_lower_bounds[j], std::min( m_upper_bounds[j], individual[j] ) );
          }
          break;

        case RANDOM_REINIT:
          for ( integer j = 0; j < m_dimension; ++j )
          {
            if ( individual[j] < m_lower_bounds[j] || individual[j] > m_upper_bounds[j] )
            {
              individual[j] = m_lower_bounds[j] +
                              m_uniform_dist( m_random_engine ) * ( m_upper_bounds[j] - m_lower_bounds[j] );
            }
          }
          break;

        case COSINE_TRANSFORM:
          for ( integer j = 0; j < m_dimension; ++j )
          {
            real_type & value = individual[j];
            if ( value < m_lower_bounds[j] || value > m_upper_bounds[j] )
            {
              // Transform to [-pi/2, pi/2] range
              real_type transformed = std::asin(
                  2.0 * ( value - m_lower_bounds[j] ) / ( m_upper_bounds[j] - m_lower_bounds[j] ) - 1.0 );
              // Add small random perturbation
              transformed += ( m_uniform_dist( m_random_engine ) - 0.5 ) * 0.1;
              // Transform back
              value = m_lower_bounds[j] +
                      ( std::sin( transformed ) + 1.0 ) * 0.5 * ( m_upper_bounds[j] - m_lower_bounds[j] );
            }
          }
          break;

        case HYBRID:
          // Use bounce-back for first half, random reinit for second half
          for ( integer j = 0; j < m_dimension; ++j )
          {
            real_type & value = individual[j];
            if ( value < m_lower_bounds[j] )
            {
              if ( j < m_dimension / 2 )
              {
                value = m_lower_bounds[j] + m_uniform_dist( m_random_engine ) * ( original[j] - m_lower_bounds[j] );
              }
              else
              {
                value = m_lower_bounds[j] +
                        m_uniform_dist( m_random_engine ) * ( m_upper_bounds[j] - m_lower_bounds[j] );
              }
            }
            else if ( value > m_upper_bounds[j] )
            {
              if ( j < m_dimension / 2 )
              {
                value = m_upper_bounds[j] + m_uniform_dist( m_random_engine ) * ( original[j] - m_upper_bounds[j] );
              }
              else
              {
                value = m_lower_bounds[j] +
                        m_uniform_dist( m_random_engine ) * ( m_upper_bounds[j] - m_lower_bounds[j] );
              }
            }
          }
          break;
      }
    }

    // Calculate population diversity
    real_type
    calculate_diversity() const
    {
      if ( m_population_size <= 1 ) return 0.0;

      Vector    mean      = m_population.rowwise().mean();
      real_type diversity = 0.0;

      for ( integer i = 0; i < m_population_size; ++i ) { diversity += ( m_population.col( i ) - mean ).squaredNorm(); }

      return std::sqrt( diversity / m_population_size );
    }

    // Adapt parameters based on success (JADE style)
    void
    adapt_parameters( real_type F, real_type CR, bool success )
    {
      if ( !m_use_parameter_adaptation ) return;

      if ( success )
      {
        m_successful_F.push_back( F );
        m_successful_CR.push_back( CR );

        // Keep memory bounded
        if ( m_successful_F.size() > static_cast<size_t>( m_memory_size ) )
        {
          m_successful_F.erase( m_successful_F.begin() );
          m_successful_CR.erase( m_successful_CR.begin() );
        }
      }

      // Update memory values using Lehmer mean for F, arithmetic mean for CR
      if ( !m_successful_F.empty() )
      {
        real_type sum_F    = 0.0;
        real_type sum_F_sq = 0.0;

        for ( real_type f : m_successful_F )
        {
          sum_F += f;
          sum_F_sq += f * f;
        }

        m_memory_F = sum_F_sq / sum_F;  // Lehmer mean

        m_memory_CR = std::accumulate( m_successful_CR.begin(), m_successful_CR.end(), 0.0 ) / m_successful_CR.size();
      }
    }

    // Generate F and CR values (adapted or fixed)
    void
    generate_parameters( real_type & F, real_type & CR )
    {
      if ( m_use_parameter_adaptation )
      {
        // Generate from Cauchy and Normal distributions around memory values
        std::cauchy_distribution<real_type> cauchy_dist( m_memory_F, 0.1 );
        std::normal_distribution<real_type> normal_dist( m_memory_CR, 0.1 );

        F  = cauchy_dist( m_random_engine );
        CR = normal_dist( m_random_engine );

        // Clamp to valid ranges
        F  = std::max<real_type>( 0.1, std::min<real_type>( 1.0, F ) );
        CR = std::max<real_type>( 0.0, std::min<real_type>( 1.0, CR ) );
      }
      else
      {
        F  = m_weight;
        CR = m_crossover_rate;
      }
    }

    // Reduce population size (L-SHADE style)
    void
    reduce_population()
    {
      if ( !m_use_population_reduction ) return;

      if ( m_current_iteration > m_max_iterations / 2 )
      {
        integer new_size = std::max( m_min_population_size,
                                     static_cast<integer>(
                                         m_initial_population_size *
                                         ( 1.0 - static_cast<real_type>( m_current_iteration ) / m_max_iterations ) ) );

        if ( new_size < m_population_size )
        {
          // Sort individuals by fitness and keep the best ones
          std::vector<integer> indices( m_population_size );
          std::iota( indices.begin(), indices.end(), 0 );

          std::sort( indices.begin(), indices.end(),
                     [this]( integer a, integer b ) { return m_fitness[a] < m_fitness[b]; } );

          Matrix new_population( m_dimension, new_size );
          Vector new_fitness( new_size );

          for ( integer i = 0; i < new_size; ++i )
          {
            new_population.col( i ) = m_population.col( indices[i] );
            new_fitness[i]          = m_fitness[indices[i]];
          }

          m_population      = std::move( new_population );
          m_fitness         = std::move( new_fitness );
          m_population_size = new_size;
        }
      }
    }

    // Add to archive
    void
    add_to_archive( const Vector & individual )
    {
      if ( !m_use_archive || m_current_archive_size >= m_archive_capacity ) return;

      if ( m_current_archive_size < m_archive.cols() ) { m_archive.col( m_current_archive_size ) = individual; }
      else
      {
        // Resize archive if needed
        m_archive.conservativeResize( m_dimension, m_archive.cols() + m_archive_capacity );
        m_archive.col( m_current_archive_size ) = individual;
      }
      m_current_archive_size++;
    }

    // Simple local search (gradient-free)
    template <typename ObjectiveFunction>
    void
    local_search( Vector & solution, real_type & fitness, ObjectiveFunction && objective, integer max_iterations = 100 )
    {
      real_type step_size    = 0.01;
      Vector    best         = solution;
      real_type best_fitness = fitness;

      for ( integer iter = 0; iter < max_iterations; ++iter )
      {
        bool improved = false;

        for ( integer j = 0; j < m_dimension; ++j )
        {
          // Try positive perturbation
          Vector candidate = solution;
          candidate[j] += step_size;
          apply_boundary_constraints( candidate, solution );
          real_type cand_fitness = objective( candidate );
          m_function_evaluations++;

          if ( cand_fitness < best_fitness )
          {
            best         = candidate;
            best_fitness = cand_fitness;
            improved     = true;
          }

          // Try negative perturbation
          candidate = solution;
          candidate[j] -= step_size;
          apply_boundary_constraints( candidate, solution );
          cand_fitness = objective( candidate );
          m_function_evaluations++;

          if ( cand_fitness < best_fitness )
          {
            best         = candidate;
            best_fitness = cand_fitness;
            improved     = true;
          }
        }

        if ( improved )
        {
          solution = best;
          fitness  = best_fitness;
          step_size *= 1.2;  // Increase step if improving
        }
        else
        {
          step_size *= 0.5;  // Decrease step if not improving
        }

        if ( step_size < 1e-10 ) break;
      }
    }

  public:
    // =========================================================================
    // Constructors
    // =========================================================================

    DifferentialEvolution()
    {
      std::random_device rd;
      m_random_engine.seed( rd() );
      m_initial_population_size = m_population_size;
    }

    explicit DifferentialEvolution( integer dimension ) : m_dimension( dimension )
    {
      std::random_device rd;
      m_random_engine.seed( rd() );
      m_initial_population_size = m_population_size;
    }

    // =========================================================================
    // Parameter setters
    // =========================================================================

    void
    set_dimension( integer dim )
    {
      m_dimension = dim;
    }

    void
    set_population_size( integer np )
    {
      m_population_size         = std::max( np, static_cast<integer>( 5 ) );
      m_initial_population_size = m_population_size;
    }

    void
    set_max_iterations( integer max_iter )
    {
      m_max_iterations = max_iter;
    }

    void
    set_strategy( Strategy strategy )
    {
      m_strategy = strategy;
    }

    void
    set_weight( real_type F )
    {
      m_weight = std::max<real_type>( 0.0, std::min<real_type>( 2.0, F ) );
    }

    void
    set_crossover_rate( real_type CR )
    {
      m_crossover_rate = std::max<real_type>( 0.0, std::min<real_type>( 1.0, CR ) );
    }

    void
    set_bounds( const Vector & lower, const Vector & upper )
    {
      if ( lower.size() != m_dimension || upper.size() != m_dimension )
      {
        throw std::invalid_argument( "Bound vectors must match dimension" );
      }
      m_lower_bounds = lower;
      m_upper_bounds = upper;
      m_use_bounds   = true;
    }

    void
    set_constraint_method( ConstraintMethod method )
    {
      m_constraint_method = method;
    }

    void
    set_tolerance( real_type tol )
    {
      m_tolerance = tol;
    }

    void
    set_value_to_reach( real_type vtr )
    {
      m_value_to_reach = vtr;
    }

    void
    set_verbose( bool verbose )
    {
      m_verbose = verbose;
    }

    void
    set_print_interval( integer interval )
    {
      m_print_interval = interval;
    }

    void
    set_seed( unsigned int seed )
    {
      m_random_engine.seed( seed );
    }

    void
    set_dither_range( real_type min, real_type max )
    {
      m_dither_min = min;
      m_dither_max = max;
    }

    void
    enable_parameter_adaptation( bool enable = true )
    {
      m_use_parameter_adaptation = enable;
      if ( enable )
      {
        m_successful_F.reserve( m_memory_size );
        m_successful_CR.reserve( m_memory_size );
      }
    }

    void
    enable_population_reduction( bool enable = true )
    {
      m_use_population_reduction = enable;
    }

    void
    enable_archive( bool enable = true )
    {
      m_use_archive = enable;
      if ( enable ) { m_archive.resize( m_dimension, m_archive_capacity ); }
    }

    void
    enable_parallel_evaluation( bool enable = true )
    {
      m_enable_parallel = enable;
    }

    void
    enable_local_search( bool enable = true )
    {
      m_use_local_search = enable;
    }

    void
    enable_fitness_cache( bool enable = true )
    {
      m_use_cache = enable;
    }

    void
    set_diversity_threshold( real_type threshold )
    {
      m_diversity_threshold = threshold;
    }

    void
    set_iteration_callback( std::function<void( integer, real_type )> callback )
    {
      m_iteration_callback = callback;
    }

    void
    set_new_best_callback( std::function<void( const Vector &, real_type )> callback )
    {
      m_new_best_callback = callback;
    }

    // =========================================================================
    // Getters
    // =========================================================================

    integer
    get_dimension() const
    {
      return m_dimension;
    }
    integer
    get_population_size() const
    {
      return m_population_size;
    }
    integer
    get_iteration() const
    {
      return m_current_iteration;
    }
    integer
    get_function_evaluations() const
    {
      return m_function_evaluations.load();
    }
    const Vector &
    get_best_solution() const
    {
      return m_best_solution;
    }
    real_type
    get_best_fitness() const
    {
      return m_best_fitness;
    }
    bool
    has_converged() const
    {
      return m_converged;
    }
    const std::string &
    get_convergence_reason() const
    {
      return m_convergence_reason;
    }
    const std::vector<real_type> &
    get_best_history() const
    {
      return m_best_history;
    }
    real_type
    get_diversity() const
    {
      return calculate_diversity();
    }

    // =========================================================================
    // Core algorithm
    // =========================================================================

    template <typename ObjectiveFunction>
    bool
    minimize( ObjectiveFunction && objective, Vector & result )
    {
      // Reset state
      initialize();
      m_best_history.clear();
      m_converged = false;
      m_convergence_reason.clear();

      // Initialize population (Dim x NP)
      initialize_population( objective );

      // Main optimization loop
      for ( m_current_iteration = 0; m_current_iteration < m_max_iterations; ++m_current_iteration )
      {
        // Generate trial population (Dim x NP)
        Matrix trial_population = m_population;
        Vector trial_fitness( m_population_size );

        // Generate trial vectors
        generate_trial_population( trial_population );

        // Evaluate trial population (possibly in parallel)
        evaluate_trial_population( trial_population, trial_fitness, objective );

        // Selection: replace if better
        for ( integer i = 0; i < m_population_size; ++i )
        {
          if ( trial_fitness[i] < m_fitness[i] )
          {
            // Add old individual to archive
            if ( m_use_archive ) { add_to_archive( m_population.col( i ) ); }

            // Update individual
            m_population.col( i ) = trial_population.col( i );
            m_fitness[i]          = trial_fitness[i];

            // Update parameter adaptation statistics
            if ( m_use_parameter_adaptation )
            {
              real_type F_used, CR_used;
              generate_parameters( F_used, CR_used );
              adapt_parameters( F_used, CR_used, true );
            }

            // Update best solution
            if ( trial_fitness[i] < m_best_fitness )
            {
              m_best_fitness  = trial_fitness[i];
              m_best_solution = trial_population.col( i );

              // Call new best callback
              if ( m_new_best_callback ) { m_new_best_callback( m_best_solution, m_best_fitness ); }
            }
          }
        }

        // Store best fitness for convergence checking
        m_best_history.push_back( m_best_fitness );

        // Reduce population size if enabled
        reduce_population();

        // Apply local search if enabled and near convergence
        if ( m_use_local_search && m_current_iteration > m_max_iterations * 0.8 )
        {
          real_type diversity = calculate_diversity();
          if ( diversity < m_local_search_threshold )
          {
            local_search( m_best_solution, m_best_fitness, objective, m_local_search_iterations );
          }
        }

        // Print progress
        if ( m_verbose &&
             ( m_current_iteration % m_print_interval == 0 || m_current_iteration == m_max_iterations - 1 ) )
        {
          print_progress();
        }

        // Call iteration callback
        if ( m_iteration_callback ) { m_iteration_callback( m_current_iteration, m_best_fitness ); }

        // Check convergence
        if ( check_convergence() )
        {
          m_converged = true;
          break;
        }
      }

      // Final local search if enabled
      if ( m_use_local_search )
      {
        local_search( m_best_solution, m_best_fitness, objective, m_local_search_iterations * 2 );
      }

      // Return result
      result = m_best_solution;
      return m_converged;
    }

    template <typename ObjectiveFunction>
    bool
    minimize( ObjectiveFunction && objective )
    {
      Vector result( m_dimension );
      return minimize( std::forward<ObjectiveFunction>( objective ), result );
    }

    // Benchmark functions for testing
    static real_type
    sphere_function( const Vector & x )
    {
      return x.squaredNorm();
    }

    static real_type
    rastrigin_function( const Vector & x )
    {
      const real_type A   = 10.0;
      real_type       sum = A * x.size();
      for ( integer i = 0; i < x.size(); ++i ) { sum += x[i] * x[i] - A * std::cos( 2.0 * M_PI * x[i] ); }
      return sum;
    }

    static real_type
    rosenbrock_function( const Vector & x )
    {
      real_type sum = 0.0;
      for ( integer i = 0; i < x.size() - 1; ++i )
      {
        sum += 100.0 * std::pow( x[i + 1] - x[i] * x[i], 2 ) + std::pow( 1.0 - x[i], 2 );
      }
      return sum;
    }

  private:
    // =========================================================================
    // Internal helper methods
    // =========================================================================

    void
    initialize()
    {
      m_current_iteration    = 0;
      m_function_evaluations = 0;
      m_best_fitness         = std::numeric_limits<real_type>::max();
      m_population.resize( m_dimension, m_population_size );  // Dim x NP
      m_fitness.resize( m_population_size );
      m_fitness.setConstant( std::numeric_limits<real_type>::max() );
      m_best_solution.resize( m_dimension );
      m_best_history.reserve( m_max_iterations );

      if ( m_use_parameter_adaptation )
      {
        m_successful_F.clear();
        m_successful_CR.clear();
      }

      if ( m_use_archive ) { m_current_archive_size = 0; }

      if ( m_use_cache ) { m_fitness_cache.clear(); }
    }

    template <typename ObjectiveFunction>
    void
    initialize_population( ObjectiveFunction && objective )
    {
      if ( m_use_bounds )
      {
        // Parallel initialization if enabled
        if ( m_enable_parallel )
        {
          for ( integer i = 0; i < m_population_size; ++i )
          {
            Vector individual( m_dimension );
            for ( integer j = 0; j < m_dimension; ++j )
            {
              individual[j] = m_lower_bounds[j] +
                              m_uniform_dist( m_random_engine ) * ( m_upper_bounds[j] - m_lower_bounds[j] );
            }

            real_type fitness = objective( individual );
            {
              m_population.col( i ) = individual;
              m_fitness[i]          = fitness;
              m_function_evaluations++;

              if ( fitness < m_best_fitness )
              {
                m_best_fitness  = fitness;
                m_best_solution = individual;
              }
            }
          }
        }
        else
        {
          for ( integer i = 0; i < m_population_size; ++i )
          {
            for ( integer j = 0; j < m_dimension; ++j )
            {
              m_population( j, i ) = m_lower_bounds[j] +
                                     m_uniform_dist( m_random_engine ) * ( m_upper_bounds[j] - m_lower_bounds[j] );
            }

            m_fitness[i] = objective( m_population.col( i ) );
            m_function_evaluations++;

            if ( m_fitness[i] < m_best_fitness )
            {
              m_best_fitness  = m_fitness[i];
              m_best_solution = m_population.col( i );
            }
          }
        }
      }
      else
      {
        // Generate in [-10, 10] range
        for ( integer i = 0; i < m_population_size; ++i )
        {
          for ( integer j = 0; j < m_dimension; ++j )
          {
            m_population( j, i ) = m_uniform_dist( m_random_engine ) * 20.0 - 10.0;
          }

          m_fitness[i] = objective( m_population.col( i ) );
          m_function_evaluations++;

          if ( m_fitness[i] < m_best_fitness )
          {
            m_best_fitness  = m_fitness[i];
            m_best_solution = m_population.col( i );
          }
        }
      }
    }

    void
    generate_trial_population( Matrix & trial_population )
    {
      // Prepare permutations
      std::vector<integer> perm1( m_population_size );
      std::vector<integer> perm2( m_population_size );
      std::vector<integer> perm3( m_population_size );

      generate_permutation( perm1 );
      generate_permutation( perm2 );
      generate_permutation( perm3 );

      // Prepare dither factor if needed
      real_type dither_factor = 1.0;
      if ( m_strategy == RAND_1_DITHER )
      {
        dither_factor = m_dither_min + m_uniform_dist( m_random_engine ) * ( m_dither_max - m_dither_min );
      }

      // Generate trial vectors for each individual
      for ( integer i = 0; i < m_population_size; ++i )
      {
        // Get random indices (ensure they're different from i)
        integer r1 = perm1[i];
        integer r2 = perm2[i];
        integer r3 = perm3[i];

        while ( r1 == i ) r1 = ( r1 + 1 ) % m_population_size;
        while ( r2 == i ) r2 = ( r2 + 1 ) % m_population_size;
        while ( r3 == i ) r3 = ( r3 + 1 ) % m_population_size;

        // Get the individuals (columns)
        Vector ind_r1 = m_population.col( r1 );
        Vector ind_r2 = m_population.col( r2 );
        Vector ind_r3 = m_population.col( r3 );

        // Get parameters (adapted or fixed)
        real_type F, CR;
        generate_parameters( F, CR );

        // Generate mutation vector based on strategy
        Vector mutation_vector( m_dimension );
        Vector original_individual = m_population.col( i );

        switch ( m_strategy )
        {
          case RAND_1:  // DE/rand/1
            mutation_vector = ind_r1 + F * ( ind_r2 - ind_r3 );
            break;

          case BEST_1:  // DE/best/1
            mutation_vector = m_best_solution + F * ( ind_r1 - ind_r2 );
            break;

          case RAND_TO_BEST_1:  // DE/rand-to-best/1
            mutation_vector = ind_r1 + F * ( m_best_solution - ind_r1 ) + F * ( ind_r2 - ind_r3 );
            break;

          case BEST_2:  // DE/best/2
          {
            integer r4      = ( r3 + 1 ) % m_population_size;
            Vector  ind_r4  = m_population.col( r4 );
            mutation_vector = m_best_solution + F * ( ind_r1 - ind_r2 ) + F * ( ind_r3 - ind_r4 );
          }
          break;

          case RAND_2:  // DE/rand/2
          {
            integer r4      = ( r3 + 1 ) % m_population_size;
            integer r5      = ( r4 + 1 ) % m_population_size;
            Vector  ind_r4  = m_population.col( r4 );
            Vector  ind_r5  = m_population.col( r5 );
            mutation_vector = ind_r1 + F * ( ind_r2 - ind_r3 ) + F * ( ind_r4 - ind_r5 );
          }
          break;

          case RAND_1_DITHER:  // DE/rand/1 with dither
            mutation_vector = ind_r1 + dither_factor * ( ind_r2 - ind_r3 );
            break;

          case EITHER_OR:  // Either-or algorithm
            if ( m_uniform_dist( m_random_engine ) < 0.5 ) { mutation_vector = ind_r1 + F * ( ind_r2 - ind_r3 ); }
            else
            {
              mutation_vector = ind_r1 + 0.5 * ( F + 1.0 ) * ( ind_r2 + ind_r3 - 2.0 * ind_r1 );
            }
            break;

          case CURRENT_TO_PBEST:  // DE/current-to-pbest/1 (JADE)
          {
            // Select pbest individual (top p% of population)
            integer              p = std::max<integer>( 2, m_population_size / 10 );
            std::vector<integer> indices( p );
            std::iota( indices.begin(), indices.end(), 0 );
            std::partial_sort( indices.begin(), indices.begin() + p, indices.end(),
                               [this]( integer a, integer b ) { return m_fitness[a] < m_fitness[b]; } );

            integer pbest_idx = indices[static_cast<integer>( m_uniform_dist( m_random_engine ) * p )];
            Vector  ind_pbest = m_population.col( pbest_idx );

            mutation_vector = original_individual + F * ( ind_pbest - original_individual ) + F * ( ind_r1 - ind_r2 );
          }
          break;

          case CURRENT_TO_PBEST_WITH_ARCHIVE:  // JADE with archive
          {
            integer              p = std::max<integer>( 2, m_population_size / 10 );
            std::vector<integer> indices( p );
            std::iota( indices.begin(), indices.end(), 0 );
            std::partial_sort( indices.begin(), indices.begin() + p, indices.end(),
                               [this]( integer a, integer b ) { return m_fitness[a] < m_fitness[b]; } );

            integer pbest_idx = indices[static_cast<integer>( m_uniform_dist( m_random_engine ) * p )];
            Vector  ind_pbest = m_population.col( pbest_idx );

            // Select from archive if available
            Vector ind_r2_used = ind_r2;
            if ( m_use_archive && m_current_archive_size > 0 )
            {
              if ( m_uniform_dist( m_random_engine ) < 0.5 )
              {
                integer archive_idx = static_cast<integer>( m_uniform_dist( m_random_engine ) *
                                                            m_current_archive_size );
                ind_r2_used         = m_archive.col( archive_idx );
              }
            }

            mutation_vector = original_individual + F * ( ind_pbest - original_individual ) +
                              F * ( ind_r1 - ind_r2_used );
          }
          break;

          default:
            mutation_vector = ind_r1 + F * ( ind_r2 - ind_r3 );
            break;
        }

        // Apply boundary constraints to mutation vector
        if ( m_use_bounds ) { apply_boundary_constraints( mutation_vector, original_individual ); }

        // Binomial crossover
        integer j_rand = static_cast<integer>( m_uniform_dist( m_random_engine ) * m_dimension );

        for ( integer j = 0; j < m_dimension; ++j )
        {
          if ( j == j_rand || m_uniform_dist( m_random_engine ) < CR )
          {
            trial_population( j, i ) = mutation_vector[j];
          }
          else
          {
            trial_population( j, i ) = original_individual[j];
          }
        }
      }
    }

    template <typename ObjectiveFunction>
    void
    evaluate_trial_population( Matrix & trial_population, Vector & trial_fitness, ObjectiveFunction && objective )
    {
      if ( m_enable_parallel )
      {
        for ( integer i = 0; i < m_population_size; ++i )
        {
          Vector individual = trial_population.col( i );

          // Check cache first if enabled
          if ( m_use_cache )
          {
            size_t                      hash = hash_vector( individual );
            std::lock_guard<std::mutex> lock( m_cache_mutex );
            auto                        it = m_fitness_cache.find( hash );
            if ( it != m_fitness_cache.end() )
            {
              trial_fitness[i] = it->second;
              continue;
            }
          }

          real_type fitness = objective( individual );

          // Update cache if enabled
          if ( m_use_cache )
          {
            std::lock_guard<std::mutex> lock( m_cache_mutex );
            m_fitness_cache[hash_vector( individual )] = fitness;
          }

          trial_fitness[i] = fitness;
          m_function_evaluations++;
        }
      }
      else
      {
        for ( integer i = 0; i < m_population_size; ++i )
        {
          Vector individual = trial_population.col( i );

          if ( m_use_cache )
          {
            size_t hash = hash_vector( individual );
            auto   it   = m_fitness_cache.find( hash );
            if ( it != m_fitness_cache.end() )
            {
              trial_fitness[i] = it->second;
              continue;
            }
          }

          trial_fitness[i] = objective( individual );

          if ( m_use_cache ) { m_fitness_cache[hash_vector( individual )] = trial_fitness[i]; }

          m_function_evaluations++;
        }
      }
    }

    bool
    check_convergence()
    {
      // Check value to reach
      if ( m_best_fitness <= m_value_to_reach )
      {
        m_convergence_reason = "Reached target value";
        return true;
      }

      // Check stagnation
      if ( m_best_history.size() > static_cast<size_t>( m_stagnation_window ) )
      {
        real_type best_recent = *std::min_element( m_best_history.end() - m_stagnation_window, m_best_history.end() );
        real_type improvement = m_best_history[m_best_history.size() - m_stagnation_window] - best_recent;

        if ( std::abs( improvement ) < m_tolerance )
        {
          m_convergence_reason = "Stagnation detected";
          return true;
        }
      }

      // Check diversity
      real_type diversity = calculate_diversity();
      if ( diversity < m_diversity_threshold )
      {
        m_convergence_reason = fmt::format( "Low diversity: {:.2e}", diversity );
        return true;
      }

      return false;
    }

    void
    print_progress()
    {
      real_type diversity = calculate_diversity();

      fmt::print(
          "DE Iteration {:5}: Best = {:.6e}, F = {:.3f}, CR = {:.3f}, NP = "
          "{:3}, Div = {:.2e}\n",
          m_current_iteration, m_best_fitness, m_weight, m_crossover_rate, m_population_size, diversity );

      if ( m_verbose && m_current_iteration % ( m_print_interval * 5 ) == 0 )
      {
        fmt::print( "  Best solution: [" );
        for ( integer j = 0; j < std::min( m_dimension, static_cast<integer>( 5 ) ); ++j )
        {
          fmt::print( "{:.6e}", m_best_solution[j] );
          if ( j < std::min( m_dimension, static_cast<integer>( 5 ) ) - 1 ) { fmt::print( ", " ); }
        }
        if ( m_dimension > 5 ) { fmt::print( ", ..." ); }
        fmt::print( "]\n" );
      }
    }
  };

}  // namespace Utils

#endif  // UTILS_DIFFERENTIAL_EVOLUTION_dot_HH
