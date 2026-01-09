/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022-2024                                                 |
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

/**
 * \file Utils_HJPatternSearch.hh
 * \brief Header-only implementation of Hooke-Jeeves Pattern Search algorithm
 *
 * \details This file provides a modern C++ implementation of the Hooke-Jeeves
 * direct search method for derivative-free optimization. The implementation
 * is header-only, uses Eigen for linear algebra operations, and {fmt} for
 * formatted output.
 *
 * \note Unicode symbols are used for enhanced console output visualization.
 *
 * \ingroup Minimize
 *
 * \author Enrico Bertolazzi
 * \version 1.0.0
 * \date 2022-2024
 */

#ifndef UTILS_HJ_PATTERN_SEARCH_HH
#define UTILS_HJ_PATTERN_SEARCH_HH

#pragma once

#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

// Required headers
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace Utils
{

  // Forward declaration
  class Console;

  /*!
   * \addtogroup Minimize
   * \brief Optimization algorithms and methods
   * @{
   */

  /**
   * \class HJPatternSearch
   * \brief Implementation of Hooke-Jeeves Pattern Search algorithm for derivative-free optimization.
   *
   * \tparam Real Floating-point type (float, double, etc.)
   *
   * \details The Hooke-Jeeves algorithm is a direct search method that doesn't require
   * gradient information. It combines exploratory moves (searching around the current point)
   * with pattern moves (accelerating in promising directions).
   *
   * ## Algorithm Overview
   *
   * 1. **Exploratory Move**: Search in coordinate directions around the current point
   * 2. **Pattern Move**: Accelerate in the direction of improvement
   * 3. **Step Reduction**: Reduce step size when no improvement is found
   * 4. **Termination**: Stop when step size is below tolerance or max iterations reached
   *
   * ## Features
   *
   * - Header-only implementation
   * - Uses Eigen for efficient vector operations
   * - Configurable tolerance, iterations, and function evaluations
   * - Verbosity control for debugging
   * - Stagnation detection
   *
   * ## Complexity
   *
   * - Memory: O(n) where n is the problem dimension
   * - Function evaluations: Variable, typically O(n * iterations)
   *
   * ## References
   *
   * 1. Hooke, R., & Jeeves, T. A. (1961). "Direct Search" Solution of Numerical
   *    and Statistical Problems. Journal of the ACM, 8(2), 212-229.
   * 2. Lewis, R. M., Torczon, V., & Trosset, M. W. (2000). Direct search methods:
   *    then and now. Journal of Computational and Applied Mathematics, 124(1-2), 191-207.
   * 3. Kolda, T. G., Lewis, R. M., & Torczon, V. (2003). Optimization by direct search:
   *    New perspectives on some classical and modern methods. SIAM Review, 45(3), 385-482.
   *
   * ## Usage Example
   *
   * \code{.cpp}
   * #include "Utils_HJPatternSearch.hh"
   *
   * using Real = double;
   *
   * // Rosenbrock function
   * Real rosenbrock(const Real* x) {
   *     return 100.0 * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1.0 - x[0], 2);
   * }
   *
   * int main() {
   *     Utils::Console console(&std::cout, 2); // Verbosity level 2
   *     Utils::HJPatternSearch<Real> solver("RosenbrockSolver");
   *
   *     std::function<Real(Real const[])> f = rosenbrock;
   *     solver.setup(2, f, &console);
   *
   *     solver.set_tolerance(1e-6);
   *     solver.set_max_iterations(1000);
   *     solver.set_max_fun_evaluation(5000);
   *
   *     Real x0[2] = {-1.2, 1.0};
   *     solver.run(x0, 0.5);
   *
   *     Real solution[2];
   *     Real fval = solver.get_last_solution(solution);
   *
   *     fmt::print("Solution: f({:.6f}, {:.6f}) = {:.6f}\n",
   *                solution[0], solution[1], fval);
   *     solver.print_info(std::cout);
   *
   *     return 0;
   * }
   * \endcode
   */
  template <typename Real> class HJPatternSearch
  {
  public:
    using Vector         = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using MapVector      = Eigen::Map<Vector>;
    using ConstMapVector = Eigen::Map<const Vector>;
    using integer        = int;
    using Function       = std::function<Real( const Real[] )>;

  private:
    // Algorithm parameters
    std::string const m_name;                ///< Name of the solver instance
    Function          m_fun;                 ///< Objective function to minimize
    Console const *   m_console{ nullptr };  ///< Console for output messages
    integer           m_verbose{ 1 };        ///< Verbosity level (0=silent, 1=basic, 2=detailed, 3=debug)

    // Algorithm constants
    Real    m_rho{ static_cast<Real>( 0.9 ) };         ///< Step reduction factor (0 < ρ < 1)
    Real    m_h{ static_cast<Real>( 0.1 ) };           ///< Current step size
    Real    m_tolerance{ static_cast<Real>( 1e-8 ) };  ///< Convergence tolerance
    integer m_dim{ 0 };                                ///< Problem dimension

    // Counters and limits
    integer m_max_iterations{ 500 };        ///< Maximum number of iterations
    integer m_max_fun_evaluations{ 1000 };  ///< Maximum function evaluations
    integer m_max_stagnations{ 10 };        ///< Maximum stagnations before stopping
    integer m_iteration_count{ 0 };         ///< Current iteration count
    integer m_fun_evaluation_count{ 0 };    ///< Function evaluation count
    integer m_stagnation_count{ 0 };        ///< Stagnation counter

    // State variables
    bool m_stencil_failure{ false };  ///< Flag indicating stencil failure
    Real m_f_best;                    ///< Best function value found
    Real m_f_old;                     ///< Previous function value

    // Storage vectors (allocated on heap)
    std::unique_ptr<Real[]> m_storage;                    ///< Raw storage for all vectors
    MapVector               m_x_best{ nullptr, 0 };       ///< Best point found
    MapVector               m_x_old{ nullptr, 0 };        ///< Previous point
    MapVector               m_direction{ nullptr, 0 };    ///< Search direction for pattern move
    MapVector               m_search_sign{ nullptr, 0 };  ///< Sign for coordinate search directions
    MapVector               m_temp_point1{ nullptr, 0 };  ///< Temporary point 1
    MapVector               m_temp_point2{ nullptr, 0 };  ///< Temporary point 2
    MapVector               m_new_point{ nullptr, 0 };    ///< New point for pattern move

    /**
     * \brief Evaluate the objective function at a point
     *
     * \param x Point at which to evaluate function
     * \return Function value at x
     *
     * \note Increments the function evaluation counter
     */
    Real evaluate_function( const Real * x ) const
    {
      ++const_cast<HJPatternSearch *>( this )->m_fun_evaluation_count;
      return m_fun( x );
    }

    /**
     * \brief Perform exploratory move (coordinate search)
     *
     * Searches along coordinate directions around the current best point.
     * Updates the best point if improvement is found.
     */
    void exploratory_move()
    {
      m_stencil_failure = true;  // Assume failure until improvement found

      for ( integer j = 0; j < m_dim; ++j )
      {
        Real step = m_search_sign( j ) * m_h;

        // Try positive direction
        m_temp_point1 = m_x_best;
        m_temp_point1( j ) += step;
        Real f1 = evaluate_function( m_temp_point1.data() );

        // Try negative direction
        m_temp_point2 = m_x_best;
        m_temp_point2( j ) -= step;
        Real f2 = evaluate_function( m_temp_point2.data() );

        // Choose the better direction
        if ( f1 < m_f_best && f1 <= f2 )
        {
          m_x_best          = m_temp_point1;
          m_f_best          = f1;
          m_stencil_failure = false;
        }
        else if ( f2 < m_f_best )
        {
          m_x_best           = m_temp_point2;
          m_f_best           = f2;
          m_search_sign( j ) = -m_search_sign( j );  // Flip search direction
          m_stencil_failure  = false;
        }

        // Check evaluation limit
        if ( m_fun_evaluation_count >= m_max_fun_evaluations ) { break; }
      }
    }

    /**
     * \brief Perform pattern move (acceleration)
     *
     * Accelerates in the direction of improvement found by exploratory move.
     * Uses line search with backtracking.
     */
    void pattern_move()
    {
      // Compute pattern direction
      m_direction = m_x_best - m_x_old;

      if ( m_direction.norm() < m_tolerance )
      {
        return;  // No significant movement
      }

      // Line search with backtracking
      Real lambda          = static_cast<Real>( 1.0 );
      Real max_improvement = static_cast<Real>( 0.0 );

      while ( m_fun_evaluation_count < m_max_fun_evaluations && lambda > static_cast<Real>( 0.1 ) )
      {
        m_new_point = m_x_best + lambda * m_direction;
        Real f_new  = evaluate_function( m_new_point.data() );

        // Armijo-like condition
        if ( f_new < m_f_best - static_cast<Real>( 0.25 ) * lambda * max_improvement )
        {
          Real improvement = ( m_f_best - f_new ) / lambda;
          if ( improvement > max_improvement ) { max_improvement = improvement; }
          m_x_best = m_new_point;
          m_f_best = f_new;
          // Could expand step here: lambda *= 2;
        }

        lambda *= static_cast<Real>( 0.5 );  // Reduce step
      }
    }

  public:
    /**
     * \brief Construct a new HJPatternSearch object
     *
     * \param name Name identifier for the solver
     */
    explicit HJPatternSearch( std::string_view name ) : m_name( name ) {}

    /**
     * \brief Destroy the HJPatternSearch object
     */
    ~HJPatternSearch() = default;

    // Delete copy constructor and assignment
    HJPatternSearch( const HJPatternSearch & )             = delete;
    HJPatternSearch & operator=( const HJPatternSearch & ) = delete;

    /**
     * \brief Get the solver name
     *
     * \return std::string_view Name of the solver
     */
    std::string_view name() const { return m_name; }

    /**
     * \brief Set up the solver with problem dimension and objective function
     *
     * \param dim Problem dimension (number of variables)
     * \param fun Objective function to minimize
     * \param console Console for output (can be nullptr)
     */
    void setup( integer dim, Function & fun, Console const * console = nullptr )
    {
      UTILS_ASSERT( dim > 0, "HJPatternSearch::setup: dimension {} must be > 0", dim );

      m_dim     = dim;
      m_fun     = fun;
      m_console = console;

      // Allocate storage for all vectors (7 vectors of size dim)
      m_storage = std::make_unique<Real[]>( 7 * dim );

      // Map Eigen vectors to storage
      Real * data = m_storage.get();
      new ( &m_x_best ) MapVector( data, dim );
      data += dim;
      new ( &m_x_old ) MapVector( data, dim );
      data += dim;
      new ( &m_direction ) MapVector( data, dim );
      data += dim;
      new ( &m_search_sign ) MapVector( data, dim );
      data += dim;
      new ( &m_temp_point1 ) MapVector( data, dim );
      data += dim;
      new ( &m_temp_point2 ) MapVector( data, dim );
      data += dim;
      new ( &m_new_point ) MapVector( data, dim );

      // Initialize search signs to positive
      m_search_sign.setConstant( static_cast<Real>( 1.0 ) );
    }

    /**
     * \brief Change the console for output messages
     *
     * \param console New console to use
     */
    void change_console( Console const * console ) { m_console = console; }

    /**
     * \brief Set verbosity level
     *
     * \param level Verbosity level:
     *              - 0: No output
     *              - 1: Basic information (iterations, function values)
     *              - 2: Detailed information
     *              - 3: Debug information (vectors, etc.)
     */
    void set_verbose( integer level ) { m_verbose = level; }

    /**
     * \brief Set convergence tolerance
     *
     * \param tol Tolerance value (must be > 0)
     */
    void set_tolerance( Real tol )
    {
      UTILS_ASSERT( tol > 0, "HJPatternSearch::set_tolerance: tolerance {} must be > 0", tol );
      m_tolerance = tol;
    }

    /**
     * \brief Set maximum number of iterations
     *
     * \param max_iter Maximum iterations (must be > 0)
     */
    void set_max_iterations( integer max_iter )
    {
      UTILS_ASSERT( max_iter > 0, "HJPatternSearch::set_max_iterations: max_iter {} must be > 0", max_iter );
      m_max_iterations = max_iter;
    }

    /**
     * \brief Set maximum number of function evaluations
     *
     * \param max_fev Maximum function evaluations (must be > 0)
     */
    void set_max_fun_evaluations( integer max_fev )
    {
      UTILS_ASSERT( max_fev > 0, "HJPatternSearch::set_max_fun_evaluations: max_fev {} must be > 0", max_fev );
      m_max_fun_evaluations = max_fev;
    }

    /**
     * \brief Set maximum number of stagnations
     *
     * \param max_stag Maximum stagnations before termination (must be > 0)
     */
    void set_max_stagnations( integer max_stag )
    {
      UTILS_ASSERT( max_stag > 0, "HJPatternSearch::set_max_stagnations: max_stag {} must be > 0", max_stag );
      m_max_stagnations = max_stag;
    }

    /**
     * \brief Set step reduction factor
     *
     * \param rho Reduction factor (must satisfy 0 < rho < 1)
     */
    void set_rho( Real rho )
    {
      UTILS_ASSERT( rho > 0 && rho < 1, "HJPatternSearch::set_rho: rho {} must be in (0,1)", rho );
      m_rho = rho;
    }

    /**
     * \brief Get information about the optimization result
     *
     * \return std::string Formatted string with optimization statistics
     */
    std::string info() const
    {
      std::string result = "\n";
      result += "╔══════════════════════════════════════════════════════════════╗\n";
      result += fmt::format( "║ Hooke-Jeeves Pattern Search Results {:>25} ║\n", " " );
      result += "╠══════════════════════════════════════════════════════════════╣\n";

      // Termination reason
      if ( m_h <= m_tolerance ) { result += "║ ✓ Converged: Step size below tolerance\n"; }
      else if ( m_iteration_count >= m_max_iterations ) { result += "║ ⚠ Stopped: Maximum iterations reached\n"; }
      else if ( m_fun_evaluation_count >= m_max_fun_evaluations )
      {
        result += "║ ⚠ Stopped: Maximum function evaluations reached\n";
      }
      else if ( m_stagnation_count >= m_max_stagnations ) { result += "║ ⚠ Stopped: Maximum stagnations reached\n"; }
      else
      {
        result += "║ ? Unknown termination reason\n";
      }

      // Statistics
      result += fmt::format( "║   Iterations: {:5d} / {:5d}\n", m_iteration_count, m_max_iterations );
      result +=
        fmt::format( "║   Function evaluations: {:5d} / {:5d}\n", m_fun_evaluation_count, m_max_fun_evaluations );
      result += fmt::format( "║   Final step size: {:.3e} (tol = {:.1e})\n", m_h, m_tolerance );
      result += fmt::format( "║   Best function value: {:.10e}\n", m_f_best );

      // Solution (truncated if many dimensions)
      result += "║   Solution: [";
      integer show_dims = std::min( m_dim, 5 );
      for ( integer i = 0; i < show_dims; ++i )
      {
        result += fmt::format( "{:.6e}", m_x_best( i ) );
        if ( i < show_dims - 1 ) result += ", ";
      }
      if ( m_dim > show_dims ) { result += fmt::format( ", ... (+{} more)", m_dim - show_dims ); }
      result += "]\n";

      result += "╚══════════════════════════════════════════════════════════════╝\n";
      return result;
    }

    /**
     * \brief Print optimization information to output stream
     *
     * \param os Output stream
     */
    void print_info( std::ostream & os ) const { os << info(); }

    /**
     * \brief Run the Hooke-Jeeves optimization algorithm
     *
     * \param x0 Initial point (array of size dim)
     * \param h0 Initial step size
     *
     * \note The algorithm stops when:
     *       - Step size h < tolerance
     *       - Maximum iterations reached
     *       - Maximum function evaluations reached
     *       - Maximum stagnations reached
     */
    void run( const Real x0[], Real h0 )
    {
      // Reset counters
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_stagnation_count     = 0;

      // Initialize
      m_h = h0;
      std::copy_n( x0, m_dim, m_x_best.data() );
      m_f_best = evaluate_function( m_x_best.data() );
      m_search_sign.setConstant( static_cast<Real>( 1.0 ) );

      // Main optimization loop
      while ( m_h > m_tolerance && m_iteration_count < m_max_iterations &&
              m_fun_evaluation_count < m_max_fun_evaluations && m_stagnation_count < m_max_stagnations )
      {
        ++m_iteration_count;

        // Print iteration info
        if ( m_verbose > 0 && m_console )
        {
          std::string line = "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━";
          std::string msg  = fmt::format(
            "⎧ Iteration {:4d}: f = {:.6e}, feval = {:4d}, h = {:.3e}\n",
            m_iteration_count,
            m_f_best,
            m_fun_evaluation_count,
            m_h );

          if ( m_verbose > 2 )
          {
            for ( integer i = 0; i < m_dim; ++i ) { msg += fmt::format( "⎪ x[{:2d}] = {:.6e}\n", i, m_x_best( i ) ); }
          }
          msg += fmt::format( "⎩{}\n", line );
          // Assuming Console has a log method
          // m_console->log(msg);
        }

        // Save old point
        m_x_old = m_x_best;
        m_f_old = m_f_best;

        // Perform exploratory move
        exploratory_move();

        // Check for stencil failure
        while ( m_stencil_failure && m_h > m_tolerance )
        {
          m_h *= m_rho;
          exploratory_move();
        }

        // If exploratory move succeeded, try pattern move
        if ( !m_stencil_failure ) { pattern_move(); }

        // Reduce step size
        m_h *= m_rho;

        // Check stagnation
        if ( std::abs( m_f_best - m_f_old ) < m_tolerance ) { ++m_stagnation_count; }
        else
        {
          m_stagnation_count = 0;
        }
      }
    }

    /**
     * \brief Get the best solution found
     *
     * \param x Output array for the solution (must have size ≥ dim)
     * \return Real Function value at the best solution
     */
    Real get_last_solution( Real x[] ) const
    {
      std::copy_n( m_x_best.data(), m_dim, x );
      return m_f_best;
    }

    /**
     * \brief Get the best function value found
     *
     * \return Real Best function value
     */
    Real get_best_value() const { return m_f_best; }

    /**
     * \brief Get iteration count
     *
     * \return integer Number of iterations performed
     */
    integer get_iteration_count() const { return m_iteration_count; }

    /**
     * \brief Get function evaluation count
     *
     * \return integer Number of function evaluations performed
     */
    integer get_fun_evaluation_count() const { return m_fun_evaluation_count; }

    /**
     * \brief Get current step size
     *
     * \return Real Current step size h
     */
    Real get_current_step() const { return m_h; }
  };

  /*! @} */  // end of Minimize group

}  // namespace Utils

#endif  // UTILS_HJ_PATTERN_SEARCH_HH
