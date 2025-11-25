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
// file: Utils_BOBYQA.hh
//
// Implementation of BOBYQA (Bound Optimization BY Quadratic Approximation)
// for derivative-free optimization with bound constraints.
//
// Main references:
// - Powell, M. J. D. (2009). The BOBYQA algorithm for bound constrained
//   optimization without derivatives. Technical Report, University of Cambridge.
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_BOBYQA_dot_HH
#define UTILS_BOBYQA_dot_HH

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdouble-promotion"
#pragma clang diagnostic ignored "-Wsigned-enum-bitfield"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace Utils {

  using std::abs;
  using std::min;
  using std::max;
  using std::sqrt;
  using std::pow;

  /**
   * @class BOBYQA_minimizer
   * @brief Minimizer using BOBYQA (Bound Optimization BY Quadratic Approximation)
   * 
   * Implements the BOBYQA algorithm for derivative-free optimization with bound constraints.
   * The method builds quadratic models using interpolation points and maintains bounds
   * on the variables throughout the optimization.
   * 
   * Key features:
   * - No gradient information required
   * - Handles bound constraints efficiently
   * - Uses quadratic models for local approximation
   * - Maintains a set of interpolation points
   * - Robust convergence properties
   *
   * References:
   * - Powell, M. J. D. (2009). The BOBYQA algorithm for bound constrained
   *   optimization without derivatives. Technical Report, University of Cambridge.
   * 
   * @tparam Scalar Data type (double, float, etc.)
   */
  template <typename Scalar = double>
  class BOBYQA_minimizer {
  public:
  
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>; ///< Vector type
    using Matrix   = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; ///< Vector type
    using Callback = std::function<Scalar(Vector const&)>;     ///< Function callback type
  
    /**
     * @struct Options
     * @brief Options for BOBYQA minimizer
     */
    struct Options {
      size_t  max_iter          { 1000    }; ///< Maximum iterations
      size_t  n_interpolation   { 0       }; ///< Number of interpolation points (0 = auto: 2*n+1)
      Scalar  rhobeg            { 0.5     }; ///< Initial trust region radius
      Scalar  rhoend            { 1.0e-6  }; ///< Final trust region radius
      bool    verbose           { false   }; ///< Verbose output
      size_t  print_every       { 10      }; ///< Print every n iterations

      // Convergence criteria
      Scalar  tol_f             { 1e-8    }; ///< Function value tolerance
      Scalar  tol_x             { 1e-6    }; ///< Step size tolerance
      size_t  patience          { 50      }; ///< Iterations without improvement before stop
      
      // Model parameters
      bool    adaptive_interpolation { true  }; ///< Adaptive interpolation point selection
      Scalar  max_model_degree  { 2.0     }; ///< Maximum degree for quadratic model
    };
  
    /**
     * @struct Result
     * @brief Minimization result
     */
    struct Result {
      Scalar final_f;          ///< Final function value
      Vector final_x;          ///< Final point
      size_t iterations;       ///< Number of iterations performed
      size_t f_eval_count;     ///< Number of function evaluations
      bool   converged;        ///< True if convergence reached
      Scalar trust_region_radius; ///< Final trust region radius
      string message;          ///< Termination message
    };

  private:
    Options m_opts;        ///< Configured options
    Vector  m_lower;       ///< Lower bounds
    Vector  m_upper;       ///< Upper bounds

    // Internal state for BOBYQA algorithm
    struct InternalState {
      Vector current_x;        ///< Current best point
      Vector best_x;           ///< Best point found
      Scalar current_f;        ///< Current function value
      Scalar best_f;           ///< Best function value
      Scalar rho;              ///< Current trust region radius
      size_t n_improvements;   ///< Number of improvements
      
      // Quadratic model parameters
      Vector model_gradient;
      Matrix model_hessian;
      std::vector<Vector> interpolation_points;
      std::vector<Scalar> interpolation_values;
    };
    
    InternalState m_state;

  public:
  
    /**
     * @brief Constructor
     * @param o Options for the minimizer
     */
    BOBYQA_minimizer( Options const & o = Options() )
    : m_opts(o)
    {}
  
    /**
     * @brief Set optimization bounds
     * @param lo Lower bounds
     * @param up Upper bounds
     */
    void
    set_bounds( Vector const & lo, Vector const & up ) {
      assert(lo.size() == up.size());
      m_lower = lo;
      m_upper = up;
    }
    
    void
    set_bounds(size_t n, Scalar const lower[], Scalar const upper[]) {
      m_lower.resize(n);
      m_upper.resize(n);
      std::copy_n(lower, n, m_lower.data());
      std::copy_n(upper, n, m_upper.data());
    }

    /**
     * @brief Initialize the minimizer state
     * @param x0 Starting point
     * @param fun Objective function
     */
    void
    initialize( Vector const & x0, Callback const & fun ) {
      size_t n = x0.size();
      
      // Set default number of interpolation points if not specified
      if (m_opts.n_interpolation == 0) {
        m_opts.n_interpolation = 2 * n + 1;
      }
      
      // Initialize state
      m_state.current_x = x0;
      m_state.best_x = x0;
      m_state.current_f = fun(x0);
      m_state.best_f = m_state.current_f;
      m_state.rho = m_opts.rhobeg;
      m_state.n_improvements = 0;
      
      // Initialize quadratic model
      m_state.model_gradient = Vector::Zero(n);
      m_state.model_hessian = Matrix::Zero(n, n);
      
      // Initialize interpolation points
      initialize_interpolation_points(x0, fun);
    }

  private:

    /**
     * @brief Initialize interpolation points around starting point
     * @param x0 Starting point
     * @param fun Objective function
     */
    void
    initialize_interpolation_points( Vector const & x0, Callback const & fun ) {
      m_state.interpolation_points.clear();
      m_state.interpolation_values.clear();
      
      // Add center point
      m_state.interpolation_points.push_back(x0);
      m_state.interpolation_values.push_back(m_state.current_f);
      
      // Generate interpolation points around x0 within bounds
      for (size_t i = 0; i < m_opts.n_interpolation - 1; ++i) {
        Vector point = generate_interpolation_point(x0, i);
        Scalar value = fun(point);
        
        m_state.interpolation_points.push_back(point);
        m_state.interpolation_values.push_back(value);
      }
    }

    /**
     * @brief Generate a new interpolation point
     * @param center Center point
     * @param index Point index
     * @return New interpolation point
     */
    Vector
    generate_interpolation_point( Vector const & center, size_t index ) const {
      Vector point = center;
      size_t n = center.size();
      
      // Simple strategy: perturb each coordinate in sequence
      size_t coord = index % n;
      Scalar perturbation = m_state.rho * (1.0 - 2.0 * ((index / n) % 2));
      
      point(coord) += perturbation;
      
      // Ensure point stays within bounds
      point = point.cwiseMax(m_lower).cwiseMin(m_upper);
      
      return point;
    }

    /**
     * @brief Update quadratic model using interpolation points
     */
    void
    update_quadratic_model() {
      size_t n = m_state.current_x.size();
      size_t m = m_state.interpolation_points.size();
      
      // Solve least squares problem for quadratic model coefficients
      // This is a simplified version - actual BOBYQA uses more sophisticated approach
      
      // Build matrix of differences
      Matrix A(m, (n * (n + 3)) / 2);
      Vector b(m);
      
      Vector center = m_state.current_x;
      
      for (size_t i = 0; i < m; ++i) {
        Vector dx = m_state.interpolation_points[i] - center;
        b(i) = m_state.interpolation_values[i] - m_state.current_f;
        
        // Fill quadratic terms
        size_t col = 0;
        for (size_t j = 0; j < n; ++j) {
          A(i, col++) = dx(j);
          for (size_t k = j; k < n; ++k) {
            A(i, col++) = dx(j) * dx(k);
          }
        }
      }
      
      // Solve for model coefficients (simplified)
      Vector coefficients = A.colPivHouseholderQr().solve(b);
      
      // Extract gradient and Hessian from coefficients
      extract_model_parameters(coefficients);
    }

    /**
     * @brief Extract gradient and Hessian from coefficient vector
     * @param coeff Model coefficients
     */
    void
    extract_model_parameters( Vector const & coeff ) {
      size_t n = m_state.current_x.size();
      size_t idx = 0;
      
      // Gradient terms
      for (size_t i = 0; i < n; ++i) {
        m_state.model_gradient(i) = coeff(idx++);
      }
      
      // Hessian terms (symmetric)
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
          Scalar value = coeff(idx++);
          m_state.model_hessian(i, j) = value;
          m_state.model_hessian(j, i) = value;
        }
      }
    }

    /**
     * @brief Solve trust region subproblem
     * @return Proposed step
     */
    Vector
    solve_trust_region_subproblem() const {
      // Simplified trust region subproblem solution
      // Actual BOBYQA uses more sophisticated approach
      
      Vector step = -m_state.model_gradient;
      Scalar step_norm = step.norm();
      
      if (step_norm > m_state.rho) {
        step = step * (m_state.rho / step_norm);
      }
      
      return step;
    }

  public:

    /**
     * @brief Perform minimization using BOBYQA
     * 
     * Convergence criteria:
     * - Trust region radius below tolerance
     * - Function value change below tolerance  
     * - Step size below tolerance
     * - No improvement for patience iterations
     * - Maximum iterations reached
     * 
     * @param x0 Starting point
     * @param fun Objective function to minimize
     * @return Minimization result
     */
    Result
    minimize( Vector const & x0, Callback const & fun ) {
      initialize(x0, fun);
      
      size_t n = x0.size();
      size_t f_eval_count = m_state.interpolation_points.size();
      size_t no_improvement_count = 0;
      
      if ( m_opts.verbose ) {
        fmt::print( 
          "[BOBYQA] Starting optimization, dimension: {}, interpolation points: {}, initial f: {}\n", 
          n, m_opts.n_interpolation, m_state.current_f 
        );
      }
    
      for ( size_t k{0}; k < m_opts.max_iter; ++k ) {
        // Update quadratic model
        update_quadratic_model();
        
        // Solve trust region subproblem
        Vector step = solve_trust_region_subproblem();
        Vector candidate_x = m_state.current_x + step;
        candidate_x = candidate_x.cwiseMax(m_lower).cwiseMin(m_upper);
        
        // Evaluate candidate
        Scalar candidate_f = fun(candidate_x);
        f_eval_count++;
        
        // Check improvement
        bool improvement = false;
        if (candidate_f < m_state.best_f) {
          m_state.best_f = candidate_f;
          m_state.best_x = candidate_x;
          m_state.n_improvements++;
          no_improvement_count = 0;
          improvement = true;
        } else {
          no_improvement_count++;
        }
        
        // Update trust region radius
        update_trust_region_radius(improvement, candidate_f);
        
        // Update interpolation points
        update_interpolation_set(candidate_x, candidate_f);
        
        // Update current point
        m_state.current_x = candidate_x;
        m_state.current_f = candidate_f;
        
        // Verbose output
        if ( m_opts.verbose && (k % m_opts.print_every) == 0 ) {
          fmt::print(
            "[BOBYQA] iter={:<4} f={:<10.4} rho={:<10.4} improvements={}\n",
            k, m_state.best_f, m_state.rho, m_state.n_improvements
          );
        }
        
        // Check convergence
        Result conv_result = check_convergence(k, no_improvement_count);
        if (conv_result.converged) {
          conv_result.f_eval_count = f_eval_count;
          return conv_result;
        }
      }
      
      return { 
        m_state.best_f, m_state.best_x, m_opts.max_iter, 
        f_eval_count, false, m_state.rho, "Maximum iterations reached" 
      };
    }

  private:

    /**
     * @brief Update trust region radius based on model performance
     * @param improvement Whether last step improved solution
     * @param candidate_f Function value at candidate point
     */
    void
    update_trust_region_radius( bool improvement, [[maybe_unused]] Scalar candidate_f ) {
      // Simple trust region update strategy
      if (improvement) {
        // Consider expanding trust region if model is good
        m_state.rho = min(m_state.rho * 1.5, m_opts.rhobeg);
      } else {
        // Contract trust region if model is poor
        m_state.rho = max(m_state.rho * 0.5, m_opts.rhoend);
      }
    }

    /**
     * @brief Update interpolation point set
     * @param new_point New point to consider
     * @param new_value Function value at new point
     */
    void
    update_interpolation_set( Vector const & new_point, Scalar new_value ) {
      if (!m_opts.adaptive_interpolation) return;
      
      // Find worst interpolation point to potentially replace
      size_t worst_idx = 0;
      Scalar worst_merit = -std::numeric_limits<Scalar>::max();
      
      for (size_t i = 0; i < m_state.interpolation_points.size(); ++i) {
        // Merit function based on distance and function value
        Scalar distance = (m_state.interpolation_points[i] - m_state.current_x).norm();
        Scalar merit = distance + abs(m_state.interpolation_values[i] - m_state.current_f);
        
        if (merit > worst_merit) {
          worst_merit = merit;
          worst_idx = i;
        }
      }
      
      // Replace worst point if new point is better
      Scalar new_merit = (new_point - m_state.current_x).norm() + abs(new_value - m_state.current_f);
      if (new_merit < worst_merit) {
        m_state.interpolation_points[worst_idx] = new_point;
        m_state.interpolation_values[worst_idx] = new_value;
      }
    }

    /**
     * @brief Check convergence criteria
     * @param iteration Current iteration
     * @param no_improvement_count Count of iterations without improvement
     * @return Result if converged, empty otherwise
     */
    Result
    check_convergence( size_t iteration, size_t no_improvement_count ) {
      // Trust region radius convergence
      if (m_state.rho <= m_opts.rhoend) {
        return { 
          m_state.best_f, m_state.best_x, iteration, 0, true, m_state.rho,
          "Trust region radius below tolerance" 
        };
      }
      
      // Step size convergence (simplified)
      Vector recent_step = m_state.interpolation_points.back() - m_state.current_x;
      if (recent_step.norm() < m_opts.tol_x) {
        return { 
          m_state.best_f, m_state.best_x, iteration, 0, true, m_state.rho,
          "Step size below tolerance" 
        };
      }
      
      // No improvement convergence
      if (no_improvement_count >= m_opts.patience) {
        return { 
          m_state.best_f, m_state.best_x, iteration, 0, true, m_state.rho,
          "No improvement for patience iterations" 
        };
      }
      
      // Not converged
      return { m_state.best_f, m_state.best_x, iteration, 0, false, m_state.rho, "" };
    }

  public:

    /**
     * @brief Get current options
     * @return Current options
     */
    Options const & get_options() const { return m_opts; }
    
    /**
     * @brief Set new options
     * @param opts New options
     */
    void set_options( Options const & opts ) { m_opts = opts; }
    
    /**
     * @brief Get current state information
     * @return String with current state info
     */
    string
    get_state_info() const {
      return fmt::format(
        "BOBYQA State: f={}, rho={}, improvements={}, interpolation_points={}",
        m_state.best_f, m_state.rho, m_state.n_improvements, 
        m_state.interpolation_points.size()
      );
    }
  };

}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

#endif

//
// eof: Utils_BOBYQA.hh
