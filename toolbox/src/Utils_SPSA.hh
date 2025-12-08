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
// file: Utils_SPSA.hh
//
// Implementation of Stochastic Optimization Algorithms using Simultaneous
// Perturbation Stochastic Approximation (SPSA) with convergence improvements.
//
// Main references:
// - Spall, J. C. (1998). Implementation of the simultaneous perturbation
//   algorithm for stochastic optimization. IEEE Transactions on Aerospace and
//   Electronic Systems.
// - Spall, J. C. (2000). Adaptive stochastic approximation by the simultaneous
//   perturbation method. IEEE Transactions on Automatic Control.
// - Spall, J. C. (2003). Introduction to stochastic search and optimization.
//   John Wiley & Sons.
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_SPSA_dot_HH
#define UTILS_SPSA_dot_HH

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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::pow;
  using std::sqrt;

#if 0

  /**
   * @class SPSAGradientEstimator
   * @brief Gradient estimation using Simultaneous Perturbation Stochastic Approximation
   * 
   * This class implements the SPSA algorithm for estimating the gradient of a function
   * using only function evaluations, without needing analytical gradients.
   * The algorithm is particularly useful for optimizing expensive or non-differentiable functions.
   * 
   * The method is based on the formula:
   * \f[
   * g_i \approx \frac{f(x + c \odot \Delta) - f(x - c \odot \Delta)}{2 c_i \Delta_i}
   * \f]
   * where \f$\Delta\f$ is a random perturbation vector and \f$c\f$ is the amplitude vector.
   * 
   * References:
   * - Spall, J. C. (1998). Implementation of the simultaneous perturbation 
   *   algorithm for stochastic optimization. IEEE Transactions on Aerospace and Electronic Systems.
   * - Spall, J. C. (2000). Adaptive stochastic approximation by the simultaneous
   *   perturbation method. IEEE Transactions on Automatic Control.
   * 
   * @tparam Scalar Data type (double, float, etc.)
   */
  template <typename Scalar = double>
  class SPSAGradientEstimator {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>; ///< Vector type
    using Callback = std::function<Scalar(Vector const &)>;    ///< Function callback type

    /**
     * @struct Options
     * @brief Options for SPSA gradient estimator
     */
    struct Options {
      size_t   repeats            = 10;            ///< Number of repetitions for gradient averaging
      Scalar   c_base             = Scalar(1e-3);  ///< Base perturbation amplitude
      Vector   c_per_component;                    ///< Component-specific amplitudes (if empty, use c_base)
      bool     use_rademacher     = true;          ///< True for ±1 noise, false for Gaussian N(0,1)
      bool     apply_precond      = true;         ///< Apply diagonal preconditioning
      Vector   preconditioner;                     ///< Preconditioner weights vector
      unsigned int rng_seed       = std::random_device{}(); ///< Random generator seed
      Scalar   min_delta_abs      = Scalar(1e-12); ///< Minimum value to avoid division by zero
    };

  private:
    Options m_opts; ///< Configured options
    Vector  m_tmp_grad; ///< Temporary buffer for gradients

  public:

    /**
     * @brief Constructor
     * @param opts Options for the estimator
     */
    SPSAGradientEstimator(Options opts = Options())
    : m_opts(std::move(opts))
    {}

    /**
     * @brief Estimate gradient at point x using SPSA
     * 
     * The algorithm estimates the gradient using the formula:
     * \f[
     * g_i \approx \frac{f(x + c \odot \Delta) - f(x - c \odot \Delta)}{2 c_i \Delta_i}
     * \f]
     * 
     * where \f$\Delta\f$ is a random perturbation vector (±1 or Gaussian) and
     * \f$\odot\f$ denotes component-wise product.
     * 
     * The method is particularly efficient for expensive functions since it requires
     * only 2 function evaluations regardless of problem dimensionality.
     * 
     * @param x Point at which to estimate gradient
     * @param f Objective function
     * @param out_grad Estimated gradient (output)
     * @return Number of function evaluations performed
     */
    size_t
    estimate( Vector const & x, Callback const & f, Vector & out_grad ) {
      size_t const n = static_cast<size_t>(x.size());
      assert(n > 0);

      // Prepare c_i
      Vector c;
      if ( m_opts.c_per_component.size() == 0 ) {
        c = Vector::Constant( n, m_opts.c_base );
      } else {
        assert( m_opts.c_per_component.size() == n );
        c = m_opts.c_per_component;
      }

      // Preconditioner
      if ( m_opts.apply_precond ) assert( m_opts.preconditioner.size() == n );

      out_grad.setZero(n);
      if ( m_opts.repeats == 0 ) return 0;

      // RNGs
      std::mt19937 rng(m_opts.rng_seed);
      std::uniform_int_distribution<size_t> rademacher_dist(0, 1); // 0 -> -1, 1 -> +1
      std::normal_distribution<Scalar>      gaussian_dist(0.0, 1.0);

      Vector x_plus(n);
      Vector x_minus(n);
      Vector delta(n);
      Scalar y_plus, y_minus;

      size_t eval_count{0};

      for ( size_t rep{0}; rep < m_opts.repeats; ++rep ) {

        // Sample delta
        if ( m_opts.use_rademacher ) {
          for ( size_t i{0}; i < n; ++i )
            delta[i] = (rademacher_dist(rng) ? Scalar(1) : Scalar(-1));
        } else {
          for ( size_t i{0}; i < n; ++i ) {
            Scalar d = gaussian_dist(rng);
            // Prevent too small delta (rare for Gaussian, but useful)
            if ( abs(d) < m_opts.min_delta_abs )
              d = std::copysign( m_opts.min_delta_abs, d == 0 ? Scalar(1) : d );
            delta[i] = d;
          }
        }

        // Construct x± = x ± c .* delta
        x_plus  = x + c.cwiseProduct(delta);
        x_minus = x - c.cwiseProduct(delta);

        // Evaluate function (note: we pass nullptr for gradient output; callback can ignore second argument)
        y_plus  = f(x_plus);
        y_minus = f(x_minus);
        eval_count += 2;

        // Update estimate: g_i += (y+ - y-) / (2 * c_i * delta_i)
        const Scalar denom_factor = Scalar(0.5); // Multiply by (y+ - y-) * 0.5 / (c_i * delta_i)
        Scalar diff = (y_plus - y_minus);

        for (size_t i{0}; i < n; ++i) {
          Scalar di = delta[i];
          // Safety: if delta very close to zero (only possible with Gaussian), skip or use fallback
          if (std::abs(di) < m_opts.min_delta_abs) {
            // Fallback: ignore contribution of this repetition for component i
            continue;
          }
          out_grad[i] += diff * (denom_factor / (c[i] * di));
        }
      } // end repeats

      // Average
      out_grad /= static_cast<Scalar>(m_opts.repeats);

      // Apply preconditioning (if requested): simply scale component-wise
      if (m_opts.apply_precond) {
        out_grad = out_grad.cwiseProduct(m_opts.preconditioner);
      }

      return eval_count;
    }

    Options const & options() const { return m_opts; }
    Options & options() { return m_opts; }
  };

#endif

  /**
   * @class SPSA_minimizer
   * @brief Minimizer using Simultaneous Perturbation Stochastic Approximation
   *
   * Implements the SPSA algorithm for gradient-free optimization. Particularly
   * effective for high-dimensional problems or noisy objective functions.
   *
   * The main algorithm combines:
   * - Gradient estimation via SPSA
   * - Stochastic descent with adaptive learning rate
   * - Box constraints handling via clipping
   * - Robustness strategies (backtracking, restart, convergence monitoring)
   *
   * Key features:
   * - Efficient for high-dimensional problems (requires only 2 function
   * evaluations per iteration)
   * - Robust to noise in function evaluations
   * - Supports box constraints through projection
   * - Implements adaptive learning rate strategies
   * - Convergence monitoring with multiple criteria
   *
   * References:
   * - Spall, J. C. (1998). Implementation of the simultaneous perturbation
   *   algorithm for stochastic optimization. IEEE Transactions on Aerospace and
   * Electronic Systems.
   * - Spall, J. C. (2000). Adaptive stochastic approximation by the
   * simultaneous perturbation method. IEEE Transactions on Automatic Control.
   * - Spall, J. C. (2003). Introduction to stochastic search and optimization.
   *   John Wiley & Sons.
   *
   * @tparam Scalar Data type (double, float, etc.)
   */
  template <typename Scalar = double>
  class SPSA_minimizer
  {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;  ///< Vector type
    using Callback = std::function<Scalar( Vector const & )>;   ///< Function callback type

    /**
     * @struct Options
     * @brief Options for SPSA minimizer
     */
    struct Options
    {
      size_t max_iter{ 1000 };        ///< Maximum iterations
      size_t gradient_avg{ 2 };       ///< SPSA averages for gradient (reduces variance)
      Scalar a0{ 0.3 };               ///< Base learning rate
      Scalar c0{ 0.05 };              ///< Perturbation amplitude
      Scalar alpha{ 0.602 };          ///< Learning rate decay exponent (standard SPSA)
      Scalar gamma{ 0.101 };          ///< Perturbation decay exponent (standard SPSA)
      bool   use_projection{ true };  ///< Enable projection to bounds
      bool   verbose{ false };        ///< Verbose output
      size_t print_every{ 100 };

      // Convergence criteria
      Scalar tol_grad{ 1e-8 };          ///< Gradient norm tolerance
      size_t tol_grad_patience{ 100 };  ///< Iterations without improvement before stop
      Scalar tol_f{ 1e-10 };            ///< Function change tolerance
      size_t tol_f_patience{ 100 };     ///< Iterations without improvement before stop
      size_t patience{ 100 };           ///< Iterations without improvement before stop

      // Adaptive step size control
      bool   use_backtracking{ true };  ///< Enable backtracking line search
      size_t max_backtrack{ 10 };       ///< Maximum backtracking steps
      Scalar backtrack_factor{ 0.5 };   ///< Step size reduction factor
      Scalar max_step_ratio{ 10.0 };    ///< Maximum step size ratio for clipping
    };

    /**
     * @struct Result
     * @brief Minimization result
     */
    struct Result
    {
      Scalar final_f;       ///< Final function value
      Scalar grad_norm;     ///< Final gradient norm
      size_t iterations;    ///< Number of iterations performed
      Vector final_x;       ///< Final point
      bool   converged;     ///< True if convergence reached
      size_t f_eval_count;  ///< Number of function evaluations
      string message;       ///< Termination message
    };

  private:
    Options      m_opts;                         ///< Configured options
    Vector       m_lower;                        ///< Lower bounds
    Vector       m_upper;                        ///< Upper bounds
    std::mt19937 rng{ std::random_device{}() };  ///< Random generator

  public:
    /**
     * @brief Constructor
     * @param o Options for the minimizer
     */
    SPSA_minimizer( Options const & o = Options() ) : m_opts( o ) {}

    /**
     * @brief Set optimization bounds
     * @param lo Lower bounds
     * @param up Upper bounds
     */
    void
    set_bounds( Vector const & lo, Vector const & up )
    {
      assert( lo.size() == up.size() );
      m_lower = lo;
      m_upper = up;
    }

    void
    set_bounds( size_t n, Scalar const lower[], Scalar const upper[] )
    {
      m_lower.resize( n );
      m_upper.resize( n );
      std::copy_n( lower, n, m_lower.data() );
      std::copy_n( upper, n, m_upper.data() );
    }

  private:
    /**
     * @brief Project point x onto bounds
     * @param x Point to project (modified in-place)
     */
    void
    project( Vector & x ) const
    {
      if ( m_opts.use_projection ) x = x.cwiseMax( m_lower ).cwiseMin( m_upper );
    }

    /**
     * @brief Generate sparse Rademacher noise (-1, +1)
     *
     * @param d Vector to store noise (output)
     */
    void
    rademacher( Vector & d )
    {
      // Use uniform distribution for [0,1]
      std::uniform_int_distribution<> U( 0, 100 );
      // Generate sparse Rademacher noise
      d = Vector::NullaryExpr( d.size(), [&]() { return -1 + 2 * ( U( rng ) % 2 ); } );
    }

    /**
     * @brief Perform backtracking line search
     * @param x Current point
     * @param g Current gradient
     * @param a_k Current learning rate
     * @param f Current function value
     * @param fun Objective function
     * @param x_new New point (output)
     * @param f_new New function value (output)
     * @param eval_count Function evaluation counter (incremented)
     * @return True if backtracking found a better point
     */
    bool
    backtracking_line_search( Vector const &   x,
                              Vector const &   g,
                              Scalar           a_k,
                              Scalar           f,
                              Callback const & fun,
                              Vector &         x_new,
                              Scalar &         f_new,
                              size_t &         eval_count ) const
    {
      if ( !m_opts.use_backtracking ) return false;

      Scalar a_temp = a_k;
      // Try up to max_backtrack reduced step sizes
      for ( size_t backtrack{ 0 }; backtrack < m_opts.max_backtrack; ++backtrack )
      {
        a_temp *= m_opts.backtrack_factor;
        Vector x_temp = x - a_temp * g;
        project( x_temp );
        Scalar f_temp = fun( x_temp );
        eval_count++;

        if ( f_temp < f )
        {
          x_new = x_temp;
          f_new = f_temp;
          return true;
        }
      }
      return false;
    }

  public:
    /**
     * @brief Perform minimization using SPSA with robust convergence handling
     *
     * Convergence criteria:
     * - Gradient norm below tolerance
     * - Function value change below tolerance
     * - No improvement for patience iterations
     * - Maximum iterations reached
     *
     * @param x0 Starting point
     * @param fun Objective function to minimize
     * @return Minimization result
     */
    Result
    minimize( Vector const & x0, Callback const & fun )
    {
      size_t n = x0.size();
      Vector x = x0;
      Vector g( n ), delta( n );
      Vector xp( n ), xm( n );

      // Initialize tracking variables
      project( x );
      Scalar f                    = fun( x );
      Scalar best_f               = f;
      Scalar previous_f           = f;
      Vector best_x               = x;
      size_t no_improvement_count = 0;
      size_t f_eval_count         = 1;

      // Track step size adaptations
      size_t total_backtrack_success = 0;
      size_t total_consecutive_grad  = 0;
      size_t total_consecutive_func  = 0;

      if ( m_opts.verbose ) fmt::print( "[SPSA] Starting optimization, dimension: {}, initial f: {}\n", n, f );

      for ( size_t k{ 0 }; k < m_opts.max_iter; ++k )
      {
        // Calculate learning rate and perturbation with standard SPSA decay
        Scalar a_k = m_opts.a0 / pow( Scalar( k + 1 ), m_opts.alpha );
        Scalar c_k = m_opts.c0 / pow( Scalar( k + 1 ), m_opts.gamma );

        g.setZero();

        // ============================
        //   SPSA GRADIENT ESTIMATION
        // ============================
        for ( size_t rep{ 0 }; rep < m_opts.gradient_avg; ++rep )
        {
          // Generate random perturbations
          rademacher( delta );

          xp = x + c_k * delta;
          xm = x - c_k * delta;

          project( xp );
          project( xm );

          // Evaluate function at perturbed points (guaranteed within bounds)
          Scalar fp = fun( xp );
          Scalar fm = fun( xm );
          f_eval_count += 2;

          // GRADIENT CALCULATION
          g.array() += ( fp - fm ) / ( xp - xm ).array();
        }

        // Average gradient estimates
        g.array() /= Scalar( m_opts.gradient_avg );
        Scalar gnorm = g.template lpNorm<Eigen::Infinity>();

        // Verbose output (every 50 iterations to avoid clutter)
        if ( m_opts.verbose && ( k % m_opts.print_every ) == 0 )
        {
          fmt::print(
              "[SPSA] iter={:<4} f={:<10.4} |g|={:<10.4} a_k={:<8.2} "
              "c_k={:<8.2} improv_count={}\n",
              k, f, gnorm, a_k, c_k, no_improvement_count );
        }

        // CONVERGENCE CRITERIA
        bool   converged           = false;
        string convergence_message = "Maximum iterations reached";

        // 1. Gradient norm criterion
        if ( !converged )
        {
          if ( gnorm < m_opts.tol_grad )
          {
            if ( ++total_consecutive_grad > m_opts.tol_grad_patience )
            {
              converged           = true;
              convergence_message = "Gradient norm below tolerance";
            }
          }
          else
          {
            total_consecutive_grad = 0;
          }
        }

        // 2. Function change criterion
        if ( !converged )
        {
          if ( abs( f - previous_f ) < m_opts.tol_f * ( 1 + abs( previous_f ) ) )
          {
            if ( ++total_consecutive_func > m_opts.tol_f_patience )
            {
              converged           = true;
              convergence_message = "Function change below tolerance";
            }
          }
          else
          {
            total_consecutive_func = 0;
          }
        }

        // 3. Stagnation criterion
        if ( !converged )
        {
          if ( no_improvement_count > m_opts.patience )
          {
            converged           = true;
            convergence_message = "No improvement for patience iterations";
          }
        }

        if ( converged )
        {
          if ( m_opts.verbose )
          {
            fmt::print( "[SPSA] Converged at iteration {}: {}\n",
                        "Final f: {}, gradient norm: {}, function evaluations: {}\n", k, convergence_message, best_f,
                        gnorm, f_eval_count );
          }
          return { best_f, gnorm, k, best_x, true, f_eval_count, convergence_message };
        }

        // ============================
        //       PARAMETER UPDATE
        // ============================
        Vector x_new = x - a_k * g;
        project( x_new );

        Scalar f_new = fun( x_new );
        f_eval_count++;

        // TRACK BEST SOLUTION
        if ( f_new < best_f )
        {
          best_f               = f_new;
          best_x               = x_new;
          no_improvement_count = 0;  // Reset counter
        }
        else
        {
          no_improvement_count++;
        }

        // BACKTRACKING LINE SEARCH if step causes deterioration
        bool backtrack_used = false;
        if ( f_new > f && m_opts.use_backtracking )
        {
          backtrack_used = backtracking_line_search( x, g, a_k, f, fun, x_new, f_new, f_eval_count );
          if ( backtrack_used )
          {
            total_backtrack_success++;
            // Update best solution if found during backtracking
            if ( f_new < best_f )
            {
              best_f               = f_new;
              best_x               = x_new;
              no_improvement_count = 0;
            }
          }
        }

        // RESTART STRATEGY if stuck for many iterations
        if ( no_improvement_count > m_opts.patience / 2 )
        {
          if ( m_opts.verbose )
          {
            fmt::print(
                "[SPSA] Restarting from best point at iteration {} (no "
                "improvement for {} iterations)\n",
                k, no_improvement_count );
          }
          x                    = best_x;
          f                    = best_f;
          no_improvement_count = 0;
          continue;  // Skip rest of iteration
        }

        // ACCEPT NEW POINT
        previous_f = f;
        x          = x_new;
        f          = f_new;
      }

      // Return best solution found (not necessarily the last one)
      if ( m_opts.verbose )
      {
        fmt::print(
            "[SPSA] Optimization finished:\n"
            "best f:               {}\n"
            "function evaluations: {}\n"
            "backtrack successes:  {}\n",
            best_f, f_eval_count, total_backtrack_success );
      }

      return { best_f, 0, m_opts.max_iter, best_x, false, f_eval_count, "Maximum iterations reached" };
    }

    /**
     * @brief Get current options
     * @return Current options
     */
    Options const &
    get_options() const
    {
      return m_opts;
    }

    /**
     * @brief Set new options
     * @param opts New options
     */
    void
    set_options( Options const & opts )
    {
      m_opts = opts;
    }
  };
}  // namespace Utils

#endif

#endif

//
// eof: Utils_SPSA.hh
//
