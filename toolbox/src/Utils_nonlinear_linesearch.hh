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
// file: Utils_Linesearch.hh
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_LINESEARCH_dot_HH
#define UTILS_LINESEARCH_dot_HH

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
#include <cmath>
#include <limits>
#include <optional>
#include <random>
#include <utility>

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;

  // ===========================================================================
  // Line Search Helper Functions
  // ===========================================================================

  namespace detail
  {
    /**
     * @brief Robust cubic interpolation for line-search (Hager–Zhang friendly).
     *
     * This version provides:
     *  - Full finite checks (NaN/inf)
     *  - Degeneracy detection (a ≈ b, fa ≈ fb, etc.)
     *  - Hager–Zhang safeguard window: [a_lo+δΔ, a_hi−σΔ]
     *  - Fallback to bisection in all doubtful cases
     */
    template <class Scalar>
    inline Scalar
    cubic_minimizer( Scalar a,
                     Scalar fa,
                     Scalar fpa,
                     Scalar b,
                     Scalar fb,
                     Scalar fpb,
                     Scalar delta = Scalar( 0.1 ),  // Hager–Zhang recommended
                     Scalar sigma = Scalar( 0.9 )   // Hager–Zhang recommended
    )
    {
      Scalar eps{ std::numeric_limits<Scalar>::epsilon() };

      // reorder interval
      Scalar lo    = min( a, b );
      Scalar hi    = max( a, b );
      Scalar width = hi - lo;

      // fallback if degenerate interval
      if ( width <= 10 * eps * ( 1 + abs( lo ) ) ) return ( a + b ) / 2;

      // -----------------------------------------------------------------------
      // Standard More–Thuente cubic setup
      // -----------------------------------------------------------------------
      Scalar denom_ab = ( a - b );
      if ( abs( denom_ab ) < eps ) return ( a + b ) / 2;

      Scalar d1    = fpa + fpb - 3 * ( ( fa - fb ) / denom_ab );
      Scalar discr = d1 * d1 - fpa * fpb;

      // discriminant negative → no real minimizer → bisection
      if ( discr < 0 ) return ( a + b ) / 2;

      Scalar sqrt_disc = sqrt( max( discr, Scalar( 0 ) ) );

      Scalar denom = fpb - fpa + 2 * sqrt_disc;
      if ( !std::isfinite( denom ) || abs( denom ) < eps ) return ( a + b ) / 2;

      // cubic minimizer
      Scalar t = b - ( b - a ) * ( ( fpb + sqrt_disc - d1 ) / denom );

      // -----------------------------------------------------------------------
      // Hager–Zhang safeguard: force t inside:
      //   [ a_lo + δΔ , a_hi − σΔ ]
      // -----------------------------------------------------------------------
      Scalar left  = lo + delta * width;
      Scalar right = lo + sigma * width;

      if ( !std::isfinite( t ) || t <= left || t >= right ) return ( a + b ) / 2;
      return t;
    }

    /**
     * @brief Compute next trial step for More-Thuente line search
     *
     * This function implements the step selection logic from More-Thuente,
     * choosing between cubic interpolation, quadratic interpolation, and
     * bisection based on the current state of the line search.
     *
     * ### Strategy
     *
     * 1. **Primary method**: Cubic interpolation using interval endpoints
     * 2. **Safeguarding**: Switch to bisection or extrapolation if cubic fails
     * 3. **Phase-dependent logic**:
     *    - Bracketing phase: Use extrapolation if cubic invalid
     *    - Zoom phase: Use bisection if cubic invalid
     *
     * @param a_lo Lower bound of current interval with function/derivative info
     * @param phi_lo Function value at a_lo
     * @param der_lo Derivative at a_lo
     * @param a_hi Upper bound of current interval
     * @param phi_hi Function value at a_hi
     * @param der_hi Derivative at a_hi
     * @param a_prev Previous trial point
     * @param phi_prev Function value at previous point
     * @param der_prev Derivative at previous point
     * @param alpha_max Maximum allowable step length
     * @param step_max Maximum step to take in this iteration
     * @param is_bracketing true if in bracketing phase, false if in zoom phase
     *
     * @return Next trial step length
     *
     * @note This is a simplified but robust version of MT step computation
     */
    template <typename Scalar>
    Scalar
    compute_step( Scalar                  a_lo,
                  Scalar                  phi_lo,
                  Scalar                  der_lo,
                  Scalar                  a_hi,
                  Scalar                  phi_hi,
                  Scalar                  der_hi,
                  [[maybe_unused]] Scalar a_prev,
                  [[maybe_unused]] Scalar phi_prev,
                  [[maybe_unused]] Scalar der_prev,
                  Scalar                  alpha_max,
                  Scalar                  step_max,
                  bool                    is_bracketing )
    {
      Scalar a_j = cubic_minimizer( a_lo, phi_lo, der_lo, a_hi, phi_hi, der_hi );
      if ( a_j <= a_lo || a_j >= a_hi )
      {
        a_j = is_bracketing ? ( a_hi + 2 * ( a_hi - a_prev ) ) : ( ( a_lo + a_hi ) / 2 );
        if ( a_j > alpha_max ) a_j = alpha_max;
      }
      Scalar tol = 1e-6 * ( a_hi - a_lo );
      a_j        = std::clamp( a_j, a_lo + tol, a_hi - tol );
      return min( a_j, step_max );
    }

  }  // namespace detail


  // ===========================================================================
  // Line Search Algorithms
  // ===========================================================================

  /**
   * @class ArmijoLineSearch
   * @brief Simple Armijo (sufficient decrease) line search
   *
   * Implements the classical Armijo backtracking line search, which finds a
   * step length satisfying the Armijo condition:
   * \f[
   *   f(x + \alpha p) \leq f(x) + c_1 \alpha \nabla f(x)^T p
   * \f]
   *
   * ### Algorithm
   *
   * 1. Start with initial step α₀ (typically 1.0)
   * 2. Evaluate f(x + α p)
   * 3. If Armijo condition satisfied, accept α
   * 4. Otherwise, reduce α ← ρ α (ρ ≈ 0.5)
   * 5. Repeat until convergence or max iterations
   *
   * ### Characteristics
   *
   * **Advantages:**
   * - Very simple and fast per iteration
   * - Robust: always terminates if gradient is descent direction
   * - No derivative evaluations at trial points needed
   * - Works well for Newton-like methods with good curvature
   *
   * **Disadvantages:**
   * - May accept steps that are too short (inefficient progress)
   * - No guarantee of adequate step length for superlinear convergence
   * - Can lead to many small steps in poorly scaled problems
   *
   * ### When to Use
   *
   * Armijo is recommended for:
   * - Newton or quasi-Newton methods with good Hessian approximation
   * - Problems where function evaluations are cheap
   * - When simplicity and reliability are more important than efficiency
   * - Initial iterations of optimization (before switching to stronger
   * conditions)
   *
   * @tparam Scalar Floating-point type (double or float)
   *
   * @note This is the simplest practical line search. For better efficiency,
   *       consider Wolfe or Strong Wolfe conditions.
   *
   * @see Armijo, L. (1966). "Minimization of functions having Lipschitz
   *      continuous first partial derivatives". Pacific Journal of Mathematics.
   * @see Nocedal & Wright (2006), Algorithm 3.1
   */
  template <typename Scalar>
  class ArmijoLineSearch
  {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1          = 1e-4;   ///< Armijo constant (typical: 1e-4)
    Scalar m_step_reduce = 0.5;    ///< Step reduction factor when Armijo fails
    Scalar m_step_expand = 1.2;    ///< Step expansion factor (currently unused)
    size_t m_max_iters   = 50;     ///< Maximum backtracking iterations
    Scalar m_epsi        = 1e-15;  ///< Minimum step size threshold

    mutable Vector m_x_new;  ///< Workspace for trial point

  public:
    /**
     * @brief Perform Armijo backtracking line search
     *
     * @tparam Callback Function type with signature: Scalar(Vector const&,
     * Vector*)
     *
     * @param f0 Function value at current point x
     * @param Df0 Directional derivative ∇f(x)^T d (must be negative)
     * @param x Current point
     * @param d Search direction
     * @param callback Function that evaluates f (and optionally ∇f)
     * @param step0 Initial trial step length (default: 1.0)
     *
     * @return Step length if successful, std::nullopt if failed
     *
     * @note If line search fails repeatedly, returns best step found even if
     *       it doesn't strictly satisfy Armijo
     */
    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar f0, Scalar Df0, Vector const & x, Vector const & d, Callback const & callback, Scalar step0 = 1 )
        const
    {
      Scalar     step   = step0;
      Scalar     c1_Df0 = m_c1 * Df0;
      auto const n{ x.size() };
      m_x_new.resize( n );

      Scalar best_step    = step0;
      Scalar best_f       = std::numeric_limits<Scalar>::max();
      size_t shrink_count = 0;
      size_t iteration    = 0;

      for ( ; iteration < m_max_iters; ++iteration )
      {
        m_x_new.noalias() = x + step * d;
        Scalar f_new      = callback( m_x_new, nullptr );

        // Track best step found
        if ( f_new < best_f )
        {
          best_f    = f_new;
          best_step = step;
        }

        if ( f_new <= f0 + step * c1_Df0 ) return std::tuple<Scalar, size_t>( step, iteration );  // Success

        // Adaptive reduction: more aggressive after multiple failures
        Scalar reduction = m_step_reduce;
        if ( shrink_count > 2 ) reduction = 0.1;

        step *= reduction;
        shrink_count++;

        // Check for excessively small steps
        if ( step < m_epsi )
        {
          // Return best step found even if it doesn't satisfy Armijo
          if ( best_f < f0 ) return std::tuple<Scalar, size_t>( best_step, iteration );
          break;
        }
      }

      // Fallback: return best step if it improved objective
      if ( best_f < f0 ) return std::tuple<Scalar, size_t>( best_step, iteration );

      return std::nullopt;
    }
  };

  /**
   * @class WolfeLineSearch
   * @brief Wolfe conditions line search (sufficient decrease + curvature)
   *
   * Implements a line search satisfying the Wolfe conditions:
   * \f{align*}{
   *   f(x + \alpha p) &\leq f(x) + c_1 \alpha \nabla f(x)^T p \quad
   * &\text{(Armijo)}\\
   *   \nabla f(x + \alpha p)^T p &\geq c_2 \nabla f(x)^T p \quad
   * &\text{(Curvature)}
   * \f}
   * with \f$0 < c_1 < c_2 < 1\f$.
   *
   * ### Algorithm Structure
   *
   * **Phase 1 - Bracketing**: Find interval [α_lo, α_hi] containing acceptable
   * step
   * - Start with α = α₀
   * - Increase α until Armijo violated or curvature satisfied
   * - If Armijo violated, bracket found
   * - If derivative becomes positive, bracket found
   *
   * **Phase 2 - Zoom**: Refine within bracket using cubic interpolation
   * - Repeatedly subdivide interval
   * - Choose new point via cubic interpolation
   * - Update bracket endpoints
   * - Terminate when curvature condition satisfied
   *
   * ### Characteristics
   *
   * **Advantages:**
   * - Guarantees superlinear convergence for quasi-Newton methods
   * - Prevents steps that are too short (via curvature condition)
   * - Well-suited for L-BFGS and conjugate gradient methods
   * - Theoretically sound with convergence guarantees
   *
   * **Disadvantages:**
   * - Requires gradient evaluations at trial points
   * - More complex than Armijo
   * - May require several iterations to converge
   *
   * ### Parameters
   *
   * Typical values:
   * - c₁ = 1e-4 (Armijo parameter)
   * - c₂ = 0.9 for Newton/quasi-Newton, 0.1 for nonlinear CG
   *
   * @tparam Scalar Floating-point type (double or float)
   *
   * @note For L-BFGS, use c₂ = 0.9. For conjugate gradient, use c₂ = 0.1-0.4.
   *
   * @see Nocedal & Wright (2006), Algorithm 3.5 and 3.6
   * @see Wolfe, P. (1969). "Convergence Conditions for Ascent Methods"
   */

  template <typename Scalar>
  class WeakWolfeLineSearch
  {
    Scalar m_c1        = Scalar( 1e-4 );
    Scalar m_c2        = Scalar( 0.1 );  // Weak curvature - allows longer steps
    Scalar m_alpha_max = Scalar( 50.0 );
    Scalar m_alpha_min = Scalar( 1e-12 );
    Scalar w_shrink    = Scalar( 0.1 );
    size_t m_max_iters = 40;
    size_t m_max_evals = 150;

    template <typename EvalFunc>
    Scalar
    zoom( Scalar           a_lo,
          Scalar           f_lo,
          Scalar           Df_lo,
          Scalar           a_hi,
          Scalar           f_hi,
          Scalar           Df_hi,
          Scalar           f0,
          Scalar           c1_Df0,
          Scalar           c2_Df0,
          EvalFunc const & eval ) const
    {
      Scalar a_j, f_j, Df_j;

      for ( size_t iter = 0; iter < m_max_iters; ++iter )
      {
        // Prefer quadratic interpolation for weak Wolfe
        if ( m_use_quadratic_interp && Df_lo != 0 )
        {
          Scalar a_hilo = a_hi - a_lo;
          a_j           = a_lo - a_hilo * Df_lo / ( Df_lo - ( f_hi - f_lo ) / a_hilo );

          // Safeguard
          Scalar low_bound  = a_lo + w_shrink * a_hilo;
          Scalar high_bound = a_hi - w_shrink * a_hilo;

          if ( a_j <= low_bound || a_j >= high_bound ) a_j = ( a_lo + a_hi ) / 2;
        }
        else
        {
          a_j = ( a_lo + a_hi ) / 2;
        }

        eval( a_j, f_j, Df_j );

        // Update bracket
        if ( f_j > f0 + a_j * c1_Df0 || f_j >= f_lo )
        {
          a_hi  = a_j;
          f_hi  = f_j;
          Df_hi = Df_j;
        }
        else
        {
          // Armijo satisfied, check weak curvature
          if ( Df_j >= c2_Df0 ) return a_j;

          // Update lower bound
          if ( Df_j * ( a_hi - a_lo ) >= 0 ) a_hi = a_lo;
          a_lo  = a_j;
          f_lo  = f_j;
          Df_lo = Df_j;
        }

        if ( abs( a_hi - a_lo ) < m_alpha_min ) return a_j;
      }
      return 0;
    }

  public:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    WeakWolfeLineSearch() {}

    // Adaptive step control
    Scalar m_extrapolation_factor = Scalar( 2.2 );
    bool   m_use_quadratic_interp = true;

    mutable Vector m_g_new, m_x_new;

    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar           f0,
                Scalar           Df0,
                Vector const &   x,
                Vector const &   d,
                Callback const & callback,
                Scalar           alpha0 = Scalar( 1 ) ) const
    {
      if ( Df0 >= 0 ) return std::nullopt;

      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );

      size_t n_evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df )
      {
        ++n_evals;
        m_x_new.noalias() = x + a * d;
        f                 = callback( m_x_new, &m_g_new );
        df                = m_g_new.dot( d );
      };

      Scalar c1_Df0 = m_c1 * Df0;
      Scalar c2_Df0 = m_c2 * Df0;

      Scalar alpha      = alpha0;
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      Scalar phi, der;
      eval( alpha, phi, der );

      // Check weak Wolfe conditions
      if ( phi <= f0 + alpha * c1_Df0 && der >= c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );

      // Simple bracketing with more aggressive extrapolation

      for ( size_t iter{ 0 }; iter < m_max_iters && n_evals < m_max_evals; ++iter )
      {
        // Armijo condition failed or function not decreasing
        if ( phi > f0 + alpha * c1_Df0 || ( iter > 0 && phi >= phi_prev ) )
        {
          Scalar new_alpha = zoom( alpha_prev, phi_prev, der_prev, alpha, phi, der, f0, c1_Df0, c2_Df0, eval );
          if ( new_alpha > 0 ) return std::tuple<Scalar, size_t>( new_alpha, n_evals );
          return std::nullopt;
        }

        // Weak curvature satisfied
        if ( der >= c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );

        // Positive gradient - we've passed a minimum
        if ( der >= 0 )
        {
          Scalar new_alpha = zoom( alpha_prev, phi_prev, der_prev, alpha, phi, der, f0, c1_Df0, c2_Df0, eval );
          if ( new_alpha > 0 ) return std::tuple<Scalar, size_t>( new_alpha, n_evals );
          return std::nullopt;
        }

        // Extrapolate more aggressively
        alpha_prev = alpha;
        phi_prev   = phi;
        der_prev   = der;

        Scalar alpha_new = alpha * m_extrapolation_factor;
        if ( alpha_new > m_alpha_max ) alpha_new = m_alpha_max;

        // Prevent tiny increments
        if ( alpha_new - alpha < m_alpha_min )
        {
          Scalar new_alpha = zoom( alpha_prev, phi_prev, der_prev, alpha, phi, der, f0, c1_Df0, c2_Df0, eval );
          if ( new_alpha > 0 ) return std::tuple<Scalar, size_t>( new_alpha, n_evals );
          return std::nullopt;
        }

        alpha = alpha_new;
        eval( alpha, phi, der );

        // Check if we got lucky with extrapolation
        if ( phi <= f0 + alpha * c1_Df0 && der >= c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );
        ;
      }
      return std::nullopt;
    }
  };

  /**
   * @class StrongWolfeLineSearch
   * @brief Strong Wolfe conditions line search
   *
   * Implements line search satisfying the **strong Wolfe conditions**:
   * \f{align*}{
   *   f(x + \alpha p) &\leq f(x) + c_1 \alpha \nabla f(x)^T p \quad
   * &\text{(Armijo)}\\
   *   |\nabla f(x + \alpha p)^T p| &\leq c_2 |\nabla f(x)^T p| \quad
   * &\text{(Strong Curvature)}
   * \f}
   *
   * The key difference from regular Wolfe is the absolute value in the
   * curvature condition, which prevents the algorithm from accepting steps
   * where the gradient is still strongly negative (indicating the function is
   * still decreasing rapidly).
   *
   * ### Why Strong Wolfe?
   *
   * The strong curvature condition ensures:
   * 1. Step is not too short (function still decreasing rapidly)
   * 2. Step is not too long (function starting to increase)
   * 3. Better guarantee of superlinear convergence
   * 4. More stable behavior near local minima
   *
   * ### Comparison with Regular Wolfe
   *
   * |  | Regular Wolfe | Strong Wolfe |
   * |--|--------------|--------------|
   * | Curvature | φ'(α) ≥ c₂ φ'(0) | \|φ'(α)\| ≤ c₂ \|φ'(0)\| |
   * | Accepts | α with φ'(α) still negative | Only near critical points |
   * | Convergence | Superlinear | Superlinear (better guaranteed) |
   * | Difficulty | Easier to satisfy | Slightly harder |
   *
   * ### When to Use
   *
   * Strong Wolfe is preferred for:
   * - L-BFGS optimization (most common choice)
   * - When high-quality steps are needed for fast convergence
   * - Problems with well-scaled variables
   * - When gradient evaluations are not too expensive
   *
   * @tparam Scalar Floating-point type
   *
   * @see Nocedal & Wright (2006), Section 3.1
   */
  template <typename Scalar>
  class StrongWolfeLineSearch
  {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1        = 1e-4;   ///< Armijo parameter
    Scalar m_c2        = 0.9;    ///< Strong curvature parameter
    Scalar m_alpha_max = 10;     ///< Maximum step length
    Scalar m_epsi      = 1e-12;  ///< Interval convergence tolerance
    size_t m_max_iters = 50;     ///< Maximum iterations

    mutable Vector m_g_new, m_x_new;

  public:
    /**
     * @brief Perform Strong Wolfe line search
     *
     * Similar to regular Wolfe but with stricter curvature condition.
     *
     * @param f0 Function value at current point
     * @param Df0 Directional derivative (must be < 0)
     * @param x Current point
     * @param d Search direction
     * @param callback Function/gradient evaluator
     * @param alpha0 Initial step (default: 1.0)
     *
     * @return Step satisfying strong Wolfe, or std::nullopt if failed
     */
    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar           f0,
                Scalar           Df0,
                Vector const &   x,
                Vector const &   d,
                Callback const & callback,
                Scalar           alpha0 = 1 ) const
    {
      if ( !( Df0 < 0 ) ) return std::nullopt;
      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t n_evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df )
      {
        ++n_evals;
        m_x_new.noalias() = x + a * d;
        f                 = callback( m_x_new, &m_g_new );
        df                = m_g_new.dot( d );
      };

      Scalar c1_Df0{ m_c1 * Df0 };
      Scalar c2_Df0{ m_c2 * Df0 };

      // Initial trial
      Scalar alpha = alpha0;
      Scalar phi, der;
      eval( alpha, phi, der );

      // Check if acceptable (note: abs(der) for strong Wolfe)
      if ( phi <= f0 + alpha * c1_Df0 && abs( der ) <= -c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );

      Scalar alpha_lo = 0;
      Scalar phi_lo   = f0;
      Scalar der_lo   = Df0;

      Scalar alpha_hi = 0;
      Scalar phi_hi   = 0;
      Scalar der_hi   = 0;

      bool   bracketed  = false;
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      // === Bracketing ===
      while ( n_evals < m_max_iters )
      {
        if ( ( phi > f0 + alpha * c1_Df0 ) || ( n_evals > 1 && phi >= phi_prev ) )
        {
          alpha_lo = alpha_prev;
          phi_lo   = phi_prev;
          der_lo   = der_prev;

          alpha_hi = alpha;
          phi_hi   = phi;
          der_hi   = der;

          bracketed = true;
          break;
        }

        if ( abs( der ) <= -c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );

        if ( der >= 0 )
        {
          alpha_lo = alpha_prev;
          phi_lo   = phi_prev;
          der_lo   = der_prev;

          alpha_hi = alpha;
          phi_hi   = phi;
          der_hi   = der;

          bracketed = true;
          break;
        }

        Scalar new_alpha = min( 2 * alpha, m_alpha_max );
        alpha_prev       = alpha;
        phi_prev         = phi;
        der_prev         = der;
        alpha            = new_alpha;
        eval( alpha, phi, der );

        if ( phi <= f0 + alpha * c1_Df0 && abs( der ) <= -c2_Df0 ) return std::tuple<Scalar, size_t>( alpha, n_evals );
      }

      if ( !bracketed ) return std::nullopt;

      // === Zoom ===
      while ( n_evals < m_max_iters )
      {
        Scalar a_j = abs( alpha_hi - alpha_lo ) > m_epsi
                         ? detail::cubic_minimizer<Scalar>( alpha_lo, phi_lo, der_lo, alpha_hi, phi_hi, der_hi )
                         : ( alpha_lo + alpha_hi ) / 2;

        Scalar phi_j, der_j;
        eval( a_j, phi_j, der_j );

        if ( ( phi_j > f0 + a_j * c1_Df0 ) || ( phi_j >= phi_lo ) )
        {
          alpha_hi = a_j;
          phi_hi   = phi_j;
          der_hi   = der_j;
        }
        else
        {
          if ( abs( der_j ) <= -c2_Df0 ) return std::tuple<Scalar, size_t>( a_j, n_evals );

          if ( der_j * ( alpha_hi - alpha_lo ) >= 0 )
          {
            alpha_hi = alpha_lo;
            phi_hi   = phi_lo;
            der_hi   = der_lo;
          }

          alpha_lo = a_j;
          phi_lo   = phi_j;
          der_lo   = der_j;
        }

        if ( abs( alpha_hi - alpha_lo ) < m_epsi ) return std::tuple<Scalar, size_t>( a_j, n_evals );
      }

      return std::nullopt;
    }
  };

  /**
   * @class GoldsteinLineSearch
   * @brief Goldstein conditions line search
   *
   * Implements Goldstein conditions, which require the step length to satisfy:
   * \f[
   *   f(x) + (1-c) \alpha \nabla f(x)^T p \leq f(x+\alpha p) \leq f(x) + c
   * \alpha \nabla f(x)^T p
   * \f]
   * where \f$0 < c < 0.5\f$ (typically c = 0.1 to 0.3).
   *
   * ### Interpretation
   *
   * The Goldstein conditions bound the function decrease from both sides:
   * - **Upper bound** (Armijo-like): Ensures sufficient decrease
   * - **Lower bound**: Prevents steps that are too short
   *
   * This creates a "band" of acceptable function values at each step length.
   *
   * ### Comparison with Wolfe
   *
   * | Feature | Goldstein | Wolfe |
   * |---------|-----------|-------|
   * | Conditions | Two-sided function bounds | Function + derivative |
   * | Derivatives | Only at current point | At current and trial |
   * | Cost | Cheaper (no trial gradients) | More expensive |
   * | Theory | Less common in modern software | Standard choice |
   * | Convergence | Guarantees for convex functions | Stronger guarantees |
   *
   * ### Algorithm Strategy
   *
   * 1. Start with α = α₀
   * 2. Check if both Goldstein conditions are satisfied
   * 3. If f too large: reduce α (upper bound violated)
   * 4. If f too small: increase α (lower bound violated)
   * 5. Repeat until both conditions satisfied
   *
   * ### When to Use
   *
   * Goldstein is useful when:
   * - Gradient evaluations are expensive
   * - Problem is moderately well-conditioned
   * - Simple implementation is preferred
   * - Historical compatibility is needed
   *
   * However, Wolfe conditions are generally preferred in modern implementations
   * due to better theoretical properties and broader applicability.
   *
   * @tparam Scalar Floating-point type
   *
   * @note This is mainly included for completeness. For production use,
   *       prefer Strong Wolfe or More-Thuente line searches.
   *
   * @see Goldstein, A.A. (1965). "On Steepest Descent". SIAM Journal on
   * Control.
   * @see Nocedal & Wright (2006), Section 3.1 for comparison with Wolfe
   */
  template <typename Scalar>
  class GoldsteinLineSearch
  {
    using Vector         = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    Scalar m_c1          = 0.1;  // Parametro per Armijo (tipicamente 0.1-0.3)
    Scalar m_step_reduce = 0.5;  // Fattore di riduzione del passo
    Scalar m_step_expand = 1.2;  // Fattore di espansione del passo (più conservativo)
    size_t m_max_iters   = 50;

    mutable Vector m_x_new;

  public:
    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar           f0,
                Scalar           Df0,
                Vector const &   x,
                Vector const &   d,
                Callback const & callback,
                Scalar           alpha0 = 1 ) const
    {
      // Controllo direzione di discesa
      if ( Df0 >= 0 ) return std::nullopt;

      Scalar alpha = alpha0;
      m_x_new.resize( x.size() );

      // Calcola i bound di Goldstein
      Scalar armijo_bound    = f0 + m_c1 * alpha * Df0;
      Scalar goldstein_bound = f0 + ( 1 - m_c1 ) * alpha * Df0;

      size_t iteration{ 0 };
      for ( ; iteration < m_max_iters; ++iteration )
      {
        m_x_new.noalias() = x + alpha * d;
        Scalar f_new      = callback( m_x_new, nullptr );  // Solo valore funzione, non gradiente

        // Verifica condizioni di Goldstein
        bool satisfies_armijo    = ( f_new <= armijo_bound );
        bool satisfies_goldstein = ( f_new >= goldstein_bound );

        if ( satisfies_armijo && satisfies_goldstein ) { return std::tuple<Scalar, size_t>( alpha, iteration ); }

        if ( !satisfies_armijo )
        {
          // Passo troppo lungo - riduci
          alpha *= m_step_reduce;
        }
        else
        {
          // Passo troppo corto (soddisfa Armijo ma non Goldstein) - aumenta
          alpha *= m_step_expand;
        }

        // Ricalcola i bound con il nuovo alpha
        armijo_bound    = f0 + m_c1 * alpha * Df0;
        goldstein_bound = f0 + ( 1 - m_c1 ) * alpha * Df0;

        // Controllo per alpha troppo piccolo
        if ( alpha < std::numeric_limits<Scalar>::epsilon() * 10 )
        {
          // fmt::print("Goldstein FAILED: step too small\n");
          return std::nullopt;
        }
      }

      // fmt::print("Goldstein FAILED: max iterations reached\n");
      return std::nullopt;
    }
  };

  /**
   * @brief Hager–Zhang inspired line search with enhancements.
   *
   * This class implements a robust line-search procedure inspired by the
   * Hager–Zhang algorithm. It enforces the Armijo (sufficient decrease)
   * condition and a (modified) Wolfe curvature condition with small
   * tolerances to avoid numerical issues. Zoom uses cubic interpolation
   * with explicit δ/σ safeguards; a bisection fallback is used whenever
   * interpolation yields an out-of-range value or NaN/inf.
   *
   * Notes:
   *  - The callback must have signature: Scalar callback(const Vector &x,
   * Vector *g) returning f(x) and writing gradient into *g (if g != nullptr).
   *
   * @tparam Scalar float/double/etc.
   */
  template <typename Scalar>
  class HagerZhangLineSearch
  {
  public:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  private:
    /** @name Tunable parameters */
    Scalar m_c1        = Scalar( 1e-4 );   ///< Armijo parameter.
    Scalar m_c2        = Scalar( 0.9 );    ///< Curvature parameter.
    Scalar m_alpha_max = Scalar( 10.0 );   ///< Maximum allowed step length.
    Scalar m_epsi      = Scalar( 1e-12 );  ///< Absolute tolerance for interval
                                           ///< collapse (used relatively).
    size_t m_max_iters = 50;               ///< Max outer iterations (extrapolation/zoom cycles).
    size_t m_max_evals = 200;              ///< Max function/gradient evaluations total.

    // Hager–Zhang specific interpolation safeguards (now used)
    Scalar m_delta     = Scalar( 0.1 );   ///< left-side safety fraction for interpolation
    Scalar m_sigma     = Scalar( 0.9 );   ///< right-side safety fraction for interpolation
    Scalar m_epsilon_k = Scalar( 1e-6 );  ///< tolerance used to relax curvature checks slightly

    // internal mutable state (so operator() can be const)
    mutable Scalar m_f0;
    mutable Scalar m_Df0;
    mutable Scalar m_c1_Df0;
    mutable Scalar m_c2_Df0;

    mutable Vector m_g_new;
    mutable Vector m_x_new;

    /**
     * @brief Zoom phase: refine [a_lo, a_hi] to find acceptable alpha.
     *
     * Uses cubic minimization (detail::cubic_minimizer) then enforces safety:
     * a_j ∈ [a_lo + δ*(a_hi-a_lo), a_hi - σ*(a_hi-a_lo)].
     * If cubic minimizer returns out of range / NaN / inf, fallback to
     * bisection.
     *
     * @tparam EvalFunc callback type used to compute phi and derivative at a
     * trial alpha.
     */
    template <typename EvalFunc>
    Scalar
    zoom( Scalar a_lo, Scalar phi_lo, Scalar der_lo, Scalar a_hi, Scalar phi_hi, Scalar der_hi, EvalFunc const & eval )
        const
    {
      // relative collapse tolerance
      const auto rel_tol = [&]( Scalar a, Scalar b ) { return abs( a - b ) <= m_epsi * ( Scalar( 1 ) + abs( a ) ); };

      for ( size_t iteration{ 0 }; iteration < m_max_iters; ++iteration )
      {
        // Attempt cubic minimization using existing helper (assumed available).
        Scalar a_j = detail::cubic_minimizer<Scalar>( a_lo, phi_lo, der_lo, a_hi, phi_hi, der_hi );

        // define safeguarded interval for a_j
        Scalar lo_guard = a_lo + m_delta * ( a_hi - a_lo );
        Scalar hi_guard = a_hi - m_sigma * ( a_hi - a_lo );

        // Fallback to bisection if cubic gives bad value or outside guard
        // interval or non-finite
        if ( !( a_j > lo_guard && a_j < hi_guard ) || !std::isfinite( a_j ) ) { a_j = ( a_lo + a_hi ) / 2; }

        Scalar phi_j, der_j;
        eval( a_j, phi_j, der_j );

        // Check Armijo at a_j
        if ( phi_j > ( m_f0 + a_j * m_c1_Df0 ) || phi_j >= phi_lo )
        {
          // a_j is too big — tighten upper bound
          a_hi   = a_j;
          phi_hi = phi_j;
          der_hi = der_j;
        }
        else
        {
          // sufficient decrease holds at a_j — check curvature
          if ( der_j >= ( m_c2_Df0 - m_epsilon_k ) )
          {
            // derivative satisfies curvature condition -> accept a_j
            return a_j;
          }

          // If derivative has same sign as (a_hi - a_lo), move hi to lo
          // (re-orient interval)
          if ( der_j * ( a_hi - a_lo ) >= Scalar( 0 ) )
          {
            a_hi   = a_lo;
            phi_hi = phi_lo;
            der_hi = der_lo;
          }

          // move lower bound up
          a_lo   = a_j;
          phi_lo = phi_j;
          der_lo = der_j;
        }

        // If bracket collapsed sufficiently (relative check), return best point
        // so far
        if ( rel_tol( a_hi, a_lo ) )
        {
          // choose midpoint as a candidate
          Scalar a_mid = ( a_lo + a_hi ) / 2;
          // evaluate midpoint if we still have eval budget
          Scalar phi_mid, der_mid;
          eval( a_mid, phi_mid, der_mid );
          // accept midpoint if curvature satisfied
          if ( phi_mid <= ( m_f0 + a_mid * m_c1_Df0 ) && der_mid >= ( m_c2_Df0 - m_epsilon_k ) ) return a_mid;
          // otherwise return the last tested a_j as best-effort
          return a_j;
        }
      }

      return 0;
    }

  public:
    /**
     * @brief Default constructor — parameters can be tuned after construction.
     *
     * @param c1 Armijo constant (typically 1e-4).
     * @param c2 Curvature constant (typically in (0.1, 0.9), e.g. 0.9).
     */
    HagerZhangLineSearch() {}

    /**
     * @brief Perform Hager–Zhang-style line search.
     *
     * @tparam Callback signature: Scalar callback(const Vector &x, Vector *g)
     *         where callback returns f(x) and writes gradient into *g (if
     * non-null).
     *
     * @param f0        initial function value f(x).
     * @param Df0       initial directional derivative g(x)^T d (must be < 0).
     * @param x         current point.
     * @param d         search direction (descent direction).
     * @param callback  function+gradient evaluator.
     * @param alpha0    initial step (default 1).
     *
     * @return optional step length satisfying Armijo and curvature conditions,
     * or std::nullopt.
     */
    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar           f0,
                Scalar           Df0,
                Vector const &   x,
                Vector const &   d,
                Callback const & callback,
                Scalar           alpha0 = Scalar( 1 ) ) const
    {
      // require descent
      if ( !( Df0 < Scalar( 0 ) ) ) return std::nullopt;

      // initialize internal state
      m_f0     = f0;
      m_Df0    = Df0;
      m_c1_Df0 = m_c1 * Df0;
      m_c2_Df0 = m_c2 * Df0;

      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t n_evals = 0;

      // evaluation wrapper (records count)
      auto eval = [&]( Scalar a, Scalar & f, Scalar & df )
      {
        ++n_evals;
        m_x_new.noalias() = x + a * d;
        f                 = callback( m_x_new, &m_g_new );
        df                = m_g_new.dot( d );
      };

      // initial values
      Scalar alpha      = min( alpha0, m_alpha_max );
      Scalar alpha_prev = Scalar( 0 );
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      // Evaluate at initial alpha
      Scalar phi, der;
      eval( alpha, phi, der );

      // Minimum step threshold to avoid pointless tiny steps
      const Scalar alpha_min_threshold = std::numeric_limits<Scalar>::epsilon() * Scalar( 100 );

      for ( size_t iter{ 0 }; iter < m_max_iters && n_evals < m_max_evals; ++iter )
      {
        // Check Armijo (sufficient decrease) and curvature (Wolfe) conditions
        // Armijo: phi <= f0 + alpha * c1 * Df0
        // Curvature: der >= c2 * Df0
        if ( phi <= ( m_f0 + alpha * m_c1_Df0 ) && der >= ( m_c2_Df0 - m_epsilon_k ) )
        {
          return std::tuple<Scalar, size_t>( alpha, n_evals );
        }

        // Bracketing: phi is too large (doesn't satisfy Armijo) or phi is not
        // decreasing
        if ( phi > ( m_f0 + alpha * m_c1_Df0 ) || ( iter > 0 && phi >= phi_prev ) )
        {
          Scalar new_alpha = zoom( alpha_prev, phi_prev, der_prev, alpha, phi, der, eval );
          if ( new_alpha > 0 ) return std::tuple<Scalar, size_t>( new_alpha, n_evals );
          return std::nullopt;
        }

        // If derivative is already small in magnitude and satisfies modified
        // curvature: Here we check |der| <= -c2*Df0 + epsilon  equivalently der
        // >= c2*Df0 - eps
        if ( abs( der ) <= ( abs( m_c2_Df0 ) + m_epsilon_k ) && der >= ( m_c2_Df0 - m_epsilon_k ) )
        {
          // Accept alpha if derivative is close to satisfying curvature
          return std::tuple<Scalar, size_t>( alpha, n_evals );
        }

        // If derivative non-negative, bracket and zoom
        if ( der >= Scalar( 0 ) )
        {
          Scalar new_alpha = zoom( alpha, phi, der, alpha_prev, phi_prev, der_prev, eval );
          if ( new_alpha > 0 ) return std::tuple<Scalar, size_t>( new_alpha, n_evals );
          return std::nullopt;
        }

        // If we've exhausted max evaluations, break
        if ( n_evals >= m_max_evals ) break;

        // Extrapolation: increase alpha (doubling), but respect maximum allowed
        // step
        Scalar alpha_new = min( Scalar( 2 ) * alpha, m_alpha_max );

        // If step growth stalls or becomes too tiny, fail
        if ( alpha_new <= alpha_min_threshold ) return std::nullopt;

        alpha_prev = alpha;
        phi_prev   = phi;
        der_prev   = der;

        alpha = alpha_new;

        eval( alpha, phi, der );
      }

      return std::nullopt;
    }
  };

  /**
   * @class MoreThuenteLineSearch
   * @brief Implements the Moré–Thuente line search algorithm for Strong Wolfe
   * conditions.
   *
   * This line search algorithm attempts to find a step length α along a descent
   * direction d that satisfies:
   *
   * 1. **Armijo condition (sufficient decrease)**
   *    \f[ φ(α) ≤ φ(0) + c_1 α φ'(0) \f]
   *
   * 2. **Curvature condition (Strong Wolfe)**
   *    \f[ |φ'(α)| ≤ c_2 |φ'(0)| \f]
   *
   * The algorithm works in two phases:
   * - **Extrapolation phase:** tries to find an interval where the conditions
   * may hold.
   * - **Zoom phase:** refines the interval using cubic interpolation until a
   * valid α is found.
   *
   * ### Algorithm Outline
   *
   * 1. Evaluate initial trial α0.
   * 2. If Strong Wolfe conditions are satisfied → return α0.
   * 3. Otherwise:
   *    - Update bracketing interval [α_lo, α_hi].
   *    - If bracketed → enter zoom phase (iterative cubic interpolation).
   *    - If not bracketed → extrapolate α (doubling or cubic extrapolation).
   * 4. Continue until:
   *    - A step satisfying Strong Wolfe is found, or
   *    - Maximum iterations / function evaluations reached.
   *
   * ### Notes / Improvements
   *
   * - Alpha_min should be enforced to avoid numerical underflow.
   * - Zoom phase should be separated for clarity and robustness.
   * - Relative tolerance for bracket collapse is recommended.
   * - Cubic extrapolation for unbracketed phase can accelerate convergence in
   * steep directions.
   * - Robust derivative checks (using absolute values) are safer than signed
   * comparisons.
   *
   * @tparam Scalar floating point type (float, double, long double)
   */

  template <typename Scalar>
  class MoreThuenteLineSearch
  {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1        = 1e-4;  // Armijo
    Scalar m_c2        = 0.9;   // Strong Wolfe (Curvature)
    Scalar m_alpha_max = 10.0;
    Scalar m_alpha_min = 1e-12;  // Passo minimo accettabile
    Scalar m_epsi      = 1e-12;
    size_t m_max_iters = 50;

    mutable Vector m_g_new, m_x_new;

    // --- Stati interni (per chiarezza, gestiti tramite i parametri passati)
    // --- Questi potrebbero essere membri privati, ma li gestiremo nel main
    // operator() per isolamento.

  public:
    template <typename Callback>
    std::optional<std::tuple<Scalar, size_t>>
    operator()( Scalar           f0,
                Scalar           Df0,
                Vector const &   x,
                Vector const &   d,
                Callback const & callback,
                Scalar           alpha0 = 1 ) const
    {
      if ( !( Df0 < 0 ) ) return std::nullopt;

      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t n_evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df )
      {
        ++n_evals;
        m_x_new.noalias() = x + a * d;
        f                 = callback( m_x_new, &m_g_new );
        df                = m_g_new.dot( d );
      };

      Scalar c1_Df0{ m_c1 * Df0 };
      Scalar c2_Df0{ m_c2 * Df0 };

      // ----------------------------------------------------
      // FASE 1: Inizializzazione e ciclo principale (Bracketing & Zoom)
      // ----------------------------------------------------

      // Punti di lavoro
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      Scalar alpha_curr = alpha0;  // Passo di prova corrente
      Scalar phi_curr, der_curr;

      // Intervallo di bracketing iniziale
      Scalar alpha_lo = 0;
      Scalar phi_lo   = f0;
      Scalar der_lo   = Df0;

      Scalar alpha_hi = m_alpha_max;

      bool bracketed = false;

      // Esegui la prima valutazione
      eval( alpha_curr, phi_curr, der_curr );

      for ( size_t k{ 0 }; k < m_max_iters; ++k )
      {
        // Check Strong Wolfe al passo corrente (MT è un raffinamento per SW)
        if ( phi_curr <= f0 + alpha_curr * c1_Df0 && abs( der_curr ) <= -c2_Df0 )
          return std::tuple<Scalar, size_t>( alpha_curr, n_evals );  // Trovato passo accettabile

        // --- Bracketing Logic (Aggiornamento dell'intervallo [alpha_lo,
        // alpha_hi]) ---

        if ( phi_curr > f0 + alpha_curr * c1_Df0 || ( bracketed && phi_curr >= phi_lo ) )
        {
          // Caso A: Violata Armijo o funzione non decrescente
          // L'intervallo [alpha_lo, alpha_curr] contiene il minimo.
          alpha_hi = alpha_curr;
          // Mantieni alpha_lo (potrebbe essere alpha_prev o 0)
          bracketed = true;
        }
        else if ( abs( der_curr ) <= -c2_Df0 )
        {
          // Caso B: Soddisfatta Armijo ma derivata non abbastanza piatta
          // (Strong Wolfe fallita) Se la derivata è troppo negativa (troppo
          // ripida), dobbiamo andare più avanti. Se der_curr < 0, il minimo è
          // oltre alpha_curr.
          alpha_lo = alpha_curr;
          phi_lo   = phi_curr;
          der_lo   = der_curr;
        }
        else if ( der_curr >= 0 )
        {
          // Caso C: Trovata pendenza positiva, l'intervallo è [alpha_prev,
          // alpha_curr]
          alpha_hi = alpha_curr;
          // Mantieni alpha_lo
          bracketed = true;
        }

        // --- Selezione del nuovo passo ---
        Scalar alpha_new;

        if ( bracketed )
        {
          // Zoom phase: usa interpolazione cubica/safeguard tra alpha_lo e
          // alpha_hi
          alpha_new = detail::compute_step( alpha_lo, phi_lo, der_lo, alpha_hi, phi_curr,
                                            der_curr,  // Nota: usa phi_curr e der_curr per alpha_hi (il passo
                                                       // "cattivo")
                                            alpha_prev, phi_prev, der_prev, m_alpha_max, m_alpha_max,
                                            false  // is_bracketing = false
          );
        }
        else
        {
          // Extrapolation phase: usa interpolazione o raddoppia il passo

          // Opzione 1: Raddoppia
          alpha_new = 2 * alpha_curr;

          // Opzione 2: Interpolazione/Estrapolazione (più avanzata, ma
          // rischiosa) alpha_new = detail::compute_step(
          //     alpha_lo, phi_lo, der_lo,
          //     alpha_curr, phi_curr, der_curr,
          //     alpha_prev, phi_prev, der_prev,
          //     m_alpha_max, m_alpha_max,
          //     true // is_bracketing = true
          // );
        }

        // --- Salvaguardia finale ---
        if ( alpha_new <= 0 || alpha_new >= m_alpha_max ) alpha_new = min( 2 * alpha_curr, m_alpha_max );
        if ( abs( alpha_hi - alpha_lo ) < m_epsi )
          return std::tuple<Scalar, size_t>( alpha_curr, n_evals );  // Convergenza dell'intervallo

        // Aggiorna stato per la prossima iterazione
        alpha_prev = alpha_curr;
        phi_prev   = phi_curr;
        der_prev   = der_curr;

        alpha_curr = alpha_new;
        eval( alpha_curr, phi_curr, der_curr );
      }

      return std::nullopt;  // Raggiunto il massimo delle iterazioni
    }
  };

}  // namespace Utils

#endif

#endif

//
// eof: Utils_Linesearch.hh
