/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003-2025                                                 |
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
 * @file Utils_AlgoHNewton.hh
 * @brief Header-only implementation of the HNewton (Hybrid Newton) algorithm
 *        for finding roots of functions with derivative information.
 *
 * The algorithm combines:
 * 1. Newton's method when derivative is available
 * 2. Quadratic/inverse interpolation to accelerate convergence
 * 3. Secant step as fallback
 * 4. Bisection when needed for stability
 *
 * The first iteration is simply a secant step.
 * Starting with the second iteration, three steps are taken in each iteration:
 *   a. First two steps are either quadratic interpolation or cubic inverse interpolation.
 *   b. The third step is a double-size secant step.
 * If the diameter of the enclosing interval obtained after those three steps is
 * larger than (b-a)/2, then an additional bisection step will be taken.
 *
 * Finds either an exact solution or an approximate solution of the equation f(x)=0
 * in the interval [a,b].
 *
 * Template: Real (float, double, etc.)
 */

#pragma once

#ifndef UTILS_ALGO_HNEWTON_dot_HH
#define UTILS_ALGO_HNEWTON_dot_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Utils.hh"
#include "Utils_fmt.hh"

namespace Utils
{
  using std::abs;
  using std::ceil;
  using std::log2;
  using std::max;
  using std::min;
  using std::pow;
  using std::signbit;
  using std::swap;

  /*!
   * \addtogroup Zeros
   * @{
   */

  /**
   * @class AlgoHNewton_base_fun
   * @brief Abstract base class for mathematical functions used in the zero search algorithm.
   *
   * This class serves as a base interface for user-defined functions that can be evaluated.
   * It allows for the implementation of the numerical method to find the solution of the
   * one-dimensional equation f(x) = 0. Users must inherit from this class and implement
   * the virtual methods to define their specific functions.
   *
   * @tparam Real Numeric type representing the data type of the function's input and output,
   *             such as float, double, etc.
   *
   * @see AlgoHNewton
   */
  template <typename Real> class AlgoHNewton_base_fun
  {
  public:
    virtual ~AlgoHNewton_base_fun() = default;

    /**
     * @brief Evaluate the function f(x)
     * @param x Point at which to evaluate f(x)
     * @return Value of f(x)
     */
    virtual Real eval( Real x ) const = 0;

    /**
     * @brief Evaluate the derivative f'(x)
     * @param x Point at which to evaluate f'(x)
     * @return Value of f'(x)
     */
    virtual Real D( Real x ) const = 0;
  };

  /**
   * @class AlgoHNewton_fun
   * @brief Adapter class for C-style functions
   *
   * This class wraps a pair of function pointers (function and its derivative)
   * into the AlgoHNewton_base_fun interface.
   *
   * @tparam Real Numeric type
   * @tparam PFUN Type of function pointer for f(x)
   * @tparam PFUN_D Type of function pointer for f'(x)
   */
  template <typename Real, typename PFUN, typename PFUN_D> class AlgoHNewton_fun : public AlgoHNewton_base_fun<Real>
  {
    PFUN   m_fun;    ///< Function pointer for f(x)
    PFUN_D m_fun_D;  ///< Function pointer for f'(x)

  public:
    /**
     * @brief Construct from function pointers
     * @param f Function pointer for f(x)
     * @param Df Function pointer for f'(x)
     */
    explicit AlgoHNewton_fun( PFUN f, PFUN_D Df ) : m_fun( f ), m_fun_D( Df ) {}

    Real eval( Real x ) const override { return m_fun( x ); }
    Real D( Real x ) const override { return m_fun_D( x ); }
  };

  /**
   * @class AlgoHNewton
   * @brief Hybrid Newton algorithm for solving f(x)=0 with derivative information
   *
   * This class implements a hybrid algorithm combining Newton's method with
   * interpolation techniques for robust and efficient root finding.
   *
   * Usage example:
   * @code{.cpp}
   * class MyFunction : public AlgoHNewton_base_fun<double> {
   * public:
   *   double eval(double x) const override { return x*x - 2; }
   *   double D(double x) const override { return 2*x; }
   * };
   *
   * MyFunction f;
   * AlgoHNewton<double> solver;
   * double root = solver.eval(0.0, 3.0, &f);
   * @endcode
   *
   * @tparam Real Numeric type (float, double, etc.)
   */
  template <typename Real> class AlgoHNewton
  {
    using Integer = int;

  private:
    // Algorithm parameters
    Real m_tolerance{ pow( machine_eps<Real>(), Real( 2. / 3. ) ) };  ///< Convergence tolerance
    bool m_converged{ false };                                        ///< Convergence flag
    Real m_kappa{ 0.05 };                                             ///< Safety parameter for step bounds

    // Iteration state
    Real m_a{ 0 }, m_fa{ 0 };  ///< Left endpoint and function value
    Real m_b{ 0 }, m_fb{ 0 };  ///< Right endpoint and function value
    Real m_c{ 0 }, m_fc{ 0 };  ///< Intermediate point and function value
    Real m_ba{ 0 };            ///< Current interval width (b - a)

    // Function interface
    AlgoHNewton_base_fun<Real> const * m_function{ nullptr };  ///< Pointer to user function

    // Counters
    Integer         m_max_iteration{ 200 };         ///< Maximum allowed iterations
    mutable Integer m_iteration_count{ 0 };         ///< Iterations used in last computation
    mutable Integer m_fun_evaluation_count{ 0 };    ///< Function evaluations count
    mutable Integer m_fun_D_evaluation_count{ 0 };  ///< Derivative evaluations count

    /**
     * @brief Evaluate function at x with counter increment
     * @param x Evaluation point
     * @return f(x)
     */
    Real evaluate( Real x ) const
    {
      ++m_fun_evaluation_count;
      return m_function->eval( x );
    }

    /**
     * @brief Evaluate derivative at x with counter increment
     * @param x Evaluation point
     * @return f'(x)
     */
    Real evaluate_D( Real x ) const
    {
      ++m_fun_D_evaluation_count;
      return m_function->D( x );
    }

    /**
     * @brief Quadratic interpolation (direct)
     *
     * Uses quadratic interpolation of f(x) at a, b, and c to get an approximate root.
     *
     * @return Interpolated root estimate in normalized coordinates [0,1] relative to [a,b]
     */
    Real p_zero2() const;

    /**
     * @brief Quadratic inverse interpolation
     *
     * Uses quadratic inverse interpolation of f(x) at a, b, and c to get an approximate root.
     * Rewritten using divided differences for numerical stability.
     *
     * @return Interpolated root estimate in normalized coordinates [0,1] relative to [a,b]
     */
    Real invp_zero2() const;

    /**
     * @brief Core algorithm implementation
     *
     * Implements the hybrid Newton algorithm as described in the class documentation.
     *
     * @return Approximate root
     */
    Real eval();

    /**
     * @brief Initialize and run algorithm with interval [a,b]
     * @param a Left endpoint
     * @param b Right endpoint
     * @return Root approximation
     */
    Real eval( Real a, Real b );

    /**
     * @brief Initialize and run algorithm with expandable interval
     *
     * If no sign change is detected in the initial interval, tries to expand
     * it within the bounds [amin, bmax] until a sign change is found.
     *
     * @param a Initial left endpoint
     * @param b Initial right endpoint
     * @param amin Minimum allowed left endpoint
     * @param bmax Maximum allowed right endpoint
     * @return Root approximation
     */
    Real eval( Real a, Real b, Real amin, Real bmax );

  public:
    /**
     * @brief Default constructor
     */
    AlgoHNewton() = default;

    /**
     * @brief Destructor
     */
    ~AlgoHNewton() = default;

    // Public interface
    // ================

    /**
     * @brief Set maximum number of iterations
     * @param mit Maximum number of iterations (must be > 0)
     */
    void set_max_iterations( Integer mit )
    {
      UTILS_ASSERT( mit > 0, "AlgoHNewton::set_max_iterations({}) argument must be >0\n", mit );
      m_max_iteration = mit;
    }

    /**
     * @brief Set tolerance based on reference value
     *
     * The tolerance is computed as: eps + 2*|B|*eps, where eps is machine epsilon.
     * This ensures the tolerance scales appropriately with the magnitude of B.
     *
     * @param B Reference value for tolerance scaling
     */
    void set_tolerance( Real B )
    {
      Real eps{ 2 * machine_eps<Real>() };
      m_tolerance = eps + 2 * abs( B ) * eps;
    }

    /**
     * @brief Solve f(x)=0 in interval [a,b]
     * @param a Lower bound of search interval
     * @param b Upper bound of search interval
     * @param fun Pointer to user function
     * @return Root approximation
     */
    Real eval( Real a, Real b, AlgoHNewton_base_fun<Real> const * fun )
    {
      m_function = fun;
      return this->eval( a, b );
    }

    /**
     * @brief Solve f(x)=0 with automatic interval expansion
     *
     * If no sign change is found in [a,b], expands the interval within [amin, bmax]
     * until a sign change is detected.
     *
     * @param a Initial guess interval lower bound
     * @param b Initial guess interval upper bound
     * @param amin Lower bound for interval expansion
     * @param bmax Upper bound for interval expansion
     * @param fun Pointer to user function
     * @return Root approximation
     */
    Real eval( Real a, Real b, Real amin, Real bmax, AlgoHNewton_base_fun<Real> const * fun )
    {
      m_function = fun;
      return this->eval( a, b, amin, bmax );
    }

    /**
     * @brief Get number of iterations used in last computation
     * @return Iteration count
     */
    Integer used_iter() const { return m_iteration_count; }

    /**
     * @brief Get number of function evaluations in last computation
     * @return Function evaluation count
     */
    Integer num_fun_eval() const { return m_fun_evaluation_count; }

    /**
     * @brief Get number of derivative evaluations in last computation
     * @return Derivative evaluation count
     */
    Integer num_fun_D_eval() const { return m_fun_D_evaluation_count; }

    /**
     * @brief Get current tolerance
     * @return Tolerance value
     */
    Real tolerance() const { return m_tolerance; }

    /**
     * @brief Check if last computation converged
     * @return true if converged, false otherwise
     */
    bool converged() const { return m_converged; }

    /**
     * @brief Get left endpoint of final interval
     * @return Left endpoint a
     */
    Real a() const { return m_a; }

    /**
     * @brief Get right endpoint of final interval
     * @return Right endpoint b
     */
    Real b() const { return m_b; }

    /**
     * @brief Get function value at left endpoint
     * @return f(a)
     */
    Real fa() const { return m_fa; }

    /**
     * @brief Get function value at right endpoint
     * @return f(b)
     */
    Real fb() const { return m_fb; }
  };

  // =================================================================
  // TEMPLATE IMPLEMENTATIONS (HEADER-ONLY)
  // =================================================================

  /**
   * @brief Initialize and run algorithm with fixed interval [a,b]
   *
   * Checks if a solution exists in [a,b] by verifying f(a)*f(b) ≤ 0.
   * If no sign change, returns a (left endpoint).
   */
  template <typename Real> Real AlgoHNewton<Real>::eval( Real a, Real b )
  {
    m_iteration_count        = 0;
    m_fun_evaluation_count   = 0;
    m_fun_D_evaluation_count = 0;

    m_a  = a;
    m_fa = evaluate( m_a );
    m_b  = b;
    m_fb = evaluate( m_b );

    // Check if a solution can exist (sign change condition)
    if ( m_fa * m_fb > 0 ) return m_a;  // No sign change in interval

    return eval();  // Start main algorithm
  }

  /**
   * @brief Initialize and run algorithm with expandable interval
   *
   * Tries to enlarge interval if no sign change is initially detected,
   * within the limits [amin, bmax].
   */
  template <typename Real> Real AlgoHNewton<Real>::eval( Real a, Real b, Real amin, Real bmax )
  {
    m_iteration_count        = 0;
    m_fun_evaluation_count   = 0;
    m_fun_D_evaluation_count = 0;

    m_a  = a;
    m_fa = evaluate( m_a );
    m_b  = b;
    m_fb = evaluate( m_b );

    // Try to enlarge interval until sign change is found
    while ( m_fa * m_fb > 0 )
    {
      if ( is_finite( m_fa ) && m_a > amin )
      {
        // Expand leftward
        m_a -= m_b - m_a;
        m_fa = evaluate( m_a );
      }
      else if ( is_finite( m_fb ) && m_b < bmax )
      {
        // Expand rightward
        m_b += m_b - m_a;
        m_fb = evaluate( m_b );
      }
      else
      {
        break;  // Cannot expand further
      }
    }

    return eval();  // Start main algorithm
  }

  /**
   * @brief Quadratic interpolation (direct method)
   *
   * Constructs quadratic polynomial through (a,fa), (b,fb), (c,fc)
   * and finds its root in normalized coordinates [0,1] relative to [a,b].
   * Uses scaling for numerical stability.
   */
  template <typename Real> Real AlgoHNewton<Real>::p_zero2() const
  {
    // Differences for divided differences
    Real D1{ m_fb - m_fa };              // f[b,a]
    Real theta{ ( m_c - m_a ) / m_ba };  // Normalized position of c
    Real D2{ m_fc - m_fa };              // f[c,a]

    // Quadratic coefficients (scaled for stability)
    Real A{ D2 / theta - D1 };          // Quadratic coefficient
    Real B{ D1 * theta - D2 / theta };  // Linear coefficient
    Real C{ m_fa * ( theta - 1 ) };     // Constant term

    // Scale coefficients to avoid overflow/underflow
    Real M{ max( max( abs( A ), abs( B ) ), abs( C ) ) };
    if ( C < 0 ) M = -M;  // Preserve sign
    A /= M;
    B /= M;
    C /= M;

    // Solve quadratic equation: A*x² + B*x + C = 0
    Real D{ B * B - 4 * A * C };     // Discriminant
    if ( D < 0 ) return -m_fa / D1;  // Fallback to secant if no real roots

    D = sqrt( D );
    // Choose numerically stable formula based on sign of B
    return B < 0 ? 2 * C / ( D - B ) : ( D + B ) / ( 2 * A );
  }

  /**
   * @brief Quadratic inverse interpolation
   *
   * Uses divided differences to construct inverse interpolating polynomial.
   * Interpolates x as a function of f at points (fa,0), (fb,1), (fc,(c-a)/(b-a)).
   * Returns normalized root estimate in [0,1].
   */
  template <typename Real> Real AlgoHNewton<Real>::invp_zero2() const
  {
    // Function values (f coordinates)
    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ m_fc };

    // x coordinates (normalized)
    Real D0{ 0 };                     // x at f=fa (normalized a = 0)
    Real D1{ 1 };                     // x at f=fb (normalized b = 1)
    Real D2{ ( m_c - m_a ) / m_ba };  // x at f=fc (normalized c position)

    // Divided differences
    Real D01{ ( D0 - D1 ) / ( x0 - x1 ) };     // First order: f[a,b]
    Real D12{ ( D1 - D2 ) / ( x1 - x2 ) };     // First order: f[b,c]
    Real D012{ ( D01 - D12 ) / ( x0 - x2 ) };  // Second order: f[a,b,c]

    // Newton polynomial evaluation at f=0
    Real O1{ 0 - x0 };           // (0 - fa)
    Real O2{ ( 0 - x1 ) * O1 };  // (0 - fb)*(0 - fa)

    Real P0{ D0 };              // Constant term
    Real P1{ P0 + D01 * O1 };   // Linear term
    Real P2{ P1 + D012 * O2 };  // Quadratic term

    // Safety check
    UTILS_ASSERT(
      is_finite( P2 ),
      "AlgoHNewton::invp_zero2(): computed NaN or Inf\n"
      "a={} f(a)={}\n"
      "b={} f(b)={}\n"
      "c={} f(c)={}\n",
      m_a,
      m_fa,
      m_b,
      m_fb,
      m_c,
      m_fc );

    return P2;  // Normalized root estimate
  }

  /**
   * @brief Core hybrid Newton algorithm implementation
   *
   * Implements the algorithm as described in the header documentation:
   * 1. Check for trivial solutions at endpoints
   * 2. Handle infinite function values with bisection
   * 3. Main loop with Newton steps, interpolation, and interval updates
   * 4. Convergence checks and final approximation
   */
  template <typename Real> Real AlgoHNewton<Real>::eval()
  {
    // Check for exact solutions at endpoints
    m_converged = m_fa == 0;
    if ( m_converged ) return m_a;
    m_converged = m_fb == 0;
    if ( m_converged ) return m_b;

    // Initial tolerance based on machine precision and interval magnitude
    m_tolerance = machine_eps<Real>() + 2 * max( abs( m_a ), abs( m_b ) ) * machine_eps<Real>();

    // Handle infinite/NaN function values at endpoints using bisection
    while ( !( is_finite( m_fa ) && is_finite( m_fb ) ) )
    {
      ++m_iteration_count;
      m_c  = ( m_a + m_b ) / 2;  // Bisection step
      m_fc = evaluate( m_c );

      m_converged = m_fc == 0;  // Exact solution found
      if ( m_converged )
      {
        m_a = m_b = m_c;
        m_fa = m_fb = 0;
        return m_c;
      }  // collapse interval

      // Update interval based on sign change
      if ( m_fa * m_fc < 0 )  // Root in [a,c]
      {
        m_b  = m_c;
        m_fb = m_fc;
      }
      else if ( m_fb * m_fc < 0 )  // Root in [c,b]
      {
        m_a  = m_c;
        m_fa = m_fc;
      }
      else
      {
        UTILS_ERROR(
          "AlgoHNewton::eval(): cannot determine which subinterval contains root\n"
          "a={} fa={}\n"
          "b={} fb={}\n"
          "c={} fc={}\n",
          m_a,
          m_fa,
          m_b,
          m_fb,
          m_c,
          m_fc );
      }

      // Check convergence after bisection
      m_converged = ( m_b - m_a ) <= m_tolerance;
      if ( m_converged ) return abs( m_fb ) < abs( m_fa ) ? m_b : m_a;
    }

    // Main tolerance for the algorithm
    m_tolerance = 10 * machine_eps<Real>() * ( 1 + max( abs( m_a ), abs( m_b ) ) );

    Real abs_fa{ abs( m_fa ) };
    Real abs_fb{ abs( m_fb ) };

    // Main iteration loop
    while ( ++m_iteration_count < m_max_iteration )
    {
      m_ba        = m_b - m_a;            // Current interval width
      m_converged = m_ba <= m_tolerance;  // Convergence check
      if ( m_converged ) break;

      // Epsilon for checking significant differences in function values
      Real epsi{ max( Real( 1 ), max( abs_fa, abs_fb ) ) * sqrt_machine_eps<Real>() };
      auto diff = [epsi]( Real a, Real b ) -> bool { return abs( a - b ) > epsi; };

      // Newton step (choose endpoint with smaller |f| for stability)
      bool a_small = abs_fa < abs_fb;
      if ( ( m_iteration_count % 3 ) == 0 ) a_small = !a_small;
      if ( a_small )
      {
        Real fa_D{ evaluate_D( m_a ) };
        m_c = m_a - m_fa / fa_D;  // Newton iteration from a
      }
      else
      {
        Real fb_D{ evaluate_D( m_b ) };
        m_c = m_b - m_fb / fb_D;  // Newton iteration from b
      }

      // Ensure Newton step stays inside interval (with safety margin)
      Real delta{ m_ba * m_kappa };  // Safety margin
      if ( Real m_bb = m_b - delta; m_c > m_bb )
        m_c = m_bb;  // Fallback
      else if ( Real m_aa = m_a + delta; m_c < m_aa )
        m_c = m_aa;  // Fallback

      // Evaluate at new point
      m_fc        = evaluate( m_c );
      m_converged = m_fc == 0;
      if ( m_converged )
      {
        m_a = m_b = m_c;
        m_fa = m_fb = 0;
        return m_c;
      }  // Exact solution found

      // Update interval if Newton step gives proper bracketing
      if ( m_fa * m_fc < 0 )  // Root in [a,c]
      {
        if ( m_c - m_a <= m_b - m_c )  // Choose smaller subinterval
        {
          m_b    = m_c;
          m_fb   = m_fc;
          abs_fb = abs( m_fb );
          continue;  // Skip to next iteration
        }
      }
      else  // Root in [c,b]
      {
        if ( m_c - m_a >= m_b - m_c )  // Choose smaller subinterval
        {
          m_a    = m_c;
          m_fa   = m_fc;
          abs_fa = abs( m_fa );
          continue;  // Skip to next iteration
        }
      }

      // Interpolation step (use inverse if all function values are distinct)
      bool all_diff{ diff( m_fa, m_fb ) && diff( m_fa, m_fc ) && diff( m_fb, m_fc ) };
      Real x{ all_diff ? invp_zero2() : p_zero2() };  // Normalized root estimate

      // Clamp to safe region within interval
      if ( x < m_kappa )
        x = m_kappa;
      else if ( x > 1 - m_kappa )
        x = 1 - m_kappa;

      // Compute new point d using interpolation estimate
      Real d{ m_a + x * m_ba };
      Real fd{ evaluate( d ) };
      m_converged = fd == 0;
      if ( m_converged )
      {
        m_a = m_b = d;
        m_fa = m_fb = 0;
        return d;
      }  // Exact solution found

      // Ensure c < d for consistent ordering
      if ( m_c > d )
      {
        swap( m_c, d );
        swap( m_fc, fd );
      }

      // Update bracketing interval based on sign changes
      if ( m_fc * fd < 0 )  // Root in [c,d]
      {
        m_a    = m_c;
        m_fa   = m_fc;
        abs_fa = abs( m_fa );
        m_b    = d;
        m_fb   = fd;
        abs_fb = abs( fd );
      }
      else
      {
        if ( m_fa * m_fc > 0 )  // Root in [d,b]
        {
          m_a    = d;
          m_fa   = fd;
          abs_fa = abs( fd );
        }
        else  // Root in [a,c]
        {
          m_b    = m_c;
          m_fb   = m_fc;
          abs_fb = abs( m_fb );
        }
      }
    }

    // Return best approximation (endpoint with smaller |f|)
    return abs_fa < abs_fb ? m_a : m_b;
  }

  /*! @} */  // End of Zeros group

}  // namespace Utils

#endif  // UTILS_ALGO_HNEWTON_dot_HH
