/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003-2022                                                 |
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
// file: Utils_zeros.hh
//

#pragma once

#ifndef UTILS_ZEROS_dot_HH
#define UTILS_ZEROS_dot_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils
{

  using std::abs;
  using std::isfinite;
  using std::numeric_limits;
  using std::pow;

  /*!
   * \defgroup Zeros Root Finding Algorithms
   * \brief A collection of iterative methods for finding roots of nonlinear equations.
   *
   * This module provides various numerical methods for solving \f$ f(x) = 0 \f$,
   * ranging from classical methods (Newton, Halley) to high-order optimal methods
   * (up to order 32). All methods require function evaluations and derivatives.
   */

  /*!
   * \ingroup Zeros
   * \class Zeros_base_fun
   * \brief Abstract base class for defining mathematical functions for root-finding.
   *
   * This class defines the interface for user-defined functions that can be evaluated
   * and differentiated up to the third derivative. It is used by all root-finding
   * algorithms in the Zeros class.
   *
   * \tparam Real Numeric type (float, double, etc.)
   *
   * \note Users must implement all pure virtual methods in derived classes.
   *
   * \see Zeros
   *
   * \example
   * \code{.cpp}
   * class MyFunction : public Zeros_base_fun<double> {
   * public:
   *     double eval(double x) const override { return x*x - 2.0; }
   *     double eval_D(double x) const override { return 2.0*x; }
   *     double eval_DD(double x) const override { return 2.0; }
   *     double eval_DDD(double x) const override { return 0.0; }
   * };
   * \endcode
   */
  template <typename Real> class Zeros_base_fun
  {
  public:
    virtual ~Zeros_base_fun() = default;

    //! Evaluate \f$ f(x) \f$
    virtual Real eval( Real x ) const = 0;

    //! Evaluate \f$ f'(x) \f$ (first derivative)
    virtual Real eval_D( Real x ) const = 0;

    //! Evaluate \f$ f''(x) \f$ (second derivative)
    virtual Real eval_DD( Real x ) const = 0;

    //! Evaluate \f$ f'''(x) \f$ (third derivative)
    virtual Real eval_DDD( Real x ) const = 0;

    //! Shorthand for eval()
    Real operator()( Real x ) const { return this->eval( x ); }

    //! Shorthand for eval_D()
    Real D( Real x ) const { return this->eval_D( x ); }

    //! Shorthand for eval_DD()
    Real DD( Real x ) const { return this->eval_DD( x ); }

    //! Shorthand for eval_DDD()
    Real DDD( Real x ) const { return this->eval_DDD( x ); }
  };

  /*!
   * \ingroup Zeros
   * \class Zeros
   * \brief Solver for finding roots of nonlinear equations \f$ f(x) = 0 \f$.
   *
   * This class implements multiple iterative methods for root finding, including:
   * - Newton-Raphson (order 2) \cite NR
   * - Chebyshev's method (order 3) \cite Chebyshev
   * - Halley's method (order 3) \cite Halley
   * - Optimal methods of orders 4, 8, 16, and 32 \cite Varona2022
   *
   * The higher-order methods (4, 8, 16, 32) are based on the work of Juan Luis Varona
   * and provide optimal efficiency according to the Kung-Traub conjecture.
   *
   * \tparam Real Numeric type (float, double, etc.)
   *
   * \note All methods require the function and its derivatives via Zeros_base_fun.
   *
   * \references
   * \anchor NR <b>[NR]</b> Newton's method. Wikipedia. https://en.wikipedia.org/wiki/Newton%27s_method\n
   * \anchor Chebyshev <b>[Chebyshev]</b> Chebyshev's method. Wikipedia.
   * https://en.wikipedia.org/wiki/Chebyshev%27s_method\n
   * \anchor Halley <b>[Halley]</b> Halley's method. Wikipedia. https://en.wikipedia.org/wiki/Halley%27s_method\n
   * \anchor Varona2022 <b>[Varona2022]</b> J.L. Varona. "An Optimal Thirty-Second-Order Iterative Method for Solving
   * Nonlinear Equations and a Conjecture". Qualitative Theory of Dynamical Systems, 2022.
   * https://doi.org/10.1007/s12346-022-00572-3
   *
   * \example
   * \code{.cpp}
   * Zeros<double> solver;
   * solver.set_tolerance(1e-12);
   * solver.set_max_iterations(50);
   *
   * MyFunction f;
   * double root = solver.solve_Newton(1.5, &f);
   * if (solver.converged()) {
   *     fmt::print("Root found: {} after {} iterations\n",
   *                root, solver.used_iter());
   * }
   * \endcode
   */
  template <typename Real> class Zeros
  {
    using Integer = int;

    Integer m_max_fun_evaluation{ 200 };  //!< Max function evaluations
    Integer m_max_iteration{ 100 };       //!< Max iterations
    Real    m_tolerance{ pow( machine_eps<Real>(), Real( 2. / 3. ) ) };
    bool    m_converged{ false };

    mutable Integer m_iteration_count{ 0 };
    mutable Integer m_fun_evaluation_count{ 0 };

    // Static failure value
    static constexpr Real FAILURE_VALUE = numeric_limits<Real>::quiet_NaN();

    // Weight functions for optimal methods (Varona 2022) with safety checks
    static Real Q( Real t )
    {
      if ( !is_finite( t ) || abs( t ) > Real( 1e10 ) ) return Real( 1.0 );
      return 1 + 2 * t;
    }

    static Real W( Real t, Real s )
    {
      if ( !is_finite( t ) || !is_finite( s ) || abs( t ) > Real( 1e10 ) || abs( s ) > Real( 1e10 ) )
        return Real( 1.0 );

      Real t2{ t * t };
      return t2 * ( 1 - 4 * t ) + ( 4 * s + 2 ) * t + s + 1;
    }

    static Real H( Real t, Real s, Real u )
    {
      if (
        !is_finite( t ) || !is_finite( s ) || !is_finite( u ) || abs( t ) > Real( 1e5 ) || abs( s ) > Real( 1e5 ) ||
        abs( u ) > Real( 1e5 ) )
        return Real( 1.0 );

      Real t1  = t * t;
      Real t2  = t1 * t1;
      Real t8  = s * s;
      Real t17 = s * t8;
      Real t23 = 2 * u;
      return ( ( 8 * u + 6 * t2 + 4 ) * s - ( 6 * t8 + 4 * ( s + u + 1 ) ) * t1 + 2 * t8 - 4 * t17 + t23 + 2 ) * t +
             t1 * ( t8 + s + u + 1 ) + ( 1 - 3 * t2 + t23 ) * s + u - t17 + 1;
    }

    static Real J( Real t, Real s, Real u, Real v )
    {
      if (
        !is_finite( t ) || !is_finite( s ) || !is_finite( u ) || !is_finite( v ) || abs( t ) > Real( 1e3 ) ||
        abs( s ) > Real( 1e3 ) || abs( u ) > Real( 1e3 ) || abs( v ) > Real( 1e3 ) )
        return Real( 1.0 );

      Real t1  = s * s;
      Real t2  = t1 * t1;
      Real t17 = t * t;
      Real t22 = u * u;
      Real t32 = t17 * t17;
      Real t34 = t * t32;
      Real t37 = t * t17;
      Real t46 = 1 + v;
      Real t65 = u + 1 + v;
      Real t76 = ( 2 - 2 * t22 + u + 4 * v ) * u;

      constexpr Real r1_2{ static_cast<Real>( 1.0 / 2.0 ) };
      constexpr Real r1_4{ static_cast<Real>( 1.0 / 4.0 ) };
      constexpr Real r2_3{ static_cast<Real>( 2.0 / 3.0 ) };
      constexpr Real r3_2{ static_cast<Real>( 3.0 / 2.0 ) };
      constexpr Real r3_4{ static_cast<Real>( 3.0 / 4.0 ) };
      constexpr Real r3_8{ static_cast<Real>( 3.0 / 8.0 ) };
      constexpr Real r5_8{ static_cast<Real>( 5.0 / 8.0 ) };
      return ( 2 * t - 1 ) * ( 2 + 5 * t ) * u * t * t2 + ( 4 * t + 1 ) * u * s * t2 +
             ( u * t22 - 2 * u * v - u - v - 1 ) * ( 4 * t17 + 3 * t + 1 ) * ( t - 1 ) -
             8 *
               ( t22 * ( t17 / 2 - r1_4 ) +
                 u * ( t17 * t32 - r5_8 * t34 - r3_4 * t32 + r3_8 * t37 + r3_4 * t17 - t / 8 - r1_4 ) +
                 r3_4 * t46 * ( t + r1_2 ) * ( t - r2_3 ) ) *
               t * t1 +
             4 *
               ( t22 * ( -r3_2 * t - r1_4 ) + u * ( t34 - t32 - r3_2 * t37 + t17 / 4 - t - r1_4 ) -
                 t46 * ( t + r3_4 ) ) *
               s * t1 +
             ( 1 + v + t65 * t17 - 4 * t65 * t37 - 3 * t65 * t32 + 6 * t65 * t34 + t76 + 4 * ( 1 + v + t76 ) * t ) * s;
    }

    bool              m_trace = false;
    std::vector<Real> m_trace_values;

    // Helper function to check if result is valid
    bool is_result_valid( Real const result ) const { return m_converged && Utils::is_finite( result ); }

    // Common convergence check with trace
    bool check_convergence( Real const f, Real const x )
    {
      if ( m_trace ) m_trace_values.push_back( x );
      m_converged = abs( f ) < m_tolerance;
      return m_converged;
    }

    // Check iteration limits
    bool check_limits() const
    {
      return m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation;
    }

    // Safe division with check
    bool safe_divide( Real & result, Real a, Real b ) const
    {
      if ( abs( b ) < machine_eps<Real>() || !Utils::is_finite( b ) ) { return false; }
      result = a / b;
      return Utils::is_finite( result );
    }

  public:
    Zeros()  = default;
    ~Zeros() = default;

    //! Set maximum number of iterations
    void set_max_iterations( Integer mit )
    {
      UTILS_ASSERT( mit > 0, "Zeros::set_max_iterations({}) argument must be >0\n", mit );
      m_max_iteration = mit;
    }

    //! Set maximum function evaluations
    void set_max_fun_evaluation( Integer mfev )
    {
      UTILS_ASSERT( mfev > 0, "Zeros::set_max_fun_evaluation({}) argument must be >0\n", mfev );
      m_max_fun_evaluation = mfev;
    }

    //! Set tolerance \f$ \epsilon \f$ (stop when \f$ |f(x)| < \epsilon \f$)
    void set_tolerance( Real tol )
    {
      UTILS_ASSERT( tol > 0, "Zeros::set_tolerance({}) argument must be >0\n", tol );
      m_tolerance = tol;
    }

    void set_trace( bool yn ) { m_trace = yn; }

    std::vector<Real> const & trace_values() const { return m_trace_values; }

    //! Check if the last result is valid
    bool is_valid_result( Real result ) const { return is_result_valid( result ); }

    /*!
     * \brief Newton-Raphson method (order 2)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Newton( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;

        Real delta;
        if ( !safe_divide( delta, f, df ) ) break;

        x -= delta;
        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Chebyshev's method (order 3)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Chebyshev( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;
        Real ddf = fun->eval_DD( x );
        ++m_fun_evaluation_count;

        Real ratio;
        if ( !safe_divide( ratio, f, df ) ) break;

        Real adjustment;
        if ( !safe_divide( adjustment, f * ddf, 2 * df * df ) ) adjustment = 0;

        x -= ratio * ( 1 + adjustment );
        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Halley's method (order 3)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Halley( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;
        Real ddf = fun->eval_DD( x );
        ++m_fun_evaluation_count;

        Real ratio;
        if ( !safe_divide( ratio, f, df ) ) break;

        Real denominator;
        if ( !safe_divide( denominator, f * ddf, 2 * df * df ) ) denominator = 1;  // Fallback to Newton

        x -= ratio / ( 1 - denominator );
        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Optimal 4th-order method (Varona 2022)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Order4( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;

        Real delta;
        if ( !safe_divide( delta, f, df ) ) break;

        Real y = x - delta;
        if ( !Utils::is_finite( y ) ) break;

        Real fy = fun->eval( y );
        ++m_fun_evaluation_count;

        // Check convergence at y
        if ( abs( fy ) < m_tolerance )
        {
          m_converged = true;
          return y;
        }

        Real t;
        if ( !safe_divide( t, fy, f ) ) break;

        Real weight = Q( t );
        x           = y - weight * ( fy / df );

        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Optimal 8th-order method (Varona 2022)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Order8( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;

        Real delta;
        if ( !safe_divide( delta, f, df ) ) break;

        Real y = x - delta;
        if ( !Utils::is_finite( y ) ) break;

        Real fy = fun->eval( y );
        ++m_fun_evaluation_count;

        // Check convergence at y
        if ( abs( fy ) < m_tolerance )
        {
          m_converged = true;
          return y;
        }

        Real t;
        if ( !safe_divide( t, fy, f ) ) break;

        Real z = y - Q( t ) * ( fy / df );
        if ( !Utils::is_finite( z ) ) break;

        Real fz = fun->eval( z );
        ++m_fun_evaluation_count;

        // Check convergence at z
        if ( abs( fz ) < m_tolerance )
        {
          m_converged = true;
          return z;
        }

        Real s;
        if ( !safe_divide( s, fz, fy ) ) break;

        Real weight = W( t, s );
        x           = z - weight * ( fz / df );

        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Optimal 16th-order method (Varona 2022)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Order16( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;

        Real delta;
        if ( !safe_divide( delta, f, df ) ) break;

        Real y = x - delta;
        if ( !Utils::is_finite( y ) ) break;

        Real fy = fun->eval( y );
        ++m_fun_evaluation_count;

        // Check convergence at y
        if ( abs( fy ) < m_tolerance )
        {
          m_converged = true;
          return y;
        }

        Real t;
        if ( !safe_divide( t, fy, f ) ) break;

        Real z = y - Q( t ) * ( fy / df );
        if ( !Utils::is_finite( z ) ) break;

        Real fz = fun->eval( z );
        ++m_fun_evaluation_count;

        // Check convergence at z
        if ( abs( fz ) < m_tolerance )
        {
          m_converged = true;
          return z;
        }

        Real s;
        if ( !safe_divide( s, fz, fy ) ) break;

        Real w = z - W( t, s ) * ( fz / df );
        if ( !Utils::is_finite( w ) ) break;

        Real fw = fun->eval( w );
        ++m_fun_evaluation_count;

        // Check convergence at w
        if ( abs( fw ) < m_tolerance )
        {
          m_converged = true;
          return w;
        }

        Real u;
        if ( !safe_divide( u, fw, fz ) ) break;

        Real weight = H( t, s, u );
        x           = w - weight * ( fw / df );

        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    /*!
     * \brief Optimal 32nd-order method (Varona 2022)
     * \param x_guess Initial guess
     * \param fun Pointer to function object
     * \return Approximate root or NaN if failed
     */
    Real solve_Order32( Real x_guess, Zeros_base_fun<Real> * fun )
    {
      if ( m_trace ) m_trace_values.clear();
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_converged            = false;

      Real x{ x_guess };
      while ( true )
      {
        Real f = fun->eval( x );
        ++m_fun_evaluation_count;

        if ( check_convergence( f, x ) ) return x;

        ++m_iteration_count;
        if ( check_limits() ) break;

        Real df = fun->eval_D( x );
        ++m_fun_evaluation_count;

        // Check for zero derivative
        if ( abs( df ) < machine_eps<Real>() )
        {
          m_converged = false;
          return FAILURE_VALUE;
        }

        Real y = x - ( f / df );
        if ( !Utils::is_finite( y ) ) break;

        Real fy = fun->eval( y );
        ++m_fun_evaluation_count;

        // Check convergence at y
        if ( abs( fy ) < m_tolerance )
        {
          m_converged = true;
          return y;
        }

        // Check for zero f to avoid division by zero
        if ( abs( f ) < machine_eps<Real>() ) break;

        Real t = fy / f;
        if ( !Utils::is_finite( t ) ) break;

        Real z = y - Q( t ) * ( fy / df );
        if ( !Utils::is_finite( z ) ) break;

        Real fz = fun->eval( z );
        ++m_fun_evaluation_count;

        // Check convergence at z
        if ( abs( fz ) < m_tolerance )
        {
          m_converged = true;
          return z;
        }

        // Check for zero fy to avoid division by zero
        if ( abs( fy ) < machine_eps<Real>() ) break;

        Real s = fz / fy;
        if ( !Utils::is_finite( s ) ) break;

        Real w = z - W( t, s ) * ( fz / df );
        if ( !Utils::is_finite( w ) ) break;

        Real fw = fun->eval( w );
        ++m_fun_evaluation_count;

        // Check convergence at w
        if ( abs( fw ) < m_tolerance )
        {
          m_converged = true;
          return w;
        }

        // Check for zero fz to avoid division by zero
        if ( abs( fz ) < machine_eps<Real>() ) break;

        Real u = fw / fz;
        if ( !Utils::is_finite( u ) ) break;

        Real h = w - H( t, s, u ) * ( fw / df );
        if ( !Utils::is_finite( h ) ) break;

        Real fh = fun->eval( h );
        ++m_fun_evaluation_count;

        // Check convergence at h
        if ( abs( fh ) < m_tolerance )
        {
          m_converged = true;
          return h;
        }

        // Check for zero fw to avoid division by zero
        if ( abs( fw ) < machine_eps<Real>() ) break;

        Real v = fh / fw;
        if ( !Utils::is_finite( v ) ) break;

        Real weight = J( t, s, u, v );
        x           = h - weight * ( fh / df );

        if ( !Utils::is_finite( x ) ) break;
      }

      m_converged = false;
      return FAILURE_VALUE;
    }

    //! Get number of iterations used in last computation
    Integer used_iter() const { return m_iteration_count; }

    //! Get number of function evaluations used in last computation
    Integer num_fun_eval() const { return m_fun_evaluation_count; }

    //! Get current tolerance
    Real tolerance() const { return m_tolerance; }

    //! Check if last computation converged
    bool converged() const { return m_converged; }

    //! Get failure value constant
    static constexpr Real failure_value() { return FAILURE_VALUE; }
  };

}  // namespace Utils

#endif

//
// eof: Utils_zeros.hh
//
