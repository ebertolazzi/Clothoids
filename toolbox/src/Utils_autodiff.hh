/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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
 * \file Utils_autodiff.hh
 * \brief Header file for automatic differentiation utilities with C++17 optimizations.
 *
 * This file provides an optimized automatic differentiation library supporting
 * forward-mode AD with up to 6th order derivatives. It includes:
 * - Unary and binary function definitions with derivative rules
 * - Macro-based derivative computation for functions with 1-6 arguments
 * - Optimized implementations using C++17 features (fold expressions, constexpr if)
 * - Integration with fmt library for formatting
 *
 * The library defines several mathematical functions (cbrt, erfc, log1p, etc.)
 * with their corresponding derivative rules, and provides macros to automatically
 * generate derivative functions for user-defined functions.
 *
 * \author Enrico Bertolazzi
 * \date 2013
 * \version Optimized for C++17 - COMPLETE VERSION
 */

//
// file: Utils_autodiff.hh
// Optimized for C++17 - COMPLETE VERSION
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_AUTODIFF_dot_HH
#define UTILS_AUTODIFF_dot_HH

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifdef _MSC_VER
#pragma warning( disable : 4127 )  ///< Disable warning for constant conditional expressions
#endif

#include "Utils/3rd/autodiff/forward/dual.hpp"
#include "Utils/3rd/autodiff/forward/real.hpp"

namespace fmt
{
  /// \brief fmt formatter specialization for autodiff::dual1st
  template <> struct formatter<autodiff::dual1st> : ostream_formatter
  {
  };
  /// \brief fmt formatter specialization for autodiff::dual2nd
  template <> struct formatter<autodiff::dual2nd> : ostream_formatter
  {
  };
  /// \brief fmt formatter specialization for autodiff::dual3rd
  template <> struct formatter<autodiff::dual3rd> : ostream_formatter
  {
  };
  /// \brief fmt formatter specialization for autodiff::dual4th
  template <> struct formatter<autodiff::dual4th> : ostream_formatter
  {
  };
}  // namespace fmt

#include <type_traits>

// ===========================================================================
// OPTIMIZED MACRO FOR UNARY FUNCTIONS (C++17)
// ===========================================================================

/**
 * \def UTILS_AUTODIFF_ADD_UNARY_FUNCTION( FUN )
 * \brief Macro to define a unary function for automatic differentiation.
 *
 * This macro generates:
 * - A tag structure for the function operation
 * - A function wrapper that returns a UnaryExpr
 * - An apply function that computes the function value and gradient
 *
 * \param FUN Name of the function to define (e.g., sin, cos, exp)
 */
#define UTILS_AUTODIFF_ADD_UNARY_FUNCTION( FUN )                                                                    \
  inline constexpr struct FUN##Op_tag                                                                               \
  {                                                                                                                 \
  } FUN##Op;                                                                                                        \
  template <typename R, Requires<isExpr<R>> = true> [[nodiscard]] AUTODIFF_DEVICE_FUNC constexpr auto FUN( R && r ) \
    -> UnaryExpr<FUN##Op_tag, R>                                                                                    \
  {                                                                                                                 \
    return { r };                                                                                                   \
  }                                                                                                                 \
  template <typename T, typename G>                                                                                 \
  AUTODIFF_DEVICE_FUNC constexpr void apply( Dual<T, G> & self, [[maybe_unused]] FUN##Op_tag )

namespace autodiff::detail
{

  using std::erfc;

  // ==========================================================================
  // OPTIMIZED: MaxN using fold expressions (C++17)
  // ==========================================================================

  /**
   * \brief Compile-time maximum of multiple size_t values using fold expressions.
   *
   * \tparam Ns Parameter pack of size_t values
   */
  template <size_t... Ns> struct MaxN
  {
    static constexpr size_t value = std::max( { Ns... } );  ///< Maximum value among Ns
  };

  /**
   * \brief Type trait to get the corresponding Dual type for a given type T.
   *
   * \tparam T Input type
   */
  template <typename T>
  using GetDual_t = HigherOrderDual<NumberTraits<T>::Order, typename NumberTraits<T>::NumericType>;

  /**
   * \brief Convert a value to its corresponding Dual type.
   *
   * \tparam T Type of the input value
   * \param x Value to convert
   * \return Corresponding Dual type with the same order as T
   */
  template <typename T> [[nodiscard]] constexpr auto to_dual( T const & x )
  {
    return GetDual_t<T>( x );
  }

  // ==========================================================================
  // OPTIMIZED: DualOrder using fold expressions (C++17)
  // ==========================================================================

  /**
   * \brief Compile-time maximum order among multiple Dual types.
   *
   * \tparam Ts Parameter pack of types
   */
  template <typename... Ts>
  struct DualOrder
  {
    static constexpr size_t value = []() {
      size_t max_order = 0;
      ((max_order = (NumberTraits<Ts>::Order > max_order) ? NumberTraits<Ts>::Order : max_order), ...);
      return max_order;
    }();
  };

  /*
  //       _          _
  //   ___| |__  _ __| |_
  //  / __| '_ \| '__| __|
  // | (__| |_) | |  | |_
  //  \___|_.__/|_|   \__|
  */

  /// \brief Cube root function for automatic differentiation
  inline constexpr struct CbrtOp_tag
  {
  } CbrtOp;

  template <typename R> using CbrtExpr = UnaryExpr<CbrtOp_tag, R>;  ///< Expression type for cube root

  /**
   * \brief Cube root function for automatic differentiation expressions.
   *
   * \tparam R Type of the expression
   * \param r Input expression
   * \return Cube root expression wrapper
   */
  template <typename R, Requires<isExpr<R>> = true> [[nodiscard]] AUTODIFF_DEVICE_FUNC constexpr auto cbrt( R && r )
    -> CbrtExpr<R>
  {
    return { r };
  }

  /**
   * \brief Apply cube root operation to a Dual number.
   *
   * Computes both the value and gradient of cbrt(x).
   * Special handling for negative values and zero.
   *
   * \tparam T Value type
   * \tparam G Gradient type
   * \param self Dual number to modify
   * \param CbrtOp_tag Tag for cube root operation
   */
  template <typename T, typename G>
  AUTODIFF_DEVICE_FUNC constexpr void apply( Dual<T, G> & self, [[maybe_unused]] CbrtOp_tag )
  {
    if ( self.val > 0 )
    {
      self.val = cbrt( self.val );
      self.grad *= 1 / ( 3 * self.val * self.val );
    }
    else if ( self.val < 0 )
    {
      self.val = -cbrt( -self.val );
      self.grad *= 1 / ( 3 * self.val * self.val );
    }
    else
    {
      self.val  = 0;
      self.grad = Utils::Inf<NumericType<T>>();
    }
  }

  /*
  //              __
  //    ___ _ __ / _| ___
  //   / _ \ '__| |_ / __|
  //  |  __/ |  |  _| (__
  //   \___|_|  |_|  \___|
  */

  /// \brief Complementary error function for automatic differentiation
  inline constexpr struct ErfcOp_tag
  {
  } ErfcOp;

  template <typename R> using ErfcExpr = UnaryExpr<ErfcOp_tag, R>;  ///< Expression type for erfc

  /**
   * \brief Complementary error function for automatic differentiation expressions.
   *
   * \tparam R Type of the expression
   * \param r Input expression
   * \return Erfc expression wrapper
   */
  template <typename R, Requires<isExpr<R>> = true> [[nodiscard]] AUTODIFF_DEVICE_FUNC constexpr auto erfc( R && r )
    -> ErfcExpr<R>
  {
    return { r };
  }

  /**
   * \brief Apply complementary error function to a Dual number.
   *
   * Computes erfc(x) and its gradient using the formula:
   * d/dx erfc(x) = -2/√π * exp(-x²)
   *
   * \tparam T Value type
   * \tparam G Gradient type
   * \param self Dual number to modify
   * \param ErfcOp_tag Tag for erfc operation
   */
  template <typename T, typename G>
  AUTODIFF_DEVICE_FUNC constexpr void apply( Dual<T, G> & self, [[maybe_unused]] ErfcOp_tag )
  {
    constexpr NumericType<T> sqrt_pi = 1.7724538509055160272981674833411451872554456638435;
    const T                  aux     = self.val;
    self.val                         = erfc( aux );
    self.grad *= -2.0 * exp( -aux * aux ) / sqrt_pi;
  }

  /*
  //                              _
  //   _ __ ___  _   _ _ __   __| |
  //  | '__/ _ \| | | | '_ \ / _` |
  //  | | | (_) | |_| | | | | (_| |
  //  |_|  \___/ \__,_|_| |_|\__,_|
  */

  /**
   * \brief Round function for Real numbers with automatic differentiation.
   *
   * Returns the nearest integer value, with derivatives set to zero
   * since round() is piecewise constant.
   *
   * \tparam N Order of the Real number
   * \tparam T Underlying numeric type
   * \param x Input Real number
   * \return Rounded Real number with zero derivatives
   */
  template <size_t N, typename T> AUTODIFF_DEVICE_FUNC constexpr auto round( Real<N, T> const & x )
  {
    Real<N, T> res;
    res[0] = std::round( x[0] );
    if constexpr ( N > 0 )
    {
      For<1, N + 1>( [&]( auto i ) constexpr { res[i] = 0; } );
    }
    return res;
  }

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( round )
  {
    using std::round;
    self.val  = round( val( self.val ) );
    self.grad = 0;  ///< Derivative is zero for round function
  }

  /**
   * \brief Round function for floating-point types.
   *
   * \tparam T Floating-point type
   * \param x Input value
   * \return Rounded value
   */
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type round( T const & x )
  {
    return round( Real<0, T>{ x } )[0];
  }

  /*
  //    __ _
  //   / _| | ___   ___  _ __
  //  | |_| |/ _ \ / _ \| '__|
  //  |  _| | (_) | (_) | |
  //  |_| |_|\___/ \___/|_|
  */

  /**
   * \brief Floor function for Real numbers with automatic differentiation.
   *
   * Returns the largest integer not greater than x, with derivatives set to zero.
   *
   * \tparam N Order of the Real number
   * \tparam T Underlying numeric type
   * \param x Input Real number
   * \return Floor value with zero derivatives
   */
  template <size_t N, typename T> AUTODIFF_DEVICE_FUNC constexpr auto floor( Real<N, T> const & x )
  {
    Real<N, T> res;
    res[0] = std::floor( x[0] );
    if constexpr ( N > 0 )
    {
      For<1, N + 1>( [&]( auto i ) constexpr { res[i] = 0; } );
    }
    return res;
  }

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( floor )
  {
    self.val  = std::floor( val( self.val ) );
    self.grad = 0;  ///< Derivative is zero for floor function
  }

  /**
   * \brief Floor function for floating-point types.
   *
   * \tparam T Floating-point type
   * \param x Input value
   * \return Floor value
   */
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type floor( T const & x )
  {
    return floor( Real<0, T>{ x } )[0];
  }

  /*
  //            _ _
  //    ___ ___(_) |
  //   / __/ _ \ | |
  //  | (_|  __/ | |
  //   \___\___|_|_|
  */

  /**
   * \brief Ceil function for Real numbers with automatic differentiation.
   *
   * Returns the smallest integer not less than x, with derivatives set to zero.
   *
   * \tparam N Order of the Real number
   * \tparam T Underlying numeric type
   * \param x Input Real number
   * \return Ceil value with zero derivatives
   */
  template <size_t N, typename T> AUTODIFF_DEVICE_FUNC constexpr auto ceil( Real<N, T> const & x )
  {
    Real<N, T> res;
    res[0] = std::ceil( x[0] );
    if constexpr ( N > 0 )
    {
      For<1, N + 1>( [&]( auto i ) constexpr { res[i] = 0; } );
    }
    return res;
  }

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( ceil )
  {
    using std::ceil;
    self.val  = ceil( val( self.val ) );
    self.grad = 0;  ///< Derivative is zero for ceil function
  }

  /**
   * \brief Ceil function for floating-point types.
   *
   * \tparam T Floating-point type
   * \param x Input value
   * \return Ceil value
   */
  template <typename T> constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type ceil( T const & x )
  {
    return ceil( Real<0, T>{ x } )[0];
  }

  /*
  //   _             _
  //  | | ___   __ _/ |_ __
  //  | |/ _ \ / _` | | '_ \
  //  | | (_) | (_| | | |_) |
  //  |_|\___/ \__, |_| .__/
  //           |___/  |_|
  */

  /**
   * \brief log1p function for Real numbers with automatic differentiation.
   *
   * Computes log(1+x) accurately for small x.
   * Implements higher-order derivatives using Faà di Bruno's formula.
   *
   * \tparam N Order of the Real number
   * \tparam T Underlying numeric type
   * \param x Input Real number
   * \return log(1+x) with derivatives up to order N
   * \throws Assertion error if x = -1 (undefined)
   */
  template <size_t N, typename T> AUTODIFF_DEVICE_FUNC constexpr auto log1p( Real<N, T> const & x )
  {
    assert( x[0] != -1 && "autodiff::log(1+x) has undefined value and derivatives when x = -1" );
    Real<N, T> log1px;
    log1px[0]     = std::log1p( x[0] );
    T one_plus_x0 = T( 1 ) + x[0];
    For<1, N + 1>(
      [&]( auto i ) constexpr
      {
        log1px[i] = x[i] - Sum<1, i>(
                             [&]( auto j ) constexpr
                             {
                               constexpr auto c = BinomialCoefficient<i.index - 1, j.index - 1>;
                               return c * x[i - j] * log1px[j];
                             } );
        log1px[i] /= one_plus_x0;
      } );
    return log1px;
  }

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( log1p )
  {
    using std::log1p;
    const T aux = One<T>() / ( One<T>() + self.val );
    self.val    = log1p( self.val );
    self.grad *= aux;  ///< d/dx log(1+x) = 1/(1+x)
  }

  /**
   * \brief log1p function for floating-point types.
   *
   * \tparam T Floating-point type
   * \param x Input value
   * \return log(1+x)
   */
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type log1p( T const & x )
  {
    return log1p( Real<0, T>{ x } )[0];
  }

  /*
  //         _              _
  //    __ _| |_ __ _ _ __ | |__
  //   / _` | __/ _` | '_ \| '_ \
  //  | (_| | || (_| | | | | | | |
  //   \__,_|\__\__,_|_| |_|_| |_|
  */

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( atanh )
  {
    using std::atanh;
    self.grad *= One<T>() / ( One<T>() - self.val * self.val );  ///< d/dx atanh(x) = 1/(1-x²)
    self.val = atanh( self.val );
  }

  /*
  //             _       _
  //    __ _ ___(_)_ __ | |__
  //   / _` / __| | '_ \| '_ \
  //  | (_| \__ \ | | | | | | |
  //   \__,_|___/_|_| |_|_| |_|
  */

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( asinh )
  {
    using std::asinh;
    self.grad *= One<T>() / sqrt( One<T>() + self.val * self.val );  ///< d/dx asinh(x) = 1/√(1+x²)
    self.val = asinh( self.val );
  }

  /*
  //                       _
  //    __ _  ___ ___  ___| |__
  //   / _` |/ __/ _ \/ __| '_ \
  //  | (_| | (_| (_) \__ \ | | |
  //   \__,_|\___\___/|___/_| |_|
  */

  UTILS_AUTODIFF_ADD_UNARY_FUNCTION( acosh )
  {
    using std::acosh;
    self.grad *= One<T>() / sqrt( self.val * self.val - One<T>() );  ///< d/dx acosh(x) = 1/√(x²-1)
    self.val = asinh( self.val );
  }

  /*
  //
  //   _ __   _____      _____ _ __
  //  | '_ \ / _ \ \ /\ / / _ \ '__|
  //  | |_) | (_) \ V  V /  __/ |
  //  | .__/ \___/ \_/\_/ \___|_|
  //  |_|
  */

  /// \brief Compute x²
  template <typename T> inline auto power2( T const & a )
  {
    return a * a;
  }
  /// \brief Compute x³
  template <typename T> inline auto power3( T const & a )
  {
    return a * a * a;
  }
  /// \brief Compute x⁴
  template <typename T> inline auto power4( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2;
  }
  /// \brief Compute x⁵
  template <typename T> inline auto power5( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a;
  }
  /// \brief Compute x⁶
  template <typename T> inline auto power6( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a2;
  }
  /// \brief Compute x⁷
  template <typename T> inline auto power7( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a2 * a;
  }
  /// \brief Compute x⁸
  template <typename T> inline auto power8( T const & a )
  {
    auto a2{ a * a };
    auto a4{ a2 * a2 };
    return a4 * a4;
  }

  /// \brief Compute 1/x²
  template <typename T> inline auto rpower2( T const & a )
  {
    return 1 / ( a * a );
  }
  /// \brief Compute 1/x³
  template <typename T> inline auto rpower3( T const & a )
  {
    return 1 / ( a * a * a );
  }
  /// \brief Compute 1/x⁴
  template <typename T> inline auto rpower4( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 );
  }
  /// \brief Compute 1/x⁵
  template <typename T> inline auto rpower5( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a );
  }
  /// \brief Compute 1/x⁶
  template <typename T> inline auto rpower6( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a2 );
  }
  /// \brief Compute 1/x⁷
  template <typename T> inline auto rpower7( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a2 * a );
  }
  /// \brief Compute 1/x⁸
  template <typename T> inline auto rpower8( T const & a )
  {
    auto a2{ a * a };
    auto a4{ a2 * a2 };
    return 1 / ( a4 * a4 );
  }

}  // namespace autodiff::detail

namespace Utils
{
  /// \brief Bring acosh function into Utils namespace
  using autodiff::detail::acosh;
  /// \brief Bring asinh function into Utils namespace
  using autodiff::detail::asinh;
  /// \brief Bring atanh function into Utils namespace
  using autodiff::detail::atanh;
  /// \brief Bring ceil function into Utils namespace
  using autodiff::detail::ceil;
  /// \brief Bring erfc function into Utils namespace
  using autodiff::detail::erfc;
  /// \brief Bring floor function into Utils namespace
  using autodiff::detail::floor;
  /// \brief Bring log1p function into Utils namespace
  using autodiff::detail::log1p;
  /// \brief Bring power2 function into Utils namespace
  using autodiff::detail::power2;
  /// \brief Bring power3 function into Utils namespace
  using autodiff::detail::power3;
  /// \brief Bring power4 function into Utils namespace
  using autodiff::detail::power4;
  /// \brief Bring power5 function into Utils namespace
  using autodiff::detail::power5;
  /// \brief Bring power6 function into Utils namespace
  using autodiff::detail::power6;
  /// \brief Bring power7 function into Utils namespace
  using autodiff::detail::power7;
  /// \brief Bring power8 function into Utils namespace
  using autodiff::detail::power8;
  /// \brief Bring round function into Utils namespace
  using autodiff::detail::round;
  /// \brief Bring rpower2 function into Utils namespace
  using autodiff::detail::rpower2;
  /// \brief Bring rpower3 function into Utils namespace
  using autodiff::detail::rpower3;
  /// \brief Bring rpower4 function into Utils namespace
  using autodiff::detail::rpower4;
  /// \brief Bring rpower5 function into Utils namespace
  using autodiff::detail::rpower5;
  /// \brief Bring rpower6 function into Utils namespace
  using autodiff::detail::rpower6;
  /// \brief Bring rpower7 function into Utils namespace
  using autodiff::detail::rpower7;
  /// \brief Bring rpower8 function into Utils namespace
  using autodiff::detail::rpower8;


  using autodiff::detail::DualOrder;


}  // namespace Utils

// ===========================================================================
// DERIVATIVE MACROS - COMPLETE FOR ALL ARGUMENTS (1-6)
// ===========================================================================

// ---------------------------------------------------------------------------
// 1-ARG DERIVATIVES
// ---------------------------------------------------------------------------
/**
 * \def UTILS_AUTODIFF_DERIV_1ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate first and second derivatives for 1-argument functions.
 *
 * Generates:
 * - PREFIXD(x): First derivative
 * - PREFIXDD(x): Second derivative
 *
 * \param INLINE Inline specifier (inline, static, etc.)
 * \param CLASS Class qualifier (empty or class name)
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_1ARG( INLINE, CLASS, PREFIX, FUN, CONST ) \
  INLINE real_type CLASS PREFIX##D( real_type const x ) CONST          \
  {                                                                    \
    autodiff::dual1st X{ x };                                          \
    X.grad = 1;                                                        \
    autodiff::dual1st res{ FUN( X ) };                                 \
    return res.grad;                                                   \
  }                                                                    \
                                                                       \
  INLINE real_type CLASS PREFIX##DD( real_type const x ) CONST         \
  {                                                                    \
    autodiff::dual2nd X{ x };                                          \
    X.val.grad = 1;                                                    \
    X.grad.val = 1;                                                    \
    autodiff::dual2nd res{ FUN( X ) };                                 \
    return res.grad.grad;                                              \
  }

// ---------------------------------------------------------------------------
// 2-ARG DERIVATIVES
// ---------------------------------------------------------------------------
/**
 * \def UTILS_AUTODIFF_DERIV_2ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate derivatives for 2-argument functions.
 *
 * Generates:
 * - First derivatives: PREFIXD_1, PREFIXD_2
 * - Second derivatives: PREFIXD_1_1, PREFIXD_1_2, PREFIXD_2_2
 *
 * \param INLINE Inline specifier
 * \param CLASS Class qualifier
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_2ARG( INLINE, CLASS, PREFIX, FUN, CONST )               \
  INLINE real_type CLASS PREFIX##D_1( real_type const x, real_type const y ) CONST   \
  {                                                                                  \
    autodiff::dual1st X{ x };                                                        \
    X.grad = 1;                                                                      \
    autodiff::dual1st res{ FUN( X, y ) };                                            \
    return res.grad;                                                                 \
  }                                                                                  \
                                                                                     \
  INLINE real_type CLASS PREFIX##D_2( real_type const x, real_type const y ) CONST   \
  {                                                                                  \
    autodiff::dual1st Y{ y };                                                        \
    Y.grad = 1;                                                                      \
    autodiff::dual1st res{ FUN( x, Y ) };                                            \
    return res.grad;                                                                 \
  }                                                                                  \
                                                                                     \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x, real_type const y ) CONST \
  {                                                                                  \
    autodiff::dual2nd X{ x };                                                        \
    X.val.grad = 1;                                                                  \
    X.grad.val = 1;                                                                  \
    autodiff::dual2nd res{ FUN( X, y ) };                                            \
    return res.grad.grad;                                                            \
  }                                                                                  \
                                                                                     \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x, real_type const y ) CONST \
  {                                                                                  \
    autodiff::dual2nd X{ x }, Y{ y };                                                \
    X.val.grad = 1;                                                                  \
    Y.grad.val = 1;                                                                  \
    autodiff::dual2nd res{ FUN( X, Y ) };                                            \
    return res.grad.grad;                                                            \
  }                                                                                  \
                                                                                     \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x, real_type const y ) CONST \
  {                                                                                  \
    autodiff::dual2nd Y{ y };                                                        \
    Y.val.grad = 1;                                                                  \
    Y.grad.val = 1;                                                                  \
    autodiff::dual2nd res{ FUN( x, Y ) };                                            \
    return res.grad.grad;                                                            \
  }

// ---------------------------------------------------------------------------
// 3-ARG DERIVATIVES
// ---------------------------------------------------------------------------
/// \brief Parameter list for 3-argument functions
#define _UTILS_AUTODIFF_3ARG_PARAMS real_type const x, real_type const y, real_type const z

/**
 * \def UTILS_AUTODIFF_DERIV_3ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate derivatives for 3-argument functions.
 *
 * Generates:
 * - First derivatives: PREFIXD_1, PREFIXD_2, PREFIXD_3
 * - Second derivatives: All combinations PREFIXD_i_j
 *
 * \param INLINE Inline specifier
 * \param CLASS Class qualifier
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_3ARG( INLINE, CLASS, PREFIX, FUN, CONST )      \
  INLINE real_type CLASS PREFIX##D_1( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X{ x };                                               \
    X.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( X, y, z ) };                                \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st Y{ y };                                               \
    Y.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( x, Y, z ) };                                \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st Z{ z };                                               \
    Z.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( x, y, Z ) };                                \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_1( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x };                                               \
    X.val.grad = 1;                                                         \
    X.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, y, z ) };                                \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_2( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x }, Y{ y };                                       \
    X.val.grad = 1;                                                         \
    Y.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, Y, z ) };                                \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_3( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x }, Z{ z };                                       \
    X.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, y, Z ) };                                \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_2( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Y{ y };                                               \
    Y.val.grad = 1;                                                         \
    Y.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, Y, z ) };                                \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_3( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Y{ y }, Z{ z };                                       \
    Y.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, Y, Z ) };                                \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_3( _UTILS_AUTODIFF_3ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Z{ z };                                               \
    Z.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, y, Z ) };                                \
    return res.grad.grad;                                                   \
  }

// ---------------------------------------------------------------------------
// 4-ARG DERIVATIVES
// ---------------------------------------------------------------------------
/// \brief Parameter list for 4-argument functions
#define _UTILS_AUTODIFF_4ARG_PARAMS real_type const x, real_type const y, real_type const z, real_type const w

/**
 * \def UTILS_AUTODIFF_DERIV_4ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate derivatives for 4-argument functions.
 *
 * Generates all first and second partial derivatives for functions with 4 arguments.
 *
 * \param INLINE Inline specifier
 * \param CLASS Class qualifier
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_4ARG( INLINE, CLASS, PREFIX, FUN, CONST )      \
  /* First derivatives */                                                   \
  INLINE real_type CLASS PREFIX##D_1( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X{ x };                                               \
    X.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( X, y, z, w ) };                             \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st Y{ y };                                               \
    Y.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( x, Y, z, w ) };                             \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st Z{ z };                                               \
    Z.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( x, y, Z, w ) };                             \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st W{ w };                                               \
    W.grad = 1;                                                             \
    autodiff::dual1st res{ FUN( x, y, z, W ) };                             \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  /* Second derivatives */                                                  \
  INLINE real_type CLASS PREFIX##D_1_1( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x };                                               \
    X.val.grad = 1;                                                         \
    X.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, y, z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_2( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x }, Y{ y };                                       \
    X.val.grad = 1;                                                         \
    Y.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, Y, z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_3( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x }, Z{ z };                                       \
    X.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, y, Z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_4( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X{ x }, W{ w };                                       \
    X.val.grad = 1;                                                         \
    W.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( X, y, z, W ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_2( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Y{ y };                                               \
    Y.val.grad = 1;                                                         \
    Y.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, Y, z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_3( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Y{ y }, Z{ z };                                       \
    Y.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, Y, Z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_4( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Y{ y }, W{ w };                                       \
    Y.val.grad = 1;                                                         \
    W.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, Y, z, W ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_3( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Z{ z };                                               \
    Z.val.grad = 1;                                                         \
    Z.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, y, Z, w ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_4( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd Z{ z }, W{ w };                                       \
    Z.val.grad = 1;                                                         \
    W.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, y, Z, W ) };                             \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_4( _UTILS_AUTODIFF_4ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd W{ w };                                               \
    W.val.grad = 1;                                                         \
    W.grad.val = 1;                                                         \
    autodiff::dual2nd res{ FUN( x, y, z, W ) };                             \
    return res.grad.grad;                                                   \
  }

// ---------------------------------------------------------------------------
// 5-ARG DERIVATIVES (COMPLETED - ALL SECOND DERIVATIVES)
// ---------------------------------------------------------------------------

/// \brief Parameter list for 5-argument functions
#define _UTILS_AUTODIFF_5ARG_PARAMS \
  real_type const x1, real_type const x2, real_type const x3, real_type const x4, real_type const x5

/**
 * \def UTILS_AUTODIFF_DERIV_5ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate derivatives for 5-argument functions.
 *
 * Generates all first and second partial derivatives for functions with 5 arguments.
 *
 * \param INLINE Inline specifier
 * \param CLASS Class qualifier
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_5ARG( INLINE, CLASS, PREFIX, FUN, CONST )      \
  /* First derivatives */                                                   \
  INLINE real_type CLASS PREFIX##D_1( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X1{ x1 };                                             \
    X1.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( X1, x2, x3, x4, x5 ) };                     \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X2{ x2 };                                             \
    X2.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, X2, x3, x4, x5 ) };                     \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X3{ x3 };                                             \
    X3.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, X3, x4, x5 ) };                     \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X4{ x4 };                                             \
    X4.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, x3, X4, x5 ) };                     \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X5{ x5 };                                             \
    X5.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, X5 ) };                     \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  /* Second derivatives - all combinations */                               \
  INLINE real_type CLASS PREFIX##D_1_1( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 };                                             \
    X1.val.grad = 1;                                                        \
    X1.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_2( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X2{ x2 };                                   \
    X1.val.grad = 1;                                                        \
    X2.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, X2, x3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_3( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X3{ x3 };                                   \
    X1.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, X3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_4( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X4{ x4 };                                   \
    X1.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, X4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X5{ x5 };                                   \
    X1.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, X5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_2( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 };                                             \
    X2.val.grad = 1;                                                        \
    X2.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_3( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X3{ x3 };                                   \
    X2.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, X3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_4( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X4{ x4 };                                   \
    X2.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, X4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X5{ x5 };                                   \
    X2.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, X5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_3( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 };                                             \
    X3.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_4( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 }, X4{ x4 };                                   \
    X3.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, X4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 }, X5{ x5 };                                   \
    X3.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, X5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_4( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X4{ x4 };                                             \
    X4.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X4{ x4 }, X5{ x5 };                                   \
    X4.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, X5 ) };                     \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_5_5( _UTILS_AUTODIFF_5ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X5{ x5 };                                             \
    X5.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5 ) };                     \
    return res.grad.grad;                                                   \
  }

// ---------------------------------------------------------------------------
// 6-ARG DERIVATIVES (COMPLETED - ALL SECOND DERIVATIVES)
// ---------------------------------------------------------------------------

/// \brief Parameter list for 6-argument functions
#define _UTILS_AUTODIFF_6ARG_PARAMS \
  real_type const x1, real_type const x2, real_type const x3, real_type const x4, real_type const x5, real_type const x6

/**
 * \def UTILS_AUTODIFF_DERIV_6ARG( INLINE, CLASS, PREFIX, FUN, CONST )
 * \brief Macro to generate derivatives for 6-argument functions.
 *
 * Generates all first and second partial derivatives for functions with 6 arguments.
 * This includes 6 first derivatives and 21 second derivatives (6 diagonal + 15 mixed).
 *
 * \param INLINE Inline specifier
 * \param CLASS Class qualifier
 * \param PREFIX Function name prefix
 * \param FUN Function to differentiate
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_DERIV_6ARG( INLINE, CLASS, PREFIX, FUN, CONST )      \
  /* First derivatives */                                                   \
  INLINE real_type CLASS PREFIX##D_1( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X1{ x1 };                                             \
    X1.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( X1, x2, x3, x4, x5, x6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X2{ x2 };                                             \
    X2.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, X2, x3, x4, x5, x6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X3{ x3 };                                             \
    X3.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, X3, x4, x5, x6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X4{ x4 };                                             \
    X4.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, x3, X4, x5, x6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X5{ x5 };                                             \
    X5.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, X5, x6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST   \
  {                                                                         \
    autodiff::dual1st X6{ x6 };                                             \
    X6.grad = 1;                                                            \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, x5, X6 ) };                 \
    return res.grad;                                                        \
  }                                                                         \
                                                                            \
  /* Second derivatives - all 21 combinations */                            \
  INLINE real_type CLASS PREFIX##D_1_1( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 };                                             \
    X1.val.grad = 1;                                                        \
    X1.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_2( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X2{ x2 };                                   \
    X1.val.grad = 1;                                                        \
    X2.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, X2, x3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_3( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X3{ x3 };                                   \
    X1.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, X3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_4( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X4{ x4 };                                   \
    X1.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, X4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X5{ x5 };                                   \
    X1.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, X5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_1_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X1{ x1 }, X6{ x6 };                                   \
    X1.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_2( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 };                                             \
    X2.val.grad = 1;                                                        \
    X2.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_3( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X3{ x3 };                                   \
    X2.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, X3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_4( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X4{ x4 };                                   \
    X2.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, X4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X5{ x5 };                                   \
    X2.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, X5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_2_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X2{ x2 }, X6{ x6 };                                   \
    X2.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_3( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 };                                             \
    X3.val.grad = 1;                                                        \
    X3.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_4( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 }, X4{ x4 };                                   \
    X3.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, X4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 }, X5{ x5 };                                   \
    X3.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, X5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_3_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X3{ x3 }, X6{ x6 };                                   \
    X3.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_4( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X4{ x4 };                                             \
    X4.val.grad = 1;                                                        \
    X4.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X4{ x4 }, X5{ x5 };                                   \
    X4.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, X5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_4_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X4{ x4 }, X6{ x6 };                                   \
    X4.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_5_5( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X5{ x5 };                                             \
    X5.val.grad = 1;                                                        \
    X5.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5, x6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_5_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X5{ x5 }, X6{ x6 };                                   \
    X5.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }                                                                         \
                                                                            \
  INLINE real_type CLASS PREFIX##D_6_6( _UTILS_AUTODIFF_6ARG_PARAMS ) CONST \
  {                                                                         \
    autodiff::dual2nd X6{ x6 };                                             \
    X6.val.grad = 1;                                                        \
    X6.grad.val = 1;                                                        \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, x5, X6 ) };                 \
    return res.grad.grad;                                                   \
  }

// ===========================================================================
// DECLARATION MACROS FOR FUNCTION SIGNATURES
// ===========================================================================

/// \brief Parameter list for 1-argument function declarations
#define UTILS_AUTODIFF_PARAMS_1 real_type const x1
/// \brief Parameter list for 2-argument function declarations
#define UTILS_AUTODIFF_PARAMS_2 UTILS_AUTODIFF_PARAMS_1, real_type const x2
/// \brief Parameter list for 3-argument function declarations
#define UTILS_AUTODIFF_PARAMS_3 UTILS_AUTODIFF_PARAMS_2, real_type const x3
/// \brief Parameter list for 4-argument function declarations
#define UTILS_AUTODIFF_PARAMS_4 UTILS_AUTODIFF_PARAMS_3, real_type const x4
/// \brief Parameter list for 5-argument function declarations
#define UTILS_AUTODIFF_PARAMS_5 UTILS_AUTODIFF_PARAMS_4, real_type const x5
/// \brief Parameter list for 6-argument function declarations
#define UTILS_AUTODIFF_PARAMS_6 UTILS_AUTODIFF_PARAMS_5, real_type const x6

/**
 * \def UTILS_AUTODIFF_DECLARE_DERIV(PREFIX, SUFFIX, PARAMS, CONST_QUAL)
 * \brief Helper macro to declare a single derivative function.
 *
 * \param PREFIX Function name prefix
 * \param SUFFIX Derivative suffix (e.g., D_1, D_1_2)
 * \param PARAMS Parameter list
 * \param CONST_QUAL Const qualifier
 */
#define UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, SUFFIX, PARAMS, CONST_QUAL ) \
  real_type PREFIX##SUFFIX( PARAMS ) CONST_QUAL;

/**
 * \def UTILS_AUTODIFF_FUN_1_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 1-variable function.
 *
 * Declares: PREFIXD_1, PREFIXD_1_1
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_1_VARS_DECL( PREFIX, CONST )                       \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D, UTILS_AUTODIFF_PARAMS_1, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, DD, UTILS_AUTODIFF_PARAMS_1, CONST )

/**
 * \def UTILS_AUTODIFF_FUN_2_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 2-variable function.
 *
 * Declares all first and second partial derivatives.
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_2_VARS_DECL( PREFIX, CONST )                         \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1, UTILS_AUTODIFF_PARAMS_2, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2, UTILS_AUTODIFF_PARAMS_2, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_1, UTILS_AUTODIFF_PARAMS_2, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_2, UTILS_AUTODIFF_PARAMS_2, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_2, UTILS_AUTODIFF_PARAMS_2, CONST )

/**
 * \def UTILS_AUTODIFF_FUN_3_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 3-variable function.
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_3_VARS_DECL( PREFIX, CONST )                         \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1, UTILS_AUTODIFF_PARAMS_3, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2, UTILS_AUTODIFF_PARAMS_3, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3, UTILS_AUTODIFF_PARAMS_3, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_1, UTILS_AUTODIFF_PARAMS_3, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_2, UTILS_AUTODIFF_PARAMS_3, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_3, UTILS_AUTODIFF_PARAMS_3, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_2, UTILS_AUTODIFF_PARAMS_3, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_3, UTILS_AUTODIFF_PARAMS_3, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_3, UTILS_AUTODIFF_PARAMS_3, CONST )

/**
 * \def UTILS_AUTODIFF_FUN_4_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 4-variable function.
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_4_VARS_DECL( PREFIX, CONST )                         \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1, UTILS_AUTODIFF_PARAMS_4, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2, UTILS_AUTODIFF_PARAMS_4, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3, UTILS_AUTODIFF_PARAMS_4, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4, UTILS_AUTODIFF_PARAMS_4, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_1, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_2, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_3, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_4, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_2, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_3, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_4, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_3, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_4, UTILS_AUTODIFF_PARAMS_4, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_4, UTILS_AUTODIFF_PARAMS_4, CONST )

/**
 * \def UTILS_AUTODIFF_FUN_5_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 5-variable function.
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_5_VARS_DECL( PREFIX, CONST )                         \
  /* Declarations for all first derivatives */                                  \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1, UTILS_AUTODIFF_PARAMS_5, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2, UTILS_AUTODIFF_PARAMS_5, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3, UTILS_AUTODIFF_PARAMS_5, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4, UTILS_AUTODIFF_PARAMS_5, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_5, UTILS_AUTODIFF_PARAMS_5, CONST )   \
  /* Declarations for all second derivatives */                                 \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_1, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_2, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_3, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_4, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_5, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_2, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_3, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_4, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_5, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_3, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_4, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_5, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_4, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_5, UTILS_AUTODIFF_PARAMS_5, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_5_5, UTILS_AUTODIFF_PARAMS_5, CONST )

/**
 * \def UTILS_AUTODIFF_FUN_6_VARS_DECL( PREFIX, CONST )
 * \brief Declare all derivatives for a 6-variable function.
 *
 * Declares all 6 first derivatives and 21 second derivatives.
 *
 * \param PREFIX Function name prefix
 * \param CONST Const qualifier
 */
#define UTILS_AUTODIFF_FUN_6_VARS_DECL( PREFIX, CONST )                         \
  /* First derivatives */                                                       \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_5, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_6, UTILS_AUTODIFF_PARAMS_6, CONST )   \
  /* Second derivatives - all combinations */                                   \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_1, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_2, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_3, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_4, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_5, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_1_6, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_2, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_3, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_4, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_5, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_2_6, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_3, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_4, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_5, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_3_6, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_4, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_5, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_4_6, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_5_5, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_5_6, UTILS_AUTODIFF_PARAMS_6, CONST ) \
  UTILS_AUTODIFF_DECLARE_DERIV( PREFIX, D_6_6, UTILS_AUTODIFF_PARAMS_6, CONST )

#define UTILS_AUTODIFF_DERIV_1ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_1_VARS_DECL( PREFIX, CONST )
#define UTILS_AUTODIFF_DERIV_2ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_2_VARS_DECL( PREFIX, CONST )
#define UTILS_AUTODIFF_DERIV_3ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_3_VARS_DECL( PREFIX, CONST )
#define UTILS_AUTODIFF_DERIV_4ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_4_VARS_DECL( PREFIX, CONST )
#define UTILS_AUTODIFF_DERIV_5ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_5_VARS_DECL( PREFIX, CONST )
#define UTILS_AUTODIFF_DERIV_6ARG_DEF( PREFIX, CONST ) UTILS_AUTODIFF_FUN_6_VARS_DECL( PREFIX, CONST )

#endif  // UTILS_AUTODIFF_dot_HH

#endif  // DOXYGEN_SHOULD_SKIP_THIS

//
// eof: Utils_autodiff.hh
//
