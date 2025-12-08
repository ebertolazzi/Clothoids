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

//
// file: Utils_autodiff.hh
//
#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_AUTODIFF_dot_HH
#define UTILS_AUTODIFF_dot_HH

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifdef _MSC_VER
#pragma warning( disable : 4127 )
#endif

#include "Utils/3rd/autodiff/forward/dual.hpp"
#include "Utils/3rd/autodiff/forward/real.hpp"

namespace fmt
{
  template <>
  struct formatter<autodiff::dual1st> : ostream_formatter
  {
  };
  template <>
  struct formatter<autodiff::dual2nd> : ostream_formatter
  {
  };
  template <>
  struct formatter<autodiff::dual3rd> : ostream_formatter
  {
  };
  template <>
  struct formatter<autodiff::dual4th> : ostream_formatter
  {
  };
}  // namespace fmt

#include <type_traits>

#define UTILS_AUTODIFF_ADD_UNARY_FUNCTION( FUN )                                                                       \
  struct FUN##Op                                                                                                       \
  {                                                                                                                    \
  };                                                                                                                   \
  template <typename R, Requires<isExpr<R>> = true>                                                                    \
  AUTODIFF_DEVICE_FUNC constexpr auto FUN( R && r ) -> UnaryExpr<FUN##Op, R>                                           \
  {                                                                                                                    \
    return { r };                                                                                                      \
  }                                                                                                                    \
  template <typename T, typename G>                                                                                    \
  AUTODIFF_DEVICE_FUNC constexpr void apply( Dual<T, G> & self, FUN##Op )

namespace autodiff::detail
{

  using std::erfc;

  template <size_t A, size_t... Rest>
  struct MaxN
  {
    static constexpr size_t value = []()
    {
      if constexpr ( sizeof...( Rest ) == 0 ) { return A; }
      else
      {
        size_t max_rest = MaxN<Rest...>::value;
        return A > max_rest ? A : max_rest;
      }
    }();
  };

  template <typename T>
  using GetDual_t = HigherOrderDual<NumberTraits<T>::Order, typename NumberTraits<T>::NumericType>;

  template <typename T>
  constexpr auto
  to_dual( T const & x )
  {
    return GetDual_t<T>( x );
  }

  // Caso base: nessun tipo (valore 0 per default)
  template <typename... Ts>
  struct DualOrder
  {
    static constexpr size_t value = 0;
  };
  template <typename T, typename... Ts>
  struct DualOrder<T, Ts...>
  {
    static constexpr size_t value = []()
    {
      if constexpr ( sizeof...( Ts ) == 0 )
      {
        return NumberTraits<T>::Order;  // Caso singolo tipo
      }
      else
      {
        constexpr size_t current_order = NumberTraits<T>::Order;
        constexpr size_t rest_order    = DualOrder<Ts...>::value;
        return ( current_order > rest_order ) ? current_order : rest_order;
      }
    }();
  };

  /*
  //       _          _
  //   ___| |__  _ __| |_
  //  / __| '_ \| '__| __|
  // | (__| |_) | |  | |_
  //  \___|_.__/|_|   \__|
  */

  // missing code for cbrt
  struct CbrtOp
  {
  };  // CUBIC ROOT OPERATOR
  template <typename R>
  using CbrtExpr = UnaryExpr<CbrtOp, R>;
  template <typename R, Requires<isExpr<R>> = true>
  AUTODIFF_DEVICE_FUNC constexpr auto
  cbrt( R && r ) -> CbrtExpr<R>
  {
    return { r };
  }

  template <typename T, typename G>
  AUTODIFF_DEVICE_FUNC constexpr void
  apply( Dual<T, G> & self, CbrtOp )
  {
    self.val = cbrt( self.val );
    self.grad *= 1 / ( 3 * self.val * self.val );
  }

  /*
  //              __
  //    ___ _ __ / _| ___
  //   / _ \ '__| |_ / __|
  //  |  __/ |  |  _| (__
  //   \___|_|  |_|  \___|
  */

  // missing code for erfc
  struct ErfcOp
  {
  };  // ERROR FUNCTION OPERATOR
  template <typename R>
  using ErfcExpr = UnaryExpr<ErfcOp, R>;
  template <typename R, Requires<isExpr<R>> = true>
  AUTODIFF_DEVICE_FUNC constexpr auto
  erfc( R && r ) -> ErfcExpr<R>
  {
    return { r };
  }

  template <typename T, typename G>
  AUTODIFF_DEVICE_FUNC constexpr void
  apply( Dual<T, G> & self, ErfcOp )
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
  //
  */

  template <size_t N, typename T>
  AUTODIFF_DEVICE_FUNC constexpr auto
  round( Real<N, T> const & x )
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
    self.grad = 0;
  }

  // overload per tipi floating-point (double, float, ...)
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type
  round( T const & x )
  {
    return round( Real<0, T>{ x } )[0];
  }

  /*
  //    __ _
  //   / _| | ___   ___  _ __
  //  | |_| |/ _ \ / _ \| '__|
  //  |  _| | (_) | (_) | |
  //  |_| |_|\___/ \___/|_|
  //
  */

  template <size_t N, typename T>
  AUTODIFF_DEVICE_FUNC constexpr auto
  floor( Real<N, T> const & x )
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
    self.grad = 0;
  }

  // overload per tipi floating-point (double, float, ...)
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type
  floor( T const & x )
  {
    return floor( Real<0, T>{ x } )[0];
  }

  /*
  //            _ _
  //    ___ ___(_) |
  //   / __/ _ \ | |
  //  | (_|  __/ | |
  //   \___\___|_|_|
  //
  */

  template <size_t N, typename T>
  AUTODIFF_DEVICE_FUNC constexpr auto
  ceil( Real<N, T> const & x )
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
    self.grad = 0;
  }

  // overload per tipi floating-point (double, float, ...)
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type
  ceil( T const & x )
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

  template <size_t N, typename T>
  AUTODIFF_DEVICE_FUNC constexpr auto
  log1p( Real<N, T> const & x )
  {
    assert( x[0] != 1 && "autodiff::log(1+x) has undefined value and derivatives when x = -1" );
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
    const T aux = One<T>() / ( One<T>() + self.val );  // 1 / (1 + x)
    self.val    = log1p( self.val );                   // log(1 + x)
    self.grad *= aux;                                  // grad * 1/(1 + x)
  }

  // overload per tipi floating-point (double, float, ...)
  template <typename T>
  constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type
  log1p( T const & x )
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
    self.val = atanh( self.val );
    self.grad *= One<T>() / ( 1 - self.val * self.val );
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
    self.val = asinh( self.val );
    self.grad *= One<T>() / ( 1 + self.val * self.val );
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
    self.val = asinh( self.val );
    self.grad *= One<T>() / sqrt( self.val * self.val - 1 );
  }

  /*
  //
  //   _ __   _____      _____ _ __
  //  | '_ \ / _ \ \ /\ / / _ \ '__|
  //  | |_) | (_) \ V  V /  __/ |
  //  | .__/ \___/ \_/\_/ \___|_|
  //  |_|
  */

  template <typename T>
  inline auto
  power2( T const & a )
  {
    return a * a;
  }
  template <typename T>
  inline auto
  power3( T const & a )
  {
    return a * a * a;
  }
  template <typename T>
  inline auto
  power4( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2;
  }
  template <typename T>
  inline auto
  power5( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a;
  }
  template <typename T>
  inline auto
  power6( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a2;
  }
  template <typename T>
  inline auto
  power7( T const & a )
  {
    auto a2{ a * a };
    return a2 * a2 * a2 * a;
  }
  template <typename T>
  inline auto
  power8( T const & a )
  {
    auto a2{ a * a };
    auto a4{ a2 * a2 };
    return a4 * a4;
  }

  template <typename T>
  inline auto
  rpower2( T const & a )
  {
    return 1 / ( a * a );
  }
  template <typename T>
  inline auto
  rpower3( T const & a )
  {
    return 1 / ( a * a * a );
  }
  template <typename T>
  inline auto
  rpower4( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 );
  }
  template <typename T>
  inline auto
  rpower5( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a );
  }
  template <typename T>
  inline auto
  rpower6( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a2 );
  }
  template <typename T>
  inline auto
  rpower7( T const & a )
  {
    auto a2{ a * a };
    return 1 / ( a2 * a2 * a2 * a );
  }
  template <typename T>
  inline auto
  rpower8( T const & a )
  {
    auto a2{ a * a };
    auto a4{ a2 * a2 };
    return 1 / ( a4 * a4 );
  }
}  // namespace autodiff::detail

namespace Utils
{
  using autodiff::detail::acosh;
  using autodiff::detail::asinh;
  using autodiff::detail::atanh;
  using autodiff::detail::ceil;
  using autodiff::detail::erfc;
  using autodiff::detail::floor;
  using autodiff::detail::log1p;
  using autodiff::detail::power2;
  using autodiff::detail::power3;
  using autodiff::detail::power4;
  using autodiff::detail::power5;
  using autodiff::detail::power6;
  using autodiff::detail::power7;
  using autodiff::detail::power8;
  using autodiff::detail::round;
  using autodiff::detail::rpower2;
  using autodiff::detail::rpower3;
  using autodiff::detail::rpower4;
  using autodiff::detail::rpower5;
  using autodiff::detail::rpower6;
  using autodiff::detail::rpower7;
  using autodiff::detail::rpower8;
}  // namespace Utils

#define UTILS_AUTODIFF_DERIV_1ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D( real_type const x ) CONST                                                          \
  {                                                                                                                    \
    autodiff::dual1st X{ x };                                                                                          \
    X.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( X ) };                                                                                 \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##DD( real_type const x ) CONST                                                         \
  {                                                                                                                    \
    autodiff::dual2nd X{ x };                                                                                          \
    X.val.grad = 1;                                                                                                    \
    X.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X ) };                                                                                 \
    return res.grad.grad;                                                                                              \
  }


#define UTILS_AUTODIFF_DERIV_2ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D_1( real_type const x, real_type const y ) CONST                                     \
  {                                                                                                                    \
    autodiff::dual1st X{ x };                                                                                          \
    X.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( X, y ) };                                                                              \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2( real_type const x, real_type const y ) CONST                                     \
  {                                                                                                                    \
    autodiff::dual1st Y{ y };                                                                                          \
    Y.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, Y ) };                                                                              \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x, real_type const y ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual2nd X{ x };                                                                                          \
    X.val.grad = 1;                                                                                                    \
    X.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y ) };                                                                              \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x, real_type const y ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, Y{ y };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, Y ) };                                                                              \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x, real_type const y ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y };                                                                                          \
    Y.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y ) };                                                                              \
    return res.grad.grad;                                                                                              \
  }


#define UTILS_AUTODIFF_DERIV_3ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D_1( real_type const x, real_type const y, real_type const z ) CONST                  \
  {                                                                                                                    \
    autodiff::dual1st X{ x };                                                                                          \
    X.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( X, y, z ) };                                                                           \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2( real_type const x, real_type const y, real_type const z ) CONST                  \
  {                                                                                                                    \
    autodiff::dual1st Y{ y };                                                                                          \
    Y.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, Y, z ) };                                                                           \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3( real_type const x, real_type const y, real_type const z ) CONST                  \
  {                                                                                                                    \
    autodiff::dual1st Z{ z };                                                                                          \
    Z.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, y, Z ) };                                                                           \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd X{ x };                                                                                          \
    X.val.grad = 1;                                                                                                    \
    X.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y, z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, Y{ y };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, Y, z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_3( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, Z{ z };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y, Z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y };                                                                                          \
    Y.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y, z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_3( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y }, Z{ z };                                                                                  \
    Y.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y, Z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_3( real_type const x, real_type const y, real_type const z ) CONST                \
  {                                                                                                                    \
    autodiff::dual2nd Z{ z };                                                                                          \
    Z.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, y, Z ) };                                                                           \
    return res.grad.grad;                                                                                              \
  }


#define UTILS_AUTODIFF_DERIV_4ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D_1( real_type const x, real_type const y, real_type const z, real_type const w )     \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual1st X{ x };                                                                                          \
    X.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( X, y, z, w ) };                                                                        \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2( real_type const x, real_type const y, real_type const z, real_type const w )     \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual1st Y{ y };                                                                                          \
    Y.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, Y, z, w ) };                                                                        \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3( real_type const x, real_type const y, real_type const z, real_type const w )     \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual1st Z{ z };                                                                                          \
    Z.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, y, Z, w ) };                                                                        \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4( real_type const x, real_type const y, real_type const z, real_type const w )     \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual1st W{ w };                                                                                          \
    W.grad = 1;                                                                                                        \
    autodiff::dual1st res{ FUN( x, y, z, W ) };                                                                        \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd X{ x };                                                                                          \
    X.val.grad = 1;                                                                                                    \
    X.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y, z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, Y{ y };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, Y, z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_3( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, Z{ z };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y, Z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_4( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd X{ x }, W{ w };                                                                                  \
    X.val.grad = 1;                                                                                                    \
    W.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( X, y, z, W ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y };                                                                                          \
    Y.val.grad = 1;                                                                                                    \
    Y.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y, z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_3( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y }, Z{ z };                                                                                  \
    Y.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y, Z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_4( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd Y{ y }, W{ w };                                                                                  \
    Y.val.grad = 1;                                                                                                    \
    W.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, Y, z, W ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_3( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd Z{ z };                                                                                          \
    Z.val.grad = 1;                                                                                                    \
    Z.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, y, Z, w ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_4( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd Z{ z }, W{ w };                                                                                  \
    Z.val.grad = 1;                                                                                                    \
    W.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, y, Z, W ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_4( real_type const x, real_type const y, real_type const z, real_type const w )   \
      CONST                                                                                                            \
  {                                                                                                                    \
    autodiff::dual2nd W{ w };                                                                                          \
    W.val.grad = 1;                                                                                                    \
    W.grad.val = 1;                                                                                                    \
    autodiff::dual2nd res{ FUN( x, y, z, W ) };                                                                        \
    return res.grad.grad;                                                                                              \
  }


#define UTILS_AUTODIFF_DERIV_5ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5 ) CONST                                                       \
  {                                                                                                                    \
    autodiff::dual1st X1{ x1 };                                                                                        \
    X1.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( X1, x2, x3, x4, x5 ) };                                                                \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5 ) CONST                                                       \
  {                                                                                                                    \
    autodiff::dual1st X2{ x2 };                                                                                        \
    X2.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, X2, x3, x4, x5 ) };                                                                \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5 ) CONST                                                       \
  {                                                                                                                    \
    autodiff::dual1st X3{ x3 };                                                                                        \
    X3.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, X3, x4, x5 ) };                                                                \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5 ) CONST                                                       \
  {                                                                                                                    \
    autodiff::dual1st X4{ x4 };                                                                                        \
    X4.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, x3, X4, x5 ) };                                                                \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5 ) CONST                                                       \
  {                                                                                                                    \
    autodiff::dual1st X5{ x5 };                                                                                        \
    X5.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, X5 ) };                                                                \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 };                                                                                        \
    X1.val.grad = 1;                                                                                                   \
    X1.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X2{ x2 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X2.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, X2, x3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X3{ x3 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, X3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X4{ x4 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, X4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X5{ x5 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, X5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 };                                                                                        \
    X2.val.grad = 1;                                                                                                   \
    X2.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X3{ x3 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, X3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X4{ x4 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, X4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X5{ x5 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, X5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 };                                                                                        \
    X3.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 }, X4{ x4 };                                                                              \
    X3.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, X4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 }, X5{ x5 };                                                                              \
    X3.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, X5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X4{ x4 };                                                                                        \
    X4.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X4{ x4 }, X5{ x5 };                                                                              \
    X4.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, X5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_5_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5 ) CONST                                 \
  {                                                                                                                    \
    autodiff::dual2nd X5{ x5 };                                                                                        \
    X5.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5 ) };                                                                \
    return res.grad.grad;                                                                                              \
  }


#define UTILS_AUTODIFF_DERIV_6ARG( INLINE, CLASS, PREFIX, FUN, CONST )                                                 \
  INLINE real_type CLASS PREFIX##D_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X1{ x1 };                                                                                        \
    X1.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( X1, x2, x3, x4, x5, x6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X2{ x2 };                                                                                        \
    X2.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, X2, x3, x4, x5, x6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X3{ x3 };                                                                                        \
    X3.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, X3, x4, x5, x6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X4{ x4 };                                                                                        \
    X4.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, x3, X4, x5, x6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X5{ x5 };                                                                                        \
    X5.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, X5, x6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,  \
                                      real_type const x5, real_type const x6 ) CONST                                   \
  {                                                                                                                    \
    autodiff::dual1st X6{ x6 };                                                                                        \
    X6.grad = 1;                                                                                                       \
    autodiff::dual1st res{ FUN( x1, x2, x3, x4, x5, X6 ) };                                                            \
    return res.grad;                                                                                                   \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_1( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 };                                                                                        \
    X1.val.grad = 1;                                                                                                   \
    X1.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_2( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X2{ x2 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X2.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, X2, x3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X3{ x3 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, X3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X4{ x4 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, X4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X5{ x5 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, X5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_1_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X1{ x1 }, X6{ x6 };                                                                              \
    X1.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( X1, x2, x3, x4, x5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_2( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 };                                                                                        \
    X2.val.grad = 1;                                                                                                   \
    X2.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X3{ x3 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, X3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X4{ x4 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, X4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X5{ x5 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, X5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_2_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X2{ x2 }, X6{ x6 };                                                                              \
    X2.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, X2, x3, x4, x5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_3( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 };                                                                                        \
    X3.val.grad = 1;                                                                                                   \
    X3.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 }, X4{ x4 };                                                                              \
    X3.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, X4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 }, X5{ x5 };                                                                              \
    X3.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, X5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_3_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X3{ x3 }, X6{ x6 };                                                                              \
    X3.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, X3, x4, x5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_4( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X4{ x4 };                                                                                        \
    X4.val.grad = 1;                                                                                                   \
    X4.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X4{ x4 }, X5{ x5 };                                                                              \
    X4.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, X5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_4_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X4{ x4 }, X6{ x6 };                                                                              \
    X4.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, X4, x5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_5_5( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X5{ x5 };                                                                                        \
    X5.val.grad = 1;                                                                                                   \
    X5.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5, x6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_5_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X5{ x5 }, X6{ x6 };                                                                              \
    X5.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, X5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }                                                                                                                    \
                                                                                                                       \
  INLINE real_type CLASS PREFIX##D_6_6( real_type const x1, real_type const x2, real_type const x3,                    \
                                        real_type const x4, real_type const x5, real_type const x6 ) CONST             \
  {                                                                                                                    \
    autodiff::dual2nd X6{ x6 };                                                                                        \
    X6.val.grad = 1;                                                                                                   \
    X6.grad.val = 1;                                                                                                   \
    autodiff::dual2nd res{ FUN( x1, x2, x3, x4, x5, X6 ) };                                                            \
    return res.grad.grad;                                                                                              \
  }

#define UTILS_AUTODIFF_DERIV_1ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D( real_type const x ) CONST;                                                                      \
  real_type PREFIX##DD( real_type const x ) CONST;

#define UTILS_AUTODIFF_DERIV_2ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D_1( real_type const x, real_type const y ) CONST;                                                 \
  real_type PREFIX##D_2( real_type const x, real_type const y ) CONST;                                                 \
  real_type PREFIX##D_1_1( real_type const x, real_type const y ) CONST;                                               \
  real_type PREFIX##D_1_2( real_type const x, real_type const y ) CONST;                                               \
  real_type PREFIX##D_2_2( real_type const x, real_type const y ) CONST;

#define UTILS_AUTODIFF_DERIV_3ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D_1( real_type const x, real_type const y, real_type const z ) CONST;                              \
  real_type PREFIX##D_2( real_type const x, real_type const y, real_type const z ) CONST;                              \
  real_type PREFIX##D_3( real_type const x, real_type const y, real_type const z ) CONST;                              \
  real_type PREFIX##D_1_1( real_type const x, real_type const y, real_type const z ) CONST;                            \
  real_type PREFIX##D_1_2( real_type const x, real_type const y, real_type const z ) CONST;                            \
  real_type PREFIX##D_1_3( real_type const x, real_type const y, real_type const z ) CONST;                            \
  real_type PREFIX##D_2_2( real_type const x, real_type const y, real_type const z ) CONST;                            \
  real_type PREFIX##D_2_3( real_type const x, real_type const y, real_type const z ) CONST;                            \
  real_type PREFIX##D_3_3( real_type const x, real_type const y, real_type const z ) CONST;

#define UTILS_AUTODIFF_DERIV_4ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D_1( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;           \
  real_type PREFIX##D_2( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;           \
  real_type PREFIX##D_3( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;           \
  real_type PREFIX##D_4( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;           \
  real_type PREFIX##D_1_1( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_1_2( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_1_3( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_1_4( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_2_2( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_2_3( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_2_4( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_3_3( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_3_4( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;         \
  real_type PREFIX##D_4_4( real_type const x, real_type const y, real_type const z, real_type const w ) CONST;

#define UTILS_AUTODIFF_DERIV_5ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5 ) CONST;                                                                   \
  real_type PREFIX##D_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5 ) CONST;                                                                   \
  real_type PREFIX##D_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5 ) CONST;                                                                   \
  real_type PREFIX##D_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5 ) CONST;                                                                   \
  real_type PREFIX##D_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5 ) CONST;                                                                   \
  real_type PREFIX##D_1_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_1_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_1_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_1_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_1_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_2_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_2_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_2_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_2_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_3_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_3_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_3_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_4_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_4_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;                                                                 \
  real_type PREFIX##D_5_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5 ) CONST;


#define UTILS_AUTODIFF_DERIV_6ARG_DEF( PREFIX, CONST )                                                                 \
  real_type PREFIX##D_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,               \
                         real_type const x5, real_type const x6 ) CONST;                                               \
  real_type PREFIX##D_1_1( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_1_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_1_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_1_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_1_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_1_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_2_2( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_2_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_2_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_2_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_2_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_3_3( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_3_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_3_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_3_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_4_4( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_4_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_4_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_5_5( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_5_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;                                             \
  real_type PREFIX##D_6_6( real_type const x1, real_type const x2, real_type const x3, real_type const x4,             \
                           real_type const x5, real_type const x6 ) CONST;

#endif

#endif

//
// eof: Utils_autodiff.hh
//
