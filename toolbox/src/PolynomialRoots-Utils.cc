/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

#include "PolynomialRoots.hh"

namespace PolynomialRoots
{

  template <> real_type machepsiT() { return std::numeric_limits<real_type>::epsilon(); }
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <> quad_real machepsiT() { return std::numeric_limits<quad_real>::epsilon(); }
#endif
  template <> real_type toleranceT() { return 1e-12; }
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <> quad_real toleranceT() { return 1e-20; }
#endif

  using std::abs;
  using std::pow;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  template <typename T_real> T_real eval_poly( T_real const op[], integer const Degree, T_real const & x )
  {
    if ( abs( x ) > 1 )
    {
      T_real res = op[Degree];
      T_real xn  = 1;
      for ( integer i = 1; i <= Degree; ++i )
      {
        res /= x;
        res += op[Degree - i];
        xn *= x;
      }
      res *= xn;
      return res;
    }
    T_real res( op[0] );
    for ( integer i = 1; i <= Degree; ++i )
    {
      res *= x;
      res += op[i];
    }
    return res;
  }

  template double    eval_poly( double const op[], integer const Degree, double const & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_real eval_poly( quad_real const op[], integer const Degree, quad_real const & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of Newton step
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  template <typename T_real> bool Newton_step( T_real const op[], integer const Degree, T_real & x )
  {
    T_real p, dp;
    if ( abs( x ) > 1 )
    {
      T_real xn = 1;
      p         = op[Degree];
      for ( integer i{ 1 }; i <= Degree; ++i )
      {
        p /= x;
        p += op[Degree - i];
        xn *= x;
      }
      p *= xn;

      dp = op[Degree];
      xn = 1;
      for ( integer i = 2; i <= Degree; ++i )
      {
        dp /= x;
        dp += ( Degree - i ) * op[Degree - i];
        xn *= x;
      }
      dp *= xn;
    }
    else
    {
      dp = Degree * op[0];
      p  = op[0];
      for ( integer i = 1; i < Degree; ++i )
      {
        p *= x;
        p += op[i];
        dp *= x;
        dp += ( Degree - i ) * op[i];
      }
      p *= x;
      p += op[Degree];
    }
    x -= p / dp;
    return true;
  }

  template bool Newton_step( double const op[], integer const Degree, double & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template bool Newton_step( quad_real const op[], integer const Degree, quad_real & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial and its derivative
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  template <typename T_real>
  void eval_poly_Dpoly( T_real const op[], integer const Degree, T_real const & x, T_real & p, T_real & dp )
  {
    if ( abs( x ) > 1 )
    {
      T_real const   y  = T_real( 1 ) / x;
      T_real const * pc = op + Degree;

      p  = *pc--;
      dp = 0;
      for ( integer i{ 1 }; i <= Degree; ++i )
      {
        dp = dp * y + p;
        p  = p * y + *pc--;
      }

      T_real xn = 1;
      for ( integer i = 0; i < Degree; ++i ) xn *= x;

      dp = static_cast<T_real>( Degree ) * p * ( xn / x ) - dp * ( xn / ( x * x ) );
      p *= xn;
    }
    else
    {
      p  = op[0];
      dp = op[0];
      p  = p * x + op[1];
      for ( integer i = 2; i <= Degree; ++i )
      {
        dp = dp * x + p;  // c_i = c_{i-1}*x + b_{i-1}
        p  = p * x + op[i];
      }
    }
  }

  template void eval_poly_Dpoly( double const op[], integer const Degree, double const & x, double & p, double & dp );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template void eval_poly_Dpoly(
    quad_real const   op[],
    integer const     Degree,
    quad_real const & x,
    quad_real &       p,
    quad_real &       dp );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  template <typename T_real, typename T_complex>
  T_complex eval_poly_complex( T_real const op[], integer const Degree, T_complex const & x )
  {
    if ( abs( x ) > 1 )
    {
      T_complex res( op[Degree], 0 );
      T_complex xn( 1, 0 );
      for ( integer i = 1; i <= Degree; ++i )
      {
        res /= x;
        res += op[Degree - i];
        xn *= x;
      }
      res *= xn;
      return res;
    }
    T_complex res( op[0], 0 );
    for ( integer i = 1; i <= Degree; ++i )
    {
      res *= x;
      res += op[i];
    }
    return res;
  }

  template real_complex eval_poly_complex( real_type const op[], integer const Degree, real_complex const & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_complex eval_poly_complex( quad_real const op[], integer const Degree, quad_complex const & x );
#endif

}  // namespace PolynomialRoots

// EOF: PolynomialRoots-Utils.cc
