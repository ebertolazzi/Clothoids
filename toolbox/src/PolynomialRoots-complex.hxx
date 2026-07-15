/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2014                                                      |
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

#ifndef POLYNOMIAL_ROOTS_COMPLEX_HXX
#define POLYNOMIAL_ROOTS_COMPLEX_HXX

namespace PolynomialRoots
{

  //! Lightweight complex number type used by the standard-precision solvers.
  class real_complex
  {
    real_type re = 0;
    real_type im = 0;

  public:
    //! Return the real part.
    real_type const & real() const { return re; }
    //! Return the imaginary part.
    real_type const & imag() const { return im; }

    //! Build the zero complex number.
    real_complex() = default;

    //! Build a complex number from real and imaginary parts.
    real_complex( real_type real_part, real_type imag_part = 0 ) : re( real_part ), im( imag_part ) {}

    //! Return the absolute value of a real scalar.
    static inline real_type abs_real( real_type x ) { return std::fabs( x ); }

    //! Add a real scalar in place.
    inline real_complex & operator+=( real_type rhs )
    {
      re += rhs;
      return *this;
    }

    //! Add a complex number in place.
    inline real_complex & operator+=( real_complex const & rhs )
    {
      re += rhs.re;
      im += rhs.im;
      return *this;
    }

    //! Subtract a real scalar in place.
    inline real_complex & operator-=( real_type rhs )
    {
      re -= rhs;
      return *this;
    }

    //! Subtract a complex number in place.
    inline real_complex & operator-=( real_complex const & rhs )
    {
      re -= rhs.re;
      im -= rhs.im;
      return *this;
    }

    //! Multiply by a real scalar in place.
    inline real_complex & operator*=( real_type rhs )
    {
      re *= rhs;
      im *= rhs;
      return *this;
    }

    //! Multiply by a complex number in place.
    inline real_complex & operator*=( real_complex const & rhs )
    {
      long double const a = static_cast<long double>( re );
      long double const b = static_cast<long double>( im );
      long double const c = static_cast<long double>( rhs.re );
      long double const d = static_cast<long double>( rhs.im );

      long double const real = a * c - b * d;
      long double const imag = a * d + b * c;

      re = static_cast<real_type>( real );
      im = static_cast<real_type>( imag );

      return *this;
    }

    //! Divide by a real scalar in place.
    inline real_complex & operator/=( real_type rhs )
    {
      re /= rhs;
      im /= rhs;
      return *this;
    }

    //! Divide by a complex number in place using a scaled algorithm.
    inline real_complex & operator/=( real_complex const & rhs )
    {
      long double const a = static_cast<long double>( re );
      long double const b = static_cast<long double>( im );
      long double const c = static_cast<long double>( rhs.re );
      long double const d = static_cast<long double>( rhs.im );

      long double const abs_c = std::fabs( c );
      long double const abs_d = std::fabs( d );

      long double real;
      long double imag;

      if ( abs_c >= abs_d )
      {
        long double const r   = d / c;
        long double const den = c + d * r;

        real = ( a + b * r ) / den;
        imag = ( b - a * r ) / den;
      }
      else
      {
        long double const r   = c / d;
        long double const den = d + c * r;

        real = ( a * r + b ) / den;
        imag = ( b * r - a ) / den;
      }

      re = static_cast<real_type>( real );
      im = static_cast<real_type>( imag );

      return *this;
    }
  };

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  //! Lightweight complex number type used by the high-precision solvers.
  class quad_complex
  {
    quad_real re = 0;
    quad_real im = 0;

  public:
    //! Return the real part.
    quad_real const & real() const { return re; }
    //! Return the imaginary part.
    quad_real const & imag() const { return im; }

    //! Build the zero complex number.
    quad_complex() = default;

    //! Build a complex number from high-precision real and imaginary parts.
    quad_complex( quad_real const & real_part, quad_real const & imag_part ) : re( real_part ), im( imag_part ) {}

    //! Build a complex number from standard-precision parts.
    quad_complex( real_type real_part, real_type imag_part = 0 ) : re( real_part ), im( imag_part ) {}

    //! Promote a standard-precision complex number.
    explicit quad_complex( real_complex const & z ) : re( z.real() ), im( z.imag() ) {}

    //! Return the absolute value of a high-precision real scalar.
    static inline quad_real abs_real( quad_real const & x ) { return x >= quad_real{ 0 } ? x : -x; }

    //! Add a high-precision real scalar in place.
    inline quad_complex & operator+=( quad_real const & rhs )
    {
      re += rhs;
      return *this;
    }

    //! Add a complex number in place.
    inline quad_complex & operator+=( quad_complex const & rhs )
    {
      re += rhs.re;
      im += rhs.im;
      return *this;
    }

    //! Subtract a high-precision real scalar in place.
    inline quad_complex & operator-=( quad_real const & rhs )
    {
      re -= rhs;
      return *this;
    }

    //! Subtract a complex number in place.
    inline quad_complex & operator-=( quad_complex const & rhs )
    {
      re -= rhs.re;
      im -= rhs.im;
      return *this;
    }

    //! Multiply by a high-precision real scalar in place.
    inline quad_complex & operator*=( quad_real const & rhs )
    {
      re *= rhs;
      im *= rhs;
      return *this;
    }

    //! Multiply by a complex number in place.
    inline quad_complex & operator*=( quad_complex const & rhs )
    {
      quad_real const a = re;
      quad_real const b = im;
      quad_real const c = rhs.re;
      quad_real const d = rhs.im;

      quad_real const real = a * c - b * d;
      quad_real const imag = a * d + b * c;

      re = real;
      im = imag;

      return *this;
    }

    //! Divide by a high-precision real scalar in place.
    inline quad_complex & operator/=( quad_real const & rhs )
    {
      re /= rhs;
      im /= rhs;
      return *this;
    }

    //! Divide by a complex number in place using a scaled algorithm.
    inline quad_complex & operator/=( quad_complex const & rhs )
    {
      quad_real const a = re;
      quad_real const b = im;
      quad_real const c = rhs.re;
      quad_real const d = rhs.im;

      quad_real const abs_c = abs_real( c );
      quad_real const abs_d = abs_real( d );

      quad_real real;
      quad_real imag;

      if ( abs_c >= abs_d )
      {
        quad_real const r   = d / c;
        quad_real const den = c + d * r;

        real = ( a + b * r ) / den;
        imag = ( b - a * r ) / den;
      }
      else
      {
        quad_real const r   = c / d;
        quad_real const den = d + c * r;

        real = ( a * r + b ) / den;
        imag = ( b * r - a ) / den;
      }

      re = real;
      im = imag;

      return *this;
    }
  };
#endif

  //! Return the sum of two standard-precision complex values.
  inline real_complex operator+( real_complex const & a, real_complex const & b )
  { return { a.real() + b.real(), a.imag() + b.imag() }; }

  //! Return the sum of a standard-precision complex value and a real scalar.
  inline real_complex operator+( real_complex const & a, real_type b )
  { return { a.real() + b, a.imag() }; }

  //! Return the sum of a real scalar and a standard-precision complex value.
  inline real_complex operator+( real_type a, real_complex const & b )
  { return { a + b.real(), b.imag() }; }

  //! Return the difference of two standard-precision complex values.
  inline real_complex operator-( real_complex const & a, real_complex const & b )
  { return { a.real() - b.real(), a.imag() - b.imag() }; }

  //! Return the difference between a standard-precision complex value and a real scalar.
  inline real_complex operator-( real_complex const & a, real_type b )
  { return { a.real() - b, a.imag() }; }

  //! Return the difference between a real scalar and a standard-precision complex value.
  inline real_complex operator-( real_type a, real_complex const & b )
  { return { a - b.real(), -b.imag() }; }

  //! Return the additive inverse of a standard-precision complex value.
  inline real_complex operator-( real_complex const & a )
  { return { -a.real(), -a.imag() }; }

  //! Return the product of two standard-precision complex values.
  inline real_complex operator*( real_complex const & a, real_complex const & b )
  {
    real_complex res{ a };
    res *= b;
    return res;
  }

  //! Return the product of a standard-precision complex value and a real scalar.
  inline real_complex operator*( real_complex const & a, real_type b )
  {
    real_complex res{ a };
    res *= b;
    return res;
  }

  //! Return the product of a real scalar and a standard-precision complex value.
  inline real_complex operator*( real_type a, real_complex const & b )
  { return { a * b.real(), a * b.imag() }; }

  //! Return the quotient of two standard-precision complex values.
  inline real_complex operator/( real_complex const & a, real_complex const & b )
  {
    real_complex res{ a };
    res /= b;
    return res;
  }

  //! Return the quotient of a standard-precision complex value and a real scalar.
  inline real_complex operator/( real_complex const & a, real_type b )
  {
    real_complex res{ a };
    res /= b;
    return res;
  }

  //! Return the quotient of a real scalar and a standard-precision complex value.
  inline real_complex operator/( real_type a, real_complex const & b )
  {
    real_complex res{ a, 0 };
    res /= b;
    return res;
  }

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  //! Return the sum of two high-precision complex values.
  inline quad_complex operator+( quad_complex const & a, quad_complex const & b )
  { return { a.real() + b.real(), a.imag() + b.imag() }; }

  //! Return the sum of a high-precision complex value and a real scalar.
  inline quad_complex operator+( quad_complex const & a, quad_real const & b )
  { return { a.real() + b, a.imag() }; }

  //! Return the sum of a real scalar and a high-precision complex value.
  inline quad_complex operator+( quad_real const & a, quad_complex const & b )
  { return { a + b.real(), b.imag() }; }

  //! Return the difference of two high-precision complex values.
  inline quad_complex operator-( quad_complex const & a, quad_complex const & b )
  { return { a.real() - b.real(), a.imag() - b.imag() }; }

  //! Return the difference between a high-precision complex value and a real scalar.
  inline quad_complex operator-( quad_complex const & a, quad_real const & b )
  { return { a.real() - b, a.imag() }; }

  //! Return the difference between a real scalar and a high-precision complex value.
  inline quad_complex operator-( quad_real const & a, quad_complex const & b )
  { return { a - b.real(), -b.imag() }; }

  //! Return the additive inverse of a high-precision complex value.
  inline quad_complex operator-( quad_complex const & a )
  { return { -a.real(), -a.imag() }; }

  //! Return the product of two high-precision complex values.
  inline quad_complex operator*( quad_complex const & a, quad_complex const & b )
  {
    quad_complex res{ a };
    res *= b;
    return res;
  }

  //! Return the product of a high-precision complex value and a real scalar.
  inline quad_complex operator*( quad_complex const & a, quad_real const & b )
  {
    quad_complex res{ a };
    res *= b;
    return res;
  }

  //! Return the product of a real scalar and a high-precision complex value.
  inline quad_complex operator*( quad_real const & a, quad_complex const & b )
  { return { a * b.real(), a * b.imag() }; }

  //! Return the quotient of two high-precision complex values.
  inline quad_complex operator/( quad_complex const & a, quad_complex const & b )
  {
    quad_complex res{ a };
    res /= b;
    return res;
  }

  //! Return the quotient of a high-precision complex value and a real scalar.
  inline quad_complex operator/( quad_complex const & a, quad_real const & b )
  {
    quad_complex res{ a };
    res /= b;
    return res;
  }

  //! Return the quotient of a real scalar and a high-precision complex value.
  inline quad_complex operator/( quad_real const & a, quad_complex const & b )
  {
    quad_complex res{ a, quad_real{ 0 } };
    res /= b;
    return res;
  }

  inline quad_real sqrt( quad_real const & value )
  {
    using boost::multiprecision::sqrt;
    return sqrt( value );
  }

  inline quad_real cbrt( quad_real const & value )
  {
    using boost::multiprecision::cbrt;
    return cbrt( value );
  }

  inline quad_real abs2( quad_complex const & value )
  {
    quad_real const ax = abs( value.real() );
    quad_real const ay = abs( value.imag() );

    if ( ax == quad_real{ 0 } ) return ay * ay;
    if ( ay == quad_real{ 0 } ) return ax * ax;

    quad_real const scale = ax >= ay ? ax : ay;
    quad_real const x     = value.real() / scale;
    quad_real const y     = value.imag() / scale;

    return ( x * x + y * y ) * scale * scale;
  }

  inline quad_real hypot( quad_real const & x, quad_real const & y )
  {
    quad_real const ax = abs( x );
    quad_real const ay = abs( y );

    if ( ax == quad_real{ 0 } ) return ay;
    if ( ay == quad_real{ 0 } ) return ax;

    quad_real const scale = ax >= ay ? ax : ay;
    quad_real const X     = x / scale;
    quad_real const Y     = y / scale;

    return sqrt( X * X + Y * Y ) * scale;
  }

  inline quad_real abs( quad_complex const & value )
  { return boost::math::hypot( value.real(), value.imag() ); }
#endif

  inline real_type abs( real_complex const & value )
  { return std::hypot( value.real(), value.imag() ); }

  inline real_type abs2( real_complex const & value )
  { return value.real() * value.real() + value.imag() * value.imag(); }

}  // namespace PolynomialRoots

#endif
