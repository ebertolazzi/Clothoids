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

#ifndef POLYNOMIAL_ROOTS_HH
#define POLYNOMIAL_ROOTS_HH

#include <cmath>
#include <cfloat>
#include <iostream>

#include <complex>
#include <cstdint>
#include <format>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <sstream>
#include <string>

//!
//! Implementation of Flocke algorithm for roots
//! of 3rd and 4th degree polynomials.
//!
//! There are 3 classed for 2nd, 3rd and 4th degree polynomial.
//! An experimental translation to C++ of a C implementation of
//! Jenkins--Traub algorithm is available.
//!
//! **References**
//!
//! - **N.Flocke**
//!   Algorithm 954: An Accurate and Efficient Cubic and Quartic
//!   Equation Solver for Physical Applications
//!   ACM TOMS, vol 41, n.4, 2015
//!
//! - **M.A. Jenkins and J.F.Traub**
//!   A Three-Stage Algorithm for Real Polynomials Using Quadratic Iteration
//!   SIAM Journal on Numerical Analysis
//!   Vol. 7, No.4 (Dec., 1970), pp.545-566
//!

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  #include <boost/multiprecision/cpp_bin_float.hpp>
  #include <boost/math/special_functions/cbrt.hpp>
#endif

namespace PolynomialRoots
{

  //! Scalar type used by the standard-precision API.
  using real_type    = double;
  //! Integer type used for degrees, counts and indices.
  using integer      = int;
  //! Output stream type used by diagnostic methods.
  using ostream_type = std::basic_ostream<char>;
  //! Input stream type reserved for formatted input helpers.
  using istream_type = std::basic_istream<char>;

  //! Maximum degree accepted by the Jenkins-Traub entry point.
  inline constexpr integer MAXDEGREE = 100;

  //! Lightweight assertion helper used by the solvers to validate inputs.
  //!
  //! \param[in] cond  condition that must hold
  //! \param[in] fmt   format string used to build the diagnostic message
  //! \param[in] args  arguments interpolated in \p fmt
  template <typename... Args> inline void root_assert( bool cond, std::format_string<Args...> fmt, Args &&... args )
  {
    if ( !cond ) std::runtime_error( std::format( fmt, std::forward<Args>( args )... ) );
  }

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  //! High-precision scalar type used when multiprecision support is enabled.
  using quad_real = boost::multiprecision::cpp_bin_float_100;
#endif

}  // namespace PolynomialRoots

namespace std
{

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <>

  struct formatter<PolynomialRoots::quad_real, char> : formatter<std::string, char>

  {
    auto format( PolynomialRoots::quad_real const & x, format_context & ctx ) const
    {
      std::ostringstream ss;
      constexpr int      digits = std::numeric_limits<PolynomialRoots::quad_real>::max_digits10 > 0
                                    ? std::numeric_limits<PolynomialRoots::quad_real>::max_digits10
                                    : 36;
      ss << std::setprecision( digits ) << x;
      return formatter<std::string, char>::format( ss.str(), ctx );
    }
  };
#endif
}  // namespace std

#include "PolynomialRoots-complex.hxx"

namespace PolynomialRoots
{
  using std::isfinite;

  //! Return the machine epsilon associated with the scalar type.
  template <typename T_real> T_real machepsiT();
  //! Return the default residual tolerance associated with the scalar type.
  template <typename T_real> T_real toleranceT();

  //! Evaluate a real polynomial at a real point with a numerically stable scheme.
  //!
  //! \param[in] op      polynomial coefficients in descending powers
  //! \param[in] Degree  polynomial degree
  //! \param[in] x       evaluation point
  //! \return            value of the polynomial at \p x
  template <typename T_real> T_real eval_poly( T_real const op[], integer Degree, T_real const & x );

  //! Evaluate a real polynomial and its first derivative at a real point.
  //!
  //! \param[in]  op      polynomial coefficients in descending powers
  //! \param[in]  Degree  polynomial degree
  //! \param[in]  x       evaluation point
  //! \param[out] p       polynomial value at \p x
  //! \param[out] dp      derivative value at \p x
  template <typename T_real>
  void eval_poly_Dpoly( T_real const op[], integer Degree, T_real const & x, T_real & p, T_real & dp );

  //! Apply one Newton update to a real polynomial root estimate.
  //!
  //! \param[in]     op      polynomial coefficients in descending powers
  //! \param[in]     Degree  polynomial degree
  //! \param[in,out] x       current estimate, overwritten with the updated value
  //! \return                `true` after performing the update
  template <typename T_real> bool Newton_step( T_real const op[], integer Degree, T_real & x );

  //! Evaluate a real polynomial at a complex point.
  //!
  //! \param[in] op      polynomial coefficients in descending powers
  //! \param[in] Degree  polynomial degree
  //! \param[in] x       complex evaluation point
  //! \return            value of the polynomial at \p x
  template <typename T_real, typename T_complex>
  T_complex eval_poly_complex( T_real const op[], integer Degree, T_complex const & x );

  //!
  //! Find roots of a generic polynomial using Jenkins-Traub method
  //!
  //! \param[in]  op     the coefficients of the polynomial
  //! \param[in]  Degree degree of the polynomial
  //! \param[out] zeror  real part of the roots
  //! \param[out] zeroi  imaginary part of the roots
  //!
  //! \return error code, 0 OK
  //!
  [[nodiscard]] int roots( real_type const * op, integer Degree, real_type * zeror, real_type * zeroi );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
   |    ___                  _           _   _
   |   / _ \ _   _  __ _  __| |_ __ __ _| |_(_) ___
   |  | | | | | | |/ _` |/ _` | '__/ _` | __| |/ __|
   |  | |_| | |_| | (_| | (_| | | | (_| | |_| | (__
   |   \__\_\\__,_|\__,_|\__,_|_|  \__,_|\__|_|\___|
   |
   |  A * x^2 + B * x + C
  \*/
  //! Quadratic polynomial class
  //!
  //! **Constructor**
  //!
  //! ```{cpp}
  //!   double a = 1;
  //!   double b = 2;
  //!   double c = 3;
  //!   Quadratic q(a,b,c); // build an solve `a x^2 + b x + c = 0`
  //!
  //!   Quadratic q;
  //!   q.setup(a,b,c); // build an solve `a x^2 + b x + c = 0`
  //! ```
  //!
  //! **Get kind of solution**
  //!
  //! ```{cpp}
  //!   int  nroots            = q.num_roots();
  //!   bool has_complex_root  = q.complex_root();
  //!   bool has_a_double_root = q.double_root();
  //! ```
  //!
  //! **Get real roots**
  //!
  //! ```{cpp}
  //!   double r_min = 0;
  //!   double r_max = 2;
  //!   double r[2];
  //!   int nroots;
  //!   nroots = p.getRealRoots( r );
  //!   nroots = p.getPositiveRoots( r );
  //!   nroots = p.getNegativeRoots( r );
  //!   nroots = p.getRootsInRange( r_min, r_max, r );
  //!   nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //! ```
  //!
  //! **Get roots**
  //!
  //! ```{cpp}
  //!   real_type    r0 = p.real_root0();
  //!   real_type    r1 = p.real_root1();
  //!   real_complex r0 = p.root0();
  //!   real_complex r1 = p.root1();
  //!
  //!   real_complex r;
  //!   double re, im;
  //!   p.getRoot0( re, im );
  //!   p.getRoot0( r );
  //!   p.getRoot1( re, im );
  //!   p.getRoot1( r );
  //! ```
  //!
  //! **Evaluate polynomial**
  //!
  //! ```{cpp}
  //!   {double or complex} v, x;
  //!   v = p.eval( x );
  //!
  //!   p.eval( x, p, dp );
  //! ```
  //!
  //! **Information**
  //!
  //! ```{cpp}
  //!   p.info( cout );
  //!   bool ok = p.check( cout );
  //! ```
  //!
  template <typename T_real, typename T_complex> class QuadraticT
  {
    T_real  m_ABC[3]{ 0, 0, 0 };
    T_real  m_r0   = 0;
    T_real  m_r1   = 0;
    integer m_nrts = 0;
    bool    m_cplx = false;
    bool    m_dblx = false;

    void find_roots();

  public:
    using value_type   = T_real;
    using complex_type = T_complex;

    //! Build an empty quadratic solver instance.
    QuadraticT() = default;

    //!
    //! Build the object that store the roots
    //! of the quadratic polynomial
    //!
    //! \f[ q(x) = a x^2 + b x + c \f]
    //!
    //! \param[in] a leading coefficient of \f$ q(x) \f$
    //! \param[in] b coefficient of \f$ x \f$
    //! \param[in] c coefficient of \f$ x^0 \f$
    //!
    QuadraticT( T_real const & a, T_real const & b, T_real const & c )
    {
      m_ABC[0] = a;
      m_ABC[1] = b;
      m_ABC[2] = c;
      // find roots only on finite values
      root_assert(
        isfinite( a ) && isfinite( b ) && isfinite( c ),
        "QuadraticT( a={}, b={}, c={} ) arguments must be finite!",
        a, b, c
      );
      find_roots();
    }

    //!
    //! Setup the object that store the roots of the quadratic polynomial
    //!
    //! \f[ q(x) = a x^2 + b x + c \f]
    //!
    //! \param[in] a leading coefficient of \f$ q(x) \f$
    //! \param[in] b coefficient of \f$ x \f$
    //! \param[in] c coefficient of \f$ x^0 \f$
    //!
    void setup( T_real const & a, T_real const & b, T_real const & c )
    {
      m_ABC[0] = a;
      m_ABC[1] = b;
      m_ABC[2] = c;
      find_roots();
    }

    //!
    //! Return the number of computed roots of
    //!
    //! \f[ q(x) = a x^2 + b x + c \f]
    //!
    //! Normally return 2 but if e.g. \f$ a = 0 \f$
    //! return 1 or less depending on the values of \f$ a, b, c \f$
    //!
    //! \return number of computed roots
    //!
    integer num_roots() const { return m_nrts; }

    //! Alias of `num_roots()`.
    integer numRoots() const { return m_nrts; }

    //!
    //! true if roots are complex conjugated
    //!
    bool complex_root() const { return m_cplx; }

    //! Alias of `complex_root()`.
    bool complexRoot() const { return m_cplx; }

    //!
    //! true if \f$ p(x) = a x^2 + b x + c = (x-r)^2 \f$
    //!
    bool double_root() const { return m_dblx; }

    //! Alias of `double_root()`.
    bool doubleRoot() const { return m_dblx; }

    //!
    //! Get the real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1 or 2
    //!
    integer get_real_roots( T_real r[] ) const;

    //! Alias of `get_real_roots()`.
    integer getRealRoots( T_real r[] ) const { return get_real_roots( r ); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1 or 2
    //!
    integer get_positive_roots( T_real r[] ) const;

    //! Alias of `get_positive_roots()`.
    integer getPositiveRoots( T_real r[] ) const { return get_positive_roots( r ); }

    //!
    //! Get negative real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1 or 2
    //!
    integer get_negative_roots( T_real r[] ) const;

    //! Alias of `get_negative_roots()`.
    integer getNegativeRoots( T_real r[] ) const { return get_negative_roots( r ); }

    //!
    //! Get real roots in a closed range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range \f$ [a,b] \f$
    //!
    integer get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_range()`.
    integer getRootsInRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range \f$ (a,b) \f$
    //!
    integer get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_open_range()`.
    integer getRootsInOpenRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_open_range( a, b, r ); }

    //!
    //! Return the first stored real component.
    //!
    //! When the polynomial has a complex conjugate pair, this is the common
    //! real part of the pair returned by `root0()` and `root1()`.
    //!
    T_real real_root0() const { return m_r0; }

    //! Return the second stored real component.
    //!
    //! When the polynomial has a complex conjugate pair, this is the positive
    //! imaginary part paired with `real_root0()`. Otherwise it is the second
    //! real root.
    T_real real_root1() const { return m_r1; }

    //! Return the first root as a complex value.
    T_complex root0() const { return m_cplx ? T_complex( m_r0, m_r1 ) : T_complex( m_r0, 0 ); }

    //! Return the second root as a complex value.
    T_complex root1() const { return m_cplx ? T_complex( m_r0, -m_r1 ) : T_complex( m_r1, 0 ); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] re the first complex root, real part
    //! \param[out] im the first complex root, imaginary part
    //!
    void get_root0( T_real & re, T_real & im ) const
    {
      if ( m_cplx )
      {
        re = m_r0;
        im = m_r1;
      }
      else
      {
        re = m_r0;
        im = 0;
      }
    }

    //! Alias of `get_root0()`.
    void getRoot0( T_real & re, T_real & im ) const { return get_root0( re, im ); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] r the first complex root
    //!
    void get_root0( T_complex & r ) const { r = m_cplx ? T_complex( m_r0, m_r1 ) : T_complex( m_r0, 0 ); }

    //! Alias of `get_root0()`.
    void getRoot0( T_complex & r ) const { return get_root0( r ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] re the second complex root, real part
    //! \param[out] im the second complex root, imaginary part
    //!
    void get_root1( T_real & re, T_real & im ) const
    {
      if ( m_cplx )
      {
        re = m_r0;
        im = -m_r1;
      }
      else
      {
        re = m_r1;
        im = 0;
      }
    }

    //! Alias of `get_root1()`.
    void getRoot1( T_real & re, T_real & im ) const { return get_root1( re, im ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] r the second complex root
    //!
    void get_root1( T_complex & r ) const { r = m_cplx ? T_complex( m_r0, -m_r1 ) : T_complex( m_r1, 0 ); }

    //! Alias of `get_root1()`.
    void getRoot1( T_complex & r ) const { return get_root1( r ); }

    //! Return the \p i-th root as a complex value.
    //!
    //! \param[in] i root index in `[0,1]`
    //! \return      requested root, or zero for an invalid index
    T_complex root( integer const i ) const
    {
      switch ( i )
      {
        case 0: return root0();
        case 1: return root1();
      }
      return 0;
    }

    //! Store the \p i-th root in split real/imaginary form.
    //!
    //! \param[in]  i   root index in `[0,1]`
    //! \param[out] re  real part of the selected root
    //! \param[out] im  imaginary part of the selected root
    void get_root( integer const i, T_real & re, T_real & im ) const
    {
      switch ( i )
      {
        case 0: return get_root0( re, im );
        case 1: return get_root1( re, im );
      }
    }

    //!
    //! Evaluate the quadratic polynomial at a real point.
    //!
    //! \param[in] x  value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    //!
    T_real eval( T_real const & x ) const { return eval_poly<T_real>( m_ABC, 2, x ); }

    //! Evaluate the quadratic polynomial at a complex point.
    //!
    //! \param[in] x  complex value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    T_complex eval( T_complex const & x ) const { return eval_poly_complex<T_real, T_complex>( m_ABC, 2, x ); }

    //!
    //! Evaluate the polynomial with its derivative
    //!
    //! \param[in]  x   value where compute \f$ p(x) \f$
    //! \param[out] p   value \f$ p(x) \f$
    //! \param[out] dp  value \f$ p'(x) \f$
    //!
    void eval( T_real const & x, T_real & p, T_real & dp ) const { eval_poly_Dpoly<T_real>( m_ABC, 2, x, p, dp ); }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void info( ostream_type & s ) const;

    //!
    //! Check tolerance and quality of the computed roots
    //!
    [[nodiscard]] bool check( ostream_type & s ) const;
  };

  using Quadratic   = QuadraticT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using QuadraticHQ = QuadraticT<quad_real, quad_complex>;
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
   |    ____      _     _
   |   / ___|   _| |__ (_) ___
   |  | |  | | | | '_ \| |/ __|
   |  | |__| |_| | |_) | | (__
   |   \____\__,_|_.__/|_|\___|
   |
   |  A * x^3 + B * x^2 + C * x + D
  \*/
  //! Cubic polynomial class
  //!
  //! **Constructor**
  //!
  //! ```{cpp}
  //!   double a = 1;
  //!   double b = 2;
  //!   double c = 3;
  //!   double d = 3;
  //!   Cubic p(a,b,c,d); // build an solve a x^3 + b x^2 + c x + d = 0
  //!
  //!   Cubic p;
  //!   p.setup(a,b,c,d); // build an solve a x^3 + b x^2 + c x + d = 0
  //! ```
  //!
  //! **Get kind of solution**
  //!
  //! ```{cpp}
  //!   int  nroots = p.num_roots();
  //!   bool has_complex_root  = p.complex_root();
  //!   bool has_a_double_root = p.double_root();
  //!   bool has_a_triple_root = p.triple_root();
  //! ```
  //!
  //! **Get real roots**
  //!
  //! ```{cpp}
  //!   double r_min = 0;
  //!   double r_max = 2;
  //!   double r[3];
  //!   int nroots;
  //!   nroots = p.getRealRoots( r );
  //!   nroots = p.getPositiveRoots( r );
  //!   nroots = p.getNegativeRoots( r );
  //!   nroots = p.getRootsInRange( r_min, r_max, r );
  //!   nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //! ```
  //!
  //! **Get roots**
  //!
  //! ```{cpp}
  //!   real_type    r0 = p.real_root0();
  //!   real_type    r1 = p.real_root1();
  //!   real_type    r2 = p.real_root2();
  //!   real_complex r0 = p.root0();
  //!   real_complex r1 = p.root1();
  //!   real_complex r2 = p.root2();
  //!
  //!   real_complex r;
  //!   real_type    re, im;
  //!   p.getRoot0( re, im );
  //!   p.getRoot0( r );
  //!   p.getRoot1( re, im );
  //!   p.getRoot1( r );
  //!   p.getRoot2( re, im );
  //!   p.getRoot2( r );
  //! ```
  //!
  //! **Evaluate polynomial**
  //!
  //! ```{cpp}
  //!   {double or complex} v, x;
  //!   v = p.eval( x );
  //!
  //!   p.eval( x, p, dp );
  //! ```
  //!
  //! **Information**
  //!
  //! ```{cpp}
  //!   p.info( cout );
  //!   bool ok = p.check( cout );
  //! ```
  //!
  template <typename T_real, typename T_complex> class CubicT
  {
    T_real  m_ABCD[4]{ 0, 0, 0, 0 };
    T_real  m_r0   = 0;
    T_real  m_r1   = 0;
    T_real  m_r2   = 0;
    integer m_nrts = 0;
    integer m_iter = 0;
    bool    m_cplx = false;  // complex root
    bool    m_dblx = false;  // double root
    bool    m_trpx = false;  // triple root

    void find_roots();

  public:
    using value_type   = T_real;
    using complex_type = T_complex;

    //!
    //! Build an empty instance of Cubic polynomial solver
    //!
    //! Build an empty cubic solver instance.
    CubicT() = default;

    //!
    //! Build the instance of the calss and compute the roots of cubic polynomial
    //! \f$ a x^3 + b x^2 + c x + d \f$
    //!
    //! \param[in] a coefficient of \f$ x^3 \f$
    //! \param[in] b coefficient of \f$ x^2 \f$
    //! \param[in] c coefficient of \f$ x   \f$
    //! \param[in] d coefficient of \f$ x^0 \f$
    //!
    CubicT( T_real const & a, T_real const & b, T_real const & c, T_real const & d )
    {
      m_ABCD[0] = a;
      m_ABCD[1] = b;
      m_ABCD[2] = c;
      m_ABCD[3] = d;
      // find roots only on finite values
      root_assert(
        isfinite( a ) && isfinite( b ) && isfinite( c ) && isfinite( d ),
        "CubicT( a={}, b={}, c={}, d={} ) arguments must be finite!",
        a, b, c, d
      );
      find_roots();
    }

    //!
    //! Compute the roots of cubic polynomial
    //! \f$ a x^3 + b x^2 + c x + d \f$
    //!
    //! \param[in] a coefficient of \f$ x^3 \f$
    //! \param[in] b coefficient of \f$ x^2 \f$
    //! \param[in] c coefficient of \f$ x   \f$
    //! \param[in] d coefficient of \f$ x^0 \f$
    //!
    void setup( T_real const & a, T_real const & b, T_real const & c, T_real const & d )
    {
      m_ABCD[0] = a;
      m_ABCD[1] = b;
      m_ABCD[2] = c;
      m_ABCD[3] = d;
      m_nrts    = 0;
      m_iter    = 0;
      m_cplx    = false;  // complex root
      m_dblx    = false;  // double root
      m_trpx    = false;  // triple root
      root_assert(
        isfinite( a ) && isfinite( b ) && isfinite( c ) && isfinite( d ),
        "CubicT::setup( a={}, b={}, c={}, d={} ) arguments must be finite!",
        a, b, c, d
      );
      find_roots();
    }

    //!
    //! Number of found roots.
    //!
    integer num_roots() const { return m_nrts; }

    //! Alias of `num_roots()`.
    integer numRoots() const { return m_nrts; }

    //!
    //! Has complex roots?
    //!
    bool complex_root() const { return m_cplx; }

    //! Alias of `complex_root()`.
    bool complexRoot() const { return m_cplx; }

    //!
    //! Has a double root?
    //!
    bool double_root() const { return m_dblx; }

    //! Alias of `double_root()`.
    bool doubleRoot() const { return m_dblx; }

    //!
    //! Has a triple root?
    //!
    bool triple_root() const { return m_trpx; }

    //! Alias of `triple_root()`.
    bool tripleRoot() const { return m_trpx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2 or 3
    //!
    integer get_real_roots( T_real r[] ) const;

    //! Alias of `get_real_roots()`.
    integer getRealRoots( T_real r[] ) const { return get_real_roots( r ); }

    //!
    //! Get positive real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2 or 3
    //!
    integer get_positive_roots( T_real r[] ) const;

    //! Alias of `get_positive_roots()`.
    integer getPositiveRoots( T_real r[] ) const { return get_positive_roots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2 or 3
    //!
    integer get_negative_roots( T_real r[] ) const;

    //! Alias of `get_negative_roots()`.
    integer getNegativeRoots( T_real r[] ) const { return get_negative_roots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_range()`.
    integer getRootsInRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_open_range()`.
    integer getRootsInOpenRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_open_range( a, b, r ); }

    //!
    //! Return the first stored real component.
    //!
    //! When `complex_root()` is true, this is the common real part of the
    //! conjugate pair returned by `root0()` and `root1()`.
    //!
    T_real real_root0() const { return m_r0; }

    //! Return the second stored real component.
    //!
    //! When `complex_root()` is true, this is the positive imaginary part of
    //! the conjugate pair. Otherwise it is the second real root.
    T_real real_root1() const { return m_r1; }

    //! Return the third stored real root.
    T_real real_root2() const { return m_r2; }

    //!
    //! First complex or real root.
    //!
    T_complex root0() const { return m_cplx ? T_complex( m_r0, m_r1 ) : T_complex( m_r0, 0 ); }

    //!
    //! Second complex or real root.
    //!
    T_complex root1() const { return m_cplx ? T_complex( m_r0, -m_r1 ) : T_complex( m_r1, 0 ); }

    //!
    //! Third complex or real root.
    //!
    T_complex root2() const { return T_complex( m_r2, 0 ); }

    //!
    //! First complex or real root.
    //!
    void get_root0( T_real & re, T_real & im ) const
    {
      if ( m_cplx )
      {
        re = m_r0;
        im = m_r1;
      }
      else
      {
        re = m_r0;
        im = 0;
      }
    }

    //! Alias of `get_root0()`.
    void getRoot0( T_real & re, T_real & im ) const { get_root0( re, im ); }

    //!
    //! First complex or real root.
    //!
    void get_root0( T_complex & r ) const { r = m_cplx ? T_complex( m_r0, m_r1 ) : T_complex( m_r0, 0 ); }

    //! Alias of `get_root0()`.
    void getRoot0( T_complex & r ) const { get_root0( r ); }

    //!
    //! Second complex or real root.
    //!
    void get_root1( T_real & re, T_real & im ) const
    {
      if ( m_cplx )
      {
        re = m_r0;
        im = -m_r1;
      }
      else
      {
        re = m_r1;
        im = 0;
      }
    }

    //! Alias of `get_root1()`.
    void getRoot1( T_real & re, T_real & im ) const { get_root1( re, im ); }

    //!
    //! Second complex or real root.
    //!
    void get_root1( T_complex & r ) const { r = m_cplx ? T_complex( m_r0, -m_r1 ) : T_complex( m_r1, 0 ); }

    //! Alias of `get_root1()`.
    void getRoot1( T_complex & r ) const { get_root1( r ); }

    //!
    //! Third complex or real root.
    //!
    void get_root2( T_real & re, T_real & im ) const
    {
      re = m_r2;
      im = 0;
    }

    //! Alias of `get_root2()`.
    void getRoot2( T_real & re, T_real & im ) const { get_root2( re, im ); }

    //!
    //! Third complex or real root.
    //!
    void get_root2( T_complex & r ) const { r = T_complex( m_r2, 0 ); }

    //! Alias of `get_root2()`.
    void getRoot2( T_complex & r ) const { get_root2( r ); }

    //! Return the \p i-th root as a complex value.
    //!
    //! \param[in] i root index in `[0,2]`
    //! \return      requested root, or zero for an invalid index
    T_complex root( integer const i ) const
    {
      switch ( i )
      {
        case 0: return root0();
        case 1: return root1();
        case 2: return root2();
      }
      return 0;
    }

    //! Store the \p i-th root in split real/imaginary form.
    //!
    //! \param[in]  i   root index in `[0,2]`
    //! \param[out] re  real part of the selected root
    //! \param[out] im  imaginary part of the selected root
    void get_root( integer const i, T_real & re, T_real & im ) const
    {
      switch ( i )
      {
        case 0: return get_root0( re, im );
        case 1: return get_root1( re, im );
        case 2: return get_root2( re, im );
      }
    }

    //!
    //! Evaluate the cubic polynomial at a real point.
    //!
    //! \param[in] x  value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    //!
    T_real eval( T_real const & x ) const { return eval_poly<T_real>( m_ABCD, 3, x ); }

    //! Evaluate the cubic polynomial at a complex point.
    //!
    //! \param[in] x  complex value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    T_complex eval( T_complex const & x ) const { return eval_poly_complex<T_real, T_complex>( m_ABCD, 3, x ); }

    //! Evaluate the cubic polynomial and its derivative at a real point.
    //!
    //! \param[in]  x   value where compute \f$ p(x) \f$
    //! \param[out] p   value \f$ p(x) \f$
    //! \param[out] dp  value \f$ p'(x) \f$
    void eval( T_real const & x, T_real & p, T_real & dp ) const { eval_poly_Dpoly<T_real>( m_ABCD, 3, x, p, dp ); }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void info( ostream_type & s ) const;

    //!
    //! Check tolerance and quality of the computed roots.
    //!
    [[nodiscard]] bool check( ostream_type & s ) const;
  };

  using Cubic   = CubicT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using CubicHQ = CubicT<quad_real, quad_complex>;
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*\
   |    ___                   _   _
   |   / _ \ _   _  __ _ _ __| |_(_) ___
   |  | | | | | | |/ _` | '__| __| |/ __|
   |  | |_| | |_| | (_| | |  | |_| | (__
   |   \__\_\\__,_|\__,_|_|   \__|_|\___|
   |
   |  A * x^4 + B * x^3 + C * x^2 + D * x + E
  \*/
  //! Quartic polynomial class
  //!
  //! **Constructor**
  //!
  //! ```{cpp}
  //!   double a = 1;
  //!   double b = 2;
  //!   double c = 3;
  //!   double d = 3;
  //!   double e = 3;
  //!   Quartic p(a,b,c,d,e); // build an solve a x^4 + b x^3 + c x^2 + d x + e = 0
  //!
  //!   Quartic p;
  //!   p.setup(a,b,c,d,e); // build an solve a x^4 + b x^3 + c x^2 + d x + e = 0
  //! ```
  //!
  //! **Get kind of solution**
  //!
  //! ```{cpp}
  //!   int nroots = p.num_roots();
  //!   int nroots = p.num_real_roots();
  //!   int nroots = p.num_complex_root();
  //! ```
  //!
  //! **Get real roots**
  //!
  //! ```{cpp}
  //!   double r_min = 0;
  //!   double r_max = 2;
  //!   double r[4];
  //!   int nroots;
  //!   nroots = p.getRealRoots( r );
  //!   nroots = p.getPositiveRoots( r );
  //!   nroots = p.getNegativeRoots( r );
  //!   nroots = p.getRootsInRange( r_min, r_max, r );
  //!   nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //! ```
  //!
  //! **Get roots**
  //!
  //! ```{cpp}
  //!   real_type    r0 = p.real_root0();
  //!   real_type    r1 = p.real_root1();
  //!   real_type    r2 = p.real_root2();
  //!   real_type    r3 = p.real_root3();
  //!   real_complex r0 = p.root0();
  //!   real_complex r1 = p.root1();
  //!   real_complex r2 = p.root2();
  //!   real_complex r3 = p.root3();
  //!
  //!   real_complex r;
  //!   real_type    re, im;
  //!   p.getRoot0( re, im );
  //!   p.getRoot0( r );
  //!   p.getRoot1( re, im );
  //!   p.getRoot1( r );
  //!   p.getRoot2( re, im );
  //!   p.getRoot2( r );
  //!   p.getRoot3( re, im );
  //!   p.getRoot3( r );
  //! ```
  //!
  //! **Evaluate polynomial**
  //!
  //! ```{cpp}
  //!   {double or complex} v, x;
  //!   v = p.eval( x );
  //!
  //!   p.eval( x, p, dp );
  //! ```
  //!
  //! **Information**
  //!
  //! ```{cpp}
  //!   p.info( cout );
  //!   bool ok = p.check( cout );
  //! ```
  //!
  template <typename T_real, typename T_complex> class QuarticT
  {
    T_real  m_ABCDE[5]{ 0, 0, 0, 0, 0 };
    T_real  m_r0    = 0;
    T_real  m_r1    = 0;
    T_real  m_r2    = 0;
    T_real  m_r3    = 0;
    integer m_iter  = 0;
    integer m_nreal = 0;
    integer m_ncplx = 0;

    void find_roots();

    bool cplx0() const { return m_ncplx > 0; }
    bool cplx1() const { return m_ncplx > 0; }
    bool cplx2() const { return m_ncplx > 2; }
    bool cplx3() const { return m_ncplx > 2; }

  public:
    using value_type   = T_real;
    using complex_type = T_complex;

    //! Build an empty quartic solver instance.
    QuarticT() = default;

    //! Build the instance and immediately compute the quartic roots.
    //!
    //! \param[in] a coefficient of \f$ x^4 \f$
    //! \param[in] b coefficient of \f$ x^3 \f$
    //! \param[in] c coefficient of \f$ x^2 \f$
    //! \param[in] d coefficient of \f$ x \f$
    //! \param[in] e coefficient of \f$ x^0 \f$
    QuarticT( T_real const & a, T_real const & b, T_real const & c, T_real const & d, T_real const & e )
    {
      m_ABCDE[0] = a;
      m_ABCDE[1] = b;
      m_ABCDE[2] = c;
      m_ABCDE[3] = d;
      m_ABCDE[4] = e;
      // find roots only on finite values
      root_assert(
        isfinite( a ) && isfinite( b ) && isfinite( c ) && isfinite( d ) && isfinite( e ),
        "QuarticT( a={}, b={}, c={}, d={}, e={} ) arguments must be finite!",
        a, b, c, d, e
      );
      find_roots();
    }

    //!
    //! Compute the roots of quartic polynomial
    //! \f$ a x^4 + b x^3 + c x^2 + d x + e \f$
    //!
    //! \param[in] a coefficient of \f$ x^4 \f$
    //! \param[in] b coefficient of \f$ x^3 \f$
    //! \param[in] c coefficient of \f$ x^2  \f$
    //! \param[in] d coefficient of \f$ x   \f$
    //! \param[in] e coefficient of \f$ x^0 \f$
    //!
    void setup( T_real const & a, T_real const & b, T_real const & c, T_real const & d, T_real const & e )
    {
      m_ABCDE[0] = a;
      m_ABCDE[1] = b;
      m_ABCDE[2] = c;
      m_ABCDE[3] = d;
      m_ABCDE[4] = e;
      m_iter     = 0;
      m_nreal    = 0;
      m_ncplx    = 0;
      // find roots only on finite values
      root_assert(
        isfinite( a ) && isfinite( b ) && isfinite( c ) && isfinite( d ) && isfinite( e ),
        "QuarticT::setup( a={}, b={}, c={}, d={}, e={} ) arguments must be finite!",
        a, b, c, d, e
      );
      find_roots();
    }

    //!
    //! Number of roots found.
    //!
    integer num_roots() const { return m_nreal + m_ncplx; }

    //! Alias of `num_roots()`.
    integer numRoots() const { return m_nreal + m_ncplx; }

    //!
    //! Number of real roots.
    //!
    integer num_real_roots() const { return m_nreal; }

    //! Alias of `num_real_roots()`.
    integer numRealRoots() const { return m_nreal; }

    //!
    //! Number of complex roots.
    //!
    integer num_complex_roots() const { return m_ncplx; }

    //! Alias of `num_complex_roots()`.
    integer numComplexRoots() const { return m_ncplx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2, 3 or 4
    //!
    integer get_real_roots( T_real r[] ) const;

    //! Alias of `get_real_roots()`.
    integer getRealRoots( T_real r[] ) const { return get_real_roots( r ); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2, 3 or 4
    //!
    integer get_positive_roots( T_real r[] ) const;

    //! Alias of `get_positive_roots()`.
    integer getPositiveRoots( T_real r[] ) const { return get_positive_roots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2, 3 or 4
    //!
    integer get_negative_roots( T_real r[] ) const;

    //! Alias of `get_negative_roots()`.
    integer getNegativeRoots( T_real r[] ) const { return get_negative_roots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_range()`.
    integer getRootsInRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const;

    //! Alias of `get_roots_in_open_range()`.
    integer getRootsInOpenRange( T_real const & a, T_real const & b, T_real r[] ) const
    { return get_roots_in_open_range( a, b, r ); }

    //!
    //! Return the first stored real component.
    //!
    //! When `num_complex_roots() > 0`, this is the real part of the first
    //! conjugate pair returned by `root0()` and `root1()`.
    //!
    T_real real_root0() const { return m_r0; }

    //! Return the second stored real component.
    //!
    //! When `num_complex_roots() > 0`, this is the positive imaginary part of
    //! the first conjugate pair. Otherwise it is the second real root.
    T_real real_root1() const { return m_r1; }

    //! Return the third stored real component.
    //!
    //! When `num_complex_roots() > 2`, this is the real part of the second
    //! conjugate pair returned by `root2()` and `root3()`.
    T_real real_root2() const { return m_r2; }

    //! Return the fourth stored real component.
    //!
    //! When `num_complex_roots() > 2`, this is the positive imaginary part of
    //! the second conjugate pair. Otherwise it is the fourth real root.
    T_real real_root3() const { return m_r3; }

    //!
    //! Return the first root as a complex value.
    //!
    T_complex root0() const { return cplx0() ? T_complex( m_r0, m_r1 ) : T_complex( m_r0, 0 ); }

    //!
    //! Return the second root as a complex value.
    //!
    T_complex root1() const { return cplx1() ? T_complex( m_r0, -m_r1 ) : T_complex( m_r1, 0 ); }

    //!
    //! Return the third root as a complex value.
    //!
    T_complex root2() const { return cplx2() ? T_complex( m_r2, m_r3 ) : T_complex( m_r2, 0 ); }

    //!
    //! Return the fourth root as a complex value.
    //!
    T_complex root3() const { return cplx3() ? T_complex( m_r2, -m_r3 ) : T_complex( m_r3, 0 ); }

    //!
    //! First real or complex root.
    //!
    void get_root0( T_real & re, T_real & im ) const
    {
      if ( cplx0() )
      {
        re = m_r0;
        im = m_r1;
      }
      else
      {
        re = m_r0;
        im = 0;
      }
    }

    //! Alias of `get_root0()`.
    void getRoot0( T_real & re, T_real & im ) const { get_root0( re, im ); }

    //!
    //! First real or complex root.
    //!
    void get_root0( T_complex & r ) const
    {
      if ( cplx0() )
        r = T_complex( m_r0, m_r1 );
      else
        r = T_complex( m_r0, 0 );
    }

    //! Alias of `get_root0()`.
    void getRoot0( T_complex & r ) const { get_root0( r ); }

    //!
    //! Second real or complex root.
    //!
    void get_root1( T_real & re, T_real & im ) const
    {
      if ( cplx1() )
      {
        re = m_r0;
        im = -m_r1;
      }
      else
      {
        re = m_r1;
        im = 0;
      }
    }

    //! Alias of `get_root1()`.
    void getRoot1( T_real & re, T_real & im ) const { get_root1( re, im ); }

    //!
    //! Second real or complex root.
    //!
    void get_root1( T_complex & r ) const
    {
      if ( cplx1() )
        r = T_complex( m_r0, -m_r1 );
      else
        r = T_complex( m_r1, 0 );
    }

    //! Alias of `get_root1()`.
    void getRoot1( T_complex & r ) const { get_root1( r ); }

    //!
    //! Third real or complex root.
    //!
    void get_root2( T_real & re, T_real & im ) const
    {
      if ( cplx2() )
      {
        re = m_r2;
        im = m_r3;
      }
      else
      {
        re = m_r2;
        im = 0;
      }
    }

    //! Alias of `get_root2()`.
    void getRoot2( T_real & re, T_real & im ) const { get_root2( re, im ); }

    //!
    //! Third real or complex root.
    //!
    void get_root2( T_complex & r ) const
    {
      if ( cplx2() )
        r = T_complex( m_r2, m_r3 );
      else
        r = T_complex( m_r2, 0 );
    }

    //! Alias of `get_root2()`.
    void getRoot2( T_complex & r ) const { get_root2( r ); }

    //!
    //! Fourth real or complex root.
    //!
    void get_root3( T_real & re, T_real & im ) const
    {
      if ( cplx3() )
      {
        re = m_r2;
        im = -m_r3;
      }
      else
      {
        re = m_r3;
        im = 0;
      }
    }

    //! Alias of `get_root3()`.
    void getRoot3( T_real & re, T_real & im ) const { get_root3( re, im ); }

    //!
    //! Fourth real or complex root.
    //!
    void get_root3( T_complex & r ) const
    {
      if ( cplx3() )
        r = T_complex( m_r2, -m_r3 );
      else
        r = T_complex( m_r3, 0 );
    }

    //! Alias of `get_root3()`.
    void getRoot3( T_complex & r ) const { get_root3( r ); }

    //! Return the \p i-th root as a complex value.
    //!
    //! \param[in] i root index in `[0,3]`
    //! \return      requested root, or zero for an invalid index
    T_complex root( integer const i ) const
    {
      switch ( i )
      {
        case 0: return root0();
        case 1: return root1();
        case 2: return root2();
        case 3: return root3();
      }
      return 0;
    }

    //! Store the \p i-th root in split real/imaginary form.
    //!
    //! \param[in]  i   root index in `[0,3]`
    //! \param[out] re  real part of the selected root
    //! \param[out] im  imaginary part of the selected root
    void get_root( integer const i, T_real & re, T_real & im ) const
    {
      switch ( i )
      {
        case 0: return get_root0( re, im );
        case 1: return get_root1( re, im );
        case 2: return get_root2( re, im );
        case 3: return get_root3( re, im );
      }
    }

    //!
    //! Evaluate the quartic polynomial at a real point.
    //!
    //! \param[in] x  value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    //!
    T_real eval( T_real const & x ) const { return eval_poly<T_real>( m_ABCDE, 4, x ); }

    //! Evaluate the quartic polynomial at a complex point.
    //!
    //! \param[in] x  complex value where compute \f$ p(x) \f$
    //! \return       the value \f$ p(x) \f$
    T_complex eval( T_complex const & x ) const { return eval_poly_complex<T_real, T_complex>( m_ABCDE, 4, x ); }

    //! Evaluate the quartic polynomial and its derivative at a real point.
    //!
    //! \param[in]  x   value where compute \f$ p(x) \f$
    //! \param[out] p   value \f$ p(x) \f$
    //! \param[out] dp  value \f$ p'(x) \f$
    void eval( T_real const & x, T_real & p, T_real & dp ) const { eval_poly_Dpoly<T_real>( m_ABCDE, 4, x, p, dp ); }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void info( ostream_type & s ) const;

    //!
    //! Check tolerance and quality of the computed roots.
    //!
    [[nodiscard]] bool check( ostream_type & s ) const;
  };

  using Quartic   = QuarticT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using QuarticHQ = QuarticT<quad_real, quad_complex>;
#endif

  /*\
   |   _   _ _   _ _
   |  | | | | |_(_) |___
   |  | | | | __| | / __|
   |  | |_| | |_| | \__ \
   |   \___/ \__|_|_|___/
  \*/

  //! Evaluate a monic cubic polynomial.
  //!
  //! \param[in] x  evaluation point
  //! \param[in] a  coefficient of \f$ x^2 \f$
  //! \param[in] b  coefficient of \f$ x \f$
  //! \param[in] c  constant coefficient
  //! \return       value of \f$ x^3 + a x^2 + b x + c \f$
  template <typename T_real>
  inline T_real evalMonicCubic( T_real const & x, T_real const & a, T_real const & b, T_real const & c )
  {
    T_real p;
    p = x + a;
    p = p * x + b;
    p = p * x + c;
    return p;
  }

  //! Evaluate a monic cubic polynomial and its first derivative.
  template <typename T_real> inline void evalMonicCubic(
    T_real const & x,
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real &       p,
    T_real &       dp )
  {
    p  = x + a;
    dp = x + p;
    p  = p * x + b;
    dp = dp * x + p;
    p  = p * x + c;
  }

  //! Evaluate a monic cubic polynomial and its first two derivatives.
  template <typename T_real> inline void evalMonicCubic(
    T_real const & x,
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real &       p,
    T_real &       dp,
    T_real &       ddp )
  {
    p   = x + a;
    dp  = x + p;      // 2*x + a
    p   = p * x + b;  // x^2 + a * x + b
    ddp = 2 * ( x + dp );
    dp  = dp * x + p;
    p   = p * x + c;
  }

  //! Evaluate a monic quartic polynomial.
  template <typename T_real> inline T_real evalMonicQuartic(
    T_real const & x,
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real const & d )
  {
    T_real p;
    p = x + a;      // x + a
    p = p * x + b;  // x^2+ a*x + b
    p = p * x + c;  // x^3+ a*x^2 + b*x + c
    p = p * x + d;  // x^4+ a*x^3 + b*x^2 + c*x + d
    return p;
  }

  //! Evaluate a monic quartic polynomial and its first derivative.
  template <typename T_real> inline void evalMonicQuartic(
    T_real const & x,
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real const & d,
    T_real &       p,
    T_real &       dp )
  {
    p  = x + a;       // x + a
    dp = x + p;       // 2*x + a
    p  = p * x + b;   // x^2+ a*x + b
    dp = dp * x + p;  // 3*x^2 + 2*a*x + b
    p  = p * x + c;   // x^3+ a*x^2 + b*x + c
    dp = dp * x + p;  // 4*x^3 + 3*a*x^2 + 2*b*x + c
    p  = p * x + d;   // x^4+ a*x^3 + b*x^2 + c*x + d
  }

  //! Evaluate a monic quartic polynomial and its first two derivatives.
  template <typename T_real> inline void evalMonicQuartic(
    T_real const & x,
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real const & d,
    T_real &       p,
    T_real &       dp,
    T_real &       ddp )
  {
    // p_{n+1}(x)   = x * p_{n}(x) + b_{n}
    // p'_{n+1}(x)  = x * p'_{n}(x) + p_{n}(x)
    // p''_{n+1}(x) = x * p''_{n}(x) + 2*p'_{n}(x)
    // ddp = 0;
    // dp  = 1;
    p = x + a;  // x + a

    ddp = 2;
    dp  = x + p;
    p   = p * x + b;

    ddp = ddp * x + 2 * dp;
    dp  = dp * x + p;
    p   = p * x + c;

    ddp = ddp * x + 2 * dp;
    dp  = dp * x + p;
    p   = p * x + d;
  }

}  // namespace PolynomialRoots

#endif
