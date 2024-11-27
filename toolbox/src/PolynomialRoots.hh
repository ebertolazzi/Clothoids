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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef POLYNOMIAL_ROOTS_HH
#define POLYNOMIAL_ROOTS_HH

#include <cmath>
#include <cfloat>
#include <complex>
#include <iostream>

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
namespace PolynomialRoots {

  using real_type    = double;
  using integer      = int;
  using complex_type = std::complex<real_type>;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  //! check if cloating point number `x` is zero
  static
  inline
  bool
  isZero( real_type x )
  { return FP_ZERO == std::fpclassify(x); }

  //! evaluate real polynomial
  real_type
  evalPoly(
    real_type const op[],
    integer         Degree,
    real_type       x
  );

  void
  evalPolyDPoly(
    real_type const op[],
    integer         Degree,
    real_type       x,
    real_type     & p,
    real_type     & dp
  );

  bool
  NewtonStep(
    real_type const op[],
    integer         Degree,
    real_type     & x
  );

  //! evaluate real polynomial with complex value
  complex_type
  evalPolyC(
    real_type const op[],
    integer         Degree,
    complex_type    x
  );

  #endif

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
  int
  roots(
    real_type const * op,
    integer           Degree,
    real_type       * zeror,
    real_type       * zeroi
  );

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
  //!   double r0 = p.real_root0();
  //!   double r1 = p.real_root1();
  //!   complex_type r0 = p.root0();
  //!   complex_type r1 = p.root1();
  //!
  //!   complex_type r;
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
  class Quadratic {
    real_type m_ABC[3]{0,0,0};
    real_type m_r0{0}, m_r1{0};
    integer   m_nrts{0};
    bool      m_cplx{false};
    bool      m_dblx{false};

    void find_roots();

  public:

    Quadratic() {}

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
    Quadratic( real_type a, real_type b, real_type c ) {
      using std::isfinite;
      real_type & A{ m_ABC[0] };
      real_type & B{ m_ABC[1] };
      real_type & C{ m_ABC[2] };
      A = a; B = b; C = c;
      // find roots only on finite values
      if ( isfinite(a) && isfinite(b) && isfinite(c) ) {
        find_roots();
      }
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
    void
    setup( real_type a, real_type b, real_type c ) {
      real_type & A = m_ABC[0];
      real_type & B = m_ABC[1];
      real_type & C = m_ABC[2];
      A = a; B = b; C = c;
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

    //! alias of `num_roots`
    integer numRoots() const { return m_nrts; }

    //!
    //! true if roots are complex conjugated
    //!
    bool complex_root() const { return m_cplx; }

    //! alias of `complex_root`
    bool complexRoot() const { return m_cplx; }

    //!
    //! true if \f$ p(x) = a x^2 + b x + c = (x-r)^2 \f$
    //!
    bool double_root() const { return m_dblx; }

    //! alias of `double_root`
    bool doubleRoot() const { return m_dblx; }

    //!
    //! Get the real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1 or 2
    //!
    integer get_real_roots( real_type r[] ) const;

    //! alias of `get_real_roots`
    integer getRealRoots( real_type r[] ) const { return get_real_roots(r); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1 or 2
    //!
    integer get_positive_roots( real_type r[] ) const;

    //! alias of `get_positive_roots`
    integer getPositiveRoots( real_type r[] ) const { return get_positive_roots(r); }

    //!
    //! Get negative real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1 or 2
    //!
    integer get_negative_roots( real_type r[] ) const;

    //! alias of `get_negative_roots`
    integer getNegativeRoots( real_type r[] ) const { return get_negative_roots(r); }

    //!
    //! Get real roots in a closed range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range \f$ [a,b] \f$
    //!
    integer get_roots_in_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_range`
    integer
    getRootsInRange( real_type a, real_type b, real_type r[] ) const
    { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range \f$ (a,b) \f$
    //!
    integer get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_open_range`
    integer
    getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const
    { return get_roots_in_open_range( a, b, r ); }

    //!
    //! The first real root
    //!
    real_type real_root0() const { return m_r0; }

    //!
    //! The second real root
    //!
    real_type real_root1() const { return m_r1; }

    //!
    //! The first complex root
    //!
    complex_type root0() const { return m_cplx ? complex_type(m_r0,m_r1) : complex_type(m_r0,0); }

    //!
    //! The second complex root
    //!
    complex_type root1() const { return m_cplx ? complex_type(m_r0,-m_r1) : complex_type(m_r1,0); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] re the first complex root, real part
    //! \param[out] im the first complex root, imaginary part
    //!
    void
    get_root0( real_type & re, real_type & im ) const {
      if ( m_cplx ) { re = m_r0; im = m_r1; }
      else          { re = m_r0; im = 0;    }
    }

    //! alias of `get_root0`
    void
    getRoot0( real_type & re, real_type & im  ) const
    { return get_root0( re, im ); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] r the first complex root
    //!
    void
    get_root0( complex_type & r ) const
    { r = m_cplx ? complex_type(m_r0,m_r1) : complex_type(m_r0,0); }

    //! alias of `get_root0`
    void getRoot0( complex_type & r ) const { return get_root0( r ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] re the second complex root, real part
    //! \param[out] im the second complex root, imaginary part
    //!
    void
    get_root1( real_type & re, real_type & im ) const {
      if ( m_cplx ) { re = m_r0; im = -m_r1; }
      else          { re = m_r1; im = 0;     }
    }

    //! alias of `get_root1`
    void
    getRoot1( real_type & re, real_type & im  ) const
    { return get_root1( re, im ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] r the second complex root
    //!
    void
    get_root1( complex_type & r ) const
    { r = m_cplx ? complex_type(m_r0,-m_r1) : complex_type(m_r1,0); }

    //! alias of `get_root1`
    void getRoot1( complex_type & r ) const { return get_root1( r ); }

    //!
    //! Evaluate the quadratic polynomial
    //!
    //! \param x  value where compute \f$ p(x) \f$
    //! \return   the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type x ) const
    { return evalPoly( m_ABC, 2, x ); }

    //!
    //! Evalute the quadratic polynomial
    //!
    //! \param[in] x value where compute \f$ p(x)=ax^2+bx+c \f$
    //! \return      the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type x ) const
    { return evalPolyC( m_ABC, 2, x ); }

    //!
    //! Evaluate the polynomial with its derivative
    //!
    //! \param[in]  x   value where compute \f$ p(x) \f$
    //! \param[out] p   value \f$ p(x) \f$
    //! \param[out] dp  value \f$ p'(x) \f$
    //!
    void
    eval( real_type x, real_type & p, real_type & dp ) const {
      evalPolyDPoly( m_ABC, 2, x, p, dp );
    }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void
    info( ostream_type & s ) const;

    //!
    //! Check tolerance and quality of the computed roots
    //!
    bool
    check( ostream_type & s ) const;
  };

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
  //!   double r0 = p.real_root0();
  //!   double r1 = p.real_root1();
  //!   double r2 = p.real_root2();
  //!   complex_type r0 = p.root0();
  //!   complex_type r1 = p.root1();
  //!   complex_type r2 = p.root2();
  //!
  //!   complex_type r;
  //!   double re, im;
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
  class Cubic {
    real_type m_ABCD[4]{0,0,0,0};
    real_type m_r0{0}, m_r1{0}, m_r2{0};
    integer   m_nrts{0};
    integer   m_iter{0};
    bool      m_cplx{false}; // complex root
    bool      m_dblx{false}; // double root
    bool      m_trpx{false}; // triple root

    void find_roots();

  public:

    //!
    //! Build an empty instance of Cubic polynomial solver
    //!
    Cubic() {}

    //!
    //! Build the instance of the calss and compute the roots of cubic polynomial
    //! \f$ a x^3 + b x^2 + c x + d \f$
    //!
    //! \param[in] a coefficient of \f$ x^3 \f$
    //! \param[in] b coefficient of \f$ x^2 \f$
    //! \param[in] c coefficient of \f$ x   \f$
    //! \param[in] d coefficient of \f$ x^0 \f$
    //!
    Cubic( real_type a, real_type b, real_type c, real_type d ) {
      using std::isfinite;
      real_type & A = m_ABCD[0];
      real_type & B = m_ABCD[1];
      real_type & C = m_ABCD[2];
      real_type & D = m_ABCD[3];
      A = a; B = b; C = c; D = d;
      // find roots only on finite values
      if ( isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d) ) {
        find_roots();
      }
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
    void
    setup( real_type a, real_type b, real_type c, real_type d ) {
      using std::isfinite;
      real_type & A{m_ABCD[0]};
      real_type & B{m_ABCD[1]};
      real_type & C{m_ABCD[2]};
      real_type & D{m_ABCD[3]};
      m_nrts = 0;
      m_iter = 0;
      m_cplx = false; // complex root
      m_dblx = false; // double root
      m_trpx = false; // triple root
      A = a; B = b; C = c; D = d;
      if ( isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d) ) {
        find_roots();
      }
    }

    //!
    //! Number of found roots.
    //!
    integer num_roots() const { return m_nrts; }

    //! alias of `num_roots`
    integer numRoots() const { return m_nrts; }

    //!
    //! Has complex roots?
    //!
    bool complex_root() const { return m_cplx; }

    //! alias of `complex_root`
    bool complexRoot() const { return m_cplx; }

    //!
    //! Has a double root?
    //!
    bool double_root() const { return m_dblx; }

    //! alias of `double_root`
    bool doubleRoot() const { return m_dblx; }

    //!
    //! Has a triple root?
    //!
    bool triple_root() const { return m_trpx; }

    //! alias of `triple_root`
    bool tripleRoot() const { return m_trpx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2 or 3
    //!
    integer get_real_roots( real_type r[] ) const;

    //! alias of `get_real_roots`
    integer
    getRealRoots( real_type r[] ) const
    { return get_real_roots( r ); }

    //!
    //! Get positive real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2 or 3
    //!
    integer get_positive_roots( real_type r[] ) const;

    //! alias of `get_positive_roots`
    integer
    getPositiveRoots( real_type r[] ) const
    { return get_positive_roots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2 or 3
    //!
    integer get_negative_roots( real_type r[] ) const;

    //! alias of `get_negative_roots`
    integer
    getNegativeRoots( real_type r[] ) const
    { return get_negative_roots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer get_roots_in_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_range`
    integer
    getRootsInRange( real_type a, real_type b, real_type r[] ) const
    { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_open_range`
    integer
    getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const
    { return get_roots_in_open_range( a, b, r ); }

    //!
    //! First real root.
    //!
    real_type real_root0() const { return m_r0; }

    //!
    //! Second real root.
    //!
    real_type real_root1() const { return m_r1; }

    //!
    //! Third real root.
    //!
    real_type real_root2() const { return m_r2; }

    //!
    //! First complex or real root.
    //!
    complex_type
    root0() const
    { return m_cplx ? complex_type(m_r0,m_r1) : complex_type(m_r0,0); }

    //!
    //! Second complex or real root.
    //!
    complex_type
    root1() const
    { return m_cplx ? complex_type(m_r0,-m_r1) : complex_type(m_r1,0); }

    //!
    //! Third complex or real root.
    //!
    complex_type
    root2() const
    { return complex_type(m_r2,0); }

    //!
    //! First complex or real root.
    //!
    void
    get_root0( real_type & re, real_type & im ) const {
      if ( m_cplx ) { re = m_r0; im = m_r1; }
      else          { re = m_r0; im = 0;    }
    }

    //! alias of `get_root0`
    void getRoot0( real_type & re, real_type & im ) const { get_root0( re, im ); }

    //!
    //! First complex or real root.
    //!
    void
    get_root0( complex_type & r ) const
    { r = m_cplx ? complex_type(m_r0,m_r1) : complex_type(m_r0,0); }

    //! alias of `get_root0`
    void getRoot0( complex_type & r ) const { get_root0( r ); }

    //!
    //! Second complex or real root.
    //!
    void
    get_root1( real_type & re, real_type & im ) const {
      if ( m_cplx ) { re = m_r0; im = -m_r1; }
      else          { re = m_r1; im = 0;     }
    }

    //! alias of `get_root1`
    void getRoot1( real_type & re, real_type & im ) const { get_root1( re, im ); }

    //!
    //! Second complex or real root.
    //!
    void
    get_root1( complex_type & r ) const
    { r = m_cplx ? complex_type(m_r0,-m_r1) : complex_type(m_r1,0); }

    //! alias of `get_root1`
    void getRoot1( complex_type & r ) const { get_root1( r ); }

    //!
    //! Third complex or real root.
    //!
    void
    get_root2( real_type & re, real_type & im ) const
    { re = m_r2; im = 0; }

    //! alias of `get_root2`
    void getRoot2( real_type & re, real_type & im ) const { get_root2( re, im ); }

    //!
    //! Third complex or real root.
    //!
    void
    get_root2( complex_type & r ) const
    { r = complex_type(m_r2,0); }

    //! alias of `get_root2`
    void getRoot2( complex_type & r ) const { get_root2( r ); }

    //!
    //! Evalute the cubic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$
    //! \return the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type x ) const
    { return evalPoly( m_ABCD, 3, x ); }

    //!
    //! Evalute the cubic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type x ) const
    { return evalPolyC( m_ABCD, 3, x ); }

    //!
    //! Evalute the polynomial with its derivative.
    //!
    //! \param x   value where compute \f$ p(x) \f$
    //! \param p   value \f$ p(x) \f$
    //! \param dp  value \f$ p'(x) \f$
    //!
    void
    eval( real_type x, real_type & p, real_type & dp ) const {
      evalPolyDPoly( m_ABCD, 3, x, p, dp );
    }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void
    info( ostream_type & s ) const;

    //!
    //! Check tolerenace and quality of the computed roots.
    //!
    bool
    check( ostream_type & s ) const;
  };

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
  //!   double r0 = p.real_root0();
  //!   double r1 = p.real_root1();
  //!   double r2 = p.real_root2();
  //!   double r3 = p.real_root3();
  //!   complex_type r0 = p.root0();
  //!   complex_type r1 = p.root1();
  //!   complex_type r2 = p.root2();
  //!   complex_type r3 = p.root3();
  //!
  //!   complex_type r;
  //!   double re, im;
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
  class Quartic {
    real_type m_ABCDE[5]{0,0,0,0,0};
    real_type m_r0{0}, m_r1{0}, m_r2{0},m_r3{0};
    integer   m_iter{0}, m_nreal{0}, m_ncplx{0};

    void find_roots();

    bool cplx0() const { return m_ncplx > 0; }
    bool cplx1() const { return m_ncplx > 0; }
    bool cplx2() const { return m_ncplx > 2; }
    bool cplx3() const { return m_ncplx > 2; }

  public:

    Quartic() {}
    Quartic(
      real_type a,
      real_type b,
      real_type c,
      real_type d,
      real_type e
    ) {
      using std::isfinite;
      real_type & A{m_ABCDE[0]};
      real_type & B{m_ABCDE[1]};
      real_type & C{m_ABCDE[2]};
      real_type & D{m_ABCDE[3]};
      real_type & E{m_ABCDE[4]};
      A = a; B = b; C = c; D = d; E = e;
      // find roots only on finite values
      if ( isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d) && isfinite(e) ) {
        find_roots();
      }
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
    void
    setup(
      real_type a,
      real_type b,
      real_type c,
      real_type d,
      real_type e
    ) {
      using std::isfinite;
      real_type & A{m_ABCDE[0]};
      real_type & B{m_ABCDE[1]};
      real_type & C{m_ABCDE[2]};
      real_type & D{m_ABCDE[3]};
      real_type & E{m_ABCDE[4]};
      A = a; B = b; C = c; D = d; E = e;
      m_iter  = 0;
      m_nreal = 0;
      m_ncplx = 0;
      // find roots only on finite values
      if ( isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d) && isfinite(e) ) {
        find_roots();
      }
    }

    //!
    //! Number of roots found.
    //!
    integer num_roots() const { return m_nreal+m_ncplx; }

    //! alias of `num_roots`
    integer numRoots() const { return m_nreal+m_ncplx; }

    //!
    //! Number of real roots.
    //!
    integer num_real_roots() const { return m_nreal; }

    //! alias of `num_real_roots`
    integer numRealRoots() const { return m_nreal; }

    //!
    //! Number of complex roots
    //!
    integer num_complex_roots() const { return m_ncplx; }

    //! alias of `num_complex_root`
    integer numComplexRoots() const { return m_ncplx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2, 3 or 4
    //!
    integer get_real_roots( real_type r[] ) const;

    //! alias of `get_real_roots`
    integer getRealRoots( real_type r[] ) const { return get_real_roots(r); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2, 3 or 4
    //!
    integer get_positive_roots( real_type r[] ) const;

    //! alias of `get_positive_roots`
    integer getPositiveRoots( real_type r[] ) const { return get_positive_roots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2, 3 or 4
    //!
    integer get_negative_roots( real_type r[] ) const;

    //! alias of `get_negative_roots`
    integer getNegativeRoots( real_type r[] ) const { return get_negative_roots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer get_roots_in_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_range`
    integer getRootsInRange( real_type a, real_type b, real_type r[] ) const { return get_roots_in_range( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const;

    //! alias of `get_roots_in_open_range`
    integer getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const { return get_roots_in_open_range( a, b, r ); }

    //!
    //! First real root.
    //!
    real_type real_root0() const { return m_r0; }

    //!
    //! Second real root.
    //!
    real_type real_root1() const { return m_r1; }

    //!
    //! Third real root.
    //!
    real_type real_root2() const { return m_r2; }

    //!
    //! Fourth real root.
    //!
    real_type real_root3() const { return m_r3; }

    //!
    //! First real or complex root.
    //!
    complex_type
    root0() const
    { return cplx0() ? complex_type(m_r0,m_r1) : complex_type(m_r0,0); }

    //!
    //! Second real or complex root.
    //!
    complex_type
    root1() const
    { return cplx1() ? complex_type(m_r0,-m_r1) : complex_type(m_r1,0); }

    //!
    //! Third real or complex root.
    //!
    complex_type
    root2() const
    { return cplx2() ? complex_type(m_r2,m_r3) : complex_type(m_r2,0); }

    //!
    //! 4th real or complex root.
    //!
    complex_type
    root3() const
    { return cplx3() ? complex_type(m_r2,-m_r3) : complex_type(m_r3,0); }

    //!
    //! First real or complex root.
    //!
    void
    get_root0( real_type & re, real_type & im ) const {
      if ( cplx0() ) { re = m_r0; im = m_r1; }
      else           { re = m_r0; im = 0;    }
    }

    //! alias of `get_root0`
    void getRoot0( real_type & re, real_type & im ) const { get_root0( re, im ); }

    //!
    //! First real or complex root.
    //!
    void
    get_root0( complex_type & r ) const {
      if ( cplx0() ) r = complex_type(m_r0,m_r1);
      else           r = complex_type(m_r0,0);
    }

    //! alias of `get_root0`
    void getRoot0( complex_type & r ) const { get_root0( r ); }

    //!
    //! Second real or complex root.
    //!
    void
    get_root1( real_type & re, real_type & im ) const {
      if ( cplx1() ) { re = m_r0; im = -m_r1; }
      else           { re = m_r1; im = 0;     }
    }

    //! alias of `get_root1`
    void getRoot1( real_type & re, real_type & im ) const { get_root1( re, im ); }

    //!
    //! Second real or complex root.
    //!
    void
    get_root1( complex_type & r ) const {
      if ( cplx1() ) r = complex_type(m_r0,-m_r1);
      else           r = complex_type(m_r1,0);
    }

    //! alias of `get_root1`
    void getRoot1( complex_type & r ) const { get_root1( r ); }

    //!
    //! Third real or complex root.
    //!
    void
    get_root2( real_type & re, real_type & im ) const {
      if ( cplx2() ) { re = m_r2; im = m_r3; }
      else           { re = m_r2; im = 0;    }
    }

    //! alias of `get_root2`
    void getRoot2( real_type & re, real_type & im ) const { get_root2( re, im ); }

    //!
    //! Third real or complex root.
    //!
    void
    get_root2( complex_type & r ) const {
      if ( cplx2() ) r = complex_type(m_r2,m_r3);
      else           r = complex_type(m_r2,0);
    }

    //! alias of `get_root2`
    void getRoot2( complex_type & r ) const { get_root2( r ); }

    //!
    //! 4th real or complex root.
    //!
    void
    get_root3( real_type & re, real_type & im ) const {
      if ( cplx3() ) { re = m_r2; im = -m_r3; }
      else           { re = m_r3; im = 0;     }
    }

    //! alias of `get_root3`
    void getRoot3( real_type & re, real_type & im ) const { get_root3( re, im ); }

    //!
    //! 4th real or complex root.
    //!
    void
    get_root3( complex_type & r ) const {
      if ( cplx3() ) r = complex_type(m_r2,-m_r3);
      else           r = complex_type(m_r3,0);
    }

    //! alias of `get_root3`
    void getRoot3( complex_type & r ) const { get_root3( r ); }

    //!
    //! Evalute the quartic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type x ) const
    { return evalPoly( m_ABCDE, 4, x ); }

    //!
    //! Evalute the quartic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type x ) const
    { return evalPolyC( m_ABCDE, 4, x ); }

    //!
    //! Evalute the polynomial with its derivative.
    //!
    //! \param x   value where compute \f$ p(x) \f$
    //! \param p   value \f$ p(x) \f$
    //! \param dp  value \f$ p'(x) \f$
    //!
    void
    eval( real_type x, real_type & p, real_type & dp ) const {
      evalPolyDPoly( m_ABCDE, 4, x, p, dp );
    }

    //!
    //! Print info of the roots of the polynomial.
    //!
    void
    info( ostream_type & s ) const;

    //!
    //! Check tolerenace and quality of the computed roots.
    //!
    bool
    check( ostream_type & s ) const;

  };

  /*\
   |   _   _ _   _ _
   |  | | | | |_(_) |___
   |  | | | | __| | / __|
   |  | |_| | |_| | \__ \
   |   \___/ \__|_|_|___/
  \*/

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  // x^3 + a*x^2 + b*x + c
  static
  inline
  real_type
  evalMonicCubic(
    real_type x,
    real_type a,
    real_type b,
    real_type c
  ) {
    real_type p;
    p = x + a;
    p = p * x + b;
    p = p * x + c;
    return p;
  }

  static
  inline
  void
  evalMonicCubic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & p,
    real_type & dp
  ) {
    p  = x + a;
    dp = x + p;
    p  = p  * x + b;
    dp = dp * x + p;
    p  = p  * x + c;
  }

  // 3*x^2 + 2*a*x + b
  // 6*x + 2*a
  static
  inline
  void
  evalMonicCubic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & p,
    real_type & dp,
    real_type & ddp
  ) {
    p   = x + a;
    dp  = x + p;      // 2*x + a
    p   = p  * x + b; // x^2 + a * x + b
    ddp = 2*(x + dp);
    dp  = dp * x + p;
    p   = p  * x + c;
  }

  // x^4 + a*x^3 + b*x^2 + c*x + d
  static
  inline
  real_type
  evalMonicQuartic(
    real_type x,
    real_type a,
    real_type b,
    real_type c,
    real_type d
  ) {
    real_type p;
    p = x + a;     // x + a
    p = p * x + b; // x^2+ a*x + b
    p = p * x + c; // x^3+ a*x^2 + b*x + c
    p = p * x + d; // x^4+ a*x^3 + b*x^2 + c*x + d
    return p;
  }

  static
  inline
  void
  evalMonicQuartic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type   d,
    real_type & p,
    real_type & dp
  ) {
    p  = x + a;      // x + a
    dp = x + p;      // 2*x + a
    p  = p  * x + b; // x^2+ a*x + b
    dp = dp * x + p; // 3*x^2 + 2*a*x + b
    p  = p  * x + c; // x^3+ a*x^2 + b*x + c
    dp = dp * x + p; // 4*x^3 + 3*a*x^2 + 2*b*x + c
    p  = p  * x + d; // x^4+ a*x^3 + b*x^2 + c*x + d
  }

  static
  inline
  void
  evalMonicQuartic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type   d,
    real_type & p,
    real_type & dp,
    real_type & ddp
  ) {
    // p_{n+1}(x)   = x * p_{n}(x) + b_{n}
    // p'_{n+1}(x)  = x * p'_{n}(x) + p_{n}(x)
    // p''_{n+1}(x) = x * p''_{n}(x) + 2*p'_{n}(x)
    // ddp = 0;
    // dp  = 1;
    p   = x + a;     // x + a

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

  // x^3 + a*x^2 + b*x + c
  static
  inline
  real_type
  eval_monic_cubic(
    real_type x,
    real_type a,
    real_type b,
    real_type c
  ) {
    return evalMonicCubic( x, a, b, c );
  }

  static
  inline
  void
  eval_monic_cubic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & p,
    real_type & dp
  ) {
    evalMonicCubic( x, a, b, c, p, dp );
  }

  // 3*x^2 + 2*a*x + b
  // 6*x + 2*a
  static
  inline
  void
  eval_monic_cubic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type & p,
    real_type & dp,
    real_type & ddp
  ) {
    evalMonicCubic( x, a, b, c, p, dp, ddp );
  }

  // x^4 + a*x^3 + b*x^2 + c*x + d
  static
  inline
  real_type
  eval_monic_quartic(
    real_type x,
    real_type a,
    real_type b,
    real_type c,
    real_type d
  ) {
    return evalMonicQuartic( x, a, b, c, d );
  }

  static
  inline
  void
  eval_monic_quartic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type   d,
    real_type & p,
    real_type & dp
  ) {
    evalMonicQuartic( x, a, b, c, d, p, dp );
  }

  static
  inline
  void
  eval_monic_quartic(
    real_type   x,
    real_type   a,
    real_type   b,
    real_type   c,
    real_type   d,
    real_type & p,
    real_type & dp,
    real_type & ddp
  ) {
    evalMonicQuartic( x, a, b, c, d, p, dp, ddp );
  }
  #endif
}

#endif
