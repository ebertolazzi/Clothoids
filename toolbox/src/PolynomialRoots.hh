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
//! \rst
//!
//! -  | **N.Flocke**
//!    | Algorithm 954: An Accurate and Efficient Cubic and Quartic
//!    | Equation Solver for Physical Applications
//!    | ACM TOMS, vol 41, n.4, 2015
//!
//! -  | **M.A. Jenkins and J.F.Traub**
//!    | A Three-Stage Algorithm for Real Polynomials Using Quadratic
//!      Iteration
//!    | SIAM Journal on Numerical Analysis
//!    | Vol. 7, No. 4 (Dec., 1970), pp. 545-566
//!
//! \endrst
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
    real_type const      op[],
    integer              Degree,
    complex_type const & x
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
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!   double a = 1;
  //!   double b = 2;
  //!   double c = 3;
  //!   Quadratic q(a,b,c); // build an solve `a x^2 + b x + c = 0`
  //!
  //!   Quadratic q;
  //!   q.setup(a,b,c); // build an solve `a x^2 + b x + c = 0`
  //!
  //! \endrst
  //!
  //! **Get kind of solution**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!   int  nroots            = q.numRoots();
  //!   bool has_complex_roots = q.complexRoots();
  //!   bool has_a_double_root = q.doubleRoot();
  //!
  //! \endrst
  //!
  //! **Get real roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!   double r_min = 0;
  //!   double r_max = 2;
  //!   double r[2];
  //!   int nroots;
  //!   nroots = p.getRealRoots( r );
  //!   nroots = p.getPositiveRoots( r );
  //!   nroots = p.getNegativeRoots( r );
  //!   nroots = p.getRootsInRange( r_min, r_max, r );
  //!   nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //!
  //! \endrst
  //!
  //! **Get roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
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
  //!
  //! \endrst
  //!
  //! **Evaluate polynomial**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!   {double or complex} v, x;
  //!   v = p.eval( x );
  //!
  //!   p.eval( x, p, dp );
  //!
  //! \endrst
  //!
  //! **Information**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!   p.info( cout );
  //!   bool ok = p.check( cout );
  //!
  //! \endrst
  //!
  class Quadratic {
    real_type ABC[3];
    real_type r0, r1;
    integer   nrts;
    bool      cplx;
    bool      dblx;

    void findRoots();

  public:

    Quadratic() : nrts(0), cplx(false), dblx(false) {}

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
    Quadratic( real_type a, real_type b, real_type c )
    : nrts(0), cplx(false), dblx(false) {
      real_type & A = ABC[0];
      real_type & B = ABC[1];
      real_type & C = ABC[2];
      A = a; B = b; C = c;
      findRoots();
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
      real_type & A = ABC[0];
      real_type & B = ABC[1];
      real_type & C = ABC[2];
      A = a; B = b; C = c;
      findRoots();
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
    integer numRoots() const { return nrts; }

    //! alias of `numRoots`
    integer num_roots() const { return nrts; }

    //!
    //! true if roots are complex conjugated
    //!
    bool complexRoots() const { return cplx; }

    //! alias of `complexRoots`
    bool complex_roots() const { return cplx; }

    //!
    //! true if \f$ p(x) = a x^2 + b x + c = (x-r)^2 \f$
    //!
    bool doubleRoot() const { return dblx; }

    //! alias of `doubleRoot`
    bool double_root() const { return dblx; }

    //!
    //! Get the real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1 or 2
    //!
    integer getRealRoots( real_type r[] ) const;

    //! alias of `getRealRoots`
    integer get_real_roots( real_type r[] ) const { return getRealRoots(r); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1 or 2
    //!
    integer getPositiveRoots( real_type r[] ) const;

    //! alias of `getPositiveRoots`
    integer get_positive_roots( real_type r[] ) const { return getPositiveRoots(r); }

    //!
    //! Get negative real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1 or 2
    //!
    integer getNegativeRoots( real_type r[] ) const;

    //! alias of `getNegativeRoots`
    integer get_negative_roots( real_type r[] ) const { return getNegativeRoots(r); }

    //!
    //! Get real roots in a closed range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range \f$ [a,b] \f$
    //!
    integer getRootsInRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInRange`
    integer
    get_roots_in_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInRange( a, b, r ); }

    //!
    //! Get real roots in an open range
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range \f$ (a,b) \f$
    //!
    integer getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInOpenRange`
    integer
    get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInOpenRange( a, b, r ); }

    //!
    //! The first real root
    //!
    real_type real_root0() const { return r0; }

    //!
    //! The second real root
    //!
    real_type real_root1() const { return r1; }

    //!
    //! The first complex root
    //!
    complex_type
    root0() const
    { return cplx ? complex_type(r0,r1) : complex_type(r0,0); }

    //!
    //! The second complex root
    //!
    complex_type
    root1() const
    { return cplx ? complex_type(r0,-r1) : complex_type(r1,0); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] re the first complex root, real part
    //! \param[out] im the first complex root, imaginary part
    //!
    void
    getRoot0( real_type & re, real_type & im ) const {
      if ( cplx ) { re = r0; im = r1; }
      else        { re = r0; im = 0;  }
    }

    //! alias of `getRoot0`
    void
    get_root0( real_type & re, real_type & im  ) const
    { return getRoot0( re, im ); }

    //!
    //! Get the first root (complex or real)
    //!
    //! \param[out] r the first complex root
    //!
    void
    getRoot0( complex_type & r ) const
    { r = cplx ? complex_type(r0,r1) : complex_type(r0,0); }

    //! alias of `getRoot0`
    void get_root0( complex_type & r ) const { return getRoot0( r ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] re the second complex root, real part
    //! \param[out] im the second complex root, imaginary part
    //!
    void
    getRoot1( real_type & re, real_type & im ) const {
      if ( cplx ) { re = r0; im = -r1; }
      else        { re = r1; im = 0;   }
    }

    //! alias of `getRoot1`
    void
    get_root1( real_type & re, real_type & im  ) const
    { return getRoot1( re, im ); }

    //!
    //! Get the second root (complex or real)
    //!
    //! \param[out] r the second complex root
    //!
    void
    getRoot1( complex_type & r ) const
    { r = cplx ? complex_type(r0,-r1) : complex_type(r1,0); }

    //! alias of `getRoot1`
    void get_root1( complex_type & r ) const { return getRoot1( r ); }

    //!
    //! Evaluate the quadratic polynomial
    //!
    //! \param x  value where compute \f$ p(x) \f$
    //! \return   the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type const & x ) const
    { return evalPoly( ABC, 2, x ); }

    //!
    //! Evalute the quadratic polynomial
    //!
    //! \param[in] x value where compute \f$ p(x)=ax^2+bx+c \f$
    //! \return      the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type const & x ) const
    { return evalPolyC( ABC, 2, x ); }

    //!
    //! Evaluate the polynomial with its derivative
    //!
    //! \param[in]  x   value where compute \f$ p(x) \f$
    //! \param[out] p   value \f$ p(x) \f$
    //! \param[out] dp  value \f$ p'(x) \f$
    //!
    void eval( real_type x, real_type & p, real_type & dp ) const;

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
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double a = 1;
  //!      double b = 2;
  //!      double c = 3;
  //!      double d = 3;
  //!      Cubic p(a,b,c,d); // build an solve ``a x^3 + b x^2 + c x + d = 0``
  //!
  //!
  //!      Cubic p;
  //!      p.setup(a,b,c,d); // build an solve ``a x^3 + b x^2 + c x + d = 0``
  //!
  //! \endrst
  //!
  //! **Get kind of solution**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      int  nroots = p.numRoots();
  //!      bool has_complex_roots = p.complexRoots();
  //!      bool has_a_double_root = p.doubleRoot();
  //!      bool has_a_triple_root = p.tripleRoot();
  //!
  //! \endrst
  //!
  //! **Get real roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double r_min = 0;
  //!      double r_max = 2;
  //!      double r[3];
  //!      int nroots;
  //!      nroots = p.getRealRoots( r );
  //!      nroots = p.getPositiveRoots( r );
  //!      nroots = p.getNegativeRoots( r );
  //!      nroots = p.getRootsInRange( r_min, r_max, r );
  //!      nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //!
  //! \endrst
  //!
  //! **Get roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double r0 = p.real_root0();
  //!      double r1 = p.real_root1();
  //!      double r2 = p.real_root2();
  //!      complex_type r0 = p.root0();
  //!      complex_type r1 = p.root1();
  //!      complex_type r2 = p.root2();
  //!
  //!      complex_type r;
  //!      double re, im;
  //!      p.getRoot0( re, im );
  //!      p.getRoot0( r );
  //!      p.getRoot1( re, im );
  //!      p.getRoot1( r );
  //!      p.getRoot2( re, im );
  //!      p.getRoot2( r );
  //!
  //! \endrst
  //!
  //! **Evaluate polynomial**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      {double or complex} v, x;
  //!      v = p.eval( x );
  //!
  //!      p.eval( x, p, dp );
  //!
  //! \endrst
  //!
  //! **Information**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      p.info( cout );
  //!      bool ok = p.check( cout );
  //!
  //! \endrst
  //!
  class Cubic {
    real_type ABCD[4];
    real_type r0, r1, r2;
    integer   nrts, iter;
    bool      cplx; // complex root
    bool      dblx; // double root
    bool      trpx; // triple root

    void findRoots();

  public:

    Cubic() : nrts(0), iter(0), cplx(false), trpx(false) {}
    Cubic( real_type _a, real_type _b, real_type _c, real_type _d )
    : nrts(0), iter(0), cplx(false), trpx(false) {
      real_type & A = ABCD[0];
      real_type & B = ABCD[1];
      real_type & C = ABCD[2];
      real_type & D = ABCD[3];
      A = _a; B = _b; C = _c; D = _d;
      findRoots();
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
      real_type & A = ABCD[0];
      real_type & B = ABCD[1];
      real_type & C = ABCD[2];
      real_type & D = ABCD[3];
      A = a; B = b; C = c; D = d;
      findRoots();
    }

    //!
    //! Number of found roots.
    //!
    integer numRoots() const { return nrts; }

    //! alias of `numRoots`
    integer num_roots() const { return nrts; }

    //!
    //! Has complex roots?
    //!
    bool complexRoots() const { return cplx; }

    //! alias of `complexRoots`
    bool complex_roots() const { return cplx; }

    //!
    //! Has a double root?
    //!
    bool doubleRoot() const { return dblx; }

    //! alias of `doubleRoot`
    bool double_root() const { return dblx; }

    //!
    //! Has a triple root?
    //!
    bool tripleRoot() const { return trpx; }

    //! alias of `tripleRoot`
    bool triple_root() const { return trpx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2 or 3
    //!
    integer getRealRoots( real_type r[] ) const;

    //! alias of `getRealRoots`
    integer
    get_real_roots( real_type r[] ) const
    { return getRealRoots( r ); }

    //!
    //! Get positive real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2 or 3
    //!
    integer getPositiveRoots( real_type r[] ) const;

    //! alias of `getPositiveRoots`
    integer
    get_positive_roots( real_type r[] ) const
    { return getPositiveRoots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2 or 3
    //!
    integer getNegativeRoots( real_type r[] ) const;

    //! alias of `getPositiveRoots`
    integer
    get_negative_roots( real_type r[] ) const
    { return getNegativeRoots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer getRootsInRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInRange`
    integer
    get_roots_in_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInRange( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInOpenRange`
    integer
    get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInOpenRange( a, b, r ); }

    //!
    //! First real root.
    //!
    real_type real_root0() const { return r0; }

    //!
    //! Second real root.
    //!
    real_type real_root1() const { return r1; }

    //!
    //! Third real root.
    //!
    real_type real_root2() const { return r2; }

    //!
    //! First complex or real root.
    //!
    complex_type
    root0() const
    { return cplx ? complex_type(r0,r1) : complex_type(r0,0); }

    //!
    //! Second complex or real root.
    //!
    complex_type
    root1() const
    { return cplx ? complex_type(r0,-r1) : complex_type(r1,0); }

    //!
    //! Third complex or real root.
    //!
    complex_type
    root2() const
    { return complex_type(r2,0); }

    //!
    //! First complex or real root.
    //!
    void
    getRoot0( real_type & re, real_type & im ) const {
      if ( cplx ) { re = r0; im = r1; }
      else        { re = r0; im = 0;  }
    }

    //! alias of `getRoot0`
    void
    get_root0( real_type & re, real_type & im ) const
    { getRoot0( re, im ); }

    //!
    //! First complex or real root.
    //!
    void
    getRoot0( complex_type & r ) const
    { r = cplx ? complex_type(r0,r1) : complex_type(r0,0); }

    //! alias of `getRoot0`
    void get_root0( complex_type & r ) const { getRoot0( r ); }

    //!
    //! Second complex or real root.
    //!
    void
    getRoot1( real_type & re, real_type & im ) const {
      if ( cplx ) { re = r0; im = -r1; }
      else        { re = r1; im = 0;   }
    }

    //! alias of `getRoot1`
    void
    get_root1( real_type & re, real_type & im ) const
    { getRoot1( re, im ); }

    //!
    //! Second complex or real root.
    //!
    void
    getRoot1( complex_type & r ) const
    { r = cplx ? complex_type(r0,-r1) : complex_type(r1,0); }

    //! alias of `getRoot1`
    void get_root1( complex_type & r ) const { getRoot1( r ); }

    //!
    //! Third complex or real root.
    //!
    void
    getRoot2( real_type & re, real_type & im ) const
    { re = r2; im = 0; }

    //! alias of `getRoot1`
    void
    get_root2( real_type & re, real_type & im ) const
    { getRoot2( re, im ); }

    //!
    //! Third complex or real root.
    //!
    void
    getRoot2( complex_type & r ) const
    { r = complex_type(r2,0); }

    //! alias of `getRoot2`
    void get_root2( complex_type & r ) const { getRoot2( r ); }

    //!
    //! Evalute the cubic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$
    //! \return the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type const & x ) const
    { return evalPoly( ABCD, 3, x ); }

    //!
    //! Evalute the cubic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type const & x ) const
    { return evalPolyC( ABCD, 3, x ); }

    //!
    //! Evalute the polynomial with its derivative.
    //!
    //! \param x   value where compute \f$ p(x) \f$
    //! \param p   value \f$ p(x) \f$
    //! \param dp  value \f$ p'(x) \f$
    //!
    void eval( real_type x, real_type & p, real_type & dp ) const;

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
   |  A * x^3 + B * x^2 + C * x + D
  \*/
  //! Quartic polynomial class
  //!
  //! **Constructor**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double a = 1;
  //!      double b = 2;
  //!      double c = 3;
  //!      double d = 3;
  //!      double e = 3;
  //!      Quartic p(a,b,c,d,e); // build an solve ``a x^4 + b x^3 + c x^2 + d x + e = 0``
  //!
  //!
  //!      Quartic p;
  //!      p.setup(a,b,c,d,e); // build an solve ``a x^4 + b x^3 + c x^2 + d x + e = 0``
  //!
  //! \endrst
  //!
  //! **Get kind of solution**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      int nroots = p.numRoots();
  //!      int nroots = p.numRealRoots();
  //!      int nroots = p.numComplexRoots();
  //!
  //! \endrst
  //!
  //! **Get real roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double r_min = 0;
  //!      double r_max = 2;
  //!      double r[4];
  //!      int nroots;
  //!      nroots = p.getRealRoots( r );
  //!      nroots = p.getPositiveRoots( r );
  //!      nroots = p.getNegativeRoots( r );
  //!      nroots = p.getRootsInRange( r_min, r_max, r );
  //!      nroots = p.getRootsInOpenRange( r_min, r_max, r );
  //!
  //! \endrst
  //!
  //! **Get roots**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      double r0 = p.real_root0();
  //!      double r1 = p.real_root1();
  //!      double r2 = p.real_root2();
  //!      double r3 = p.real_root3();
  //!      complex_type r0 = p.root0();
  //!      complex_type r1 = p.root1();
  //!      complex_type r2 = p.root2();
  //!      complex_type r3 = p.root3();
  //!
  //!      complex_type r;
  //!      double re, im;
  //!      p.getRoot0( re, im );
  //!      p.getRoot0( r );
  //!      p.getRoot1( re, im );
  //!      p.getRoot1( r );
  //!      p.getRoot2( re, im );
  //!      p.getRoot2( r );
  //!      p.getRoot3( re, im );
  //!      p.getRoot3( r );
  //!
  //! \endrst
  //!
  //! **Evaluate polynomial**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      {double or complex} v, x;
  //!      v = p.eval( x );
  //!
  //!      p.eval( x, p, dp );
  //!
  //! \endrst
  //!
  //! **Information**
  //!
  //! \rst
  //!
  //! .. code-block:: cpp
  //!
  //!      p.info( cout );
  //!      bool ok = p.check( cout );
  //!
  //! \endrst
  //!
  class Quartic {
    real_type ABCDE[5];
    real_type r0, r1, r2, r3;
    integer   iter, nreal, ncplx;

    void findRoots();

    bool cplx0() const { return ncplx > 0; }
    bool cplx1() const { return ncplx > 0; }
    bool cplx2() const { return ncplx > 2; }
    bool cplx3() const { return ncplx > 2; }

  public:

    Quartic() : iter(0), nreal(0), ncplx(0) {}
    Quartic(
      real_type _a,
      real_type _b,
      real_type _c,
      real_type _d,
      real_type _e
    )
    : iter(0), nreal(0), ncplx(0) {
      real_type & A = ABCDE[0];
      real_type & B = ABCDE[1];
      real_type & C = ABCDE[2];
      real_type & D = ABCDE[3];
      real_type & E = ABCDE[4];
      A = _a; B = _b; C = _c; D = _d; E = _e;
      findRoots();
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
      real_type & A = ABCDE[0];
      real_type & B = ABCDE[1];
      real_type & C = ABCDE[2];
      real_type & D = ABCDE[3];
      real_type & E = ABCDE[4];
      A = a; B = b; C = c; D = d; E = e;
      findRoots();
    }

    //!
    //! Number of roots found.
    //!
    integer numRoots() const { return nreal+ncplx; }

    //! alias of `numRoots`
    integer num_roots() const { return nreal+ncplx; }

    //!
    //! Number of real roots.
    //!
    integer numRealRoots() const { return nreal; }

    //! alias of `numRealRoots`
    integer num_real_roots() const { return nreal; }

    //!
    //! Number of complex roots
    //!
    integer numComplexRoots() const { return ncplx; }

    //! alias of `numComplexRoots`
    integer num_complex_roots() const { return ncplx; }

    //!
    //! Get the real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots, 0, 1, 2, 3 or 4
    //!
    integer getRealRoots( real_type r[] ) const;

    //! alias of `getRealRoots`
    integer
    get_real_roots( real_type r[] ) const
    { return getRealRoots(r); }

    //!
    //! Get positive real roots
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of positive real roots, 0, 1, 2, 3 or 4
    //!
    integer getPositiveRoots( real_type r[] ) const;

    //! alias of `getPositiveRoots`
    integer
    get_positive_roots( real_type r[] ) const
    { return getPositiveRoots( r ); }

    //!
    //! Get negative real roots.
    //!
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of negative real roots, 0, 1, 2, 3 or 4
    //!
    integer getNegativeRoots( real_type r[] ) const;

    //! alias of `getNegativeRoots`
    integer
    get_negative_roots( real_type r[] ) const
    { return getNegativeRoots( r ); }

    //!
    //! Get real roots in a closed range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the range [a,b]
    //!
    integer getRootsInRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInRange`
    integer
    get_roots_in_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInRange( a, b, r ); }

    //!
    //! Get real roots in an open range.
    //!
    //! \param[in]  a left side of the range
    //! \param[in]  b right side of the range
    //! \param[out] r vector that will be filled with the real roots
    //! \return the total number of real roots in the open range (a,b)
    //!
    integer getRootsInOpenRange( real_type a, real_type b, real_type r[] ) const;

    //! alias of `getRootsInRange`
    integer
    get_roots_in_open_range( real_type a, real_type b, real_type r[] ) const
    { return getRootsInOpenRange( a, b, r ); }

    //!
    //! First real root.
    //!
    real_type real_root0() const { return r0; }

    //!
    //! Second real root.
    //!
    real_type real_root1() const { return r1; }

    //!
    //! Third real root.
    //!
    real_type real_root2() const { return r2; }

    //!
    //! Fourth real root.
    //!
    real_type real_root3() const { return r3; }

    //!
    //! First real or complex root.
    //!
    complex_type
    root0() const
    { return cplx0() ? complex_type(r0,r1) : complex_type(r0,0); }

    //!
    //! Second real or complex root.
    //!
    complex_type
    root1() const
    { return cplx1() ? complex_type(r0,-r1) : complex_type(r1,0); }

    //!
    //! Third real or complex root.
    //!
    complex_type
    root2() const
    { return cplx2() ? complex_type(r2,r3) : complex_type(r2,0); }

    //!
    //! 4th real or complex root.
    //!
    complex_type
    root3() const
    { return cplx3() ? complex_type(r2,-r3) : complex_type(r3,0); }

    //!
    //! First real or complex root.
    //!
    void
    getRoot0( real_type & re, real_type & im ) const {
      if ( cplx0() ) { re = r0; im = r1; }
      else           { re = r0; im = 0;  }
    }

    //! alias of `getRoot0`
    void
    get_root0( real_type & re, real_type & im ) const {
      getRoot0( re, im );
    }

    //!
    //! First real or complex root.
    //!
    void
    getRoot0( complex_type & r ) const {
      if ( cplx0() ) r = complex_type(r0,r1);
      else           r = complex_type(r0,0);
    }

    //! alias of `getRoot0`
    void get_root0( complex_type & r ) const { getRoot0( r ); }

    //!
    //! Second real or complex root.
    //!
    void
    getRoot1( real_type & re, real_type & im ) const {
      if ( cplx1() ) { re = r0; im = -r1; }
      else           { re = r1; im = 0;   }
    }

    //! alias of `getRoot1`
    void
    get_root1( real_type & re, real_type & im ) const {
      getRoot1( re, im );
    }

    //!
    //! Second real or complex root.
    //!
    void
    getRoot1( complex_type & r ) const {
      if ( cplx1() ) r = complex_type(r0,-r1);
      else           r = complex_type(r1,0);
    }

    //! alias of `getRoot1`
    void get_root1( complex_type & r ) const { getRoot1( r ); }

    //!
    //! Third real or complex root.
    //!
    void
    getRoot2( real_type & re, real_type & im ) const {
      if ( cplx2() ) { re = r2; im = r3; }
      else           { re = r2; im = 0;  }
    }

    //! alias of `getRoot2`
    void
    get_root2( real_type & re, real_type & im ) const {
      getRoot2( re, im );
    }

    //!
    //! Third real or complex root.
    //!
    void
    getRoot2( complex_type & r ) const {
      if ( cplx2() ) r = complex_type(r2,r3);
      else           r = complex_type(r2,0);
    }

    //! alias of `getRoot2`
    void get_root2( complex_type & r ) const { getRoot2( r ); }

    //!
    //! 4th real or complex root.
    //!
    void
    getRoot3( real_type & re, real_type & im ) const {
      if ( cplx3() ) { re = r2; im = -r3; }
      else           { re = r3; im = 0;   }
    }

    //! alias of `getRoot3`
    void
    get_root3( real_type & re, real_type & im ) const {
      getRoot3( re, im );
    }

    //!
    //! 4th real or complex root.
    //!
    void
    getRoot3( complex_type & r ) const {
      if ( cplx3() ) r = complex_type(r2,-r3);
      else           r = complex_type(r3,0);
    }

    //! alias of `getRoot3`
    void get_root3( complex_type & r ) const { getRoot3( r ); }

    //!
    //! Evalute the quartic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    real_type
    eval( real_type const & x ) const
    { return evalPoly( ABCDE, 4, x ); }

    //!
    //! Evalute the quartic polynomial.
    //!
    //! \param x   value where compute \f$ p(x) \f$, x complex
    //! \return the value \f$ p(x) \f$
    //!
    complex_type
    eval( complex_type const & x ) const
    { return evalPolyC( ABCDE, 4, x ); }

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
