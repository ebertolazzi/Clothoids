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

#ifndef RPOLY_HH
#define RPOLY_HH

#include <utility>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

/*
..
.. N. FLOCKE
.. Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
.. for Physical Applications
.. ACM TOMS, Vol. 41, No. 4, 2015.
.. DOI: http://dx.doi.org/10.1145/2699468
..
*/

#include <cstdint>

namespace PolynomialRoots {

  using real_type    = double;
  using integer      = int;
  using complex_type = std::complex<real_type>;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  static int       const bitsValueType = std::numeric_limits<real_type>::digits;
  static real_type const splitFactor   = static_cast<real_type>((std::uint64_t(1)<<(bitsValueType-2))+1); // one extra digit is implicitly 1

  /*
  ||         _   _ _
  ||   _   _| |_(_) |___
  ||  | | | | __| | / __|
  ||  | |_| | |_| | \__ \
  ||   \__,_|\__|_|_|___/
  */
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // a + b = x + err
  static
  inline
  void
  TwoSum(
    real_type   a,
    real_type   b,
    real_type & x,
    real_type & err
  ) {
    x = a+b;
    real_type z = x-a;
    err = (a-(x-z))+(b-z);
  }

  static
  inline
  void
  TwoSum(
    complex_type   a,
    complex_type   b,
    complex_type & x,
    complex_type & err
  ) {
    real_type s1, e1, s2, e2;
    TwoSum( a.real(), b.real(), s1, e1 );
    TwoSum( a.imag(), b.imag(), s2, e2 );
    x   = complex_type(s1,s2);
    err = complex_type(e1,e2);
  }

  // a = x + y
  static
  inline
  void
  Split( real_type a, real_type & x, real_type & y ) {
    real_type c = splitFactor*a;
    x = c-(c-a);
    y = a-x;
  }

  // a * b = x + err
  static
  inline
  void
  TwoProduct(
    real_type   a,
    real_type   b,
    real_type & x,
    real_type & err
  ) {
    real_type a1, a2, b1, b2;
    Split( a, a1, a2 );
    Split( b, b1, b2 );
    x   = a*b;
    err = a2*b2-(((x-a1*b1)-a2*b1)-a1*b2);
  }

  static
  inline
  void
  TwoProduct(
    complex_type   a,
    complex_type   b,
    complex_type & p,
    complex_type & e,
    complex_type & f,
    complex_type & g
  ) {
    real_type z1, z2, z3, z4, z5, z6, h1, h2, h3, h4, h5, h6;
    TwoProduct(a.real(), b.real(), z1, h1 );
    TwoProduct(a.imag(), b.imag(), z2, h2 );
    TwoProduct(a.real(), b.imag(), z3, h3 );
    TwoProduct(a.imag(), b.real(), z4, h4 );
    TwoSum(z1, -z2, z5, h5);
    TwoSum(z3, z4, z6, h6);
    p = complex_type(z5,z6);
    e = complex_type(h1,h3);
    f = complex_type(-h2,h4);
    g = complex_type(h5,h6);
  }

  #endif

}

#endif
