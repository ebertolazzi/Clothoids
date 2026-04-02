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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-function"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wvla-extension"
#pragma clang diagnostic ignored "-Wvla"
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include "PolynomialRoots.hh"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace PolynomialRoots {

  using  std::abs;
  static constexpr real_type machepsi{ std::numeric_limits<real_type>::epsilon() };

  integer
  Quadratic::get_real_roots( real_type r[] ) const {
    integer nr{0};
    if ( !m_cplx ) {
      r[nr++] = m_r0;
      if ( m_nrts > 1 ) r[nr++] = m_r1;
    }
    return nr;
  }

  integer
  Quadratic::get_positive_roots( real_type r[] ) const {
    integer nr{0};
    if ( !m_cplx ) {
      if ( m_nrts > 0 && m_r0 > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > 0 ) r[nr++] = m_r1;
    }
    return nr;
  }

  integer
  Quadratic::get_negative_roots( real_type r[] ) const {
    integer nr{0};
    if ( !m_cplx ) {
      if ( m_nrts > 0 && m_r0 < 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 < 0 ) r[nr++] = m_r1;
    }
    return nr;
  }

  integer
  Quadratic::get_roots_in_range( real_type const a, real_type const b, real_type r[] ) const {
    integer nr{0};
    if ( !m_cplx ) {
      if ( m_nrts > 0 && m_r0 >= a && m_r0 <= b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 >= a && m_r1 <= b ) r[nr++] = m_r1;
    }
    return nr;
  }

  integer
  Quadratic::get_roots_in_open_range( real_type const a, real_type const b, real_type r[] ) const {
    integer nr{0};
    if ( !m_cplx ) {
      if ( m_nrts > 0 && m_r0 > a && m_r0 < b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > a && m_r1 < b ) r[nr++] = m_r1;
    }
    return nr;
  }

  // void
  // Quadratic::eval( real_type x, real_type & p, real_type & dp ) const {
  //   real_type const & A = ABC[0];
  //   real_type const & B = ABC[1];
  //   real_type const & C = ABC[2];
  //   if ( std::abs(x) > 1 ) {
  //     real_type z  = 1/x;
  //     p = (((C*z+B)*z+A)/x)/x;
  //   } else {
  //     p = (A*x+B)*x+C;
  //   }
  //   dp = 2*A*x+B;
  // }

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/

  void
  Quadratic::find_roots() {
    real_type const & A{ m_ABC[0] };
    real_type const & B{ m_ABC[1] };
    real_type const & C{ m_ABC[2] };

    m_r0 = m_r1 = 0;
    m_nrts = 0;
    m_cplx = m_dblx = false;

    if ( isZero(A) ) { // less than two roots b*z + c = 0
      if ( !isZero(B) ) { m_nrts = 1; m_r0 = -C/B; }
    } else if ( isZero(C) ) { // a*z^2 + b*z  = 0
      m_nrts = 2;
      m_dblx = isZero(B);
      if ( !m_dblx ) {
        m_r0 = -B/A;
        if ( m_r0 < 0 ) std::swap(m_r0,m_r1);
      }
    } else { // Compute discriminant avoiding overflow.
      real_type const hb    { B/2          }; // b now b/2
      real_type const abs_b { std::abs(hb) };
      real_type const abs_c { std::abs(C)  };
      real_type e, d;
      if ( abs_b < abs_c ) {
        e = C < 0 ? -A : A;
        e = (hb*hb)-e*abs_c;
        d = std::sqrt(std::abs(e));
      } else {
        e = 1 - (A/hb)*(C/hb);
        d = std::sqrt(std::abs(e))*abs_b;
      }
      m_nrts = 2;
      m_cplx = e < 0;
      if ( m_cplx ) {
        // complex conjugate zeros
        m_r0 = -hb/A;         // real part
        m_r1 = std::abs(d/A); // immaginary part
      } else {
        // real zeros
        m_dblx = isZero(d);
        if ( m_dblx ) {
          m_r0 = m_r1 = -hb/A;
        } else {
          if ( hb >= 0 ) d = -d;
          m_r0 = (d-hb)/A;
          //r1 = (-d-hb)/a;
          if ( !isZero(m_r0) ) m_r1 = (C/m_r0)/A;
          if ( m_r0 > m_r1 ) std::swap(m_r0,m_r1); // order roots
        }
      }
    }
  }

  void
  Quadratic::info( ostream_type & s ) const {
    real_type const & A{ m_ABC[0] };
    real_type const & B{ m_ABC[1] };
    real_type const & C{ m_ABC[2] };
    s << "\npoly A=" << A << " B=" << B << " C=" << C
      << "\nn. roots = " << m_nrts
      << "\ncomplex  = " << (m_cplx?"YES":"NO")
      << "\ndouble   = " << (m_dblx?"YES":"NO");
    if ( m_cplx ) {
      s << "\nx0 = (" << m_r0 << "," <<  m_r1 << ')'
        << "\nx1 = (" << m_r0 << "," << -m_r1 << ')';
    } else if ( m_dblx ) {
      s << "\nx0 = x1 = " << m_r0;
    } else if ( m_nrts == 1 ) {
      s << "\nx0 = " << m_r0;
    } else if ( m_nrts == 2 ) {
      s << "\nx0 = " << m_r0
        << "\nx1 = " << m_r1;
    }
    s << '\n';
  }

  bool
  Quadratic::check( ostream_type & s ) const {
    real_type const & A{ m_ABC[0] };
    real_type const & B{ m_ABC[1] };
    real_type const & C{ m_ABC[2] };
    bool ok{ true };
    real_type const epsi{ 10 * ( std::abs(A) +
                                 std::abs(B) +
                                 std::abs(C) ) * machepsi };
    if ( m_cplx ) {
      real_type const z0{ std::abs(eval( root0() )) };
      real_type const z1{ std::abs(eval( root1() )) };
      s << "|p(r0)| = " << z0
        << "\n|p(r1)| = " << z1
        << '\n';
      ok = z0 < epsi && z1 < epsi;
    } else if ( m_nrts == 1 ) {
      real_type const z0{ eval( real_root0() ) };
      s << "p(r0) = " << z0  << '\n';
      ok = std::abs(z0) < epsi;
    } else if ( m_nrts == 2 ) {
      real_type const z0{ eval( real_root0() ) };
      real_type const z1{ eval( real_root1() ) };
      s << "p(r0) = " << z0
        << "\np(r1) = " << z1
        << '\n';
      ok = std::abs(z0) < epsi && std::abs(z1) < epsi;
    }
    return ok;
  }

}

#endif

// EOF: PolynomialRoots-1-Quadratic.cc
