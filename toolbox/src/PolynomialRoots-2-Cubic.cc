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
 |      Universita` degli Studi di Trento                                   |
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

#define MAX_ITER_SAFETY 50

namespace PolynomialRoots {

  using std::abs;
  static constexpr real_type machepsi { std::numeric_limits<real_type>::epsilon() };
  static constexpr real_type third    { 1./3. };
  static constexpr real_type one27th  { 1./27. };
  static constexpr real_type two27th  { 2./27. };

  integer
  Cubic::get_real_roots( real_type r[] ) const {
    integer nr{0};
    if ( m_cplx ) {
      if ( m_nrts > 2 ) r[nr++] = m_r2;
    } else {
      if ( m_nrts > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 ) r[nr++] = m_r1;
      if ( m_nrts > 2 ) r[nr++] = m_r2;
    }
    return nr;
  }

  integer
  Cubic::get_positive_roots( real_type r[] ) const {
    integer nr{0};
    if ( m_cplx ) {
      if ( m_nrts > 2 && m_r2 > 0  ) r[nr++] = m_r2;
    } else {
      if ( m_nrts > 0 && m_r0 > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > 0 ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 > 0 ) r[nr++] = m_r2;
    }
    return nr;
  }

  integer
  Cubic::get_negative_roots( real_type r[] ) const {
    integer nr{0};
    if ( m_cplx ) {
      if ( m_nrts > 2 && m_r2 < 0  ) r[nr++] = m_r2;
    } else {
      if ( m_nrts > 0 && m_r0 < 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 < 0 ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 < 0 ) r[nr++] = m_r2;
    }
    return nr;
  }

  integer
  Cubic::get_roots_in_range( real_type const a, real_type const b, real_type r[] ) const {
    integer nr{0};
    if ( m_cplx ) {
      if ( m_nrts > 2 && m_r2 >= a && m_r2 <= b ) r[nr++] = m_r2;
    } else {
      if ( m_nrts > 0 && m_r0 >= a && m_r0 <= b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 >= a && m_r1 <= b ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 >= a && m_r2 <= b ) r[nr++] = m_r2;
    }
    return nr;
  }

  integer
  Cubic::get_roots_in_open_range( real_type const a, real_type const b, real_type r[] ) const {
    integer nr{0};
    if ( m_cplx ) {
      if ( m_nrts > 2 && m_r2 > a && m_r2 < b ) r[nr++] = m_r2;
    } else {
      if ( m_nrts > 0 && m_r0 > a && m_r0 < b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > a && m_r1 < b ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 > a && m_r2 < b ) r[nr++] = m_r2;
    }
    return nr;
  }

  // void
  // Cubic::eval( real_type x, real_type & p, real_type & dp ) const {
  //   real_type const & A = ABCD[0];
  //   real_type const & B = ABCD[1];
  //   real_type const & C = ABCD[2];
  //   real_type const & D = ABCD[3];
  //   if ( std::abs(x) > 1 ) {
  //     real_type x2 = x*x;
  //     real_type x3 = x2*x;
  //     p  = (((D/x+C)/x+B)/x+A)*x3;
  //     dp = ((C/x+2*B)/x+3*A)*x2;
  //   } else {
  //     p  = ((A*x+B)*x+C)*x+D;
  //     dp = (3*A*x+2*B)*x+C;
  //   }
  // }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  static
  real_type
  guess1( real_type const a[3] ) {
    constexpr real_type p{  1.09574    };
    constexpr real_type q{ -3.239E-1   };
    constexpr real_type r{ -3.239E-1   };
    constexpr real_type s{  9.57439E-2 };
    return p+q*a[1]+r*a[2]+s*a[1]*a[2];
  }

  static
  real_type
  guess2( real_type const a[3] ) {
    constexpr real_type p{ -1.09574    };
    constexpr real_type q{  3.239E-1   };
    constexpr real_type r{ -3.239E-1   };
    constexpr real_type s{  9.57439E-2 };
    return p+q*a[1]+r*a[2]+s*a[1]*a[2];
  }

  static
  real_type
  guess3( real_type const a[3] ) {
    constexpr real_type p{  1.14413    };
    constexpr real_type q{ -2.75509E-1 };
    constexpr real_type r{ -4.45578E-1 };
    constexpr real_type s{ -2.59342E-2 };
    if ( real_type const t{a[2]/3}; a[0] < t*(2*t*t-1) ) return  p+q*a[0]+r*a[2]+s*a[0]*a[2];
    return -p+q*a[0]+r*a[2]-s*a[0]*a[2];
  }

  static
  real_type
  guess4( real_type const a[3] ) {
    constexpr real_type q{ -7.71845E-1 };
    constexpr real_type s{ -2.28155E-1 };
    if ( a[0] > 0 ) return (q+s*a[2])*a[0];
    return (q-s*a[2])*a[0];
  }

  static
  real_type
  guess5( real_type const a[3] ) {
    real_type p, q, r, s;
    real_type const tmp{ two27th-a[1]/3 };
    if ( a[1] <= third ) {
      if ( a[0] < tmp ) {
        p =  8.78558E-1;
        q = -5.71888E-1;
        r = -7.11154E-1;
        s = -3.22313E-1;
      } else {
        p = -1.92823E-1;
        q = -5.66324E-1;
        r = +5.05734E-1;
        s = -2.64881E-1;
      }
    } else {
      if ( a[0] < tmp ) {
        p = 1.19748;
        q = -2.83772E-1;
        r = -8.37476E-1;
        s = -3.56228E-1;
      } else {
        p = -3.45219E-1;
        q = -4.01231E-1;
        r =  2.07216E-1;
        s = -4.45532E-3;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1];
  }

  static
  real_type
  guess6( real_type const a[3] ) {
    real_type p, q, r, s;
    real_type const tmp{ a[1]/3-two27th };
    if ( a[1] <= third ) {
      if ( a[0] > tmp ) {
        p = -8.78558E-1;
        q = -5.71888E-1;
        r =  7.11154E-1;
        s = -3.22313E-1;
      } else {
        p =  1.92823E-1;
        q = -5.66324E-1;
        r = -5.05734E-1;
        s = -2.64881E-1;
      }
    } else {
      if ( a[0] > tmp ) {
        p = -1.19748;
        q = -2.83772E-1;
        r =  8.37476E-1;
        s = -3.56228E-1;
      } else {
        p =  3.45219E-1;
        q = -4.01231E-1;
        r = -2.07216E-1;
        s = -4.45532E-3;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1];
  }

  /*
  ||   _   _               _              ____  _               _   _
  ||  | \ | | _____      _| |_ ___  _ __ | __ )(_)___  ___  ___| |_(_) ___  _ __
  ||  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \|  _ \| / __|/ _ \/ __| __| |/ _ \| '_ \
  ||  | |\  |  __/\ V  V /| || (_) | | | | |_) | \__ \  __/ (__| |_| | (_) | | | |
  ||  |_| \_|\___| \_/\_/  \__\___/|_| |_|____/|_|___/\___|\___|\__|_|\___/|_| |_|
  */

  // x^3 + a * x^2 + b * x + c
  static
  integer
  NewtonBisection(
    real_type const a,
    real_type const b,
    real_type const c,
    real_type &     x
  ) {
    real_type p, dp;
    evalMonicCubic( x, a, b, c, p, dp );
    real_type t{p}; // save p(x) for sign comparison
    x -= p/dp; // 1st improved root

    integer iter       { 1     };
    integer oscillate  { 0     };
    integer nconverged { 0     };
    bool    bisection  { false };
    bool    converged  { false };
    real_type s(0), u(0); // to mute warning
    while ( ! ( nconverged > 1 || bisection ) && iter < MAX_ITER_SAFETY ) {
      ++iter;
      evalMonicCubic( x, a, b, c, p, dp );
      if ( p*t < 0 ) { // does Newton start oscillating ?
        if ( p < 0 ) {
          ++oscillate; // increment oscillation counter
          s = x;       // save lower bisection bound
        } else {
          u = x; // save upper bisection bound
        }
        t = p; // save current p(x)
      }
      dp = p/dp; // Newton correction
      x -= dp;   // new Newton root
      bisection = oscillate > 2; // activate bisection
      converged = std::abs(dp) <= (1+std::abs(x)) * machepsi; // Newton convergence indicator
      if ( converged ) ++nconverged; else nconverged = 0;
    }
    if ( bisection ) {
      t = u - s; // initial bisection interval
      while ( std::abs(t) > (1+std::abs(x)) * machepsi && iter < MAX_ITER_SAFETY ) { // bisection iterates
        ++iter;
        p = evalMonicCubic( x, a, b, c );
        if ( p < 0 ) s = x;
        else         u = x; // keep bracket on root
        t = (u-s)/2; // new bisection interval
        x = s + t;   // new bisection root
      }
    }
    return iter;
  }

  /*\
   *  Calculate the zeros of the cubic a*z^3 + b*z^2 + c*z + d.
  \*/

  void
  Cubic::find_roots() {
    real_type const & A{ m_ABCD[0] };
    real_type const & B{ m_ABCD[1] };
    real_type const & C{ m_ABCD[2] };
    real_type const & D{ m_ABCD[3] };
    m_nrts = m_iter = 0;
    m_cplx = m_dblx = m_trpx = false;
    // special cases
    if ( isZero(A) ) {
      Quadratic const qsolve( B, C, D );
      m_nrts = qsolve.num_roots();
      m_cplx = qsolve.complex_root();
      m_dblx = qsolve.double_root();
      m_r0   = qsolve.real_root0();
      m_r1   = qsolve.real_root1();
      return;
    }
    if ( isZero(D) ) {
      Quadratic const qsolve( A, B, C );
      m_nrts = qsolve.num_roots()+1;
      m_cplx = qsolve.complex_root();
      m_r0   = qsolve.real_root0();
      m_r1   = qsolve.real_root1();
      m_r2   = 0;
      if ( !m_cplx ) { // reorder
        if ( m_r2 < m_r1 ) std::swap( m_r1, m_r2 );
        if ( m_r1 < m_r0 ) std::swap( m_r0, m_r1 );
        if ( m_r2 < m_r1 ) std::swap( m_r1, m_r2 );
      }
      return;
    }
    /*              _
    ||  ___ __ __ _| |___
    || (_-</ _/ _` | / -_)
    || /__/\__\__,_|_\___|
    */
    // x^3 + aa * x^2 + bb * x + cc
    real_type const aa{ B/A };
    real_type const bb{ C/A };
    real_type const cc{ D/A };
    // scale Cubic Monic Polynomial
    real_type const absa{ std::abs(aa)            };
    real_type const absb{ std::sqrt(std::abs(bb)) };
    real_type const absc{ std::cbrt(std::abs(cc)) };

    integer i_case{ 0 }; // c MAX
    if ( absa < absb ) {
      if ( absc < absb ) i_case = 1; // |a| < |b| and |b| < |c| --> b MAX
      // |a| < |b| <= |c| --> c MAX
    } else {
      if ( absc < absa ) i_case = 2; // |b| <= |a| and |c| < |a| --> a MAX
      // |b| <= |a| < |c| --> c MAX
    }

    real_type scale(0), a[3];
    switch ( i_case ) {
    case 0:
      scale = absc;
      a[2]  = aa/absc;
      a[1]  = (bb/absc)/absc;
      a[0]  = cc > 0 ? 1 : -1;
      break;
    case 1:
      scale = absb;
      a[2]  = aa/absb;
      a[1]  = bb > 0 ? 1 : -1;
      a[0]  = ((cc/absb)/absb)/absb;
      break;
    case 2:
      scale = absa;
      a[2]  = aa > 0 ? 1 : -1;
      a[1]  = (bb/absa)/absa;
      a[0]  = ((cc/absa)/absa)/absa;
      break;
    }

    /*
    ||   __ _ _  _ ___ ______
    ||  / _` | || / -_|_-<_-<
    ||  \__, |\_,_\___/__/__/
    ||  |___/
    */
    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    integer iclass{-1};
    switch ( i_case ) {
      case 0: iclass = a[0] > 0 ? 2 : 1; break;
      case 1: iclass = a[1] > 0 ? 4 : 3; break;
      case 2: iclass = a[2] > 0 ? 6 : 5; break;
    }
    bool use_shifted{false};
    m_trpx = false;
    switch ( iclass ) {
    case 1: m_r2 = guess1(a); break;
    case 2: m_r2 = guess2(a); break;
    case 3: m_r2 = guess3(a); break;
    case 4: m_r2 = guess4(a); break;
    case 5:
      m_r0   = a[1]-third;
      m_r1   = a[0]+one27th;
      m_trpx = std::abs(m_r0) <= machepsi && std::abs(m_r1) <= machepsi; // check for triple root
      if ( m_trpx ) { m_r0 = m_r1 = m_r2 = third * scale; m_nrts = 3; return; }
      use_shifted = std::abs(m_r0) <= 0.01 && std::abs(m_r1) <= 0.01;
      m_r2 = guess5(a);
      break;
    case 6:
      m_r0   = a[1]-third;
      m_r1   = a[0]-one27th;
      m_trpx = std::abs(m_r0) <= machepsi && std::abs(m_r1) <= machepsi; // check for triple root
      if ( m_trpx ) { m_r0 = m_r1 = m_r2 = -third * scale; m_nrts = 3; return; }
      use_shifted = std::abs(m_r0) <= 0.01 && std::abs(m_r1) <= 0.01;
      m_r2 = guess6(a);
      break;
    }

    /*
    ||          _
    ||  ___ ___| |_ _____
    || (_-</ _ \ \ V / -_)
    || /__/\___/_|\_/\___|
    */
    m_iter = 0;
    if ( use_shifted ) {
      if ( iclass == 5 ) {
        // a[2] == -1
        // y^3 + (a[1]-1/3)* y + (a[0]+a[1]/3-2/27), x = y+1/3
        m_r2 -= third; // shift guess
        m_iter = NewtonBisection( 0, m_r0, a[0]+a[1]/3-two27th, m_r2 );
        m_r2 += third; // unshift solution
      } else {
        // a[2] == 1
        // y^3 + (a[1]-1/3)* y + (a[0]-a[1]/3+2/27), x = y+1/3
        m_r2 += third; // shift guess
        m_r1 -= a[1]/3-one27th;
        //if ( std::abs(r3) < machepsi ) r3 = 0;
        m_iter = NewtonBisection( 0, m_r0, a[0]-a[1]/3+two27th, m_r2 );
        m_r2 -= third; // unshift solution
      }
    } else {
      m_iter = NewtonBisection( a[2], a[1], a[0], m_r2 );
    }

    // scale root
    m_r2 *= scale;

    /*
    // deflate
    // x^3 + aa*x^2 + bb*x + cc
    //    = (x-r2)*(x^2+b1*x+b0)
    //    = x^3 + x^2 * ( b1 - r2 ) + x * ( b0 - b1*r2 ) - r2 * b0
    //
    //    aa == b1 - r2
    //    bb == b0 - b1*r2
    //    cc == -r2 * b0
    //
    //  Solve the overdetermined linear system:
    //
    //  / 0    1  \            / aa + r2 \
    //  |         |  / b0 \    |         |
    //  | 1   -r2 |  |    |  = |   bb    |
    //  |         |  \ b1 /    |         |
    //  \ -r2  0  /            \   cc    /
    //
    //  if |r2| < 1 then solve
    //
    //  / 0    1  \  / b0 \    / aa + r2 \
    //  |         |  |    | =  |         |
    //  \ -r2  0  /  \ b1 /    \   cc    /
    //
    //  otherwise solve
    //
    //  / 1   -r2 \  / b0 \   / bb \
    //  |         |  |    | = |    |
    //  \ -r2  0  /  \ b1 /   \ cc /
    */
    real_type const b0 { -cc/m_r2 };
    real_type const b1 { std::abs(m_r2) < 1 ? aa+m_r2 : -(cc/m_r2+bb)/m_r2 };

    // solve quadratic polynomial
    Quadratic const qsolve( 1.0, b1, b0 );
    m_nrts = qsolve.num_roots()+1;
    m_cplx = qsolve.complex_root();
    m_dblx = qsolve.double_root();
    m_r0   = qsolve.real_root0();
    m_r1   = qsolve.real_root1();

    if ( !m_cplx ) { // if real roots sort it!
      if ( m_r1 > m_r2 ) std::swap(m_r1,m_r2);
      if ( m_r0 > m_r1 ) std::swap(m_r0,m_r1);
      if ( m_r1 > m_r2 ) std::swap(m_r1,m_r2);
    }

  }

  void
  Cubic::info( ostream_type & s ) const {
    real_type const & A{ m_ABCD[0] };
    real_type const & B{ m_ABCD[1] };
    real_type const & C{ m_ABCD[2] };
    real_type const & D{ m_ABCD[3] };
    s << "\npoly a=" << A << " b=" << B << " c=" << C << " d=" << D
      << "\nn. roots = " << m_nrts
      << "\ncomplex  = " << (m_cplx?"YES":"NO")
      << "\ntriple   = " << (m_trpx?"YES":"NO")
      << "\ndouble   = " << (m_dblx?"YES":"NO");
    if ( m_cplx ) {
      s << "\nx0 = (" << m_r0 << "," <<  m_r1 << ')'
        << "\nx1 = (" << m_r0 << "," << -m_r1 << ')';
      if ( m_nrts > 2 ) s << "\nx3 = " << m_r2;
    } else {
      if ( m_nrts > 0 ) s << "\nx0 = " << m_r0;
      if ( m_nrts > 1 ) s << "\nx1 = " << m_r1;
      if ( m_nrts > 2 ) s << "\nx2 = " << m_r2;
    }
    s << '\n';
  }

  bool
  Cubic::check( ostream_type & s ) const {
    real_type const & A{ m_ABCD[0] };
    real_type const & B{ m_ABCD[1] };
    real_type const & C{ m_ABCD[2] };
    real_type const & D{ m_ABCD[3] };
    bool ok{ true };
    real_type const epsi{ 10 * ( std::abs(A) +
                                 std::abs(B) +
                                 std::abs(C) +
                                 std::abs(D) ) * machepsi };
    if ( m_cplx ) {
      real_type const z0{ std::abs(eval( root0() )) };
      real_type const z1{ std::abs(eval( root1() )) };
      real_type const z2{ std::abs(eval( root2() )) };
      real_type const zr{ eval( real_root2() )      };
      s << "|p(r0)| = " << z0
        << "\n|p(r1)| = " << z1
        << "\n|p(r2)| = " << z2
        << "\np(real_part(r2)) = " << zr
        << '\n';
      ok = z0 < epsi && z1 < epsi && z2 < epsi;
    } else if ( m_nrts == 1 ) {
      real_type const z0{ eval( real_root0() ) };
      s << "p(r0) = " << z0  << '\n';
      ok = std::abs(z0) < epsi;
    } else if ( m_nrts == 2 ) {
      real_type const z0{ std::abs(eval( root0() )) };
      real_type const z1{ std::abs(eval( root1() )) };
      s << "p(r0) = " << z0
        << "\np(r1) = " << z1
        << '\n';
      ok = std::abs(z0) < epsi && std::abs(z1) < epsi;
    } else if ( m_nrts == 3 ) {
      real_type const z0{ eval( real_root0() ) };
      real_type const z1{ eval( real_root1() ) };
      real_type const z2{ eval( real_root2() ) };
      s << "p(r0) = " << z0
        << "\np(r1) = " << z1
        << "\np(r2) = " << z2
        << '\n';
      ok = std::abs(z0) < epsi && std::abs(z1) < epsi && std::abs(z2) < epsi;
    }
    return ok;
  }

}

#endif

// EOF: PolynomialRoots-2-Cubic.cc
