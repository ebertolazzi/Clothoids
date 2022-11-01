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
#include "PolynomialRoots-Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace PolynomialRoots {

  // static real_type const machepsi = std::numeric_limits<real_type>::epsilon();

  using std::pair;
  using std::abs;
  using std::pow;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  real_type
  evalPoly(
    real_type const op[],
    integer         Degree,
    real_type       x
  ) {
    bool reverse = std::abs(x) > 1;
    if ( reverse ) {
      real_type res(op[Degree]);
      real_type xn(1);
      for ( integer i = 1; i <= Degree; ++i ) {
        res /= x;
        res += op[Degree-i];
        xn  *= x;
      }
      res *= xn;
      return res;
    } else {
      real_type res(op[0]);
      for ( integer i = 1; i <= Degree; ++i ) {
        res *= x;
        res += op[i];
      }
      return res;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of Newton step
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  bool
  NewtonStep(
    real_type const op[],
    integer         Degree,
    real_type     & x
  ) {
    real_type p, dp, xn;
    bool reverse = std::abs(x) > 1;
    if ( reverse ) {
      p  = op[Degree];
      xn = 1;
      for ( integer i = 1; i <= Degree; ++i ) {
        p  /= x;
        p  += op[Degree-i];
        xn *= x;
      }
      p *= xn;

      dp = op[Degree];
      xn = 1;
      for ( integer i = 2; i <= Degree; ++i ) {
        dp /= x;
        dp += (Degree-i)*op[Degree-i];
        xn *= x;
      }
      dp *= xn;
    } else {
      dp = Degree*op[0];
      p  = op[0];
      for ( integer i = 1; i < Degree; ++i ) {
        p  *= x;
        p  += op[i];
        dp *= x;
        dp += (Degree-i)*op[i];
      }
      p *= x;
      p += op[Degree];
    }
    x -= p/dp;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial and its derivative
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  void
  evalPolyDPoly(
    real_type const op[],
    integer         Degree,
    real_type       x,
    real_type     & p,
    real_type     & dp
  ) {
    bool reverse = std::abs(x) > 1;
    if ( reverse ) {
      real_type const * pc = op+Degree;
      real_type rD = Degree;
      p  = *pc--;
      dp = rD*p;
      real_type xn(1);
      for ( integer i = 1; i < Degree; ++i ) {
        p  /= x;
        p  += *pc;
        dp /= x;
        dp += rD*(*pc--);
        xn *= x;
      }
      dp *= xn;
      p   = p/x + *pc;
      xn *= x;
      p  *= xn;
    } else {
      p = op[0]*x+op[1];
      for ( integer i = 2; i <= Degree; ++i ) {
        p  *= x;
        p  += op[i];
        dp *= x;
        dp += i*op[i];
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  complex_type
  evalPolyC(
    real_type const      op[],
    integer              Degree,
    complex_type const & x
  ) {
    bool reverse = std::abs(x) > 1;
    if ( reverse ) {
      complex_type res(op[Degree]);
      complex_type xn(1,0);
      for ( integer i = 1; i <= Degree; ++i ) {
        res /= x;
        res += op[Degree-i];
        xn  *= x;
      }
      res *= xn;
      return res;
    } else {
      complex_type res(op[0]);
      for ( integer i = 1; i <= Degree; ++i ) {
        res *= x;
        res += op[i];
      }
      return res;
    }
  }

  //============================================================================

  /*
  ..  scale roots by pair.second, after scaling the polinomial has the form
  ..
  ..  p[0] + p[1]*x + ... + p[n]*x^n
  ..
  ..  pair.first is the index such that p[pair.first] = +-1
  ..
  */

  static
  pair<integer,real_type>
  scalePolynomial(
    integer         n, // degree
    real_type const p[],
    real_type       ps[]
  ) {
    integer   i_max = n;
    real_type an    = p[n];
    real_type scale = -1;
    integer   i     = n;
    while ( --i >= 0 ) {
      ps[i] = p[i]/an;
      real_type scale1 = std::pow( std::abs(ps[i]), 1.0/(n-i) );
      if ( scale1 > scale ) { scale = scale1; i_max = i; }
    }
    // scale coeffs
    pair<integer,real_type> res(i_max,scale);
    real_type s = scale;
    for ( i = 1; i <= n; ++i, s *= scale ) ps[n-i] /= s;
    ps[n] = 1;
    return res;
  }

  //============================================================================

  /*
  .. divide a(x)  by (x-r) so that a(x) = (x-r) * b(x)
  */

  static
  void
  deflatePolynomial(
    integer         n, // degree
    real_type const a[],
    real_type       r,
    real_type       b[]
  ) {
    // crossover index for forward/backward deflation
    // G. Peters and J. H. Wilkinson.
    // Practical problems arising in the solution of polynomial equations.
    // J. Inst. Math. Appl. 8 (1971), 16–35.
    integer   i_cross = 0;
    real_type v_cross = std::abs(a[0]);
    real_type r1      = r;
    for ( integer i = 1; i < n; ++i, r1 *= r ) {
      real_type v_cross1 = std::abs(a[i]*r1);
      if ( v_cross1 > v_cross )
        { v_cross = v_cross1; i_cross = i; }
    }
    b[n-1] = 1;
    if ( isZero(a[n]-1) ) {
      if ( i_cross > 0 ) {
        b[0] = -a[0] / r;
        for ( integer j = 1; j < i_cross; ++j )
          b[j] = (a[j]-b[j-1]) / r;
      }
      for ( integer j = n-2; j >= i_cross; --j )
        b[j] = a[j+1]+r*b[j+1];
    } else {
      real_type an = a[n];
      if ( i_cross > 0 ) {
        b[0] = -(a[0]/an) / r;
        for ( integer j = 1; j < i_cross; ++j )
          b[j] = (a[j]/an-b[j-1]) / r;
      }
      for ( integer j = n-2; j >= i_cross; --j )
        b[j] = a[j+1]/an+r*b[j+1];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  // p0 + p1*x + p2*x^2 + ... + pn*x^n
  //
  // (p0/x^n + p1/x^(n-1) + p2/(x^(n-2) + ... + pn)*x^n
  //
  #if 0
  // UNUSED
  real_type
  CompHorner(
    real_type const p[],
    integer         Degree,
    real_type       x
  ) {
    real_type xabs = std::abs(x);
    bool reverse = xabs > 1;
    //bool reverse = false;
    if ( reverse ) x = real_type(1)/x;
    integer ii0 = reverse ? Degree : 0;
    real_type res(p[ii0]);
    real_type c = 0;
    for ( integer i = 1; i <= Degree; ++i ) {
      integer ii = reverse ? Degree-i : i;
      real_type tmp, pi, sigma;
      //TwoProduct( res, x, tmp, pi );
      //TwoSum( tmp, p[ii], res, sigma );
      res = res * x + p[ii];
      //c = c * x + (pi+sigma);
    }
    //res += c;
    if ( reverse ) res *= std::pow(x,Degree);
    return res;
  }
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^3 + A x^2 + B x + C
  static
  inline
  void
  scaleCubicMonicPolynomial(
    real_type   A,
    real_type   B,
    real_type   C,
    real_type & AS,
    real_type & BS,
    real_type & CS,
    integer   & i_case,
    real_type & scale
  ) {

    real_type a = std::abs(A);
    real_type b = std::sqrt(abs(B));
    real_type c = std::cbrt(abs(C));

    if ( a < b ) {
      if ( b < c ) i_case = 0; // a < b < c --> c MAX
      else         i_case = 1; // a < b and c <= b --> b MAX
    } else {
      if ( a < c ) i_case = 0; // b <= a < c --> c MAX
      else         i_case = 2; // b <= a  and c <= a --> a MAX
    }

    switch ( i_case ) {
    case 0:
      scale = c;
      AS    = A/c;
      BS    = (B/c)/c;
      CS    = C > 0 ? 1 : -1;
      break;
    case 1:
      scale = b;
      AS    = A/b;
      BS    = B > 0 ? 1 : -1;
      CS    = ((C/b)/b)/b;
      break;
    case 2:
      scale = a;
      AS    = A > 0 ? 1 : -1;
      BS    = (B/a)/a;
      CS    = ((C/a)/a)/a;
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a3*x^2+b1*x+b0)
  static
  void
  deflateCubicPolynomial(
    real_type   a3,
    real_type   a2,
    real_type   a1,
    real_type   a0,
    real_type   r,
    real_type & b1,
    real_type & b0
  ) {
    integer   i_cross  = 0;
    real_type r2       = r*r;
    real_type v_cross  = std::abs(a0);
    real_type v_cross1 = std::abs(a1*r);
    if ( v_cross1 > v_cross ) { v_cross = v_cross1; i_cross = 1; }
      v_cross1 = std::abs(a2*r2);
      if ( v_cross1 > v_cross ) { v_cross = v_cross1; i_cross = 2; }
      v_cross1 = std::abs(a3*r*r2);
      if ( v_cross1 > v_cross ) i_cross = 3;
      switch ( i_cross ) {
      case 0: b1 = a2+a3*r; b0 = a1+r*b1; break;
      case 1: b1 = a2+a3*r; b0 = -a0/r;   break;
      case 2:
      case 3: b0 = -a0/r; b1 = (b0-a1)/r; break;
    }
  }

}

#endif

// EOF: PolynomialRoots-Utils.cc
