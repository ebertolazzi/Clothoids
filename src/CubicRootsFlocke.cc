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

#include "CubicRootsFlocke.hh"

namespace G2lib {

  using std::abs;
  using std::swap;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^3 + A x^2 + B x + C
  static
  inline
  void
  scaleCubicMonicPolynomial( real_type   A,
                             real_type   B,
                             real_type   C,
                             real_type & AS,
                             real_type & BS,
                             real_type & CS,
                             int_type  & i_case,
                             real_type & scale ) {

    real_type a = abs(A);
    real_type b = sqrt(abs(B));
    real_type c = cbrt(abs(C));

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
  deflateCubicPolynomial( real_type   a3,
                          real_type   a2,
                          real_type   a1,
                          real_type   a0,
                          real_type   r,
                          real_type & b1,
                          real_type & b0 ) {
    int_type  i_cross  = 0;
    real_type r2       = r*r;
    real_type v_cross  = abs(a0);
    real_type v_cross1 = abs(a1*r);
    if ( v_cross1 > v_cross ) { v_cross = v_cross1; i_cross = 1; }
    v_cross1 = abs(a2*r2);
    if ( v_cross1 > v_cross ) { v_cross = v_cross1; i_cross = 2; }
    v_cross1 = abs(a3*r*r2);
    if ( v_cross1 > v_cross ) i_cross = 3;
    switch ( i_cross ) {
      case 0: b1 = a2+a3*r; b0 = a1+r*b1; break;
      case 1: b1 = a2+a3*r; b0 = -a0/r; break;
      case 2:
      case 3: b0 = -a0/r; b1 = (b0-a1)/r; break;
    }
  }

  // x^3 + a*x^2 + b*x + c
  static
  inline
  real_type
  evalMonicCubic( real_type x,
                  real_type a,
                  real_type b,
                  real_type c ) {
    return ((x+a)*x+b)*x+c;
  }

  static
  inline
  void
  evalMonicCubic( real_type   x,
                  real_type   a,
                  real_type   b,
                  real_type   c,
                  real_type & p,
                  real_type & dp ) {
    p  = x + a;
    dp = x + p;
    p  = p  * x + b;
    dp = dp * x + p;
    p  = p  * x + c;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  static
  int_type
  zeroCubicByNewtonBisection( real_type const a,
                              real_type const b,
                              real_type const c,
                              real_type     & x ) {

    real_type p, dp;
    evalMonicCubic( x, a, b, c, p, dp );
    real_type t = p; // save p(x) for sign comparison
    x -= p/dp; // 1st improved root

    int_type  iter      = 1;
    int_type  oscillate = 0;
    bool      bisection = false;
    bool      converged = false;
    real_type s(0), u(0); // to mute warning
    while ( ! (converged||bisection) ) {
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
      converged = abs(dp) <= abs(x) * machepsi; // Newton convergence indicator
    }
    if ( bisection ) {
      t = u - s; // initial bisection interval
      while ( abs(t) > abs(x) * machepsi ) { // bisection iterates
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/
  void
  solveQuadratic( real_type   a,
                  real_type   b,
                  real_type   c,
                  real_type & r1,
                  real_type & r2,
                  int_type  & nr,
                  int_type  & nc ) {
    r1 = r2 = 0;
    nr = nc = 0;
    if ( isZero(a) ) { // less than two roots b*z + c = 0
      if ( !isZero(b) ) { nr = 1; r1 = -c/b; }
    } else if ( isZero(c) ) { // a*z^2 + b*z  = 0
      nr = 2;
      r1 = -b/a;
      if ( r1 > 0 ) swap(r1,r2);
    } else { // Compute discriminant avoiding overflow.
      b /= 2; // b now b/2
      real_type abs_b = abs(b);
      real_type abs_c = abs(c);
      real_type e, d;
      if ( abs_b < abs_c ) {
        e = c < 0 ? -a : a;
        e = b*(b/abs_c) - e;
        d = sqrt(abs(e))*sqrt(abs_c);
      } else {
        e = 1 - (a/b)*(c/b);
        d = sqrt(abs(e))*abs_b;
      }
      bool real_root = e >= 0;
      if ( real_root ) {       // complex conjugate zeros
        if ( b >= 0 ) d = -d; // real zeros
        r1 = (d-b)/a;
        if ( !isZero(r1) ) {
          r2 = (c/r1)/a;
          if ( r1 > r2 ) swap(r1,r2); // order roots
        }
        nr = 2;
      } else {
        r1 = -b/a;          // real part
        r2 = abs(d/a); // immaginary part
        nc = 2;
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  static
  inline
  real_type
  guess1( real_type const a[3] ) {
    real_type const p =  1.09574;
    real_type const q = -3.239E-1;
    real_type const r = -3.239E-1;
    real_type const s =  9.57439E-2;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2];
  }

  static
  inline
  real_type
  guess2( real_type const a[3] ) {
    real_type const p = -1.09574;
    real_type const q =  3.239E-1;
    real_type const r = -3.239E-1;
    real_type const s =  9.57439E-2;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2];
  }

  static
  inline
  real_type
  guess3( real_type const a[3] ) {
    real_type const p =  1.14413;
    real_type const q = -2.75509E-1;
    real_type const r = -4.45578E-1;
    real_type const s = -2.59342E-2;
    real_type t = a[2]/3;
    if ( a[0] < t*(2*t*t-1) ) return  p+q*a[0]+r*a[2]+s*a[0]*a[2];
    else                      return -p+q*a[0]+r*a[2]-s*a[0]*a[2];
  }

  static
  inline
  real_type
  guess4( real_type const a[3] ) {
    real_type const q = -7.71845E-1;
    real_type const s = -2.28155E-1;
    if ( a[0] > 0 ) return (q+s*a[2])*a[0];
    else            return (q-s*a[2])*a[0];
  }

  static
  inline
  real_type
  guess5( real_type const a[3] ) {
    real_type p, q, r, s;
    real_type tmp = two27th-a[1]/3;
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
  inline
  real_type
  guess6( real_type const a[3] ) {
    real_type p, q, r, s;
    real_type tmp = a[1]/3-two27th;
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

  /*\
  ... Calculate the zeros of the cubic A*z^3 + B*z^2 + C*z + D.
  ...
  ... N. FLOCKE, Flash Center for Computational Science, University of Chicago
  ... Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  ... for Physical Applications
  ... ACM Transactions on Mathematical Software, Vol. 41, No. 4, 2015.
  ... DOI: http://dx.doi.org/10.1145/2699468
  \*/

  int_type
  solveCubic( real_type   A,
              real_type   B,
              real_type   C,
              real_type   D,
              real_type & r1,
              real_type & r2,
              real_type & r3,
              int_type  & nr,
              int_type  & nc ) {

    // special cases
    if ( isZero(A) ) {
      solveQuadratic( B, C, D, r1, r2, nr, nc );
      return 0;
    }
    if ( isZero(D) ) {
      r1 = 0;
      solveQuadratic( A, B, C, r2, r3, nr, nc );
      if ( nr == 1 ) { // caso degenere
        if ( r1 > r2 ) swap(r1,r2);
      } else if ( nr == 2 ) {
        if ( r1 > r2 ) swap(r1,r2);
        if ( r2 > r3 ) swap(r2,r3);
      }
      ++nr;
      return 0;
    }

    real_type scale, a[3];
    int_type  i_case;
    scaleCubicMonicPolynomial( B/A, C/A, D/A, a[2], a[1], a[0], i_case, scale );

    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    int_type iclass = -1;
    switch ( i_case ) {
      case 0: iclass = a[0] > 0 ? 2 : 1; break;
      case 1: iclass = a[1] > 0 ? 4 : 3; break;
      case 2: iclass = a[2] > 0 ? 6 : 5; break;
    }
    bool use_shifted = false;
    bool triple_root = false;
    switch ( iclass ) {
      case 1: r1 = guess1(a); break;
      case 2: r1 = guess2(a); break;
      case 3: r1 = guess3(a); break;
      case 4: r1 = guess4(a); break;
      case 5:
        r2 = a[1]-third;
        r3 = a[0]+one27th;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi;
        r1 = guess5(a);
        break;
      case 6:
        r2 = a[1]-third;
        r3 = a[0]-one27th;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi;
        r1 = guess6(a);
        break;
    }
    int_type iter = 0;
    if ( triple_root ) {
      nr = 3;
      if ( iclass == 5 ) r1 = r2 = r3 = -third * scale;
      else               r1 = r2 = r3 =  third * scale;
      return iter;
    } else if ( use_shifted ) {
      if ( iclass == 5 ) {
        // y^3 + A * y + (B+A/3), y = x-1/3
        // B = a[0]+1./27.;
        r1 -= third;
        r3 += third * r2;
        //if ( abs(r3) < machepsi ) r3 = 0;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 );
        r1 += third;
      } else {
        // y^3 + A * y + (B-A/3), y = x+1/3
        // B = a[0]-1./27.;
        r1 += third;
        r3 -= third * r2;
        //if ( abs(r3) < machepsi ) r3 = 0;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 );
        r1 -= third;
      }
    } else {
      iter = zeroCubicByNewtonBisection( a[2], a[1], a[0], r1 );
    }
    // scale
    r1 *= scale;
    
    real_type p  = ((A*r1+B)*r1+C)*r1+D;
    real_type dp = ((4*A*r1+3*B)*r1+2*C)*r1;
    
    r1 -= p/dp;
/*
    real_type const pp[]  = { A, B, C, D };
    real_type pH = CompHorner( pp, 3, r1, true );
    cout << "pH = " << pH << "\n";

    // una extra correzione con Newton dopo riscalatura
    real_type const dpp[] = { 3*A, 2*B, C };
    for ( int k = 0; k < 10; ++k ) {
      real_type pH  = CompHorner( pp, 3, r1, true );
      real_type dpH = (3*A*r1+2*B)*r1+C;
      r1 -= pH/dpH;
    }

    pH = CompHorner( pp, 3, r1, true );
    cout << "pH = " << pH << "\n";
*/
    // deflate
    real_type b0, b1;
    deflateCubicPolynomial( A, B, C, D, r1, b1, b0 );
    solveQuadratic( A, b1, b0, r2, r3, nr, nc );
    if ( nr == 2 ) { // if real roots sort it!
      if ( r1 > r2 ) swap(r1,r2);
      if ( r2 > r3 ) swap(r2,r3);
    }
    ++nr; // one more real root
    return iter;
  }

}

// EOF: PolynomialRoots.cc
