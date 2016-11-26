/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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

namespace PolynomialRoots {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^3 + A x^2 + B x + C
  static
  inline
  void
  scaleCubicMonicPolynomial( valueType   A,
                             valueType   B,
                             valueType   C,
                             valueType & AS,
                             valueType & BS,
                             valueType & CS,
                             indexType & i_case,
                             valueType & scale ) {

    valueType a = abs(A) ;
    valueType b = sqrt(abs(B)) ;
    valueType c = cbrt(abs(C)) ;

    if ( a < b ) {
      if ( b < c ) i_case = 0 ; // a < b < c --> c MAX
      else         i_case = 1 ; // a < b and c <= b --> b MAX
    } else {
      if ( a < c ) i_case = 0 ; // b <= a < c --> c MAX
      else         i_case = 2 ; // b <= a  and c <= a --> a MAX
    }

    switch ( i_case ) {
      case 0:
        scale = c ;
        AS    = A/c ;
        BS    = (B/c)/c ;
        CS    = C > 0 ? 1 : -1 ;
      break ;
      case 1:
        scale = b ;
        AS    = A/b ;
        BS    = B > 0 ? 1 : -1 ;
        CS    = ((C/b)/b)/b ;
      break ;
      case 2:
        scale = a ;
        AS    = A > 0 ? 1 : -1 ;
        BS    = (B/a)/a ;
        CS    = ((C/a)/a)/a ;
      break ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a3*x^2+b1*x+b0)
  static
  void
  deflateCubicPolynomial( valueType   a3,
                          valueType   a2,
                          valueType   a1,
                          valueType   a0,
                          valueType   r,
                          valueType & b1,
                          valueType & b0 ) {
    indexType i_cross  = 0 ;
    valueType r2       = r*r ;
    valueType v_cross  = abs(a0) ;
    valueType v_cross1 = abs(a1*r) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 1 ; }
    v_cross1 = abs(a2*r2) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 2 ; }
    v_cross1 = abs(a3*r*r2) ;
    if ( v_cross1 > v_cross ) i_cross = 3 ;
    switch ( i_cross ) {
      case 0: b1 = a2+a3*r ; b0 = a1+r*b1 ; break;
      case 1: b1 = a2+a3*r ; b0 = -a0/r   ; break;
      case 2:
      case 3: b0 = -a0/r ; b1 = (b0-a1)/r ; break;
    }
  }

  // x^3 + a*x^2 + b*x + c
  static
  inline
  valueType
  evalMonicCubic( valueType x,
                  valueType a,
                  valueType b,
                  valueType c ) {
    return ((x+a)*x+b)*x+c ;
  }

  static
  inline
  void
  evalMonicCubic( valueType   x,
                  valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & p,
                  valueType & dp ) {
    p  = x + a ;
    dp = x + p ;
    p  = p  * x + b ;
    dp = dp * x + p ;
    p  = p  * x + c ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  static
  indexType
  zeroCubicByNewtonBisection( valueType const a,
                              valueType const b,
                              valueType const c,
                              valueType     & x ) {

    valueType p, dp ;
    evalMonicCubic( x, a, b, c, p, dp ) ;
    valueType t = p ; // save p(x) for sign comparison
    x -= p/dp ; // 1st improved root

    indexType iter      = 1 ;
    indexType oscillate = 0 ;
    bool      bisection = false ;
    bool      converged = false ;
    valueType s(0), u(0) ; // to mute warning
    while ( ! (converged||bisection) ) {
      ++iter ;
      evalMonicCubic( x, a, b, c, p, dp ) ;
      if ( p*t < 0 ) { // does Newton start oscillating ?
        if ( p < 0 ) {
          ++oscillate ; // increment oscillation counter
          s = x ;       // save lower bisection bound
        } else {
          u = x ; // save upper bisection bound
        }
        t = p ; // save current p(x)
      }
      dp = p/dp ; // Newton correction
      x -= dp ;   // new Newton root
      bisection = oscillate > 2 ; // activate bisection
      converged = abs(dp) <= abs(x) * machepsi ; // Newton convergence indicator
    }
    if ( bisection ) {
      t = u - s ; // initial bisection interval
      while ( abs(t) > abs(x) * machepsi ) { // bisection iterates
        ++iter ;
        p = evalMonicCubic( x, a, b, c ) ;
        if ( p < 0 ) s = x ;
        else         u = x ; // keep bracket on root
        t = (u-s)/2 ; // new bisection interval
        x = s + t ;   // new bisection root
      }
    }
    return iter ;
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
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2,
                  indexType & nr,
                  indexType & nc ) {
    r1 = r2 = 0 ;
    nr = nc = 0 ;
    if ( a == 0 ) { // less than two roots b*z + c = 0
      if ( b != 0 ) { nr = 1 ; r1 = -c/b ; }
    } else if ( c == 0 ) { // a*z^2 + b*z  = 0
      nr = 2 ;
      r1 = -b/a ;
      if ( r1 > 0 ) std::swap(r1,r2) ;
    } else { // Compute discriminant avoiding overflow.
      b /= 2 ; // b now b/2
      valueType abs_b = abs(b) ;
      valueType abs_c = abs(c) ;
      valueType e, d ;
      if ( abs_b < abs_c ) {
        e = c < 0 ? -a : a ;
        e = b*(b/abs_c) - e ;
        d = sqrt(abs(e))*sqrt(abs_c);
      } else {
        e = 1 - (a/b)*(c/b);
        d = sqrt(abs(e))*abs_b ;
      }
      bool real_root = e >= 0 ;
      if ( real_root ) {       // complex conjugate zeros
        if ( b >= 0 ) d = -d ; // real zeros
        r1 = (d-b)/a;
        if ( r1 != 0 ) {
          r2 = (c/r1)/a ;
          if ( r1 > r2 ) std::swap(r1,r2) ; // order roots
        }
        nr = 2 ;
      } else {
        r1 = -b/a ;          // real part
        r2 = std::abs(d/a) ; // immaginary part
        nc = 2 ;
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  static
  inline
  valueType
  guess1( valueType const a[3] ) {
    valueType const p =  1.09574 ;
    valueType const q = -3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess2( valueType const a[3] ) {
    valueType const p = -1.09574 ;
    valueType const q =  3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess3( valueType const a[3] ) {
    valueType const p =  1.14413    ;
    valueType const q = -2.75509E-1 ;
    valueType const r = -4.45578E-1 ;
    valueType const s = -2.59342E-2 ;
    valueType t = a[2]/3 ;
    if ( a[0] < t*(2*t*t-1) ) return  p+q*a[0]+r*a[2]+s*a[0]*a[2] ;
    else                      return -p+q*a[0]+r*a[2]-s*a[0]*a[2] ;
  }

  static
  inline
  valueType
  guess4( valueType const a[3] ) {
    valueType const q = -7.71845E-1 ;
    valueType const s = -2.28155E-1 ;
    if ( a[0] > 0 ) return (q+s*a[2])*a[0] ;
    else            return (q-s*a[2])*a[0] ;
  }

  static
  inline
  valueType
  guess5( valueType const a[3] ) {
    valueType p, q, r, s ;
    valueType tmp = two27th-a[1]/3 ;
    if ( a[1] <= third ) {
      if ( a[0] < tmp ) {
        p =  8.78558E-1 ;
        q = -5.71888E-1 ;
        r = -7.11154E-1 ;
        s = -3.22313E-1 ;
      } else {
        p = -1.92823E-1 ;
        q = -5.66324E-1 ;
        r = +5.05734E-1 ;
        s = -2.64881E-1 ;
      }
    } else {
      if ( a[0] < tmp ) {
        p = 1.19748 ;
        q = -2.83772E-1 ;
        r = -8.37476E-1 ;
        s = -3.56228E-1 ;
      } else {
        p = -3.45219E-1 ;
        q = -4.01231E-1 ;
        r =  2.07216E-1 ;
        s = -4.45532E-3 ;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1] ;
  }

  static
  inline
  valueType
  guess6( valueType const a[3] ) {
    valueType p, q, r, s ;
    valueType tmp = a[1]/3-two27th ;
    if ( a[1] <= third ) {
      if ( a[0] > tmp ) {
        p = -8.78558E-1 ;
        q = -5.71888E-1 ;
        r =  7.11154E-1 ;
        s = -3.22313E-1 ;
      } else {
        p =  1.92823E-1 ;
        q = -5.66324E-1 ;
        r = -5.05734E-1 ;
        s = -2.64881E-1 ;
      }
    } else {
      if ( a[0] > tmp ) {
        p = -1.19748 ;
        q = -2.83772E-1 ;
        r =  8.37476E-1 ;
        s = -3.56228E-1 ;
      } else {
        p =  3.45219E-1 ;
        q = -4.01231E-1 ;
        r = -2.07216E-1 ;
        s = -4.45532E-3 ;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1] ;
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

  indexType
  solveCubic( valueType   A,
              valueType   B,
              valueType   C,
              valueType   D,
              valueType & r1,
              valueType & r2,
              valueType & r3,
              indexType & nr,
              indexType & nc ) {

    // special cases
    if ( A == 0 ) {
      solveQuadratic( B, C, D, r1, r2, nr, nc ) ;
      return 0 ;
    }
    if ( D == 0 ) {
      r1 = 0 ;
      solveQuadratic( A, B, C, r2, r3, nr, nc ) ;
      if ( nr == 1 ) { // caso degenere
        if ( r1 > r2 ) std::swap(r1,r2) ;
      } else if ( nr == 2 ) {
        if ( r1 > r2 ) std::swap(r1,r2) ;
        if ( r2 > r3 ) std::swap(r2,r3) ;
      }
      ++nr ;
      return 0 ;
    }

    valueType scale, a[3] ;
    indexType i_case ;
    scaleCubicMonicPolynomial( B/A, C/A, D/A, a[2], a[1], a[0], i_case, scale ) ;

    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    indexType iclass = -1 ;
    switch ( i_case ) {
      case 0: iclass = a[0] > 0 ? 2 : 1 ; break ;
      case 1: iclass = a[1] > 0 ? 4 : 3 ; break ;
      case 2: iclass = a[2] > 0 ? 6 : 5 ; break ;
    }
    bool use_shifted = false ;
    bool triple_root = false ;
    switch ( iclass ) {
      case 1: r1 = guess1(a) ; break ;
      case 2: r1 = guess2(a) ; break ;
      case 3: r1 = guess3(a) ; break ;
      case 4: r1 = guess4(a) ; break ;
      case 5:
        r2 = a[1]-third ;
        r3 = a[0]+one27th ;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01 ;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi ;
        r1 = guess5(a) ;
        break ;
      case 6:
        r2 = a[1]-third ;
        r3 = a[0]-one27th ;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01 ;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi ;
        r1 = guess6(a) ;
        break ;
    }
    indexType iter = 0 ;
    if ( triple_root ) {
      nr = 3 ;
      if ( iclass == 5 ) r1 = r2 = r3 = -third * scale ;
      else               r1 = r2 = r3 =  third * scale ;
      return iter ;
    } else if ( use_shifted ) {
      if ( iclass == 5 ) {
        // y^3 + A * y + (B+A/3), y = x-1/3
        // B = a[0]+1./27. ;
        r1 -= third ;
        r3 += third * r2 ;
        //if ( abs(r3) < machepsi ) r3 = 0 ;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 ) ;
        r1 += third ;
      } else {
        // y^3 + A * y + (B-A/3), y = x+1/3
        // B = a[0]-1./27. ;
        r1 += third ;
        r3 -= third * r2 ;
        //if ( abs(r3) < machepsi ) r3 = 0 ;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 ) ;
        r1 -= third ;
      }
    } else {
      iter = zeroCubicByNewtonBisection( a[2], a[1], a[0], r1 ) ;
    }
    // scale
    r1 *= scale ;
    
    valueType p  = ((A*r1+B)*r1+C)*r1+D ;
    valueType dp = ((4*A*r1+3*B)*r1+2*C)*r1 ;
    
    r1 -= p/dp ;
/*
    valueType const pp[]  = { A, B, C, D } ;
    valueType pH = CompHorner( pp, 3, r1, true ) ;
    std::cout << "pH = " << pH << "\n" ;

    // una extra correzione con Newton dopo riscalatura
    valueType const dpp[] = { 3*A, 2*B, C } ;
    for ( int k = 0 ; k < 10 ; ++k ) {
      valueType pH  = CompHorner( pp, 3, r1, true ) ;
      valueType dpH = (3*A*r1+2*B)*r1+C ;
      r1 -= pH/dpH ;
    }

    pH = CompHorner( pp, 3, r1, true ) ;
    std::cout << "pH = " << pH << "\n" ;
*/
    // deflate
    valueType b0, b1 ;
    deflateCubicPolynomial( A, B, C, D, r1, b1, b0 ) ;
    solveQuadratic( A, b1, b0, r2, r3, nr, nc ) ;
    if ( nr == 2 ) { // if real roots sort it!
      if ( r1 > r2 ) std::swap(r1,r2) ;
      if ( r2 > r3 ) std::swap(r2,r3) ;
    }
    ++nr ; // one more real root
    return iter ;
  }

}

// EOF: PolynomialRoots.cc
