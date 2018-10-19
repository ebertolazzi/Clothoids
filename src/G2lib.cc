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

#include "G2lib.hh"

#include <algorithm>

namespace G2lib {

  real_type const machepsi     = std::numeric_limits<real_type>::epsilon();
  real_type const machepsi10   = 10*machepsi;
  real_type const machepsi100  = 100*machepsi;
  real_type const machepsi1000 = 1000*machepsi;
  real_type const m_pi         = 3.14159265358979323846264338328;  // pi
  real_type const m_pi_2       = 1.57079632679489661923132169164;  // pi/2
  real_type const m_2pi        = 6.28318530717958647692528676656;  // 2*pi
  real_type const m_1_pi       = 0.318309886183790671537767526745; // 1/pi
  real_type const m_1_sqrt_pi  = 0.564189583547756286948079451561; // 1/sqrt(pi)

  void
  rangeSymm( real_type & ang ) {
    ang = fmod( ang, m_2pi );
    while ( ang < -m_pi ) ang += m_2pi;
    while ( ang >  m_pi ) ang -= m_2pi;
  }

  static
  inline
  real_type
  power2( real_type a )
  { return a*a; }

  static
  inline
  real_type
  power3( real_type a )
  { return a*a*a; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // sin(x)/x
  */
  real_type
  Sinc( real_type x ) {
    if ( std::abs(x) < 0.02 ) {
      real_type x2 = x*x;
      return 1-(x2/6)*(1-(x2/20)*(1-x2/42));
    } else {
      return sin(x)/x;
    }
  }

  real_type
  Sinc_D( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.04 ) return -(x/3)*(1-(x2/10)*(1-(x2/28)*(1-(x2/54))));
    else                      return (cos(x)-sin(x)/x)/x;
  }

  real_type
  Sinc_DD( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.02 ) return -1./3.+x2*(0.1-x2*((1.0/168.0)-(x2/6480)));
    else                      return ((2/x2-1)*sin(x)-2*cos(x)/x)/x;
  }

  real_type
  Sinc_DDD( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.009 ) return (1.0/5.0+(-1.0/42.0+(1.0/1080.0)*x2)*x2)*x;
    else                       return ((6/x2-1)*cos(x)+(3-6/x2)*sin(x)/x)/x;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // (1-cos(x))/x
  */
  real_type
  Cosc( real_type x ) {
    if ( std::abs(x) < 0.04 ) {
      real_type x2 = x*x;
      return (x/2)*(1-(x2/12)*(1-(x2/30)*(1-x2/56)));
    } else {
      return (1-cos(x))/x;
    }
  }

  real_type
  Cosc_D( real_type x ) {
    if ( std::abs(x) < 0.02 ) {
      real_type x2  = x*x;
      return 0.5*(1-(x2/4)*(1-(x2/18)*(1-(x2/40))));
    } else {
      return (sin(x)+(cos(x)-1)/x)/x;
    }
  }

  real_type
  Cosc_DD( real_type x ) {
    real_type x2  = x*x;
    if ( std::abs(x) < 0.04 ) return -(x/4)*(1-(x2/9)*(1-((3.0/80.0)*x2)*(1-((2.0/105.0)*x2))));
    else                      return ((1-2/x2)*cos(x)+(2/x-sin(x))/x)/x;
  }

  real_type
  Cosc_DDD( real_type x ) {
    real_type x2  = x*x;
    if ( std::abs(x) < 0.02 ) return -(1-(x2/3)*(1-(x2/16)*(1-(2.0/75.0)*x2)))/4.0;
    else                      return ((6/x2-1)*sin(x)+((6/x2-3)*cos(x)-6/x2)/x)/x;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // atan(x)/x
  */
  real_type
  Atanc( real_type x ) {
    if ( std::abs(x) < 0.03 ) {
      real_type x2 = x*x;
      return 1-x2*((1./3.)-x2*((1./5.)-x2*((1./7.)-x2*((1./9.)-(x2/11)))));
    } else {
      return atan(x)/x;
    }
  }

  real_type
  Atanc_D( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.03 ) {
      return x*( -(2./3.) + x2*( (4./5.) + x2*( -(6./7.) + x2*( (8./9.) + x2*( -(10./11.) + x2*(12./13.))))));
    } else {
      return (1/(1+x2)-atan(x)/x)/x;
    }
  }

  real_type
  Atanc_DD( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.02 ) {
      return -2./3.+ x2*( (12./5.) + (-(30./7.) + x2 * ( (56./9.) + x2*( -(90./11.) + x2 * (132./13.)))));
    } else {
      return (2*atan(x)/x-(4*x2+2)/power2(1+x2))/x2;
    }
  }

  real_type
  Atanc_DDD( real_type x ) {
    real_type x2 = x*x;
    if ( std::abs(x) < 0.02 ) {
      return x*(24./5.+x2*(-120./7. + x2 * (112./3. + x2 * (-720./11. + x2*(1320./13. - x2*728./5.)))));
    } else {
      return ( ((18*x2+16)*x2+6)/power3(x2+1)-6*atan(x)/x )/(x2*x);
    }
  }

  /*\
   |   ____        _           ____       ____
   |  / ___|  ___ | |_   _____|___ \__  _|___ \
   |  \___ \ / _ \| \ \ / / _ \ __) \ \/ / __) |
   |   ___) | (_) | |\ V /  __// __/ >  < / __/
   |  |____/ \___/|_| \_/ \___|_____/_/\_\_____|
  \*/

  bool
  Solve2x2::factorize( real_type A[2][2] ) {
    // full pivoting
    real_type Amax = std::abs(A[0][0]);
    real_type tmp  = std::abs(A[0][1]);
    int_type ij = 0;
    if ( tmp > Amax ) { ij = 1; Amax = tmp; }
    tmp = std::abs(A[1][0]);
    if ( tmp > Amax ) { ij = 2; Amax = tmp; }
    tmp = std::abs(A[1][1]);
    if ( tmp > Amax ) { ij = 3; Amax = tmp; }
    if ( isZero(Amax) ) return false;
    if ( (ij&0x01) == 0x01 ) { j[0] = 1; j[1] = 0; }
    else                     { j[0] = 0; j[1] = 1; }
    if ( (ij&0x02) == 0x02 ) { i[0] = 1; i[1] = 0; }
    else                     { i[0] = 0; i[1] = 1; }
    // apply factorization
    LU[0][0] = A[i[0]][j[0]];
    LU[0][1] = A[i[0]][j[1]];
    LU[1][0] = A[i[1]][j[0]];
    LU[1][1] = A[i[1]][j[1]];

    LU[1][0] /= LU[0][0];
    LU[1][1] -= LU[1][0]*LU[0][1];
    // check for singularity
    singular = std::abs( LU[1][1] ) < epsi;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Solve2x2::solve( real_type const b[2], real_type x[2] ) const {
    if ( singular ) {
      // L^+ Pb
      real_type tmp = (b[i[0]] + LU[1][0]*b[i[1]]) /
                      ( (1+power2(LU[1][0]) ) * ( power2(LU[0][0])+power2(LU[0][1]) ) );
      x[j[0]] = tmp*LU[0][0];
      x[j[1]] = tmp*LU[0][1];
      // check consistency
      tmp = (LU[0][0]*x[j[0]]+LU[0][1]*x[j[1]]);
      return hypot( b[i[0]]-tmp, b[i[1]]+tmp*LU[1][0] ) < hypot(b[0],b[1])*epsi;
    } else { // non singular
      // L^(-1) Pb
      x[j[0]] = b[i[0]];
      x[j[1]] = b[i[1]]-LU[1][0]*x[j[0]];
      // U^(-1) x
      x[j[1]] /= LU[1][1];
      x[j[0]]  = (x[j[0]]-LU[0][1]*x[j[1]])/LU[0][0];
      return FP_INFINITE != std::fpclassify(x[0]) &&
             FP_NAN      != std::fpclassify(x[0]) &&
             FP_INFINITE != std::fpclassify(x[1]) &&
             FP_NAN      != std::fpclassify(x[1]);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  isCounterClockwise( real_type const P1[2],
                      real_type const P2[2],
                      real_type const P3[2] ) {
    real_type dx1 = P2[0] - P1[0];
    real_type dy1 = P2[1] - P1[1];
    real_type dx2 = P3[0] - P1[0];
    real_type dy2 = P3[1] - P1[1];
    real_type tol = machepsi10*(hypot(dx1,dy1)*hypot(dx2,dy2));
    real_type det = dx1*dy2 - dy1*dx2;
    if      ( det >  tol ) return  1;
    else if ( det < -tol ) return -1;
    return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  isPointInTriangle( real_type const point[2],
                     real_type const p1[2],
                     real_type const p2[2],
                     real_type const p3[2] ) {
    int_type d = isCounterClockwise(p1, p2, p3);
    int_type a = isCounterClockwise(p1, p2, point);
    int_type b = isCounterClockwise(p2, p3, point);
    int_type c = isCounterClockwise(p3, p1, point);
    if ( d < 0) { a = -a; b = -b; c = -c; }
    if ( a < 0 ) return -1;
    if ( b < 0 ) return -1;
    if ( c < 0 ) return -1;
    if ( a+b+c == 3 ) return 1;
    return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  updateInterval( int_type      & lastInterval,
                  real_type       x,
                  real_type const Xvec[],
                  int_type        npts ) {

    if ( npts <= 2 ) { lastInterval = 0; return; } // nothing to search

    // optimized interval search
    real_type const * XL = Xvec + lastInterval;
    if ( XL[1] <= x ) { // x on the right
      if ( x >= Xvec[npts-2] ) { // x in [X[npt-2],X[npts-1]]
        lastInterval = npts-2; // last interval
      } else if ( x < XL[2] ) { // x in [XL[1],XL[2])
        ++lastInterval;
      } else { // x >= XL[2] search the right interval
        real_type const * XE = Xvec+npts;
        lastInterval += int_type(std::lower_bound( XL, XE, x )-XL);
        if ( Xvec[lastInterval] > x ) --lastInterval;
      }
    } else if ( x < XL[0] ) { // on the left
      if ( x < Xvec[1] ) { // x in [X[0],X[1])
        lastInterval = 0; // first interval
      } else if ( XL[-1] <= x ) { // x in [XL[-1],XL[0])
        --lastInterval;
      } else {
        lastInterval = int_type(std::lower_bound( Xvec, XL, x )-Xvec);
        if ( Xvec[lastInterval] > x ) --lastInterval;
      }
    } else {
      // x in the interval [XL[0],XL[1]) nothing to do
    }
  }

}

// EOF: G2lib.cc
