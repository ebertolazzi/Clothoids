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

namespace G2lib {

  valueType const machepsi    = std::numeric_limits<valueType>::epsilon() ;
  valueType const m_pi        = 3.14159265358979323846264338328  ; // pi
  valueType const m_pi_2      = 1.57079632679489661923132169164  ; // pi/2
  valueType const m_2pi       = 6.28318530717958647692528676656  ; // 2*pi
  valueType const m_1_pi      = 0.318309886183790671537767526745 ; // 1/pi
  valueType const m_1_sqrt_pi = 0.564189583547756286948079451561 ; // 1/sqrt(pi)

  void
  rangeSymm( valueType & ang ) {
    ang = fmod( ang, m_2pi ) ;
    while ( ang < -m_pi ) ang += m_2pi ;
    while ( ang >  m_pi ) ang -= m_2pi ;
  }

  static
  inline
  valueType
  power2( valueType a )
  { return a*a ; }

  static
  inline
  valueType
  power3( valueType a )
  { return a*a*a ; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // sin(x)/x
  */
  valueType
  Sinc( valueType x ) {
    if ( std::abs(x) < 0.02 ) {
      valueType x2 = x*x ;
      return 1-(x2/6)*(1-(x2/20)*(1-x2/42));
    } else {
      return sin(x)/x;
    }
  }

  valueType
  Sinc_D( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.04 ) return -(x/3)*(1-(x2/10)*(1-(x2/28)*(1-(x2/54)))) ;
    else                      return (cos(x)-sin(x)/x)/x ;
  }

  valueType
  Sinc_DD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) return -1./3.+x2*(0.1-x2*((1.0/168.0)-(x2/6480))) ;
    else                      return ((2/x2-1)*sin(x)-2*cos(x)/x)/x ;
  }

  valueType
  Sinc_DDD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.009 ) return (1.0/5.0+(-1.0/42.0+(1.0/1080.0)*x2)*x2)*x;
    else                       return ((6/x2-1)*cos(x)+(3-6/x2)*sin(x)/x)/x ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // (1-cos(x))/x
  */
  valueType
  Cosc( valueType x ) {
    if ( std::abs(x) < 0.04 ) {
      valueType x2 = x*x ;
      return (x/2)*(1-(x2/12)*(1-(x2/30)*(1-x2/56)));
    } else {
      return (1-cos(x))/x;
    }
  }

  valueType
  Cosc_D( valueType x ) {
    if ( std::abs(x) < 0.02 ) {
      valueType x2  = x*x;
      return 0.5*(1-(x2/4)*(1-(x2/18)*(1-(x2/40))));
    } else {
      return (sin(x)+(cos(x)-1)/x)/x ;
    }
  }

  valueType
  Cosc_DD( valueType x ) {
    valueType x2  = x*x;
    if ( std::abs(x) < 0.04 ) return -(x/4)*(1-(x2/9)*(1-((3.0/80.0)*x2)*(1-((2.0/105.0)*x2))));
    else                      return ((1-2/x2)*cos(x)+(2/x-sin(x))/x)/x ;
  }

  valueType
  Cosc_DDD( valueType x ) {
    valueType x2  = x*x;
    if ( std::abs(x) < 0.02 ) return -(1-(x2/3)*(1-(x2/16)*(1-(2.0/75.0)*x2)))/4.0 ;
    else                      return ((6/x2-1)*sin(x)+((6/x2-3)*cos(x)-6/x2)/x)/x ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  // atan(x)/x
  */
  valueType
  Atanc( valueType x ) {
    if ( std::abs(x) < 0.03 ) {
      valueType x2 = x*x ;
      return 1-x2*((1./3.)-x2*((1./5.)-x2*((1./7.)-x2*((1./9.)-(x2/11))))) ;
    } else {
      return atan(x)/x;
    }
  }
  valueType
  Atanc_D( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.03 ) {
      return x*( -(2./3.) + x2*( (4./5.) + x2*( -(6./7.) + x2*( (8./9.) + x2*( -(10./11.) + x2*(12./13.)))))) ;
    } else {
      return (1/(1+x2)-atan(x)/x)/x ;
    }
  }

  valueType
  Atanc_DD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) {
      return -2./3.+ x2*( (12./5.) + (-(30./7.) + x2 * ( (56./9.) + x2*( -(90./11.) + x2 * (132./13.))))) ;
    } else {
      return (2*atan(x)/x-(4*x2+2)/power2(1+x2))/x2 ;
    }
  }

  valueType
  Atanc_DDD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) {
      return x*(24./5.+x2*(-120./7. + x2 * (112./3. + x2 * (-720./11. + x2*(1320./13. - x2*728./5.))))) ;
    } else {
      return ( ((18*x2+16)*x2+6)/power3(x2+1)-6*atan(x)/x )/(x2*x) ;
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
  Solve2x2::factorize( valueType A[2][2] ) {
    // full pivoting
    valueType Amax = std::abs(A[0][0]) ;
    valueType tmp  = std::abs(A[0][1]) ;
    indexType ij = 0 ;
    if ( tmp > Amax ) { ij = 1 ; Amax = tmp ; }
    tmp = std::abs(A[1][0]) ;
    if ( tmp > Amax ) { ij = 2 ; Amax = tmp ; }
    tmp = std::abs(A[1][1]) ;
    if ( tmp > Amax ) { ij = 3 ; Amax = tmp ; }
    if ( Amax == 0 ) return false ;
    if ( (ij&0x01) == 0x01 ) { j[0] = 1 ; j[1] = 0 ; }
    else                     { j[0] = 0 ; j[1] = 1 ; }
    if ( (ij&0x02) == 0x02 ) { i[0] = 1 ; i[1] = 0 ; }
    else                     { i[0] = 0 ; i[1] = 1 ; }
    // apply factorization
    LU[0][0] = A[i[0]][j[0]] ;
    LU[0][1] = A[i[0]][j[1]] ;
    LU[1][0] = A[i[1]][j[0]] ;
    LU[1][1] = A[i[1]][j[1]] ;

    LU[1][0] /= LU[0][0] ;
    LU[1][1] -= LU[1][0]*LU[0][1] ;
    // check for singularity
    singular = std::abs( LU[1][1] ) < epsi ;
    return true ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Solve2x2::solve( valueType const b[2], valueType x[2] ) const {
    if ( singular ) {
      // L^+ Pb
      valueType tmp = (b[i[0]] + LU[1][0]*b[i[1]]) /
                      ( (1+power2(LU[1][0]) ) * ( power2(LU[0][0])+power2(LU[0][1]) ) ) ;
      x[j[0]] = tmp*LU[0][0] ;
      x[j[1]] = tmp*LU[0][1] ;
      // check consistency
      tmp = (LU[0][0]*x[j[0]]+LU[0][1]*x[j[1]]) ;
      return hypot( b[i[0]]-tmp, b[i[1]]+tmp*LU[1][0] ) < hypot(b[0],b[1])*epsi ;
    } else { // non singular
      // L^(-1) Pb
      x[j[0]] = b[i[0]] ;
      x[j[1]] = b[i[1]]-LU[1][0]*x[j[0]] ;
      // U^(-1) x
      x[j[1]] /= LU[1][1] ;
      x[j[0]]  = (x[j[0]]-LU[0][1]*x[j[1]])/LU[0][0] ;
      return FP_INFINITE != std::fpclassify(x[0]) &&
             FP_NAN      != std::fpclassify(x[0]) &&
             FP_INFINITE != std::fpclassify(x[1]) &&
             FP_NAN      != std::fpclassify(x[1]);
    }
  }

}

// EOF: G2lib.cc
