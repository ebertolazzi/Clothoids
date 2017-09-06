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

  /*
  // sin(x)/x
  */
  valueType
  Sinc( valueType x ) {
    if ( std::abs(x) < 0.002 ) return 1+x*(1/6-x*x/20);
    else                       return sin(x)/x;
  }

  valueType
  Sinc_D( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) return (-1.0/3.0+(1.0/30.0-(1.0/840.0)*x2)*x2)*x;
    else                      return (cos(x)*x-sin(x))/x2 ;
  }

  valueType
  Sinc_DD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.05 ) return -1.0/3.0+(1.0/10.0+(-1.0/168.0+(1.0/6480.0)*x2)*x2)*x2;
    else                      return ((2-x2)*sin(x)/x-2*cos(x))/x2 ;
  }

  valueType
  Sinc_DDD( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) return (1.0/5.0+(-1.0/42.0+(1.0/1080.0)*x2)*x2)*x;
    else                      return ((3*x2-6)*sin(x)+x*(6-x2)*cos(x))/(x2*x2) ;
  }

  /*
  // (1-cos(x))/x
  */
  valueType
  Cosc( valueType x ) {
    valueType x2 = x*x ;
    if ( std::abs(x) < 0.02 ) {
      return (1.0/2.0+(-.01/24.0+(1.0/720.0)*x2)*x2)*x;
    } else {
      return (1-cos(x))/x;
    }
  }

  valueType
  Cosc_D( valueType x ) {
    valueType x2  = x*x;
    if ( std::abs(x) < 0.05 ) return 1.0/2.0+(-1.0/8.0+(1.0/144.0-(1.0/5760.0)*x2)*x2)*x2;
    else                      return (sin(x)*x+cos(x)-1)/x2 ;
  }

  valueType
  Cosc_DD( valueType x ) {
    valueType x2  = x*x;
    if ( std::abs(x) < 0.02 ) return (-1.0/4.0+(1.0/36.0-(1.0/960.0)*x2)*x2)*x;
    else                      return ((2+(x2-2)*cos(x))/x-2*sin(x))/x2 ;
  }

  valueType
  Cosc_DDD( valueType x ) {
    valueType x2  = x*x;
    if ( std::abs(x) < 0.05 ) return -1.0/4.0+(1.0/12.0+(-1.0/192.0+(1.0/7200.0)*x2)*x2)*x2 ;
    else                      return ((6-x2)*sin(x)-6+(6-3*x2)*cos(x)/x)/(x2*x) ;
  }

  static
  inline
  valueType
  power2( valueType a )
  { return a*a ; }

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
