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

///
/// file: Arcs.cc
///

#include "Clothoid.hh"
#include "CubicRootsFlocke.hh"

#include <cmath>
#include <cfloat>
#include <sstream>
#include <stdexcept>

#ifndef CLOTHOID_ASSERT
  #define CLOTHOID_ASSERT(COND,MSG)         \
    if ( !(COND) ) {                        \
      std::ostringstream ost ;              \
      ost << "On line: " << __LINE__        \
          << " file: " << __FILE__          \
          << '\n' << MSG << '\n' ;          \
      throw std::runtime_error(ost.str()) ; \
    }
#endif

namespace Clothoid {

  using namespace std ;

  static const valueType m_pi = 3.14159265358979323846264338328  ; // pi

  /*\
   |    ____ _          _
   |   / ___(_)_ __ ___| | ___  ___
   |  | |   | | '__/ __| |/ _ \/ __|
   |  | |___| | | | (__| |  __/\__ \
   |   \____|_|_|  \___|_|\___||___/
  \*/

  static
  inline
  valueType
  powersub( valueType a, valueType b)
  { return (a+b)*(a-b); }

  /*
    http://www.lucidarme.me/?p=490
  */
  void
  CircleTangentPoints( valueType PA[2],
                       valueType rA,
                       valueType PB[2],
                       valueType rB,
                       bool &    external_tangents,
                       valueType PTE0[2][2],
                       valueType PTE1[2][2],
                       bool &    internal_tangents,
                       valueType PTI0[2][2],
                       valueType PTI1[2][2] ) {

    // Compute distance between circle centers
    valueType D = hypot( PB[0]-PA[0], PB[1]-PA[1] );

    // First case : process external tangents
    valueType disc1 = powersub(D,rB-rA) ;
    external_tangents = disc1 >= 0 ;

    if ( external_tangents ) {
      // Compute the lenght of the tangents
      valueType L = sqrt(disc1);

      // Compute the parameters
      valueType R1     = hypot(L,rB);
      valueType Sigma1 = 0.25 * sqrt( powersub( D+rA, R1 )* powersub( R1, D-rA ) );
      valueType R2     = hypot(L,rA);
      valueType Sigma2 = 0.25 * sqrt( powersub( D+rB, R2 )* powersub( R2, D-rB ) );

      valueType xM  = (PA[0]+PB[0])/2 ;
      valueType yM  = (PA[1]+PB[1])/2 ;
      valueType Tx  = 2*(PA[0]-PB[0])/(D*D) ;
      valueType Ty  = 2*(PA[1]-PB[1])/(D*D) ;
      valueType bf0 = powersub(rA,R1)/4;
      valueType bf1 = powersub(rB,R2)/4;

      valueType TxS1 = Tx*Sigma1 ;
      valueType TxS2 = Tx*Sigma2 ;
      valueType TyS1 = Ty*Sigma1 ;
      valueType TyS2 = Ty*Sigma2 ;

      valueType xM0 = xM - Tx*bf0 ;
      valueType yM0 = yM - Ty*bf0 ;
      valueType xM1 = xM + Tx*bf1 ;
      valueType yM1 = yM + Ty*bf1 ;

      // Compute the first tangent
      PTE0[0][0] = xM0 + TyS1;
      PTE0[0][1] = yM0 - TxS1;
      PTE0[1][0] = xM1 + TyS2;
      PTE0[1][1] = yM1 - TxS2;

      // Compute second tangent
      PTE1[0][0] = xM0 - TyS1;
      PTE1[0][1] = yM0 + TxS1;
      PTE1[1][0] = xM1 - TyS2;
      PTE1[1][1] = yM1 + TxS2;
    }

    // Second case : process internal tangents
    valueType disc2 = powersub(D,rB+rA) ;
    internal_tangents = disc2 >= 0 ;

    if ( internal_tangents ) {
      // Compute the lenght of the tangents
      valueType L = sqrt(disc2);

      // Compute the parameters
      valueType R1     = hypot(L,rB);
      valueType Sigma1 = 0.25 * sqrt ( powersub(D+rA,R1)*powersub(R1,D-rA) );
      valueType R2     = hypot(L,rA);
      valueType Sigma2 = 0.25 * sqrt ( powersub(D+rB,R2)*powersub(R2,D-rB) );

      valueType xM  = (PA[0]+PB[0])/2 ;
      valueType yM  = (PA[1]+PB[1])/2 ;
      valueType Tx  = 2*(PA[0]-PB[0])/(D*D) ;
      valueType Ty  = 2*(PA[1]-PB[1])/(D*D) ;
      valueType bf0 = powersub(rA,R1)/4;
      valueType bf1 = powersub(rB,R2)/4;

      valueType TxS1 = Tx*Sigma1;
      valueType TxS2 = Tx*Sigma2;
      valueType TyS1 = Ty*Sigma1;
      valueType TyS2 = Ty*Sigma2;

      valueType xM0 = xM - Tx*bf0;
      valueType yM0 = yM - Ty*bf0;
      valueType xM1 = xM + Tx*bf1;
      valueType yM1 = yM + Ty*bf1;

      // Compute the first tangent
      PTI0[0][0] = xM0 + TyS1;
      PTI0[0][1] = yM0 - TxS1;
      PTI0[1][0] = xM1 - TyS2;
      PTI0[1][1] = yM1 + TxS2;
    
      // Compute second tangent
      PTI0[0][0] = xM0 - TyS1;
      PTI0[0][1] = yM0 + TxS1;
      PTI0[1][0] = xM1 + TyS2;
      PTI0[1][1] = yM1 - TxS2;
    }
  }

  bool
  CircleLineTransition( valueType C[2],
                        valueType r,
                        valueType P[2],
                        valueType theta,
                        valueType C0[2],
                        valueType C1[2] ) {
    valueType Nx =  sin(theta) ;
    valueType Ny = -cos(theta) ;
    valueType Dx = C[0] - P[0] ;
    valueType Dy = C[1] - P[1] ;
    valueType delta = (Dx*Dx+Dy*Dy-r*r)/2 ;
    if ( delta <= 0 ) return false ;
    valueType s0 = delta/(Dx*Nx+Dy*Ny-r) ;
    C0[0] = P[0]+s0*Nx ;
    C0[1] = P[1]+s0*Ny ;
    valueType s1 = delta/(Dx*Nx+Dy*Ny+r) ;
    C1[0] = P[0]+s1*Nx ;
    C1[1] = P[1]+s1*Ny ;
    return true ;
  }

  indexType
  ArcToNURBS( valueType x0,
              valueType y0,
              valueType theta0,
              valueType x1,
              valueType y1,
              valueType knots[7],
              valueType Poly[4][3] ) {
    valueType dx = x1 - x0;
    valueType dy = y1 - y0;
    valueType alpha  = atan2(dy,dx);
    valueType L      = hypot(dy,dx);
    valueType th     = theta0-alpha ;
    valueType theta1 = alpha-th ;
    if ( abs(th) > m_pi/3 ) {
      valueType w   = cos(th/2);
      valueType s   = L/(4*w*w);
      valueType xm0 = x0+s*cos(theta0) ;
      valueType ym0 = y0+s*sin(theta0) ;
      valueType xm1 = x1-s*cos(theta1) ;
      valueType ym1 = y1-s*sin(theta1) ;
      knots[0] = knots[1] = knots[2] = 0 ;
      knots[3] = 1 ;
      knots[4] = knots[5] = knots[6] = 2 ;
      Poly[0][0] = x0 ; Poly[1][0] = w*xm0 ; Poly[2][0] = w*xm1 ; Poly[3][0] = x1 ;
      Poly[0][1] = y0 ; Poly[1][1] = w*ym0 ; Poly[2][1] = w*ym1 ; Poly[3][1] = y1 ;
      Poly[0][2] = 1  ; Poly[1][2] = w     ; Poly[2][2] = w ;     Poly[3][2] = 1  ;
      return 4 ;
    } else {
      valueType w  = cos(th);
      valueType s  = L/(2*w);
      valueType xm = x0+s*cos(theta0) ;
      valueType ym = y0+s*sin(theta0) ;
      knots[0] = knots[1] = knots[2] = 0 ;
      knots[3] = knots[4] = knots[5] = 1 ;
      Poly[0][0] = x0 ; Poly[1][0] = w*xm ; Poly[2][0] = x1 ;
      Poly[0][1] = y0 ; Poly[1][1] = w*ym ; Poly[2][1] = y1 ;
      Poly[0][2] = 1  ; Poly[1][2] = w    ; Poly[2][2] = 1  ;
      return 3 ;
    }
  }

  /*\
   |   ____  _
   |  | __ )(_) __ _ _ __ ___
   |  |  _ \| |/ _` | '__/ __|
   |  | |_) | | (_| | | | (__
   |  |____/|_|\__,_|_|  \___|
  \*/

  bool
  Biarc( valueType   x0,
         valueType   y0,
         valueType   theta0,
         valueType   x1,
         valueType   y1,
         valueType   theta1,
         valueType & xs,
         valueType & ys,
         valueType & thstar,
         valueType & L0,
         valueType & L1 ) {

    valueType b[2] = { x1-x0, y1-y0 } ;
    valueType alpha = atan2(b[1],b[0]);
    // put in range
    rangeSymm(theta0);
    rangeSymm(theta1);
    valueType c0      = cos(theta0);
    valueType s0      = sin(theta0);
    valueType c1      = cos(theta1);
    valueType s1      = sin(theta1);
    valueType thave   = (theta0+theta1)/2;
              thstar  = 2*alpha-thave;
    valueType thstar0 = thstar-theta0;
    valueType thstar1 = thstar-theta1;
    valueType Sinc0   = Sinc(thstar0);
    valueType Cosc0   = Cosc(thstar0);
    valueType Sinc1   = Sinc(thstar1);
    valueType Cosc1   = Cosc(thstar1);

    Solve2x2 solver;
    valueType A[2][2] = {
      { c0*Sinc0-s0*Cosc0, c1*Sinc1-s1*Cosc1 },
      { s0*Sinc0+c0*Cosc0, s1*Sinc1+c1*Cosc1 }
    };
    solver.factorize(A);
    valueType st[2] ;
    solver.solve( b, st );
    valueType epsi = 100*hypot(b[0],b[1])*std::numeric_limits<valueType>::epsilon();
    bool ok = FP_INFINITE != std::fpclassify(st[0]) &&
              FP_NAN      != std::fpclassify(st[0]) &&
              FP_INFINITE != std::fpclassify(st[1]) &&
              FP_NAN      != std::fpclassify(st[1]) &&
              st[0] > -epsi && st[0] > -epsi ;

    if ( ok ) {
      L0 = st[0] ;
      L1 = st[1] ;
      xs = x0 + L0*A[0][0];
      ys = y0 + L0*A[1][0];
    }

    return ok ;
  }

}

///
/// eof: Circles.cc
///
