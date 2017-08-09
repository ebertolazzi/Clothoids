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
  ArcToNURBS( valueType theta0,
              valueType x0,
              valueType y0,
              valueType c0,
              valueType s0,
              valueType L,
              valueType kappa,
              valueType knots[12],
              valueType Poly[9][3] ) {

    valueType dtheta = L*kappa ;
    indexType ns     = std::floor(3*std::abs(dtheta)/m_pi) ;
    if      ( ns < 1 ) ns = 1 ;
    else if ( ns > 4 ) ns = 4 ;
    valueType th = dtheta/(2*ns) ;
    valueType D  = std::abs(tan(th)/kappa);
    valueType w  = cos(th) ;
    knots[0] = knots[1] = knots[2] = 0 ;
    Poly[0][0] = x0 ;
    Poly[0][1] = y0 ;
    Poly[0][2] = 1  ;
    indexType kk = 0 ;
    for ( indexType k = 0 ; k < ns ; ++k ) {
      valueType tth = theta0+kk*th ;
      valueType x1 = Poly[kk][0]+D*cos(tth) ;
      valueType y1 = Poly[kk][1]+D*sin(tth) ;
      valueType kL = ((k+1)*dtheta)/ns ;
      valueType s  = Sinc(kL) ;
      valueType c  = Cosc(kL) ;
      valueType Lk = ((k+1)*L)/ns ;
      valueType x2 = x0+Lk*(c0*s-s0*c) ;
      valueType y2 = y0+Lk*(s0*s+c0*c) ;
      ++kk;
      Poly[kk][0] = w*x1 ;
      Poly[kk][1] = w*y1 ;
      Poly[kk][2] = w  ;
      ++kk;
      Poly[kk][0] = x2 ;
      Poly[kk][1] = y2 ;
      Poly[kk][2] = 1  ;
      knots[kk+1] = k+1 ;
      knots[kk+2] = k+1 ;
    }
    knots[kk+3] = ns ;
    return 1+2*ns;
  }

  /*\
   |   ____  _
   |  | __ )(_) __ _ _ __ ___
   |  |  _ \| |/ _` | '__/ __|
   |  | |_) | | (_| | | | (__
   |  |____/|_|\__,_|_|  \___|
  \*/

  bool
  Biarc( BiarcData & BiData, bool compute_thstar ) {

    valueType dx = BiData.x1-BiData.x0 ;
    valueType dy = BiData.y1-BiData.y0 ;

    BiData.alpha = atan2(dy,dx);
    valueType d = hypot(dy,dx);

    // put in range
    valueType theta0 = BiData.theta0-BiData.alpha;
    valueType theta1 = BiData.theta1-BiData.alpha;

    rangeSymm(theta0);
    rangeSymm(theta1);

    valueType thstar = compute_thstar ? -(theta0+theta1)/2 : BiData.thetas + BiData.alpha ;

    valueType c0 = cos(theta0);
    valueType s0 = sin(theta0);
    valueType c1 = cos(theta1);
    valueType s1 = sin(theta1);

    valueType thstar0 = thstar-theta0;
    valueType thstar1 = thstar-theta1;

    valueType Sinc0 = Sinc(thstar0);
    valueType Cosc0 = Cosc(thstar0);

    valueType Sinc1 = Sinc(thstar1);
    valueType Cosc1 = Cosc(thstar1);

    valueType A[2][2] = {
      { c0*Sinc0-s0*Cosc0, c1*Sinc1-s1*Cosc1 },
      { s0*Sinc0+c0*Cosc0, s1*Sinc1+c1*Cosc1 }
    };

    valueType b[2] = { 1, 0 };

    Solve2x2 solver;
    solver.factorize(A);
    valueType st[2] ;
    bool ok = solver.solve( b, st );
    if ( ok ) {
      valueType epsi = 100*std::numeric_limits<valueType>::epsilon();
      ok = st[0] > epsi && st[1] > epsi ; // NO ZERO LENGHT SOLUTION
    }
    if ( ok ) {

      BiData.L0     = d*st[0];
      BiData.L1     = d*st[1];
      BiData.kappa0 = thstar0/BiData.L0;
      BiData.kappa1 = -thstar1/BiData.L1;

      valueType ca  = cos(BiData.alpha);
      valueType sa  = sin(BiData.alpha);
      BiData.xs = BiData.x0 + BiData.L0*(A[0][0]*ca-A[1][0]*sa);
      BiData.ys = BiData.y0 + BiData.L0*(A[0][0]*sa+A[1][0]*ca);

      BiData.thetas = thstar+BiData.alpha;
      BiData.cs     = cos(BiData.thetas);
      BiData.ss     = sin(BiData.thetas);

      BiData.c0 = cos(BiData.theta0);
      BiData.s0 = sin(BiData.theta0);
      BiData.c1 = cos(BiData.theta1);
      BiData.s1 = sin(BiData.theta1);
    }

    return ok ;
  }

}

///
/// eof: Circles.cc
///
