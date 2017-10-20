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

///
/// file: Biarc.cc
///

#include "Biarc.hh"

#include <limits>

namespace Biarc {

  using namespace std ;

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

  /*\
   |   ____  _
   |  | __ )(_) __ _ _ __ ___
   |  |  _ \| |/ _` | '__/ __|
   |  | |_) | | (_| | | | (__
   |  |____/|_|\__,_|_|  \___|
  \*/

  bool
  Biarc::build( valueType x0,
                valueType y0,
                valueType theta0,
                valueType x1,
                valueType y1,
                valueType theta1 ) {

    valueType dx = x1-x0 ;
    valueType dy = y1-y0 ;

    alpha = atan2(dy,dx);
    valueType d = hypot(dy,dx);

    // put in range
    valueType th0 = theta0-alpha;
    valueType th1 = theta1-alpha;

    rangeSymm(th0);
    rangeSymm(th1);

    //valueType thstar = compute_thstar ? -(th0+th1)/2 : BiData.thetas + alpha ;

    valueType thstar = -(th0+th1)/2 ;

    valueType c0 = cos(th0);
    valueType s0 = sin(th0);
    valueType c1 = cos(th1);
    valueType s1 = sin(th1);

    valueType thstar0 = thstar-th0;
    valueType thstar1 = thstar-th1;

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

      valueType L0     = d*st[0];
      valueType L1     = d*st[1];
      valueType kappa0 = thstar0/L0;
      valueType kappa1 = -thstar1/L1;

      C0.build( x0, y0, theta0, kappa0, L0 );

      valueType ca = cos(alpha);
      valueType sa = sin(alpha);

      xs     = x0 + L0*(A[0][0]*ca-A[1][0]*sa);
      ys     = y0 + L0*(A[0][0]*sa+A[1][0]*ca);
      thetas = thstar+alpha;
      cs     = cos(thetas);
      ss     = sin(thetas);

      C1.build( xs, ys, thetas, kappa1, L1 );
    }

    return ok ;
  }

  valueType
  Biarc::X( valueType s ) const {
    if ( s < C0.getL() ) return C0.X(s);
    else                 return C1.X(s-C0.getL());
  }

  valueType
  Biarc::Y( valueType s ) const {
    if ( s < C0.getL() ) return C0.Y(s);
    else                 return C1.Y(s-C0.getL());
  }

  valueType
  Biarc::theta( valueType s ) const {
    if ( s < C0.getL() ) return C0.theta(s);
    else                 return C1.theta(s-C0.getL());
  }

  std::ostream &
  operator << ( std::ostream & stream, Biarc const & bi ) {
    stream <<   "Biarc"
           << "\nC0 = " << bi.C0
           << "\nC1 = " << bi.C1
           << "\nxs     = " << bi.xs
           << "\nys     = " << bi.ys
           << "\nthetas = " << bi.thetas
           << "\n" ;
    return stream ;
  }

}

///
/// eof: Biarc.cc
///
