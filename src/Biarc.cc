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

namespace G2lib {

  using namespace std ;

  /*\
   |    ____ _          _
   |   / ___(_)_ __ ___| | ___  ___
   |  | |   | | '__/ __| |/ _ \/ __|
   |  | |___| | | | (__| |  __/\__ \
   |   \____|_|_|  \___|_|\___||___/
  \*/

  #if 0

  static
  inline
  valueType
  powersub( valueType a, valueType b)
  { return (a+b)*(a-b); }

  /*
  //  http://www.lucidarme.me/?p=490
  */

  static // unused for the moment
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
      PTI1[0][0] = xM0 - TyS1;
      PTI1[0][1] = yM0 + TxS1;
      PTI1[1][0] = xM1 + TyS2;
      PTI1[1][1] = yM1 - TxS2;
    }
  }

  static // unused for the moment
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

  #endif

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

      valueType xs     = x0 + L0*(A[0][0]*ca-A[1][0]*sa);
      valueType ys     = y0 + L0*(A[0][0]*sa+A[1][0]*ca);
      valueType thetas = thstar+alpha;
      //valueType cs     = cos(thetas);
      //valueType ss     = sin(thetas);

      C1.build( xs, ys, thetas, kappa1, L1 );
    }

    return ok ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Biarc::build_3P( valueType x0,
                   valueType y0,
                   valueType x1,
                   valueType y1,
                   valueType x2,
                   valueType y2 ) {

    valueType dxa   = x1-x0 ;
    valueType dya   = y1-y0 ;
    valueType dxb   = x2-x1 ;
    valueType dyb   = y2-y1 ;
    valueType La    = hypot(dya,dxa) ;
    valueType Lb    = hypot(dyb,dxb) ;
    valueType arg   = (dxa*dxb + dya * dyb)/(La*Lb) ;
    if      ( arg >  1 ) arg = 1 ;
    else if ( arg < -1 ) arg = -1 ;
    valueType omega = acos(arg) ;

    valueType at = (La/(La+Lb))*omega;
    valueType bt = (Lb/(La+Lb))*omega;
    // find solution using Halley
    valueType Delta = 0 ;
    bool found = false ;
    for ( indexType iter = 0 ; iter < 10 && !found ; ++iter ) {
      valueType ga[3], gb[3] ;
      gfun( at+Delta, ga );
      gfun( bt-Delta, gb );
      valueType f   = ga[0]/La - gb[0]/Lb ;
      valueType df  = ga[1]/La + gb[1]/Lb ;
      valueType ddf = ga[2]/La - gb[2]/Lb ;
      valueType h   = (df*f)/(df*df-0.5*f*ddf) ;
      Delta -= h ;
      found = abs(h) < 1e-10 && abs(f) < 1e-10 ;
    }

    if ( found ) {
      at += Delta ; bt -= Delta ;
      valueType tha = atan2(dya,dxa) ;
      valueType thb = atan2(dyb,dxb) ;
      if ( dxa*dyb < dya*dxb ) {
        tha += at ;
        thb += bt ;
      } else {
        tha -= at ;
        thb -= bt ;
      }
      C0.build_G1( x0, y0, tha, x1, y1 );
      C1.build_G1( x1, y1, thb, x2, y2 );
    }

    return found ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  Biarc::X( valueType s ) const {
    if ( s < C0.length() ) return C0.X(s);
    else                   return C1.X(s-C0.length());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  Biarc::Y( valueType s ) const {
    if ( s < C0.length() ) return C0.Y(s);
    else                   return C1.Y(s-C0.length());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  Biarc::theta( valueType s ) const {
    if ( s < C0.length() ) return C0.theta(s);
    else                   return C1.theta(s-C0.length());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  Biarc::kappa( valueType s ) const {
    if ( s < C0.length() ) return C0.kappa();
    else                   return C1.kappa();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  Biarc::eval( valueType   s,
               valueType & th,
               valueType & k,
               valueType & x,
               valueType & y ) const {
    if ( s < C0.length() ) {
      th = C0.theta(s) ;
      k  = C0.kappa() ;
      C0.eval(s,x,y);
    } else {
      s -= C0.length() ;
      th = C1.theta(s) ;
      k  = C1.kappa() ;
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  Biarc::eval( valueType s, valueType & x, valueType & y ) const {
    if ( s < C0.length() ) {
      C0.eval(s,x,y);
    } else {
      s -= C0.length() ;
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  Biarc::eval_D( valueType s, valueType & x_D, valueType & y_D ) const {
    if ( s < C0.length() ) {
      C0.eval_D(s,x_D,y_D);
    } else {
      s -= C0.length() ;
      C1.eval_D(s,x_D,y_D);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  Biarc::eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const {
    if ( s < C0.length() ) {
      C0.eval_DD(s,x_DD,y_DD);
    } else {
      s -= C0.length() ;
      C1.eval_DD(s,x_DD,y_DD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  void
  Biarc::eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const {
    if ( s < C0.length() ) {
      C0.eval_DDD(s,x_DDD,y_DDD);
    } else {
      s -= C0.length() ;
      C1.eval_DDD(s,x_DDD,y_DDD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  Biarc::closestPoint( valueType   x,
                       valueType   y,
                       valueType & X,
                       valueType & Y,
                       valueType & S ) const {
    valueType dst0 = C0.closestPoint( x, y, X, Y, S );
    valueType X1, Y1, S1 ;
    valueType dst1 = C1.closestPoint( x, y, X1, Y1, S1 );
    if ( dst0 <= dst1 ) return dst0 ;
    X = X1 ; Y= Y1 ; S = S1 + C0.length() ;
    return dst1 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::changeOrigin( valueType newx0, valueType newy0 ) {
    C1.translate( newx0-C0.xBegin(), newy0-C0.yBegin() );
    C0.changeOrigin( newx0, newy0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::translate( valueType tx, valueType ty ) {
    C0.translate( tx, ty );
    C1.translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::rotate( valueType angle, valueType cx, valueType cy ) {
    C0.rotate( angle, cx, cy );
    C1.rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::reverse() {
    CircleArc tmp = C0 ;
    C0 = C1 ;
    C1 = tmp ;
    C0.reverse();
    C1.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::scale( valueType scl ) {
    valueType newx0 = C0.xBegin() + scl*(C1.xBegin()-C0.xBegin()) ;
    valueType newy0 = C0.yBegin() + scl*(C1.yBegin()-C0.yBegin()) ;
    C1.changeOrigin( newx0, newy0 );
    C1.scale( scl );
    C0.scale( scl );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::ostream &
  operator << ( std::ostream & stream, Biarc const & bi ) {
    stream <<   "Biarc"
           << "\nC0 = " << bi.C0
           << "\nC1 = " << bi.C1
           << "\n" ;
    return stream ;
  }

}

///
/// eof: Biarc.cc
///
