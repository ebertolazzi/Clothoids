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

  using namespace std;

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
  real_type
  powersub( real_type a, real_type b)
  { return (a+b)*(a-b); }

  /*
  //  http://www.lucidarme.me/?p=490
  */

  static // unused for the moment
  void
  CircleTangentPoints( real_type PA[2],
                       real_type rA,
                       real_type PB[2],
                       real_type rB,
                       bool &    external_tangents,
                       real_type PTE0[2][2],
                       real_type PTE1[2][2],
                       bool &    internal_tangents,
                       real_type PTI0[2][2],
                       real_type PTI1[2][2] ) {

    // Compute distance between circle centers
    real_type D = hypot( PB[0]-PA[0], PB[1]-PA[1] );

    // First case : process external tangents
    real_type disc1 = powersub(D,rB-rA);
    external_tangents = disc1 >= 0;

    if ( external_tangents ) {
      // Compute the lenght of the tangents
      real_type L = sqrt(disc1);

      // Compute the parameters
      real_type R1     = hypot(L,rB);
      real_type Sigma1 = 0.25 * sqrt( powersub( D+rA, R1 )* powersub( R1, D-rA ) );
      real_type R2     = hypot(L,rA);
      real_type Sigma2 = 0.25 * sqrt( powersub( D+rB, R2 )* powersub( R2, D-rB ) );

      real_type xM  = (PA[0]+PB[0])/2;
      real_type yM  = (PA[1]+PB[1])/2;
      real_type Tx  = 2*(PA[0]-PB[0])/(D*D);
      real_type Ty  = 2*(PA[1]-PB[1])/(D*D);
      real_type bf0 = powersub(rA,R1)/4;
      real_type bf1 = powersub(rB,R2)/4;

      real_type TxS1 = Tx*Sigma1;
      real_type TxS2 = Tx*Sigma2;
      real_type TyS1 = Ty*Sigma1;
      real_type TyS2 = Ty*Sigma2;

      real_type xM0 = xM - Tx*bf0;
      real_type yM0 = yM - Ty*bf0;
      real_type xM1 = xM + Tx*bf1;
      real_type yM1 = yM + Ty*bf1;

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
    real_type disc2 = powersub(D,rB+rA);
    internal_tangents = disc2 >= 0;

    if ( internal_tangents ) {
      // Compute the lenght of the tangents
      real_type L = sqrt(disc2);

      // Compute the parameters
      real_type R1     = hypot(L,rB);
      real_type Sigma1 = 0.25 * sqrt ( powersub(D+rA,R1)*powersub(R1,D-rA) );
      real_type R2     = hypot(L,rA);
      real_type Sigma2 = 0.25 * sqrt ( powersub(D+rB,R2)*powersub(R2,D-rB) );

      real_type xM  = (PA[0]+PB[0])/2;
      real_type yM  = (PA[1]+PB[1])/2;
      real_type Tx  = 2*(PA[0]-PB[0])/(D*D);
      real_type Ty  = 2*(PA[1]-PB[1])/(D*D);
      real_type bf0 = powersub(rA,R1)/4;
      real_type bf1 = powersub(rB,R2)/4;

      real_type TxS1 = Tx*Sigma1;
      real_type TxS2 = Tx*Sigma2;
      real_type TyS1 = Ty*Sigma1;
      real_type TyS2 = Ty*Sigma2;

      real_type xM0 = xM - Tx*bf0;
      real_type yM0 = yM - Ty*bf0;
      real_type xM1 = xM + Tx*bf1;
      real_type yM1 = yM + Ty*bf1;

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
  CircleLineTransition( real_type C[2],
                        real_type r,
                        real_type P[2],
                        real_type theta,
                        real_type C0[2],
                        real_type C1[2] ) {
    real_type Nx =  sin(theta);
    real_type Ny = -cos(theta);
    real_type Dx = C[0] - P[0];
    real_type Dy = C[1] - P[1];
    real_type delta = (Dx*Dx+Dy*Dy-r*r)/2;
    if ( delta <= 0 ) return false;
    real_type s0 = delta/(Dx*Nx+Dy*Ny-r);
    C0[0] = P[0]+s0*Nx;
    C0[1] = P[1]+s0*Ny;
    real_type s1 = delta/(Dx*Nx+Dy*Ny+r);
    C1[0] = P[0]+s1*Nx;
    C1[1] = P[1]+s1*Ny;
    return true;
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
  Biarc::build( real_type x0,
                real_type y0,
                real_type theta0,
                real_type x1,
                real_type y1,
                real_type theta1 ) {

    real_type dx = x1-x0;
    real_type dy = y1-y0;
    real_type d  = hypot(dy,dx);

    omega = atan2(dy,dx);

    // put in range
    real_type th0 = theta0 - omega;
    real_type th1 = theta1 - omega;

    rangeSymm(th0);
    rangeSymm(th1);

    real_type thstar = - (th0+th1)/2;

    real_type dth  = th1 - th0;
    real_type dth0 = thstar - th0;
    real_type dth1 = thstar - th1;

    real_type t  = d * (Sinc(dth/4) / Sinc(dth/2) );
    real_type l0 = t/(2*Sinc( dth0/2 ));
    real_type l1 = t/(2*Sinc( dth1/2 ));

    real_type epsi = 100*std::numeric_limits<real_type>::epsilon();

    if ( l0 > epsi && l1 > epsi ) {

      real_type k0 = dth0/l0;
      real_type k1 = -dth1/l1;

      C0.build( x0, y0, theta0, k0, l0 );

      real_type an = omega+(thstar+th0)/2;
      real_type xs = x0 + (t/2)*cos(an);
      real_type ys = y0 + (t/2)*sin(an);

      C1.build( xs, ys, omega+thstar, k1, l1 );
      return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Biarc::build_3P( real_type x0,
                   real_type y0,
                   real_type x1,
                   real_type y1,
                   real_type x2,
                   real_type y2 ) {

    real_type dxa   = x1-x0;
    real_type dya   = y1-y0;
    real_type dxb   = x2-x1;
    real_type dyb   = y2-y1;
    real_type La    = hypot(dya,dxa);
    real_type Lb    = hypot(dyb,dxb);
    real_type arg   = (dxa*dxb + dya * dyb)/(La*Lb);
    if      ( arg >  1 ) arg = 1;
    else if ( arg < -1 ) arg = -1;
    real_type om = acos(arg);

    real_type at = (La/(La+Lb))*om;
    real_type bt = (Lb/(La+Lb))*om;
    // find solution using Halley
    real_type Delta = 0;
    bool found = false;
    for ( int_type iter = 0; iter < 10 && !found; ++iter ) {
      real_type ga[3], gb[3];
      gfun( at+Delta, ga );
      gfun( bt-Delta, gb );
      real_type f   = ga[0]/La - gb[0]/Lb;
      real_type df  = ga[1]/La + gb[1]/Lb;
      real_type ddf = ga[2]/La - gb[2]/Lb;
      real_type h   = (df*f)/(df*df-0.5*f*ddf);
      Delta -= h;
      found = abs(h) < 1e-10 && abs(f) < 1e-10;
    }

    if ( found ) {
      at += Delta; bt -= Delta;
      real_type tha = atan2(dya,dxa);
      real_type thb = atan2(dyb,dxb);
      if ( dxa*dyb < dya*dxb ) {
        tha += at;
        thb += bt;
      } else {
        tha -= at;
        thb -= bt;
      }
      C0.build_G1( x0, y0, tha, x1, y1 );
      C1.build_G1( x1, y1, thb, x2, y2 );
    }

    return found;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::theta( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.theta(s);
    else          return C1.theta(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::theta_D( real_type s ) const
  { return kappa(s); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::theta_DD( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::theta_DDD( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::kappa( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.kappa();
    else          return C1.kappa();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::kappa_D( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::kappa_DD( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::kappa_DDD( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X(s);
    else          return C1.X(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_D( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_D(s);
    else          return C1.X_D(s-L0);
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_DD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_DD(s);
    else          return C1.X_DD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_DDD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_DDD(s);
    else          return C1.X_DDD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y(s);
    else          return C1.Y(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_D( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_D(s);
    else          return C1.Y_D(s-L0);
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_DD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_DD(s);
    else          return C1.Y_DD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_DDD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_DDD(s);
    else          return C1.Y_DDD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::XY( real_type s, real_type & x, real_type & y ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.XY(s,x,y);
    else          return C1.XY(s-L0,x,y);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::XY( real_type s, real_type t, real_type & x, real_type & y ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.XY(s,t,x,y);
    else          return C1.XY(s-L0,t,x,y);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::TG( real_type s, real_type & tx, real_type & ty ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.TG(s,tx,ty);
    else          return C1.TG(s-L0,tx,ty);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::NOR( real_type s, real_type & nx, real_type & ny ) const{
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.NOR(s,nx,ny);
    else          return C1.NOR(s-L0,nx,ny);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::NOR_D( real_type s, real_type & nx_D, real_type & ny_D ) const{
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.NOR_D(s,nx_D,ny_D);
    else          return C1.NOR_D(s-L0,nx_D,ny_D);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::NOR_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const{
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.NOR_DD(s,nx_DD,ny_DD);
    else          return C1.NOR_DD(s-L0,nx_DD,ny_DD);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::NOR_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const{
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.NOR_DDD(s,nx_DDD,ny_DDD);
    else          return C1.NOR_DDD(s-L0,nx_DDD,ny_DDD);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::findST( real_type   x,
                 real_type   y,
                 real_type & s,
                 real_type & t ) const {
    real_type S, T;
    C0.findST( x, y, s, t );
    C1.findST( x, y, S, T );
    real_type L0 = C0.length();
    bool ok0 = s >= 0 && s <= L0;
    bool ok1 = S >= 0 && S <= C1.length();
    if ( ok0 ) {
      if ( ok1 && std::abs(t) > std::abs(T) ) {
        s = S + L0;
        t = T;
      }
    } else if ( ok1 ) {
      s = S + L0;
      t = T;
    } else {
      s = t = 0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval( real_type   s,
               real_type & th,
               real_type & k,
               real_type & x,
               real_type & y ) const {
    if ( s < C0.length() ) {
      th = C0.theta(s);
      k  = C0.kappa();
      C0.eval(s,x,y);
    } else {
      s -= C0.length();
      th = C1.theta(s);
      k  = C1.kappa();
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval( real_type   s,
               real_type & x,
               real_type & y ) const {
    if ( s < C0.length() ) {
      C0.eval(s,x,y);
    } else {
      s -= C0.length();
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_D( real_type   s,
                 real_type & x_D,
                 real_type & y_D ) const {
    if ( s < C0.length() ) {
      C0.eval_D(s,x_D,y_D);
    } else {
      s -= C0.length();
      C1.eval_D(s,x_D,y_D);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DD( real_type   s,
                  real_type & x_DD,
                  real_type & y_DD ) const {
    if ( s < C0.length() ) {
      C0.eval_DD(s,x_DD,y_DD);
    } else {
      s -= C0.length();
      C1.eval_DD(s,x_DD,y_DD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DDD( real_type   s,
                   real_type & x_DDD,
                   real_type & y_DDD ) const {
    if ( s < C0.length() ) {
      C0.eval_DDD(s,x_DDD,y_DDD);
    } else {
      s -= C0.length();
      C1.eval_DDD(s,x_DDD,y_DDD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval( real_type   s,
               real_type   t,
               real_type & x,
               real_type & y ) const {
    if ( s < C0.length() ) {
      C0.eval(s,t,x,y);
    } else {
      s -= C0.length();
      C1.eval(s,t,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_D( real_type   s,
                 real_type   t,
                 real_type & x_D,
                 real_type & y_D ) const {
    if ( s < C0.length() ) {
      C0.eval_D(s,t,x_D,y_D);
    } else {
      s -= C0.length();
      C1.eval_D(s,t,x_D,y_D);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DD( real_type   s,
                  real_type   t,
                  real_type & x_DD,
                  real_type & y_DD ) const {
    if ( s < C0.length() ) {
      C0.eval_DD(s,t,x_DD,y_DD);
    } else {
      s -= C0.length();
      C1.eval_DD(s,t,x_DD,y_DD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DDD( real_type   s,
                   real_type   t,
                   real_type & x_DDD,
                   real_type & y_DDD ) const {
    if ( s < C0.length() ) {
      C0.eval_DDD(s,t,x_DDD,y_DDD);
    } else {
      s -= C0.length();
      C1.eval_DDD(s,t,x_DDD,y_DDD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::closestPoint( real_type   x,
                       real_type   y,
                       real_type & X,
                       real_type & Y,
                       real_type & S ) const {
    real_type dst0 = C0.closestPoint( x, y, X, Y, S );
    real_type X1, Y1, S1;
    real_type dst1 = C1.closestPoint( x, y, X1, Y1, S1 );
    if ( dst0 <= dst1 ) return dst0;
    X = X1; Y= Y1; S = S1 + C0.length();
    return dst1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::changeOrigin( real_type newx0, real_type newy0 ) {
    C1.translate( newx0-C0.xBegin(), newy0-C0.yBegin() );
    C0.changeOrigin( newx0, newy0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::translate( real_type tx, real_type ty ) {
    C0.translate( tx, ty );
    C1.translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::rotate( real_type angle, real_type cx, real_type cy ) {
    C0.rotate( angle, cx, cy );
    C1.rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::reverse() {
    CircleArc tmp = C0;
    C0 = C1;
    C1 = tmp;
    C0.reverse();
    C1.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::scale( real_type scl ) {
    real_type newx0 = C0.xBegin() + scl*(C1.xBegin()-C0.xBegin());
    real_type newy0 = C0.yBegin() + scl*(C1.yBegin()-C0.yBegin());
    C1.changeOrigin( newx0, newy0 );
    C1.scale( scl );
    C0.scale( scl );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, Biarc const & bi ) {
    stream << "Biarc\n"
           << "C0\n" << bi.C0
           << "C1\n" << bi.C1
           << "\n";
    return stream;
  }

}

///
/// eof: Biarc.cc
///
