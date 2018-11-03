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

  using std::numeric_limits;
  using std::abs;

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
  CircleTangentPoints(
    real_type PA[2],
    real_type rA,
    real_type PB[2],
    real_type rB,
    bool &    external_tangents,
    real_type PTE0[2][2],
    real_type PTE1[2][2],
    bool &    internal_tangents,
    real_type PTI0[2][2],
    real_type PTI1[2][2]
  ) {

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
  CircleLineTransition(
    real_type C[2],
    real_type r,
    real_type P[2],
    real_type theta,
    real_type C0[2],
    real_type C1[2]
  ) {
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

  Biarc::Biarc( BaseCurve const & C )
  : BaseCurve(G2LIB_BIARC)
  {
    switch ( C.type() ) {
    case G2LIB_BIARC:
      *this = *static_cast<Biarc const *>(&C);
      break;
    case G2LIB_LINE:
    case G2LIB_CIRCLE:
    case G2LIB_CLOTHOID:
    case G2LIB_CLOTHOID_LIST:
    case G2LIB_POLYLINE:
      G2LIB_ASSERT( false,
                    "Biarc constructor cannot convert from: " <<
                    CurveType_name[C.type()] );
    }
  }

  bool
  Biarc::build(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1
  ) {

    real_type dx = x1-x0;
    real_type dy = y1-y0;
    real_type d  = hypot(dy,dx);

    real_type omega = atan2(dy,dx);

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

    real_type epsi = 100*numeric_limits<real_type>::epsilon();

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
  Biarc::build_3P(
    real_type x0,
    real_type y0,
    real_type x1,
    real_type y1,
    real_type x2,
    real_type y2
  ) {

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

  void
  Biarc::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    C0.bbox( xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    C1.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::bbox(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    C0.bbox( offs, xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    C1.bbox( offs, xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  /*\
   |  _                        __
   | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
   | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
   | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
   |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
  \*/

  void
  Biarc::reverse() {
    CircleArc TMP(C0);
    C0 = C1;  C0.reverse();
    C1 = TMP; C1.reverse();
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

  void
  Biarc::changeOrigin( real_type newx0, real_type newy0 ) {
    C0.changeOrigin(newx0,newy0);
    C1.changeOrigin(C0.xEnd(),C0.yEnd());
  }

  void
  Biarc::trim( real_type s_begin, real_type s_end ) {
    G2LIB_ASSERT( s_end > s_begin,
                  "Biarc::trim(begin=" << s_begin <<
                  ", s_end=" << s_end << ") s_end must be > s_begin" );
    real_type L0 = C0.length();
    if ( s_end <= L0 ) {
      C0.trim( s_begin, s_end );
      C1 = C0;
      real_type ss = C0.length();
      C0.trim( 0, ss/2 );
      C1.trim( ss/2, ss );
    } else if ( s_begin >= L0 ) {
      C1.trim( s_begin-L0, s_end-L0 );
      C0 = C1;
      real_type ss = C0.length();
      C0.trim( 0, ss/2 );
      C1.trim( ss/2, ss );
    } else {
      C0.trim( s_begin, L0 );
      C1.trim( 0, s_end-L0 );
    }
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
  Biarc::theta_D( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.k;
    else          return C1.k;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/*
  real_type
  Biarc::kappa( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.k;
    else          return C1.k;
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
  Biarc::tx( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tx(s);
    else          return C1.tx(s-L0);
  }
*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::tx( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tx(s);
    else          return C1.tx(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::tx_D( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tx_D(s);
    else          return C1.tx_D(s-L0);
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::tx_DD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tx_DD(s);
    else          return C1.tx_DD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::tx_DDD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tx_DDD(s);
    else          return C1.tx_DDD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::ty( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.ty(s);
    else          return C1.ty(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::ty_D( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.ty_D(s);
    else          return C1.ty_D(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::ty_DD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.ty_DD(s);
    else          return C1.ty_DD(s-L0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::ty_DDD( real_type s ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.ty_DDD(s);
    else          return C1.ty_DDD(s-L0);
  }

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

  real_type
  Biarc::X( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X(s,offs);
    else          return C1.X(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_D( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_D(s,offs);
    else          return C1.X_D(s-L0,offs);
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_DD( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_DD(s,offs);
    else          return C1.X_DD(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::X_DDD( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.X_DDD(s,offs);
    else          return C1.X_DDD(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y(s,offs);
    else          return C1.Y(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_D( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_D(s,offs);
    else          return C1.Y_D(s-L0,offs);
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_DD( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_DD(s,offs);
    else          return C1.Y_DD(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Biarc::Y_DDD( real_type s, real_type offs ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.Y_DDD(s,offs);
    else          return C1.Y_DDD(s-L0,offs);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::tg( real_type s, real_type & tx, real_type & ty ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tg(s,tx,ty);
    else          return C1.tg(s-L0,tx,ty);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::tg_D( real_type s, real_type & tx_D, real_type & ty_D ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tg_D(s,tx_D,ty_D);
    else          return C1.tg_D(s-L0,tx_D,ty_D);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::tg_DD( real_type s, real_type & tx_DD, real_type & ty_DD ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tg_DD(s,tx_DD,ty_DD);
    else          return C1.tg_DD(s-L0,tx_DD,ty_DD);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::tg_DDD( real_type s, real_type & tx_DDD, real_type & ty_DDD ) const {
    real_type L0 = C0.length();
    if ( s < L0 ) return C0.tg_DDD(s,tx_DDD,ty_DDD);
    else          return C1.tg_DDD(s-L0,tx_DDD,ty_DDD);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::evaluate(
    real_type   s,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    if ( s < C0.length() ) {
      th = C0.theta(s);
      k  = C0.curvature();
      C0.eval(s,x,y);
    } else {
      s -= C0.length();
      th = C1.theta(s);
      k  = C1.curvature();
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    if ( s < C0.length() ) {
      C0.eval(s,x,y);
    } else {
      s -= C0.length();
      C1.eval(s,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( s < C0.length() ) {
      C0.eval_D(s,x_D,y_D);
    } else {
      s -= C0.length();
      C1.eval_D(s,x_D,y_D);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( s < C0.length() ) {
      C0.eval_DD(s,x_DD,y_DD);
    } else {
      s -= C0.length();
      C1.eval_DD(s,x_DD,y_DD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( s < C0.length() ) {
      C0.eval_DDD(s,x_DDD,y_DDD);
    } else {
      s -= C0.length();
      C1.eval_DDD(s,x_DDD,y_DDD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval(
    real_type   s,
    real_type   t,
    real_type & x,
    real_type & y
  ) const {
    if ( s < C0.length() ) {
      C0.eval(s,t,x,y);
    } else {
      s -= C0.length();
      C1.eval(s,t,x,y);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_D(
    real_type   s,
    real_type   t,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( s < C0.length() ) {
      C0.eval_D(s,t,x_D,y_D);
    } else {
      s -= C0.length();
      C1.eval_D(s,t,x_D,y_D);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DD(
    real_type   s,
    real_type   t,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( s < C0.length() ) {
      C0.eval_DD(s,t,x_DD,y_DD);
    } else {
      s -= C0.length();
      C1.eval_DD(s,t,x_DD,y_DD);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::eval_DDD(
    real_type   s,
    real_type   t,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( s < C0.length() ) {
      C0.eval_DDD(s,t,x_DDD,y_DDD);
    } else {
      s -= C0.length();
      C1.eval_DDD(s,t,x_DDD,y_DDD);
    }
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::intersect(
    Biarc const   & B,
    IntersectList & ilist,
    bool            swap_s_vals
  ) const {
    IntersectList ilist00, ilist01, ilist10, ilist11;
    C0.intersect( B.C0, ilist00, false );
    C0.intersect( B.C1, ilist01, false );
    C1.intersect( B.C0, ilist10, false );
    C1.intersect( B.C1, ilist11, false );
    real_type L  = C0.length();
    real_type LB = B.C0.length();
    IntersectList::iterator it;
    ilist.reserve( ilist.size() +
                   ilist00.size() +
                   ilist01.size() +
                   ilist10.size() +
                   ilist11.size() );
    for ( it = ilist01.begin(); it != ilist01.end(); ++it ) it->second += LB;
    for ( it = ilist10.begin(); it != ilist10.end(); ++it ) it->first  += L;
    for ( it = ilist11.begin(); it != ilist11.end(); ++it )
      { it->first += L; it->second += LB; }

    if ( swap_s_vals ) {
      for ( it = ilist00.begin(); it != ilist00.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist01.begin(); it != ilist01.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist10.begin(); it != ilist10.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist11.begin(); it != ilist11.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
    } else {
      for ( it = ilist00.begin(); it != ilist00.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist01.begin(); it != ilist01.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist10.begin(); it != ilist10.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist11.begin(); it != ilist11.end(); ++it ) ilist.push_back( *it );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Biarc::intersect(
    real_type       offs,
    Biarc const   & B,
    real_type       offs_B,
    IntersectList & ilist,
    bool            swap_s_vals
  ) const {
    IntersectList ilist00, ilist01, ilist10, ilist11;
    C0.intersect( offs, B.C0, offs_B, ilist00, false );
    C0.intersect( offs, B.C1, offs_B, ilist01, false );
    C1.intersect( offs, B.C0, offs_B, ilist10, false );
    C1.intersect( offs, B.C1, offs_B, ilist11, false );
    real_type L  = C0.length();
    real_type LB = B.C0.length();
    IntersectList::iterator it;
    ilist.reserve( ilist.size() +
                   ilist00.size() +
                   ilist01.size() +
                   ilist10.size() +
                   ilist11.size() );
    for ( it = ilist01.begin(); it != ilist01.end(); ++it ) it->second += LB;
    for ( it = ilist10.begin(); it != ilist10.end(); ++it ) it->first  += L;
    for ( it = ilist11.begin(); it != ilist11.end(); ++it )
      { it->first += L; it->second += LB; }

    if ( swap_s_vals ) {
      for ( it = ilist00.begin(); it != ilist00.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist01.begin(); it != ilist01.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist10.begin(); it != ilist10.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
      for ( it = ilist11.begin(); it != ilist11.end(); ++it )
        ilist.push_back( Ipair(it->second,it->first) );
    } else {
      for ( it = ilist00.begin(); it != ilist00.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist01.begin(); it != ilist01.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist10.begin(); it != ilist10.end(); ++it ) ilist.push_back( *it );
      for ( it = ilist11.begin(); it != ilist11.end(); ++it ) ilist.push_back( *it );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  Biarc::closestPoint(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    int_type res  = C0.closestPoint( qx, qy, x,  y,  s,  t,  dst  );
    int_type res1 = C1.closestPoint( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1;
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  Biarc::closestPoint(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    int_type res  = C0.closestPoint( qx, qy, offs, x,  y,  s,  t,  dst  );
    int_type res1 = C1.closestPoint( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1;
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, Biarc const & bi ) {
    stream << "C0\n" << bi.C0
           << "C1\n" << bi.C1
           << "\n";
    return stream;
  }

}

///
/// eof: Biarc.cc
///
