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

#include "Clothoid.hh"
#include "Biarc.hh"
//#include "CubicRootsFlocke.hh"

#include <cmath>
#include <cfloat>

namespace Clothoid {

  using namespace std ;

  static const valueType m_pi = 3.14159265358979323846264338328  ; // pi

  /*\
   |    ____ ____     _       _
   |   / ___|___ \ __| | __ _| |_ __ _
   |  | |  _  __) / _` |/ _` | __/ _` |
   |  | |_| |/ __/ (_| | (_| | || (_| |
   |   \____|_____\__,_|\__,_|\__\__,_|
  \*/

  void
  G2data::setup( valueType _x0,
                 valueType _y0,
                 valueType _theta0,
                 valueType _kappa0,
                 valueType _x1,
                 valueType _y1,
                 valueType _theta1,
                 valueType _kappa1 ) {

    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = _theta0;
    kappa0 = _kappa0 ;
    x1     = _x1 ;
    y1     = _y1 ;
    theta1 = _theta1 ;
    kappa1 = _kappa1 ;

    // scale problem
    valueType dx = x1 - x0 ;
    valueType dy = y1 - y0 ;
    phi    = atan2( dy, dx ) ;
    lambda = hypot( dx, dy ) ;

    valueType C = dx/lambda ;
    valueType S = dy/lambda ;
    lambda /= 2 ;

    xbar = -(x0*C+y0*S+lambda) ;
    ybar = x0*S-y0*C ;

    th0 = theta0 - phi ;
    th1 = theta1 - phi ;

    k0 = kappa0*lambda ;
    k1 = kappa1*lambda ;

    DeltaK     = k1 - k0 ;
    DeltaTheta = th1 - th0 ;
  }

  void
  G2data::setTolerance( valueType tol ) {
    G2LIB_ASSERT( tol > 0 && tol <= 0.1,
                  "setTolerance, tolerance = " << tol << " must be in (0,0.1]" ) ;
    tolerance = tol ;
  }

  void
  G2data::setMaxIter( int miter ) {
    G2LIB_ASSERT( miter > 0 && miter <= 1000,
                  "setMaxIter, maxIter = " << miter << " must be in [1,1000]" ) ;
    maxIter = miter ;
  }


  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/

  void
  G2solve2arc::evalA( valueType   alpha,
                      valueType   L,
                      valueType & A,
                      valueType & A_1,
                      valueType & A_2 ) const {
    valueType K  = k0+k1 ;
    valueType aK = alpha*DeltaK ;
    A   = alpha*(L*(aK-K)+2*DeltaTheta) ;
    A_1 = (2*aK-K)*L+2*DeltaTheta;
    A_2 = alpha*(aK-K) ;
  }

  void
  G2solve2arc::evalG( valueType alpha,
                      valueType L,
                      valueType th,
                      valueType k,
                      valueType G[2],
                      valueType G_1[2],
                      valueType G_2[2] ) const {

    valueType A, A_1, A_2, X[3], Y[3] ;
    evalA( alpha, L, A, A_1, A_2 ) ;
    valueType ak = alpha*k ;
    valueType Lk = L*k ;
    GeneralizedFresnelCS( 3, A, ak*L, th, X, Y );

    G[0]   = alpha*X[0] ;
    G_1[0] = X[0]-alpha*(Y[2]*A_1/2+Y[1]*Lk) ;
    G_2[0] =     -alpha*(Y[2]*A_2/2+Y[1]*ak) ;

    G[1]   = alpha*Y[0] ;
    G_1[1] = Y[0]+alpha*(X[2]*A_1/2+X[1]*Lk) ;
    G_2[1] =      alpha*(X[2]*A_2/2+X[1]*ak) ;

  }

  void
  G2solve2arc::evalFJ( valueType const vars[2],
                       valueType       F[2],
                       valueType       J[2][2] ) const {

    valueType alpha = vars[0] ;
    valueType L     = vars[1] ;
    valueType G[2], G_1[2], G_2[2] ;

    evalG( alpha, L, th0, k0, G, G_1, G_2 ) ;

    F[0]    = G[0] - 2/L ;       F[1]    = G[1] ;
    J[0][0] = G_1[0] ;           J[1][0] = G_1[1] ;
    J[0][1] = G_2[0] + 2/(L*L) ; J[1][1] = G_2[1] ;

    evalG( alpha-1, L, th1, k1, G, G_1, G_2 ) ;
    F[0]    -= G[0] ;   F[1]    -= G[1] ;
    J[0][0] -= G_1[0] ; J[1][0] -= G_1[1] ;
    J[0][1] -= G_2[0] ; J[1][1] -= G_2[1] ;
  }
  
  // ---------------------------------------------------------------------------

  int
  G2solve2arc::solve() {
    Solve2x2 solver ;
    valueType X[2] = { 0.5, 2 } ;
    int iter = 0 ;
    bool converged = false ;
    do {
      valueType F[2], J[2][2], d[2] ;
      evalFJ( X, F, J ) ;
      if ( !solver.factorize( J ) ) break ;
      solver.solve( F, d ) ;
      valueType lenF = hypot(F[0],F[1]) ;
      X[0] -= d[0] ;
      X[1] -= d[1] ;
      converged = lenF < tolerance ;
    } while ( ++iter < maxIter && !converged ) ;
    if ( converged ) converged = X[1] > 0 && X[0] > 0 && X[0] < 1 ;
    if ( converged ) buildSolution( X[0], X[1] ) ;
    return converged ? iter : -1 ;
  }

  // **************************************************************************

  void
  G2solve2arc::buildSolution( valueType alpha, valueType L ) {
    valueType beta = 1-alpha ;
    valueType LL   = L*lambda ;
    valueType s0   = LL*alpha ;
    valueType s1   = LL*beta ;

    valueType tmp = k0*alpha+k1*beta-2*DeltaTheta/L ;

    valueType dk0 = -(k0+tmp)/(alpha*s0) ;
    valueType dk1 =  (k1+tmp)/(alpha*s1) ;

    S0.setup( x0, y0, theta0, kappa0, dk0,  0, s0 ) ;
    S1.setup( x1, y1, theta1, kappa1, dk1, -s1, 0 ) ;
    S1.change_origin( -s1 ) ;
  }

  /*\
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/

  void
  G2solve3arc::setTolerance( valueType tol ) {
    G2LIB_ASSERT( tol > 0 && tol <= 0.1,
                  "G2solve3arc::setTolerance, tolerance = " << tol <<
                  " must be in (0,0.1]" ) ;
    tolerance = tol ;
  }

  void
  G2solve3arc::setMaxIter( int miter ) {
    G2LIB_ASSERT( miter > 0 && miter <= 1000,
                  "G2solve3arc::setMaxIter, maxIter = " << miter <<
                  " must be in [1,1000]" ) ;
    maxIter = miter ;
  }

  int
  G2solve3arc::build( valueType _x0,
                      valueType _y0,
                      valueType _theta0,
                      valueType _kappa0,
                      valueType _x1,
                      valueType _y1,
                      valueType _theta1,
                      valueType _kappa1,
                      valueType L0,
                      valueType L1 ) {

    // use G1 for guess and length estimation
    Biarc::Biarc bi( _x0, _y0, _theta0, _x1, _y1, _theta1 ) ;
    valueType L3      = std::min(m_pi,bi.getC0().getL()+bi.getC1().getL())/3;
    valueType dth_max = m_pi ;
    if ( L0 <= 0 ) {
      valueType kappa = std::abs( _kappa0 + bi.getC0().getKappa() ) ;
      L0 = L3*kappa > dth_max ? dth_max / kappa : L3 ;
    }
    if ( L1 <= 0 ) {
      valueType kappa = std::abs( _kappa1 + bi.getC1().getKappa() ) ;
      L1 = L3*kappa > dth_max ? dth_max / kappa : L3 ;
      if ( L1*kappa > dth_max ) L1 = dth_max / kappa ;
    }

    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = _theta0;
    kappa0 = _kappa0 ;
    x1     = _x1 ;
    y1     = _y1 ;
    theta1 = _theta1 ;
    kappa1 = _kappa1 ;

    // scale problem
    valueType dx = x1 - x0 ;
    valueType dy = y1 - y0 ;
    phi    = atan2( dy, dx ) ;
    Lscale = 2/hypot( dx, dy ) ;

    th0 = theta0 - phi ;
    th1 = theta1 - phi ;

    valueType k0 = kappa0/Lscale ;
    valueType k1 = kappa1/Lscale ;

    s0 = L0*Lscale ;
    s1 = L1*Lscale ;

    K0 = k0*s0 ;
    K1 = k1*s1 ;

    valueType t0 = 2*th0+K0;
    valueType t1 = 2*th1-K1;

    c0  = s0*s1;
    c1  = 2 * s0;
    c2  = 0.25*((K1-6*K0-6*th0-2*th1)*s0 - 3*K0*s1);
    c3  = -c0 * (K0 + th0);
    c4  = 2 * s1;
    c5  = 0.25*((6*K1-K0-2*th0-6*th1)*s1 + 3*K1*s0);
    c6  = c0 * (K1 - th1);
    c7  = -0.5*(s0 + s1);
    c8  = th0 + th1 + 0.5*(K0 - K1);
    c9  = 0.25*(t1*s0 + t0*s1);
    c10 = 0.5*(s1 - s0);
    c11 = 0.5*(th1 - th0) - 0.25*(K0 + K1);
    c12 = 0.25*(t1*s0 - t0*s1);
    c13 = 0.5*s0*s1;
    c14 = 0.75*(s0 + s1);

    //return solve( 1.5*L3+0.5*(L0-L1), SM.theta(1.5*L3) ) ;
    //return solve( L3, SM.theta(1.5*L3), 3*L3 ) ;
    //return solve( L3, (theta0+theta1)/2, 3*L3 ) ;
    //return solve( L3, SM.theta(1.5*L3) ) ;
    return solve( 2, 0 ) ;
  }

  // **************************************************************************

  void
  G2solve3arc::evalFJ( valueType const vars[2],
                       valueType       F[2],
                       valueType       J[2][2]) const {

    valueType sM  = vars[0];
    valueType thM = vars[1];

    valueType dsM = 1.0 / (c13+(c14+sM)*sM);
    valueType sM2 = sM*sM;
    valueType stM = sM*thM;

    valueType dK0 = dsM*(c0*thM + c1*stM - K0*sM2 + c2*sM + c3);
    valueType dK1 = dsM*(c0*thM + c4*stM + K1*sM2 + c5*sM + c6);
    valueType dKM = dsM*sM*(c7*thM - 2*stM + c8*sM + c9);
    valueType KM  = dsM*sM*(c10*thM + c11*sM + c12);

    valueType X0[3],  Y0[3],
              X1[3],  Y1[3],
              XMp[3], YMp[3],
              XMm[3], YMm[3];
    GeneralizedFresnelCS( 3, dK0,  K0, th0, X0,  Y0);
    GeneralizedFresnelCS( 3, dK1, -K1, th1, X1,  Y1);
    GeneralizedFresnelCS( 3, dKM,  KM, thM, XMp, YMp);
    GeneralizedFresnelCS( 3, dKM, -KM, thM, XMm, YMm);

    // in the standard problem dx = 2, dy = 0
    //valueType dx  = x1 - x0;
    //valueType dy  = y1 - y0;
    F[0] = s0*X0[0] + s1*X1[0] + sM*(XMm[0] + XMp[0]) - 2 ;
    F[1] = s0*Y0[0] + s1*Y1[0] + sM*(YMm[0] + YMp[0]) - 0 ;

    // calcolo J(F)
    valueType dsM2 = dsM*dsM;
    valueType g0 = -(2 * sM + c14)*dsM2;
    valueType g1 = (c13 - sM2)*dsM2;
    valueType g2 = sM*(sM*c14+2*c13)*dsM2;

    valueType dK0_sM  = (c0*thM+c3)*g0 + (c1*thM+c2)*g1 - K0*g2 ;
    valueType dK1_sM  = (c0*thM+c6)*g0 + (c4*thM+c5)*g1 + K1*g2 ;
    valueType dKM_sM  = (c7*thM+c9)*g1 + (c8-2*thM)*g2;
    valueType KM_sM   = (c10*thM+c12)*g1 + c11*g2;

    valueType dK0_thM = (c0+c1*sM)*dsM;
    valueType dK1_thM = (c0+c4*sM)*dsM;
    valueType dKM_thM = (c7-2*sM)*sM*dsM;
    valueType KM_thM  = c10*sM*dsM;

    // coeff fresnel per f_j per lo jacobiano
    valueType f0 = -0.5*s0*Y0[2];
    valueType f1 = -0.5*s1*Y1[2];
    valueType f2 = -0.5*sM*(YMm[2] + YMp[2]);
    valueType f3 = sM*(YMm[1] - YMp[1]);
    valueType f4 = 0.5*s0*X0[2];
    valueType f5 = 0.5*s1*X1[2];
    valueType f6 = 0.5*sM*(XMm[2] + XMp[2]);
    valueType f7 = sM*(XMp[1] - XMm[1]);
    valueType f8 = XMp[0]+XMm[0];
    valueType f9 = YMp[0]+YMm[0];

    J[0][0] = f0 * dK0_sM  + f1 * dK1_sM  + f2 * dKM_sM  + f3 * KM_sM  + f8 ;
    J[0][1] = f0 * dK0_thM + f1 * dK1_thM + f2 * dKM_thM + f3 * KM_thM - sM * f9 ;
    J[1][0] = f4 * dK0_sM  + f5 * dK1_sM  + f6 * dKM_sM  + f7 * KM_sM  + f9 ;
    J[1][1] = f4 * dK0_thM + f5 * dK1_thM + f6 * dKM_thM + f7 * KM_thM + sM * f8 ;
  }

  // **************************************************************************

  int
  G2solve3arc::solve( valueType sM_guess,
                      valueType thM_guess ) {

    Solve2x2 solver;
    valueType F[2], J[2][2], d[2], X[2];
    X[0] = sM_guess ;
    X[1] = thM_guess ;

    int iter = 0;
    bool converged = false;
    try {
      do {
        evalFJ(X, F, J);
        if (!solver.factorize(J)) break;
        solver.solve(F, d);
        X[0] -= d[0];
        X[1] -= d[1];
        valueType lenF = hypot(F[0], F[1]);
        converged = lenF < tolerance;
      } while ( ++iter < maxIter && !converged );

      // re-check solution
      if ( converged )
        converged = FP_INFINITE != std::fpclassify(X[0]) &&
                    FP_NAN      != std::fpclassify(X[0]) &&
                    FP_INFINITE != std::fpclassify(X[1]) &&
                    FP_NAN      != std::fpclassify(X[1]) &&
                    X[0] > 0;
    }
    catch (...) {
      // nothing to do
    }
    if ( converged ) buildSolution(X[0], X[1]); // costruisco comunque soluzione
    return converged ? iter : -1;
  }

  // **************************************************************************
  static
  inline
  valueType
  power2( valueType a )
  { return a*a ; }

  void 
  G2solve3arc::buildSolution( valueType sM, valueType thM ) {

    // ricostruzione dati clotoidi trasformati
    valueType dsM = 1.0 / (c13+(c14+sM)*sM);
    valueType sM2 = sM*sM;
    valueType stM = sM*thM;

    valueType dK0 = dsM*(c0*thM + c1*stM - K0*sM2 + c2*sM + c3);
    valueType dK1 = dsM*(c0*thM + c4*stM + K1*sM2 + c5*sM + c6);
    valueType dKM = dsM*sM*(c7*thM - 2*stM + c8*sM + c9);
    valueType KM  = dsM*sM*(c10*thM + c11*sM + c12);

    valueType xa, ya, xmL, ymL;
    GeneralizedFresnelCS( dK0,  K0, th0, xa,  ya  );
    GeneralizedFresnelCS( dKM, -KM, thM, xmL, ymL );

    valueType xM = s0 * xa + sM * xmL - 1;
    valueType yM = s0 * ya + sM * ymL;

    // rovescia trasformazione standard
    valueType L0 = s0/Lscale;
    valueType L1 = s1/Lscale;
    valueType LM = sM/Lscale;

    dK0 *= power2(Lscale/s0) ;
    dK1 *= power2(Lscale/s1) ;
    dKM *= power2(Lscale/sM) ;
    KM  *= Lscale/sM ;

    S0.setup( x0, y0, theta0, kappa0, dK0,   0, L0 );
    S1.setup( x1, y1, theta1, kappa1, dK1, -L1, 0  );

    // la trasformazione inversa da [-1,1] a (x0,y0)-(x1,y1)
    // g(x,y) = RotInv(phi)*(1/lambda*[X;Y] - [xbar;ybar]) = [x;y]

    valueType C  = cos(phi);
    valueType S  = sin(phi);
    valueType dx = (xM + 1) / Lscale;
    valueType dy = yM / Lscale;
    SM.setup( x0 + C * dx - S * dy,
              y0 + C * dy + S * dx,
              thM + phi, KM, dKM, -LM, LM );

    S1.change_origin(-L1);
    SM.change_origin(-LM);
  }

  // **************************************************************************

  valueType
  G2solve3arc::thetaMinMax( valueType & thMin, valueType & thMax ) const {
    valueType thMin1, thMax1 ;
    S0.thetaMinMax( thMin,  thMax ) ;
    S1.thetaMinMax( thMin1, thMax1 ) ;
    if ( thMin > thMin1 ) thMin = thMin1 ;
    if ( thMax < thMax1 ) thMax = thMax1 ;
    SM.thetaMinMax( thMin1, thMax1 ) ;
    if ( thMin > thMin1 ) thMin = thMin1 ;
    if ( thMax < thMax1 ) thMax = thMax1 ;
    return thMax-thMin ;
  }

  // **************************************************************************

  valueType
  G2solve3arc::curvatureMinMax( valueType & kMin, valueType & kMax ) const {
    valueType kMin1, kMax1 ;
    S0.thetaMinMax( kMin,  kMax ) ;
    S1.thetaMinMax( kMin1, kMax1 ) ;
    if ( kMin > kMin1 ) kMin = kMin1 ;
    if ( kMax < kMax1 ) kMax = kMax1 ;
    SM.thetaMinMax( kMin1, kMax1 ) ;
    if ( kMin > kMin1 ) kMin = kMin1 ;
    if ( kMax < kMax1 ) kMax = kMax1 ;
    return kMax-kMin ;
  }

}

// EOF: ClothoidG2.cc
