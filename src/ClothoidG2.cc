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

namespace G2lib {

  using namespace std ;

  inline
  valueType
  power2( valueType a )
  { return a*a ; }

  inline
  valueType
  power3( valueType a )
  { return a*a*a ; }

  inline
  valueType
  power4( valueType a )
  { valueType a2 = a*a ; return a2*a2 ; }

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/

  int
  G2solve2arc::build( valueType _x0,
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

    DeltaK     = k1  - k0 ;
    DeltaTheta = th1 - th0 ;

    return solve();
  }

  void
  G2solve2arc::setTolerance( valueType tol ) {
    G2LIB_ASSERT( tol > 0 && tol <= 0.1,
                  "G2solve2arc::setTolerance, tolerance = " << tol << " must be in (0,0.1]" ) ;
    tolerance = tol ;
  }

  void
  G2solve2arc::setMaxIter( int miter ) {
    G2LIB_ASSERT( miter > 0 && miter <= 1000,
                  "G2solve2arc::setMaxIter, maxIter = " << miter << " must be in [1,1000]" ) ;
    maxIter = miter ;
  }

  void
  G2solve2arc::evalA( valueType   alpha,
                      valueType   L,
                      valueType & A ) const {
    valueType K  = k0+k1 ;
    valueType aK = alpha*DeltaK ;
    A = alpha*(L*(aK-K)+2*DeltaTheta) ;
  }

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
                      valueType G[2] ) const {
    valueType A, X, Y ;
    evalA( alpha, L, A ) ;
    valueType ak = alpha*k ;
    GeneralizedFresnelCS( A, ak*L, th, X, Y );
    G[0] = alpha*X ;
    G[1] = alpha*Y ;
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
  G2solve2arc::evalF( valueType const vars[2], valueType F[2] ) const {
    valueType alpha = vars[0] ;
    valueType L     = vars[1] ;
    valueType G[2] ;
    evalG( alpha, L, th0, k0, G ) ;
    F[0] = G[0] - 2/L ;  F[1] = G[1] ;
    evalG( alpha-1, L, th1, k1, G ) ;
    F[0] -= G[0] ; F[1] -= G[1] ;
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
      #if 0
      X[0] -= d[0];
      X[1] -= d[1];
      #else
      valueType FF[2], dd[2], XX[2];
      // Affine invariant Newton solver
      valueType nd = hypot( d[0], d[1] ) ;
      bool step_found = false ;
      valueType tau = 2 ;
      do {
        tau  /= 2 ;
        XX[0] = X[0]-tau*d[0];
        XX[1] = X[1]-tau*d[1];
        evalF(XX, FF);
        solver.solve(FF, dd);
        step_found = hypot( dd[0], dd[1] ) <= (1-tau/2)*nd + 1e-6
                     && XX[0] > 0 && XX[0] < 1 && XX[1] > 0 ;
      } while ( tau > 1e-6 && !step_found );
      if ( !step_found ) break ;
      X[0] = XX[0];
      X[1] = XX[1];
      #endif
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
    valueType s0   = L*alpha ;
    valueType s1   = L*beta ;
    valueType tmp  = 2*DeltaTheta-L*(k0+k1);
    valueType A0   = alpha*(s0*DeltaK+tmp) ;
    valueType A1   = beta*(s1*DeltaK-tmp) ;

    valueType dk0  = A0/(s0*s0) ;
    valueType dk1  = A1/(s1*s1) ;

    // transform solution from (-1,0)--(1,0) to (x0,y0)--(x1,y1)
    //S0.build( -1, 0, th0, k0, dk0, s0 ) ;
    //S1.build(  1, 0, th1, k1, dk1, s1 ) ;
    //S1.changeCurvilinearOrigin( -s1, s1 ) ;
    s0  *= lambda ;
    s1  *= lambda ;
    dk0 /= lambda*lambda ;
    dk1 /= lambda*lambda ;

    S0.build( x0, y0, theta0, kappa0, dk0, s0 ) ;
    S1.build( x1, y1, theta1, kappa1, dk1, s1 ) ;
    S1.changeCurvilinearOrigin( -s1, s1 ) ;
  }

  /*\
   |    ____ ____            _            ____ _     ____
   |   / ___|___ \ ___  ___ | |_   _____ / ___| |   / ___|
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |   | |  | |
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/ |___| |__| |___
   |   \____|_____|___/\___/|_| \_/ \___|\____|_____\____|
  \*/

  int
  G2solveCLC::build( valueType _x0,
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

    return solve();
  }

  void
  G2solveCLC::setTolerance( valueType tol ) {
    G2LIB_ASSERT( tol > 0 && tol <= 0.1,
                  "G2solveCLC::setTolerance, tolerance = " << tol << " must be in (0,0.1]" ) ;
    tolerance = tol ;
  }

  void
  G2solveCLC::setMaxIter( int miter ) {
    G2LIB_ASSERT( miter > 0 && miter <= 1000,
                  "G2solveCLC::setMaxIter, maxIter = " << miter << " must be in [1,1000]" ) ;
    maxIter = miter ;
  }

  // ---------------------------------------------------------------------------

  int
  G2solveCLC::solve() {
    valueType X0[3], Y0[3], X1[3], Y1[3] ;
      valueType thM = 0, sM = 0.0 ;
    int iter = 0 ;
    bool converged = false ;
    do {
      valueType D0 = thM - th0 ;
      valueType D1 = thM - th1 ;

      GeneralizedFresnelCS( 3, 2*D0, -2*D0, D0, X0, Y0 ) ;
      GeneralizedFresnelCS( 3, 2*D1, -2*D1, D1, X1, Y1 ) ;

      valueType F  = D0*k1*Y0[0]-D1*k0*Y1[0] - k0*k1*sin(thM) ;
      valueType dF = D0*k1*(X0[2]-2*X0[1]+X0[0])-D1*k0*(X1[2]-2*X1[1]+X1[0]) - k0*k1*cos(thM) + k1*Y0[0]-k0*Y1[0];

      if ( std::abs(dF) < 1e-10 ) break ;
      valueType d = F/dF ;
      #if 0
      thM -= d;
      #else
      valueType FF, dd, thM1 ;
      // Affine invariant Newton solver
      bool step_found = false ;
      valueType tau = 2 ;
      do {
        tau  /= 2 ;
        thM1 = thM-tau*d;
        D0 = thM1 - th0 ;
        D1 = thM1 - th1 ;
        GeneralizedFresnelCS( 1, 2*D0, -2*D0, D0, X0, Y0 ) ;
        GeneralizedFresnelCS( 1, 2*D1, -2*D1, D1, X1, Y1 ) ;
        FF = D0*k1*Y0[0]-D1*k0*Y1[0] - k0*k1*sin(thM1) ;
        dd = FF/dF ;
        step_found = std::abs( dd ) <= (1-tau/2)*std::abs(d) + 1e-6 ;
      } while ( tau > 1e-6 && !step_found );
      if ( !step_found ) break ;
      thM = thM1 ;
      #endif
      converged = std::abs(d) < tolerance ;
    } while ( ++iter < maxIter && !converged ) ;
    if ( converged ) {
      valueType D0 = thM - th0 ;
      valueType D1 = thM - th1 ;
      GeneralizedFresnelCS( 1, 2*D0, -2*D0, D0, X0, Y0 ) ;
      GeneralizedFresnelCS( 1, 2*D1, -2*D1, D1, X1, Y1 ) ;
      sM = cos(thM) + D1*X1[0]/k1 - D0*X0[0]/k0 ;
      converged = sM > 0 && sM < 1e100 ;
    }
    if ( converged ) converged = buildSolution( sM, thM ) ;
    return converged ? iter : -1 ;
  }

  // **************************************************************************

  bool
  G2solveCLC::buildSolution( valueType sM, valueType thM ) {
    valueType dk0 = 0.5*power2(k0/lambda)/(th0-thM) ;
    valueType dk1 = 0.5*power2(k1/lambda)/(th1-thM) ;
    valueType L0  = 2*lambda*(thM-th0)/k0 ;
    valueType L1  = 2*lambda*(th1-thM)/k1 ;

    if ( ! ( L0 > 0 && L1 > 0 ) ) return false ;

    S0.build( x0, y0, theta0, kappa0, dk0, L0 ) ;
    S1.build( x1, y1, theta1, kappa1, dk1, L1 ) ;
    S1.changeCurvilinearOrigin( -L1, L1 ) ;
    SM.build( S0.xEnd(), S0.yEnd(), S0.thetaEnd(), 0, 0, 2*sM*lambda ) ;

    return true ;
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
                      valueType Dmax,
                      valueType dmax ) {
    try {
      // save data
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0;
      kappa0 = _kappa0 ;
      x1     = _x1 ;
      y1     = _y1 ;
      theta1 = _theta1 ;
      kappa1 = _kappa1 ;

      // transform to reference frame
      valueType dx = x1 - x0 ;
      valueType dy = y1 - y0 ;
      phi    = atan2( dy, dx ) ;
      Lscale = 2/hypot( dx, dy ) ;

      th0 = theta0 - phi ;
      th1 = theta1 - phi ;

      // put in range
      rangeSymm(th0);
      rangeSymm(th1);

      K0 = (kappa0/Lscale) ; // k0
      K1 = (kappa1/Lscale) ; // k1

      if ( Dmax <= 0 ) Dmax = m_pi ;
      if ( dmax <= 0 ) dmax = m_pi/8 ;

      if ( Dmax > m_2pi  ) Dmax = m_2pi ;
      if ( dmax > m_pi/4 ) dmax = m_pi/4 ;

      // compute guess G1
      ClothoidCurve SG ;
      SG.build_G1( -1, 0, th0, 1, 0, th1 ) ;

      valueType kA = SG.kappaBegin() ;
      valueType kB = SG.kappaEnd() ;
      valueType dk = std::abs(SG.kappa_D()) ;
      valueType L3 = SG.length()/3 ;

      valueType tmp = 0.5*std::abs(K0-kA)/dmax;
      s0 = L3;
      if ( tmp*s0 > 1 ) s0 = 1/tmp ;
      tmp = (std::abs(K0+kA)+s0*dk)/(2*Dmax) ;
      if ( tmp*s0 > 1 ) s0 = 1/tmp ;

      tmp = 0.5*std::abs(K1-kB)/dmax ;
      s1 = L3;
      if ( tmp*s1 > 1 ) s1 = 1/tmp ;
      tmp = (std::abs(K1+kB)+s1*dk)/(2*Dmax) ;
      if ( tmp*s1 > 1 ) s1 = 1/tmp ;

      valueType dth   = std::abs(th0-th1) / m_2pi ;
      valueType scale = power3(cos( power4(dth)*m_pi_2 )) ;
      s0 *= scale ;
      s1 *= scale ;

      valueType L   = (3*L3-s0-s1)/2 ;
      valueType thM = SG.theta(s0+L) ;
      th0 = SG.thetaBegin() ;
      th1 = SG.thetaEnd() ;

      // setup

      K0 *= s0;
      K1 *= s1;

      valueType t0 = 2*th0+K0;
      valueType t1 = 2*th1-K1;

      c0  = s0*s1;
      c1  = 2 * s0;
      c2  = 0.25*((K1-6*(K0+th0)-2*th1)*s0 - 3*K0*s1);
      c3  = -c0 * (K0 + th0);
      c4  = 2 * s1;
      c5  = 0.25*((6*(K1-th1)-K0-2*th0)*s1 + 3*K1*s0);
      c6  = c0 * (K1 - th1);
      c7  = -0.5*(s0 + s1);
      c8  = th0 + th1 + 0.5*(K0 - K1);
      c9  = 0.25*(t1*s0 + t0*s1);
      c10 = 0.5*(s1 - s0);
      c11 = 0.5*(th1 - th0) - 0.25*(K0 + K1);
      c12 = 0.25*(t1*s0 - t0*s1);
      c13 = 0.5*s0*s1;
      c14 = 0.75*(s0 + s1);
      return solve( L, thM ) ;
    } catch (...) {
      return -1 ;
      // nothing to do
    }
  }

  // **************************************************************************

  void
  G2solve3arc::evalF( valueType const vars[2], valueType F[2] ) const {

    valueType sM  = vars[0];
    valueType thM = vars[1];

    valueType dsM = 1.0 / (c13+(c14+sM)*sM);
    valueType dK0 = dsM*(c0*thM + sM*(c1*thM - K0*sM + c2) + c3);
    valueType dK1 = dsM*(c0*thM + sM*(c4*thM + K1*sM + c5) + c6);
    valueType dKM = dsM*sM*( thM*(c7-2*sM) + c8*sM + c9);
    valueType KM  = dsM*sM*(c10*thM + c11*sM + c12);

    valueType X0, Y0, X1, Y1, XMp, YMp, XMm, YMm ;
    GeneralizedFresnelCS( dK0,  K0, th0, X0,  Y0);
    GeneralizedFresnelCS( dK1, -K1, th1, X1,  Y1);
    GeneralizedFresnelCS( dKM,  KM, thM, XMp, YMp);
    GeneralizedFresnelCS( dKM, -KM, thM, XMm, YMm);

    // in the standard problem dx = 2, dy = 0
    //valueType dx  = x1 - x0;
    //valueType dy  = y1 - y0;
    F[0] = s0*X0 + s1*X1 + sM*(XMm + XMp) - 2 ;
    F[1] = s0*Y0 + s1*Y1 + sM*(YMm + YMp) - 0 ;
  }


  // **************************************************************************

  void
  G2solve3arc::evalFJ( valueType const vars[2],
                       valueType       F[2],
                       valueType       J[2][2]) const {

    valueType sM  = vars[0];
    valueType thM = vars[1];

    valueType dsM   = 1.0 / (c13+(c14+sM)*sM);
    valueType dsMsM = dsM*sM;
    valueType dK0   = dsM*(c0*thM + sM*(c1*thM + c2 - sM*K0) + c3);
    valueType dK1   = dsM*(c0*thM + sM*(c4*thM + c5 + sM*K1) + c6);
    valueType dKM   = dsMsM*(thM*(c7-2*sM) + c8*sM + c9);
    valueType KM    = dsMsM*(c10*thM + c11*sM + c12);

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
    valueType t0 = XMp[0]+XMm[0];
    valueType t1 = YMp[0]+YMm[0];
    F[0] = s0*X0[0] + s1*X1[0] + sM*t0 - 2 ;
    F[1] = s0*Y0[0] + s1*Y1[0] + sM*t1 - 0 ;

    // calcolo J(F)
    valueType dsM2 = dsM*dsM;
    valueType g0   = -(2 * sM + c14)*dsM2;
    valueType g1   = (c13 - sM*sM)*dsM2;
    valueType g2   = sM*(sM*c14+2*c13)*dsM2;

    valueType dK0_sM  = (c0*thM+c3)*g0 + (c1*thM+c2)*g1 - K0*g2 ;
    valueType dK1_sM  = (c0*thM+c6)*g0 + (c4*thM+c5)*g1 + K1*g2 ;
    valueType dKM_sM  = (c7*thM+c9)*g1 + (c8-2*thM)*g2;
    valueType KM_sM   = (c10*thM+c12)*g1 + c11*g2;

    valueType dK0_thM = (c0+c1*sM)*dsM;
    valueType dK1_thM = (c0+c4*sM)*dsM;
    valueType dKM_thM = (c7-2*sM)*dsMsM;
    valueType KM_thM  = c10*dsMsM;

    // coeff fresnel per f_j per lo jacobiano
    valueType f0 = -0.5*s0*Y0[2];
    valueType f1 = -0.5*s1*Y1[2];
    valueType f2 = -0.5*sM*(YMm[2] + YMp[2]);
    valueType f3 = sM*(YMm[1] - YMp[1]);
    valueType f4 = 0.5*s0*X0[2];
    valueType f5 = 0.5*s1*X1[2];
    valueType f6 = 0.5*sM*(XMm[2] + XMp[2]);
    valueType f7 = sM*(XMp[1] - XMm[1]);

    J[0][0] = f0 * dK0_sM  + f1 * dK1_sM  + f2 * dKM_sM  + f3 * KM_sM  + t0 ;
    J[0][1] = f0 * dK0_thM + f1 * dK1_thM + f2 * dKM_thM + f3 * KM_thM - sM * t1 ;
    J[1][0] = f4 * dK0_sM  + f5 * dK1_sM  + f6 * dKM_sM  + f7 * KM_sM  + t1 ;
    J[1][1] = f4 * dK0_thM + f5 * dK1_thM + f6 * dKM_thM + f7 * KM_thM + sM * t0 ;
  }

  // **************************************************************************

  int
  G2solve3arc::solve( valueType sM_guess,
                      valueType thM_guess ) {

    Solve2x2 solver;
    valueType F[2], d[2], X[2], J[2][2];
    X[0] = sM_guess ;
    X[1] = thM_guess ;

    //valueType thmin = min(th0,th1)-2*m_2pi ;
    //valueType thmax = max(th0,th1)+2*m_2pi ;

    int iter = 0;
    bool converged = false;
    try {
      do {
        evalFJ(X, F, J);
        valueType lenF = hypot(F[0], F[1]);
        converged = lenF < tolerance;
        if ( converged || !solver.factorize(J) ) break;
        solver.solve(F, d);
        #if 1
        X[0] -= d[0];
        X[1] -= d[1];
        #else
        valueType FF[2], dd[2], XX[2];
        // Affine invariant Newton solver
        valueType nd = hypot( d[0], d[1] ) ;
        bool step_found = false ;
        valueType tau = 2 ;
        do {
          tau  /= 2 ;
          XX[0] = X[0]-tau*d[0];
          XX[1] = X[1]-tau*d[1];
          evalF(XX, FF);
          solver.solve(FF, dd);
          step_found = hypot( dd[0], dd[1] ) <= (1-tau/2)*nd + 1e-6 ;
                       //&& XX[0] > 0 ; // && XX[0] > X[0]/4 && XX[0] < 4*X[0] ;
                       //&& XX[1] > thmin && XX[1] < thmax ;
        } while ( tau > 1e-6 && !step_found );
        if ( !step_found ) break ;
        X[0] = XX[0];
        X[1] = XX[1];
        #endif
      } while ( ++iter < maxIter );

      // re-check solution
      if ( converged )
        converged = FP_INFINITE != std::fpclassify(X[0]) &&
                    FP_NAN      != std::fpclassify(X[0]) &&
                    FP_INFINITE != std::fpclassify(X[1]) &&
                    FP_NAN      != std::fpclassify(X[1]) ;
    }
    catch (...) {
      cout << "PASSA\n" ;
      // nothing to do
    }
    if ( converged ) buildSolution(X[0], X[1]); // costruisco comunque soluzione
    return converged ? iter : -1;
  }

  // **************************************************************************

  void 
  G2solve3arc::buildSolution( valueType sM, valueType thM ) {
    // soluzione nel frame di riferimento
    /* valueType k0 = K0
     S0.build( -1, 0, th0, k0, dK0,   0, L0 );
     S1.build( x1, y1, phi+th1, kappa1, dK1, -L1, 0  );
     S1.change_origin(-L1);
    */

    // ricostruzione dati clotoidi trasformati
    valueType dsM = 1.0 / (c13+(c14+sM)*sM);
    valueType dK0 = dsM*(c0*thM + sM*(c1*thM - K0*sM + c2) + c3);
    valueType dK1 = dsM*(c0*thM + sM*(c4*thM + K1*sM + c5) + c6);
    valueType dKM = dsM*sM*(c7*thM + sM*(c8 - 2*thM) + c9);
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

    //th0 = theta0 - phi ;
    //th1 = theta1 - phi ;
    S0.build( x0, y0, phi+th0, kappa0, dK0, L0 );
    S1.build( x1, y1, phi+th1, kappa1, dK1, L1 );
    S1.changeCurvilinearOrigin( -L1, L1 );

    // la trasformazione inversa da [-1,1] a (x0,y0)-(x1,y1)
    // g(x,y) = RotInv(phi)*(1/lambda*[X;Y] - [xbar;ybar]) = [x;y]

    valueType C  = cos(phi);
    valueType S  = sin(phi);
    valueType dx = (xM + 1) / Lscale;
    valueType dy = yM / Lscale;
    SM.build( x0 + C * dx - S * dy,
              y0 + C * dy + S * dx,
              thM + phi, KM, dKM, 2*LM );
    SM.changeCurvilinearOrigin( -LM, 2*LM );
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
    S0.curvatureMinMax( kMin,  kMax ) ;
    S1.curvatureMinMax( kMin1, kMax1 ) ;
    if ( kMin > kMin1 ) kMin = kMin1 ;
    if ( kMax < kMax1 ) kMax = kMax1 ;
    SM.curvatureMinMax( kMin1, kMax1 ) ;
    if ( kMin > kMin1 ) kMin = kMin1 ;
    if ( kMax < kMax1 ) kMax = kMax1 ;
    return kMax-kMin ;
  }

  valueType
  G2solve3arc::theta( valueType s ) const {
    if ( s < S0.length() ) return S0.theta(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.theta(s) ;
    s -= S0.length() ;
    return S1.theta(s) ;
  }

  valueType
  G2solve3arc::theta_D( valueType s ) const {
    if ( s < S0.length() ) return S0.theta_D(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.theta_D(s) ;
    s -= S0.length() ;
    return S1.theta_D(s) ;
  }

  valueType
  G2solve3arc::theta_DD( valueType s ) const {
    if ( s < S0.length() ) return S0.theta_DD(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.theta_DD(s) ;
    s -= S0.length() ;
    return S1.theta_DD(s) ;
  }

  valueType
  G2solve3arc::theta_DDD( valueType s ) const {
    if ( s < S0.length() ) return S0.theta_DDD(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.theta_DDD(s) ;
    s -= S0.length() ;
    return S1.theta_DDD(s) ;
  }

  valueType
  G2solve3arc::X( valueType s ) const {
    if ( s < S0.length() ) return S0.X(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.X(s) ;
    s -= S0.length() ;
    return S1.X(s) ;
  }

  valueType
  G2solve3arc::Y( valueType s ) const {
    if ( s < S0.length() ) return S0.Y(s) ;
    s -= S0.length() ;
    if ( s < SM.length() ) return SM.Y(s) ;
    s -= S0.length() ;
    return S1.Y(s) ;
  }

  void
  G2solve3arc::eval( valueType   s,
                     valueType & theta,
                     valueType & kappa,
                     valueType & x,
                     valueType & y ) const {
    if ( s < S0.length() ) {
      S0.eval(s, theta, kappa, x, y ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval(s, theta, kappa, x, y ) ;
      } else {
        s -= SM.length() ;
        S1.eval(s, theta, kappa, x, y ) ;
      }
    }
  }

  void
  G2solve3arc::eval( valueType s, valueType & x, valueType & y ) const {
    if ( s < S0.length() ) {
      S0.eval(s, x, y ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval(s, x, y ) ;
      } else {
        s -= SM.length() ;
        S1.eval(s, x, y ) ;
      }
    }
  }

  void
  G2solve3arc::eval_D( valueType s, valueType & x_D, valueType & y_D ) const {
    if ( s < S0.length() ) {
      S0.eval_D(s, x_D, y_D ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_D(s, x_D, y_D ) ;
      } else {
        s -= SM.length() ;
        S1.eval_D(s, x_D, y_D ) ;
      }
    }
  }

  void
  G2solve3arc::eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const {
    if ( s < S0.length() ) {
      S0.eval_DD(s, x_DD, y_DD ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_DD(s, x_DD, y_DD ) ;
      } else {
        s -= SM.length() ;
        S1.eval_DD(s, x_DD, y_DD ) ;
      }
    }
  }

  void
  G2solve3arc::eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const {
    if ( s < S0.length() ) {
      S0.eval_DDD(s, x_DDD, y_DDD ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_DDD(s, x_DDD, y_DDD ) ;
      } else {
        s -= SM.length() ;
        S1.eval_DDD(s, x_DDD, y_DDD ) ;
      }
    }
  }

  // offset curve
  void
  G2solve3arc::eval( valueType s, valueType offs, valueType & x, valueType & y ) const {
    if ( s < S0.length() ) {
      S0.eval(s, offs, x, y ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval(s, offs, x, y ) ;
      } else {
        s -= SM.length() ;
        S1.eval(s, offs, x, y ) ;
      }
    }
  }

  void
  G2solve3arc::eval_D( valueType s, valueType offs, valueType & x_D, valueType & y_D ) const {
    if ( s < S0.length() ) {
      S0.eval_D(s, offs, x_D, y_D ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_D(s, offs, x_D, y_D ) ;
      } else {
        s -= SM.length() ;
        S1.eval_D(s, offs, x_D, y_D ) ;
      }
    }
  }

  void
  G2solve3arc::eval_DD( valueType s, valueType offs, valueType & x_DD, valueType & y_DD ) const {
    if ( s < S0.length() ) {
      S0.eval_DD(s, offs, x_DD, y_DD ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_DD(s, offs, x_DD, y_DD ) ;
      } else {
        s -= SM.length() ;
        S1.eval_DD(s, offs, x_DD, y_DD ) ;
      }
    }
  }

  void
  G2solve3arc::eval_DDD( valueType s, valueType offs, valueType & x_DDD, valueType & y_DDD ) const {
    if ( s < S0.length() ) {
      S0.eval_DDD(s, offs, x_DDD, y_DDD ) ;
    } else {
      s -= S0.length() ;
      if ( s < SM.length() ) {
        SM.eval_DDD(s, offs, x_DDD, y_DDD ) ;
      } else {
        s -= SM.length() ;
        S1.eval_DDD(s, offs, x_DDD, y_DDD ) ;
      }
    }
  }

}

// EOF: ClothoidG2.cc
