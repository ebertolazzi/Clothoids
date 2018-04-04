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

#include "Circle.hh"
#include "Clothoid.hh"
#include "CubicRootsFlocke.hh"

#include <cmath>
#include <cfloat>

namespace G2lib {

  using namespace std ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  bool
  closest4( valueType            epsi,
            ClothoidData const & CD,
            valueType            L,
            valueType            qx,
            valueType            qy,
            valueType          & X,
            valueType          & Y,
            valueType          & S ) {

    // S = GUESS
    int nb = 0 ;
    valueType theta, kappa, dS ;
    for ( int iter = 0 ; iter < 20 && nb < 2 ; ++iter ) {
      // approx clothoid with a circle
      CD.eval( S, theta, kappa, X, Y );

      valueType CS = cos(theta) ;
      valueType SS = sin(theta) ;

      valueType dx  = X - qx ;
      valueType dy  = Y - qy ;
      valueType a0  = CS * dy - SS * dx ;
      valueType b0  = SS * dy + CS * dx ;
      valueType tmp = a0*kappa ;

      if ( 1+2*tmp > 0 ) {

        tmp = b0/(1+tmp) ;
        dS = -tmp*Atanc(tmp*kappa) ;

      } else {

        valueType om = atan2( b0, a0+1/kappa ) ;
        if ( kappa < 0 ) {
          if ( om < 0 ) om += m_pi ;
          else          om -= m_pi ;
        }

        dS = -om/kappa ;
      }

      S += dS ;
      if ( abs( dS ) < epsi ) {
        if ( S < 0 || S > L ) break ;
        CD.eval( S, X, Y );
        return true ;
      }

      // check divergence
      if ( S < 0 || S > L ) ++nb ; else nb = 0 ;

    }
    return false ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closest3( valueType            epsi,
            ClothoidData const & CD,
            valueType            L,
            valueType            qx,
            valueType            qy,
            valueType          & X,
            valueType          & Y,
            valueType          & S ) {
    // S=0
    valueType dx0   = CD.x0-qx ;
    valueType dy0   = CD.y0-qy ;
    valueType tphi0 = CD.theta0 - atan2( dy0, dx0 ) ;
    valueType d0, X0, Y0, S0 = 0 ;
    bool min0 = cos(tphi0) > 0 ; // distance increasing
    if ( !min0 ) min0 = !closest4( epsi, CD, L, qx, qy, X0, Y0, S0 ) ;
    if ( min0 ) {
      S0 = 0 ; d0 = hypot( dx0, dy0 ) ;
    } else {
      d0 = hypot( X0-qx, Y0-qy ) ;
    }

    // S=L
    valueType thetaL, kappaL ;
    CD.eval( L, thetaL, kappaL, X, Y );
    valueType dxL   = X-qx ;
    valueType dyL   = Y-qy ;
    valueType tphiL = thetaL - atan2( dyL, dxL ) ;
    valueType dL, XL, YL, SL = L ;
    bool minL = cos(tphiL) < 0 ; // distance increasing
    if ( !minL ) minL = !closest4( epsi, CD, L, qx, qy, XL, YL, SL ) ;
    if ( minL ) {
      SL = L ; dL = hypot( dxL, dyL ) ;
    } else {
      dL = hypot( XL-qx, YL-qy ) ;
    }

    if ( min0 && minL ) {
      S = L/2 ;
      if ( closest4( epsi, CD, L, qx, qy, X, Y, S ) ) {
        valueType dx = X-qx ;
        valueType dy = Y-qy ;
        valueType d  = hypot( dx, dy ) ;
        if ( d < d0 && d < dL ) return d ;
      }
    }

    if ( dL < d0 ) {
      S = SL ; X = XL ; Y = YL ; return dL ;
    } else {
      S = S0 ; X = X0 ; Y = Y0 ; return d0 ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closest2( valueType            epsi,
            ClothoidData const & CD,
            valueType            L,
            valueType            qx,
            valueType            qy,
            valueType          & X,
            valueType          & Y,
            valueType          & S ) {

    valueType theta1 = CD.theta(L) ;
    if ( abs(theta1-CD.theta0) <= m_2pi )
      return closest3( epsi, CD, L, qx, qy, X, Y, S ) ;

    valueType Xi, Yi ; // point at infinity
    CD.Pinfinity( Xi, Yi, true ) ;
    valueType di   = hypot( qx-Xi, qy-Yi ) ;
    valueType d0   = hypot( qx-CD.x0, qy-CD.x0 ) ;
    valueType pi4  = 4*m_pi ;
    valueType adk  = abs(CD.dk) ;
    valueType pidk = pi4*adk;

    if ( di >= d0 ) {
      valueType LS = pi4/(abs(CD.k0)+sqrt(pidk+CD.k0*CD.k0)) ;
      return closest3( epsi, CD, LS, qx, qy, X, Y, S ) ;
    }

    valueType x1, y1 ;
    CD.eval( L, x1, y1 ) ;
    valueType d1 = hypot( qx-x1, qy-y1 ) ;
    if ( di <= d1 ) {
      valueType kappa1 = CD.kappa(L) ;
      valueType LS     = (abs(kappa1)+sqrt(pidk+kappa1*kappa1))/adk ;
      ClothoidData CD1 ;
      CD1.x0     = x1 ;
      CD1.y0     = y1 ;
      CD1.theta0 = theta1+m_pi ;
      CD1.k0     = -kappa1 ;
      CD1.dk     = CD.dk ;
      valueType d = closest3( epsi, CD1, LS, qx, qy, X, Y, S ) ;
      S = L-S ;
      return d ;
    }

    valueType ss = L/2, thetas(0), kappas(0), xs(0), ys(0) ;
    for ( int iter = 0 ; iter < 20 ; ++iter ) {
      CD.eval( ss, thetas, kappas, xs, ys ) ;
      valueType rhox = qx - Xi ;
      valueType rhoy = qy - Yi ;
      valueType rho  = hypot( rhox, rhoy ) ;
      valueType f    = rho - di ;
      if ( abs(f) < epsi ) break ;
      valueType phi   = atan2( rhoy, rhox ) ;
      valueType drho  = cos(thetas - phi) ;
      valueType t     = sin(thetas - phi) ;
      valueType ddrho = t*(kappas-t/rho) ;
      ss -= (f*drho)/((drho*drho)-f*ddrho/2) ;
    }

    valueType t = abs(kappas)+sqrt(pidk+kappas*kappas) ;
    ClothoidData CP, CM ;

    CP.x0 = CM.x0 = xs ;
    CP.y0 = CM.y0 = ys ;
    CP.dk = CM.dk = CD.dk ;
    CP.theta0 = thetas ;      CP.k0 =  kappas ;
    CM.theta0 = thetas+m_pi ; CM.k0 = -kappas ;

    valueType LP = min(L-ss,pi4/t);
    valueType LM = min(ss,t/adk);

    valueType dp = closest3( epsi, CP, LP, qx, qy,  X,  Y, S ) ;
    valueType dm = closest3( epsi, CM, LM, qx, qy, x1, y1, t ) ;
    if ( dp <= dm ) { S += ss ; return dp ; }
    X = x1 ;
    Y = y1 ;
    S = ss - t ;
    return dm ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::closestPoint( valueType   qx,
                               valueType   qy,
                               valueType & X,
                               valueType & Y,
                               valueType & S ) const {
    valueType epsi = 1e-10 ;
    // check if flex is inside curve, if so then split
    if ( CD.k0*CD.dk >= 0 ) { // flex on the left
      return closest2( epsi, CD, L, qx, qy, X, Y, S );
    } else if ( CD.dk*CD.kappa(L) <= 0 ) { // flex on the right, reverse curve
      ClothoidData CD1 ;
      CD.eval( L, CD1.theta0, CD1.k0, CD1.x0, CD1.y0 );
      CD1.theta0 += m_pi ;
      CD1.k0      = -CD1.k0 ;
      CD1.dk      = CD.dk ;
      valueType d = closest2( epsi, CD1, L, qx, qy, X, Y, S );
      S = L-S ;
      return d ;
    }

    // flex inside, split clothoid
    ClothoidData C0, C1 ;
    valueType sflex = -CD.k0/CD.dk ;
    valueType kappa1 ;
    eval( sflex, C0.theta0, kappa1, C0.x0, C0.y0 );
    C1.x0     = C0.x0 ;
    C1.y0     = C0.y0 ;
    C1.theta0 = C0.theta0+m_pi ;
    C0.k0     = C1.k0 = 0 ;
    C0.dk     = C1.dk = CD.dk ;

    valueType d0 = closest2( epsi, C0, L-sflex, qx, qy, X, Y, S  );
    valueType x1, y1, s1 ;
    valueType d1 = closest2( epsi, C1, sflex, qx, qy, x1, y1, s1 );

    if ( d1 < d0 ) {
      S = sflex - s1 ; X = x1 ; Y = y1 ;
      return d1 ;
    }
    S += sflex ;
    return d0 ;
  }

/*
+      // approx clothoid with a circle
       valueType theta, kappa ;
       eval( S, theta, kappa, X, Y );
-      valueType dx   = X-x ;
-      valueType dy   = Y-y ;
-      valueType d    = hypot( dx, dy );
-      valueType tphi = theta - atan2( dy, dx ) ;
-      valueType f    = d*cos(tphi) ;
-      valueType g    = d*sin(tphi) ;
-      valueType df   = 1-kappa*g ;
-      valueType ddf  = -kappa*f*(dk+kappa) ;
-
*/

}

// EOF: ClothoidDistance.cc
