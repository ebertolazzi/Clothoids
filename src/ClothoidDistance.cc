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
  valueType
  closest_by_sample( ClothoidData const & CD,
                     valueType            L,
                     valueType            qx,
                     valueType            qy,
                     valueType          & X,
                     valueType          & Y,
                     valueType          & S ) {
    S = 0 ;
    X = CD.x0 ;
    Y = CD.y0 ;
    valueType DS  = std::min(0.1,L/100) ;
    valueType DST = hypot( X-qx, Y-qy ) ;
    valueType SSS = DS ;
    while ( SSS <= L ) {
      valueType theta, kappa, XS, YS ;
      CD.eval( SSS, theta, kappa, XS, YS );
      valueType dst = hypot( XS-qx, YS-qy ) ;
      if ( dst < DST ) {
        DST = dst ;
        S   = SSS ;
        X   = XS ;
        Y   = YS ;
      }
      SSS += DS ;
    }
    return DST ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  bool
  closestPointQC2( valueType            epsi,
                   ClothoidData const & CD,
                   valueType            L,
                   valueType            qx,
                   valueType            qy,
                   valueType          & S ) {

    //valueType X, Y ;
    //closest_by_sample( CD, L, qx, qy, X, Y, S );
    //return true ;

    // S = GUESS
    int nb = 0 ;
    valueType theta, kappa, dS, XS, YS ;
    valueType SS = S ;
    for ( int iter = 0 ; iter < 20 && nb < 2 ; ++iter ) {
      // approx clothoid with a circle
      CD.eval( SS, theta, kappa, XS, YS );

      valueType CSS = cos(theta) ;
      valueType SSS = sin(theta) ;

      valueType dx  = XS - qx ;
      valueType dy  = YS - qy ;
      valueType a0  = CSS * dy - SSS * dx ;
      valueType b0  = SSS * dy + CSS * dx ;
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

      SS += dS ;
      if ( std::abs( dS ) < epsi ) {
        if ( SS < -epsi || SS > L+epsi ) return false ;
        S = SS ;
        return true ;
      }

      // check divergence
      if ( SS < 0 ) {
        SS = 0 ; ++nb ;
      } else if ( SS > L ) {
        SS = L ; ++nb ;
      } else {
        nb = 0 ;
      }

    }
    return false ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closestPointQC1( valueType            epsi,
                   ClothoidData const & CD,
                   valueType            L,
                   valueType            qx,
                   valueType            qy,
                   valueType          & X,
                   valueType          & Y,
                   valueType          & S ) {

    valueType phi0 = CD.theta0 - atan2( CD.y0 - qy, CD.x0 - qx ) ;
    bool ok0 = cos(phi0) < 0 ; // distanza decrescente

    valueType theta1, kappa1, x1, y1 ;
    CD.eval( L, theta1, kappa1, x1, y1 ) ;
    valueType phi1 = theta1 - atan2( y1 - qy, x1 - qx ) ;
    bool ok1 = cos(phi1) > 0 ; // distanza crescente

    valueType s0 = 0, x0 = CD.x0, y0 = CD.y0 ;
    if ( ok0 ) ok0 = closestPointQC2( epsi, CD, L, qx, qy, s0 ) ;
    if ( ok0 ) CD.eval( s0, x0, y0 ) ;
    valueType d0 = hypot( x0-qx, y0-qy ) ;

    valueType s1 = L ;
    if ( ok1 ) ok1 = closestPointQC2( epsi, CD, L, qx, qy, s1 ) ;
    if ( ok1 ) CD.eval( s1, x1, y1 ) ;
    valueType d1 = hypot( x1-qx, y1-qy ) ;

    if ( s1 - s0 > 2 * epsi ) { // buoni entrambi estremi
      S = (s0+s1)/2 ;
      bool okm = closestPointQC2( epsi, CD, L, qx, qy, S ) ;
      if ( okm ) {
        CD.eval( S, X, Y ) ;
        valueType dm = hypot( X-qx, Y-qy ) ;
        if ( dm < d0 && dm < d1 ) return dm ;
      }
    }

    if ( d0 < d1 ) {
      S = s0 ; X = x0 ; Y = y0 ; return d0 ;
    } else {
      S = s1 ; X = x1 ; Y = y1 ; return d1 ;
    }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closestPointQC( valueType            epsi,
                  ClothoidData const & CD,
                  valueType            L,
                  valueType            qx,
                  valueType            qy,
                  valueType          & X,
                  valueType          & Y,
                  valueType          & S ) {

    valueType DTheta = std::abs( CD.theta(L) - CD.theta0 ) ;
    if ( DTheta <= m_2pi )
      return closestPointQC1( epsi, CD, L, qx, qy, X, Y, S );

    valueType cx = CD.c0x() ;
    valueType cy = CD.c0y() ;

    if ( hypot( CD.x0 - cx, CD.y0 - cy ) <= hypot( qx - cx, qy - cy ) ) {
      valueType ell = CD.aplus(m_2pi) ;
      return closestPointQC1( epsi, CD, ell, qx, qy, X, Y, S );
    }

    ClothoidData CD1 ;
    CD.reverse( L, CD1 );
    cx = CD1.c0x() ;
    cy = CD1.c0y() ;

    if ( hypot( CD1.x0 - cx, CD1.y0 - cy ) >= hypot( qx - cx, qy - cy ) ) {
      valueType ell = CD1.aplus(m_2pi) ;
      valueType d = closestPointQC1( epsi, CD1, ell, qx, qy, X, Y, S );
      S = L - S ;
      return d ;
    }

    valueType ell = CD.aplus(DTheta/2) ;
    valueType d0  = closestPointQC( epsi, CD, ell, qx, qy, X, Y, S );

    CD.eval( ell, CD1 );

    valueType X1, Y1, S1 ;
    valueType d1 = closestPointQC( epsi, CD1, L-ell, qx, qy, X1, Y1, S1 );

    if ( d1 < d0 ) { S = ell+S1 ; X = X1 ; Y = Y1 ; return d1 ; }

    return d0 ; // closest_by_sample( CD, L, qx, qy, X, Y, S );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  bool
  closestPointStandard3( valueType   epsi,
                         valueType   a,
                         valueType   b,
                         valueType   qx,
                         valueType   qy,
                         valueType & S ) {
    // S = GUESS
    int nb = 0 ;
    valueType SS = S, dS ;
    for ( int iter = 0 ; iter < 20 && nb < 2 ; ++iter ) {
      // approx clothoid with a circle
      valueType theta = m_pi_2 * (SS*SS) ;
      valueType kappa = m_pi * SS ;
      valueType XS, YS ;
      FresnelCS( SS, XS, YS ) ;

      valueType CSS = cos(theta) ;
      valueType SSS = sin(theta) ;

      valueType dx  = XS - qx ;
      valueType dy  = YS - qy ;
      valueType a0  = CSS * dy - SSS * dx ;
      valueType b0  = SSS * dy + CSS * dx ;
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

      SS += dS ;
      if ( std::abs( dS ) < epsi ) {
        if ( SS < a-epsi|| SS > b+epsi ) break ;
        S = SS ;
        return false ; // failed = false
      }

      // check divergence
      if ( SS < a || SS > b ) ++nb ; else nb = 0 ;

    }
    return true ; // failed = true
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closestPointStandard2( valueType   epsi,
                         valueType   a,
                         valueType   b,
                         valueType   qx,
                         valueType   qy,
                         valueType & S ) {

    valueType x, y ;
    FresnelCS( a, x, y ) ;
    valueType phia = m_pi_2 * (a*a) - atan2( y - qy, x - qx ) ;
    bool ok0 = cos(phia) < 0 ; // distanza decrescente

    FresnelCS( b, x, y ) ;
    valueType phib = m_pi_2 * (b*b) - atan2( y - qy, x - qx ) ;
    bool ok1 = cos(phib) > 0 ; // distanza crescente

    valueType s0 = a ;
    if ( ok0 ) closestPointStandard3( epsi, a, b, qx, qy, s0 ) ;
    FresnelCS( s0, x, y ) ;
    valueType d0 = hypot( x-qx, y-qy ) ;

    valueType s1 = b ;
    if ( ok1 ) closestPointStandard3( epsi, a, b, qx, qy, s1 ) ;
    FresnelCS( s1, x, y ) ;
    valueType d1 = hypot( x-qx, y-qy ) ;

    if ( s1 - s0 > 2 * epsi ) { // buoni entrambi estremi
      S = (s0+s1)/2 ;
      bool failm = closestPointStandard3( epsi, a, b, qx, qy, S ) ;
      if ( !failm ) {
        FresnelCS( S, x, y ) ;
        valueType dm = hypot( x-qx, y-qy ) ;
        if ( dm < d0 && dm < d1 ) return dm ;
      }
    }

    if ( d0 < d1 ) {
      S = s0 ; return d0 ;
    } else {
      S = s1 ; return d1 ;
    }

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closestPointStandard( valueType            epsi,
                        ClothoidData const & CD,
                        valueType            L,
                        valueType            qx,
                        valueType            qy,
                        valueType          & S ) {

    // transform to standard clothoid
    valueType sflex  = -CD.kappa0/CD.dk ;

    G2LIB_ASSERT( sflex <= 0,
                  " bad sflex = " << sflex ) ;

    valueType thflex = CD.theta0 + 0.5*CD.kappa0*sflex ;
    valueType ssf    = sin(thflex) ;
    valueType csf    = cos(thflex) ;
    valueType gamma  = sqrt(std::abs(CD.dk)/m_pi) ;
    valueType a      = -sflex*gamma ;
    valueType b      = (L-sflex)*gamma ;
    valueType xflex, yflex ;
    CD.eval( sflex, xflex, yflex ) ;
    valueType xx = qx - xflex ;
    valueType yy = qy - yflex ;
    // gamma * R^(-1)
    valueType qqx = gamma * (  csf * xx + ssf * yy ) ;
    valueType qqy = gamma * ( -ssf * xx + csf * yy ) ;
    // M^(-1)
    if ( CD.dk < 0 ) qqy = -qqy ;

    // now in standard form
    if ( b*b-a*a <= 4 ) {
      valueType d = closestPointStandard2( epsi, a, b, qqx, qqy, S ) ;
      S = sflex + S/gamma ;
      return d/gamma ;
    }

    FresnelCS( a, xx, yy ) ;
    valueType di = hypot(qqx-0.5,qqy-0.5);
    valueType da = hypot(xx-0.5,yy-0.5);

    if ( di >= da ) {
      valueType La = 4/(a+sqrt(a*a+4)) ;
      valueType d = closestPointStandard2( epsi, a, a+La, qqx, qqy, S ) ;
      S = sflex + S/gamma ;
      return d/gamma ;
    }

    FresnelCS( b, xx, yy ) ;
    valueType db = hypot(xx-0.5,yy-0.5);

    if ( di <= db ) {
      valueType Lb = 4/(b+sqrt(b*b-4)) ;
      G2LIB_ASSERT( Lb <= b-a, " bad Lb = " << Lb ) ;
      valueType d = closestPointStandard2( epsi, b-Lb, b, qqx, qqy, S ) ;
      S = sflex + S/gamma ;
      return d/gamma ;
    }

    valueType ss = a;
    bool converged = false ;
    for ( int iter = 0 ; iter < 20 && !converged ; ++iter ) {
      FresnelCS( ss, xx, yy ) ;
      valueType kappa = m_pi * ss;
      valueType theta = m_pi_2 * (ss*ss) ;
      valueType rhox  = xx - 0.5 ;
      valueType rhoy  = yy - 0.5 ;
      valueType rho   = hypot( rhox, rhoy ) ;
      valueType f     = rho - di ;
      //if ( std::abs(f) < epsi ) break ;
      valueType tphi  = theta - atan2( rhoy, rhox );
      valueType df    = cos( tphi );
      valueType t     = sin( tphi );
      valueType ddf   = t*(kappa-t/rho);
      valueType ds    = (f*df)/((df*df)-f*ddf/2) ;
      ss -= ds ;
      converged = std::abs(ds) < epsi ;
    }

    G2LIB_ASSERT( converged,
                  " closestPointStandard not converged " ) ;

    valueType Lp = min( b-ss, 4/(ss+sqrt(ss*ss+4)) ) ;
    valueType Lm = min( ss-a, 4/(ss+sqrt(ss*ss-4)) ) ;

    valueType sp, sm ;
    valueType dp = closestPointStandard2( epsi, ss, ss+Lp, qqx, qqy, sp ) ;
    valueType dm = closestPointStandard2( epsi, ss-Lm, ss, qqx, qqy, sm ) ;

    if ( dp < dm ) {
      S = sflex + sp/gamma ;
      return dp/gamma ;
    } else {
      S = sflex + sm/gamma ;
      return dm/gamma ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  closestPoint1( valueType            epsi,
                 ClothoidData const & CD,
                 valueType            L,
                 valueType            qx,
                 valueType            qy,
                 valueType          & X,
                 valueType          & Y,
                 valueType          & S ) {

    valueType DK = sqrt(10*m_2pi*std::abs(CD.dk)) ;
    if ( std::abs(CD.kappa0) >= DK ) {
      return closestPointQC( epsi, CD, L, qx, qy, X, Y, S );
    }

    if ( std::abs(CD.kappa0)+std::abs(CD.dk)*L <= DK ) {
      valueType d = closestPointStandard( epsi, CD, L, qx, qy, S );
      CD.eval( S, X, Y ) ;
      return d ;
    }

    valueType ell = (DK-std::abs(CD.kappa0))/std::abs(CD.dk) ;

    G2LIB_ASSERT( ell > 0 && ell < L,
                  " bad ell = " << ell << " L = " << L ) ;

    ClothoidData CDS ;
    CD.eval( ell, CDS );

    valueType S0 ;
    valueType d0 = closestPointStandard( epsi, CD, ell, qx, qy, S0 );
    valueType d1 = closestPointQC( epsi, CDS, L-ell, qx, qy, X, Y, S );
    if ( d0 < d1 ) {
      S = S0 ;
      CD.eval( S, X, Y ) ;
      return d0 ;
    } else {
      S += ell ;
      return d1 ;
    }
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

    if ( CD.kappa0*CD.dk >= 0 ) { // flex on the left
      return closestPoint1( epsi, CD, L, qx, qy, X, Y, S );
    }

    if ( CD.dk*CD.kappa(L) <= 0 ) { // flex on the right, reverse curve
      ClothoidData CD1 ;
      CD.reverse( L, CD1 ) ;
      valueType d = closestPoint1( epsi, CD1, L, qx, qy, X, Y, S );
      S = L-S ;
      return d ;
    }

    // flex inside, split clothoid
    ClothoidData C0, C1 ;
    valueType sflex = CD.split_at_flex( C0, C1 ) ;

    valueType d0 = closestPoint1( epsi, C0, L-sflex, qx, qy, X, Y, S  );
    valueType x1, y1, s1 ;
    valueType d1 = closestPoint1( epsi, C1, sflex, qx, qy, x1, y1, s1 );

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
