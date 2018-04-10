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

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------

  static valueType const CF[] = { 2.989696028701907,  0.716228953608281,
                                 -0.458969738821509, -0.502821153340377,
                                  0.261062141752652, -0.045854475238709 } ;
  int
  buildClothoid( valueType   x0,
                 valueType   y0,
                 valueType   theta0,
                 valueType   x1,
                 valueType   y1,
                 valueType   theta1,
                 valueType & k,
                 valueType & dk,
                 valueType & L ) {

    // traslazione in (0,0)
    valueType dx  = x1 - x0 ;
    valueType dy  = y1 - y0 ;
    valueType r   = hypot( dx, dy ) ;
    valueType phi = atan2( dy, dx ) ;

    valueType phi0 = theta0 - phi ;
    valueType phi1 = theta1 - phi ;
    
    phi0 -= m_2pi*round(phi0/m_2pi) ;
    phi1 -= m_2pi*round(phi1/m_2pi) ;

    if ( phi0 >  m_pi ) phi0 -= m_2pi ;
    if ( phi0 < -m_pi ) phi0 += m_2pi ;
    if ( phi1 >  m_pi ) phi1 -= m_2pi ;
    if ( phi1 < -m_pi ) phi1 += m_2pi ;

    valueType delta = phi1 - phi0 ;

    // punto iniziale
    valueType X  = phi0*m_1_pi ;
    valueType Y  = phi1*m_1_pi ;
    valueType xy = X*Y ;
    Y *= Y ; X *= X ;
    valueType A  = (phi0+phi1)*(CF[0]+xy*(CF[1]+xy*CF[2])+(CF[3]+xy*CF[4])*(X+Y)+CF[5]*(X*X+Y*Y)) ;

    // newton
    valueType g=0, dg, intC[3], intS[3] ;
    indexType niter = 0 ;
    do {
      GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS ) ;
      g   = intS[0] ;
      dg  = intC[2] - intC[1] ;
      A  -= g / dg ;
    } while ( ++niter <= 10 && std::abs(g) > 1E-12 ) ;

    G2LIB_ASSERT( std::abs(g) < 1E-8, "Newton do not converge, g = " << g << " niter = " << niter ) ;
    GeneralizedFresnelCS( 2*A, delta-A, phi0, intC[0], intS[0] ) ;
    L = r/intC[0] ;

    G2LIB_ASSERT( L > 0, "Negative length L = " << L ) ;
    k  = (delta-A)/L ;
    dk = 2*A/L/L ;
    
    return niter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int
  buildClothoid( valueType   x0,
                 valueType   y0,
                 valueType   theta0,
                 valueType   x1,
                 valueType   y1,
                 valueType   theta1,
                 valueType & k,
                 valueType & dk,
                 valueType & L,
                 valueType & k_1,
                 valueType & dk_1,
                 valueType & L_1,
                 valueType & k_2,
                 valueType & dk_2,
                 valueType & L_2 ) {

    // traslazione in (0,0)
    valueType dx  = x1 - x0 ;
    valueType dy  = y1 - y0 ;
    valueType r   = hypot( dx, dy ) ;
    valueType phi = atan2( dy, dx ) ;

    valueType phi0 = theta0 - phi ;
    valueType phi1 = theta1 - phi ;
    
    phi0 -= m_2pi*round(phi0/m_2pi) ;
    phi1 -= m_2pi*round(phi1/m_2pi) ;

    if ( phi0 >  m_pi ) phi0 -= m_2pi ;
    if ( phi0 < -m_pi ) phi0 += m_2pi ;
    if ( phi1 >  m_pi ) phi1 -= m_2pi ;
    if ( phi1 < -m_pi ) phi1 += m_2pi ;

    valueType delta = phi1 - phi0 ;

    // punto iniziale
    valueType X  = phi0*m_1_pi ;
    valueType Y  = phi1*m_1_pi ;
    valueType xy = X*Y ;
    Y *= Y ; X *= X ;
    valueType A = (phi0+phi1)*(CF[0]+xy*(CF[1]+xy*CF[2])+(CF[3]+xy*CF[4])*(X+Y)+CF[5]*(X*X+Y*Y)) ;

    // newton
    valueType g=0, dg, intC[3], intS[3] ;
    indexType niter = 0 ;
    do {
      GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS ) ;
      g   = intS[0] ;
      dg  = intC[2] - intC[1] ;
      A  -= g / dg ;
    } while ( ++niter <= 10 && std::abs(g) > 1E-12 ) ;

    G2LIB_ASSERT( std::abs(g) < 1E-8, "Newton do not converge, g = " << g << " niter = " << niter ) ;
    GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS ) ;
    L = r/intC[0] ;

    G2LIB_ASSERT( L > 0, "Negative length L = " << L ) ;
    k  = (delta-A)/L ;
    dk = 2*A/L/L ;

    valueType alpha = intC[0]*intC[1] + intS[0]*intS[1] ;
    valueType beta  = intC[0]*intC[2] + intS[0]*intS[2] ;
    valueType gamma = intC[0]*intC[0] + intS[0]*intS[0] ;
    valueType tx    = intC[1]-intC[2] ;
    valueType ty    = intS[1]-intS[2] ;
    valueType txy   = L*(intC[1]*intS[2]-intC[2]*intS[1]) ;
    valueType omega = L*(intS[0]*tx-intC[0]*ty) - txy ;

    delta = intC[0]*tx + intS[0]*ty ;

    L_1  = omega/delta ;
    L_2  = txy/delta ;

    delta *= L ;
    k_1  = (beta-gamma-k*omega)/delta ;
    k_2  = -(beta+k*txy)/delta ;

    delta *= L/2 ;
    dk_1 = (gamma-alpha-dk*omega*L)/delta ;
    dk_2 = (alpha-dk*txy*L)/delta ;

    return niter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int
  ClothoidCurve::build( valueType x1,
                        valueType y1,
                        valueType theta1 ) {

    // traslazione in (0,0)
    valueType dx  = x1 - CD.x0 ;
    valueType dy  = y1 - CD.y0 ;
    valueType r   = hypot( dx, dy ) ;
    valueType phi = atan2( dy, dx ) ;

    valueType phi0 = CD.theta0 - phi ;
    valueType phi1 = theta1 - phi ;
    
    phi0 -= m_2pi*round(phi0/m_2pi) ;
    phi1 -= m_2pi*round(phi1/m_2pi) ;

    if ( phi0 >  m_pi ) phi0 -= m_2pi ;
    if ( phi0 < -m_pi ) phi0 += m_2pi ;
    if ( phi1 >  m_pi ) phi1 -= m_2pi ;
    if ( phi1 < -m_pi ) phi1 += m_2pi ;

    valueType delta = phi1 - phi0 ;

    // punto iniziale
    valueType X  = phi0*m_1_pi ;
    valueType Y  = phi1*m_1_pi ;
    valueType xy = X*Y ;
    Y *= Y ; X *= X ;
    valueType A  = (phi0+phi1)*(CF[0]+xy*(CF[1]+xy*CF[2])+(CF[3]+xy*CF[4])*(X+Y)+CF[5]*(X*X+Y*Y)) ;

    // newton
    valueType g=0, dg, intC[3], intS[3] ;
    indexType niter = 0 ;
    do {
      GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS ) ;
      g   = intS[0] ;
      dg  = intC[2] - intC[1] ;
      A  -= g / dg ;
    } while ( ++niter <= 10 && std::abs(g) > 1E-12 ) ;

    G2LIB_ASSERT( std::abs(g) < 1E-8, "Newton do not converge, g = " << g << " niter = " << niter ) ;
    GeneralizedFresnelCS( 2*A, delta-A, phi0, intC[0], intS[0] ) ;
    L = r/intC[0] ;

    G2LIB_ASSERT( L > 0, "Negative length L = " << L ) ;
    CD.kappa0 = (delta-A)/L ;
    CD.dk     = 2*A/L/L ;

    return niter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  kappa( valueType theta0, valueType theta ) {
    valueType x = theta0*theta0 ;
    valueType a = -3.714 + x * 0.178 ;
    valueType b = -1.913 - x * 0.0753 ;
    valueType c =  0.999 + x * 0.03475 ;
    valueType d =  0.191 - x * 0.00703 ;
    valueType e =  0.500 - x * -0.00172 ;
    valueType t = d*theta0+e*theta ;
    return a*theta0+b*theta+c*(t*t*t) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  theta_guess( valueType theta0, valueType k0, bool & ok ) {
    valueType x   = theta0*theta0 ;
    valueType a   = -3.714 + x * 0.178 ;
    valueType b   = -1.913 - x * 0.0753 ;
    valueType c   =  0.999 + x * 0.03475 ;
    valueType d   =  0.191 - x * 0.00703 ;
    valueType e   =  0.500 - x * -0.00172 ;
    valueType e2  = e*e ;
    valueType dt  = d*theta0 ;
    valueType dt2 = dt*dt ;
    valueType A   = c*e*e2 ;
    valueType B   = 3*(c*d*e2*theta0) ;
    valueType C   = 3*c*e*dt2 + b ;
    valueType D   = c*(dt*dt2) + a*theta0 - k0 ;

    valueType r[3] ;
    indexType nr, nc ;
    PolynomialRoots::solveCubic( A, B, C, D, r[0], r[1], r[2], nr, nc ) ;
    // cerco radice reale piu vicina
    valueType theta ;
    switch ( nr ) {
    case 0:
    default:
      ok = false ;
      return 0 ;
    case 1:
      theta = r[0] ;
      break ;
    case 2:
      if ( abs(r[0]-theta0) < abs(r[1]-theta0) ) theta = r[0] ;
      else                                       theta = r[1] ;
      break ;
    case 3:
      theta = r[0] ;
      for ( indexType i = 1 ; i < 3 ; ++i ) {
        if ( abs(theta-theta0) > abs(r[i]-theta0) )
          theta = r[i] ;
      }
      break ;
    }
    ok = abs(theta-theta0) < m_pi ;
    return theta ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::build_forward( valueType _x0,
                                valueType _y0,
                                valueType _theta0,
                                valueType _k,
                                valueType _x1,
                                valueType _y1,
                                valueType tol ) {

    CD.x0     = _x0 ;
    CD.y0     = _y0 ;
    CD.theta0 = _theta0 ;
    CD.kappa0 = _k ;

    // Compute guess angles
    valueType len  = hypot( _y1-_y0, _x1-_x0 ) ;
    valueType arot = atan2( _y1-_y0, _x1-_x0 ) ;
    valueType th0  = CD.theta0 - arot ;
    // normalize angle
    while ( th0 >  m_pi ) th0 -= m_2pi ;
    while ( th0 < -m_pi ) th0 += m_2pi ;

    // solve the problem from (0,0) to (1,0)
    valueType k0    = CD.kappa0*len ;
    valueType alpha = 2.6 ;
    valueType thmin = max(-m_pi,-CD.theta0/2-alpha) ;
    valueType thmax = min( m_pi,-CD.theta0/2+alpha) ;
    valueType Kmin  = kappa( th0, thmax ) ;
    valueType Kmax  = kappa( th0, thmin ) ;
    bool ok ;
    valueType th = theta_guess( th0, max(min(k0,Kmax),Kmin), ok ) ;
    if ( ok ) {
      for ( indexType iter = 0 ; iter < 10 ; ++iter ) {
        valueType dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ;
        buildClothoid( 0, 0, th0,
                       1, 0, th,
                       CD.kappa0, CD.dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ) ;
        valueType f   = CD.kappa0 - k0 ;
        valueType df  = k_2 ;
        valueType dth = f/df ;
        th -= dth ;
        if ( abs(dth) < tol && abs(f) < tol ) {
          // transform solution
          buildClothoid( CD.x0, CD.y0, CD.theta0,
                         _x1, _y1, arot + th,
                         _k, dk, L ) ;
          return true ;
        }
      }
    }
    return false ;
  }

}

// EOF: ClothoidG1.cc
