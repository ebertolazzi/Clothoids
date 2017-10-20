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
#include "CubicRootsFlocke.hh"

#include <cmath>
#include <cfloat>

namespace Clothoid {

  using namespace std ;

  static const valueType m_pi   = 3.14159265358979323846264338328  ; // pi
  static const valueType m_pi_2 = 1.57079632679489661923132169164  ; // pi/2
  static const valueType m_2pi  = 6.28318530717958647692528676656  ; // 2*pi
  static const valueType m_1_pi = 0.318309886183790671537767526745 ; // 1/pi

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

  // ---------------------------------------------------------------------------
  int
  ClothoidCurve::build( valueType x1,
                        valueType y1,
                        valueType theta1 ) {

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
    valueType L = r/intC[0] ;

    G2LIB_ASSERT( L > 0, "Negative length L = " << L ) ;
    k  = (delta-A)/L ;
    dk = 2*A/L/L ;

    s_min = 0 ;
    s_max = L ;

    return niter ;
  }

  valueType
  ClothoidCurve::X( valueType s ) const {
    valueType C, S ;
    GeneralizedFresnelCS( dk*s*s, k*s, theta0, C, S ) ;
    return x0 + s*C ;
  }
  
  valueType
  ClothoidCurve::Y( valueType s ) const {
    valueType C, S ;
    GeneralizedFresnelCS( dk*s*s, k*s, theta0, C, S ) ;
    return y0 + s*S ;
  }

  void
  ClothoidCurve::eval( valueType   s,
                       valueType & theta,
                       valueType & kappa,
                       valueType & x,
                       valueType & y ) const {
    valueType C, S ;
    GeneralizedFresnelCS( dk*s*s, k*s, theta0, C, S ) ;
    x     = x0 + s*C ;
    y     = y0 + s*S ;
    theta = theta0 + s*(k+s*(dk/2)) ;
    kappa = k + s*dk ;
  }

  void
  ClothoidCurve::eval( valueType s, valueType & x, valueType & y ) const {
    valueType C, S ;
    GeneralizedFresnelCS( dk*s*s, k*s, theta0, C, S ) ;
    x = x0 + s*C ;
    y = y0 + s*S ;
  }

  void
  ClothoidCurve::eval_D( valueType s, valueType & x_D, valueType & y_D ) const {
    valueType theta = theta0 + s*(k+s*(dk/2)) ;
    x_D = cos(theta) ;
    y_D = sin(theta) ;
  }

  void
  ClothoidCurve::eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const {
    valueType theta   = theta0 + s*(k+s*(dk/2)) ;
    valueType theta_D = k+s*dk ;
    x_DD = -sin(theta)*theta_D ;
    y_DD =  cos(theta)*theta_D  ;
  }

  void
  ClothoidCurve::eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const {
    valueType theta   = theta0 + s*(k+s*(dk/2)) ;
    valueType theta_D = k+s*dk ;
    valueType C       = cos(theta) ;
    valueType S       = sin(theta) ;
    valueType th2     = theta_D*theta_D ;
    x_DDD = -C*th2-S*dk ;
    y_DDD = -S*th2+C*dk  ;
  }

  // offset curve
  void
  ClothoidCurve::eval( valueType s, valueType offs, valueType & x, valueType & y ) const {
    valueType C, S ;
    GeneralizedFresnelCS( dk*s*s, k*s, theta0, C, S ) ;
    valueType theta = theta0 + s*(k+s*(dk/2)) ;
    valueType nx    = -sin(theta) ;
    valueType ny    =  cos(theta) ;
    x = x0 + s*C + offs * nx ;
    y = y0 + s*S + offs * ny ;
  }

  void
  ClothoidCurve::eval_D( valueType s, valueType offs, valueType & x_D, valueType & y_D ) const {
    valueType theta   = theta0 + s*(k+s*(dk/2)) ;
    valueType theta_D = k+s*dk ;
    valueType scale   = 1-offs*theta_D ;
    x_D = cos(theta)*scale ;
    y_D = sin(theta)*scale ;
  }

  void
  ClothoidCurve::eval_DD( valueType s, valueType offs, valueType & x_DD, valueType & y_DD ) const {
    valueType theta   = theta0 + s*(k+s*(dk/2)) ;
    valueType theta_D = k+s*dk ;
    valueType C       = cos(theta) ;
    valueType S       = sin(theta) ;
    valueType tmp1    = theta_D*(1-theta_D*offs) ;
    valueType tmp2    = offs*dk ;
    x_DD = -tmp1*S - C*tmp2 ;
    y_DD =  tmp1*C - S*tmp2 ;
  }

  void
  ClothoidCurve::eval_DDD( valueType s, valueType offs, valueType & x_DDD, valueType & y_DDD ) const {
    valueType theta   = theta0 + s*(k+s*(dk/2)) ;
    valueType theta_D = k+s*dk ;
    valueType C       = cos(theta) ;
    valueType S       = sin(theta) ;
    valueType tmp1    = theta_D*theta_D*(theta_D*offs-1) ;
    valueType tmp2    = dk*(1-3*theta_D*offs) ;
    x_DDD = tmp1*C-tmp2*S ;
    y_DDD = tmp1*S+tmp2*C ;
  }

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

  bool
  ClothoidCurve::build_forward( valueType _x0,
                                valueType _y0,
                                valueType _theta0,
                                valueType _k,
                                valueType _x1,
                                valueType _y1,
                                valueType tol ) {

    x0     = _x0 ;
    y0     = _y0 ;
    theta0 = _theta0 ;
    k      = _k ;
    s_min  = 0 ;

    // Compute guess angles
    valueType len  = hypot( _y1-_y0, _x1-_x0 ) ;
    valueType arot = atan2( _y1-_y0, _x1-_x0 ) ;
    valueType th0  = theta0 - arot ;
    // normalize angle
    while ( th0 >  m_pi ) th0 -= m_2pi ;
    while ( th0 < -m_pi ) th0 += m_2pi ;

    // solve the problem from (0,0) to (1,0)
    valueType k0    = k*len ;
    valueType alpha = 2.6 ;
    valueType thmin = max(-m_pi,-theta0/2-alpha) ;
    valueType thmax = min( m_pi,-theta0/2+alpha) ;
    valueType Kmin  = kappa( th0, thmax ) ;
    valueType Kmax  = kappa( th0, thmin ) ;
    bool ok ;
    valueType th = theta_guess( th0, max(min(k0,Kmax),Kmin), ok ) ;
    if ( ok ) {
      for ( indexType iter = 0 ; iter < 10 ; ++iter ) {
        valueType dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ;
        buildClothoid( 0, 0, th0,
                       1, 0, th,
                       k, dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ) ;
        valueType f   = k - k0 ;
        valueType df  = k_2 ;
        valueType dth = f/df ;
        th -= dth ;
        if ( abs(dth) < tol && abs(f) < tol ) {
          // transform solution
          buildClothoid( x0, y0, theta0,
                         _x1, _y1, arot + th,
                         _k, dk, s_max ) ;
          return true ;
        }
      }
    }
    return false ;
  }

  void
  ClothoidCurve::change_origin( valueType s0 ) {
    valueType new_theta, new_kappa, new_x0, new_y0 ;
    eval( s0, new_theta, new_kappa, new_x0, new_y0 ) ;
    x0     = new_x0 ;
    y0     = new_y0 ;
    theta0 = new_theta ;
    k      = new_kappa ;
    s_min -= s0 ;
    s_max -= s0 ;
  }

  bool
  ClothoidCurve::bbTriangle( valueType offs,
                             valueType p0[2],
                             valueType p1[2],
                             valueType p2[2] ) const {
    valueType theta_max = theta( s_max ) ;
    valueType theta_min = theta( s_min ) ;
    valueType dtheta    = std::abs( theta_max-theta_min ) ;
    if ( dtheta < m_pi_2 ) {
      valueType alpha, t0[2] ;
      eval( s_min, offs, p0[0], p0[1] ) ;
      eval_D( s_min, t0[0], t0[1] ) ; // no offset
      if ( dtheta > 0.0001 * m_pi_2 ) {
        valueType t1[2] ;
        eval( s_max, offs, p1[0], p1[1] ) ;
        eval_D( s_max, t1[0], t1[1] ) ; // no offset
        // risolvo il sistema
        // p0 + alpha * t0 = p1 + beta * t1
        // alpha * t0 - beta * t1 = p1 - p0
        valueType det = t1[0]*t0[1]-t0[0]*t1[1] ;
        alpha = ((p1[1]-p0[1])*t1[0] - (p1[0]-p0[0])*t1[1])/det ;
      } else {
        // se angolo troppo piccolo uso approx piu rozza
        alpha = s_max - s_min ;
      }
      p2[0] = p0[0] + alpha*t0[0] ;
      p2[1] = p0[1] + alpha*t0[1] ;
      return true ;
    } else {
      return false ;
    }
  }

  void
  ClothoidCurve::bbSplit(
    valueType               split_angle,
    valueType               split_size,
    valueType               split_offs,
    vector<ClothoidCurve> & c,
    vector<T2D>           & t
  ) const {

    // step 0: controllo se curvatura passa per 0
    valueType k_min = theta_D( s_min ) ;
    valueType k_max = theta_D( s_max ) ;
    c.clear() ;
    t.clear() ;
    if ( k_min * k_max < 0 ) {
      // risolvo (s-s_min)*dk+k_min = 0 --> s = s_min-k_min/dk
      valueType s_med = s_min-k_min/dk ;
      ClothoidCurve tmp(*this) ;
      tmp.trim(s_min,s_med) ;
      tmp.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
      tmp.trim(s_med,s_max) ;
      tmp.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    } else {
      bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    }
  }

  static
  valueType
  abs2pi( valueType a ) {
    a = std::abs(a) ;
    while ( a > m_pi ) a -= m_2pi ;
    return std::abs(a) ;
  }

  void
  ClothoidCurve::bbSplit_internal(
    valueType               split_angle,
    valueType               split_size,
    valueType               split_offs,
    vector<ClothoidCurve> & c,
    vector<T2D>           & t
  ) const {

    valueType theta_min, kappa_min, x_min, y_min,
              theta_max, kappa_max, x_max, y_max ;

    eval( s_min, theta_min, kappa_min, x_min, y_min ) ;
    eval( s_max, theta_max, kappa_max, x_max, y_max ) ;

    valueType dtheta = std::abs( theta_max - theta_min ) ;
    valueType dx     = x_max - x_min ;
    valueType dy     = y_max - y_min ;
    valueType len    = hypot( dy, dx ) ;
    valueType dangle = abs2pi(atan2( dy, dx )-theta_min) ;
    if ( dtheta <= split_angle && len*tan(dangle) <= split_size ) {
      T2D tt ;
      this->bbTriangle(split_offs,tt) ;
      c.push_back(*this) ;
      t.push_back(tt) ;
    } else {
      ClothoidCurve cc(*this) ;
      valueType s_med = (s_min+s_max)/2 ;
      cc.trim(s_min,s_med) ;
      cc.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
      cc.trim(s_med,s_max) ;
      cc.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    }
  }

  bool
  ClothoidCurve::intersect_internal( ClothoidCurve & c1,
                                     valueType       c1_offs,
                                     valueType     & s1,
                                     ClothoidCurve & c2,
                                     valueType       c2_offs,
                                     valueType     & s2,
                                     indexType       max_iter,
                                     valueType       tolerance ) const {
    valueType angle1a = c1.theta(c1.s_min) ;
    valueType angle1b = c1.theta(c1.s_max) ;
    valueType angle2a = c2.theta(c2.s_min) ;
    valueType angle2b = c2.theta(c2.s_max) ;
    // cerca angoli migliori per partire
    valueType dmax = abs2pi(angle1a-angle2a) ;
    valueType dab  = abs2pi(angle1a-angle2b) ;
    valueType dba  = abs2pi(angle1b-angle2a) ;
    valueType dbb  = abs2pi(angle1b-angle2b) ;
    s1 = c1.s_min ; s2 = c2.s_min ;
    if ( dmax < dab ) { dmax = dab ; s2 = c2.s_max ; }
    if ( dmax < dba ) { dmax = dba ; s1 = c1.s_min ; s2 = c2.s_min ; }
    if ( dmax < dbb ) {              s1 = c1.s_min ; s2 = c2.s_max ; }
    for ( indexType i = 0 ; i < max_iter ; ++i ) {
      valueType t1[2], t2[2], p1[2], p2[2] ;
      c1.eval( s1, c1_offs, p1[0], p1[1] ) ;
      c1.eval_D( s1, c1_offs, t1[0], t1[1] ) ;
      c2.eval( s2, c2_offs, p2[0], p2[1] ) ;
      c2.eval_D( s2, c2_offs, t2[0], t2[1] ) ;
      /*
      // risolvo il sistema
      // p1 + alpha * t1 = p2 + beta * t2
      // alpha * t1 - beta * t2 = p2 - p1
      //
      //  / t1[0] -t2[0] \ / alpha \ = / p2[0] - p1[0] \
      //  \ t1[1] -t2[1] / \ beta  /   \ p2[1] - p1[1] /
      */
      valueType det = t2[0]*t1[1]-t1[0]*t2[1] ;
      valueType px  = p2[0]-p1[0] ;
      valueType py  = p2[1]-p1[1] ;
      s1 += (py*t2[0] - px*t2[1])/det ;
      s2 += (t1[0]*py - t1[1]*px)/det ;
      if ( s1 <= c1.s_min || s1 >= c1.s_max ||
           s2 <= c2.s_min || s2 >= c2.s_max ) break ;
      if ( std::abs(px) <= tolerance ||
           std::abs(py) <= tolerance ) return true ;
    }
    return false ;
  }

  void
  ClothoidCurve::intersect( valueType             offs,
                            ClothoidCurve const & clot,
                            valueType             clot_offs,
                            vector<valueType>   & s1,
                            vector<valueType>   & s2,
                            indexType             max_iter,
                            valueType             tolerance ) const {
    vector<ClothoidCurve> c0, c1 ;
    vector<T2D>           t0, t1 ;
    bbSplit( m_pi/50, (s_max-s_min)/3, offs, c0, t0 ) ;
    clot.bbSplit( m_pi/50, (clot.s_max-clot.s_min)/3, clot_offs, c1, t1 ) ;
    s1.clear() ;
    s2.clear() ;
    for ( indexType i = 0 ; i < indexType(c0.size()) ; ++i ) {
      for ( indexType j = 0 ; j < indexType(c1.size()) ; ++j ) {
        if ( t0[i].overlap(t1[j]) ) {
          // uso newton per cercare intersezione
          valueType tmp_s1, tmp_s2 ;
          bool ok = intersect_internal( c0[i], offs,      tmp_s1,
                                        c1[j], clot_offs, tmp_s2,
                                        max_iter, tolerance ) ;
          if ( ok ) {
            s1.push_back(tmp_s1) ;
            s2.push_back(tmp_s2) ;
          }
        }
      }
    }
  }
  
  // collision detection
  bool
  ClothoidCurve::approsimate_collision( valueType             offs,
                                        ClothoidCurve const & clot,
                                        valueType             clot_offs,
                                        valueType             max_angle,
                                        valueType             max_size ) const {
    vector<ClothoidCurve> c0, c1 ;
    vector<T2D>           t0, t1 ;
    bbSplit( max_angle, max_size, offs, c0, t0 ) ;
    clot.bbSplit( max_angle, max_size, clot_offs, c1, t1 ) ;
    for ( indexType i = 0 ; i < indexType(c0.size()) ; ++i ) {
      for ( indexType j = 0 ; j < indexType(c1.size()) ; ++j ) {
        if ( t0[i].overlap(t1[j]) ) return true ;
      }
    }
    return false ;
  }

  void
  ClothoidCurve::rotate( valueType angle, valueType cx, valueType cy ) {
    valueType dx  = x0 - cx ;
    valueType dy  = y0 - cy ;
    valueType C   = cos(angle) ;
    valueType S   = sin(angle) ;
    valueType ndx = C*dx - S*dy ;
    valueType ndy = C*dy + S*dx ;
    x0      = cx + ndx ;
    y0      = cy + ndy ;
    theta0 += angle ;
  }

  void
  ClothoidCurve::scale( valueType s ) {
    k     /= s ;
    dk    /= s*s ;
    s_min *= s ;
    s_max *= s ;
  }

  void
  ClothoidCurve::reverse() {
    theta0 = theta0 + m_pi ;
    if ( theta0 > m_pi ) theta0 -= 2*m_pi ;
    k     = -k ;
    valueType tmp = s_max ;
    s_max = -s_min ;
    s_min = -tmp ;
  }

  valueType
  ClothoidCurve::thetaTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType kL  = k+dk*s_min ;
    valueType kR  = k+dk*s_max ;
    valueType thL = s_min*(k+dk*s_min/2) ;
    valueType thR = s_max*(k+dk*s_max/2) ;
    if ( kL*kR < 0 ) {
      valueType root = -k/dk ;
      if ( root > s_min && root < s_max ) {
        valueType thM  = root*(k+dk*root/2) ;
        return std::abs( thR - thM ) + std::abs( thM - thL ) ;
      }
    }
    return std::abs( thR - thL ) ;
  }

  valueType
  ClothoidCurve::thetaMinMax( valueType & thMin, valueType & thMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType kL  = k+dk*s_min ;
    valueType kR  = k+dk*s_max ;
    valueType thL = s_min*(k+dk*s_min/2) ;
    valueType thR = s_max*(k+dk*s_max/2) ;
    if ( thL < thR ) { thMin = thL ; thMax = thR ; }
    else             { thMin = thR ; thMax = thL ; }
    if ( kL*kR < 0 ) {
      valueType root = -k/dk ;
      if ( root > s_min && root < s_max ) {
        valueType thM = root*(k+dk*root/2) ;
        if      ( thM < thMin ) thMin = thM ;
        else if ( thM > thMax ) thMax = thM ;
      }
    }
    return thMax - thMin ;
  }

  valueType
  ClothoidCurve::curvatureMinMax( valueType & kMin, valueType & kMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    kMin = k+dk*s_min ;
    kMax = k+dk*s_max ;
    if ( kMax < kMin ) std::swap( kMax, kMin ) ;
    return kMax - kMin ;
  }

  valueType
  ClothoidCurve::curvatureTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType km = k+s_min*dk ;
    valueType kp = k+s_max*dk ;
    return std::abs(kp-km) ;
  }

  valueType
  ClothoidCurve::integralCurvature2() const {
    return (s_max-s_min)*( k*(k+(s_max+s_min)*dk) +
                          (s_max*s_max+s_max*s_min+s_min*s_min)*dk*dk/3 ) ;
  }

  valueType
  ClothoidCurve::integralJerk2() const {
    valueType s_min2 = s_min*s_min ;
    valueType s_min3 = s_min*s_min2 ;
    valueType s_min4 = s_min2*s_min2 ;
    valueType k2     = k*k ;
    valueType k3     = k*k2 ;
    valueType k4     = k2*k2 ;
    valueType t1     = s_max+s_min ;
    valueType t2     = s_max*t1+s_min2 ;
    valueType t3     = s_max*t2+s_min3 ;
    valueType t4     = s_max*t3+s_min4 ;
    return ((((t4/5*dk+t3*k)*dk+(1+2*t2)*k2)*dk+2*t1*k3)*dk+k4)*(s_max-s_min) ;
  }

  valueType
  ClothoidCurve::integralSnap2() const {
    valueType s_min2 = s_min*s_min  ;
    valueType s_min3 = s_min*s_min2 ;
    valueType s_min4 = s_min3*s_min ;
    valueType s_min5 = s_min4*s_min ;
    valueType s_min6 = s_min5*s_min ;
    valueType k2     = k*k ;
    valueType k3     = k*k2 ;
    valueType k4     = k2*k2 ;
    valueType k5     = k4*k ;
    valueType k6     = k4*k2 ;
    valueType dk2    = dk*dk ;
    valueType dk3    = dk*dk2 ;
    valueType dk4    = dk2*dk2 ;
    valueType dk5    = dk4*dk ;
    valueType dk6    = dk4*dk2 ;
    valueType t2     = s_max+s_min ;
    valueType t3     = s_max*t2+s_min2 ;
    valueType t4     = s_max*t3+s_min3 ;
    valueType t5     = s_max*t4+s_min4 ;
    valueType t6     = s_max*t5+s_min5 ;
    valueType t7     = s_max*t6+s_min6 ;

    return ( (t7/7)*dk6 + dk5*k*t6 + 3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 + 3*dk*k5*t2 + 9*dk2*k*t2 +
             k6+9*k2*dk ) * ( s_max - s_min ) ;
  }

  std::ostream &
  operator << ( std::ostream & stream, ClothoidCurve const & c ) {
    stream <<   "x0     = " << c.x0
           << "\ny0     = " << c.y0
           << "\ntheta0 = " << c.theta0
           << "\nk      = " << c.k
           << "\ndk     = " << c.dk
           << "\nL      = " << c.s_max-c.s_min
           << "\ns_min  = " << c.s_min
           << "\ns_max  = " << c.s_max
           << "\n" ;
    return stream ;
  }

}

// EOF: ClothoidG1.cc
