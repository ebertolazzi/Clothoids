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

  void
  ClothoidCurve::changeCurvilinearOrigin( valueType s0, valueType newL ) {
    valueType new_theta, new_kappa, new_x0, new_y0 ;
    eval( s0, new_theta, new_kappa, new_x0, new_y0 ) ;
    CD.x0     = new_x0 ;
    CD.y0     = new_y0 ;
    CD.theta0 = new_theta ;
    CD.k0     = new_kappa ;
    L         = newL ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::rotate( valueType angle, valueType cx, valueType cy ) {
    valueType dx  = CD.x0 - cx ;
    valueType dy  = CD.y0 - cy ;
    valueType C   = cos(angle) ;
    valueType S   = sin(angle) ;
    valueType ndx = C*dx - S*dy ;
    valueType ndy = C*dy + S*dx ;
    CD.x0      = cx + ndx ;
    CD.y0      = cy + ndy ;
    CD.theta0 += angle ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::scale( valueType s ) {
    CD.k0 /= s ;
    CD.dk /= s*s ;
    L     *= s ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::reverse() {
    CD.theta0 += m_pi ;
    if ( CD.theta0 > m_pi ) CD.theta0 -= 2*m_pi ;
    CD.k0 = -CD.k0 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::thetaTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType kL  = CD.k0 ;
    valueType kR  = CD.k0+CD.dk*L ;
    valueType thL = 0 ;
    valueType thR = L*(CD.k0+0.5*CD.dk*L) ;
    if ( kL*kR < 0 ) {
      valueType root = -CD.k0/CD.dk ;
      if ( root > 0 && root < L ) {
        valueType thM  = root*(CD.k0+0.5*CD.dk*root) ;
        return std::abs( thR - thM ) + std::abs( thM - thL ) ;
      }
    }
    return std::abs( thR - thL ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::thetaMinMax( valueType & thMin, valueType & thMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType kL  = CD.k0 ;
    valueType kR  = CD.k0+CD.dk*L ;
    valueType thL = 0 ;
    valueType thR = L*(CD.k0+0.5*CD.dk*L) ;
    if ( thL < thR ) { thMin = thL ; thMax = thR ; }
    else             { thMin = thR ; thMax = thL ; }
    if ( kL*kR < 0 ) {
      valueType root = -CD.k0/CD.dk ;
      if ( root > 0 && root < L ) {
        valueType thM = CD.deltaTheta(root) ;
        if      ( thM < thMin ) thMin = thM ;
        else if ( thM > thMax ) thMax = thM ;
      }
    }
    return thMax - thMin ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::curvatureMinMax( valueType & kMin, valueType & kMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    kMin = CD.k0 ;
    kMax = CD.kappa(L);
    if ( kMax < kMin ) std::swap( kMax, kMin ) ;
    return kMax - kMin ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::curvatureTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType km = CD.k0 ;
    valueType kp = CD.kappa(L) ;
    return std::abs(kp-km) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralCurvature2() const {
    return L*( CD.k0*(CD.k0+L*CD.dk) + (L*L)*CD.dk*CD.dk/3 ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralJerk2() const {
    valueType k2 = CD.k0*CD.k0 ;
    valueType k3 = CD.k0*k2 ;
    valueType k4 = k2*k2 ;
    valueType t1 = L ;
    valueType t2 = L*t1 ;
    valueType t3 = L*t2 ;
    valueType t4 = L*t3 ;
    return ((((t4/5*CD.dk+t3*CD.k0)*CD.dk+(1+2*t2)*k2)*CD.dk+2*t1*k3)*CD.dk+k4)*L ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralSnap2() const {
    valueType k2  = CD.k0*CD.k0 ;
    valueType k3  = CD.k0*k2 ;
    valueType k4  = k2*k2 ;
    valueType k5  = k4*CD.k0 ;
    valueType k6  = k4*k2 ;
    valueType dk2 = CD.dk*CD.dk ;
    valueType dk3 = CD.dk*dk2 ;
    valueType dk4 = dk2*dk2 ;
    valueType dk5 = dk4*CD.dk ;
    valueType dk6 = dk4*dk2 ;
    valueType t2  = L ;
    valueType t3  = L*t2 ;
    valueType t4  = L*t3 ;
    valueType t5  = L*t4 ;
    valueType t6  = L*t5 ;
    valueType t7  = L*t6 ;

    return ( (t7/7)*dk6 + dk5*CD.k0*t6 + 3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 + 3*CD.dk*k5*t2 + 9*dk2*CD.k0*t2 +
             k6+9*k2*CD.dk ) * L ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::bbTriangle( valueType offs,
                             valueType p0[2],
                             valueType p1[2],
                             valueType p2[2] ) const {
    valueType theta_max = CD.theta( L ) ;
    valueType theta_min = CD.theta0 ;
    valueType dtheta    = std::abs( theta_max-theta_min ) ;
    if ( dtheta < m_pi_2 ) {
      valueType alpha, t0[2] ;
      eval( 0, offs, p0[0], p0[1] ) ;
      eval_D( 0, t0[0], t0[1] ) ; // no offset
      if ( dtheta > 0.0001 * m_pi_2 ) {
        valueType t1[2] ;
        eval( L, offs, p1[0], p1[1] ) ;
        eval_D( L, t1[0], t1[1] ) ; // no offset
        // risolvo il sistema
        // p0 + alpha * t0 = p1 + beta * t1
        // alpha * t0 - beta * t1 = p1 - p0
        valueType det = t1[0]*t0[1]-t0[0]*t1[1] ;
        alpha = ((p1[1]-p0[1])*t1[0] - (p1[0]-p0[0])*t1[1])/det ;
      } else {
        // se angolo troppo piccolo uso approx piu rozza
        alpha = L ;
      }
      p2[0] = p0[0] + alpha*t0[0] ;
      p2[1] = p0[1] + alpha*t0[1] ;
      return true ;
    } else {
      return false ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::bbSplit(
    valueType               split_angle,
    valueType               split_size,
    valueType               split_offs,
    vector<ClothoidCurve> & c,
    vector<T2D>           & t
  ) const {

    // step 0: controllo se curvatura passa per 0
    valueType k_min = theta_D( 0 ) ;
    valueType k_max = theta_D( L ) ;
    c.clear() ;
    t.clear() ;
    if ( k_min * k_max < 0 ) {
      // risolvo (s-s_min)*dk+k_min = 0 --> s = s_min-k_min/dk
      valueType s_med = -k_min/CD.dk ;
      ClothoidCurve tmp(*this) ;
      tmp.trim(0,s_med) ;
      tmp.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
      tmp.trim(s_med,L) ;
      tmp.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    } else {
      bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  valueType
  abs2pi( valueType a ) {
    a = std::abs(a) ;
    while ( a > m_pi ) a -= m_2pi ;
    return std::abs(a) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

    eval( 0, theta_min, kappa_min, x_min, y_min ) ;
    eval( L, theta_max, kappa_max, x_max, y_max ) ;

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
      valueType s_med = L/2 ;
      cc.trim(0,s_med) ;
      cc.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
      cc.trim(s_med,L) ;
      cc.bbSplit_internal( split_angle, split_size, split_offs, c, t ) ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::intersect_internal( ClothoidCurve & c1,
                                     valueType       c1_offs,
                                     valueType     & s1,
                                     ClothoidCurve & c2,
                                     valueType       c2_offs,
                                     valueType     & s2,
                                     indexType       max_iter,
                                     valueType       tolerance ) const {
    valueType angle1a = c1.theta(0) ;
    valueType angle1b = c1.theta(c1.L) ;
    valueType angle2a = c2.theta(0) ;
    valueType angle2b = c2.theta(c2.L) ;
    // cerca angoli migliori per partire
    valueType dmax = abs2pi(angle1a-angle2a) ;
    valueType dab  = abs2pi(angle1a-angle2b) ;
    valueType dba  = abs2pi(angle1b-angle2a) ;
    valueType dbb  = abs2pi(angle1b-angle2b) ;
    s1 = s2 = 0 ;
    if ( dmax < dab ) { dmax = dab ; s2 = c2.L ; }
    if ( dmax < dba ) { dmax = dba ; s1 = 0 ; s2 = 0 ; }
    if ( dmax < dbb ) {              s1 = 0 ; s2 = c2.L ; }
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
      if ( s1 <= 0 || s1 >= c1.L ||
           s2 <= 0 || s2 >= c2.L ) break ;
      if ( std::abs(px) <= tolerance ||
           std::abs(py) <= tolerance ) return true ;
    }
    return false ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    bbSplit( m_pi/50, L/3, offs, c0, t0 ) ;
    clot.bbSplit( m_pi/50, clot.L/3, clot_offs, c1, t1 ) ;
    s1.clear() ;
    s2.clear() ;
    for ( unsigned i = 0 ; i < unsigned(c0.size()) ; ++i ) {
      for ( unsigned j = 0 ; j < unsigned(c1.size()) ; ++j ) {
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // collision detection
  bool
  ClothoidCurve::approximate_collision( valueType             offs,
                                        ClothoidCurve const & clot,
                                        valueType             clot_offs,
                                        valueType             max_angle,
                                        valueType             max_size ) const {
    vector<ClothoidCurve> c0, c1 ;
    vector<T2D>           t0, t1 ;
    bbSplit( max_angle, max_size, offs, c0, t0 ) ;
    clot.bbSplit( max_angle, max_size, clot_offs, c1, t1 ) ;
    for ( unsigned i = 0 ; i < unsigned(c0.size()) ; ++i ) {
      for ( unsigned j = 0 ; j < unsigned(c1.size()) ; ++j ) {
        if ( t0[i].overlap(t1[j]) ) return true ;
      }
    }
    return false ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::ostream &
  operator << ( std::ostream & stream, ClothoidCurve const & c ) {
    stream <<   "x0     = " << c.CD.x0
           << "\ny0     = " << c.CD.y0
           << "\ntheta0 = " << c.CD.theta0
           << "\nk      = " << c.CD.k0
           << "\ndk     = " << c.CD.dk
           << "\nL      = " << c.L
           << "\n" ;
    return stream ;
  }

}

// EOF: Clothoid.cc
