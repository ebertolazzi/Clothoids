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
    CD.kappa0 = new_kappa ;
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
    CD.kappa0 /= s ;
    CD.dk     /= s*s ;
    L         *= s ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::reverse() {
    CD.theta0 += m_pi ;
    if ( CD.theta0 > m_pi ) CD.theta0 -= 2*m_pi ;
    CD.kappa0 = -CD.kappa0 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::thetaTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType kL  = CD.kappa0 ;
    valueType kR  = CD.kappa(L) ;
    valueType thL = 0 ;
    valueType thR = CD.deltaTheta(L) ;
    if ( kL*kR < 0 ) {
      valueType root = -CD.kappa0/CD.dk ;
      if ( root > 0 && root < L ) {
        valueType thM  = CD.deltaTheta(root) ;
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
    valueType kL  = CD.kappa0 ;
    valueType kR  = CD.kappa(L) ;
    valueType thL = 0 ;
    valueType thR = CD.deltaTheta(L) ;
    if ( thL < thR ) { thMin = thL ; thMax = thR ; }
    else             { thMin = thR ; thMax = thL ; }
    if ( kL*kR < 0 ) {
      valueType root = -CD.kappa0/CD.dk ;
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
    kMin = CD.kappa0 ;
    kMax = CD.kappa(L);
    if ( kMax < kMin ) std::swap( kMax, kMin ) ;
    return kMax - kMin ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::curvatureTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk ;
    valueType km = CD.kappa0 ;
    valueType kp = CD.kappa(L) ;
    return std::abs(kp-km) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralCurvature2() const {
    return L*( CD.kappa0*(CD.kappa0+L*CD.dk) + (L*L)*CD.dk*CD.dk/3 ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralJerk2() const {
    valueType k2 = CD.kappa0*CD.kappa0 ;
    valueType k3 = CD.kappa0*k2 ;
    valueType k4 = k2*k2 ;
    valueType t1 = L ;
    valueType t2 = L*t1 ;
    valueType t3 = L*t2 ;
    valueType t4 = L*t3 ;
    return ((((t4/5*CD.dk+t3*CD.kappa0)*CD.dk+(1+2*t2)*k2)*CD.dk+2*t1*k3)*CD.dk+k4)*L ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidCurve::integralSnap2() const {
    valueType k2  = CD.kappa0*CD.kappa0 ;
    valueType k3  = CD.kappa0*k2 ;
    valueType k4  = k2*k2 ;
    valueType k5  = k4*CD.kappa0 ;
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

    return ( (t7/7)*dk6 + dk5*CD.kappa0*t6 + 3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 + 3*CD.dk*k5*t2 + 9*dk2*CD.kappa0*t2 +
             k6+9*k2*CD.dk ) * L ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::bbSplit(
    valueType        split_angle,
    valueType        split_size,
    valueType        split_offs,
    vector<bbData> & bb
  ) const {

    // step 0: controllo se curvatura passa per 0
    valueType k_min = theta_D( 0 ) ;
    valueType k_max = theta_D( L ) ;

    bb.clear() ;

    bbData2 data ;
    data.split_angle = split_angle ;
    data.split_size  = split_size  ;
    data.split_offs  = split_offs  ;
    data.cd          = this->CD ;
    data.s0          = 0 ;

    if ( k_min * k_max < 0 ) {
      // risolvo (s-s_min)*dk+k_min = 0 --> s = s_min-k_min/dk
      valueType s_med = -k_min/CD.dk ;
      data.L  = s_med ;
      bbSplit_internal( data, bb ) ;
      // trim
      CD.eval( s_med,
               data.cd.theta0, data.cd.kappa0,
               data.cd.x0, data.cd.y0 ) ;
      data.s0 = s_med ;
      data.L  = this->L - s_med ;
      bbSplit_internal( data, bb ) ;
    } else {
      data.L  = this->L ;
      bbSplit_internal( data, bb ) ;
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
    bbData2 const  & data,
    vector<bbData> & bbV
  ) const {

    valueType theta_min, kappa_min, x_min, y_min,
              theta_max, kappa_max, x_max, y_max ;

    data.cd.eval( 0,      theta_min, kappa_min, x_min, y_min ) ;
    data.cd.eval( data.L, theta_max, kappa_max, x_max, y_max ) ;

    valueType dtheta = std::abs( theta_max - theta_min ) ;
    valueType dx     = x_max - x_min ;
    valueType dy     = y_max - y_min ;
    valueType len    = hypot( dy, dx ) ;
    valueType dangle = abs2pi(atan2( dy, dx )-theta_min) ;
    if ( dtheta          <= data.split_angle &&
         len*tan(dangle) <= data.split_size ) {
      bbData bb ;
      valueType p0[2], p1[2], p2[2] ;
      bool ok = data.cd.bbTriangle( data.L, data.split_offs, p0, p1, p2 ) ;
      G2LIB_ASSERT( ok, "ClothoidCurve::bbSplit_internal, bad bounding box" ) ;
      bb.t.setup( p0, p1, p2 ) ;
      bb.s0 = data.s0 ;
      bb.L  = data.L ;
      bb.cd = data.cd ;
      bbV.push_back(bb) ;
    } else {
      bbData2 d ;
      valueType Lh = data.L / 2 ;
      d.split_angle = data.split_angle ;
      d.split_size  = data.split_size ;
      d.split_offs  = data.split_offs ;
      d.s0          = data.s0 ;
      d.L           = Lh ;
      d.cd          = data.cd ;
      bbSplit_internal( d, bbV ) ;

      // trim
      data.cd.eval( Lh,
                    d.cd.theta0, d.cd.kappa0,
                    d.cd.x0, d.cd.y0 ) ;
      d.s0 = data.s0 + Lh  ;
      d.L  = Lh ;
      bbSplit_internal( d, bbV ) ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::intersect_internal( bbData const & c1,
                                     valueType      c1_offs,
                                     valueType    & s1,
                                     bbData const & c2,
                                     valueType      c2_offs,
                                     valueType    & s2,
                                     indexType      max_iter,
                                     valueType      tolerance ) const {
    valueType angle1a = c1.cd.theta(0) ;
    valueType angle1b = c1.cd.theta(c1.L) ;
    valueType angle2a = c2.cd.theta(0) ;
    valueType angle2b = c2.cd.theta(c2.L) ;
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
      c1.cd.eval  ( s1, c1_offs, p1[0], p1[1] ) ;
      c1.cd.eval_D( s1, c1_offs, t1[0], t1[1] ) ;
      c2.cd.eval  ( s2, c2_offs, p2[0], p2[1] ) ;
      c2.cd.eval_D( s2, c2_offs, t2[0], t2[1] ) ;
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
    vector<bbData> bbV0,  bbV1 ;
    bbSplit( m_pi/50, L/3, offs, bbV0 ) ;
    clot.bbSplit( m_pi/50, clot.L/3, clot_offs, bbV1 ) ;
    s1.clear() ;
    s2.clear() ;
    for ( unsigned i = 0 ; i < unsigned(bbV0.size()) ; ++i ) {
      bbData const & bbi = bbV0[i] ;
      for ( unsigned j = 0 ; j < unsigned(bbV1.size()) ; ++j ) {
        bbData const & bbj = bbV1[j] ;
        if ( bbi.t.overlap(bbj.t) ) {
          // uso newton per cercare intersezione
          valueType tmp_s1, tmp_s2 ;
          bool ok = intersect_internal( bbi, offs,      tmp_s1,
                                        bbj, clot_offs, tmp_s2,
                                        max_iter, tolerance ) ;
          if ( ok ) {
            s1.push_back(bbi.s0+tmp_s1) ;
            s2.push_back(bbj.s0+tmp_s2) ;
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
    vector<bbData> bbV0, bbV1 ;
    bbSplit( max_angle, max_size, offs, bbV0 ) ;
    clot.bbSplit( max_angle, max_size, clot_offs, bbV1 ) ;
    for ( unsigned i = 0 ; i < unsigned(bbV0.size()) ; ++i ) {
      bbData & bbi = bbV0[i] ;
      for ( unsigned j = 0 ; j < unsigned(bbV1.size()) ; ++j ) {
        bbData & bbj = bbV1[j] ;
        if ( bbi.t.overlap(bbj.t) ) return true ;
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
           << "\nkappa0 = " << c.CD.kappa0
           << "\ndk     = " << c.CD.dk
           << "\nL      = " << c.L
           << "\n" ;
    return stream ;
  }

}

// EOF: Clothoid.cc
