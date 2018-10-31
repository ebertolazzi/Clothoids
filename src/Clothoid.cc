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

  using std::vector;
  using std::abs;
  using std::swap;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::bbox( real_type & xmin,
                       real_type & ymin,
                       real_type & xmax,
                       real_type & ymax ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::bbox" ) ;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::bbox( real_type   offs,
                       real_type & xmin,
                       real_type & ymin,
                       real_type & xmax,
                       real_type & ymax ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::bbox" ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::collision( BaseCurve const & obj ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::collision" );
    return false;
  }

  bool
  ClothoidCurve::collision( real_type         offs,
                            BaseCurve const & obj,
                            real_type         offs_obj ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::collision" ) ;
    return false;
  }

  void
  ClothoidCurve::intersect( BaseCurve const & obj,
                            IntersectList   & ilist,
                            bool              swap_s_vals ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::intersect" ) ;
  }

  void
  ClothoidCurve::intersect( real_type         offs,
                            BaseCurve const & obj,
                            real_type         offs_obj,
                            IntersectList   & ilist,
                            bool              swap_s_vals ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::intersect" ) ;
  }

  int_type
  ClothoidCurve::projection( real_type   qx,
                             real_type   qy,
                             real_type & x,
                             real_type & y,
                             real_type & s ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::projection" ) ;
    return 0;
  }

  int_type // true if projection is unique and orthogonal
  ClothoidCurve::projection( real_type   qx,
                             real_type   qy,
                             real_type   offs,
                             real_type & x,
                             real_type & y,
                             real_type & s ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::projection" ) ;
    return 0;
  }

  real_type
  ClothoidCurve::closestPoint( real_type   qx,
                               real_type   qy,
                               real_type   offs,
                               real_type & x,
                               real_type & y,
                               real_type & s ) const {
    G2LIB_ASSERT( false, "DA FARE ClothoidCurve::closestPoint" ) ;
    return 0;
  }

  static real_type const one_degree = m_pi/180;
  static real_type const n90_degree = 90*one_degree;
  static int_type  const max_level  = 10;

  void
  ClothoidCurve::bbTriangles( real_type              offs,
                              vector<Triangle2D>   & tvec,
                              bbDataForSplit const & data,
                              real_type              max_angle,
                              real_type              max_size,
                              int_type               level ) const {

    real_type dtheta = abs(data.theta0 - data.theta1);
    real_type dx     = data.x1-data.x0;
    real_type dy     = data.y1-data.y0;
    if ( level >= max_level || (hypot(dx,dy) <= max_size && dtheta <= max_angle) ) {
      real_type tx0   = cos(data.theta0);
      real_type ty0   = sin(data.theta0);
      real_type alpha = data.s1-data.s0; // se angolo troppo piccolo uso approx piu rozza
      if ( dtheta > one_degree && dtheta < n90_degree ) {
        real_type tx1 = cos(data.theta1);
        real_type ty1 = sin(data.theta1);
        real_type det = tx1*ty0-tx0*ty1;
        alpha = (dy*tx1 - dx*ty1)/det;
      }
      real_type x2 = data.x0 + alpha*tx0;
      real_type y2 = data.y0 + alpha*ty0;
      //real_type x2 = (data.x0 + data.x1)/2;
      //real_type y2 = (data.y0 + data.y1)/2;
      Triangle2D t( data.x0, data.y0,
                    x2, y2,
                    data.x1, data.y1 );
      tvec.push_back( t );
      return;
    }

    real_type sM     = (data.s0+data.s1)/2;
    real_type thetaM = CD.theta(sM);
    real_type xM, yM;
    CD.eval( sM, offs, xM, yM );

    bbDataForSplit dataNew;
    dataNew.s0     = data.s0;
    dataNew.s1     = sM;
    dataNew.theta0 = data.theta0;
    dataNew.theta1 = thetaM;
    dataNew.x0     = data.x0;
    dataNew.y0     = data.y0;
    dataNew.x1     = xM;
    dataNew.y1     = yM;
    bbTriangles( offs, tvec, dataNew, max_angle, max_size, level+1 );

    dataNew.s0     = sM;
    dataNew.s1     = data.s1;
    dataNew.theta0 = thetaM;
    dataNew.theta1 = data.theta1;
    dataNew.x0     = xM;
    dataNew.y0     = yM;
    dataNew.x1     = data.x1;
    dataNew.y1     = data.y1;
    bbTriangles( offs, tvec, dataNew, max_angle, max_size, level+1 );
  }

  void
  ClothoidCurve::bbTriangles( real_type            offs,
                              vector<Triangle2D> & tvec,
                              real_type            max_angle,
                              real_type            max_size ) const {

    bbDataForSplit data;
    data.s0     = 0;
    data.theta0 = CD.theta0;
    CD.eval( 0, offs, data.x0, data.y0 ) ;

    if ( CD.kappa0*CD.dk >= 0 || CD.dk*CD.kappa(L) <= 0 ) {
      data.s1     = L;
      data.theta1 = CD.theta(L);
      CD.eval( L, offs, data.x1, data.y1 ) ;
      bbTriangles( offs, tvec, data, max_angle, max_size, 0 );
    } else {
      // flex inside, split clothoid
      real_type sflex = -CD.kappa0/CD.dk;
      data.s1     = sflex;
      data.theta1 = CD.theta(sflex);
      CD.eval( sflex, offs, data.x1, data.y1 ) ;
      bbTriangles( offs, tvec, data, max_angle, max_size, 0 );
      data.s0     = data.s1;
      data.theta0 = data.theta1;
      data.x0     = data.x1;
      data.y0     = data.y1;
      data.s1     = L;
      data.theta1 = CD.theta(L);
      CD.eval( L, offs, data.x1, data.y1 ) ;
      bbTriangles( offs, tvec, data, max_angle, max_size, 0 );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::thetaTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = CD.kappa0;
    real_type kR  = CD.kappa(L);
    real_type thL = 0;
    real_type thR = CD.deltaTheta(L);
    if ( kL*kR < 0 ) {
      real_type root = -CD.kappa0/CD.dk;
      if ( root > 0 && root < L ) {
        real_type thM  = CD.deltaTheta(root);
        return abs( thR - thM ) + abs( thM - thL );
      }
    }
    return abs( thR - thL );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::thetaMinMax( real_type & thMin, real_type & thMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = CD.kappa0;
    real_type kR  = CD.kappa(L);
    real_type thL = 0;
    real_type thR = CD.deltaTheta(L);
    if ( thL < thR ) { thMin = thL; thMax = thR; }
    else             { thMin = thR; thMax = thL; }
    if ( kL*kR < 0 ) {
      real_type root = -CD.kappa0/CD.dk;
      if ( root > 0 && root < L ) {
        real_type thM = CD.deltaTheta(root);
        if      ( thM < thMin ) thMin = thM;
        else if ( thM > thMax ) thMax = thM;
      }
    }
    return thMax - thMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvatureMinMax( real_type & kMin, real_type & kMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    kMin = CD.kappa0;
    kMax = CD.kappa(L);
    if ( kMax < kMin ) swap( kMax, kMin );
    return kMax - kMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvatureTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type km = CD.kappa0;
    real_type kp = CD.kappa(L);
    return abs(kp-km);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralCurvature2() const {
    return L*( CD.kappa0*(CD.kappa0+L*CD.dk) + (L*L)*CD.dk*CD.dk/3 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralJerk2() const {
    real_type k2 = CD.kappa0*CD.kappa0;
    real_type k3 = CD.kappa0*k2;
    real_type k4 = k2*k2;
    real_type t1 = L;
    real_type t2 = L*t1;
    real_type t3 = L*t2;
    real_type t4 = L*t3;
    return ((((t4/5*CD.dk+t3*CD.kappa0)*CD.dk+(1+2*t2)*k2)*CD.dk+2*t1*k3)*CD.dk+k4)*L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralSnap2() const {
    real_type k2  = CD.kappa0*CD.kappa0;
    real_type k3  = CD.kappa0*k2;
    real_type k4  = k2*k2;
    real_type k5  = k4*CD.kappa0;
    real_type k6  = k4*k2;
    real_type dk2 = CD.dk*CD.dk;
    real_type dk3 = CD.dk*dk2;
    real_type dk4 = dk2*dk2;
    real_type dk5 = dk4*CD.dk;
    real_type dk6 = dk4*dk2;
    real_type t2  = L;
    real_type t3  = L*t2;
    real_type t4  = L*t3;
    real_type t5  = L*t4;
    real_type t6  = L*t5;
    real_type t7  = L*t6;

    return ( (t7/7)*dk6 + dk5*CD.kappa0*t6 + 3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 + 3*CD.dk*k5*t2 + 9*dk2*CD.kappa0*t2 +
             k6+9*k2*CD.dk ) * L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::bbSplit(
    real_type        split_angle,
    real_type        split_size,
    real_type        split_offs,
    vector<bbData> & bb,
    bool             reset_bb
  ) const {

    // step 0: controllo se curvatura passa per 0
    real_type k_min = theta_D( 0 );
    real_type k_max = theta_D( L );

    if ( reset_bb ) bb.clear();

    bbData2 data;
    data.split_angle = split_angle;
    data.split_size  = split_size;
    data.split_offs  = split_offs;
    data.cd          = this->CD;
    data.s0          = 0;

    if ( k_min * k_max < 0 ) {
      // risolvo (s-s_min)*dk+k_min = 0 --> s = s_min-k_min/dk
      real_type s_med = -k_min/CD.dk;
      data.L  = s_med;
      bbSplit_internal( data, bb );
      // trim
      CD.evaluate( s_med,
                   data.cd.theta0, data.cd.kappa0,
                   data.cd.x0, data.cd.y0 );
      data.s0 = s_med;
      data.L  = this->L - s_med;
      bbSplit_internal( data, bb );
    } else {
      data.L  = this->L;
      bbSplit_internal( data, bb );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  real_type
  abs2pi( real_type a ) {
    a = abs(a);
    while ( a > m_pi ) a -= m_2pi;
    return abs(a);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::bbSplit_internal(
    bbData2 const  & data,
    vector<bbData> & bbV
  ) const {

    real_type theta_min, kappa_min, x_min, y_min,
              theta_max, kappa_max, x_max, y_max;

    data.cd.evaluate( 0,      theta_min, kappa_min, x_min, y_min );
    data.cd.evaluate( data.L, theta_max, kappa_max, x_max, y_max );

    real_type dtheta = abs( theta_max - theta_min );
    real_type dx     = x_max - x_min;
    real_type dy     = y_max - y_min;
    real_type len    = hypot( dy, dx );
    real_type dangle = abs2pi(atan2( dy, dx )-theta_min);
    if ( dtheta          <= data.split_angle &&
         len*tan(dangle) <= data.split_size ) {
      bbData bb;
      real_type x0, y0, x1, y1, x2, y2;
      bool ok = data.cd.bbTriangle( data.L, data.split_offs,
                                    x0, y0, x1, y1, x2, y2 );
      G2LIB_ASSERT( ok, "ClothoidCurve::bbSplit_internal, bad bounding box" );
      bb.t.build( x0, y0, x1, y1, x2, y2 );
      bb.s0 = data.s0;
      bb.L  = data.L;
      bb.cd = data.cd;
      bbV.push_back(bb);
    } else {
      bbData2 d;
      real_type Lh = data.L / 2;
      d.split_angle = data.split_angle;
      d.split_size  = data.split_size;
      d.split_offs  = data.split_offs;
      d.s0          = data.s0;
      d.L           = Lh;
      d.cd          = data.cd;
      bbSplit_internal( d, bbV );

      // trim
      data.cd.evaluate( Lh,
                        d.cd.theta0, d.cd.kappa0,
                        d.cd.x0, d.cd.y0 );
      d.s0 = data.s0 + Lh;
      d.L  = Lh;
      bbSplit_internal( d, bbV );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::intersect_internal(
    bbData const & c1,
    real_type      c1_offs,
    real_type    & s1,
    bbData const & c2,
    real_type      c2_offs,
    real_type    & s2,
    int_type       max_iter,
    real_type      tolerance
  ) const {
    real_type angle1a = c1.cd.theta(0);
    real_type angle1b = c1.cd.theta(c1.L);
    real_type angle2a = c2.cd.theta(0);
    real_type angle2b = c2.cd.theta(c2.L);
    // cerca angoli migliori per partire
    real_type dmax = abs2pi(angle1a-angle2a);
    real_type dab  = abs2pi(angle1a-angle2b);
    real_type dba  = abs2pi(angle1b-angle2a);
    real_type dbb  = abs2pi(angle1b-angle2b);
    s1 = s2 = 0;
    if ( dmax < dab ) { dmax = dab; s2 = c2.L; }
    if ( dmax < dba ) { dmax = dba; s1 = 0; s2 = 0; }
    if ( dmax < dbb ) {              s1 = 0; s2 = c2.L; }
    int_type nout = 0;
    for ( int_type i = 0; i < max_iter; ++i ) {
      real_type t1[2], t2[2], p1[2], p2[2];
      c1.cd.eval  ( s1, c1_offs, p1[0], p1[1] );
      c1.cd.eval_D( s1, c1_offs, t1[0], t1[1] );
      c2.cd.eval  ( s2, c2_offs, p2[0], p2[1] );
      c2.cd.eval_D( s2, c2_offs, t2[0], t2[1] );
      /*
      // risolvo il sistema
      // p1 + alpha * t1 = p2 + beta * t2
      // alpha * t1 - beta * t2 = p2 - p1
      //
      //  / t1[0] -t2[0] \ / alpha \ = / p2[0] - p1[0] \
      //  \ t1[1] -t2[1] / \ beta  /   \ p2[1] - p1[1] /
      */
      real_type det = t2[0]*t1[1]-t1[0]*t2[1];
      real_type px  = p2[0]-p1[0];
      real_type py  = p2[1]-p1[1];
      s1 += (py*t2[0] - px*t2[1])/det;
      s2 += (t1[0]*py - t1[1]*px)/det;
      if ( ! ( isfinite(s1) && isfinite(s1) ) ) break;
      bool out = false;
      if      ( s1 <= 0    ) { out = true; s1 = 0; }
      else if ( s1 >= c1.L ) { out = true; s1 = c1.L; }
      if      ( s2 <= 0    ) { out = true; s2 = 0; }
      else if ( s2 >= c2.L ) { out = true; s2 = c2.L; }
      if ( out ) {
        if ( ++nout > 3 ) break;
      } else {
        if ( abs(px) <= tolerance &&
             abs(py) <= tolerance ) return true;
      }
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::intersect(
    real_type             offs,
    ClothoidCurve const & clot,
    real_type             clot_offs,
    vector<real_type>   & s1,
    vector<real_type>   & s2,
    int_type              max_iter,
    real_type             tolerance
  ) const {
    vector<bbData> bbV0, bbV1;
    bbSplit( m_pi/50, L/3, offs, bbV0 );
    clot.bbSplit( m_pi/50, clot.L/3, clot_offs, bbV1 );
    s1.clear();
    s2.clear();
    for ( unsigned i = 0; i < unsigned(bbV0.size()); ++i ) {
      bbData const & bbi = bbV0[i];
      for ( unsigned j = 0; j < unsigned(bbV1.size()); ++j ) {
        bbData const & bbj = bbV1[j];
        if ( bbi.t.overlap(bbj.t) ) {
          // uso newton per cercare intersezione
          real_type tmp_s1, tmp_s2;
          bool ok = intersect_internal( bbi, offs,      tmp_s1,
                                        bbj, clot_offs, tmp_s2,
                                        max_iter, tolerance );
          if ( ok ) {
            s1.push_back(bbi.s0+tmp_s1);
            s2.push_back(bbj.s0+tmp_s2);
          }
        }
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // collision detection
  bool
  ClothoidCurve::approximate_collision(
    real_type             offs,
    ClothoidCurve const & clot,
    real_type             clot_offs,
    real_type             max_angle,
    real_type             max_size
  ) const {
    vector<bbData> bbV0, bbV1;
    bbSplit( max_angle, max_size, offs, bbV0 );
    clot.bbSplit( max_angle, max_size, clot_offs, bbV1 );
    for ( unsigned i = 0; i < unsigned(bbV0.size()); ++i ) {
      bbData & bbi = bbV0[i];
      for ( unsigned j = 0; j < unsigned(bbV1.size()); ++j ) {
        bbData & bbj = bbV1[j];
        if ( bbi.t.overlap(bbj.t) ) return true;
      }
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::findST( real_type   x,
                         real_type   y,
                         real_type & s,
                         real_type & t ) const {
    real_type X, Y, nx, ny;
    real_type d = closestPoint( x, y, X, Y, s );
    nor( s, nx, ny );
    t = nx*(x-X) + ny*(y-Y);
    // check if projection is orthogonal on the curve
    #if 0
      real_type abst = abs(t);
      return abs(d-abst) <= machepsi1000*(1+abst);
    #else
      eval( s, t, X, Y );
      real_type err = hypot( x-X, y-Y );
      return err < 1e-8*(1+d);
    #endif
    //return abs(d-abst) <= 1e-3*(1+abst);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, ClothoidCurve const & c ) {
    stream <<   "x0     = " << c.CD.x0
           << "\ny0     = " << c.CD.y0
           << "\ntheta0 = " << c.CD.theta0
           << "\nkappa0 = " << c.CD.kappa0
           << "\ndk     = " << c.CD.dk
           << "\nL      = " << c.L
           << "\n";
    return stream;
  }

}

// EOF: Clothoid.cc
