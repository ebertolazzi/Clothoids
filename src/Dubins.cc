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
/// file: Dubins.cc
///

#include "Clothoids.hh"
#include "PolynomialRoots.hh"

namespace G2lib {

  using PolynomialRoots::Quadratic;

  /*\
   |   ____        _     _
   |  |  _ \ _   _| |__ (_)_ __  ___
   |  | | | | | | | '_ \| | '_ \/ __|
   |  | |_| | |_| | |_) | | | | \__ \
   |  |____/ \__,_|_.__/|_|_| |_|___/
  \*/

  static
  void
  into_0_2pi( real_type & a ) {
    while ( a < 0            ) a += Utils::m_2pi;
    while ( a > Utils::m_2pi ) a -= Utils::m_2pi;
  }

  static
  void
  minus_pi_pi( real_type & a ) {
    while ( a < -Utils::m_pi ) a += Utils::m_2pi;
    while ( a >  Utils::m_pi ) a -= Utils::m_2pi;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::setup( GenericContainer const & gc ) {
    // @@@@@@@@@@@ DA FARE @@@@@@@@@@@@@@
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Dubins::build( LineSegment const & )   { UTILS_ERROR("can convert from LineSegment to Dubins\n"); }
  void Dubins::build( CircleArc const & )      { UTILS_ERROR("can convert from CircleArc to Dubins\n"); }
  void Dubins::build( Biarc const & )          { UTILS_ERROR("can convert from Biarc to Dubins\n"); }
  void Dubins::build( ClothoidCurve const & )  { UTILS_ERROR("can convert from ClothoidCurve to Dubins\n"); }
  void Dubins::build( PolyLine const & )       { UTILS_ERROR("can convert from PolyLine to Dubins\n"); }
  void Dubins::build( BiarcList const & )      { UTILS_ERROR("can convert from BiarcList to Dubins\n"); }
  void Dubins::build( ClothoidList const & )   { UTILS_ERROR("can convert from ClothoidList to Dubins\n"); }

  void
  Dubins::build( Dubins const & DB ) {
    m_C0 = DB.m_C0;
    m_C1 = DB.m_C1;
    m_C2 = DB.m_C2;
    m_solution_type = DB.m_solution_type;
  }

  bool
  Dubins::build(
    real_type x0,
    real_type y0,
    real_type theta0,

    real_type x3,
    real_type y3,
    real_type theta3,

    real_type k_max
  ) {

    // setup in standard form
    real_type dx    = x3-x0;
    real_type dy    = y3-y0;
    real_type th    = atan2( dy, dx );
    real_type alpha = theta0-th;
    real_type beta  = theta3-th;
    real_type L     = k_max * hypot( dx, dy );

    // metto angolo in -Pi, Pi
    minus_pi_pi( alpha );
    minus_pi_pi( beta );

    // common constants
    real_type Ca  = cos(alpha);
    real_type Sa  = sin(alpha);
    real_type Cb  = cos(beta);
    real_type Sb  = sin(beta);
    real_type ab  = alpha-beta;
    real_type Cab = cos(ab);
    real_type Sab = sin(ab);
    real_type LCa = L*Ca;
    real_type LSa = L*Sa;

    real_type min_len = Utils::Inf<real_type>();
    real_type L0{0}, L1{0}, L2{0}, s0{0}, s1{0}, s2{0};

    // LSL, RSR, LSR, RSL
    real_type d0[4]{1,-1,1,-1};
    real_type d2[4]{1,-1,-1,1};
    DubinsType stype[4]{DubinsType::LSL, DubinsType::RSR, DubinsType::LSR, DubinsType::RSL };

    for ( integer i = 0; i < 4; ++i ) {
      real_type d__0 = d0[i];
      real_type d__2 = d2[i];
      Quadratic Q( (Cab-1)*d__2 - LSa,
                   2*d__0*(LCa+d__2*Sab),
                   LSa + 2*d__0 - (Cab+1)*d__2 );
      real_type Z[2];
      integer nr = Q.get_real_roots( Z );
      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type z  = Z[ir];
        real_type l0 = Utils::m_pi+2*atan(z);
        real_type z2 = z*z;
        real_type C0 = (z2-1)/(z2+1);
        real_type S0 = -2*z/(z2+1);
        real_type l1 = (LCa + Sab*d__2)*C0 + ((d__2*Cab-LSa)*d__0-1)*S0;
        if ( l1 < 0 ) continue;
        real_type l2 = -d__2*(ab+d__0*l0);
        into_0_2pi(l2);
        real_type len = l0+l1+l2;
        if ( len < min_len ) {
          min_len = len;
          L0 = l0;   L1 = l1; L2 = l2;
          s0 = d__0; s1 = 0;  s2 = d__2;
          m_solution_type = stype[i];
        }
      }
    }

    if ( L < 4 ) {
      real_type dpm[2]{1,-1};
      DubinsType stype[2]{ DubinsType::LRL, DubinsType::RLR };
      for ( integer i = 0; i < 2; ++i ) {
        real_type d = dpm[i];
        Quadratic Q( (L/2 + d*(Sa+Sb))*L + Cab - 1,
                     4*(L*Cb + d*Sab),
                     (L/2 + (Sa - 3*Sb)*d)*L - 3*(Cab-1) );
        real_type Z[2];
        integer nr = Q.get_real_roots( Z );
        for ( integer ir = 0; ir < nr; ++ir ) {
          real_type z  = Z[ir];
          real_type l2 = Utils::m_pi+2*atan(z);
          real_type z2 = z*z;
          real_type C0 = ((LSa*d + Cab + 1)*z2 + 4*z*d*Sab + LSa*d - 3*Cab + 1); // /(2*(z2+1));
          real_type S0 = ((LCa - d*Sab)*z2 + 4*z*Cab + LCa + 3*d*Sab);           // /(2*(z2+1));
          real_type l0 = atan2( S0, C0 );
          if ( l0 < 0 ) l0 += Utils::m_2pi;
          real_type l1 = l0+l2+d*ab;
          into_0_2pi(l1);
          real_type len = l0+l1+l2;
          if ( len < min_len ) {
            min_len = len;
            L0 = l0; L1 = l1; L2 = l2;
            s0 = d;  s1 = -d; s2 = d;
            m_solution_type = stype[i];
          }
        }
      }
    }

    if ( min_len == Utils::Inf<real_type>() ) return false;
    m_C0.build( x0,           y0,           theta0,           s0*k_max, L0/k_max );
    m_C1.build( m_C0.x_end(), m_C0.y_end(), m_C0.theta_end(), s1*k_max, L1/k_max );
    m_C2.build( m_C1.x_end(), m_C1.y_end(), m_C1.theta_end(), s2*k_max, L2/k_max );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::length() const {
    return m_C0.length()+m_C1.length()+m_C2.length();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::length_ISO( real_type offs ) const {
    return m_C0.length_ISO( offs )+m_C1.length_ISO( offs )+m_C2.length_ISO( offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::theta( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.theta(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.theta(s);
    s -= m_C2.length();
    return m_C2.theta(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::theta_D( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.theta_D(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.theta_D(s);
    s -= m_C1.length();
    return m_C2.theta_D(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::X( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.X(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.X(s);
    s -= m_C1.length();
    return m_C2.X(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::Y( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.Y(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.Y(s);
    s -= m_C1.length();
    return m_C2.Y(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::X_D( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.X_D(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.X_D(s);
    s -= m_C1.length();
    return m_C2.X_D(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::Y_D( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.Y_D(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.Y_D(s);
    s -= m_C1.length();
    return m_C2.Y_D(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::X_DD( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.X_DD(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.X_DD(s);
    s -= m_C1.length();
    return m_C2.X_DD(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::Y_DD( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.Y_DD(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.Y_DD(s);
    s -= m_C1.length();
    return m_C2.Y_DD(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::X_DDD( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.X_DDD(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.X_DDD(s);
    s -= m_C1.length();
    return m_C2.X_DDD(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  Dubins::Y_DDD( real_type s ) const {
    if ( s < m_C0.length() ) return m_C0.Y_DDD(s);
    s -= m_C0.length();
    if ( s < m_C1.length() ) return m_C1.Y_DDD(s);
    s -= m_C1.length();
    return m_C2.Y_DDD(s);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval(
    real_type   s,
    real_type & theta,
    real_type & kappa,
    real_type & x,
    real_type & y
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.evaluate( s, theta, kappa, x, y );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.evaluate( s, theta, kappa, x, y );
      } else {
        s -= m_C1.length();
        m_C2.evaluate( s, theta, kappa, x, y );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval(s, x, y );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval(s, x, y );
      } else {
        s -= m_C1.length();
        m_C2.eval(s, x, y );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_D(s, x_D, y_D );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_D(s, x_D, y_D );
      } else {
        s -= m_C1.length();
        m_C2.eval_D(s, x_D, y_D );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_DD(s, x_DD, y_DD );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_DD(s, x_DD, y_DD );
      } else {
        s -= m_C1.length();
        m_C2.eval_DD(s, x_DD, y_DD );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_DDD(s, x_DDD, y_DDD );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_DDD(s, x_DDD, y_DDD );
      } else {
        s -= m_C1.length();
        m_C2.eval_DDD(s, x_DDD, y_DDD );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void
  Dubins::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_ISO( s, offs, x, y );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_ISO( s, offs, x, y );
      } else {
        s -= m_C1.length();
        m_C2.eval_ISO( s, offs, x, y );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_ISO_D( s, offs, x_D, y_D );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_ISO_D( s, offs, x_D, y_D );
      } else {
        s -= m_C1.length();
        m_C2.eval_ISO_D( s, offs, x_D, y_D );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_ISO_DD( s, offs, x_DD, y_DD );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_ISO_DD( s, offs, x_DD, y_DD );
      } else {
        s -= m_C1.length();
        m_C2.eval_ISO_DD( s, offs, x_DD, y_DD );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    if ( s < m_C0.length() ) {
      m_C0.eval_ISO_DDD( s, offs, x_DDD, y_DDD );
    } else {
      s -= m_C0.length();
      if ( s < m_C1.length() ) {
        m_C1.eval_ISO_DDD( s, offs, x_DDD, y_DDD );
      } else {
        s -= m_C1.length();
        m_C2.eval_ISO_DDD( s, offs, x_DDD, y_DDD );
      }
    }
  }

  void
  Dubins::reverse() {
    CircleArc TMP(m_C0); m_C0 = m_C2; m_C2 = TMP;
    m_C0.reverse();
    m_C1.reverse();
    m_C2.reverse();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::scale( real_type scl ) {
    m_C0.scale( scl );
    m_C1.scale( scl );
    m_C1.change_origin( m_C0.x_end(), m_C0.y_end() );
  }

  void
  Dubins::change_origin( real_type newx0, real_type newy0 ) {
    m_C0.change_origin(newx0,newy0);
    m_C1.change_origin(m_C0.x_end(),m_C0.y_end());
    m_C2.change_origin(m_C1.x_end(),m_C1.y_end());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_C0.bbox( xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_C1.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
    m_C2.bbox( xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::bbox_ISO(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    m_C0.bbox_ISO( offs, xmin, ymin, xmax, ymax );
    real_type xmi1, ymi1, xma1, yma1;
    m_C1.bbox_ISO( offs, xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
    m_C2.bbox_ISO( offs, xmi1, ymi1, xma1, yma1 );
    if ( xmi1 < xmin ) xmin = xmi1;
    if ( xma1 > xmax ) xmax = xma1;
    if ( ymi1 < ymin ) ymin = ymi1;
    if ( yma1 > ymax ) ymax = yma1;
  }

  /*\
   |             _ _ _     _
   |    ___ ___ | | (_)___(_) ___  _ __
   |   / __/ _ \| | | / __| |/ _ \| '_ \
   |  | (_| (_) | | | \__ \ | (_) | | | |
   |   \___\___/|_|_|_|___/_|\___/|_| |_|
  \*/

  bool
  Dubins::collision( Dubins const & B ) const {
    return m_C0.collision( B.m_C0 ) ||
           m_C0.collision( B.m_C1 ) ||
           m_C0.collision( B.m_C2 ) ||

           m_C1.collision( B.m_C0 ) ||
           m_C1.collision( B.m_C1 ) ||
           m_C1.collision( B.m_C2 ) ||

           m_C2.collision( B.m_C0 ) ||
           m_C2.collision( B.m_C1 ) ||
           m_C2.collision( B.m_C2 );
  }

  bool
  Dubins::collision_ISO(
    real_type      offs,
    Dubins const & B,
    real_type      offs_B
  ) const {
    return m_C0.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C0.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C0.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C1.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C1.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C2, offs_B ) ||

           m_C2.collision_ISO( offs, B.m_C0, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C1, offs_B ) ||
           m_C2.collision_ISO( offs, B.m_C2, offs_B );
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Dubins::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      return this->collision( C );
    } else {
      return G2lib::collision( this, pC );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Dubins::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      return G2lib::collision_ISO( this, offs, pC, offs_C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::intersect(
    Dubins const  & B,
    IntersectList & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist02,
                  ilist10, ilist11, ilist12,
                  ilist20, ilist21, ilist22;

    m_C0.intersect( B.m_C0, ilist00 );
    m_C0.intersect( B.m_C1, ilist01 );
    m_C0.intersect( B.m_C2, ilist02 );

    m_C1.intersect( B.m_C0, ilist10 );
    m_C1.intersect( B.m_C1, ilist11 );
    m_C1.intersect( B.m_C2, ilist12 );

    m_C2.intersect( B.m_C0, ilist20 );
    m_C2.intersect( B.m_C1, ilist21 );
    m_C2.intersect( B.m_C2, ilist22 );

    real_type L0  = m_C0.length();
    real_type L1  = L0+m_C1.length();
    real_type LB0 = B.m_C0.length();
    real_type LB1 = LB0+B.m_C1.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +
                   ilist02.size() +

                   ilist10.size() +
                   ilist11.size() +
                   ilist12.size() +

                   ilist20.size() +
                   ilist21.size() +
                   ilist22.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist02 ) { it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist12 ) { it.first += L0; it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist20 ) { it.first += L1; ilist.push_back( it ); }
    for ( auto & it : ilist21 ) { it.first += L1; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist22 ) { it.first += L1; it.second += LB1; ilist.push_back( it ); }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Dubins::intersect_ISO(
    real_type       offs,
    Dubins const  & B,
    real_type       offs_B,
    IntersectList & ilist
  ) const {

    IntersectList ilist00, ilist01, ilist02,
                  ilist10, ilist11, ilist12,
                  ilist20, ilist21, ilist22;

    m_C0.intersect_ISO( offs, B.m_C0, offs_B, ilist00 );
    m_C0.intersect_ISO( offs, B.m_C1, offs_B, ilist01 );
    m_C0.intersect_ISO( offs, B.m_C2, offs_B, ilist02 );

    m_C1.intersect_ISO( offs, B.m_C0, offs_B, ilist10 );
    m_C1.intersect_ISO( offs, B.m_C1, offs_B, ilist11 );
    m_C1.intersect_ISO( offs, B.m_C2, offs_B, ilist12 );

    m_C2.intersect_ISO( offs, B.m_C0, offs_B, ilist20 );
    m_C2.intersect_ISO( offs, B.m_C1, offs_B, ilist21 );
    m_C2.intersect_ISO( offs, B.m_C2, offs_B, ilist22 );

    real_type L0  = m_C0.length();
    real_type L1  = L0+m_C1.length();
    real_type LB0 = B.m_C0.length();
    real_type LB1 = LB0+B.m_C1.length();

    ilist.reserve( ilist.size() +

                   ilist00.size() +
                   ilist01.size() +
                   ilist02.size() +

                   ilist10.size() +
                   ilist11.size() +
                   ilist12.size() +

                   ilist20.size() +
                   ilist21.size() +
                   ilist22.size() );

    for ( auto & it : ilist00 ) ilist.push_back( it );
    for ( auto & it : ilist01 ) { it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist02 ) { it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist10 ) { it.first += L0; ilist.push_back( it ); }
    for ( auto & it : ilist11 ) { it.first += L0; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist12 ) { it.first += L0; it.second += LB1; ilist.push_back( it ); }

    for ( auto & it : ilist20 ) { it.first += L1; ilist.push_back( it ); }
    for ( auto & it : ilist21 ) { it.first += L1; it.second += LB0; ilist.push_back( it ); }
    for ( auto & it : ilist22 ) { it.first += L1; it.second += LB1; ilist.push_back( it ); }

 }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    real_type x1, y1, s1, t1, dst1;
    integer res  = m_C0.closest_point_ISO( qx, qy, x,  y,  s,  t,  dst  );
    integer res1 = m_C1.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length()+m_C2.length();
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  void
  Dubins::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      this->intersect( C, ilist );
    } else {
      G2lib::intersect( this, pC, ilist );
    }
  }

  void
  Dubins::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::DUBINS ) {
      Dubins const & C = *static_cast<Dubins const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  Dubins::closest_point_ISO(
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
    integer res  = m_C0.closest_point_ISO( qx, qy, offs, x,  y,  s,  t,  dst  );
    integer res1 = m_C1.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length_ISO(offs);
      t   = t1;
      dst = dst1;
      res = res1;
    }
    res1 = m_C2.closest_point_ISO( qx, qy, offs, x1, y1, s1, t1, dst1 );
    if ( dst1 < dst ) {
      x   = x1;
      y   = y1;
      s   = s1+m_C0.length_ISO(offs)+m_C2.length_ISO(offs);
      t   = t1;
      dst = dst1;
      res = res1;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, Dubins const & bi ) {
    stream
      << "C0\n" << bi.m_C0
      << "C1\n" << bi.m_C1
      << "C2\n" << bi.m_C2
      << "\n";
    return stream;
  }

}

///
/// eof: Dubins.cc
///
