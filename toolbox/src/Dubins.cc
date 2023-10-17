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

}

///
/// eof: Dubins.cc
///
