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

  #define PRINT_DEBUG(S0,S1,S2) \
  fmt::print( "LENGTH = {}\n", l0+l1+l2 ); \
  m_C0.build( x0,           y0,           theta0,           S0*k_max, l0/k_max ); \
  m_C1.build( m_C0.x_end(), m_C0.y_end(), m_C0.theta_end(), S1*k_max, l1/k_max ); \
  m_C2.build( m_C1.x_end(), m_C1.y_end(), m_C1.theta_end(), S2*k_max, l2/k_max ); \
  fmt::print( "C0\n{}\n", m_C0 ); \
  fmt::print( "C1\n{}\n", m_C1 ); \
  fmt::print( "C2\n{}\n", m_C2 )

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
    real_type d     = k_max * hypot( dx, dy );
    real_type d2    = d*d;

    // metto angolo in -Pi, Pi
    minus_pi_pi( alpha );
    minus_pi_pi( beta );

    // common constants
    real_type Ca  = cos(alpha);
    real_type Sa  = sin(alpha);
    real_type Cad = Ca*d;
    real_type Sad = Sa*d;
    real_type Cb  = cos(beta);
    real_type Sb  = sin(beta);
    real_type Cbd = Cb*d;
    real_type Sbd = Sb*d;
    real_type Cab = cos(alpha-beta);
    real_type Sab = sin(alpha-beta);

    real_type CF1  = 1-Cab-Sad;
    real_type CF2  = 2*(Sab-Cad);
    real_type CF3  = 3+Cab+Sad;

    real_type CF4  = Cab-Sad-1;
    real_type CF5  = 2*(Sab+Cad);
    real_type CF6  = Sad-Cab-3;

    real_type CF7  = d2+2*(Cab-Sad-Sbd-1);
    real_type CF8  = -8*(Sab+Cad);
    real_type CF9  = d2+6*(1-Cab+Sad)-2*Sbd;

    real_type CF10 = d2+2*(Cab+Sad+Sbd-1);
    real_type CF11 = 8*(Sab-Cad);
    real_type CF12 = d*d+6*(1-Cab-Sad)+2*Sbd;

    real_type min_len = Utils::Inf<real_type>();
    real_type L0{0}, L1{0}, L2{0}, s0{0}, s1{0}, s2{0};

    { // caso LSL
      real_type t1 = atan2(Cb-Ca, d-Sb+Sa);
      real_type l0 = t1-alpha;
      real_type l2 = beta-t1;
      if ( l0 < 0 ) l0 += Utils::m_2pi;
      if ( l2 < 0 ) l2 += Utils::m_2pi;
      //real_type l1 = d*cos(alpha+l0)-sin(l0)-sin(l2);
      real_type l1 = (Cab-Sad-1)*sin(l0) + (Cad+Sab)*cos(l0);
      if ( l1 >= 0 ) {
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(1,0,1);
        if ( len < min_len ) {
          min_len = len;
          L0 = l0; L1 = l1; L2 = l2;
          s0 = 1;  s1 = 0;  s2 = 1;
        }
      }
    }

    { // caso RSR
      real_type t1 = atan2(Ca-Cb, d+Sb-Sa);
      real_type l0 = alpha-t1;
      real_type l2 = t1-beta;
      if ( l0 < 0 ) l0 += Utils::m_2pi;
      if ( l2 < 0 ) l2 += Utils::m_2pi;
      //real_type l1 = d*cos(alpha-l0)-sin(l2)-sin(l0);
      real_type l1 = (Sad+Cab-1)*sin(l0) + (Cad-Sab)*cos(l0);
      if ( l1 >= 0 ) {
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(-1,0,-1);
        if ( len < min_len ) {
          min_len = len;
          L0 = l0; L1 = l1; L2 = l2;
          s0 = -1; s1 =  0; s2 = -1;
        }
      }
    }

    { // caso LSR
      Quadratic Q( CF3, CF2, CF1 );
      integer   nr;
      real_type Z[2];
      nr = Q.get_real_roots( Z );
      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type z = Z[ir];
        real_type l0 = 2*atan(z);
        if ( l0 < 0 ) l0 += Utils::m_2pi;
        real_type S0 = sin(l0);
        real_type C0 = cos(l0);
        real_type l1 = (Cad-Sab)*C0 - (Sad+Cab+1)*S0;
        if ( l1 < 0 ) continue;
        real_type l2 = alpha-beta+l0;
        into_0_2pi( l2 );
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(1,0,-1);
        if ( len >= min_len ) continue;
        min_len = len;
        L0 = l0; L1 = l1; L2 = l2;
        s0 =  1; s1 =  0; s2 = -1;
      }
    }

    { // caso RSL
      Quadratic Q( CF6, CF5, CF4 );
      integer   nr;
      real_type Z[2];
      nr = Q.get_real_roots( Z );
      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type z  = Z[ir];
        real_type l0 = 2*atan(z);
        if ( l0 < 0 ) l0 += Utils::m_2pi;
        real_type S0 = sin(l0);
        real_type C0 = cos(l0);
        real_type l1 = (Sad-Cab-1)*S0 + (Cad+Sab)*C0;
        if ( l1 < 0 ) continue;
        real_type l2 = beta-alpha+l0;
        into_0_2pi( l2 );
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(-1,0,1);
        if ( len >= min_len ) continue;
        min_len = len;
        L0 = l0; L1 = l1; L2 = l2;
        s0 = -1; s1 =  0; s2 =  1;
      }
    }

    { // caso LRL
      Quadratic Q( CF9, CF8, CF7 );
      integer   nr;
      real_type Z[2];
      nr = Q.get_real_roots( Z );
      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type z = Z[ir];
        real_type l0 = 2*atan(z);
        if ( l0 < 0 ) l0 += Utils::m_2pi;
        real_type S0 = sin(l0);
        real_type C0 = cos(l0);
        real_type l1 = atan2( (Cab-Sad-1)*S0 + (Cad + Sab)*C0, 2 - (Cad+Sab)*S0 - (Sad-Cab+1)*C0 );
        if ( l1 < 0 ) l1 += Utils::m_2pi; // @@@@@@@@@@@@@
        real_type l2 = beta-alpha-l0+l1;
        into_0_2pi( l2 );
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(1,-1,1);
        if ( len >= min_len ) continue;
        min_len = len;
        L0 = l0; L1 = l1; L2 = l2;
        s0 =  1; s1 = -1; s2 =  1;
      }
    }

    { // caso RLR
      Quadratic Q( CF12, CF11, CF10 );
      integer   nr;
      real_type Z[2];
      nr = Q.get_real_roots( Z );
      for ( integer ir = 0; ir < nr; ++ir ) {
        real_type z = Z[ir];
        real_type l0 = 2*atan(z);
        if ( l0 < 0 ) l0 += Utils::m_2pi;
        real_type S0 = sin(l0);
        real_type C0 = cos(l0);
        real_type l1 = atan2( (Sad+Cab-1)*S0 + (Cad-Sab)*C0, 2 - (Cad-Sab)*S0 + (Sad+Cab-1)*C0 );
        if ( l1 < 0 ) l1 += Utils::m_2pi; // @@@@@@@@@@@@@
        real_type l2 = alpha-beta-l0+l1;
        into_0_2pi( l2 );
        real_type len = l0+l1+l2;
        //PRINT_DEBUG(-1,1,-1);
        if ( len >= min_len ) continue;
        min_len = len;
        L0 = l0; L1 = l1; L2 = l2;
        s0 = -1; s1 =  1; s2 = -1;
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
