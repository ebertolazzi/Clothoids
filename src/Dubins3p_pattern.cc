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
/// file: Dubins3p_pattern.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_AlgoBracket.hh"

namespace G2lib {

  void
  Dubins3p::get_sample_angles(
    real_type const xi,
    real_type const yi,
    real_type const thetai,
    real_type const xm,
    real_type const ym,
    real_type const xf,
    real_type const yf,
    real_type const thetaf,
    real_type const k_max,
    real_type const tolerance,
    vector<real_type> & angles
  ) const {

    // select angles
    integer const NSEG{ static_cast<integer>(std::floor(Utils::m_2pi / m_sample_angle)) };
    angles.clear();
    angles.reserve( 2*NSEG + 12 );
    {
      real_type ang[12];
      integer npts{ this->get_range_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, ang ) };
      if ( npts > 0 ) {
        for ( integer i{0}; i < npts; ++i ) {
          real_type a{ i == 0 ? ang[npts-1]-Utils::m_2pi : ang[i-1] };
          //real_type b{ i == npts ? ang[0]+Utils::m_2pi      : ang[i]   };
          real_type const b{ ang[i] };
          real_type const delta{ std::min( (b-a)/2.99999, m_sample_angle ) };
          while ( a < b ) {
            real_type aa{ a };
            if      ( aa < 0            ) aa += Utils::m_2pi;
            else if ( aa > Utils::m_2pi ) aa -= Utils::m_2pi;
            angles.push_back( aa );
            a += delta;
          }
        }
      } else {
        real_type a{ 0 };
        while ( a < Utils::m_2pi ) { angles.push_back( a ); a += m_sample_angle; }
      }
    }
    std::sort( angles.begin(), angles.end() );
    // remove duplicates
    integer i{static_cast<integer>(angles.size())};
    for ( --i; i > 0; --i ) {
      if ( std::abs(angles[i]-angles[i-1]) < tolerance ) {
        angles.erase(angles.begin()+i);
      }
    }
  }

  bool
  Dubins3p::build_pattern_search(
    real_type       xi,
    real_type       yi,
    real_type       thetai,
    real_type       xm,
    real_type       ym,
    real_type       xf,
    real_type       yf,
    real_type       thetaf,
    real_type       k_max,
    real_type const tolerance,
    bool      const use_trichotomy,
    bool      const use_bracket
  ) {

    m_evaluation = 0;

    typedef struct Dubins3p_data {
      Dubins    D0{"temporary Dubins A"};
      Dubins    D1{"temporary Dubins B"};
      real_type thetam{0};
      real_type len{0};
      real_type len_D{0};

      void
      copy( Dubins3p_data const & rhs ) {
        D0.copy(rhs.D0);
        D1.copy(rhs.D1);
        thetam = rhs.thetam;
        len    = rhs.len;
        len_D  = rhs.len_D;
      }

    } Dubins3p_data;

    auto check_change_sign = []( Dubins3p_data const & A, Dubins3p_data const & B ) -> bool {
      bool ok{ A.D0.icode() == B.D0.icode() &&
               A.D1.icode() == B.D1.icode() };
      if ( ok ) ok = A.len_D * B.len_D <= 0;
      return ok;
    };

    auto eval3p = [this,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max,tolerance,use_trichotomy,use_bracket]( Dubins3p_data & D3P, real_type const thetam ) -> void {
      D3P.thetam = thetam;
      bool ok{ D3P.D0.build( xi, yi, thetai, xm, ym, thetam, k_max ) };
      if ( ok ) ok = D3P.D1.build( xm, ym, thetam, xf, yf, thetaf,     k_max );
      UTILS_ASSERT(
        ok,
        "Dubins3p::build_pattern_search(\n"
        "  xi             = {}\n"
        "  yi             = {}\n"
        "  thetai         = {}\n"
        "  xm             = {}\n"
        "  ym             = {}\n"
        "  xf             = {}\n"
        "  yf             = {}\n"
        "  thetaf         = {}\n"
        "  k_max          = {}\n"
        "  tolerance      = {}\n"
        "  use_trichotomy = {}\n"
        "  use_bracket    = {}\n"
        ") failed\n",
        xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, tolerance, use_trichotomy, use_bracket
      );
      D3P.len   = D3P.D0.length()       + D3P.D1.length();
      D3P.len_D = D3P.D0.length_Dbeta() + D3P.D1.length_Dalpha();
      ++m_evaluation;
    };

    auto eval_for_bracket = [this,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max,tolerance,use_trichotomy,use_bracket]( real_type const theta ) -> real_type {
      Dubins D0{"temporary Dubins A"};
      Dubins D1{"temporary Dubins B"};
      bool ok{ D0.build( xi, yi, thetai, xm, ym, theta,  k_max ) };
      if ( ok ) ok = D1.build( xm, ym, theta,  xf, yf, thetaf, k_max );
      UTILS_ASSERT(
        ok,
        "Dubins3p::build_pattern_search(\n"
        "  xi             = {}\n"
        "  yi             = {}\n"
        "  thetai         = {}\n"
        "  xm             = {}\n"
        "  ym             = {}\n"
        "  xf             = {}\n"
        "  yf             = {}\n"
        "  thetaf         = {}\n"
        "  k_max          = {}\n"
        "  tolerance      = {}\n"
        "  use_trichotomy = {}\n"
        "  use_bracket    = {}\n"
        ") failed\n",
        xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, tolerance, use_trichotomy, use_bracket
      );
      ++m_evaluation;
      return D0.length_Dbeta() + D1.length_Dalpha();
    };

    auto do_bracket = [&check_change_sign,&eval3p,&eval_for_bracket]( Dubins3p_data & L, Dubins3p_data & C, Dubins3p_data & R ) -> bool {
      Utils::AlgoBracket<real_type> algoBracket;
      real_type theta{0};
      bool ok{ false };
      if ( check_change_sign(L,C) ) {
        theta = algoBracket.eval3( L.thetam, C.thetam, L.len_D, C.len_D, eval_for_bracket );
        ok    = algoBracket.converged();
      } else if ( check_change_sign(C,R) ) {
        theta = algoBracket.eval3( C.thetam, R.thetam, C.len_D, R.len_D, eval_for_bracket );
        ok    = algoBracket.converged();
      }
      // if ok collapse to one point
      if ( ok ) { eval3p( C, theta ); L.copy(C); R.copy(C); }
      return ok;
    };

    auto bracketing = [&eval3p,&do_bracket,&use_bracket]( Dubins3p_data & A, Dubins3p_data & P3, Dubins3p_data & B ) -> void {
      bool ok{ use_bracket };
      if ( ok ) ok = do_bracket( A, P3, B );
      if ( !ok ) {
        Dubins3p_data P2;
        eval3p( P2, (A.thetam+2*P3.thetam)/3 );
        if ( P2.len <= P3.len ) {
          Dubins3p_data P1;
          eval3p( P1, (2*A.thetam+P3.thetam)/3 );
          if ( P1.len <= P2.len ) { P3.copy(P1); B.copy(P2); }
          else                    { A.copy(P1);  B.copy(P3); P3.copy(P2); }
        } else {
          Dubins3p_data P4;
          eval3p( P4, (B.thetam+2*P3.thetam)/3 );
          if ( P4.len <= P3.len ) {
            Dubins3p_data P5;
            eval3p( P5, (2*B.thetam+P3.thetam)/3 );
            if ( P5.len <= P4.len ) { A.copy(P4); P3.copy(P5); }
            else                    { A.copy(P3); P3.copy(P4); B.copy(P5); }
          } else {
            A.copy(P2);
            B.copy(P4);
          }
        }
      }
    };

    auto simple_search = [&eval3p,&do_bracket,use_bracket]( Dubins3p_data & L, Dubins3p_data & C, Dubins3p_data & R ) -> void {
      bool ok{ use_bracket };
      if ( ok ) ok = do_bracket( L, C, R );
      if ( !ok ) {
        Dubins3p_data LL, RR;
        eval3p( LL, (C.thetam+L.thetam)/2 );
        eval3p( RR, (C.thetam+R.thetam)/2 );
        if ( LL.len < RR.len ) {
          if ( LL.len < C.len ) { R.copy(C);  C.copy(LL); }
          else                  { L.copy(LL); R.copy(RR); }
        } else {
          if ( RR.len < C.len ) { L.copy(C);  C.copy(RR); }
          else                  { L.copy(LL); R.copy(RR); }
        }
      }
    };

    vector<real_type> angles;
    this->get_sample_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, tolerance, angles );

    vector<Dubins3p_data> DB(angles.size());
    Dubins3p_data L, C, R, BEST;

    // initialize and find min
    real_type lmin{Utils::Inf<real_type>()};
    integer const NSEG{static_cast<integer>(angles.size())};
    eval3p( C, angles[NSEG-2] );
    eval3p( R, angles[NSEG-1] );

    for ( real_type const a : angles ) {
      L.copy( C );
      C.copy( R );
      eval3p( R, a );

      if ( C.len <= L.len && C.len <= R.len ) {
        if ( C.len < 1.5*lmin ) {
          Dubins3p_data LL, CC, RR;
          LL.copy(L);
          CC.copy(C);
          RR.copy(R);
          // make angles monotone increasing
          if ( LL.thetam > CC.thetam ) LL.thetam -= Utils::m_2pi;
          if ( RR.thetam < CC.thetam ) RR.thetam += Utils::m_2pi;
          if ( use_trichotomy ) {
            while ( RR.thetam-LL.thetam > tolerance && m_evaluation < m_max_evaluation )
              bracketing( LL, CC, RR );
          } else {
            while ( RR.thetam-LL.thetam > tolerance && m_evaluation < m_max_evaluation )
              simple_search( LL, CC, RR );
          }
          if ( CC.len < lmin ) {
            lmin = CC.len;
            BEST.copy(CC);
          }
        }
      }
    }

    m_Dubins0.copy(BEST.D0);
    m_Dubins1.copy(BEST.D1);

    return true;
  }

}

///
/// eof: Dubins3p_pattern.cc
///
