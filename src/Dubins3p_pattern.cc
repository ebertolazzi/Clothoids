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

namespace G2lib {

  void
  Dubins3p::get_sample_angles(
    real_type xi,
    real_type yi,
    real_type thetai,
    real_type xm,
    real_type ym,
    real_type xf,
    real_type yf,
    real_type thetaf,
    real_type k_max,
    real_type tolerance,
    vector<real_type> & angles
  ) const {

    // select angles
    integer NSEG{ integer(std::floor(Utils::m_2pi / m_sample_angle)) };
    angles.clear();
    angles.reserve( 2*NSEG + 12 );
    {
      real_type ang[12];
      integer npts = this->get_range_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, ang );
      if ( npts > 0 ) {
        for ( integer i{0}; i < npts; ++i ) {
          real_type a{ i == 0 ? ang[npts-1]-Utils::m_2pi : ang[i-1] };
          //real_type b{ i == npts ? ang[0]+Utils::m_2pi      : ang[i]   };
          real_type b{ ang[i] };
          real_type delta{ std::min( (b-a)/3, m_sample_angle ) };
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
    integer i{integer(angles.size())};
    for ( --i; i > 0; --i ) {
      if ( std::abs(angles[i]-angles[i-1]) < tolerance ) {
        angles.erase(angles.begin()+i);
      }
    }
  }

  bool
  Dubins3p::build_pattern_search(
    real_type xi,
    real_type yi,
    real_type thetai,
    real_type xm,
    real_type ym,
    real_type xf,
    real_type yf,
    real_type thetaf,
    real_type k_max,
    real_type tolerance,
    bool      use_trichotomy
  ) {

    m_evaluation = 0;

    typedef struct Dubins3p_data {
      Dubins    D0{"temporary Dubins A"};
      Dubins    D1{"temporary Dubins B"};
      real_type thetam{0};
      real_type len{0};

      void
      copy( Dubins3p_data const & rhs ) {
        D0.copy(rhs.D0);
        D1.copy(rhs.D1);
        thetam = rhs.thetam;
        len    = rhs.len;
      }

      //bool
      //compare( Dubins3p_data const & D ) const {
      //  return D0.solution_type() == D.D0.solution_type() &&
      //         D1.solution_type() == D.D1.solution_type();
      //}

    } Dubins3p_data;

    auto eval3p = [this,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max]( Dubins3p_data & D3P ) -> void {
      D3P.D0.build( xi, yi, thetai, xm, ym, D3P.thetam, k_max );
      D3P.D1.build( xm, ym, D3P.thetam, xf, yf, thetaf, k_max );
      D3P.len = D3P.D0.length() + D3P.D1.length();
      ++m_evaluation;
    };

    auto bracketing = [&eval3p]( Dubins3p_data & A, Dubins3p_data & P3, Dubins3p_data & B ) -> void {
      Dubins3p_data P1, P2, P4, P5;
      P2.thetam = (A.thetam+2*P3.thetam)/3;
      eval3p(P2);
      if ( P2.len <= P3.len ) {
        P1.thetam = (2*A.thetam+P3.thetam)/3;
        eval3p(P1);
        if ( P1.len <= P2.len ) { P3.copy(P1); B.copy(P2); }
        else                    { A.copy(P1); B.copy(P3); P3.copy(P2); }
      } else {
        P4.thetam = (B.thetam+2*P3.thetam)/3;
        eval3p(P4);
        if ( P4.len <= P3.len ) {
          P5.thetam = (2*B.thetam+P3.thetam)/3;
          eval3p(P5);
          if ( P5.len <= P4.len ) { A.copy(P4); P3.copy(P5); }
          else                    { A.copy(P3); P3.copy(P4); B.copy(P5); }
        } else {
          A.copy(P2);
          B.copy(P4);
        }
      }
    };

    auto simple_search = [&eval3p]( Dubins3p_data & L, Dubins3p_data & C, Dubins3p_data & R ) -> void {
      Dubins3p_data LL, RR;
      LL.thetam = (C.thetam+L.thetam)/2; eval3p( LL );
      RR.thetam = (C.thetam+R.thetam)/2; eval3p( RR );
      if ( LL.len < RR.len ) {
        if ( LL.len < C.len ) { R.copy(C);  C.copy(LL); }
        else                  { L.copy(LL); R.copy(RR); }
      } else {
        if ( RR.len < C.len ) { L.copy(C);  C.copy(RR); }
        else                  { L.copy(LL); R.copy(RR); }
      }
    };

    vector<real_type> angles;
    this->get_sample_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, tolerance, angles );

    vector<Dubins3p_data> DB(angles.size());
    Dubins3p_data L, C, R, BEST;

    // initialize and find min
    real_type lmin{Utils::Inf<real_type>()};
    integer NSEG{integer(angles.size())};
    C.thetam = angles[NSEG-2]; eval3p( C );
    R.thetam = angles[NSEG-1]; eval3p( R );
    for ( real_type a : angles ) {
      L.copy( C );
      C.copy( R );
      R.thetam = a;
      eval3p( R );
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
