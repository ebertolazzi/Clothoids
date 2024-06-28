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

    // select angles
    integer NSEG{ integer(std::floor(Utils::m_2pi / m_sample_angle)) };
    vector<real_type> angles;
    angles.reserve( NSEG + 12 );
    for ( integer i{0}; i < NSEG; ++i ) angles.push_back( i*m_sample_angle );
    {
      real_type ang[12];
      integer npts = this->get_range_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, ang );
      for ( integer i{0}; i < npts; ++i ) {
        // Find the position to insert the new value using lower_bound
        auto it = std::lower_bound(angles.begin(), angles.end(), ang[i]);
        angles.insert(it,ang[i]);
      }
    }

    vector<Dubins3p_data> DB(angles.size());
    Dubins3p_data L, C, R;

    // initialize and find min
    integer   imin{0};
    real_type lmin{Utils::Inf<real_type>()};
    NSEG = 0;
    for ( real_type a : angles ) {
      Dubins3p_data & db{ DB[NSEG] };
      db.thetam = a;
      eval3p( db );
      if ( db.len < lmin ) { imin = NSEG; lmin = db.len; }
      ++NSEG;
    }

    // select interval
    L.copy( DB[(imin+NSEG-1)%NSEG] );
    C.copy( DB[imin]               );
    R.copy( DB[(imin+1)%NSEG]      );

    // make angles monotone increasing
    if ( imin == 0      ) L.thetam -= Utils::m_2pi;
    if ( imin == NSEG-1 ) R.thetam += Utils::m_2pi;

    if ( use_trichotomy ) {
      while ( R.thetam-L.thetam > tolerance && m_evaluation < m_max_evaluation )
        bracketing( L, C, R );
    } else {
      while ( R.thetam-L.thetam > tolerance && m_evaluation < m_max_evaluation )
        simple_search( L, C, R );
    }

    m_Dubins0.copy(C.D0);
    m_Dubins1.copy(C.D1);

    return true;
  }

}

///
/// eof: Dubins3p_pattern.cc
///
