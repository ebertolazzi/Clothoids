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

#include "PolyLine.hh"
#include "Line.hh"
#include "Circle.hh"
#include "Biarc.hh"
#include "Clothoid.hh"

namespace G2lib {

  /*\
   |  ____       _       _     _
   | |  _ \ ___ | |_   _| |   (_)_ __   ___
   | | |_) / _ \| | | | | |   | | '_ \ / _ \
   | |  __/ (_) | | |_| | |___| | | | |  __/
   | |_|   \___/|_|\__, |_____|_|_| |_|\___|
   |               |___/
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::init( real_type x0, real_type y0 ) {
    xe = x0;
    ye = y0;
    s0.clear();
    s0.push_back(0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( real_type x, real_type y ) {
    LineSegment s;
    s.build_2P( xe, ye, x, y );
    lvec.push_back( s );
    s0.push_back( s0.back() + s.length() );
    xe = x;
    ye = y;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( LineSegment const & C ) {
    lvec.push_back( C );
    LineSegment & S = lvec.back();
    S.changeOrigin( xe, ye );
    s0.push_back( s0.back() + S.length() );
    xe = S.xEnd();
    ye = S.yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( CircleArc const & C, real_type tol ) {
    real_type L  = C.length();
    int_type  ns = std::ceil( L / C.lenTolerance( tol ) );
    real_type tx = xe - C.xBegin();
    real_type ty = ye - C.yBegin();
    for ( int_type i = 1; i < ns; ++i ) {
      real_type s = (i*L)/ns;
      push_back( tx + C.X(s), ty + C.Y(s) );
    }
    push_back( tx + C.xEnd(), ty + C.yEnd() );
    xe = tx + C.xEnd();
    ye = ty + C.yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( Biarc const & B, real_type tol ) {
    CircleArc const & C0 = B.getC0();
    CircleArc const & C1 = B.getC1();
    real_type L0  = C0.length();
    real_type L1  = C1.length();
    int_type  ns0 = std::ceil( L0 / C0.lenTolerance( tol ) );
    int_type  ns1 = std::ceil( L1 / C1.lenTolerance( tol ) );

    real_type tx = xe - C0.xBegin();
    real_type ty = ye - C0.yBegin();

    for ( int_type i = 1; i < ns0; ++i ) {
      real_type s = (i*L0)/ns0;
      push_back( tx + C0.X(s), ty + C0.Y(s) );
    }
    push_back( tx + C1.xBegin(), ty + C1.yBegin() );
    for ( int_type i = 1; i < ns1; ++i ) {
      real_type s = (i*L1)/ns1;
      push_back( tx + C1.X(s), ty + C1.Y(s) );
    }
    push_back( tx + C1.xEnd(), ty + C1.yEnd() );
    xe = tx + C1.xEnd();
    ye = ty + C1.yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidCurve const & C, real_type tol ) {
    real_type tx = xe - C.xBegin();
    real_type ty = ye - C.yBegin();
    real_type s  = 0;

    real_type ds, absk, tmp;
    while ( true ) {
      absk = std::abs(C.kappa(s));
      tmp  = absk*tol;
      ds   = tmp > 0 ? 2*(m_pi-acos(tmp-1))/absk : C.length();
      absk += ds*std::abs(C.dkappa());
      tmp  = absk*tol;
      ds   = tmp > 0 ? 2*(m_pi-acos(tmp-1))/absk : C.length();

      s += ds;
      if ( s >= C.length() ) {
        push_back( tx + C.xEnd(), ty + C.yEnd() );
        break;
      }
      push_back( tx + C.X(s), ty + C.Y(s) );
    }
    xe = tx + C.xEnd();
    ye = ty + C.yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidList const & L, real_type tol ) {
    int_type ns = L.numSegment();
    for ( int_type idx = 0; idx < ns; ++idx ) {
      ClothoidCurve const & C = L.get( idx );
      push_back( C, tol );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( LineSegment const & C ) {
    init( C.xBegin(), C.yBegin() );
    push_back( C.xEnd(), C.yEnd() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( CircleArc const & C, real_type tol ) {
    init( C.xBegin(), C.yBegin() );
    push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( Biarc const & C, real_type tol ) {
    init( C.xBegin0(), C.yBegin0() );
    push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( ClothoidCurve const & C, real_type tol ) {
    init( C.xBegin(), C.yBegin() );
    push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( ClothoidList const & L, real_type tol ) {
    init( L.xBegin(), L.yBegin() );
    push_back( L, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::closestPoint( real_type   x,
                          real_type   y,
                          real_type & X,
                          real_type & Y,
                          real_type & S ) const{

    G2LIB_ASSERT( !lvec.empty(), "PolyLine::closestPoint, empty list" );
    std::vector<LineSegment>::const_iterator ic = lvec.begin();
    std::vector<real_type>::const_iterator   is = s0.begin();
    real_type DST = ic->closestPoint( x, y, X, Y, S );
    for ( ++ic, ++is; ic != lvec.end(); ++ic, ++is ) {
      real_type X1, Y1, S1;
      real_type DST1 = ic->closestPoint( x, y, X1, Y1, S1 );
      if ( DST1 < DST ) {
        DST = DST1;
        X   = X1;
        Y   = Y1;
        S   = *is + S1;
      }
    }
    return DST;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::intersect( PolyLine const         & pl,
                       std::vector<real_type> & ss0,
                       std::vector<real_type> & ss1 ) const {
    G2LIB_ASSERT( !lvec.empty(), "PolyLine::intersect, empty list" );
    G2LIB_ASSERT( !pl.lvec.empty(), "PolyLine::intersect, empty secondary list" );
    ss0.clear();
    ss1.clear();
    std::vector<LineSegment>::const_iterator ic0 = lvec.begin();
    std::vector<real_type>::const_iterator   is0 = s0.begin();
    while ( ic0 != lvec.end() ) {
      std::vector<LineSegment>::const_iterator ic1 = pl.lvec.begin();
      std::vector<real_type>::const_iterator   is1 = pl.s0.begin();
      while ( ic1 != pl.lvec.end() ) {
        real_type a0, a1;
        bool ok = ic0->intersect( *ic1, a0, a1 );
        if ( ok ) {
          ss0.push_back( (*is0) + a0 );
          ss1.push_back( (*is1) + a1 );
        }
        ++ic1;
        ++is1;
      }
      ++ic0;
      ++is0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

// EOF: PolyLine.cc
