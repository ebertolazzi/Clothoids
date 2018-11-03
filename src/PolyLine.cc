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
#include "ClothoidList.hh"

// Workaround for Visual Studio
#ifdef min
  #undef min
#endif

#ifdef max
  #undef max
#endif

#include <algorithm>

namespace G2lib {

  using std::min;
  using std::max;
  using std::abs;
  using std::cout;
  using std::vector;
  using std::ceil;

  /*\
   |  ____       _       _     _
   | |  _ \ ___ | |_   _| |   (_)_ __   ___
   | | |_) / _ \| | | | | |   | | '_ \ / _ \
   | |  __/ (_) | | |_| | |___| | | | |  __/
   | |_|   \___/|_|\__, |_____|_|_| |_|\___|
   |               |___/
  \*/

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  PolyLine::PolyLine( BaseCurve const & C )
  //: BaseCurve(G2LIB_CLOTHOID)
  : aabb_done(false)
  {
    switch ( C.type() ) {
    case G2LIB_LINE:
      build( *static_cast<LineSegment const *>(&C) );
      break;
    case G2LIB_POLYLINE:
      //copy( *static_cast<PolyLine const *>(&C) );
      //break;
    case G2LIB_CIRCLE:
    case G2LIB_CLOTHOID:
    case G2LIB_BIARC:
    case G2LIB_CLOTHOID_LIST:
      G2LIB_ASSERT( false,
                    "PolyLine constructor cannot convert from: " <<
                    CurveType_name[C.type()] );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::search( real_type s ) const {

    G2LIB_ASSERT( !s0.empty(),
                  "PolyLine::search(" << s << ") empty PolyLine" );

    int_type  npts = int_type(s0.size());
    real_type sl   = s0.front();
    real_type sr   = s0.back();
    G2LIB_ASSERT( s >= sl && s <= sr,
                  "PolyLine::search( " << s <<
                  " ) out of range: [" << sl << ", " << sr << "]" );

    if      ( isegment < 0      ) isegment = 0;
    else if ( isegment > npts-2 ) isegment = npts-2;

    updateInterval( isegment, s, &s0.front(), npts );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LineSegment const &
  PolyLine::getSegment( int_type n ) const {
    G2LIB_ASSERT( !lvec.empty(),
                  "PolyLine::getSegment(...) empty PolyLine" );
    G2LIB_ASSERT( n >= 0 && n < int_type(lvec.size()),
                  "PolyLine::getSegment( " << n <<
                  " ) out of range [0," << lvec.size()-1 << "]" );
    return lvec[size_t(n)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::polygon( real_type x[], real_type y[]) const {
    int_type n = int_type(lvec.size());
    for ( int_type k = 0; k < n; ++k ) {
      x[size_t(k)] = lvec[size_t(k)].xBegin();
      y[size_t(k)] = lvec[size_t(k)].yBegin();
    }
    x[size_t(n)] = lvec[size_t(n-1)].xEnd();
    y[size_t(n)] = lvec[size_t(n-1)].yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {

    G2LIB_ASSERT( !lvec.empty(), "PolyLine::bbox, empty list" );

    if ( aabb_done ) {
      aabb.bbox( xmin, ymin, xmax, ymax );
    } else {
      vector<LineSegment>::const_iterator ic = lvec.begin();
      xmin = xmax = ic->xBegin();
      ymin = ymax = ic->yBegin();
      for ( ++ic; ic != lvec.end(); ++ic ) {
        real_type x = ic->xBegin();
        real_type y = ic->yBegin();
        if      ( x < xmin ) xmin = x;
        else if ( x > xmax ) xmax = x;
        if      ( y < ymin ) ymin = y;
        else if ( y > ymax ) ymax = y;
      }
      --ic;
      real_type x = ic->xEnd();
      real_type y = ic->yEnd();
      if      ( x < xmin ) xmin = x;
      else if ( x > xmax ) xmax = x;
      if      ( y < ymin ) ymin = y;
      else if ( y > ymax ) ymax = y;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build_AABBtree( AABBtree & aabbtree ) const {
    #ifdef G2LIB_USE_CXX11
    vector<shared_ptr<BBox const> > bboxes;
    #else
    vector<BBox const *> bboxes;
    #endif
    bboxes.reserve(lvec.size());
    vector<LineSegment>::const_iterator it;
    int_type ipos = 0;
    for ( it = lvec.begin(); it != lvec.end(); ++it, ++ipos ) {
      real_type xmin, ymin, xmax, ymax;
      it->bbox( xmin, ymin, xmax, ymax );
      #ifdef G2LIB_USE_CXX11
      bboxes.push_back( make_shared<BBox const>(
        xmin, ymin, xmax, ymax, G2LIB_LINE, ipos
      ) );
      #else
      bboxes.push_back( new BBox( xmin, ymin, xmax, ymax, G2LIB_LINE, ipos ) );
      #endif
    }
    aabbtree.build(bboxes);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::init( real_type x0, real_type y0 ) {
    xe = x0;
    ye = y0;
    s0.clear();
    s0.push_back(0);
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( real_type x, real_type y ) {
    LineSegment s;
    s.build_2P( xe, ye, x, y );
    lvec.push_back( s );
    real_type slast = s0.back() + s.length();
    s0.push_back( slast );
    xe = x;
    ye = y;
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( LineSegment const & C ) {
    lvec.push_back( C );
    LineSegment & S = lvec.back();
    S.changeOrigin( xe, ye );
    real_type slast = s0.back() + S.length();
    s0.push_back( slast );
    xe = S.xEnd();
    ye = S.yEnd();
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( CircleArc const & C, real_type tol ) {
    real_type L  = C.length();
    int_type  ns = int_type(ceil( L / C.lenTolerance( tol ) ));
    real_type tx = xe - C.xBegin();
    real_type ty = ye - C.yBegin();
    for ( int_type i = 1; i < ns; ++i ) {
      real_type s = (i*L)/ns;
      push_back( tx + C.X(s), ty + C.Y(s) );
    }
    push_back( tx + C.xEnd(), ty + C.yEnd() );
    xe = tx + C.xEnd();
    ye = ty + C.yEnd();
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( Biarc const & B, real_type tol ) {
    CircleArc const & C0 = B.getC0();
    CircleArc const & C1 = B.getC1();
    real_type L0  = C0.length();
    real_type L1  = C1.length();
    int_type  ns0 = int_type(ceil( L0 / C0.lenTolerance( tol ) ));
    int_type  ns1 = int_type(ceil( L1 / C1.lenTolerance( tol ) ));

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
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidCurve const & C, real_type tol ) {

    real_type L    = C.length();
    real_type absk = max(abs(C.kappaBegin()), abs(C.kappaEnd()));
    real_type tmp  = absk*tol - 1;
    int_type ns = 1;
    if ( tmp > -1 ) ns = int_type( ceil( L*absk/(2*(m_pi-acos(tmp))) ) );

    real_type tx = xe - C.xBegin();
    real_type ty = ye - C.yBegin();
    for ( int_type i = 1; i < ns; ++i ) {
      real_type s = (i*L)/ns;
      push_back( tx + C.X(s), ty + C.Y(s) );
    }

    push_back( tx + C.xEnd(), ty + C.yEnd() );
    xe = tx + C.xEnd();
    ye = ty + C.yEnd();
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidList const & L, real_type tol ) {
    int_type ns = L.numSegment();
    for ( int_type idx = 0; idx < ns; ++idx ) {
      ClothoidCurve const & C = L.get( idx );
      push_back( C, tol );
    }
    // aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( real_type const x[],
                   real_type const y[],
                   int_type npts ) {
    init( x[0], y[0] );
    for ( int_type k = 1; k < npts; ++k )
      push_back( x[k], y[k] );
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
    init( C.xBegin(), C.yBegin() );
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

  void
  PolyLine::reverse() {
    vector<LineSegment>::iterator il;
    for ( il = lvec.begin(); il != lvec.end(); ++il )
      il->reverse();
    std::reverse(lvec.begin(),lvec.end());
    int_type k = 1;
    for ( il = lvec.begin(); il != lvec.end(); ++il, ++k )
      s0[size_t(k)] = s0[size_t(k-1)] + il->length();
    aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::closestPoint(
    real_type   x,
    real_type   y,
    real_type & X,
    real_type & Y,
    real_type & S
  ) const{

    G2LIB_ASSERT( !lvec.empty(), "PolyLine::closestPoint, empty list" );
    vector<LineSegment>::const_iterator ic = lvec.begin();
    vector<real_type>::const_iterator   is = s0.begin();
    real_type DST, T;
    ic->closestPoint( x, y, X, Y, S, T, DST );
    for ( ++ic, ++is; ic != lvec.end(); ++ic, ++is ) {
      real_type X1, Y1, S1, T1, DST1;
      ic->closestPoint( x, y, X1, Y1, S1, T1, DST1 );
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
  PolyLine::intersect(
    PolyLine const    & pl,
    vector<real_type> & ss0,
    vector<real_type> & ss1
  ) const {
    G2LIB_ASSERT( !lvec.empty(),
                  "PolyLine::intersect, empty list" );
    G2LIB_ASSERT( !pl.lvec.empty(),
                  "PolyLine::intersect, empty secondary list" );

#if 1
    build_AABBtree();
    pl.build_AABBtree();
    AABBtree::VecPairPtrBBox intersectionList;
    aabb.intersect( pl.aabb, intersectionList );
    AABBtree::VecPairPtrBBox::const_iterator ip;
    for ( ip = intersectionList.begin(); ip != intersectionList.end(); ++ip ) {
      size_t ipos0 = size_t(ip->first->Ipos());
      size_t ipos1 = size_t(ip->second->Ipos());
      G2LIB_ASSERT( ipos0 < lvec.size(),    "Bad ipos0 = " << ipos0 );
      G2LIB_ASSERT( ipos1 < pl.lvec.size(), "Bad ipos1 = " << ipos1 );
      real_type sss0, sss1;
      bool ok = lvec[ipos0].intersect(pl.lvec[ipos1],sss0,sss1);
      if ( ok ) {
        ss0.push_back(sss0+s0[ipos0]);
        ss1.push_back(sss1+pl.s0[ipos1]);
      }
    }

#else
    ss0.clear();
    ss1.clear();
    vector<LineSegment>::const_iterator ic0 = lvec.begin();
    vector<real_type>::const_iterator   is0 = s0.begin();
    while ( ic0 != lvec.end() ) {
      vector<LineSegment>::const_iterator ic1 = pl.lvec.begin();
      vector<real_type>::const_iterator   is1 = pl.s0.begin();
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
#endif

  }

  bool
  PolyLine::intersect( PolyLine const & pl ) const {
    G2LIB_ASSERT( !lvec.empty(),
                  "PolyLine::intersect, empty list" );
    G2LIB_ASSERT( !pl.lvec.empty(),
                  "PolyLine::intersect, empty secondary list" );
    vector<LineSegment>::const_iterator ic0 = lvec.begin();
    vector<real_type>::const_iterator   is0 = s0.begin();
    while ( ic0 != lvec.end() ) {
      vector<LineSegment>::const_iterator ic1 = pl.lvec.begin();
      vector<real_type>::const_iterator   is1 = pl.s0.begin();
      while ( ic1 != pl.lvec.end() ) {
        if ( ic0->collision( *ic1 ) ) return true;
        ++ic1;
        ++is1;
      }
      ++ic0;
      ++is0;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, PolyLine const & P ) {
    stream <<   "nseg   = " << P.numSegment()
           << "\nxBegin = " << P.xBegin()
           << "\nybegin = " << P.yBegin()
           << "\nxEnd   = " << P.xEnd()
           << "\nyEnd   = " << P.yEnd()
           << "\nlength = " << P.length()
           << "\n";
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

// EOF: PolyLine.cc
