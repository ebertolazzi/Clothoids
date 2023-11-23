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

#include "Clothoids.hh"

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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using LS_dist_type = vector<LineSegment>::difference_type;
  #endif

  void PolyLine::build( CircleArc const & )     { UTILS_ERROR("can convert from CircleArc to PolyLine\n"); }
  void PolyLine::build( ClothoidCurve const & ) { UTILS_ERROR("can convert from Clothoid to PolyLine\n"); }
  void PolyLine::build( Biarc const & )         { UTILS_ERROR("can convert from Biarc to PolyLine\n"); }
  void PolyLine::build( BiarcList const & )     { UTILS_ERROR("can convert from BiarcList to PolyLine\n"); }
  void PolyLine::build( PolyLine const & PL )   { *this = PL; }
  void PolyLine::build( ClothoidList const & )  { UTILS_ERROR("can convert from ClothoidList to PolyLine\n"); }

  /*\
   |  ____       _       _     _
   | |  _ \ ___ | |_   _| |   (_)_ __   ___
   | | |_) / _ \| | | | | |   | | '_ \ / _ \
   | |  __/ (_) | | |_| | |___| | | | |  __/
   | |_|   \___/|_|\__, |_____|_|_| |_|\___|
   |               |___/
  \*/

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  PolyLine::PolyLine( BaseCurve const * pC ) {

    G2LIB_DEBUG_MESSAGE( "PolyLine convert: {}\n", pC->type_name() );

    this->resetLastInterval();
    switch ( pC->type() ) {
    case CurveType::LINE:
      G2LIB_DEBUG_MESSAGE( "to -> LineSegment\n" );
      this->build( *static_cast<LineSegment const *>(pC) );
      break;
    case CurveType::POLYLINE:
      G2LIB_DEBUG_MESSAGE( "to -> PolyLine\n" );
      this->copy( *static_cast<PolyLine const *>(pC) );
      break;
    default:
      UTILS_ERROR(
        "PolyLine constructor cannot convert from: {}\n",
        pC->type_name()
      );
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PolyLine::PolyLine( LineSegment const & LS ) {
    this->resetLastInterval();
    this->init( LS.x_begin(), LS.y_begin() );
    this->push_back( LS );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PolyLine::PolyLine( CircleArc const & C, real_type tol ) {
    this->resetLastInterval();
    this->init( C.x_begin(), C.y_begin() );
    this->push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PolyLine::PolyLine( Biarc const & B, real_type tol ) {
    this->resetLastInterval();
    this->init( B.x_begin(), B.y_begin() );
    this->push_back( B, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PolyLine::PolyLine( ClothoidCurve const & C, real_type tol ) {
    this->resetLastInterval();
    this->init( C.x_begin(), C.y_begin() );
    this->push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PolyLine::PolyLine( ClothoidList const & PL, real_type tol ) {
    this->resetLastInterval();
    this->init( PL.x_begin(), PL.y_begin() );
    this->push_back( PL, tol );
  }

  real_type
  PolyLine::length_ISO( real_type ) const {
    UTILS_ERROR0( "PolyLine::length( offs ) not available!\n" );
    //return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  PolyLine::find_at_s( real_type & s ) const {
    #ifdef CLOTHOIDS_USE_THREADS
    bool ok;
    integer & lastInterval = *m_lastInterval.search( std::this_thread::get_id(), ok );
    #else
    integer & lastInterval = m_lastInterval;
    #endif
    Utils::searchInterval<integer,real_type>(
      static_cast<integer>(m_s0.size()),
      m_s0.data(), s, lastInterval, false, true
    );
    return lastInterval;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::init() {
    m_s0.clear();
    m_polylineList.clear();
    this->resetLastInterval();
    m_aabb_done = false;
    this->resetLastInterval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::copy( PolyLine const & PL ) {
    this->init();
    m_polylineList.reserve( PL.m_polylineList.size() );
    std::copy(
      PL.m_polylineList.begin(),
      PL.m_polylineList.end(),
      back_inserter(m_polylineList)
    );
    m_s0.reserve( PL.m_s0.size() );
    std::copy( PL.m_s0.begin(), PL.m_s0.end(), back_inserter(m_s0) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LineSegment const &
  PolyLine::getSegment( integer n ) const {
    UTILS_ASSERT0(
      !m_polylineList.empty(),
      "PolyLine::getSegment(...) empty PolyLine\n"
    );
    UTILS_ASSERT(
      n >= 0 && n < integer(m_polylineList.size()),
      "PolyLine::getSegment( {} ) out of range [0,{}]\n",
      n, m_polylineList.size()-1
    );
    return m_polylineList[size_t(n)];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::polygon( real_type * x, real_type * y) const {
    integer n = integer(m_polylineList.size());
    for ( size_t k = 0; k < size_t(n); ++k ) {
      x[k] = m_polylineList[k].x_begin();
      y[k] = m_polylineList[k].y_begin();
    }
    x[size_t(n)] = m_polylineList[size_t(n-1)].x_end();
    y[size_t(n)] = m_polylineList[size_t(n-1)].y_end();
  }

  /*\
   |   _     _
   |  | |__ | |__   _____  __
   |  | '_ \| '_ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |_.__/|_.__/ \___/_/\_\
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::bbox(
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {

    UTILS_ASSERT0( !m_polylineList.empty(), "PolyLine::bbox, empty list\n" );

    if ( m_aabb_done ) {
      real_type bb_min[2], bb_max[2];
      m_aabb_tree.get_root_bbox( bb_min, bb_max );
      xmin = bb_min[0]; ymin = bb_min[1];
      xmax = bb_max[0]; ymax = bb_max[1];
    } else {
      vector<LineSegment>::const_iterator ic = m_polylineList.begin();
      xmin = xmax = ic->x_begin();
      ymin = ymax = ic->y_begin();
      for ( ++ic; ic != m_polylineList.end(); ++ic ) {
        real_type x = ic->x_begin();
        real_type y = ic->y_begin();
        if      ( x < xmin ) xmin = x;
        else if ( x > xmax ) xmax = x;
        if      ( y < ymin ) ymin = y;
        else if ( y > ymax ) ymax = y;
      }
      --ic;
      real_type x = ic->x_end();
      real_type y = ic->y_end();
      if      ( x < xmin ) xmin = x;
      else if ( x > xmax ) xmax = x;
      if      ( y < ymin ) ymin = y;
      else if ( y > ymax ) ymax = y;
    }
  }

  /*\
   |  _    _   _____    _                _
   | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
   | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
   | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
   |                               |___/
  \*/

  void
  PolyLine::bbTriangles(
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<LineSegment>::const_iterator ic = m_polylineList.begin();
    for ( integer ipos = icurve; ic != m_polylineList.end(); ++ic, ++ipos )
      ic->bbTriangles( tvec, max_angle, max_size, ipos );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  PolyLine::bbTriangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<LineSegment>::const_iterator ic = m_polylineList.begin();
    for ( integer ipos = icurve; ic != m_polylineList.end(); ++ic, ++ipos )
      ic->bbTriangles_ISO( offs, tvec, max_angle, max_size, ipos );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::theta( real_type s ) const {
    integer idx = this->find_at_s( s );
    return m_polylineList[idx].m_theta0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::theta_D( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::theta_DD( real_type ) const
  { return 0; }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  PolyLine::theta_DDD( real_type ) const
  { return 0; }

  /*\
   |  _                        __
   | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
   | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
   | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
   |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::scale( real_type sfactor ) {
    vector<LineSegment>::iterator ic = m_polylineList.begin();
    real_type newx0 = ic->x_begin();
    real_type newy0 = ic->y_begin();
    m_s0[0] = 0;
    for ( size_t k=0; ic != m_polylineList.end(); ++ic, ++k ) {
      ic->scale( sfactor );
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::reverse() {
    std::reverse( m_polylineList.begin(), m_polylineList.end() );
    vector<LineSegment>::iterator ic = m_polylineList.begin();
    ic->reverse();
    real_type newx0 = ic->x_end();
    real_type newy0 = ic->y_end();
    m_s0[0] = 0;
    m_s0[1] = ic->length();
    size_t k = 1;
    for ( ++ic; ic != m_polylineList.end(); ++ic, ++k ) {
      ic->reverse();
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::change_origin( real_type newx0, real_type newy0 ) {
    for ( auto & L : m_polylineList ) {
      L.change_origin( newx0, newy0 );
      newx0 = L.x_end();
      newy0 = L.y_end();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::trim( real_type s_begin, real_type s_end ) {
    UTILS_ASSERT(
      s_begin >= m_s0.front() && s_end <= m_s0.back() && s_end > s_begin,
      "void::trim( s_begin={}, s_end={} ) bad range, must be in [{},{}]\n",
      s_begin, s_end, m_s0.front(), m_s0.back()
    );

    size_t i_begin = size_t(find_at_s(s_begin));
    size_t i_end   = size_t(find_at_s(s_end));
    m_polylineList[i_begin].trim( s_begin-m_s0[i_begin], m_s0[i_begin+1] );
    m_polylineList[i_end].trim( m_s0[i_end], s_end-m_s0[i_end] );
    m_polylineList.erase( m_polylineList.begin()+LS_dist_type(i_end+1), m_polylineList.end() );
    m_polylineList.erase( m_polylineList.begin(), m_polylineList.begin()+LS_dist_type(i_begin) );
    vector<LineSegment>::iterator ic = m_polylineList.begin();
    m_s0[0] = 0;
    size_t k{0};
    for (; ic != m_polylineList.end(); ++ic, ++k )
      m_s0[k+1] = m_s0[k] + ic->length();
    this->resetLastInterval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::trim( real_type s_begin, real_type s_end, PolyLine & newPL ) const {

    newPL.init();

    if ( m_polylineList.empty() ) return;

    // put in range
    real_type L = this->length();
    while ( s_begin > L ) s_begin -= L;
    while ( s_begin < 0 ) s_begin += L;
    while ( s_end   > L ) s_end   -= L;
    while ( s_end   < 0 ) s_end   += L;

    integer n_seg   = integer( m_polylineList.size() );
    integer i_begin = find_at_s( s_begin );
    integer i_end   = find_at_s( s_end );

    if ( s_begin < s_end ) {
      // get initial and final segment
      if ( i_begin == i_end ) { // stesso segmento
        real_type   ss0 = m_s0[i_begin];
        LineSegment LL  = m_polylineList[i_begin];
        LL.trim( s_begin-ss0, s_end-ss0 );
        newPL.push_back( LL );
      } else {
        LineSegment L0 = m_polylineList[i_begin];
        L0.trim( s_begin - m_s0[i_begin], L0.length() );
        newPL.push_back( L0 );

        for ( ++i_begin; i_begin < i_end; ++i_begin )
          newPL.push_back( m_polylineList[i_begin] );

        LineSegment L1 = m_polylineList[i_end];
        L1.trim( 0, s_end - m_s0[i_end] );
        newPL.push_back( L1 );
      }
    } else {
      LineSegment L0 = m_polylineList[i_begin];
      L0.trim( s_begin - m_s0[i_begin], L0.length() );
      newPL.push_back( L0 );

      for ( ++i_begin; i_begin < n_seg; ++i_begin )
        newPL.push_back( m_polylineList[i_begin] );

      for ( i_begin = 0; i_begin < i_end; ++i_begin )
        newPL.push_back( m_polylineList[i_begin] );

      LineSegment L1 = m_polylineList[i_end];
      L1.trim( 0, s_end - m_s0[i_end] );
      newPL.push_back( L1 );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build_AABBtree() const {

    #ifdef CLOTHOIDS_USE_THREADS
    std::lock_guard<std::mutex> lock(m_aabb_mutex);
    #endif

    if ( m_aabb_done ) return;

    integer ipos{0};
    integer nobj{ integer(m_polylineList.size()) };
    m_aabb_tree.set_max_num_objects_per_node( G2LIB_AABB_CUT );
    m_aabb_tree.allocate( nobj, 2 ); // nbox, space dimension
    real_type bbox_min[2], bbox_max[2];
    for ( auto const & line : m_polylineList ) {
      line.bbox( bbox_min[0], bbox_min[1], bbox_max[0], bbox_max[1] );
      m_aabb_tree.replace_bbox( bbox_min, bbox_max, ipos );
      ++ipos;
    }
    m_aabb_tree.build();
    m_aabb_done = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::init( real_type x0, real_type y0 ) {
    m_xe = x0;
    m_ye = y0;
    m_polylineList.clear();
    m_s0.clear();
    m_s0.emplace_back(0);
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( real_type x, real_type y ) {
    LineSegment s;
    s.build_2P( m_xe, m_ye, x, y );
    m_polylineList.emplace_back( s );
    real_type slast = m_s0.back() + s.length();
    m_s0.emplace_back( slast );
    m_xe = x;
    m_ye = y;
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( LineSegment const & C ) {
    m_polylineList.emplace_back( C );
    LineSegment & S = m_polylineList.back();
    S.change_origin( m_xe, m_ye );
    real_type slast = m_s0.back() + S.length();
    m_s0.emplace_back( slast );
    m_xe = S.x_end();
    m_ye = S.y_end();
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( CircleArc const & C, real_type tol ) {
    real_type L  = C.length();
    integer   ns = integer(ceil( L / C.len_tolerance( tol ) ));
    real_type tx = m_xe - C.x_begin();
    real_type ty = m_ye - C.y_begin();
    for ( integer i = 1; i < ns; ++i ) {
      real_type s = (i*L)/ns;
      this->push_back( tx + C.X(s), ty + C.Y(s) );
    }
    this->push_back( tx + C.x_end(), ty + C.y_end() );
    m_xe = tx + C.x_end();
    m_ye = ty + C.y_end();
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( Biarc const & B, real_type tol ) {
    CircleArc const & C0 = B.C0();
    CircleArc const & C1 = B.C1();
    real_type L0  = C0.length();
    real_type L1  = C1.length();
    integer   ns0 = integer(ceil( L0 / C0.len_tolerance( tol ) ));
    integer   ns1 = integer(ceil( L1 / C1.len_tolerance( tol ) ));

    real_type tx = m_xe - C0.x_begin();
    real_type ty = m_ye - C0.y_begin();

    for ( integer i = 1; i < ns0; ++i ) {
      real_type s = (i*L0)/ns0;
      this->push_back( tx + C0.X(s), ty + C0.Y(s) );
    }
    this->push_back( tx + C1.x_begin(), ty + C1.y_begin() );
    for ( integer i = 1; i < ns1; ++i ) {
      real_type s = (i*L1)/ns1;
      this->push_back( tx + C1.X(s), ty + C1.Y(s) );
    }
    this->push_back( tx + C1.x_end(), ty + C1.y_end() );
    m_xe = tx + C1.x_end();
    m_ye = ty + C1.y_end();
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidCurve const & C, real_type tol ) {

    real_type L    = C.length();
    real_type absk = max(abs(C.kappa_begin()), abs(C.kappa_end()));
    real_type tmp  = absk*tol - 1;
    integer   ns   = 1;
    if ( tmp > -1 ) ns = integer( ceil( L*absk/(2*(Utils::m_pi-acos(tmp))) ) );

    real_type tx = m_xe - C.x_begin();
    real_type ty = m_ye - C.y_begin();
    for ( integer i = 1; i < ns; ++i ) {
      real_type s = (i*L)/ns;
      this->push_back( tx + C.X(s), ty + C.Y(s) );
    }

    this->push_back( tx + C.x_end(), ty + C.y_end() );
    m_xe = tx + C.x_end();
    m_ye = ty + C.y_end();
    m_aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::push_back( ClothoidList const & L, real_type tol ) {
    integer ns = L.num_segments();
    for ( integer idx = 0; idx < ns; ++idx ) {
      ClothoidCurve const & C = L.get( idx );
      this->push_back( C, tol );
    }
    // aabb_done = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build(
    real_type const * x,
    real_type const * y,
    integer           npts
  ) {
    init( x[0], y[0] );
    for ( integer k = 1; k < npts; ++k )
      this->push_back( x[k], y[k] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( LineSegment const & C ) {
    init( C.x_begin(), C.y_begin() );
    this->push_back( C.x_end(), C.y_end() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( CircleArc const & C, real_type tol ) {
    init( C.x_begin(), C.y_begin() );
    this->push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( Biarc const & C, real_type tol ) {
    init( C.x_begin(), C.y_begin() );
    this->push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( ClothoidCurve const & C, real_type tol ) {
    init( C.x_begin(), C.y_begin() );
    this->push_back( C, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  PolyLine::build( ClothoidList const & L, real_type tol ) {
    init( L.x_begin(), L.y_begin() );
    this->push_back( L, tol );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  PolyLine::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::POLYLINE ) {
      PolyLine const & C = *static_cast<PolyLine const *>(pC);
      return this->collision( C );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::POLYLINE ) {
        PolyLine C(pC);
        return this->collision( C );
      } else {
        return G2lib::collision( this, pC );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  PolyLine::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::POLYLINE ) {
      PolyLine const & C = *static_cast<PolyLine const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::POLYLINE ) {
        PolyLine C(pC);
        return this->collision_ISO( offs, C, offs_C );
      } else {
        return G2lib::collision_ISO( this, offs, pC, offs_C );
      }
    }
  }

  void
  PolyLine::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::POLYLINE ) {
      PolyLine const & C = *static_cast<PolyLine const *>(pC);
      this->intersect( C, ilist );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::POLYLINE ) {
        PolyLine C(pC);
        this->intersect( C, ilist );
      } else {
        G2lib::intersect( this, pC, ilist );
      }
    }
  }

  void
  PolyLine::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::POLYLINE ) {
      PolyLine const & C = *static_cast<PolyLine const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::POLYLINE ) {
        PolyLine C(pC);
        this->intersect_ISO( offs, C, offs_C, ilist );
      } else {
        G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  PolyLine::closest_point_ISO(
    real_type   x,
    real_type   y,
    real_type & X,
    real_type & Y,
    real_type & S,
    real_type & T,
    real_type & DST
  ) const{
    UTILS_ASSERT0(
      !m_polylineList.empty(),
      "PolyLine::closest_point_ISO, empty list\n"
    );
    integer ipos{0};
    real_type X1, Y1, S1, T1, DST1;

    this->build_AABBtree();
    DST = Utils::Inf<real_type>();

    if ( m_aabb_tree.num_tree_nodes() > G2LIB_AABB_MIN_NODES && intersect_with_AABBtree ) {
      AABB_SET candidateList;
      real_type xy[2] = { x, y };
      m_aabb_tree.min_distance_candidates( xy, candidateList );
      UTILS_ASSERT(
        !candidateList.empty(),
        "PolyLine::closest_point_ISO, empty candidate list, #{}\n{}\n",
        candidateList.size(),
        m_aabb_tree.info()
      );
      for ( auto i : candidateList ) {
        LineSegment const & LS = m_polylineList[i];
        LS.closest_point_ISO( x, y, X1, Y1, S1, T1, DST1 );
        if ( DST1 < DST ) {
          DST  = DST1;
          X    = X1;
          Y    = Y1;
          S    = m_s0[i] + S1;
          T    = T1;
          ipos = i;
        }
      }
    } else {
      integer i{0};
      for ( LineSegment const & LS : m_polylineList ) {
        LS.closest_point_ISO( x, y, X1, Y1, S1, T1, DST1 );
        if ( DST1 < DST ) {
          DST  = DST1;
          X    = X1;
          Y    = Y1;
          S    = m_s0[ipos] + S1;
          T    = T1;
          ipos = i;
        }
        ++i;
      }
    }
    m_polylineList[size_t(ipos)].eval_ISO( S - m_s0[ipos], T, X1, Y1 );
    real_type err = hypot( x - X1, y - Y1 );
    real_type tol = (DST > 1 ? DST*machepsi1000 : machepsi1000);
    if ( err > tol ) return -(ipos+1);
    return ipos;
  }

  /*\
   |             _ _ _     _
   |    ___ ___ | | (_)___(_) ___  _ __
   |   / __/ _ \| | | / __| |/ _ \| '_ \
   |  | (_| (_) | | | \__ \ | (_) | | | |
   |   \___\___/|_|_|_|___/_|\___/|_| |_|
  \*/

  bool
  PolyLine::collision( PolyLine const & PL ) const {
    this->build_AABBtree();
    PL.build_AABBtree();

    AABB_MAP intersectList;
    m_aabb_tree.intersect_and_refine( PL.m_aabb_tree, intersectList );
    for ( auto const & I : intersectList ) {
      integer i = I.first;
      UTILS_ASSERT_DEBUG(
        i >= 0 && i < integer(m_polylineList.size()),
        "PolyLine::collision( PL ) i={} out of range [0,{})\n",
        i, m_polylineList.size()
      );
      LineSegment const & LS1 = m_polylineList[i];
      for ( auto const & j : I.second ) {
        UTILS_ASSERT_DEBUG(
          j >= 0 && j < integer(PL.m_polylineList.size()),
          "PolyLine::collision( PL ) j={} out of range [0,{})\n",
          j, PL.m_polylineList.size()
        );
        LineSegment const & LS2 = PL.m_polylineList[j];
        bool collide = LS1.collision( LS2 );
        if ( collide ) return true;
      }
    }
    return false;
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  void
  PolyLine::intersect(
    PolyLine const    & PL,
    vector<real_type> & ss0,
    vector<real_type> & ss1
  ) const {
    UTILS_ASSERT0( !m_polylineList.empty(),    "PolyLine::intersect, empty list\n" );
    UTILS_ASSERT0( !PL.m_polylineList.empty(), "PolyLine::intersect, empty secondary list\n" );
    if ( intersect_with_AABBtree ) {
      build_AABBtree();
      PL.build_AABBtree();
      AABB_MAP intersectList;
      m_aabb_tree.intersect_and_refine( PL.m_aabb_tree, intersectList );
      for ( auto const & ilist : intersectList ) {
        integer ipos0 = ilist.first;
        UTILS_ASSERT_DEBUG( ipos0 < integer(m_polylineList.size()), "PolyLine::intersect, bad ipos0 = {}\n", ipos0 );
        LineSegment const & LS0 = m_polylineList[ipos0];
        for ( auto const & ipos1 : ilist.second ) {
          UTILS_ASSERT_DEBUG( ipos1 < integer(PL.m_polylineList.size()), "PolyLine::intersect, bad ipos1 = {}\n", ipos1 );
          LineSegment const & LS1 = PL.m_polylineList[ipos1];
          real_type sss0, sss1;
          bool ok = LS0.intersect( LS1, sss0, sss1 );
          if ( ok ) {
            ss0.emplace_back( sss0 + m_s0[ipos0] );
            ss1.emplace_back( sss1 + PL.m_s0[ipos1] );
          }
        }
      }
    } else {
      ss0.clear();
      ss1.clear();
      integer ipos0{0};
      for ( auto const & LS0 : m_polylineList ) {
        integer ipos1{0};
        for ( auto const & LS1 : PL.m_polylineList ) {
          real_type sss0, sss1;
          bool ok = LS0.intersect( LS1, sss0, sss1 );
          if ( ok ) {
            ss0.emplace_back( sss0 + m_s0[ipos0] );
            ss1.emplace_back( sss1 + PL.m_s0[ipos1] );
          }
        }
      }
    }
  }

  void
  PolyLine::intersect(
    PolyLine const & pl,
    IntersectList  & ilist
  ) const {
    vector<real_type> s1, s2;
    this->intersect( pl, s1, s2 );
    ilist.reserve( ilist.size() + s1.size() );
    for ( size_t i=0; i < s1.size(); ++i )
      ilist.emplace_back( s1[i], s2[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, PolyLine const & P ) {
    fmt::print( stream,
      "nseg    = {}\n"
      "x_begin = {}\n"
      "y_begin = {}\n"
      "x_end   = {}\n"
      "y_end   = {}\n"
      "length  = {}\n",
      P.num_segments(),
      P.x_begin(),
      P.y_begin(),
      P.x_end(),
      P.y_end(),
      P.length()
    );
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

// EOF: PolyLine.cc
