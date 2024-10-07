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
#include "Clothoids_fmt.hh"

#include <cfloat>
#include <limits>

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace G2lib {

  using std::lower_bound;
  using std::vector;
  using std::swap;
  using std::abs;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::setup( GenericContainer const & gc ) {
    string cwhere{ fmt::format("BiarcList[{}]::setup( gc ):", this->name() ) };
    char const * where{ cwhere.c_str() };
    GenericContainer::vec_real_type const & x = gc.get_map_vec_real("x", where );
    GenericContainer::vec_real_type const & y = gc.get_map_vec_real("y", where );
    integer n{ integer(x.size()) };
    UTILS_ASSERT(
      n == integer( y.size() ),
      "BiarcList[{}]::setup( gc ) (size(x)={}) != (size(y)={})\n",
      this->name(), x.size(), y.size()
    );
    bool ok{true};
    if ( gc.exists("theta") ) {
      GenericContainer::vec_real_type const & theta = gc.get_map_vec_real("theta", where );
      UTILS_ASSERT(
        n == integer( theta.size() ),
        "BiarcList[{}]::setup( gc ) (size(x)={}) != (size(theta)={})\n",
        this->name(), x.size(), theta.size()
      );
      ok = this->build_G1( n, x.data(), y.data(), theta.data() );
    } else {
      ok = this->build_G1( n, x.data(), y.data() );
    }
    UTILS_ASSERT( ok, "BiarcList[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( LineSegment const & LS ) {
    this->reset_last_interval();
    this->init();
    this->push_back( LS );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( CircleArc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( ClothoidCurve const & ) {
    UTILS_ERROR("can convert from ClothoidCurve to BiarcList\n");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( Biarc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( BiarcList const & C ) {
    *this = C;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( PolyLine const & ) {
    UTILS_ERROR("can convert from PolyLine to BiarcList\n");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( ClothoidList const & ) {
    UTILS_ERROR("can convert from ClothoidList to BiarcList\n");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( Dubins const & ) {
    UTILS_ERROR("can convert from Dubins to BiarcList\n");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::build( Dubins3p const & ) {
    UTILS_ERROR("can convert from Dubins3p to BiarcList\n");
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/

  BiarcList::BiarcList( LineSegment const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  BiarcList::BiarcList( CircleArc const & C )   : BaseCurve( C.name() ) { this->build( C ); }
  BiarcList::BiarcList( Biarc const & C )       : BaseCurve( C.name() ) { this->build( C ); }
  BiarcList::BiarcList( PolyLine const & C )    : BaseCurve( C.name() ) { this->build( C ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  BiarcList::BiarcList( BaseCurve const * pC ) : BiarcList( pC->name() ) {

    G2LIB_DEBUG_MESSAGE( "BiarcList convert: {}\n", pC->type_name() );

    this->init();
    switch ( pC->type() ) {
    case CurveType::LINE:
      G2LIB_DEBUG_MESSAGE( "LineSegment -> Biarc\n" );
      this->push_back( *static_cast<LineSegment const *>(pC) );
      break;
    case CurveType::CIRCLE:
      G2LIB_DEBUG_MESSAGE( "CircleArc -> Biarc\n" );
      this->push_back( *static_cast<CircleArc const *>(pC) );
      break;
    case CurveType::BIARC:
      G2LIB_DEBUG_MESSAGE( "to -> Biarc\n" );
      this->push_back( *static_cast<Biarc const *>(pC) );
      break;
    case CurveType::POLYLINE:
      G2LIB_DEBUG_MESSAGE( "to -> PolyLine\n" );
      this->push_back( *static_cast<PolyLine const *>(pC) );
      break;
    case CurveType::BIARC_LIST:
      G2LIB_DEBUG_MESSAGE( "to -> BiarcList\n" );
      this->copy( *static_cast<BiarcList const *>(pC) );
      break;
    default:
      UTILS_ERROR(
        "BiarcList constructor cannot convert from: {}\n",
        pC->type_name()
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::init() {
    m_s0.clear();
    m_biarc_list.clear();
    this->reset_last_interval();
    m_aabb_done = false;
    m_aabb_triangles.clear();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::copy( BiarcList const & L ) {
    this->init();
    m_biarc_list.reserve(L.m_biarc_list.size());
    std::copy( L.m_biarc_list.begin(),
               L.m_biarc_list.end(),
               back_inserter(m_biarc_list) );
    m_s0.reserve(L.m_s0.size());
    std::copy( L.m_s0.begin(), L.m_s0.end(), back_inserter(m_s0) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  BiarcList::find_at_s( real_type & s ) const {
    #ifdef CLOTHOIDS_USE_THREADS
    std::unique_lock<std::mutex> lock(m_last_interval_mutex);
    auto id = std::this_thread::get_id();
    auto it = m_last_interval.find(id);
    if ( it == m_last_interval.end() ) {
      it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
      *it->second.get() = 0;
    }
    integer & last_interval{ *it->second.get() };
    lock.unlock();
    #else
    integer & last_interval = m_last_interval;
    #endif
    Utils::search_interval<integer,real_type>(
      static_cast<integer>(m_s0.size()),
      &m_s0.front(), s, last_interval, false, true
    );
    return last_interval;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::reserve( integer n ) {
    m_s0.reserve(size_t(n+1));
    m_biarc_list.reserve(size_t(n));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back( LineSegment const & LS ) {
    if ( m_biarc_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(LS.length());
    } else {
      m_s0.emplace_back(m_s0.back()+LS.length());
    }
    Biarc tmp(&LS);
    m_biarc_list.push_back(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back( CircleArc const & C ) {
    if ( m_biarc_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(C.length());
    } else {
      m_s0.emplace_back(m_s0.back()+C.length());
    }
    Biarc tmp(&C);
    m_biarc_list.push_back(tmp);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back( Biarc const & c ) {
    if ( m_biarc_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(c.length());
    } else {
      m_s0.emplace_back(m_s0.back()+c.length());
    }
    m_biarc_list.push_back(c);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back( PolyLine const & c ) {
    m_s0.reserve( m_s0.size() + c.m_polyline_list.size() + 1 );
    m_biarc_list.reserve( m_biarc_list.size() + c.m_polyline_list.size() );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    for ( LineSegment const & LS : c.m_polyline_list ) {
      m_s0.emplace_back(m_s0.back()+LS.length());
      Biarc B(&LS);
      m_biarc_list.push_back(B);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back_G1(
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    UTILS_ASSERT0( !m_biarc_list.empty(), "BiarcList::push_back_G1(...) empty list!\n" );
    Biarc c{"BiarcList::push_back_G1 temporary c"};
    real_type x0     = m_biarc_list.back().x_end();
    real_type y0     = m_biarc_list.back().y_end();
    real_type theta0 = m_biarc_list.back().theta_end();
    c.build( x0, y0, theta0, x1, y1, theta1 );
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::push_back_G1(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    Biarc c{"BiarcList::push_back_G1 temporary c"};
    c.build( x0, y0, theta0, x1, y1, theta1 );
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  BiarcList::build_G1(
    integer         n,
    real_type const x[],
    real_type const y[],
    real_type const theta[]
  ) {
    UTILS_ASSERT0(
      n > 1,
      "BiarcList::build_G1, at least 2 points are necessary\n"
    );
    init();
    reserve( n-1 );
    Biarc c{"BiarcList::build_G1 temporary c"};
    for ( integer k = 1; k < n; ++k ) {
      c.build( x[k-1], y[k-1], theta[k-1], x[k], y[k], theta[k] );
      this->push_back(c);
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  BiarcList::build_G1(
    integer         n,
    real_type const x[],
    real_type const y[]
  ) {
    size_t nn = size_t(n);
    Utils::Malloc<real_type> mem( "BiarcList::build_G1" );
    mem.allocate( 5 * nn );
    real_type * theta     = mem( nn );
    real_type * theta_min = mem( nn );
    real_type * theta_max = mem( nn );
    real_type * omega     = mem( nn );
    real_type * len       = mem( nn );
    G2lib::xy_to_guess_angle(
      n, x, y, theta, theta_min, theta_max, omega, len
    );
    return this->build_G1( n, x, y, theta );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Biarc const &
  BiarcList::get( integer idx ) const {
    UTILS_ASSERT(
      !m_biarc_list.empty(),
      "BiarcList::get( {} ) empty list\n", idx
    );
    try {
      return m_biarc_list.at(idx);
    } catch ( std::exception & exc ) {
      UTILS_ERROR( "BiarcList::get( {} ): {}\n", idx, exc.what() );
    } catch ( ... ) {
      UTILS_ERROR( "BiarcList::get( {} ): unknown error\n", idx );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Biarc const &
  BiarcList::get_at_s( real_type s ) const
  { return this->get(find_at_s(s)); }

  /*\
   |   _                  _   _
   |  | | ___ _ __   __ _| |_| |__
   |  | |/ _ \ '_ \ / _` | __| '_ \
   |  | |  __/ | | | (_| | |_| | | |
   |  |_|\___|_| |_|\__, |\__|_| |_|
   |                |___/
  \*/

  real_type
  BiarcList::length() const {
    return m_s0.back() - m_s0.front();
  }

  real_type
  BiarcList::length_ISO( real_type offs ) const {
    real_type L{0};
    for ( auto const & b : m_biarc_list ) L += b.length_ISO( offs );
    return L;
  }

  real_type
  BiarcList::segment_length( integer nseg ) const {
    Biarc const & c = this->get( nseg );
    return c.length();
  }

  real_type
  BiarcList::segment_length_ISO( integer nseg, real_type offs ) const {
    Biarc const & c = this->get( nseg );
    return c.length_ISO( offs );
  }

  /*\
   |   _    _   _____    _                _
   |  | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
   |  | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
   |  |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
   |                                |___/
  \*/

  void
  BiarcList::bb_triangles(
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<Biarc>::const_iterator ic = m_biarc_list.begin();
    for ( integer ipos = icurve; ic != m_biarc_list.end(); ++ic, ++ipos )
      ic->bb_triangles( tvec, max_angle, max_size, ipos );
  }

  void
  BiarcList::bb_triangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<Biarc>::const_iterator ic = m_biarc_list.begin();
    for ( integer ipos = icurve; ic != m_biarc_list.end(); ++ic, ++ipos )
      ic->bb_triangles_ISO( offs, tvec, max_angle, max_size, ipos );
  }

  void
  BiarcList::bb_triangles_SAE(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<Biarc>::const_iterator ic = m_biarc_list.begin();
    for ( integer ipos = icurve; ic != m_biarc_list.end(); ++ic, ++ipos )
      ic->bb_triangles_SAE( offs, tvec, max_angle, max_size, ipos );
  }

  /*\
   |   _     _
   |  | |__ | |__   _____  __
   |  | '_ \| '_ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |_.__/|_.__/ \___/_/\_\
  \*/

  void
  BiarcList::bbox_ISO(
    real_type   offs,
    real_type & xmin,
    real_type & ymin,
    real_type & xmax,
    real_type & ymax
  ) const {
    vector<Triangle2D> tvec;
    bb_triangles_ISO( offs, tvec, Utils::m_pi/18, 1e100 );
    xmin = ymin = Utils::Inf<real_type>();
    xmax = ymax = -xmin;
    for ( auto const & t : tvec ) {
      // - - - - - - - - - - - - - - - - - - - -
      if      ( t.x1() < xmin ) xmin = t.x1();
      else if ( t.x1() > xmax ) xmax = t.x1();
      if      ( t.x2() < xmin ) xmin = t.x2();
      else if ( t.x2() > xmax ) xmax = t.x2();
      if      ( t.x3() < xmin ) xmin = t.x3();
      else if ( t.x3() > xmax ) xmax = t.x3();
      // - - - - - - - - - - - - - - - - - - - -
      if      ( t.y1() < ymin ) ymin = t.y1();
      else if ( t.y1() > ymax ) ymax = t.y1();
      if      ( t.y2() < ymin ) ymin = t.y2();
      else if ( t.y2() > ymax ) ymax = t.y2();
      if      ( t.y3() < ymin ) ymin = t.y3();
      else if ( t.y3() > ymax ) ymax = t.y3();
    }
  }

  /*\
   |  _   _          _
   | | |_| |__   ___| |_ __ _
   | | __| '_ \ / _ \ __/ _` |
   | | |_| | | |  __/ || (_| |
   |  \__|_| |_|\___|\__\__,_|
  \*/

  real_type
  BiarcList::theta( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.theta( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::theta_D( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.theta_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::theta_DD( real_type s ) const  {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.theta_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::theta_DDD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.theta_DDD( s - m_s0[idx] );
  }

  /*\
   |  _____                   _   _   _
   | |_   _|   __ _ _ __   __| | | \ | |
   |   | |    / _` | '_ \ / _` | |  \| |
   |   | |   | (_| | | | | (_| | | |\  |
   |   |_|    \__,_|_| |_|\__,_| |_| \_|
  \*/

  real_type
  BiarcList::tx( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tx( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::ty( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.ty( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::tx_D( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tx_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::ty_D( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.ty_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::tx_DD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tx_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::ty_DD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.ty_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::tx_DDD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tx_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::ty_DDD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.ty_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::tg(
    real_type   s,
    real_type & tg_x,
    real_type & tg_y
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tg( s - m_s0[idx], tg_x, tg_y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::tg_D(
    real_type   s,
    real_type & tg_x_D,
    real_type & tg_y_D
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tg_D( s - m_s0[idx], tg_x_D, tg_y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::tg_DD(
    real_type   s,
    real_type & tg_x_DD,
    real_type & tg_y_DD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tg_DD( s - m_s0[idx], tg_x_DD, tg_y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::tg_DDD(
    real_type   s,
    real_type & tg_x_DDD,
    real_type & tg_y_DDD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.tg_DDD( s - m_s0[idx], tg_x_DDD, tg_y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::evaluate(
    real_type   s,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    c.evaluate( s - m_s0[idx], th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::evaluate_ISO(
    real_type   s,
    real_type   offs,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    c.evaluate_ISO( s - m_s0[idx], offs, th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y( real_type s ) const  {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_D( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_D( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_DD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_DD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_DDD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_DDD( real_type s ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval( s - m_s0[idx], x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_D( s - m_s0[idx], x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_DD( s - m_s0[idx], x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_DDD( s - m_s0[idx], x_DDD, y_DDD );
  }

  /*\
   |         __  __          _
   |   ___  / _|/ _|___  ___| |_
   |  / _ \| |_| |_/ __|/ _ \ __|
   | | (_) |  _|  _\__ \  __/ |_
   |  \___/|_| |_| |___/\___|\__|
  \*/

  real_type
  BiarcList::X_ISO( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_ISO( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_ISO( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_ISO( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_ISO_D( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_ISO_D( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_ISO_D( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_ISO_D( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_ISO_DD( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_ISO_DD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_ISO_DD( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_ISO_DD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::X_ISO_DDD( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.X_ISO_DDD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  BiarcList::Y_ISO_DDD( real_type s, real_type offs ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.Y_ISO_DDD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_ISO( s - m_s0[idx], offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_ISO_D( s - m_s0[idx], offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_ISO_DD( s - m_s0[idx], offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    integer idx = this->find_at_s( s );
    Biarc const & c = this->get( idx );
    return c.eval_ISO_DDD( s - m_s0[idx], offs, x_DDD, y_DDD );
  }

  /*\
   |  _                        __
   | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
   | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
   | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
   |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
  \*/

  void
  BiarcList::translate( real_type tx, real_type ty ) {
    for ( auto & B : m_biarc_list ) B.translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::rotate( real_type angle, real_type cx, real_type cy ) {
    for ( auto & B : m_biarc_list ) B.rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::scale( real_type sfactor ) {
    vector<Biarc>::iterator ic = m_biarc_list.begin();
    real_type newx0 = ic->x_begin();
    real_type newy0 = ic->y_begin();
    m_s0[0] = 0;
    for ( size_t k=0; ic != m_biarc_list.end(); ++ic, ++k ) {
      ic->scale( sfactor );
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::reverse() {
    std::reverse( m_biarc_list.begin(), m_biarc_list.end() );
    vector<Biarc>::iterator ic = m_biarc_list.begin();
    ic->reverse();
    real_type newx0 = ic->x_end();
    real_type newy0 = ic->y_end();
    m_s0[0] = 0;
    m_s0[1] = ic->length();
    size_t k = 1;
    for ( ++ic; ic != m_biarc_list.end(); ++ic, ++k ) {
      ic->reverse();
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::change_origin( real_type newx0, real_type newy0 ) {
    for ( auto & B : m_biarc_list ) {
      B.change_origin( newx0, newy0 );
      newx0 = B.x_end();
      newy0 = B.y_end();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::trim( real_type s_begin, real_type s_end ) {
    UTILS_ASSERT(
      s_begin >= m_s0.front() && s_end <= m_s0.back() && s_end > s_begin,
      "BiarcList::trim( s_begin={}, s_end={} ) bad range, must be in [ {}, {} ]\n",
      s_begin, s_end, m_s0.front(), m_s0.back()
    );

    size_t i_begin = size_t( find_at_s( s_begin ) );
    size_t i_end   = size_t( find_at_s( s_end ) );
    if ( i_begin == i_end ) {
      m_biarc_list[i_begin].trim( s_begin-m_s0[i_begin], s_end-m_s0[i_begin] );
    } else {
      m_biarc_list[i_begin].trim( s_begin-m_s0[i_begin], m_s0[i_begin+1]-m_s0[i_begin] );
      m_biarc_list[i_end].trim( 0, s_end-m_s0[i_end] );
    }
    m_biarc_list.erase( m_biarc_list.begin()+i_end+1, m_biarc_list.end() );
    m_biarc_list.erase( m_biarc_list.begin(), m_biarc_list.begin()+i_begin );
    if ( m_biarc_list.back().length() <= machepsi100 ) m_biarc_list.pop_back();
    vector<Biarc>::iterator ic = m_biarc_list.begin();
    m_s0.resize( m_biarc_list.size() + 1 );
    m_s0[0] = 0;
    size_t k{0};
    for ( ++ic; ic != m_biarc_list.end(); ++ic, ++k )
      m_s0[k+1] = m_s0[k] + ic->length();
    this->reset_last_interval();
  }

  /*\
   |     _        _    ____  ____  _
   |    / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |   / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |  / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   | /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  void
  BiarcList::build_AABBtree_ISO(
    real_type offs,
    real_type max_angle,
    real_type max_size
  ) const {

    #ifdef CLOTHOIDS_USE_THREADS
    std::lock_guard<std::mutex> lock(m_aabb_mutex);
    #endif

    if ( m_aabb_done &&
         Utils::is_zero( offs-m_aabb_offs ) &&
         Utils::is_zero( max_angle-m_aabb_max_angle ) &&
         Utils::is_zero( max_size-m_aabb_max_size ) ) return;

    bb_triangles_ISO( offs, m_aabb_triangles, max_angle, max_size );

    integer ipos{0};
    integer nobj{ integer(m_aabb_triangles.size()) };
    m_aabb_tree.set_max_num_objects_per_node( G2LIB_AABB_CUT );
    m_aabb_tree.allocate( nobj, 2 ); // nbox, space dimension
    real_type bbox_min[2], bbox_max[2];
    for ( Triangle2D const & T : m_aabb_triangles ) {
      T.bbox( bbox_min[0], bbox_min[1], bbox_max[0], bbox_max[1] );
      m_aabb_tree.replace_bbox( bbox_min, bbox_max, ipos );
      ++ipos;
    }

    m_aabb_tree.build();
    m_aabb_done      = true;
    m_aabb_offs      = offs;
    m_aabb_max_angle = max_angle;
    m_aabb_max_size  = max_size;
  }
  #endif

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  BiarcList::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::BIARC_LIST ) {
      BiarcList const & C = *static_cast<BiarcList const *>(pC);
      return this->collision( C );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::BIARC_LIST ) {
        BiarcList C(pC);
        return this->collision( C );
      } else {
        return G2lib::collision( this, pC );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  BiarcList::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::BIARC_LIST ) {
      BiarcList const & C = *static_cast<BiarcList const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::BIARC_LIST ) {
        BiarcList C(pC);
        return this->collision_ISO( offs, C, offs_C );
      } else {
        return G2lib::collision_ISO( this, offs, pC, offs_C );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  BiarcList::collision_ISO(
    real_type         offs,
    BiarcList const & BL,
    real_type         offs_BL
  ) const {
    this->build_AABBtree_ISO( offs );
    BL.build_AABBtree_ISO( offs_BL );
    AABB_MAP intersectList;
    m_aabb_tree.intersect_and_refine( BL.m_aabb_tree, intersectList );
    for ( auto const & I : intersectList ) {
      integer i = I.first;
      UTILS_ASSERT_DEBUG(
        i >= 0 && i < integer(m_aabb_triangles.size()),
        "BiarcList::collision_ISO( offs={}, BL, offs_BL={} ) i={} out of range [0,{})\n",
        offs, offs_BL, i, m_aabb_triangles.size()
      );
      Triangle2D const & T1  = m_aabb_triangles.at(i);
      Biarc      const & BA1 = m_biarc_list.at(T1.Icurve());
      for ( auto const & j : I.second ) {
        UTILS_ASSERT_DEBUG(
          j >= 0 && j < integer(BL.m_aabb_triangles.size()),
          "BiarcList::collision_ISO( offs={}, BL, offs_BL={} ) j={} out of range [0,{})\n",
          offs, offs_BL, j, BL.m_aabb_triangles.size()
        );
        Triangle2D const & T2  = BL.m_aabb_triangles.at(j);
        Biarc      const & BA2 = BL.m_biarc_list.at(T2.Icurve());
        bool collide = BA1.collision_ISO( offs, BA2, offs_BL );
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
   |
  \*/

  void
  BiarcList::intersect_ISO(
    real_type         offs,
    BiarcList const & BL,
    real_type         offs_BL,
    IntersectList   & ilist
  ) const {

    if ( intersect_with_AABBtree ) {

      this->build_AABBtree_ISO( offs );
      BL.build_AABBtree_ISO( offs_BL );
      AABB_MAP intersectList;
      m_aabb_tree.intersect_and_refine( BL.m_aabb_tree, intersectList );
      for ( auto const & I : intersectList ) {
        integer i = I.first;
        UTILS_ASSERT_DEBUG(
          i >= 0 && i < integer(m_aabb_triangles.size()),
          "BiarcList::intersect_ISO( offs={}, BL, offs_BL={}, ilist ) i={} out of range [0,{})\n",
          offs, offs_BL, i, m_aabb_triangles.size()
        );
        Triangle2D const & T1  = m_aabb_triangles.at(i);
        Biarc      const & BA1 = m_biarc_list.at(T1.Icurve());

        for ( integer j : I.second ) {
          UTILS_ASSERT_DEBUG(
            j >= 0 && j < integer(BL.m_aabb_triangles.size()),
            "BiarcList::intersect_ISO( offs={}, BL, offs_BL={}, ilist ) j={} out of range [0,{})\n",
            offs, offs_BL, j, BL.m_aabb_triangles.size()
          );
          Triangle2D const & T2  = BL.m_aabb_triangles.at(j);
          Biarc      const & BA2 = BL.m_biarc_list.at(T2.Icurve());

          IntersectList ilist1;
          BA1.intersect_ISO( offs, BA2, offs_BL, ilist1 );

          for ( auto const & it : ilist1 ) {
            real_type ss1 = it.first  + m_s0.at(T1.Icurve());
            real_type ss2 = it.second + BL.m_s0.at(T2.Icurve());
            ilist.emplace_back( ss1, ss2 );
          }
        }
      }

    } else {

      bb_triangles_ISO( offs, m_aabb_triangles, Utils::m_pi/18, 1e100 );
      BL.bb_triangles_ISO( offs_BL, BL.m_aabb_triangles, Utils::m_pi/18, 1e100 );

      for ( Triangle2D const & T1 : m_aabb_triangles ) {
        Biarc const & BA1 = m_biarc_list.at(T1.Icurve());

        for ( Triangle2D const & T2 : BL.m_aabb_triangles ) {
          Biarc const & BA2 = BL.m_biarc_list.at(T2.Icurve());

          IntersectList ilist1;
          BA1.intersect_ISO( offs, BA2, offs_BL, ilist1 );

          for ( auto const & it : ilist1 ) {
            real_type ss1 = it.first  + m_s0.at(T1.Icurve());
            real_type ss2 = it.second + BL.m_s0.at(T2.Icurve());
            ilist.emplace_back( ss1, ss2 );
          }
        }
      }
    }
  }

  void
  BiarcList::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::BIARC_LIST ) {
      BiarcList const & C = *static_cast<BiarcList const *>(pC);
      this->intersect( C, ilist );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::BIARC_LIST ) {
        BiarcList C(pC);
        this->intersect( C, ilist );
      } else {
        G2lib::intersect( this, pC, ilist );
      }
    }
  }

  void
  BiarcList::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::BIARC_LIST ) {
      BiarcList const & C = *static_cast<BiarcList const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::BIARC_LIST ) {
        BiarcList C(pC);
        this->intersect_ISO( offs, C, offs_C, ilist );
      } else {
        G2lib::intersect_ISO( this, offs, pC, offs_C, ilist );
      }
    }
  }

  /*\
   |      _ _     _
   |   __| (_)___| |_ __ _ _ __   ___ ___
   |  / _` | / __| __/ _` | '_ \ / __/ _ \
   | | (_| | \__ \ || (_| | | | | (_|  __/
   |  \__,_|_|___/\__\__,_|_| |_|\___\___|
  \*/

  integer
  BiarcList::closest_point_internal(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & DST
  ) const {

    this->build_AABBtree_ISO( offs );

    integer icurve{0};
    DST = Utils::Inf<real_type>();

    if ( m_aabb_tree.num_tree_nodes() > G2LIB_AABB_MIN_NODES && intersect_with_AABBtree ) {

      AABB_SET candidateList;
      real_type xy[2] = { qx, qy };
      m_aabb_tree.min_distance_candidates( xy, candidateList );
      UTILS_ASSERT0(
         candidateList.size() > 0,
        "BiarcList::closest_point_internal no candidate\n"
      );
      for ( integer ipos : candidateList ) {
        Triangle2D const & T = m_aabb_triangles.at(ipos);
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss, tt;
          m_biarc_list.at(T.Icurve()).closest_point_ISO( qx, qy, offs, xx, yy, ss, tt, dst );
          if ( dst < DST ) {
            DST    = dst;
            s      = ss + m_s0.at(T.Icurve());
            x      = xx;
            y      = yy;
            icurve = T.Icurve();
          }
        }
      }

    } else {

      for ( Triangle2D const & T : m_aabb_triangles ) {
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss, tt;
          m_biarc_list.at(T.Icurve()).closest_point_ISO( qx, qy, offs, xx, yy, ss, tt, dst );
          if ( dst < DST ) {
            DST    = dst;
            s      = ss + m_s0.at(T.Icurve());
            x      = xx;
            y      = yy;
            icurve = T.Icurve();
          }
        }
      }
    }
    return icurve;
  }

  integer
  BiarcList::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & DST
  ) const {

    integer icurve = this->closest_point_internal( qx, qy, offs, x, y, s, DST );

    // check if projection is orthogonal
    real_type nx, ny;
    m_biarc_list.at(icurve).nor_ISO( s - m_s0.at(icurve), nx, ny );
    real_type qxx = qx - x;
    real_type qyy = qy - y;
    t = qxx * nx + qyy * ny - offs; // signed distance
    real_type pt = abs(qxx * ny - qyy * nx);
    G2LIB_DEBUG_MESSAGE(
      "BiarcList::closest_point_ISO\n"
      "||P-P0|| = {} and {}, |(P-P0).T| = {}\n",
      DST, hypot(qxx,qyy), pt
    );
    return pt > GLIB2_TOL_ANGLE*hypot(qxx,qyy) ? -(icurve+1) : icurve;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  BiarcList::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst
  ) const {
    return closest_point_ISO( qx, qy, 0, x, y, s, t, dst );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::get_STK(
    real_type s[],
    real_type theta[],
    real_type kappa[]
  ) const {
    vector<Biarc>::const_iterator ic = m_biarc_list.begin();
    integer   k{0};
    real_type ss{0};
    while ( ic != m_biarc_list.end() ) {
      s[k]     = ss;
      theta[k] = ic->theta_begin();
      kappa[k] = ic->kappa_begin();
      ss      += ic->length();
      ++k;
      ++ic;
    }
    --ic;
    s[k]     = ss;
    theta[k] = ic->theta_end();
    kappa[k] = ic->kappa_end();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  BiarcList::get_XY( real_type x[], real_type y[] ) const {
    vector<Biarc>::const_iterator ic = m_biarc_list.begin();
    integer k{0};
    while ( ic != m_biarc_list.end() ) {
      x[k] = ic->x_begin();
      y[k] = ic->y_begin();
      ++k; ++ic;
    }
    --ic;
    x[k] = ic->x_end();
    y[k] = ic->y_end();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  BiarcList::findST1(
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    UTILS_ASSERT0( !m_biarc_list.empty(), "BiarcList::findST, empty list\n" );
    vector<Biarc>::const_iterator     ic = m_biarc_list.begin();
    vector<real_type>::const_iterator is = m_s0.begin();

    s = t = 0;
    integer ipos{0};
    integer iseg{0};
    real_type S, T;
    bool ok = ic->findST_ISO( x, y, S, T );
    if ( ok ) {
      s    = *is + S;
      t    = T;
      iseg = 0;
    }

    for ( ++ic, ++is, ++ipos;
          ic != m_biarc_list.end();
          ++ic, ++is, ++ipos ) {
      bool ok1 = ic->findST_ISO( x, y, S, T );
      if ( ok && ok1 ) ok1 = abs(T) < abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = *is + S;
        t    = T;
        iseg = ipos;
      }
    }

    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  BiarcList::findST1(
    integer     ibegin,
    integer     iend,
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    UTILS_ASSERT0(
      !m_biarc_list.empty(), "BiarcList::findST, empty list\n"
    );
    UTILS_ASSERT(
      ibegin >= 0 && ibegin <= iend &&
      iend < integer(m_biarc_list.size()),
      "BiarcList::findST( ibegin={}, iend={}, x, y, s, t )\n"
      "bad range not in [0,{}]\n",
      ibegin, iend, m_biarc_list.size()-1
    );
    s = t = 0;
    integer iseg = 0;
    bool ok = false;
    for ( integer k = ibegin; k <= iend; ++k ) {
      Biarc const & ck = m_biarc_list[k];
      real_type S, T;
      bool ok1 = ck.findST_ISO( x, y, S, T );
      if ( ok && ok1 ) ok1 = abs(T) < abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = m_s0[k] + S;
        t    = T;
        iseg = k;
      }
    }
    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  BiarcList::info() const
  { return fmt::format( "BiarcList\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `BiarcList` object
  //!
  //!  \param stream the output stream
  //!  \param CL     an instance of `BiarcList` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, BiarcList const & CL ) {
    for ( auto const & b : CL.m_biarc_list )
      stream << b << '\n';
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

// EOF: BiarcList.cc
