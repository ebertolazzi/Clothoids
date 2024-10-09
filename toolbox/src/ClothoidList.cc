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
#include "Utils_Algo748.hh"

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

  #ifdef G2LIB_DEBUG
    #define G2LIB_DEBUG_TICTOC Utils::TicToc tictoc
    #define G2LIB_DEBUG_TIC    tictoc.tic()
    #define G2LIB_DEBUG_TOC    tictoc.toc()
  #else
    #define G2LIB_DEBUG_TICTOC
    #define G2LIB_DEBUG_TIC
    #define G2LIB_DEBUG_TOC
  #endif

  using std::lower_bound;
  using std::vector;
  using std::swap;
  using std::abs;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::setup( GenericContainer const & gc ) {
    string cwhere{ fmt::format("ClothoidList[{}]::setup( gc ):", this->name() ) };
    char const * where{ cwhere.c_str() };
    GenericContainer::vec_real_type const & x = gc.get_map_vec_real("x", where );
    GenericContainer::vec_real_type const & y = gc.get_map_vec_real("y", where );
    integer n{ integer(x.size()) };
    UTILS_ASSERT(
      n == integer( y.size() ),
      "ClothoidList[{}]::setup( gc ) (size(x)={}) != (size(y)={})\n",
      this->name(), x.size(), y.size()
    );
    bool ok{true};
    if ( gc.exists("theta") ) {
      GenericContainer::vec_real_type const & theta = gc.get_map_vec_real("theta", where );
      UTILS_ASSERT(
        n == integer( theta.size() ),
        "ClothoidList[{}]::setup( gc ) (size(x)={}) != (size(theta)={})\n",
        this->name(), x.size(), theta.size()
      );
      ok = this->build_G1( n, x.data(), y.data(), theta.data() );
    } else {
      ok = this->build_G1( n, x.data(), y.data() );
    }
    UTILS_ASSERT( ok, "ClothoidList[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( LineSegment const & LS ) {
    this->reset_last_interval();
    this->init();
    this->push_back( LS );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( CircleArc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( Biarc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C.C0() );
    this->push_back( C.C1() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( BiarcList const & c ) {
    this->reset_last_interval();
    this->init();
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( ClothoidCurve const & c ) {
    this->reset_last_interval();
    this->init();
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( PolyLine const & pl ) {
    this->reset_last_interval();
    this->init();
    this->push_back( pl );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( ClothoidList const & pl ) {
    this->reset_last_interval();
    this->init();
    this->push_back( pl );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( G2solve2arc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( G2solve3arc const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( G2solveCLC const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::build( Dubins const & C ) {
    this->reset_last_interval();
    this->init();
    this->push_back( C );
  }

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/

  ClothoidList::ClothoidList( LineSegment   const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  ClothoidList::ClothoidList( CircleArc     const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  ClothoidList::ClothoidList( Biarc         const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  ClothoidList::ClothoidList( BiarcList     const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  ClothoidList::ClothoidList( ClothoidCurve const & C ) : BaseCurve( C.name() ) { this->build( C ); }
  ClothoidList::ClothoidList( PolyLine      const & C ) : BaseCurve( C.name() ) { this->build( C ); }

  ClothoidList::ClothoidList( G2solve2arc const & C, string const & name ) : BaseCurve( name ) { this->build( C ); }
  ClothoidList::ClothoidList( G2solve3arc const & C, string const & name ) : BaseCurve( name ) { this->build( C ); }
  ClothoidList::ClothoidList( G2solveCLC  const & C, string const & name ) : BaseCurve( name ) { this->build( C ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidList::ClothoidList( BaseCurve const * pC ) : BaseCurve( pC->name() )  {

    G2LIB_DEBUG_MESSAGE( "ClothoidList convert: {}\n", pC->type_name() );

    this->reset_last_interval();
    this->init();
    switch ( pC->type() ) {
    case CurveType::LINE:
      G2LIB_DEBUG_MESSAGE( "LineSegment -> ClothoidList\n" );
      this->push_back( *static_cast<LineSegment const *>(pC) );
      break;
    case CurveType::CIRCLE:
      G2LIB_DEBUG_MESSAGE( "CircleArc -> ClothoidList\n" );
      this->push_back( *static_cast<CircleArc const *>(pC) );
      break;
    case CurveType::BIARC:
      G2LIB_DEBUG_MESSAGE( "Biarc -> ClothoidList\n" );
      this->push_back( *static_cast<Biarc const *>(pC) );
      break;
    case CurveType::CLOTHOID:
      G2LIB_DEBUG_MESSAGE( "ClothoidCurve -> ClothoidList\n" );
      this->push_back( *static_cast<ClothoidCurve const *>(pC) );
      break;
    case CurveType::POLYLINE:
      G2LIB_DEBUG_MESSAGE( "PolyLine -> ClothoidList\n" );
      this->push_back( *static_cast<PolyLine const *>(pC) );
      break;
    case CurveType::BIARC_LIST:
      G2LIB_DEBUG_MESSAGE( "BiarcList -> ClothoidList\n" );
      this->push_back( *static_cast<BiarcList const *>(pC) );
      break;
    case CurveType::CLOTHOID_LIST:
      this->copy( *static_cast<ClothoidList const *>(pC) );
      break;
    case CurveType::DUBINS:
      this->push_back( *static_cast<Dubins const *>(pC) );
      break;
    case CurveType::DUBINS3P:
      this->push_back( *static_cast<Dubins3p const *>(pC) );
      break;
    //default:
    //  UTILS_ERROR(
    //    "ClothoidList::ClothoidList, missing conversion for type {}",
    //    pC->type_name()
    //  );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::find_at_s( real_type & s ) const {
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
      m_s0.data(), s, last_interval, m_curve_is_closed, true
    );
    return last_interval;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::wrap_in_range( real_type & s ) const {
    real_type a = m_s0.front();
    real_type b = m_s0.back();
    real_type L = b-a;
    s -= a;
    s  = fmod( s, L );
    if ( s < 0 ) s += L;
    // while ( s < 0 ) s += L;
    // while ( s > L ) s -= L;
    s += a;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::init() {
    m_s0.clear();
    m_clothoid_list.clear();
    m_aabb_done = false;
    m_aabb_triangles.clear();
    this->reset_last_interval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::copy( ClothoidList const & L ) {
    this->init();
    m_clothoid_list.reserve( L.m_clothoid_list.size() );
    std::copy(
      L.m_clothoid_list.begin(), L.m_clothoid_list.end(), back_inserter(m_clothoid_list)
    );
    m_s0.reserve( L.m_s0.size() );
    std::copy( L.m_s0.begin(), L.m_s0.end(), back_inserter(m_s0) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reserve( integer n ) {
    m_s0.reserve(size_t(n+1));
    m_clothoid_list.reserve(size_t(n));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( LineSegment const & LS ) {
    if ( m_clothoid_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(LS.length());
    } else {
      m_s0.emplace_back(m_s0.back()+LS.length());
    }
    m_clothoid_list.emplace_back(LS);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( CircleArc const & C ) {
    if ( m_clothoid_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(C.length());
    } else {
      m_s0.emplace_back(m_s0.back()+C.length());
    }
    m_clothoid_list.emplace_back(C);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( Biarc const & c ) {
    if ( m_clothoid_list.empty() ) m_s0.emplace_back( 0 );
    CircleArc const & C0 = c.C0();
    CircleArc const & C1 = c.C1();
    m_s0.emplace_back( m_s0.back()+C0.length() );
    m_s0.emplace_back( m_s0.back()+C1.length() );
    m_clothoid_list.emplace_back( C0 );
    m_clothoid_list.emplace_back( C1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( ClothoidCurve const & c ) {
    if ( m_clothoid_list.empty() ) {
      m_s0.emplace_back(0);
      m_s0.emplace_back(c.length());
    } else {
      m_s0.emplace_back(m_s0.back()+c.length());
    }
    m_clothoid_list.emplace_back(c);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( BiarcList const & c ) {
    m_s0.reserve( m_s0.size() + c.m_biarc_list.size() + 1 );
    m_clothoid_list.reserve( m_clothoid_list.size() + 2*c.m_biarc_list.size() );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    for ( auto const & b : c.m_biarc_list ) {
      m_s0.emplace_back(m_s0.back()+b.length());
      m_clothoid_list.emplace_back(b.C0());
      m_clothoid_list.emplace_back(b.C1());
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( PolyLine const & c ) {
    m_s0.reserve( m_s0.size() + c.m_polyline_list.size() + 1 );
    m_clothoid_list.reserve( m_clothoid_list.size() + c.m_polyline_list.size() );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    for ( auto const & L : c.m_polyline_list ) {
      m_s0.emplace_back(m_s0.back()+L.length());
      ClothoidCurve C(&L);
      m_clothoid_list.emplace_back(C);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( ClothoidList const & c ) {
    m_s0.reserve( m_s0.size() + c.m_clothoid_list.size() + 1 );
    m_clothoid_list.reserve( m_clothoid_list.size() + c.m_clothoid_list.size() );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    for ( auto const & C : c.m_clothoid_list ) {
      m_s0.emplace_back(m_s0.back()+C.length());
      m_clothoid_list.emplace_back(C);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( G2solve2arc const & C ) {
    m_s0.reserve( m_s0.size() + 2 );
    m_clothoid_list.reserve( m_clothoid_list.size() + 2 );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    m_s0.emplace_back(m_s0.back()+C.S0().length());
    m_clothoid_list.emplace_back(C.S0());

    m_s0.emplace_back(m_s0.back()+C.S1().length());
    m_clothoid_list.emplace_back(C.S1());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( G2solve3arc const & C ) {
    m_s0.reserve( m_s0.size() + 3 );
    m_clothoid_list.reserve( m_clothoid_list.size() + 3 );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    m_s0.emplace_back(m_s0.back()+C.S0().length());
    m_clothoid_list.emplace_back(C.S0());

    m_s0.emplace_back(m_s0.back()+C.SM().length());
    m_clothoid_list.emplace_back(C.SM());

    m_s0.emplace_back(m_s0.back()+C.S1().length());
    m_clothoid_list.emplace_back(C.S1());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( G2solveCLC const & C ) {
    m_s0.reserve( m_s0.size() + 3 );
    m_clothoid_list.reserve( m_clothoid_list.size() + 3 );

    if ( m_s0.empty() ) m_s0.emplace_back(0);

    m_s0.emplace_back(m_s0.back()+C.S0().length());
    m_clothoid_list.emplace_back(C.S0());

    m_s0.emplace_back(m_s0.back()+C.SM().length());
    m_clothoid_list.emplace_back(C.SM());

    m_s0.emplace_back(m_s0.back()+C.S1().length());
    m_clothoid_list.emplace_back(C.S1());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( Dubins const & c ) {
    if ( m_clothoid_list.empty() ) m_s0.emplace_back( 0 );
    CircleArc const & C0 = c.C0();
    CircleArc const & C1 = c.C1();
    CircleArc const & C2 = c.C2();
    m_s0.emplace_back( m_s0.back()+C0.length() );
    m_s0.emplace_back( m_s0.back()+C1.length() );
    m_s0.emplace_back( m_s0.back()+C2.length() );
    m_clothoid_list.emplace_back( C0 );
    m_clothoid_list.emplace_back( C1 );
    m_clothoid_list.emplace_back( C2 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( Dubins3p const & c3p ) {
    if ( m_clothoid_list.empty() ) m_s0.emplace_back( 0 );
    CircleArc const & C0 = c3p.C0();
    CircleArc const & C1 = c3p.C1();
    CircleArc const & C2 = c3p.C2();
    CircleArc const & C3 = c3p.C3();
    CircleArc const & C4 = c3p.C4();
    CircleArc const & C5 = c3p.C5();
    m_s0.emplace_back( m_s0.back()+C0.length() );
    m_s0.emplace_back( m_s0.back()+C1.length() );
    m_s0.emplace_back( m_s0.back()+C2.length() );
    m_s0.emplace_back( m_s0.back()+C3.length() );
    m_s0.emplace_back( m_s0.back()+C4.length() );
    m_s0.emplace_back( m_s0.back()+C5.length() );
    m_clothoid_list.emplace_back( C0 );
    m_clothoid_list.emplace_back( C1 );
    m_clothoid_list.emplace_back( C2 );
    m_clothoid_list.emplace_back( C3 );
    m_clothoid_list.emplace_back( C4 );
    m_clothoid_list.emplace_back( C5 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type kappa0,
    real_type dkappa,
    real_type L
  ) {
    UTILS_ASSERT0( !m_clothoid_list.empty(), "ClothoidList::push_back_G1(...) empty list!\n" );
    ClothoidCurve c{"ClothoidList::push_back temporary c"};
    real_type x0     = m_clothoid_list.back().x_end();
    real_type y0     = m_clothoid_list.back().y_end();
    real_type theta0 = m_clothoid_list.back().theta_end();
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type kappa0,
    real_type dkappa,
    real_type L
  ) {
    ClothoidCurve c{"ClothoidList::push_back temporary c"};
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    push_back( c );
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    UTILS_ASSERT0( !m_clothoid_list.empty(), "ClothoidList::push_back_G1(...) empty list!\n" );
    ClothoidCurve c{"ClothoidList::push_back_G1 temporary c"};
    real_type x0     = m_clothoid_list.back().x_end();
    real_type y0     = m_clothoid_list.back().y_end();
    real_type theta0 = m_clothoid_list.back().theta_end();
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1
  ) {
    ClothoidCurve c{"ClothoidList::push_back_G1 temporary c"};
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    this->push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1(
    integer         n,
    real_type const x[],
    real_type const y[]
  ) {
    UTILS_ASSERT0( n > 1, "ClothoidList::build_G1, at least 2 points are necessary\n" );

    ClothoidCurve c{"ClothoidList::build_G1 temporary c"};

    init();
    reserve( n-1 );

    if ( n == 2 ) {

      real_type theta = atan2( y[1] - y[0], x[1] - x[0] );
      c.build_G1( x[0], y[0], theta, x[1], y[1], theta );
      this->push_back(c);

    } else {

      Biarc b("build_guess_theta temporary b");
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10;
      real_type thetaC(0);
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] );
        UTILS_ASSERT0( ok, "ClothoidList::build_G1, failed\n" );
        thetaC = b.theta_middle();
      }
      ok = b.build_3P( x[0], y[0], x[1], y[1], x[2], y[2] );
      UTILS_ASSERT0( ok, "ClothoidList::build_G1, failed\n" );
      real_type theta0 = ciclic ? thetaC : b.theta_begin();
      real_type theta1 = b.theta_middle();
      c.build_G1( x[0], y[0], theta0, x[1], y[1], theta1 );
      this->push_back(c);
      for ( integer k = 2; k < n-1; ++k ) {
        theta0 = theta1;
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] );
        UTILS_ASSERT0( ok, "ClothoidList::build_G1, failed\n" );
        theta1 = b.theta_middle();
        c.build_G1( x[k-1], y[k-1], theta0, x[k], y[k], theta1 );
        this->push_back(c);
      }
      theta0 = theta1;
      theta1 = ciclic ? thetaC : b.theta_end();
      c.build_G1( x[n-2], y[n-2], theta0, x[n-1], y[n-1], theta1 );
      this->push_back(c);
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1(
    integer         n,
    real_type const x[],
    real_type const y[],
    real_type const theta[]
  ) {
    UTILS_ASSERT0( n > 1, "ClothoidList::build_G1, at least 2 points are necessary\n" );

    ClothoidCurve c{"ClothoidList::build_G1 temporary c"};

    init();
    reserve( n-1 );

    for ( integer k = 1; k < n; ++k ) {
      c.build_G1( x[k-1], y[k-1], theta[k-1], x[k], y[k], theta[k] );
      this->push_back(c);
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::smooth_quasi_G2(
    integer     max_iter,
    real_type   epsi,
    real_type & max_dK
  ) {

    integer n{integer(m_clothoid_list.size())};

    ClothoidCurve CL("smooth_quasi_G2-left"), CR("smooth_quasi_G2-right");

    auto do_smooth = []( ClothoidCurve * CL, ClothoidCurve * CR ) -> real_type {
      real_type theta0 = CL->theta_begin();
      real_type x0     = CL->x_begin();
      real_type y0     = CL->y_begin();
      real_type theta1 = CL->theta_end();
      real_type x1     = CL->x_end();
      real_type y1     = CL->y_end();
      real_type theta2 = CR->theta_end();
      real_type x2     = CR->x_end();
      real_type y2     = CR->y_end();
      real_type dK     = CL->kappa_end()-CR->kappa_begin();
      auto fun = [CL,CR,x0,y0,theta0,x1,y1,x2,y2,theta2]( real_type theta ) {
        CL->build_G1( x0, y0, theta0, x1, y1, theta  );
        CR->build_G1( x1, y1, theta,  x2, y2, theta2 );
        return CL->kappa_end()-CR->kappa_begin();
      };
      Utils::Algo748<real_type> solver;
      real_type theta = solver.eval2( theta1-Utils::m_pi/20,theta1+Utils::m_pi/20,-Utils::m_pi/2,Utils::m_pi/2,fun);
      fun(theta);
      return std::abs(dK);
    };

    for ( integer iter{0}; iter < max_iter; ++iter ) {
      max_dK = 0;
      for ( integer i{1}; i < n; i +=2 ) {
        real_type dK = do_smooth( &m_clothoid_list[i-1], &m_clothoid_list[i] );
        if ( dK > max_dK ) max_dK = dK;
      }
      for ( integer i{2}; i < n; i +=2 ) {
        real_type dK = do_smooth( &m_clothoid_list[i-1], &m_clothoid_list[i] );
        if ( dK > max_dK ) max_dK = dK;
      }
      if ( m_curve_is_closed ) {
        real_type dK = do_smooth( &m_clothoid_list.front(), &m_clothoid_list.back() );
        if ( dK > max_dK ) max_dK = dK;
      }
      if ( max_dK < epsi ) break;
    }
    return max_dK < epsi;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build(
    real_type       x0,
    real_type       y0,
    real_type       theta0,
    integer         n,
    real_type const s[],
    real_type const kappa[]
  ) {
    if ( n < 2 ) return false;
    real_type tol = abs(s[n-1]-s[0])*machepsi10; // minimum admissible length

    init();
    real_type k  = kappa[0];
    real_type L  = s[1]-s[0];
    real_type dk = (kappa[1]-k)/L;
    UTILS_ASSERT(
      Utils::is_finite( k ) && Utils::is_finite( L ) && Utils::is_finite( dk ),
      "ClothoidList::build, failed first segment found\n"
      "L = {} k = {} dk = {}\n",
      L, k, dk
    );
    this->push_back( x0, y0, theta0, k, dk, L );
    for ( integer i = 2; i < n; ++i ) {
      k  = kappa[i-1];
      L  = s[i]-s[i-1];
      if ( abs(L) < tol ) {
        fmt::print( "ClothoidList::build, skipping segment N.{}\n", i);
        continue; // skip too small segment
      }
      dk = (kappa[i]-k)/L;
      UTILS_ASSERT(
        Utils::is_finite( k ) && Utils::is_finite( L ) && Utils::is_finite( dk ),
        "ClothoidList::build, failed at segment N.{} found\n"
        "L = {} k = {} dk = {}\n",
        i, L, k, dk
      );
      this->push_back( k, dk, L );
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_raw(
    integer         n,
    real_type const x[],
    real_type const y[],
    real_type const abscissa[],
    real_type const theta[],
    real_type const kappa[]
  ) {
    if ( n < 2 ) return false;
    init();
    m_clothoid_list.reserve(size_t(n-1));
    real_type const * px = x;
    real_type const * py = y;
    real_type const * pa = abscissa;
    real_type const * pt = theta;
    real_type const * pk = kappa;
    for ( integer i = 1; i < n; ++i, ++px, ++py, ++pa, ++pt, ++pk ) {
      real_type dk = pk[1]-pk[0];
      real_type L  = pa[1]-pa[0];
      this->push_back( *px, *py, *pt, *pk, dk, L );
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::get( integer idx ) const {
    UTILS_ASSERT(
      !m_clothoid_list.empty(), "ClothoidList::get( {} ) empty list\n", idx
    );
    try {
      return m_clothoid_list.at(idx);
    } catch ( std::exception & exc ) {
      UTILS_ERROR( "ClothoidList::get( {} ): {}\n", idx, exc.what() );
    } catch ( ... ) {
      UTILS_ERROR( "ClothoidList::get( {} ): unknown error\n", idx );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::get_at_s( real_type s ) const {
    integer idx = this->find_at_s(s);
    return get( idx );
  }

  /*\
   |   _                  _   _
   |  | | ___ _ __   __ _| |_| |__
   |  | |/ _ \ '_ \ / _` | __| '_ \
   |  | |  __/ | | | (_| | |_| | | |
   |  |_|\___|_| |_|\__, |\__|_| |_|
   |                |___/
  \*/

  real_type
  ClothoidList::length() const {
    if ( m_clothoid_list.empty() ) return real_type(0);
    return m_s0.back() - m_s0.front();
  }

  real_type
  ClothoidList::length_ISO( real_type offs ) const {
    real_type L{0};
    for ( ClothoidCurve const & C : m_clothoid_list ) L += C.length_ISO( offs );
    return L;
  }

  real_type
  ClothoidList::segment_length( integer nseg ) const {
    ClothoidCurve const & c = get( nseg );
    return c.length();
  }

  real_type
  ClothoidList::segment_length_ISO( integer nseg, real_type offs ) const {
    ClothoidCurve const & c = get( nseg );
    return c.length_ISO( offs );
  }

  /*\
   |  _    _   _____    _                _
   | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
   | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
   | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
   |                               |___/
  \*/

  void
  ClothoidList::bb_triangles(
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    for ( integer ipos = icurve; ic != m_clothoid_list.end(); ++ic, ++ipos )
      ic->bb_triangles( tvec, max_angle, max_size, ipos );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidList::bb_triangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    for ( integer ipos = icurve; ic != m_clothoid_list.end(); ++ic, ++ipos )
      ic->bb_triangles_ISO( offs, tvec, max_angle, max_size, ipos );
  }

  /*\
   |   _     _
   |  | |__ | |__   _____  __
   |  | '_ \| '_ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |_.__/|_.__/ \___/_/\_\
  \*/

  void
  ClothoidList::bbox_ISO(
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
    for ( auto const & T : tvec ) {
      // - - - - - - - - - - - - - - - - - - - -
      if      ( T.x1() < xmin ) xmin = T.x1();
      else if ( T.x1() > xmax ) xmax = T.x1();
      if      ( T.x2() < xmin ) xmin = T.x2();
      else if ( T.x2() > xmax ) xmax = T.x2();
      if      ( T.x3() < xmin ) xmin = T.x3();
      else if ( T.x3() > xmax ) xmax = T.x3();
      // - - - - - - - - - - - - - - - - - - - -
      if      ( T.y1() < ymin ) ymin = T.y1();
      else if ( T.y1() > ymax ) ymax = T.y1();
      if      ( T.y2() < ymin ) ymin = T.y2();
      else if ( T.y2() > ymax ) ymax = T.y2();
      if      ( T.y3() < ymin ) ymin = T.y3();
      else if ( T.y3() > ymax ) ymax = T.y3();
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
  ClothoidList::theta( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.theta( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_D( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_DD( real_type s ) const  {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_DDD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.theta_DDD( s - m_s0[idx] );
  }

  /*\
   |  _____                   _    _   _
   | |_   _|   __ _ _ __   __| |  | \ | |
   |   | |    / _` | '_ \ / _` |  |  \| |
   |   | |   | (_| | | | | (_| |  | |\  |
   |   |_|    \__,_|_| |_|\__,_|  |_| \_|
  \*/

  real_type
  ClothoidList::tx( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tx( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.ty( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_D( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_D( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_DD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_DD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::tx_DDD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tx_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::ty_DDD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.ty_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg(
    real_type   s,
    real_type & tg_x,
    real_type & tg_y
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tg( s - m_s0[idx], tg_x, tg_y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_D(
    real_type   s,
    real_type & tg_x_D,
    real_type & tg_y_D
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_D( s - m_s0[idx], tg_x_D, tg_y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_DD(
    real_type   s,
    real_type & tg_x_DD,
    real_type & tg_y_DD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_DD( s - m_s0[idx], tg_x_DD, tg_y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::tg_DDD(
    real_type   s,
    real_type & tg_x_DDD,
    real_type & tg_y_DDD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.tg_DDD( s - m_s0[idx], tg_x_DDD, tg_y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::evaluate(
    real_type   s,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    c.evaluate( s - m_s0[idx], th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::evaluate_ISO(
    real_type   s,
    real_type   offs,
    real_type & th,
    real_type & k,
    real_type & x,
    real_type & y
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    c.evaluate_ISO( s - m_s0[idx], offs, th, k, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y( real_type s ) const  {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_D( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_D( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_D( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_DD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_DD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_DD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_DDD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_DDD( real_type s ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_DDD( s - m_s0[idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval(
    real_type   s,
    real_type & x,
    real_type & y
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval( s - m_s0[idx], x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D(
    real_type   s,
    real_type & x_D,
    real_type & y_D
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_D( s - m_s0[idx], x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD(
    real_type   s,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_DD( s - m_s0[idx], x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD(
    real_type   s,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
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
  ClothoidList::X_ISO( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_D( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_D( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_D( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_D( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_DD( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_DD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_DD( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_DD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X_ISO_DDD( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.X_ISO_DDD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y_ISO_DDD( real_type s, real_type offs ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.Y_ISO_DDD( s - m_s0[idx], offs );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO(
    real_type   s,
    real_type   offs,
    real_type & x,
    real_type & y
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO( s - m_s0[idx], offs, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_D(
    real_type   s,
    real_type   offs,
    real_type & x_D,
    real_type & y_D
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO_D( s - m_s0[idx], offs, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_DD(
    real_type   s,
    real_type   offs,
    real_type & x_DD,
    real_type & y_DD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
    return c.eval_ISO_DD( s - m_s0[idx], offs, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_ISO_DDD(
    real_type   s,
    real_type   offs,
    real_type & x_DDD,
    real_type & y_DDD
  ) const {
    integer idx = find_at_s( s );
    ClothoidCurve const & c = get( idx );
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
  ClothoidList::translate( real_type tx, real_type ty ) {
    for ( ClothoidCurve & C : m_clothoid_list ) C.translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::rotate( real_type angle, real_type cx, real_type cy ) {
    for ( ClothoidCurve & C : m_clothoid_list ) C.rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::scale( real_type sfactor ) {
    vector<ClothoidCurve>::iterator ic = m_clothoid_list.begin();
    real_type newx0 = ic->x_begin();
    real_type newy0 = ic->y_begin();
    m_s0[0] = 0;
    for ( size_t k=0; ic != m_clothoid_list.end(); ++ic, ++k ) {
      ic->scale( sfactor );
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reverse() {
    std::reverse( m_clothoid_list.begin(), m_clothoid_list.end() );
    vector<ClothoidCurve>::iterator ic = m_clothoid_list.begin();
    ic->reverse();
    real_type newx0 = ic->x_end();
    real_type newy0 = ic->y_end();
    m_s0[0] = 0;
    m_s0[1] = ic->length();
    size_t k = 1;
    for ( ++ic; ic != m_clothoid_list.end(); ++ic, ++k ) {
      ic->reverse();
      ic->change_origin( newx0, newy0 );
      newx0     = ic->x_end();
      newy0     = ic->y_end();
      m_s0[k+1] = m_s0[k] + ic->length();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::change_origin( real_type newx0, real_type newy0 ) {
    for ( ClothoidCurve & C : m_clothoid_list ) {
      C.change_origin( newx0, newy0 );
      newx0 = C.x_end();
      newy0 = C.y_end();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::trim( real_type s_begin, real_type s_end ) {

    ClothoidList newCL{"ClothoidList::trim temporary newCL"};
    this->trim( s_begin, s_end, newCL );
    this->copy( newCL );

#if 0
    UTILS_ASSERT(
      s_begin >= m_s0.front() && s_end <= m_s0.back() && s_end > s_begin,
      "ClothoidList::trim( s_begin={}, s_end={} ) bad range, must be in [{},{}]\n",
      s_begin, s_end, m_s0.front(), m_s0.back()
    );

    size_t i_begin = size_t( find_at_s( s_begin ) );
    size_t i_end   = size_t( find_at_s( s_end   ) );
    if ( i_begin == i_end ) {
      m_clothoid_list[i_begin].trim( s_begin-m_s0[i_begin], s_end-m_s0[i_begin] );
    } else {
      m_clothoid_list[i_begin].trim( s_begin-m_s0[i_begin], m_s0[i_begin+1]-m_s0[i_begin] );
      m_clothoid_list[i_end].trim( 0, s_end-m_s0[i_end] );
    }
    m_clothoid_list.erase( m_clothoid_list.begin()+i_end+1, m_clothoid_list.end() );
    m_clothoid_list.erase( m_clothoid_list.begin(), m_clothoid_list.begin()+i_begin );
    if ( m_clothoid_list.back().m_L <= machepsi100 ) m_clothoid_list.pop_back();
    vector<ClothoidCurve>::iterator ic = m_clothoid_list.begin();
    m_s0.resize( m_clothoid_list.size() + 1 );
    m_s0[0] = 0;
    size_t k{0};
    for (; ic != m_clothoid_list.end(); ++ic, ++k )
      m_s0[k+1] = m_s0[k] + ic->length();
    this->reset_last_interval();
#endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::trim( real_type s_begin, real_type s_end, ClothoidList & newCL ) const {

    newCL.init();

    if ( m_clothoid_list.empty() ) return;

    // put in range
    real_type L = this->length();
    while ( s_begin > L ) s_begin -= L;
    while ( s_begin < 0 ) s_begin += L;
    while ( s_end   > L ) s_end   -= L;
    while ( s_end   < 0 ) s_end   += L;

    integer n_seg   = integer( m_clothoid_list.size() );
    integer i_begin = find_at_s( s_begin );
    integer i_end   = find_at_s( s_end );

    if ( s_begin < s_end ) {
      // get initial and final segment
      if ( i_begin == i_end ) { // stesso segmento
        real_type   ss0 = m_s0[i_begin];
        ClothoidCurve C = m_clothoid_list[i_begin];
        C.trim( s_begin-ss0, s_end-ss0 );
        newCL.push_back( C );
      } else {
        ClothoidCurve C0 = m_clothoid_list[i_begin];
        C0.trim( s_begin - m_s0[i_begin], C0.length() );
        newCL.push_back( C0 );

        for ( ++i_begin; i_begin < i_end; ++i_begin )
          newCL.push_back( m_clothoid_list[i_begin] );

        ClothoidCurve C1 = m_clothoid_list[i_end];
        C1.trim( 0, s_end - m_s0[i_end] );
        newCL.push_back( C1 );
      }
    } else {
      ClothoidCurve C0 = m_clothoid_list[i_begin];
      C0.trim( s_begin - m_s0[i_begin], C0.length() );
      newCL.push_back( C0 );

      for ( ++i_begin; i_begin < n_seg; ++i_begin )
        newCL.push_back( m_clothoid_list[i_begin] );

      for ( i_begin = 0; i_begin < i_end; ++i_begin )
        newCL.push_back( m_clothoid_list[i_begin] );

      ClothoidCurve C1 = m_clothoid_list[i_end];
      C1.trim( 0, s_end - m_s0[i_end] );
      newCL.push_back( C1 );
    }
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
  ClothoidList::build_AABBtree_ISO(
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
  ClothoidList::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::CLOTHOID_LIST ) {
      ClothoidList const & C = *static_cast<ClothoidList const *>(pC);
      return this->collision( C );
    } else {
      ClothoidList C(pC);
      return this->collision( C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::CLOTHOID_LIST ) {
      ClothoidList const & C = *static_cast<ClothoidList const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      ClothoidList C(pC);
      return this->collision_ISO( offs, C, offs_C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::collision_ISO(
    real_type            offs,
    ClothoidList const & CL,
    real_type            offs_CL
  ) const {

    G2LIB_DEBUG_TICTOC;

    G2LIB_DEBUG_TIC;
    this->build_AABBtree_ISO( offs );
    CL.build_AABBtree_ISO( offs_CL );
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::collision_ISO( offs={}, CL, offs_L={} ) build_AABBtree_ISO elapsed {}ms\n",
      offs, offs_CL, tictoc.elapsed_ms()
    );

    G2LIB_DEBUG_TIC;
    AABB_MAP intersectList;
    m_aabb_tree.intersect_and_refine( CL.m_aabb_tree, intersectList );
    G2LIB_DEBUG_TOC;
    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::collision_ISO intersect_and_refine elapsed {}ms candidated #{}\n",
      tictoc.elapsed_ms(), intersectList.size()
    );

    G2LIB_DEBUG_TIC;
    bool collide = false;
    for ( auto const & I : intersectList ) {
      integer i = I.first;
      UTILS_ASSERT_DEBUG(
        i >= 0 && i < integer(m_aabb_triangles.size()),
        "ClothoidList::collision_ISO( offs={}, C, offs_CL={} ) i={} out of range [0,{})\n",
        offs, offs_CL, i, m_aabb_triangles.size()
      );
      Triangle2D    const & T1 = m_aabb_triangles.at(i);
      ClothoidCurve const & C1 = m_clothoid_list.at(T1.Icurve());
      for ( integer j : I.second ) {
        UTILS_ASSERT_DEBUG(
          j >= 0 && j < integer(CL.m_aabb_triangles.size()),
          "ClothoidList::collision_ISO( offs={}, CL, offs_CL={} ) j={} out of range [0,{})\n",
          offs, offs_CL, j, CL.m_aabb_triangles.size()
        );
        Triangle2D    const & T2 = CL.m_aabb_triangles.at(j);
        ClothoidCurve const & C2 = CL.m_clothoid_list.at(T2.Icurve());
        collide = C1.collision_ISO( offs, C2, offs_CL );
        if ( collide ) break;
      }
      if ( collide ) break;
    }

    G2LIB_DEBUG_TOC;
    G2LIB_DEBUG_MESSAGE(
      "ClothoidList collision_ISO: collisions elapsed {}ms\n",
      tictoc.elapsed_ms()
    );

    return collide;
  }

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
   |
  \*/


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::intersect_ISO(
    real_type            offs,
    ClothoidList const & CL,
    real_type            offs_CL,
    IntersectList      & ilist
  ) const {

    G2LIB_DEBUG_TICTOC;

    if ( intersect_with_AABBtree ) {

      G2LIB_DEBUG_TIC;
      this->build_AABBtree_ISO( offs );
      CL.build_AABBtree_ISO( offs_CL );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::intersect_ISO( offs={}, CL, offs_CL={}, ilist ) build_AABBtree_ISO elapsed {}ms\n",
        offs, offs_CL, tictoc.elapsed_ms()
      );

      G2LIB_DEBUG_TIC;
      AABB_MAP intersectList;
      m_aabb_tree.intersect_and_refine( CL.m_aabb_tree, intersectList );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::intersect_ISO intersect_and_refine elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), intersectList.size()
      );

      G2LIB_DEBUG_TIC;
      for ( auto const & I : intersectList ) {
        integer i = I.first;
        UTILS_ASSERT_DEBUG(
          i >= 0 && i < integer(m_aabb_triangles.size()),
          "ClothoidList::intersect_ISO( offs={}, CL, offs_CL={}, ilist ) i={} out of range [0,{})\n",
          offs, offs_CL, i, m_aabb_triangles.size()
        );
        Triangle2D    const & T1 = m_aabb_triangles.at(i);
        ClothoidCurve const & C1 = m_clothoid_list.at(T1.Icurve());

        for ( integer j : I.second ) {
          UTILS_ASSERT_DEBUG(
            j >= 0 && j < integer(CL.m_aabb_triangles.size()),
            "ClothoidList::intersect_ISO( offs={}, CL, offs_CL={}, ilist ) j={} out of range [0,{})\n",
            offs, offs_CL, j, CL.m_aabb_triangles.size()
          );
          Triangle2D    const & T2 = CL.m_aabb_triangles.at(j);
          ClothoidCurve const & C2 = CL.m_clothoid_list.at(T2.Icurve());

          real_type ss1, ss2;
          bool converged = C1.aabb_intersect_ISO( T1, offs, &C2, T2, offs_CL, ss1, ss2 );

          if ( converged ) {
            ss1 += m_s0.at(T1.Icurve());
            ss2 += CL.m_s0.at(T2.Icurve());
            ilist.emplace_back( ss1, ss2 );
          }
        }
      }

      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::intersect_ISO: intersect elapsed {}ms\n",
        tictoc.elapsed_ms()
      );

    } else {

      G2LIB_DEBUG_TIC;
      bb_triangles_ISO( offs, m_aabb_triangles, Utils::m_pi/18, 1e100 );
      CL.bb_triangles_ISO( offs_CL, CL.m_aabb_triangles, Utils::m_pi/18, 1e100 );

      for ( Triangle2D const & T1 : m_aabb_triangles ) {
        for ( Triangle2D const & T2 : CL.m_aabb_triangles ) {

          ClothoidCurve const & C1 = m_clothoid_list.at(T1.Icurve());
          ClothoidCurve const & C2 = CL.m_clothoid_list.at(T2.Icurve());

          real_type ss1, ss2;
          bool converged = C1.aabb_intersect_ISO( T1, offs, &C2, T2, offs_CL, ss1, ss2 );

          if ( converged ) {
            ss1 += m_s0.at(T1.Icurve());
            ss2 += CL.m_s0.at(T2.Icurve());
            ilist.emplace_back( ss1, ss2 );
          }
        }
      }

      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve::intersect_ISO: collisions elapsed {}ms noAABB\n",
        tictoc.elapsed_ms()
      );

    }
  }

  void
  ClothoidList::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CLOTHOID_LIST ) {
      ClothoidList const & C = *static_cast<ClothoidList const *>(pC);
      this->intersect( C, ilist );
    } else {
      ClothoidList C(pC);
      this->intersect( C, ilist );
    }
  }

  void
  ClothoidList::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CLOTHOID_LIST ) {
      ClothoidList const & C = *static_cast<ClothoidList const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      ClothoidList C(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
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
  ClothoidList::closest_point_internal(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & DST
  ) const {
    G2LIB_DEBUG_TICTOC;

    G2LIB_DEBUG_TIC;
    this->build_AABBtree_ISO( offs );
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::closest_point_internal build_AABBtree_ISO elapsed {}ms\n",
      tictoc.elapsed_ms()
    );

    integer icurve{0};
    DST = Utils::Inf<real_type>();

    if ( m_aabb_tree.num_tree_nodes() > G2LIB_AABB_MIN_NODES && intersect_with_AABBtree ) {

      G2LIB_DEBUG_TIC;
      AABB_SET candidateList;
      real_type xy[2] = { qx, qy };
      m_aabb_tree.min_distance_candidates( xy, candidateList );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_internal min_distance_candidates elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), candidateList.size()
      );

      UTILS_ASSERT0(
         candidateList.size() > 0,
        "ClothoidList::closest_point_internal no candidate\n"
      );

      G2LIB_DEBUG_TIC;
      for ( integer ipos : candidateList ) {
        Triangle2D const & T = m_aabb_triangles.at(ipos);
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss;
          ClothoidCurve const & C = m_clothoid_list.at(T.Icurve());
          C.closest_point_internal( T.S0(), T.S1(), qx, qy, offs, xx, yy, ss, dst );
          if ( dst < DST ) {
            DST    = dst;
            s      = ss + m_s0.at(T.Icurve());
            x      = xx;
            y      = yy;
            icurve = T.Icurve();
          }
        }
      }
      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_internal minimize elapsed {}ms\n",
        tictoc.elapsed_ms()
      );

    } else {

      G2LIB_DEBUG_TIC;
      for ( Triangle2D const & T : m_aabb_triangles ) {
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss, tt;
          m_clothoid_list.at(T.Icurve()).closest_point_ISO( qx, qy, offs, xx, yy, ss, tt, dst );
          if ( dst < DST ) {
            DST    = dst;
            s      = ss + m_s0.at(T.Icurve());
            x      = xx;
            y      = yy;
            icurve = T.Icurve();
          }
        }
      }
      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_internal minimize elapsed {}ms noAABB\n",
        tictoc.elapsed_ms()
      );

    }
    return icurve;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::closest_point_ISO(
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
    UTILS_ASSERT(
      icurve >= 0 && icurve < integer(m_clothoid_list.size()),
      "ClothoidList::closest_point_ISO\n"
      "call to closest_point_internal return icurve = {}\n"
      "icurve must be in [0,{})\n",
      icurve, m_clothoid_list.size()
    );

    // check if projection is orthogonal
    real_type nx, ny;
    m_clothoid_list.at(icurve).nor_ISO( s - m_s0.at(icurve), nx, ny );
    real_type qxx = qx - x;
    real_type qyy = qy - y;
    t = qxx * nx + qyy * ny - offs; // signed distance
    real_type pt = abs(qxx * ny - qyy * nx);
    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::closest_point_ISO\n"
      "||P-P0|| = {} and {}, |(P-P0).T| = {}\n",
      DST, hypot(qxx,qyy), pt
    );
    return pt > GLIB2_TOL_ANGLE*hypot(qxx,qyy) ? -(icurve+1) : icurve;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::closest_point_ISO(
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

  /*\
   |      _ _     _
   |   __| (_)___| |_ __ _ _ __   ___ ___
   |  / _` | / __| __/ _` | '_ \ / __/ _ \
   | | (_| | \__ \ || (_| | | | | (_|  __/
   |  \__,_|_|___/\__\__,_|_| |_|\___\___|
  \*/

  integer
  ClothoidList::closest_segment( real_type qx, real_type qy ) const {

    G2LIB_DEBUG_TICTOC;

    G2LIB_DEBUG_TIC;
    this->build_AABBtree_ISO( 0 );
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::closest_segment build_AABBtree_ISO elapsed {}ms\n",
      tictoc.elapsed_ms()
    );

    integer   icurve{0};
    real_type DST{ Utils::Inf<real_type>() };

    if ( m_aabb_tree.num_tree_nodes() > G2LIB_AABB_MIN_NODES && intersect_with_AABBtree ) {

      G2LIB_DEBUG_TIC;
      AABB_SET candidateList;
      real_type xy[2] = { qx, qy };
      m_aabb_tree.min_distance_candidates( xy, candidateList );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_segment min_distance_candidates elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), candidateList.size()
      );

      UTILS_ASSERT0(
         candidateList.size() > 0,
        "ClothoidList::closest_segment no candidate\n"
      );

      G2LIB_DEBUG_TIC;
      for ( integer ipos : candidateList ) {
        Triangle2D const & T = m_aabb_triangles.at(ipos);
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss, tt;
          m_clothoid_list.at(T.Icurve()).closest_point_ISO( qx, qy, 0, xx, yy, ss, tt, dst );
          if ( dst < DST ) {
            DST    = dst;
            icurve = T.Icurve();
          }
        }
      }
      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_segment elapsed {}ms\n",
        tictoc.elapsed_ms()
      );
    } else {

      G2LIB_DEBUG_TIC;
      for ( Triangle2D const & T : m_aabb_triangles ) {
        real_type dst = T.dist_min( qx, qy ); // distanza approssimata con triangolo
        if ( dst < DST ) {
          // refine distance
          real_type xx, yy, ss, tt;
          m_clothoid_list.at(T.Icurve()).closest_point_ISO( qx, qy, 0, xx, yy, ss, tt, dst );
          if ( dst < DST ) {
            DST    = dst;
            icurve = T.Icurve();
          }
        }
      }
      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_segment elapsed {}ms noAABB\n",
        tictoc.elapsed_ms()
      );
    }
    return icurve;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::closest_point_in_range_ISO(
    real_type   qx,
    real_type   qy,
    integer     icurve_begin,
    integer     icurve_end,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst,
    integer   & icurve
  ) const {
    UTILS_ASSERT0(
      !m_clothoid_list.empty(),
      "ClothoidList::closest_point_in_range_ISO, empty list\n"
    );
    integer nsegs = this->num_segments();
    if ( nsegs == 1 ) { // only 1 segment to check
      icurve = 0;
      integer res = m_clothoid_list.front().closest_point_ISO( qx, qy, x, y, s, t, dst );
      s += m_s0[0];
      return res;
    }

    integer ib = icurve_begin % nsegs; // to avoid infinite loop in case of bad input
    integer ie = icurve_end   % nsegs; // to avoid infinite loop in case of bad input
    if ( ib < 0 ) ib += nsegs;
    if ( ie < 0 ) ie += nsegs;
    UTILS_ASSERT(
      ib >= 0 && ie >= 0,
      "ClothoidList::closest_point_in_range_ISO, ib = {} ie = {}\n",
      ib, ie
    );

    icurve = ib;
    integer res = m_clothoid_list.at(icurve).closest_point_ISO( qx, qy, x, y, s, t, dst );
    s += m_s0.at(icurve);

    G2LIB_DEBUG_MESSAGE(
      "ClothoidList::closest_point_in_range_ISO\n"
      "first segment #{} dst = {} res = {}\n",
      icurve, dst, res
    );

    if ( ib == ie ) return res; // only one segment to check

    integer iseg = ib;
    do {
      if ( ++iseg >= nsegs ) iseg -= nsegs; // next segment
      real_type C_x, C_y, C_s, C_t, C_dst;
      integer C_res = m_clothoid_list.at(iseg).closest_point_ISO(
        qx, qy, C_x, C_y, C_s, C_t, C_dst
      );
      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_in_range_ISO: segment #{} dst = {} res = {}\n",
        iseg, C_dst, C_res
      );
      if ( C_dst < dst ) {
        dst    = C_dst;
        x      = C_x;
        y      = C_y;
        s      = C_s + m_s0[iseg];
        t      = C_t;
        icurve = iseg;
        res    = C_res;
        G2LIB_DEBUG_MESSAGE(
          "ClothoidList::closest_point_in_range_ISO, new min at s = {}, res = {}\n", s, res
        );
      }
    } while ( iseg != ie );
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::closest_point_in_s_range_ISO(
    real_type   qx,
    real_type   qy,
    real_type   s_begin,
    real_type   s_end,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & dst,
    integer   & icurve
  ) const {
    UTILS_ASSERT0(
      !m_clothoid_list.empty(),
      "ClothoidList::closest_point_in_s_range_ISO, empty list\n"
    );
    // put in range
    while ( s_begin < 0              ) s_begin += this->length();
    while ( s_begin > this->length() ) s_begin -= this->length();
    while ( s_end   < 0              ) s_end   += this->length();
    while ( s_end   > this->length() ) s_end   -= this->length();

    // get initial and final segment
    integer i_begin = find_at_s( s_begin );
    integer i_end   = find_at_s( s_end );
    integer res     = 0;
    if ( i_begin == i_end ) {
      // stesso segmento
      real_type     ss0 = m_s0[i_begin];
      ClothoidCurve C   = m_clothoid_list[i_begin]; // crea copia
      C.trim( s_begin-ss0, s_end-ss0 );
      res = C.closest_point_ISO( qx, qy, x, y, s, t, dst );
      s += s_begin;
    } else {
      // segmenti consecutivi
      integer   res1;
      real_type x1, y1, s1, t1, dst1;

      real_type     ss0 = m_s0[i_begin];
      real_type     ss1 = m_s0[i_end];
      ClothoidCurve C0  = m_clothoid_list[i_begin]; // crea copia
      ClothoidCurve C1  = m_clothoid_list[i_end];   // crea copia

      // taglia il segmento
      C0.trim( s_begin-ss0, C0.length() );

      // calcolo closest point
      res = C0.closest_point_ISO( qx, qy, x, y, s, t, dst );
      s  += s_begin;
      icurve = i_begin;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_in_s_range_ISO: first segment {} dst = {} res = {}\n",
        i_begin, dst, res
      );

      C1.trim( 0, s_end-ss1 );
      res1 = C1.closest_point_ISO( qx, qy, x1, y1, s1, t1, dst1 );
      s1 += ss1;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidList::closest_point_in_s_range_ISO: last segment {} dst = {} res = {}\n",
        i_end, dst1, res1
      );

      if ( dst1 < dst ) {
        x = x1; y = y1; s = s1; t = t1;
        dst = dst1; res = res1; icurve = i_end;
      }

      // ci sono altri segmenti?
      if ( i_end < i_begin ) i_end += integer(m_clothoid_list.size());
      ++i_begin;
      if ( i_begin < i_end ) {
        integer icurve1;
        res1 = this->closest_point_in_range_ISO(
          qx, qy, i_begin, i_end, x1, y1, s1, t1, dst1, icurve1
        );
        G2LIB_DEBUG_MESSAGE(
          "ClothoidList::closest_point_in_s_range_ISO: range [{},{}] dst = {} res = {}\n",
          i_begin, i_end, dst1, res1
        );
        if ( dst1 < dst ) {
          x      = x1;
          y      = y1;
          s      = s1;
          t      = t1;
          dst    = dst1;
          res    = res1;
          icurve = icurve1;
          G2LIB_DEBUG_MESSAGE(
            "ClothoidList::closest_point_in_s_range_ISO: new min at s = {}, res = {}\n", s, res
          );
        }
      }
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::closest_point_by_sample(
    real_type   ds,
    real_type   qx,
    real_type   qy,
    real_type & X,
    real_type & Y,
    real_type & S
  ) const {
    real_type dst{Utils::Inf<real_type>()};
    integer   i{0};
    for ( ClothoidCurve const & C : m_clothoid_list ) {
      real_type xx, yy, ss;
      real_type dd = C.closest_point_by_sample( ds, qx, qy, xx, yy, ss );
      if ( dd < dst ) {
        dst = dd;
        X   = xx;
        Y   = yy;
        S   = ss+m_s0[i];
      }
      ++i;
    }
    return dst;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::get_SK( real_type s[], real_type kappa[] ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    integer   k{0};
    real_type ss{0};
    while ( ic != m_clothoid_list.end() ) {
      s[k]     = ss;
      kappa[k] = ic->kappa_begin();
      ss      += ic->length();
      ++k;
      ++ic;
    }
    --ic; // last element
    s[k]     = ss;
    kappa[k] = ic->kappa_end();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::get_STK(
    real_type s[],
    real_type theta[],
    real_type kappa[]
  ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    integer   k{0};
    real_type ss{0};
    while ( ic != m_clothoid_list.end() ) {
      s[k]     = ss;
      theta[k] = ic->theta_begin();
      kappa[k] = ic->kappa_begin();
      ss      += ic->length();
      ++k;
      ++ic;
    }
    --ic; // last element
    s[k]     = ss;
    theta[k] = ic->theta_end();
    kappa[k] = ic->kappa_end();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::get_XY( real_type x[], real_type y[] ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    integer k{0};
    while ( ic != m_clothoid_list.end() ) {
      x[k] = ic->x_begin();
      y[k] = ic->y_begin();
      ++k; ++ic;
    }
    --ic;
    x[k] = ic->x_end();
    y[k] = ic->y_end();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::get_delta_theta( real_type delta_theta[] ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    integer k{0};
    for ( ++ic; ic != m_clothoid_list.end(); ++ic, ++k ) {
      real_type tmp = ic->theta_begin()-ic[-1].theta_end();
      if      ( tmp >  Utils::m_pi ) tmp -= Utils::m_2pi;
      else if ( tmp < -Utils::m_pi ) tmp += Utils::m_2pi;
      delta_theta[k] = tmp;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::get_delta_kappa( real_type deltaKappa[] ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    integer k{0};
    for ( ++ic; ic != m_clothoid_list.end(); ++ic, ++k  )
      deltaKappa[k] = ic->kappa_begin()-ic[-1].kappa_end();

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidList::findST1(
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    UTILS_ASSERT0( !m_clothoid_list.empty(), "ClothoidList::findST, empty list\n" );
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    vector<real_type>::const_iterator     is = m_s0.begin();

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
          ic != m_clothoid_list.end();
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
  ClothoidList::findST1(
    integer     ibegin,
    integer     iend,
    real_type   x,
    real_type   y,
    real_type & s,
    real_type & t
  ) const {

    UTILS_ASSERT0(
      !m_clothoid_list.empty(),
      "ClothoidList::findST, empty list\n"
    );
    UTILS_ASSERT(
      ibegin >= 0 && ibegin <= iend && iend < integer(m_clothoid_list.size()),
      "ClothoidList::findST( ibegin={}, iend={}, x, y, s, t ) bad range not in [0,{}]\n",
      ibegin, iend, m_clothoid_list.size()-1
    );
    s = t = 0;
    integer iseg{0};
    bool ok = false;
    for ( integer k = ibegin; k <= iend; ++k ) {
      ClothoidCurve const & ck = m_clothoid_list[k];
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

  void
  ClothoidList::export_table( ostream_type & stream ) const {
    stream << "x\ty\ttheta0\tkappa0\tdkappa\tL\n";
    for ( auto const & C : m_clothoid_list )
      fmt::print( stream,
        "{}\t{}\t{}\t{}\t{}\t{}\n",
        C.x_begin(),
        C.y_begin(),
        C.theta_begin(),
        C.kappa_begin(),
        C.dkappa(),
        C.length()
      );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_ruby( ostream_type & stream ) const {
    stream << "data = {\n";
    for ( auto const & C : m_clothoid_list )
      fmt::print( stream,
        "{}\t{}\t{}\t{}\t{}\t{}\n",
        C.x_begin(),
        C.y_begin(),
        C.theta_begin(),
        C.kappa_begin(),
        C.dkappa(),
        C.length()
      );
    stream << "}\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `ClothoidList` object
  //!
  //!  \param stream the output stream
  //!  \param CL     an instance of `ClothoidList` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, ClothoidList const & CL ) {
    for ( auto const & C : CL.m_clothoid_list )
      stream << C << '\n';
    return stream;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static
  void
  save_segment( ostream_type & stream, ClothoidCurve const & c ) {
    fmt::print( stream,
      "{:<24}\t{:<24}\t{:<24}\t{:<24}\n"
      "{:<24}\t{:<24}\t{:<24}\t{:<24}\n",
      //------------------
      fmt::format("{:.20}",c.x_begin()),
      fmt::format("{:.20}",c.y_begin()),
      fmt::format("{:.20}",c.theta_begin()),
      fmt::format("{:.20}",c.kappa_begin()),
      //------------------
      fmt::format("{:.20}",c.x_end()),
      fmt::format("{:.20}",c.y_end()),
      fmt::format("{:.20}",c.theta_end()),
      fmt::format("{:.20}",c.kappa_end())
    );
  }
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static
  bool
  load_segment(
    istream_type  & stream,
    ClothoidCurve & c,
    real_type       epsi
  ) {
    string line1, line2;
    while ( stream.good() ) {
      if ( !getline(stream,line1) ) return false;
      if ( line1[0] != '#' ) break;
    }
    if ( !stream.good() ) return false;
    while ( stream.good() ) {
      if ( !getline(stream,line2) ) return false;
      if ( line2[0] != '#' ) break;
    }
    if ( !stream.good() ) return false;
    std::istringstream iss1(line1);
    std::istringstream iss2(line2);
    real_type x0, y0, x1, y1, theta0, theta1, kappa0, kappa1;
    iss1 >> x0 >> y0 >> theta0 >> kappa0;
    iss2 >> x1 >> y1 >> theta1 >> kappa1;
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    // check segment
    real_type err1 = std::abs( kappa0 - c.kappa_begin() ) * c.length();
    real_type err2 = std::abs( kappa1 - c.kappa_end() ) * c.length();
    UTILS_ASSERT(
      err1 < epsi && err2 < epsi,
      "load_segment, failed tolerance on curvature\n"
      "begin error = {}, end error = {}\n",
      err1, err2
    );
    return true;
  }
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::save( ostream_type & stream ) const {
    vector<ClothoidCurve>::const_iterator ic = m_clothoid_list.begin();
    stream << "# x y theta kappa\n";
    for ( integer nseg = 1; ic != m_clothoid_list.end(); ++ic, ++nseg ) {
      stream << "# segment n." << nseg << '\n';
      save_segment( stream, *ic );
    }
    stream << "# EOF\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::load( istream_type & stream, real_type epsi ) {
    this->init();
    while ( stream.good() ) {
      ClothoidCurve c{"ClothoidList::load temporary c"};
      bool ok = load_segment( stream, c, epsi );
      if ( !ok ) break;
      this->push_back( c );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  ClothoidList::info() const
  { return fmt::format( "ClothoidList\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  G2solveCLC::save( ostream_type & stream ) const {
    stream << "# x y theta kappa\n";
    save_segment( stream, m_S0 );
    save_segment( stream, m_SM );
    save_segment( stream, m_S1 );
    stream << "# EOF\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  G2solve3arc::save( ostream_type & stream ) const {
    stream << "# x y theta kappa\n";
    save_segment( stream, m_S0 );
    save_segment( stream, m_SM );
    save_segment( stream, m_S1 );
    stream << "# EOF\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  ClothoidSplineG2::info() const
  { return fmt::format( "ClothoidSplineG2\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}

// EOF: ClothoidList.cc
