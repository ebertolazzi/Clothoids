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

// workaround for windows that defines max and min as macros!
#ifdef max
  #undef max
#endif
#ifdef min
  #undef min
#endif

#include <cmath>
#include <cfloat>
#include <algorithm>

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

  using std::vector;
  using std::abs;
  using std::min;
  using std::max;
  using std::swap;
  using std::ceil;
  using std::floor;
  using std::isfinite;

  integer   ClothoidCurve::m_max_iter  = 10;
  real_type ClothoidCurve::m_tolerance = 0.01745329252; // 1 degree

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::setup( GenericContainer const & gc ) {
    string cwhere{ fmt::format("ClothoidCurve[{}]::setup( gc ):", this->name() ) };
    char const * where{ cwhere.c_str() };
    real_type x0     = gc.get_map_number("x0",     where );
    real_type y0     = gc.get_map_number("y0",     where );
    real_type theta0 = gc.get_map_number("theta0", where );
    real_type x1     = gc.get_map_number("x1",     where );
    real_type y1     = gc.get_map_number("y1",     where );
    real_type theta1 = gc.get_map_number("theta1", where );
    bool ok = this->build_G1( x0, y0, theta0, x1, y1, theta1 );
    UTILS_ASSERT( ok, "ClothoidCurve[{}]::setup( gc ) failed\n", this->name() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve::ClothoidCurve( string const & name )
  : BaseCurve( name )
  {
    m_CD.m_x0     = 0;
    m_CD.m_y0     = 0;
    m_CD.m_theta0 = 0;
    m_CD.m_kappa0 = 0;
    m_CD.m_dk     = 0;
    m_L           = 0;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve( ClothoidCurve const & s )
  : BaseCurve( s.name() )
  { this->copy(s); }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type k,
    real_type dk,
    real_type L,
    string const & name
  ) : BaseCurve( name ) {
    m_CD.m_x0     = x0;
    m_CD.m_y0     = y0;
    m_CD.m_theta0 = theta0;
    m_CD.m_kappa0 = k;
    m_CD.m_dk     = dk;
    m_L           = L;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve(
    real_type const   P0[],
    real_type         theta0,
    real_type const   P1[],
    real_type         theta1,
    string    const & name
  ) : BaseCurve( name ) {
    build_G1( P0[0], P0[1], theta0, P1[0], P1[1], theta1 );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::copy( ClothoidCurve const & c ) {
    m_CD        = c.m_CD;
    m_L         = c.m_L;
    m_aabb_done = false;
    m_aabb_triangles.clear();
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::build( LineSegment const & LS ) {
    m_CD.m_x0     = LS.m_x0;
    m_CD.m_y0     = LS.m_y0;
    m_CD.m_theta0 = LS.m_theta0;
    m_CD.m_kappa0 = 0;
    m_CD.m_dk     = 0;
    m_L           = LS.m_L;
    m_aabb_done   = false;
    m_aabb_triangles.clear();
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::build( CircleArc const & C ) {
    m_CD.m_x0     = C.m_x0;
    m_CD.m_y0     = C.m_y0;
    m_CD.m_theta0 = C.m_theta0;
    m_CD.m_kappa0 = C.m_k;
    m_CD.m_dk     = 0;
    m_L           = C.m_L;
    m_aabb_done   = false;
    m_aabb_triangles.clear();
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::build( ClothoidCurve const & C )
  { this->copy(C); }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void ClothoidCurve::build( Biarc const & )        { UTILS_ERROR("can convert from Biarc to ClothoidCurve\n"); }
  void ClothoidCurve::build( PolyLine const & )     { UTILS_ERROR("can convert from PolyLine to ClothoidCurve\n"); }
  void ClothoidCurve::build( BiarcList const & )    { UTILS_ERROR("can convert from BiarcList to ClothoidCurve\n"); }
  void ClothoidCurve::build( ClothoidList const & ) { UTILS_ERROR("can convert from ClothoidList to ClothoidCurve\n"); }
  void ClothoidCurve::build( Dubins const & )       { UTILS_ERROR("can convert from Dubins to ClothoidCurve\n"); }
  void ClothoidCurve::build( Dubins3p const & )     { UTILS_ERROR("can convert from Dubins3p to ClothoidCurve\n"); }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::build(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type k,
    real_type dk,
    real_type L
  ) {
    UTILS_ASSERT(
      L > 0,
      "ClothoidCurve::build( x0={}, y0={}, theta0={}, k={}, dk={}, L={} )\n"
      "L must be positive!\n",
      x0, y0, theta0, k, dk, L
    );
    m_CD.m_x0     = x0;
    m_CD.m_y0     = y0;
    m_CD.m_theta0 = theta0;
    m_CD.m_kappa0 = k;
    m_CD.m_dk     = dk;
    m_L           = L;
    m_aabb_done   = false;
    m_aabb_triangles.clear();
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  int
  ClothoidCurve::build_G1(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1,
    real_type tol
  ) {
    m_aabb_done = false;
    m_aabb_triangles.clear();
    return m_CD.build_G1( x0, y0, theta0, x1, y1, theta1, tol, m_L );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  int
  ClothoidCurve::build_G1_D(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type x1,
    real_type y1,
    real_type theta1,
    real_type L_D[2],
    real_type k_D[2],
    real_type dk_D[2],
    real_type tol
  ) {
    m_aabb_done = false;
    m_aabb_triangles.clear();
    return m_CD.build_G1( x0, y0, theta0, x1, y1, theta1, tol, m_L, true, L_D, k_D, dk_D );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  bool
  ClothoidCurve::build_forward(
    real_type x0,
    real_type y0,
    real_type theta0,
    real_type kappa0,
    real_type x1,
    real_type y1,
    real_type tol
  ) {
    m_aabb_done = false;
    m_aabb_triangles.clear();
    return m_CD.build_forward( x0, y0, theta0, kappa0, x1, y1, tol, m_L );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve( LineSegment const & LS )
  : BaseCurve( LS.name() ) {
    m_CD.m_x0     = LS.m_x0;
    m_CD.m_y0     = LS.m_y0;
    m_CD.m_theta0 = LS.m_theta0;
    m_CD.m_kappa0 = 0;
    m_CD.m_dk     = 0;
    m_L           = LS.m_L;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve( CircleArc const & C )
  : BaseCurve( C.name() ) {
    m_CD.m_x0     = C.m_x0;
    m_CD.m_y0     = C.m_y0;
    m_CD.m_theta0 = C.m_theta0;
    m_CD.m_kappa0 = C.m_k;
    m_CD.m_dk     = 0;
    m_L           = C.m_L;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


  ClothoidCurve::ClothoidCurve( BaseCurve const * pC )
  : ClothoidCurve( pC->name() ) {

    G2LIB_DEBUG_MESSAGE( "ClothoidCurve convert: {}\n", pC->type_name() );

    switch ( pC->type() ) {
    case CurveType::LINE:
      G2LIB_DEBUG_MESSAGE( "LineSegment -> ClothoidCurve\n" );
      this->build( *static_cast<LineSegment const *>(pC) );
      break;
    case CurveType::CIRCLE:
      G2LIB_DEBUG_MESSAGE( "CircleArc -> ClothoidCurve\n" );
      this->build( *static_cast<CircleArc const *>(pC) );
      break;
    case CurveType::CLOTHOID:
      G2LIB_DEBUG_MESSAGE( "ClothoidCurve -> ClothoidCurve\n" );
      this->copy( *static_cast<ClothoidCurve const *>(pC) );
      break;
    default:
      UTILS_ERROR(
        "ClothoidList constructor cannot convert from: {}\n",
        pC->type_name()
      );
      break;
    }
  }

  real_type
  ClothoidCurve::length_ISO( real_type ) const {
    UTILS_ERROR0( "Offset length not available for Clothoids\n" );
    //return 0;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::optimized_sample_internal_ISO(
    real_type           s_begin,
    real_type           s_end,
    real_type           offs,
    real_type           ds,
    real_type           max_angle,
    vector<real_type> & s
  ) const {
    real_type ss  = s_begin;
    real_type thh = theta(s_begin);
    for ( integer npts = 0; ss < s_end; ++npts ) {
      UTILS_ASSERT0(
        npts < 100000000,
        "ClothoidCurve::optimized_sample_internal "
        "is generating too much points (>100000000)\n"
        "something is going wrong or parameters are not well set\n"
      );
      // estimate angle variation and compute step accodingly
      real_type k   = m_CD.kappa( ss );
      real_type dss = ds/(1+k*offs); // scale length with offset
      real_type sss = ss + dss;
      if ( sss > s_end ) {
        sss = s_end;
        dss = s_end-ss;
      }
      if ( abs(k*dss) > max_angle ) {
        dss = abs(max_angle/k);
        sss = ss + dss;
      }
      // check and recompute if necessary
      real_type thhh = theta(sss);
      if ( abs(thh-thhh) > max_angle ) {
        k    = m_CD.kappa( sss );
        dss  = abs(max_angle/k);
        sss  = ss + dss;
        thhh = theta(sss);
      }
      ss  = sss;
      thh = thhh;
      s.emplace_back(ss);
    }
    s.back() = s_end;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::optimized_sample_ISO(
    real_type           offs,
    integer             npts,
    real_type           max_angle,
    vector<real_type> & s
  ) const {
    s.clear();
    s.reserve( size_t(npts) );
    s.emplace_back(0);

    real_type ds = m_L/npts;
    if ( m_CD.m_kappa0*m_CD.m_dk >= 0 || m_CD.kappa(m_L)*m_CD.m_dk <= 0 ) {
      optimized_sample_internal_ISO( 0, m_L, offs, ds, max_angle, s );
    } else {
      // flex inside, split clothoid
      real_type sflex = -m_CD.m_kappa0/m_CD.m_dk;
      optimized_sample_internal_ISO( 0,   sflex, offs, ds, max_angle, s );
      optimized_sample_internal_ISO( sflex, m_L, offs, ds, max_angle, s );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  /*\
   |  _    _   _____    _                _
   | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
   | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
   | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
   |                               |___/
  \*/

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::bb_triangles_internal_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            s_begin,
    real_type            s_end,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {

    static real_type const one_degree = Utils::m_pi/180;

    real_type ss  = s_begin;
    real_type thh = m_CD.theta(ss);
    real_type MX  = min( m_L, max_size );
    for ( integer npts = 0; ss < s_end; ++npts ) {
      UTILS_ASSERT0(
        npts < 100000000,
        "ClothoidCurve::bb_triangles_internal "
        "is generating too much triangles (>100000000)\n"
        "something is going wrong or parameters are not well set\n"
      );

      // estimate angle variation and compute step accodingly
      real_type k   = m_CD.kappa( ss );
      real_type dss = MX/(1+k*offs); // scale length with offset
      real_type sss = ss + dss;
      if ( sss > s_end ) {
        sss = s_end;
        dss = s_end-ss;
      }
      if ( abs(k*dss) > max_angle ) {
        dss = abs(max_angle/k);
        sss = ss + dss;
      }
      // check and recompute if necessary
      real_type thhh = theta(sss);
      if ( abs(thh-thhh) > max_angle ) {
        k    = m_CD.kappa( sss );
        dss  = abs(max_angle/k);
        sss  = ss + dss;
        thhh = theta(sss);
      }

      real_type x0, y0, x1, y1;
      m_CD.eval_ISO( ss,  offs, x0, y0 );
      m_CD.eval_ISO( sss, offs, x1, y1 );

      real_type tx0   = cos(thh);
      real_type ty0   = sin(thh);
      real_type alpha = sss-ss; // se angolo troppo piccolo uso approx piu rozza
      if ( abs(thh-thhh) > one_degree ) {
        real_type tx1 = cos(thhh);
        real_type ty1 = sin(thhh);
        real_type det = tx1 * ty0 - tx0 * ty1;
        real_type dx  = x1-x0;
        real_type dy  = y1-y0;
        alpha = (dy*tx1 - dx*ty1)/det;
      }

      real_type x2 = x0 + alpha*tx0;
      real_type y2 = y0 + alpha*ty0;
      Triangle2D t( x0, y0, x2, y2, x1, y1, ss, sss, icurve );

      tvec.emplace_back( t );

      ss  = sss;
      thh = thhh;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::bb_triangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    integer              icurve
  ) const {
    if ( m_CD.m_kappa0*m_CD.m_dk >= 0 || m_CD.kappa(m_L)*m_CD.m_dk <= 0 ) {
      bb_triangles_internal_ISO( offs, tvec, 0, m_L, max_angle, max_size, icurve );
    } else {
      // flex inside, split clothoid
      real_type sflex = -m_CD.m_kappa0/m_CD.m_dk;
      bb_triangles_internal_ISO( offs, tvec, 0,   sflex, max_angle, max_size, icurve );
      bb_triangles_internal_ISO( offs, tvec, sflex, m_L, max_angle, max_size, icurve );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  /*\
   |  ___ ___
   | | _ ) _ ) _____ __
   | | _ \ _ \/ _ \ \ /
   | |___/___/\___/_\_\
  \*/

  void
  ClothoidCurve::bbox_ISO(
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

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  /*\
   |     _        _    ____  ____  _
   |    / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |   / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |  / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   | /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  void
  ClothoidCurve::build_AABBtree_ISO(
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
    for ( auto const & clot : m_aabb_triangles ) {
      clot.bbox( bbox_min[0], bbox_min[1], bbox_max[0], bbox_max[1] );
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |            _ _ _     _
   |   ___ ___ | | (_)___(_) ___  _ __
   |  / __/ _ \| | | / __| |/ _ \| '_ \
   | | (_| (_) | | | \__ \ | (_) | | | |
   |  \___\___/|_|_|_|___/_|\___/|_| |_|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::collision_ISO(
    real_type             offs,
    ClothoidCurve const & C,
    real_type             offs_C
  ) const {

    G2LIB_DEBUG_MESSAGE(
      "ClothoidCurve::collision_ISO( offs={}, C, offs_C={}\n",
      offs, offs_C
    );

    G2LIB_DEBUG_TICTOC;

    G2LIB_DEBUG_TIC;
    this->build_AABBtree_ISO( offs );
    C.build_AABBtree_ISO( offs_C );
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidCurve::collision_ISO( offs={}, C, offs_C={} ) build_AABBtree_ISO elapsed {}ms\n",
      offs, offs_C, tictoc.elapsed_ms()
    );

    G2LIB_DEBUG_TIC;
    AABB_MAP intersectList;
    m_aabb_tree.intersect_and_refine( C.m_aabb_tree, intersectList );
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidCurve::collision_ISO intersect_and_refine elapsed {}ms candidated #{}\n",
      tictoc.elapsed_ms(), intersectList.size()
    );

    G2LIB_DEBUG_TIC;
    bool collide = false;
    for ( auto const & I : intersectList ) {
      integer i = I.first;
      UTILS_ASSERT_DEBUG(
        i >= 0 && i < integer(m_aabb_triangles.size()),
        "ClothoidCurve::collision_ISO( offs={}, C, offs_C={} ) i={} out of range [0,{})\n",
        offs, offs_C, i, m_aabb_triangles.size()
      );
      Triangle2D const & T1 = m_aabb_triangles.at(i);
      for ( integer j : I.second ) {
        UTILS_ASSERT_DEBUG(
          j >= 0 && j < integer(C.m_aabb_triangles.size()),
          "ClothoidCurve::collision_ISO( offs={}, C, offs_C={} ) j={} out of range [0,{})\n",
          offs, offs_C, j, C.m_aabb_triangles.size()
        );
        Triangle2D const & T2 = C.m_aabb_triangles.at(j);
        real_type ss1, ss2;
        collide = this->aabb_intersect_ISO( T1, offs, &C, T2, offs_C, ss1, ss2 );
        if ( collide ) break;
      }
      if ( collide ) break;
    }

    G2LIB_DEBUG_TOC;
    G2LIB_DEBUG_MESSAGE(
      "ClothoidCurve collision_ISO: collisions elapsed {}ms\n",
      tictoc.elapsed_ms()
    );

    return collide;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Collision detection
  //!
  //! \param[in] offs      curve offset
  //! \param[in] C         curve to compare for collision detection
  //! \param[in] offs_C    curve offset
  //! \param[in] max_angle maximum angle variation
  //! \param[in] max_size  if the segment is larger then this parameter is split
  //!
  bool
  ClothoidCurve::approximate_collision_ISO(
    real_type             offs,
    ClothoidCurve const & C,
    real_type             offs_C,
    real_type             max_angle,
    real_type             max_size
  ) const {

    this->build_AABBtree_ISO( offs, max_angle, max_size );
    C.build_AABBtree_ISO( offs_C, max_angle, max_size );

    AABB_MAP intersectList;
    m_aabb_tree.intersect_and_refine( C.m_aabb_tree, intersectList );

    for ( auto const & I : intersectList ) {
      integer i = I.first;
      UTILS_ASSERT_DEBUG(
        i >= 0 && i < integer(m_aabb_triangles.size()),
        "ClothoidCurve::collision_ISO( offs={}, C, offs_C={} ) i={} out of range [0,{})\n",
        offs, offs_C, i, m_aabb_triangles.size()
      );
      Triangle2D const & T1 = m_aabb_triangles.at(i);
      for ( integer j : I.second ) {
        UTILS_ASSERT_DEBUG(
          j >= 0 && j < integer(C.m_aabb_triangles.size()),
          "ClothoidCurve::collision_ISO( offs={}, C, offs_C={} ) j={} out of range [0,{})\n",
          offs, offs_C, j, C.m_aabb_triangles.size()
        );
        Triangle2D const & T2 = C.m_aabb_triangles.at(j);
        bool collide = T1.overlap(T2);
        if ( collide ) return true;
      }
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::collision( BaseCurve const * pC ) const {
    if ( pC->type() == CurveType::CLOTHOID ) {
      ClothoidCurve const & C = *static_cast<ClothoidCurve const *>(pC);
      return this->collision( C );
    } else {
      ClothoidCurve C(pC);
      return this->collision( C );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::collision_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C
  ) const {
    if ( pC->type() == CurveType::CLOTHOID ) {
      ClothoidCurve const & C = *static_cast<ClothoidCurve const *>(pC);
      return this->collision_ISO( offs, C, offs_C );
    } else {
      CurveType CT = curve_promote( this->type(), pC->type() );
      if ( CT == CurveType::CLOTHOID ) {
        ClothoidCurve C(pC);
        return this->collision_ISO( offs, C, offs_C );
      } else {
        return G2lib::collision_ISO( this, offs, pC, offs_C );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |  _       _                          _
   | (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   | | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   | | | | | | ||  __/ |  \__ \  __/ (__| |_
   | |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::aabb_intersect_ISO(
    Triangle2D    const & T1,
    real_type             offs,
    ClothoidCurve const * pC,
    Triangle2D    const & T2,
    real_type             offs_C,
    real_type           & ss1,
    real_type           & ss2
  ) const {
    real_type eps1   = machepsi1000*m_L;
    real_type eps2   = machepsi1000*pC->m_L;
    real_type s1_min = T1.S0()-eps1;
    real_type s1_max = T1.S1()+eps1;
    real_type s2_min = T2.S0()-eps2;
    real_type s2_max = T2.S1()+eps2;
    integer   nout   = 0;
    bool converged   = false;

    ss1 = (s1_min+s1_max)/2;
    ss2 = (s2_min+s2_max)/2;
    for ( integer i = 0; i < m_max_iter && !converged; ++i ) {
      real_type t1[2], t2[2], p1[2], p2[2];
      m_CD.eval_ISO  ( ss1, offs, p1[0], p1[1] );
      m_CD.eval_ISO_D( ss1, offs, t1[0], t1[1] );
      pC->m_CD.eval_ISO  ( ss2, offs_C, p2[0], p2[1] );
      pC->m_CD.eval_ISO_D( ss2, offs_C, t2[0], t2[1] );
      /*
      // risolvo il sistema
      // p1 + alpha * t1 = p2 + beta * t2
      // alpha * t1 - beta * t2 = p2 - p1
      //
      //  / t1[0] -t2[0] \ / alpha \ = / p2[0] - p1[0] \
      //  \ t1[1] -t2[1] / \ beta  /   \ p2[1] - p1[1] /
      */
      real_type det = t2[0]*t1[1]-t1[0]*t2[1];
      real_type px  = p2[0]-p1[0];
      real_type py  = p2[1]-p1[1];
      ss1 += (py*t2[0] - px*t2[1])/det;
      ss2 += (t1[0]*py - t1[1]*px)/det;
      if ( ! ( isfinite(ss1) && isfinite(ss1) ) ) break;
      bool out = false;
      if      ( ss1 < s1_min ) { out = true; ss1 = s1_min; }
      else if ( ss1 > s1_max ) { out = true; ss1 = s1_max; }
      if      ( ss2 < s2_min ) { out = true; ss2 = s2_min; }
      else if ( ss2 > s2_max ) { out = true; ss2 = s2_max; }
      if ( out ) {
        if ( ++nout > 3 ) break;
      } else {
        converged = abs(px) <= m_tolerance && abs(py) <= m_tolerance;
      }
    }
    if ( converged ) {
      if      ( ss1 < T1.S0() ) ss1 = T1.S0();
      else if ( ss1 > T1.S1() ) ss1 = T1.S1();
      if      ( ss2 < T2.S0() ) ss2 = T2.S0();
      else if ( ss2 > T2.S1() ) ss2 = T2.S1();
    }
    return converged;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::intersect_ISO(
    real_type             offs,
    ClothoidCurve const & C,
    real_type             offs_C,
    IntersectList       & ilist
  ) const {
    real_type ss1, ss2;
    G2LIB_DEBUG_TICTOC;

    if ( intersect_with_AABBtree ) {
      G2LIB_DEBUG_TIC;

      this->build_AABBtree_ISO( offs );
      C.build_AABBtree_ISO( offs_C );

      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve intersect_ISO: build_AABBtree_ISO elapsed {}ms\n",
        tictoc.elapsed_ms()
      );

      G2LIB_DEBUG_TIC;
      AABB_MAP intersectList;
      m_aabb_tree.intersect_and_refine( C.m_aabb_tree, intersectList );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve intersect_ISO: intersect_and_refine elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), intersectList.size()
      );

      G2LIB_DEBUG_TIC;
      for ( auto const & I : intersectList ) {
        integer i = I.first;
        UTILS_ASSERT_DEBUG(
          i >= 0 && i < integer(m_aabb_triangles.size()),
          "ClothoidCurve::intersect_ISO( offs={}, C, offs_C={}, ilist ) i={} out of range [0,{})\n",
          offs, offs_C, i, m_aabb_triangles.size()
        );
        Triangle2D const & T1 = m_aabb_triangles.at(i);
        for ( integer j : I.second ) {
          UTILS_ASSERT_DEBUG(
            j >= 0 && j < integer(C.m_aabb_triangles.size()),
            "ClothoidCurve::intersect_ISO( offs={}, C, offs_C={}, ilist ) j={} out of range [0,{})\n",
            offs, offs_C, j, C.m_aabb_triangles.size()
          );
          Triangle2D const & T2 = C.m_aabb_triangles.at(j);
          bool converged = aabb_intersect_ISO( T1, offs, &C, T2, offs_C, ss1, ss2 );
          if ( converged ) ilist.emplace_back( ss1, ss2 );
        }
      }

      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve intersect_ISO: intersections elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), intersectList.size()
      );

    } else {

      G2LIB_DEBUG_TIC;
      bb_triangles_ISO( offs, m_aabb_triangles, Utils::m_pi/18, 1e100 );
      C.bb_triangles_ISO( offs_C, C.m_aabb_triangles, Utils::m_pi/18, 1e100 );

      for ( Triangle2D const & T1 : m_aabb_triangles ) {
        for ( Triangle2D const & T2 : C.m_aabb_triangles ) {
          bool converged = aabb_intersect_ISO( T1, offs, &C, T2, offs_C, ss1, ss2 );
          if ( converged ) ilist.emplace_back( ss1, ss2 );
        }
      }

      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve intersect_ISO: intersections elapsed {}ms noAABB tree\n",
        tictoc.elapsed_ms()
      );
    }
  }

  void
  ClothoidCurve::intersect(
    BaseCurve const * pC,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CLOTHOID ) {
      ClothoidCurve const & C = *static_cast<ClothoidCurve const *>(pC);
      this->intersect( C, ilist );
    } else {
      ClothoidCurve C(pC);
      this->intersect( C, ilist );
    }
  }

  void
  ClothoidCurve::intersect_ISO(
    real_type         offs,
    BaseCurve const * pC,
    real_type         offs_C,
    IntersectList   & ilist
  ) const {
    if ( pC->type() == CurveType::CLOTHOID ) {
      ClothoidCurve const & C = *static_cast<ClothoidCurve const *>(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    } else {
      ClothoidCurve C(pC);
      this->intersect_ISO( offs, C, offs_C, ilist );
    }
  }

  /*\
   |        _                     _   ____       _       _
   |    ___| | ___  ___  ___  ___| |_|  _ \ ___ (_)_ __ | |_
   |   / __| |/ _ \/ __|/ _ \/ __| __| |_) / _ \| | '_ \| __|
   |  | (__| | (_) \__ \  __/\__ \ |_|  __/ (_) | | | | | |_
   |   \___|_|\___/|___/\___||___/\__|_|   \___/|_|_| |_|\__|
  \*/

  void
  ClothoidCurve::closest_point_internal(
    real_type   s_begin,
    real_type   s_end,
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & dst
  ) const {
    #if 1
    // minimize using circle approximation
    s = (s_begin + s_end)/2;
    integer nout{0};
    integer n_ok{0};
    for ( integer iter = 0; iter < m_max_iter; ++iter ) {
      // osculating circle
      m_CD.eval_ISO( s, offs, x, y );
      real_type th = m_CD.theta( s );
      real_type kk = m_CD.kappa( s );
      real_type sc = 1+kk*offs;
      real_type ds = projectPointOnCircle( x, y, th, kk/sc, qx, qy )/sc;

      s += ds;

      bool out = false;
      if      ( s <= s_begin ) { out = true; s = s_begin; }
      else if ( s >= s_end   ) { out = true; s = s_end; }

      if ( out ) {
        if ( ++nout > 3 ) break;
      } else {
        // force one more itaration to improve accuracy
        if ( abs(ds) <= m_tolerance && ++n_ok >= 2 ) break;
      }
    }
    dst = hypot( qx-x, qy-y );
    #else
    real_type ds = (s_end-s_begin)/10;
    for ( integer iter = 0; iter <= 10; ++iter ) {
      real_type ss = s_begin + iter * ds;
      real_type xx, yy;
      CD.eval( ss, offs, xx, yy );
      real_type dx = xx-qx;
      real_type dy = yy-qy;
      real_type dst1 = hypot( dx, dy );
      if ( dst1 < dst ) {
        s   = ss;
        x   = xx;
        y   = yy;
        dst = dst1;
      }
    }
    #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidCurve::closest_point_internal(
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
    real_type xx, yy, ss;
    this->build_AABBtree_ISO( offs );
    DST = Utils::Inf<real_type>();
    G2LIB_DEBUG_TOC;

    G2LIB_DEBUG_MESSAGE(
      "ClothoidCurve::closest_point_internal: build_AABBtree_ISO elapsed {}ms\n",
      tictoc.elapsed_ms()
    );

    if ( m_aabb_tree.num_tree_nodes() > G2LIB_AABB_MIN_NODES && intersect_with_AABBtree ) {

      G2LIB_DEBUG_TIC;
      AABB_SET candidateList;
      real_type xy[2] = { qx, qy };
      m_aabb_tree.min_distance_candidates( xy, candidateList );
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve::closest_point_internal min_distance_candidates elapsed {}ms candidated #{}\n",
        tictoc.elapsed_ms(), candidateList.size()
      );

      UTILS_ASSERT0(
         candidateList.size() > 0,
        "ClothoidCurve::closest_point_internal no candidate\n"
      );

      G2LIB_DEBUG_TIC;
      for ( integer ipos : candidateList ) {
        UTILS_ASSERT_DEBUG(
          ipos >= 0 && ipos < integer(m_aabb_triangles.size()),
          "ClothoidCurve::closest_point_internal( qx={}, qy={}, offs={}, x, y, s, DST ) ipos={} out of range [0,{})\n",
          qx, qy, offs, ipos, m_aabb_triangles.size()
        );
        Triangle2D const & T = m_aabb_triangles.at(ipos);
        real_type dst = T.dist_min( qx, qy );
        if ( dst < DST ) {
          // refine distance
          closest_point_internal( T.S0(), T.S1(), qx, qy, offs, xx, yy, ss, dst );
          if ( dst < DST ) {
            DST = dst;
            s   = ss;
            x   = xx;
            y   = yy;
          }
        }
      }
      G2LIB_DEBUG_TOC;

      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve::closest_point_internal minimize elapsed {}ms\n",
        tictoc.elapsed_ms()
      );

    } else {

      G2LIB_DEBUG_TIC;
      for ( Triangle2D const & T : m_aabb_triangles ) {
        real_type dst = T.dist_min( qx, qy );
        if ( dst < DST ) {
          // refine distance
          closest_point_internal( T.S0(), T.S1(), qx, qy, offs, xx, yy, ss, dst );
          if ( dst < DST ) {
            DST = dst;
            s   = ss;
            x   = xx;
            y   = yy;
          }
        }
      }
      G2LIB_DEBUG_TOC;
      G2LIB_DEBUG_MESSAGE(
        "ClothoidCurve::closest_point_internal minimize elapsed {}ms (no AABB tree)\n",
        tictoc.elapsed_ms()
      );

    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  integer
  ClothoidCurve::closest_point_ISO(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & DST
  ) const {

    this->closest_point_internal( qx, qy, offs, x, y, s, DST );

    // check if projection is orthogonal
    real_type nx, ny;
    nor_ISO( s, nx, ny );
    real_type qxx = qx - x;
    real_type qyy = qy - y;
    t = qxx * nx + qyy * ny - offs; // signed distance
    real_type pt = abs(qxx * ny - qyy * nx);
    G2LIB_DEBUG_MESSAGE(
      "Clothoid::closest_point_ISO\n"
      "||P-P0|| = {} and {}, |(P-P0).T| = {}\n",
      DST, hypot(qxx,qyy), pt
    );
    return pt > GLIB2_TOL_ANGLE*hypot(qxx,qyy) ? -1 : 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::theta_total_variation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = m_CD.m_kappa0;
    real_type kR  = m_CD.kappa(m_L);
    real_type thL = 0;
    real_type thR = m_CD.delta_theta(m_L);
    if ( kL*kR < 0 ) {
      real_type root = -m_CD.m_kappa0/m_CD.m_dk;
      if ( root > 0 && root < m_L ) {
        real_type thM  = m_CD.delta_theta(root);
        return abs( thR - thM ) + abs( thM - thL );
      }
    }
    return abs( thR - thL );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::theta_min_max( real_type & thMin, real_type & thMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = m_CD.m_kappa0;
    real_type kR  = m_CD.kappa(m_L);
    real_type thL = 0;
    real_type thR = m_CD.delta_theta(m_L);
    if ( thL < thR ) { thMin = thL; thMax = thR; }
    else             { thMin = thR; thMax = thL; }
    if ( kL*kR < 0 ) {
      real_type root = -m_CD.m_kappa0/m_CD.m_dk;
      if ( root > 0 && root < m_L ) {
        real_type thM = m_CD.delta_theta(root);
        if      ( thM < thMin ) thMin = thM;
        else if ( thM > thMax ) thMax = thM;
      }
    }
    return thMax - thMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvature_min_max( real_type & kMin, real_type & kMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    kMin = m_CD.m_kappa0;
    kMax = m_CD.kappa(m_L);
    if ( kMax < kMin ) swap( kMax, kMin );
    return kMax - kMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvature_total_variation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type km = m_CD.m_kappa0;
    real_type kp = m_CD.kappa(m_L);
    return abs(kp-km);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integral_curvature2() const {
    return m_L*( m_CD.m_kappa0*(m_CD.m_kappa0+m_L*m_CD.m_dk) +
                 (m_L*m_L)*m_CD.m_dk*m_CD.m_dk/3 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integral_jerk2() const {
    real_type k2 = m_CD.m_kappa0*m_CD.m_kappa0;
    real_type k3 = m_CD.m_kappa0*k2;
    real_type k4 = k2*k2;
    real_type t1 = m_L;
    real_type t2 = m_L*t1;
    real_type t3 = m_L*t2;
    real_type t4 = m_L*t3;
    return ((((t4/5*m_CD.m_dk+t3*m_CD.m_kappa0)*m_CD.m_dk+(1+2*t2)*k2)*m_CD.m_dk+2*t1*k3)*m_CD.m_dk+k4)*m_L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integral_snap2() const {
    real_type k2  = m_CD.m_kappa0*m_CD.m_kappa0;
    real_type k3  = m_CD.m_kappa0*k2;
    real_type k4  = k2*k2;
    real_type k5  = k4*m_CD.m_kappa0;
    real_type k6  = k4*k2;
    real_type dk2 = m_CD.m_dk*m_CD.m_dk;
    real_type dk3 = m_CD.m_dk*dk2;
    real_type dk4 = dk2*dk2;
    real_type dk5 = dk4*m_CD.m_dk;
    real_type dk6 = dk4*dk2;
    real_type t2  = m_L;
    real_type t3  = m_L*t2;
    real_type t4  = m_L*t3;
    real_type t5  = m_L*t4;
    real_type t6  = m_L*t5;
    real_type t7  = m_L*t6;

    return ( (t7/7)*dk6 +
             dk5*m_CD.m_kappa0*t6 +
             3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 +
             3*m_CD.m_dk*k5*t2 +
             9*dk2*m_CD.m_kappa0*t2 +
             k6+9*k2*m_CD.m_dk ) * m_L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  string
  ClothoidCurve::info() const
  { return fmt::format( "Clothoid\n{}\n", *this ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!
  //!  Print on strem the `ClothoidCurve` object
  //!
  //!  \param stream the output stream
  //!  \param c      an instance of `ClothoidCurve` object
  //!  \return the output stream
  //!
  ostream_type &
  operator << ( ostream_type & stream, ClothoidCurve const & c ) {
    fmt::print( stream,
      "x0     = {:<12} x1     = {:<12}\n"
      "y0     = {:<12} y1     = {:<12}\n"
      "theta0 = {:<12} theta1 = {:<12}\n"
      "kappa0 = {:<12} kappa1 = {:<12}\n"
      "dk     = {:<12} L      = {:<12}\n",
      fmt::format("{:.6}",c.x_begin()),
      fmt::format("{:.6}",c.x_end()),
      fmt::format("{:.6}",c.y_begin()),
      fmt::format("{:.6}",c.y_end()),
      fmt::format("{:.6}",c.theta_begin()),
      fmt::format("{:.6}",c.theta_end()),
      fmt::format("{:.6}",c.kappa_begin()),
      fmt::format("{:.6}",c.kappa_end()),
      fmt::format("{:.6}",c.m_CD.m_dk),
      fmt::format("{:.6}",c.m_L)
    );
    return stream;
  }

}

// EOF: Clothoid.cc
