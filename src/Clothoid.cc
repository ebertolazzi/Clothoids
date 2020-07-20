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

#include "Line.hh"
#include "Circle.hh"
#include "Biarc.hh"
#include "Clothoid.hh"
#include "ClothoidList.hh"
#include "PolyLine.hh"

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

  using std::vector;
  using std::abs;
  using std::min;
  using std::max;
  using std::swap;
  using std::ceil;
  using std::floor;
  using std::isfinite;
  using std::numeric_limits;

  int_type  ClothoidCurve::max_iter  = 10;
  real_type ClothoidCurve::tolerance = 1e-9;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ClothoidCurve::ClothoidCurve( BaseCurve const & C )
  : BaseCurve(G2LIB_CLOTHOID)
  , aabb_done(false)
  {
    switch ( C.type() ) {
    case G2LIB_LINE:
      build( *static_cast<LineSegment const *>(&C) );
      break;
    case G2LIB_CIRCLE:
      build( *static_cast<CircleArc const *>(&C) );
      break;
    case G2LIB_CLOTHOID:
      copy( *static_cast<ClothoidCurve const *>(&C) );
      break;
    case G2LIB_BIARC:
    case G2LIB_BIARC_LIST:
    case G2LIB_CLOTHOID_LIST:
    case G2LIB_POLYLINE:
      G2LIB_DO_ERROR(
        "ClothoidList constructor cannot convert from: " <<
        CurveType_name[C.type()]
      )
    }
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
    for ( int_type npts = 0; ss < s_end; ++npts ) {
      G2LIB_ASSERT(
        npts < 100000000,
        "ClothoidCurve::optimized_sample_internal " <<
        "is generating too much points (>100000000)\n" <<
        "something is going wrong or parameters are not well set"
      )
      // estimate angle variation and compute step accodingly
      real_type k   = CD.kappa( ss );
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
        k    = CD.kappa( sss );
        dss  = abs(max_angle/k);
        sss  = ss + dss;
        thhh = theta(sss);
      }
      ss  = sss;
      thh = thhh;
      s.push_back(ss);
    }
    s.back() = s_end;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::optimized_sample_ISO(
    real_type           offs,
    int_type            npts,
    real_type           max_angle,
    vector<real_type> & s
  ) const {
    s.clear();
    s.reserve( size_t(npts) );
    s.push_back(0);

    real_type ds = L/npts;
    if ( CD.kappa0*CD.dk >= 0 || CD.kappa(L)*CD.dk <= 0 ) {
      optimized_sample_internal_ISO( 0, L, offs, ds, max_angle, s );
    } else {
      // flex inside, split clothoid
      real_type sflex = -CD.kappa0/CD.dk;
      optimized_sample_internal_ISO( 0, sflex, offs, ds, max_angle, s );
      optimized_sample_internal_ISO( sflex, L, offs, ds, max_angle, s );
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
  ClothoidCurve::bbTriangles_internal_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            s_begin,
    real_type            s_end,
    real_type            max_angle,
    real_type            max_size,
    int_type             icurve
  ) const {

    static real_type const one_degree = m_pi/180;

    real_type ss  = s_begin;
    real_type thh = CD.theta(ss);
    real_type MX  = min( L, max_size );
    for ( int_type npts = 0; ss < s_end; ++npts ) {
      G2LIB_ASSERT(
        npts < 100000000,
        "ClothoidCurve::bbTriangles_internal " <<
        "is generating too much triangles (>100000000)\n" <<
        "something is going wrong or parameters are not well set"
      )

      // estimate angle variation and compute step accodingly
      real_type k   = CD.kappa( ss );
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
        k    = CD.kappa( sss );
        dss  = abs(max_angle/k);
        sss  = ss + dss;
        thhh = theta(sss);
      }

      real_type x0, y0, x1, y1;
      CD.eval_ISO( ss,  offs, x0, y0 );
      CD.eval_ISO( sss, offs, x1, y1 );

      real_type tx0    = cos(thh);
      real_type ty0    = sin(thh);
      real_type alpha  = sss-ss; // se angolo troppo piccolo uso approx piu rozza
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

      tvec.push_back( t );

      ss  = sss;
      thh = thhh;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  ClothoidCurve::bbTriangles_ISO(
    real_type            offs,
    vector<Triangle2D> & tvec,
    real_type            max_angle,
    real_type            max_size,
    int_type             icurve
  ) const {
    if ( CD.kappa0*CD.dk >= 0 || CD.kappa(L)*CD.dk <= 0 ) {
      bbTriangles_internal_ISO( offs, tvec, 0, L, max_angle, max_size, icurve );
    } else {
      // flex inside, split clothoid
      real_type sflex = -CD.kappa0/CD.dk;
      bbTriangles_internal_ISO( offs, tvec, 0, sflex, max_angle, max_size, icurve );
      bbTriangles_internal_ISO( offs, tvec, sflex, L, max_angle, max_size, icurve );
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
    bbTriangles_ISO( offs, tvec, m_pi/18, 1e100 );
    xmin = ymin = numeric_limits<real_type>::infinity();
    xmax = ymax = -xmin;
    vector<Triangle2D>::const_iterator it;
    for ( it = tvec.begin(); it != tvec.end(); ++it ) {
      // - - - - - - - - - - - - - - - - - - - -
      if      ( it->x1() < xmin ) xmin = it->x1();
      else if ( it->x1() > xmax ) xmax = it->x1();
      if      ( it->x2() < xmin ) xmin = it->x2();
      else if ( it->x2() > xmax ) xmax = it->x2();
      if      ( it->x3() < xmin ) xmin = it->x3();
      else if ( it->x3() > xmax ) xmax = it->x3();
      // - - - - - - - - - - - - - - - - - - - -
      if      ( it->y1() < ymin ) ymin = it->y1();
      else if ( it->y1() > ymax ) ymax = it->y1();
      if      ( it->y2() < ymin ) ymin = it->y2();
      else if ( it->y2() > ymax ) ymax = it->y2();
      if      ( it->y3() < ymin ) ymin = it->y3();
      else if ( it->y3() > ymax ) ymax = it->y3();
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
  void
  ClothoidCurve::build_AABBtree_ISO(
    real_type offs,
    real_type max_angle,
    real_type max_size
  ) const {

    if ( aabb_done &&
         isZero( offs-aabb_offs ) &&
         isZero( max_angle-aabb_max_angle ) &&
         isZero( max_size-aabb_max_size ) ) return;

    #ifdef G2LIB_USE_CXX11
    vector<shared_ptr<BBox const> > bboxes;
    #else
    vector<BBox const *> bboxes;
    #endif

    bbTriangles_ISO( offs, aabb_tri, max_angle, max_size );
    bboxes.reserve(aabb_tri.size());
    vector<Triangle2D>::const_iterator it;
    int_type ipos = 0;
    for ( it = aabb_tri.begin(); it != aabb_tri.end(); ++it, ++ipos ) {
      real_type xmin, ymin, xmax, ymax;
      it->bbox( xmin, ymin, xmax, ymax );
      #ifdef G2LIB_USE_CXX11
      bboxes.push_back( make_shared<BBox const>(
        xmin, ymin, xmax, ymax, G2LIB_CLOTHOID, ipos
      ) );
      #else
      bboxes.push_back(
        new BBox( xmin, ymin, xmax, ymax, G2LIB_CLOTHOID, ipos )
      );
      #endif
    }
    aabb_tree.build(bboxes);
    aabb_done      = true;
    aabb_offs      = offs;
    aabb_max_angle = max_angle;
    aabb_max_size  = max_size;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |            _ _ _     _
   |   ___ ___ | | (_)___(_) ___  _ __
   |  / __/ _ \| | | / __| |/ _ \| '_ \
   | | (_| (_) | | | \__ \ | (_) | | | |
   |  \___\___/|_|_|_|___/_|\___/|_| |_|
  \*/

  bool
  ClothoidCurve::collision( ClothoidCurve const & C ) const {
    this->build_AABBtree_ISO( 0 );
    C.build_AABBtree_ISO( 0 );
    T2D_collision_ISO fun( this, 0, &C, 0 );
    return aabb_tree.collision( C.aabb_tree, fun, false );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidCurve::collision_ISO(
    real_type             offs,
    ClothoidCurve const & C,
    real_type             offs_C
  ) const {
    this->build_AABBtree_ISO( offs );
    C.build_AABBtree_ISO( offs_C );
    T2D_collision_ISO fun( this, offs, &C, offs_C );
    return aabb_tree.collision( C.aabb_tree, fun, false );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // collision detection
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
    T2D_approximate_collision fun( this, &C );
    return aabb_tree.collision( C.aabb_tree, fun, false );
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
    real_type eps1   = machepsi1000*L;
    real_type eps2   = machepsi1000*pC->L;
    real_type s1_min = T1.S0()-eps1;
    real_type s1_max = T1.S1()+eps1;
    real_type s2_min = T2.S0()-eps2;
    real_type s2_max = T2.S1()+eps2;
    int_type  nout   = 0;
    bool converged   = false;

    ss1 = (s1_min+s1_max)/2;
    ss2 = (s2_min+s2_max)/2;
    for ( int_type i = 0; i < max_iter && !converged; ++i ) {
      real_type t1[2], t2[2], p1[2], p2[2];
      CD.eval_ISO  ( ss1, offs, p1[0], p1[1] );
      CD.eval_ISO_D( ss1, offs, t1[0], t1[1] );
      pC->CD.eval_ISO  ( ss2, offs_C, p2[0], p2[1] );
      pC->CD.eval_ISO_D( ss2, offs_C, t2[0], t2[1] );
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
        converged = abs(px) <= tolerance && abs(py) <= tolerance;
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
    IntersectList       & ilist,
    bool                  swap_s_vals
  ) const {
    if ( intersect_with_AABBtree ) {
      this->build_AABBtree_ISO( offs );
      C.build_AABBtree_ISO( offs_C );
      AABBtree::VecPairPtrBBox iList;
      aabb_tree.intersect( C.aabb_tree, iList );
      AABBtree::VecPairPtrBBox::const_iterator ip;

      for ( ip = iList.begin(); ip != iList.end(); ++ip ) {
        size_t ipos1 = size_t(ip->first->Ipos());
        size_t ipos2 = size_t(ip->second->Ipos());

        Triangle2D const & T1 = aabb_tri[ipos1];
        Triangle2D const & T2 = C.aabb_tri[ipos2];

        real_type ss1, ss2;
        bool converged = aabb_intersect_ISO( T1, offs, &C, T2, offs_C, ss1, ss2 );

        if ( converged ) {
          if ( swap_s_vals ) swap( ss1, ss2 );
          ilist.push_back( Ipair( ss1, ss2 ) );
        }
      }
    } else {
      bbTriangles_ISO( offs, aabb_tri, m_pi/18, 1e100 );
      C.bbTriangles_ISO( offs_C, C.aabb_tri, m_pi/18, 1e100 );
      for ( vector<Triangle2D>::const_iterator i1 = aabb_tri.begin();
            i1 != aabb_tri.end(); ++i1 ) {
        for ( vector<Triangle2D>::const_iterator i2 = C.aabb_tri.begin();
              i2 != C.aabb_tri.end(); ++i2 ) {
          Triangle2D const & T1 = *i1;
          Triangle2D const & T2 = *i2;

          real_type ss1, ss2;
          bool converged = aabb_intersect_ISO( T1, offs, &C, T2, offs_C, ss1, ss2 );

          if ( converged ) {
            if ( swap_s_vals ) swap( ss1, ss2 );
            ilist.push_back( Ipair( ss1, ss2 ) );
          }
        }
      }
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
  ClothoidCurve::closestPoint_internal_ISO(
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
    int_type nout = 0;
    for ( int_type iter = 0; iter < max_iter; ++iter ) {
      // osculating circle
      CD.eval_ISO( s, offs, x, y );
      real_type th = CD.theta( s );
      real_type kk = CD.kappa( s );
      real_type sc = 1+kk*offs;
      real_type ds = projectPointOnCircle( x, y, th, kk/sc, qx, qy )/sc;

      s += ds;

      bool out = false;
      if      ( s <= s_begin ) { out = true; s = s_begin; }
      else if ( s >= s_end   ) { out = true; s = s_end; }

      if ( out ) {
        if ( ++nout > 3 ) break;
      } else {
        if ( abs(ds) <= tolerance ) break;
      }
    }
    dst = hypot( qx-x, qy-y );
    #else
    real_type ds = (s_end-s_begin)/10;
    for ( int_type iter = 0; iter <= 10; ++iter ) {
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

  int_type
  ClothoidCurve::closestPoint_ISO(
    real_type   qx,
    real_type   qy,
    real_type   offs,
    real_type & x,
    real_type & y,
    real_type & s,
    real_type & t,
    real_type & DST
  ) const {
    DST = numeric_limits<real_type>::infinity();
    this->build_AABBtree_ISO( offs );

    AABBtree::VecPtrBBox candidateList;
    aabb_tree.min_distance( qx, qy, candidateList );
    AABBtree::VecPtrBBox::const_iterator ic;
    G2LIB_ASSERT(
      candidateList.size() > 0,
      "ClothoidCurve::closestPoint no candidate"
    )
    for ( ic = candidateList.begin(); ic != candidateList.end(); ++ic ) {
      size_t ipos = size_t((*ic)->Ipos());
      Triangle2D const & T = aabb_tri[ipos];
      real_type dst = T.distMin( qx, qy );
      if ( dst < DST ) {
        // refine distance
        real_type xx, yy, ss;
        closestPoint_internal_ISO(
          T.S0(), T.S1(), qx, qy, offs, xx, yy, ss, dst
        );
        if ( dst < DST ) {
          DST = dst;
          s   = ss;
          x   = xx;
          y   = yy;
        }
      }
    }
    real_type nx, ny;
    nor_ISO( s, nx, ny );
    t = (qx-x) * nx + (qy-y) * ny - offs;
    real_type err = abs( abs(t) - DST );
    if ( err > DST*machepsi1000 ) return -1;
    return 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::thetaTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = CD.kappa0;
    real_type kR  = CD.kappa(L);
    real_type thL = 0;
    real_type thR = CD.deltaTheta(L);
    if ( kL*kR < 0 ) {
      real_type root = -CD.kappa0/CD.dk;
      if ( root > 0 && root < L ) {
        real_type thM  = CD.deltaTheta(root);
        return abs( thR - thM ) + abs( thM - thL );
      }
    }
    return abs( thR - thL );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::thetaMinMax( real_type & thMin, real_type & thMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type kL  = CD.kappa0;
    real_type kR  = CD.kappa(L);
    real_type thL = 0;
    real_type thR = CD.deltaTheta(L);
    if ( thL < thR ) { thMin = thL; thMax = thR; }
    else             { thMin = thR; thMax = thL; }
    if ( kL*kR < 0 ) {
      real_type root = -CD.kappa0/CD.dk;
      if ( root > 0 && root < L ) {
        real_type thM = CD.deltaTheta(root);
        if      ( thM < thMin ) thMin = thM;
        else if ( thM > thMax ) thMax = thM;
      }
    }
    return thMax - thMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvatureMinMax( real_type & kMin, real_type & kMax ) const {
    // cerco punto minimo parabola
    // root = -k/dk;
    kMin = CD.kappa0;
    kMax = CD.kappa(L);
    if ( kMax < kMin ) swap( kMax, kMin );
    return kMax - kMin;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::curvatureTotalVariation() const {
    // cerco punto minimo parabola
    // root = -k/dk;
    real_type km = CD.kappa0;
    real_type kp = CD.kappa(L);
    return abs(kp-km);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralCurvature2() const {
    return L*( CD.kappa0*(CD.kappa0+L*CD.dk) + (L*L)*CD.dk*CD.dk/3 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralJerk2() const {
    real_type k2 = CD.kappa0*CD.kappa0;
    real_type k3 = CD.kappa0*k2;
    real_type k4 = k2*k2;
    real_type t1 = L;
    real_type t2 = L*t1;
    real_type t3 = L*t2;
    real_type t4 = L*t3;
    return ((((t4/5*CD.dk+t3*CD.kappa0)*CD.dk+(1+2*t2)*k2)*CD.dk+2*t1*k3)*CD.dk+k4)*L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidCurve::integralSnap2() const {
    real_type k2  = CD.kappa0*CD.kappa0;
    real_type k3  = CD.kappa0*k2;
    real_type k4  = k2*k2;
    real_type k5  = k4*CD.kappa0;
    real_type k6  = k4*k2;
    real_type dk2 = CD.dk*CD.dk;
    real_type dk3 = CD.dk*dk2;
    real_type dk4 = dk2*dk2;
    real_type dk5 = dk4*CD.dk;
    real_type dk6 = dk4*dk2;
    real_type t2  = L;
    real_type t3  = L*t2;
    real_type t4  = L*t3;
    real_type t5  = L*t4;
    real_type t6  = L*t5;
    real_type t7  = L*t6;

    return ( (t7/7)*dk6 + dk5*CD.kappa0*t6 + 3*dk4*k2*t5 + 5*dk3*k3*t4 +
             5*dk2*k4*t3 + 3*dk3*t3 + 3*CD.dk*k5*t2 + 9*dk2*CD.kappa0*t2 +
             k6+9*k2*CD.dk ) * L;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ostream_type &
  operator << ( ostream_type & stream, ClothoidCurve const & c ) {
    stream <<   "x0     = " << c.CD.x0
           << "\ny0     = " << c.CD.y0
           << "\ntheta0 = " << c.CD.theta0
           << "\nkappa0 = " << c.CD.kappa0
           << "\ndk     = " << c.CD.dk
           << "\nL      = " << c.L
           << "\n";
    return stream;
  }

}

// EOF: Clothoid.cc
