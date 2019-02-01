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

#include "G2lib.hh"
#include "Line.hh"
#include "PolyLine.hh"
#include "Circle.hh"
#include "Clothoid.hh"
#include "ClothoidList.hh"

#include <map>
#include <algorithm>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

namespace G2lib {

  using std::map;
  using std::pair;

  using std::numeric_limits;
  using std::fpclassify;
  using std::lower_bound;
  using std::abs;
  using std::sqrt;

  /*\
   |   _       _                          _
   |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
   |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
   |  | | | | | ||  __/ |  \__ \  __/ (__| |_
   |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
  \*/

  typedef pair<CurveType,CurveType> Ppair;

  static map<Ppair,CurveType> const promote_map = {
    {Ppair( G2LIB_LINE, G2LIB_LINE ),          G2LIB_LINE},
    {Ppair( G2LIB_LINE, G2LIB_CIRCLE ),        G2LIB_CIRCLE},
    {Ppair( G2LIB_LINE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_LINE, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_LINE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_LINE, G2LIB_POLYLINE ),      G2LIB_POLYLINE},

    {Ppair( G2LIB_CIRCLE, G2LIB_LINE ),          G2LIB_CIRCLE},
    {Ppair( G2LIB_CIRCLE, G2LIB_CIRCLE ),        G2LIB_CIRCLE},
    {Ppair( G2LIB_CIRCLE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_CIRCLE, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CIRCLE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CIRCLE, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_BIARC, G2LIB_LINE ),          G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_BIARC ),         G2LIB_BIARC},
    {Ppair( G2LIB_BIARC, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_CLOTHOID, G2LIB_LINE ),          G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CIRCLE ),        G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_LINE ),          G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_POLYLINE, G2LIB_LINE ),          G2LIB_POLYLINE},
    {Ppair( G2LIB_POLYLINE, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_POLYLINE ),      G2LIB_POLYLINE}
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  collision( BaseCurve const & obj1,
             BaseCurve const & obj2 ) {
    #ifdef DEBUG
    std::cout << "collision " << CurveType_name[obj1.type()]
              << " with " << CurveType_name[obj2.type()]
              << " using " << CurveType_name[promote_map.at(Ppair(obj1.type(),obj2.type()))]
              << '\n';
    #endif

    bool ok = false;
    switch ( promote_map.at(Ppair(obj1.type(),obj2.type())) ) {
    case G2LIB_LINE:
      {
        LineSegment L1( obj1 );
        LineSegment L2( obj2 );
        ok = L1.collision( L2 );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C1( obj1 );
        CircleArc C2( obj2 );
        ok = C1.collision( C2 );
      }
      break;
    case G2LIB_CLOTHOID:
      {
        ClothoidCurve C1( obj1 );
        ClothoidCurve C2( obj2 );
        ok = C1.collision( C2 );
      }
      break;
    case G2LIB_BIARC:
      {
        Biarc const & B1 = *static_cast<Biarc const *>(&obj1);
        Biarc const & B2 = *static_cast<Biarc const *>(&obj2);
        ok = B1.collision( B2 );
      }
      break;
    case G2LIB_CLOTHOID_LIST:
      {
        ClothoidList C1( obj1 );
        ClothoidList C2( obj2 );
        ok = C1.collision(C2);
      }
      break;
    case G2LIB_POLYLINE:
      {
        PolyLine C1( obj1 );
        PolyLine C2( obj2 );
        ok = C1.collision( C2 );
      }
      break;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  collision( BaseCurve const & obj1,
             real_type         offs1,
             BaseCurve const & obj2,
             real_type         offs2 ) {
    #ifdef DEBUG
    std::cout << "collision (offs) " << CurveType_name[obj1.type()]
              << " with " << CurveType_name[obj2.type()]
              << " using " << CurveType_name[promote_map.at(Ppair(obj1.type(),obj2.type()))]
              << '\n';
    #endif

    bool ok = false;
    switch ( promote_map.at(Ppair(obj1.type(),obj2.type())) ) {
    case G2LIB_LINE:
      {
        LineSegment L1( obj1 );
        LineSegment L2( obj2 );
        ok = L1.collision( offs1, L2, offs2 );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C1( obj1 );
        CircleArc C2( obj2 );
        ok = C1.collision( offs1, C2, offs2 );
      }
      break;
    case G2LIB_CLOTHOID:
      {
        ClothoidCurve C1( obj1 );
        ClothoidCurve C2( obj2 );
        ok = C1.collision( offs1, C2, offs2 );
      }
      break;
    case G2LIB_BIARC:
      {
        Biarc const & B1 = *static_cast<Biarc const *>(&obj1);
        Biarc const & B2 = *static_cast<Biarc const *>(&obj2);
        ok = B1.collision( offs1, B2, offs2 );
      }
      break;
    case G2LIB_CLOTHOID_LIST:
      {
        ClothoidList C1( obj1 );
        ClothoidList C2( obj2 );
        ok = C1.collision( offs1, C2, offs2 );
      }
      break;
    case G2LIB_POLYLINE:
      {
        PolyLine C1( obj1 );
        PolyLine C2( obj2 );
        ok = C1.collision( offs1, C2, offs2 );
      }
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  intersect( BaseCurve const & obj1,
             BaseCurve const & obj2,
             IntersectList   & ilist,
             bool              swap_s_vals ) {

    #ifdef DEBUG
    std::cout << "intersect " << CurveType_name[obj1.type()]
              << " with " << CurveType_name[obj2.type()]
              << " using " << CurveType_name[promote_map.at(Ppair(obj1.type(),obj2.type()))]
              << '\n';
    #endif

    switch ( promote_map.at(Ppair(obj1.type(),obj2.type())) ) {
    case G2LIB_LINE:
      {
        LineSegment L1( obj1 );
        LineSegment L2( obj2 );
        L1.intersect( L2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C1( obj1 );
        CircleArc C2( obj2 );
        C1.intersect( C2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CLOTHOID:
      {
        ClothoidCurve C1( obj1 );
        ClothoidCurve C2( obj2 );
        C1.intersect( C2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_BIARC:
      {
        Biarc const & B1 = *static_cast<Biarc const *>(&obj1);
        Biarc const & B2 = *static_cast<Biarc const *>(&obj2);
        B1.intersect( B2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CLOTHOID_LIST:
      {
        ClothoidList C1( obj1 );
        ClothoidList C2( obj2 );
        C1.intersect( C2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_POLYLINE:
      {
        PolyLine C1( obj1 );
        PolyLine C2( obj2 );
        C1.intersect( C2, ilist, swap_s_vals );
      }
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  intersect( BaseCurve const & obj1,
             real_type         offs1,
             BaseCurve const & obj2,
             real_type         offs2,
             IntersectList   & ilist,
             bool              swap_s_vals ) {

    #ifdef DEBUG
    std::cout << "intersect (offs) " << CurveType_name[obj1.type()]
              << " with " << CurveType_name[obj2.type()]
              << " using " << CurveType_name[promote_map.at(Ppair(obj1.type(),obj2.type()))]
              << '\n';
    #endif

    switch ( promote_map.at(Ppair(obj1.type(),obj2.type())) ) {
    case G2LIB_LINE:
      {
        LineSegment L1( obj1 );
        LineSegment L2( obj2 );
        L1.intersect( offs1, L2, offs2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CIRCLE:
      {
        CircleArc C1( obj1 );
        CircleArc C2( obj2 );
        C1.intersect( offs1, C2, offs2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CLOTHOID:
      {
        ClothoidCurve C1( obj1 );
        ClothoidCurve C2( obj2 );
        C1.intersect( offs1, C2, offs2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_BIARC:
      {
        Biarc const & B1 = *static_cast<Biarc const *>(&obj1);
        Biarc const & B2 = *static_cast<Biarc const *>(&obj2);
        B1.intersect( offs1, B2, offs2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_CLOTHOID_LIST:
      {
        ClothoidList C1( obj1 );
        ClothoidList C2( obj2 );
        C1.intersect( offs1, C2, offs2, ilist, swap_s_vals );
      }
      break;
    case G2LIB_POLYLINE:
      {
        PolyLine C1( obj1 );
        PolyLine C2( obj2 );
        C1.intersect( offs1, C2, offs2, ilist, swap_s_vals );
      }
      break;
    }
  }
}

// EOF: G2lib_intersect.cc
