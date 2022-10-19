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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

namespace G2lib {

  using std::map;
  using std::pair;

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

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  using Ppair = pair<CurveType,CurveType>;

  static map<Ppair,CurveType> const promote_map = {

    {Ppair( G2LIB_LINE, G2LIB_LINE ),          G2LIB_LINE},
    {Ppair( G2LIB_LINE, G2LIB_CIRCLE ),        G2LIB_CIRCLE},
    {Ppair( G2LIB_LINE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_LINE, G2LIB_BIARC ),         G2LIB_BIARC_LIST},
    {Ppair( G2LIB_LINE, G2LIB_BIARC_LIST ),    G2LIB_BIARC_LIST},
    {Ppair( G2LIB_LINE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_LINE, G2LIB_POLYLINE ),      G2LIB_POLYLINE},

    {Ppair( G2LIB_CIRCLE, G2LIB_LINE ),          G2LIB_CIRCLE},
    {Ppair( G2LIB_CIRCLE, G2LIB_CIRCLE ),        G2LIB_CIRCLE},
    {Ppair( G2LIB_CIRCLE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_CIRCLE, G2LIB_BIARC ),         G2LIB_BIARC_LIST},
    {Ppair( G2LIB_CIRCLE, G2LIB_BIARC_LIST ),    G2LIB_BIARC_LIST},
    {Ppair( G2LIB_CIRCLE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CIRCLE, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_BIARC, G2LIB_LINE ),          G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_BIARC ),         G2LIB_BIARC},
    {Ppair( G2LIB_BIARC, G2LIB_BIARC_LIST ),    G2LIB_BIARC_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_BIARC, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_CLOTHOID, G2LIB_LINE ),          G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CIRCLE ),        G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID},
    {Ppair( G2LIB_CLOTHOID, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID, G2LIB_BIARC_LIST ),    G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_LINE ),          G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_BIARC_LIST ),    G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_CLOTHOID_LIST, G2LIB_POLYLINE ),      G2LIB_CLOTHOID_LIST},

    {Ppair( G2LIB_POLYLINE, G2LIB_LINE ),          G2LIB_POLYLINE},
    {Ppair( G2LIB_POLYLINE, G2LIB_CIRCLE ),        G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_CLOTHOID ),      G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_BIARC ),         G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_BIARC_LIST ),    G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_CLOTHOID_LIST ), G2LIB_CLOTHOID_LIST},
    {Ppair( G2LIB_POLYLINE, G2LIB_POLYLINE ),      G2LIB_POLYLINE}
  };

  CurveType curve_promote( CurveType A, CurveType B ) {
    return promote_map.at(Ppair(A,B));
  }

  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  collision( BaseCurve const * pC1, BaseCurve const * pC2 ) {

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::collision {} vs {} ADDRS: {}, {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

    bool ok = false;

    try {

      CurveType CT = curve_promote( pC1->type(), pC2->type() );
      switch ( CT ) {
      case G2LIB_LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          ok = L1.collision( L2 );
        }
        break;
      case G2LIB_CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          ok = C1.collision( C2 );
        }
        break;
      case G2LIB_BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          ok = B1.collision( B2 );
        }
        break;
      case G2LIB_CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          ok = C1.collision( C2 );
        }
        break;
      case G2LIB_POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          ok = PL1.collision( PL2 );
        }
        break;
      case G2LIB_BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          ok = BL1.collision( BL2 );
        }
        break;
      case G2LIB_CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          ok = CL1.collision( CL2 );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::collision, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {

      std::cerr << "G2lib::collision error: " << e.what() << '\n';
      throw;

    } catch (...) {
      std::cerr << "G2lib::collision unknown error!\n";
      throw;
    }

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::collision {} vs {} ADDRS: {}, {} ok = {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2), ok
    );

    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  collision_ISO(
    BaseCurve const * pC1,
    real_type         offs1,
    BaseCurve const * pC2,
    real_type         offs2
  ) {

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::collision_ISO {} vs {} ADDRS: {}, {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

    bool ok = false;

    try {

      CurveType CT = curve_promote( pC1->type(), pC2->type() );
      switch ( CT ) {
      case G2LIB_LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          ok = L1.collision_ISO( offs1, L2, offs2 );
        }
        break;
      case G2LIB_CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          ok = C1.collision_ISO( offs1, C2, offs2 );
        }
        break;
      case G2LIB_BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          ok = B1.collision_ISO( offs1, B2, offs2 );
        }
        break;
      case G2LIB_CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          ok = C1.collision_ISO( offs1, C2, offs2 );
        }
        break;
      case G2LIB_POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          ok = PL1.collision_ISO( offs1, PL2, offs2 );
        }
        break;
      case G2LIB_BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          ok = BL1.collision_ISO( offs1, BL2, offs2 );
        }
        break;
      case G2LIB_CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          ok = CL1.collision_ISO( offs1, CL2, offs2 );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::collision_ISO, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {

      std::cerr << "G2lib::collision_ISO error: " << e.what() << '\n';
      throw;

    } catch (...) {
      std::cerr << "G2lib::collision_ISO unknown error!\n";
      throw;
    }

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::collision_ISO {} vs {} ADDRS: {}, {} ok = {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2), ok
    );

    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  intersect(
    BaseCurve const * pC1,
    BaseCurve const * pC2,
    IntersectList   & ilist
  ) {

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::intersect {} vs {} ADDRS: {}, {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

    try {

      CurveType CT = curve_promote( pC1->type(), pC2->type() );
      switch ( CT ) {
      case G2LIB_LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          L1.intersect( L2, ilist );
        }
        break;
      case G2LIB_CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          C1.intersect( C2, ilist );
        }
        break;
      case G2LIB_BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          B1.intersect( B2, ilist );
        }
        break;
      case G2LIB_CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          C1.intersect( C2, ilist );
        }
        break;
      case G2LIB_POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          PL1.intersect( PL2, ilist );
        }
        break;
      case G2LIB_BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          BL1.intersect( BL2, ilist  );
        }
        break;
      case G2LIB_CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          CL1.intersect( CL2, ilist );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::intersect, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {

      std::cerr << "G2lib::intersect error: " << e.what() << '\n';
      throw;

    } catch (...) {
      std::cerr << "G2lib::intersect unknown error!\n";
      throw;
    }

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::intersect {} vs {} ADDRS: {}, {} DONE\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  intersect_ISO(
    BaseCurve const * pC1, real_type offs1,
    BaseCurve const * pC2, real_type offs2,
    IntersectList   & ilist
  ) {

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::intersect_ISO {} vs {} ADDRS: {}, {}\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

    try {

      CurveType CT = curve_promote( pC1->type(), pC2->type() );
      switch ( CT ) {
      case G2LIB_LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          L1.intersect_ISO( offs1, L2, offs2, ilist );
        }
        break;
      case G2LIB_CIRCLE:
         G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          C1.intersect_ISO( offs1, C2, offs2, ilist );
        }
        break;
      case G2LIB_BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          B1.intersect_ISO( offs1, B2, offs2, ilist );
        }
        break;
      case G2LIB_CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          C1.intersect_ISO( offs1, C2, offs2, ilist );
        }
        break;
      case G2LIB_POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          PL1.intersect_ISO( offs1, PL2, offs2, ilist );
        }
        break;
      case G2LIB_BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          BL1.intersect_ISO( offs1, BL2, offs2, ilist );
        }
        break;
      case G2LIB_CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          CL1.intersect_ISO( offs1, CL2, offs2, ilist );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::intersect_ISO, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {

      std::cerr << "G2lib::intersect_ISO error: " << e.what() << '\n';
      throw;

    } catch (...) {
      std::cerr << "G2lib::intersect_ISO unknown error!\n";
      throw;
    }

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::intersect_ISO {} vs {} ADDRS: {}, {} DONE\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

  }
}

// EOF: G2lib_intersect.cc
