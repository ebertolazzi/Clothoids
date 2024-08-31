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

  using Ppair = pair<CurveType,CurveType>; //!< Pair of curve type

  static map<Ppair,CurveType> const promote_map = {

    {Ppair( CurveType::LINE, CurveType::LINE ),          CurveType::LINE},
    {Ppair( CurveType::LINE, CurveType::CIRCLE ),        CurveType::CIRCLE},
    {Ppair( CurveType::LINE, CurveType::CLOTHOID ),      CurveType::CLOTHOID},
    {Ppair( CurveType::LINE, CurveType::BIARC ),         CurveType::BIARC_LIST},
    {Ppair( CurveType::LINE, CurveType::BIARC_LIST ),    CurveType::BIARC_LIST},
    {Ppair( CurveType::LINE, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::LINE, CurveType::POLYLINE ),      CurveType::POLYLINE},
    {Ppair( CurveType::LINE, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::CIRCLE, CurveType::LINE ),          CurveType::CIRCLE},
    {Ppair( CurveType::CIRCLE, CurveType::CIRCLE ),        CurveType::CIRCLE},
    {Ppair( CurveType::CIRCLE, CurveType::CLOTHOID ),      CurveType::CLOTHOID},
    {Ppair( CurveType::CIRCLE, CurveType::BIARC ),         CurveType::BIARC_LIST},
    {Ppair( CurveType::CIRCLE, CurveType::BIARC_LIST ),    CurveType::BIARC_LIST},
    {Ppair( CurveType::CIRCLE, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CIRCLE, CurveType::POLYLINE ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CIRCLE, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::BIARC, CurveType::LINE ),          CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::BIARC, CurveType::CIRCLE ),        CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::BIARC, CurveType::CLOTHOID ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::BIARC, CurveType::BIARC ),         CurveType::BIARC},
    {Ppair( CurveType::BIARC, CurveType::BIARC_LIST ),    CurveType::BIARC_LIST},
    {Ppair( CurveType::BIARC, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::BIARC, CurveType::POLYLINE ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::BIARC, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::CLOTHOID, CurveType::LINE ),          CurveType::CLOTHOID},
    {Ppair( CurveType::CLOTHOID, CurveType::CIRCLE ),        CurveType::CLOTHOID},
    {Ppair( CurveType::CLOTHOID, CurveType::CLOTHOID ),      CurveType::CLOTHOID},
    {Ppair( CurveType::CLOTHOID, CurveType::BIARC ),         CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID, CurveType::BIARC_LIST ),    CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID, CurveType::POLYLINE ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::CLOTHOID_LIST, CurveType::LINE ),          CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::CIRCLE ),        CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::CLOTHOID ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::BIARC ),         CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::BIARC_LIST ),    CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::POLYLINE ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::CLOTHOID_LIST, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::POLYLINE, CurveType::LINE ),          CurveType::POLYLINE},
    {Ppair( CurveType::POLYLINE, CurveType::CIRCLE ),        CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::POLYLINE, CurveType::CLOTHOID ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::POLYLINE, CurveType::BIARC ),         CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::POLYLINE, CurveType::BIARC_LIST ),    CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::POLYLINE, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::POLYLINE, CurveType::POLYLINE ),      CurveType::POLYLINE},
    {Ppair( CurveType::POLYLINE, CurveType::DUBINS ),        CurveType::CLOTHOID_LIST},

    {Ppair( CurveType::DUBINS, CurveType::LINE ),          CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::CIRCLE ),        CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::CLOTHOID ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::BIARC ),         CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::BIARC_LIST ),    CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::CLOTHOID_LIST ), CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::POLYLINE ),      CurveType::CLOTHOID_LIST},
    {Ppair( CurveType::DUBINS, CurveType::DUBINS ),        CurveType::DUBINS}
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
      case CurveType::LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          ok = L1.collision( L2 );
        }
        break;
      case CurveType::CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          ok = C1.collision( C2 );
        }
        break;
      case CurveType::BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          ok = B1.collision( B2 );
        }
        break;
      case CurveType::CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          ok = C1.collision( C2 );
        }
        break;
      case CurveType::POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          ok = PL1.collision( PL2 );
        }
        break;
      case CurveType::BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          ok = BL1.collision( BL2 );
        }
        break;
      case CurveType::CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          ok = CL1.collision( CL2 );
        }
        break;
      case CurveType::DUBINS:
        {
          Dubins const & DB1{*static_cast<Dubins const *>(pC1)};
          Dubins const & DB2{*static_cast<Dubins const *>(pC2)};
          ok = DB1.collision( DB2 );
        }
        break;
      case CurveType::DUBINS3P:
        {
          Dubins3p const & DB1{*static_cast<Dubins3p const *>(pC1)};
          Dubins3p const & DB2{*static_cast<Dubins3p const *>(pC2)};
          ok = DB1.collision( DB2 );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::collision, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {
      G2LIB_DEBUG_MESSAGE( "G2lib::collision error: {}\n", e.what() );
      UTILS_ERROR( "G2lib::collision error: {}\n", e.what() );
    } catch (...) {
      G2LIB_DEBUG_MESSAGE( "G2lib::collision unknown error!\n" );
      UTILS_ERROR( "G2lib::collision unknown error\n" );
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
      case CurveType::LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          ok = L1.collision_ISO( offs1, L2, offs2 );
        }
        break;
      case CurveType::CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          ok = C1.collision_ISO( offs1, C2, offs2 );
        }
        break;
      case CurveType::BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          ok = B1.collision_ISO( offs1, B2, offs2 );
        }
        break;
      case CurveType::CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          ok = C1.collision_ISO( offs1, C2, offs2 );
        }
        break;
      case CurveType::POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          ok = PL1.collision_ISO( offs1, PL2, offs2 );
        }
        break;
      case CurveType::BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          ok = BL1.collision_ISO( offs1, BL2, offs2 );
        }
        break;
      case CurveType::CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          ok = CL1.collision_ISO( offs1, CL2, offs2 );
        }
        break;
      case CurveType::DUBINS:
        {
          Dubins const & DB1{*static_cast<Dubins const *>(pC1)};
          Dubins const & DB2{*static_cast<Dubins const *>(pC2)};
          ok = DB1.collision_ISO( offs1, DB2, offs2 );
        }
        break;
      case CurveType::DUBINS3P:
        {
          Dubins3p const & DB1{*static_cast<Dubins3p const *>(pC1)};
          Dubins3p const & DB2{*static_cast<Dubins3p const *>(pC2)};
          ok = DB1.collision_ISO( offs1, DB2, offs2 );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::collision_ISO, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {
      G2LIB_DEBUG_MESSAGE( "G2lib::collision_ISO error: {}\n", e.what() );
      UTILS_ERROR( "G2lib::collision_ISO error: {}\n", e.what() );
    } catch (...) {
      G2LIB_DEBUG_MESSAGE( "G2lib::collision_ISO unknown error!\n" );
      UTILS_ERROR( "G2lib::collision_ISO unknown error\n" );
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
      case CurveType::LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          L1.intersect( L2, ilist );
        }
        break;
      case CurveType::CIRCLE:
        G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          C1.intersect( C2, ilist );
        }
        break;
      case CurveType::BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          B1.intersect( B2, ilist );
        }
        break;
      case CurveType::CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          C1.intersect( C2, ilist );
        }
        break;
      case CurveType::POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          PL1.intersect( PL2, ilist );
        }
        break;
      case CurveType::BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          BL1.intersect( BL2, ilist  );
        }
        break;
      case CurveType::CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          CL1.intersect( CL2, ilist );
        }
        break;
      case CurveType::DUBINS:
        {
          Dubins const & DB1{*static_cast<Dubins const *>(pC1)};
          Dubins const & DB2{*static_cast<Dubins const *>(pC2)};
          DB1.intersect( DB2, ilist );
        }
        break;
      case CurveType::DUBINS3P:
        {
          Dubins3p const & DB1{*static_cast<Dubins3p const *>(pC1)};
          Dubins3p const & DB2{*static_cast<Dubins3p const *>(pC2)};
          DB1.intersect( DB2, ilist );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::intersect, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {
      G2LIB_DEBUG_MESSAGE( "G2lib::intersect error: {}\n", e.what() );
      UTILS_ERROR( "G2lib::intersect error: {}\n", e.what() );
    } catch (...) {
      G2LIB_DEBUG_MESSAGE( "G2lib::intersect unknown error!\n" );
      UTILS_ERROR( "G2lib::intersect unknown error\n" );
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
      case CurveType::LINE:
        G2LIB_DEBUG_MESSAGE( "promote -> LineSegment\n" );
        {
          LineSegment L1( pC1 );
          LineSegment L2( pC2 );
          L1.intersect_ISO( offs1, L2, offs2, ilist );
        }
        break;
      case CurveType::CIRCLE:
         G2LIB_DEBUG_MESSAGE( "promote -> CircleArc\n" );
        {
          CircleArc C1( pC1 );
          CircleArc C2( pC2 );
          C1.intersect_ISO( offs1, C2, offs2, ilist );
        }
        break;
      case CurveType::BIARC:
        G2LIB_DEBUG_MESSAGE( "promote -> Biarc\n" );
        {
          Biarc B1( pC1 );
          Biarc B2( pC2 );
          B1.intersect_ISO( offs1, B2, offs2, ilist );
        }
        break;
      case CurveType::CLOTHOID:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidCurve\n" );
        {
          ClothoidCurve C1( pC1 );
          ClothoidCurve C2( pC2 );
          C1.intersect_ISO( offs1, C2, offs2, ilist );
        }
        break;
      case CurveType::POLYLINE:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          PolyLine PL1( pC1 );
          PolyLine PL2( pC2 );
          PL1.intersect_ISO( offs1, PL2, offs2, ilist );
        }
        break;
      case CurveType::BIARC_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> PolyLine\n" );
        {
          BiarcList BL1( pC1 );
          BiarcList BL2( pC2 );
          BL1.intersect_ISO( offs1, BL2, offs2, ilist );
        }
        break;
      case CurveType::CLOTHOID_LIST:
        G2LIB_DEBUG_MESSAGE( "promote -> ClothoidList\n" );
        {
          ClothoidList CL1( pC1 );
          ClothoidList CL2( pC2 );
          CL1.intersect_ISO( offs1, CL2, offs2, ilist );
        }
        break;
      case CurveType::DUBINS:
        {
          Dubins const & DB1{*static_cast<Dubins const *>(pC1)};
          Dubins const & DB2{*static_cast<Dubins const *>(pC2)};
          DB1.intersect_ISO( offs1, DB2, offs2, ilist );
        }
        break;
      case CurveType::DUBINS3P:
        {
          Dubins3p const & DB1{*static_cast<Dubins3p const *>(pC1)};
          Dubins3p const & DB2{*static_cast<Dubins3p const *>(pC2)};
          DB1.intersect_ISO( offs1, DB2, offs2, ilist );
        }
        break;
      //default:
      //  UTILS_ERROR0( "G2lib::intersect_ISO, missing curve type" );
      //  break;
      }

    } catch ( std::exception const & e ) {
      G2LIB_DEBUG_MESSAGE( "G2lib::intersect_ISO error: {}\n", e.what() );
      UTILS_ERROR( "G2lib::intersect_ISO error: {}\n", e.what() );
    } catch (...) {
      G2LIB_DEBUG_MESSAGE( "G2lib::intersect_ISO unknown error!\n" );
      UTILS_ERROR( "G2lib::intersect_ISO unknown error\n" );
    }

    G2LIB_DEBUG_MESSAGE(
      "G2Lib::intersect_ISO {} vs {} ADDRS: {}, {} DONE\n",
      pC1->type_name(), pC2->type_name(), fmt::ptr(pC1), fmt::ptr(pC2)
    );

  }
}

// EOF: G2lib_intersect.cc
