/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2018                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Paolo Bevilacqua and Enrico Bertolazzi                              |
 |                                                                          |
 |      (1) Dipartimento di Ingegneria e Scienza dell'Informazione          |
 |      (2) Dipartimento di Ingegneria Industriale                          |
 |                                                                          |
 |      Universita` degli Studi di Trento                                   |
 |      email: paolo.bevilacqua@unitn.it                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Clothoids.hh
///

#pragma once

#ifndef CLOTHOIDS_dot_HH
#define CLOTHOIDS_dot_HH

// comment to disable threads support
#define CLOTHOIDS_USE_THREADS 1

#include "Utils.hh"
#include "Utils_AABB_tree.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iterator>

#include <memory>  // shared_ptr

#ifdef G2LIB_DEBUG
  #define G2LIB_DEBUG_MESSAGE(...) std::cout << fmt::format(__VA_ARGS__) << std::flush
#else
  #define G2LIB_DEBUG_MESSAGE(...)
#endif

#ifndef GLIB2_TOL_ANGLE
  #define GLIB2_TOL_ANGLE 1e-8
#endif

namespace G2lib {

  using std::string;
  using std::vector;
  using std::map;
  using std::set;

  using istream_type = std::basic_istream<char>;
  using ostream_type = std::basic_ostream<char>;
  using real_type    = double;
  using integer      = int;
  using AABB_TREE    = Utils::AABBtree<real_type>;
  using AABB_SET     = Utils::AABBtree<real_type>::AABB_SET;
  using AABB_MAP     = Utils::AABBtree<real_type>::AABB_MAP;

  using CurveType = enum class CurveType : integer {
    LINE,
    POLYLINE,
    CIRCLE,
    BIARC,
    BIARC_LIST,
    CLOTHOID,
    CLOTHOID_LIST
  };

  inline
  string
  to_string( CurveType n ) {
    string res = "";
    switch ( n ) {
    case CurveType::LINE:          res = "LINE";          break;
    case CurveType::POLYLINE:      res = "POLYLINE";      break;
    case CurveType::CIRCLE:        res = "CIRCLE";        break;
    case CurveType::BIARC:         res = "BIARC";         break;
    case CurveType::BIARC_LIST:    res = "BIARC_LIST";    break;
    case CurveType::CLOTHOID:      res = "CLOTHOID";      break;
    case CurveType::CLOTHOID_LIST: res = "CLOTHOID_LIST"; break;
    }
    return res;
  };

  extern CurveType curve_promote( CurveType, CurveType );

  class LineSegment;
  class CircleArc;
  class Biarc;
  class ClothoidCurve;
  class PolyLine;
  class BiarcList;
  class ClothoidList;
}

#include "Clothoids/G2lib.hxx"
#include "Clothoids/Triangle2D.hxx"
#include "Clothoids/BBox.hxx"
#include "Clothoids/BaseCurve.hxx"
#include "Clothoids/Fresnel.hxx"
#include "Clothoids/Line.hxx"
#include "Clothoids/Circle.hxx"
#include "Clothoids/Biarc.hxx"
#include "Clothoids/Clothoid.hxx"
#include "Clothoids/PolyLine.hxx"
#include "Clothoids/BiarcList.hxx"
#include "Clothoids/ClothoidList.hxx"
#include "Clothoids/ClothoidAsyPlot.hxx"
#include "Clothoids/Dubins.hxx"

namespace fmt {
  template <> struct formatter<G2lib::Biarc>            : ostream_formatter {};
  template <> struct formatter<G2lib::BiarcList>        : ostream_formatter {};
  template <> struct formatter<G2lib::BBox>             : ostream_formatter {};
  template <> struct formatter<G2lib::CircleArc>        : ostream_formatter {};
  template <> struct formatter<G2lib::ClothoidCurve>    : ostream_formatter {};
  template <> struct formatter<G2lib::ClothoidSplineG2> : ostream_formatter {};
  template <> struct formatter<G2lib::ClothoidList>     : ostream_formatter {};
  template <> struct formatter<G2lib::LineSegment>      : ostream_formatter {};
  template <> struct formatter<G2lib::PolyLine>         : ostream_formatter {};
  template <> struct formatter<G2lib::Triangle2D>       : ostream_formatter {};
}

namespace G2lib {

  using std::string;
  using std::vector;
  using std::map;
  using std::set;

  using istream_type = std::basic_istream<char>;
  using ostream_type = std::basic_ostream<char>;
  using real_type    = double;
  using integer      = int;
  using AABB_TREE    = Utils::AABBtree<real_type>;
  using AABB_SET     = Utils::AABBtree<real_type>::AABB_SET;
  using AABB_MAP     = Utils::AABBtree<real_type>::AABB_MAP;

  extern CurveType curve_promote( CurveType, CurveType );

  class LineSegment;
  class CircleArc;
  class Biarc;
  class ClothoidCurve;
  class PolyLine;
  class BiarcList;
  class ClothoidList;
}

#endif

///
/// eof: Clothoids.hh
///
