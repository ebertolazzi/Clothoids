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

#ifdef NO_SYSTEM_UTILS
  #include "Utils.hh"
  #include "Utils_AABB_tree.hh"
#else
  #include <Utils.hh>
  #include <Utils_AABB_tree.hh>
#endif

#include "GenericContainer/GenericContainer.hh"

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

  using istream_type     = std::basic_istream<char>;             //!< input streaming
  using ostream_type     = std::basic_ostream<char>;             //!< output streaming
  using real_type        = double;                               //!< real type number
  using integer          = int;                                  //!< integer type number
  using AABB_TREE        = Utils::AABBtree<real_type>;           //!< `AABB` tree type
  using AABB_SET         = Utils::AABBtree<real_type>::AABB_SET; //!< Set type used in `AABB` tree object
  using AABB_MAP         = Utils::AABBtree<real_type>::AABB_MAP; //!< Map type used in `AABB` tree object
  using GenericContainer = GC_namespace::GenericContainer;       //!< Generic container object

  //!
  //! Enumeration type for curve type
  //!
  using CurveType = enum class CurveType : integer {
    LINE,
    POLYLINE,
    CIRCLE,
    BIARC,
    BIARC_LIST,
    CLOTHOID,
    CLOTHOID_LIST,
    DUBINS,
    DUBINS3P
  };

  //!
  //! Convert curve type to a string
  //!
  inline
  string
  to_string( CurveType n ) {
    string res{""};
    switch ( n ) {
    case CurveType::LINE:          res = "LINE";          break;
    case CurveType::POLYLINE:      res = "POLYLINE";      break;
    case CurveType::CIRCLE:        res = "CIRCLE";        break;
    case CurveType::BIARC:         res = "BIARC";         break;
    case CurveType::BIARC_LIST:    res = "BIARC_LIST";    break;
    case CurveType::CLOTHOID:      res = "CLOTHOID";      break;
    case CurveType::CLOTHOID_LIST: res = "CLOTHOID_LIST"; break;
    case CurveType::DUBINS:        res = "DUBINS";        break;
    case CurveType::DUBINS3P:      res = "DUBINS3P";      break;
    }
    return res;
  };

  //!
  //! Given two curve type determine curve type that cointain both type
  //!
  //! \param[in] A first curve type
  //! \param[in] B second curve type
  //! \return the curve type super type of both
  //!
  extern CurveType curve_promote( CurveType A, CurveType B );

  class LineSegment;
  class CircleArc;
  class Biarc;
  class ClothoidCurve;
  class PolyLine;
  class BiarcList;
  class ClothoidList;
  class Dubins;
  class Dubins3p;
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
#include "Clothoids/Dubins3p.hxx"

#endif

///
/// eof: Clothoids.hh
///
