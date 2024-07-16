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
/// file: Clothoids_fmt.hh
///

#pragma once

#ifndef CLOTHOIDS_FMT_dot_HH
#define CLOTHOIDS_FMT_dot_HH

#ifdef NO_SYSTEM_UTILS
  #include "Utils_fmt.hh"
#else
  #include <Utils_fmt.hh>
#endif

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
  template <> struct formatter<G2lib::Dubins>           : ostream_formatter {};
  template <> struct formatter<G2lib::Dubins3p>         : ostream_formatter {};
}

#endif

///
/// eof: Clothoids_fmt.hh
///
