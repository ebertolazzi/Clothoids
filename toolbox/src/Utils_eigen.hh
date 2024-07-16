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

///
/// file: Utils_eigen.hh
///
#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_EIGEN_dot_HH
#define UTILS_EIGEN_dot_HH

#include "Utils.hh"

#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated-copy-with-dtor"
#endif

#ifdef _MSC_VER
  #pragma warning( push )
  #pragma warning( disable : 4127 )
#endif

#ifndef EIGEN_DONT_PARALLELIZE
  #define EIGEN_DONT_PARALLELIZE
#endif

#ifndef EIGEN_NO_AUTOMATIC_RESIZING
  #define EIGEN_NO_AUTOMATIC_RESIZING
#endif

#include "Eigen/Core"
#include "Eigen/Dense"
#include <type_traits>

namespace fmt {
  template <typename TYPE, int ROW, int COL>
  struct formatter<Eigen::Matrix<TYPE,ROW,COL>> : ostream_formatter {};

  template <typename PlainObjectType, int MapOptions, typename StrideType>
  struct formatter<Eigen::Map<PlainObjectType,MapOptions,StrideType>> : ostream_formatter {};

  template <typename MAT>
  struct formatter<Eigen::Transpose<MAT>> : ostream_formatter {};

  template <typename EXPR>
  struct formatter<Eigen::WithFormat<EXPR>> : ostream_formatter {};
}

#ifdef __clang__
  #pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
  #pragma warning( pop )
#endif

#endif

#endif

///
/// eof: Utils_eigen.hh
///
