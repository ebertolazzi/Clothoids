/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_eigen.hh
//
#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_EIGEN_dot_HH
#define UTILS_EIGEN_dot_HH

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wunused-template"
#endif

#ifdef _MSC_VER
#pragma warning( disable : 4127 )
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#ifndef EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_NO_AUTOMATIC_RESIZING
#endif

#include <type_traits>

#include "Utils/3rd/Eigen/Core"
#include "Utils/3rd/Eigen/Dense"
#include "Utils/3rd/Eigen/QR"
#include "Utils/3rd/Eigen/Sparse"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace fmt
{
  template <typename TYPE, int ROW, int COL>
  struct formatter<Eigen::Matrix<TYPE, ROW, COL>> : ostream_formatter
  {
  };

  template <typename PlainObjectType, int MapOptions, typename StrideType>
  struct formatter<Eigen::Map<PlainObjectType, MapOptions, StrideType>> : ostream_formatter
  {
  };

  template <typename MAT>
  struct formatter<Eigen::Transpose<MAT>> : ostream_formatter
  {
  };

  template <typename EXPR>
  struct formatter<Eigen::WithFormat<EXPR>> : ostream_formatter
  {
  };
}  // namespace fmt

#endif

#endif

#endif

//
// eof: Utils_eigen.hh
//
