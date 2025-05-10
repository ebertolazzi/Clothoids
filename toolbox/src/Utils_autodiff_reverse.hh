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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_eigen.hh
//
#pragma once

#ifndef UTILS_AUTODIFF_REVERSE_dot_HH
#define UTILS_AUTODIFF_REVERSE_dot_HH

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils/3rd/autodiff/reverse/var.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace fmt {
  template <> struct formatter<autodiff::var> : ostream_formatter {};
}
#endif

#endif

//
// eof: Utils_eigen.hh
//
