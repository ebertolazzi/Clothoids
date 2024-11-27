/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
/// file: Utils_TVD.hh
///

#pragma once

#ifndef UTILS_TVD_dot_HH
#define UTILS_TVD_dot_HH

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "Utils.hh"

namespace Utils {

  template <typename Real>
  class TVD {
    using Integer = int;
  public:

    static
    void
    denoise(
      Integer    N,
      Real const y[],
      Real       lambda,
      Real       x[]
    ) {
      denoise( N, y, 1, lambda, x, 1 );
    }

    static
    void
    denoise(
      Integer    N,
      Real const y[],
      Integer    incy,
      Real       lambda,
      Real       x[],
      Integer    incx
    );
  };

  #ifndef UTILS_OS_WINDOWS
  extern template class TVD<float>;
  extern template class TVD<double>;
  #endif

}

#endif

///
/// EOF: Utils_TVD.hh
///
