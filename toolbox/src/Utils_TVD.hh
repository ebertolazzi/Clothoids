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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_TVD.hh
//

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

  //!
  //!  \brief Class for performing Total Variation Denoising (TVD) on 1D signals.
  //!
  //!  The `TVD` class implements algorithms for denoising signals by minimizing
  //!  total variation. The primary method provided is `denoise`, which applies
  //!  total variation denoising to a given input signal.
  //!
  //!
  //!  This algorithm minimizes the objective function defined as:
  //!
  //!  \f[
  //!    \textrm{minimize}\quad \sum_k (y_k - x_k)^2 + \lambda \sum_k |x_{k+1}-x_k|
  //!  \f]
  //!
  //!  where \f$ y_k \f$ is the input signal, \f$ x_k \f$
  //!  is the denoised signal, and \f$ \lambda \f$ is
  //!  a regularization parameter that controls the trade-off
  //!  between data fidelity and smoothness of the result.
  //!
  //!  \note
  //!  This implementation is based on the work of:
  //!
  //!    | Laurent Condat.
  //!    | A Direct Algorithm for 1D Total Variation Denoising.
  //!    | IEEE Signal Processing Letters,
  //!    | Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057.
  //!    | DOI: <10.1109/LSP.2013.2278339>.
  //!    | Available at: <https://hal.science/hal-00675043v4>
  //!
  //! \tparam Real The data type used for the computation (e.g., float, double).
  //!
  template <typename Real>
  class TVD {
    using Integer = int;
  public:

    //!
    //!  \brief Denoises a 1D signal using total variation denoising.
    //!
    //!  This method is a convenience wrapper that calls the more detailed `denoise` method
    //!  with a default increment of 1 for both the input and output arrays.
    //!
    //!  \param N      The number of elements in the input signal.
    //!  \param y      A pointer to the input signal array (size N).
    //!  \param lambda The regularization parameter that controls the smoothness.
    //!  \param x      A pointer to the output denoised signal array (size N).
    //!
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

    //!
    //! \brief Performs total variation denoising on a 1D signal.
    //!
    //! This method minimizes the objective function defined as:
    //!
    //! \f[
    //!    \textrm{minimize}\quad \sum_k (y_k - x_k)^2 + \lambda \sum_k |x_{k+1}-x_k|
    //! \f]
    //!
    //! where \f$ y_k \f$ is the input signal, \f$ x_k \f$
    //! is the denoised signal, and \f$ \lambda \f$ is
    //! a regularization parameter that controls the trade-off
    //! between data fidelity and smoothness of the result.
    //!
    //! \param N      The number of elements in the input signal.
    //! \param y      A pointer to the input signal array (size N).
    //! \param incy   The increment between consecutive elements in the input array (y).
    //! \param lambda The regularization parameter that controls the smoothness.
    //! \param x      A pointer to the output denoised signal array (size N).
    //! \param incx   The increment between consecutive elements in the output array (x).
    //!
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

//
// eof: Utils_TVD.hh
//
