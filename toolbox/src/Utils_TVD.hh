/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022-2023                                                 |
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Utils.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  /*!
   * \brief Total Variation Denoising (TVD) for 1D signals
   *
   * This class implements the efficient O(N) algorithm for 1D Total Variation
   * Denoising (TVD) as described by Laurent Condat. The algorithm minimizes
   * the objective function:
   *
   * \f[
   *   \text{minimize} \quad \frac{1}{2} \sum_{i=1}^{N} (y_i - x_i)^2 +
   *   \lambda \sum_{i=1}^{N-1} |x_{i+1} - x_i|
   * \f]
   *
   * where:
   * - \f$ y \f$ is the input noisy signal
   * - \f$ x \f$ is the denoised output signal
   * - \f$ \lambda \f$ is the regularization parameter controlling smoothness
   *
   * \tparam Real Floating-point type (float, double, etc.)
   *
   * \note
   * This implementation is based on:
   * - Laurent Condat, "A Direct Algorithm for 1D Total Variation Denoising",
   *   IEEE Signal Processing Letters, 2013, 20(11), pp.1054-1057.
   * - DOI: \c 10.1109/LSP.2013.2278339
   * - Available at: https://hal.science/hal-00675043v4
   *
   * \see
   * For mathematical background on Total Variation regularization:
   * - Rudin, Osher, Fatemi, "Nonlinear total variation based noise removal
   *   algorithms", Physica D, 1992.
   * - Chambolle, "An algorithm for total variation minimization and
   *   applications", Journal of Mathematical Imaging and Vision, 2004.
   */
  template <typename Real> class TVD
  {
    using Integer = int;

  public:
    /*!
     * \brief Denoise a 1D signal with default stride
     *
     * Convenience wrapper for denoising with unit stride.
     *
     * \param[in]  N      Number of elements in the signal
     * \param[in]  y      Input noisy signal (array of size N)
     * \param[in]  lambda Regularization parameter (≥ 0)
     *                    - Larger λ → smoother output
     *                    - λ = 0 → output equals input
     * \param[out] x      Output denoised signal (array of size N)
     *
     * \pre N ≥ 1
     * \pre λ ≥ 0
     * \pre y and x must not overlap (unless y == x)
     */
    static void denoise( Integer N, Real const y[], Real lambda, Real x[] ) { denoise( N, y, 1, lambda, x, 1 ); }

    /*!
     * \brief Denoise a 1D signal with custom strides
     *
     * Implements the Condat algorithm for 1D TV denoising with O(N) complexity.
     * The algorithm maintains two solution paths (vmin, vmax) and updates
     * them according to the optimality conditions.
     *
     * \param[in]  N      Number of elements in the signal
     * \param[in]  y      Input noisy signal
     * \param[in]  incy   Stride between consecutive elements in y
     * \param[in]  lambda Regularization parameter (≥ 0)
     * \param[out] x      Output denoised signal
     * \param[in]  incx   Stride between consecutive elements in x
     *
     * \note
     * Algorithm overview:
     * 1. Initialize vmin = y₀ - λ, vmax = y₀ + λ
     * 2. Scan signal maintaining optimality conditions
     * 3. When conditions violated, output constant segment
     * 4. Repeat until entire signal processed
     *
     * \complexity O(N) time, O(1) additional space
     */
    static void denoise( Integer N, Real const y[], Integer incy, Real lambda, Real x[], Integer incx )
    {
      // Input validation (in release mode, these are no-ops if NDEBUG is defined)
      UTILS_ASSERT( N >= 1, "TVD::denoise: N must be positive, got {}", N );
      UTILS_ASSERT( lambda >= 0, "TVD::denoise: lambda must be non-negative, got {}", lambda );
      UTILS_ASSERT( incy != 0, "TVD::denoise: incy cannot be zero" );
      UTILS_ASSERT( incx != 0, "TVD::denoise: incx cannot be zero" );

      // Initialize output to zero
      for ( Integer j = 0; j < N; ++j ) x[j * incx] = 0;

      // Algorithm variables (following Condat's notation)
      Real vmin = y[0] - lambda;  // Lower solution path
      Real vmax = y[0] + lambda;  // Upper solution path
      Real umin = +lambda;        // Cumulative error for vmin
      Real umax = -lambda;        // Cumulative error for vmax

      Integer k  = 0;  // Current index
      Integer k0 = 0;  // Start of current segment
      Integer km = 0;  // Last index where umin reached +λ
      Integer kp = 0;  // Last index where umax reached -λ

      Integer const nm1 = N - 1;

      while ( true )
      {
        if ( k >= nm1 )
        {
          // End of signal reached
          x[nm1 * incx] = vmin + umin;
          break;
        }

        // Check optimality conditions
        if ( y[( k + 1 ) * incy] + umin < vmin - lambda )
        {
          // Condition (b1): vmin is too high
          for ( Integer j = k0; j <= km; ++j ) x[j * incx] = vmin;

          // Restart from next point
          k  = km + 1;
          kp = km = k0 = k;

          // Reinitialize paths
          vmin = y[k * incy];
          vmax = vmin + 2 * lambda;
          umin = +lambda;
          umax = -lambda;
        }
        else if ( y[( k + 1 ) * incy] + umax > vmax + lambda )
        {
          // Condition (b2): vmax is too low
          for ( Integer j = k0; j <= kp; ++j ) x[j * incx] = vmax;

          // Restart from next point
          k  = kp + 1;
          kp = km = k0 = k;

          // Reinitialize paths
          vmax = y[k * incy];
          vmin = vmax - 2 * lambda;
          umin = +lambda;
          umax = -lambda;
        }
        else
        {
          // Condition (b3): Continue current segment
          ++k;
          Real yk = y[k * incy];
          umin += yk - vmin;
          umax += yk - vmax;

          // Update vmin if umin exceeds +λ
          if ( umin >= lambda )
          {
            vmin += ( umin - lambda ) / ( k - k0 + 1 );
            umin = lambda;
            km   = k;
          }

          // Update vmax if umax falls below -λ
          if ( umax <= -lambda )
          {
            vmax += ( umax + lambda ) / ( k - k0 + 1 );
            umax = -lambda;
            kp   = k;
          }
        }

        // Handle end of signal
        if ( k == nm1 )
        {
          if ( umin < 0 )
          {
            // vmin is too high, need negative jump
            for ( Integer j = k0; j <= km; ++j ) x[j * incx] = vmin;

            k  = km + 1;
            km = k0 = k;
            vmin    = y[k * incy];
            umin    = lambda;
            umax    = vmin + lambda - vmax;
          }
          else if ( umax > 0 )
          {
            // vmax is too low, need positive jump
            for ( Integer j = k0; j <= kp; ++j ) x[j * incx] = vmax;

            k  = kp + 1;
            kp = k0 = k;
            vmax    = y[k * incy];
            umax    = -lambda;
            umin    = vmax - lambda - vmin;
          }
          else
          {
            // Final segment
            Real tmp = vmin + umin / ( k - k0 + 1 );
            for ( Integer j = k0; j < N; ++j ) x[j * incx] = tmp;
            break;
          }
        }
      }
    }
  };

  // Common instantiations (kept for compatibility)
  using TVDf = TVD<float>;
  using TVDd = TVD<double>;

}  // namespace Utils

#endif

//
// eof: Utils_TVD.hh
//
