 /*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2023                                                      |
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

#include "Utils_TVD.hh"

#include <vector>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils {

  /*
  // minimize norm(y - x,2)^2 + lambda * sum(abs(diff( x )))
  //
  // Laurent Condat.
  // A Direct Algorithm for 1D Total Variation Denoising.
  // IEEE Signal Processing Letters,
  // Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057.
  // <10.1109/LSP.2013.2278339>. <hal-00675043v4>
  // https://hal.science/hal-00675043v4
  */
  template <typename Real>
  void
  TVD<Real>::denoise(
    Integer    N,
    Real const y[],
    Integer    incy,
    Real       lambda,
    Real       x[],
    Integer    incx
  ) {
    for ( Integer j = 0; j < N; ++j ) x[j*incx] = 0;
    Real vmin = y[0]-lambda;
    Real vmax = y[0]+lambda;
    Real umin = +lambda;
    Real umax = -lambda;
    Integer k  = 0;
    Integer k0 = 0;
    Integer km = 0;
    Integer kp = 0;
    Integer nm1 = N-1;
    while ( true ) {
      if ( k >= nm1 ) { x[nm1*incx] = vmin+umin; break; }
      if ( y[(k+1)*incy]+umin < vmin-lambda ) {
        // ----------- (b1)
        for ( Integer j = k0; j <= km; ++j ) x[j*incx] = vmin;
        // -----------
        k = km+1; kp = km = k0 = k;
        // -----------
        vmin = y[k*incy];
        vmax = vmin+2*lambda;
        umin = +lambda;
        umax = -lambda;
      } else if ( y[(k+1)*incy]+umax > vmax+lambda ) {
        // ----------- (b2)
        for ( Integer j = k0; j <= kp; ++j ) x[j*incx] = vmax;
        // -----------
        k = kp+1; kp = km = k0 = k;
        // -----------
        vmax = y[k*incy];
        vmin = vmax-2*lambda;
        umin = +lambda;
        umax = -lambda;
      } else {
        // ----------- (b3)
        ++k;
        Real yk = y[k*incy];
        umin += yk - vmin;
        umax += yk - vmax;
        if ( umin >= lambda ) {
          vmin += (umin-lambda)/(k-k0+1);
          umin = lambda;
          km   = k;
        }
        if ( umax <= -lambda ) {
          vmax += (umax+lambda)/(k-k0+1);
          umax = -lambda;
          kp   = k;
        }
      }
      if ( k == nm1 ) {
        if ( umin < 0 ) {
          // vmin is too high and a negative jump is necessary:
          for ( Integer j = k0; j <= km; ++j ) x[j*incx] = vmin;
          k = km+1; km = k0 = k;
          vmin = y[k*incy];
          umin = lambda;
          umax = vmin+lambda-vmax;
        } else if ( umax > 0 ) {
          // vmax is too low and a positive jump isnecessary:
          for ( Integer j = k0; j <= kp; ++j ) x[j*incx] = vmax;
          k = kp+1; kp = k0 = k;
          vmax = y[k*incy];
          umax = -lambda;
          umin = vmax-lambda-vmin;
        } else {
          Real tmp = vmin+umin/(k-k0+1);
          for ( Integer j = k0; j < N; ++j ) x[j*incx] = tmp;
          break;
        }
      }
    }
  }

  template class TVD<float>;
  template class TVD<double>;
}
