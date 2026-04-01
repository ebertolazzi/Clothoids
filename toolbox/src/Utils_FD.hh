/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003-2022                                                 |
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
// file: Utils_FD.hh
//

#pragma once

#ifndef UTILS_FD_dot_HH
#define UTILS_FD_dot_HH

#include "Utils.hh"

namespace Utils
{

  /**
   * @brief First derivative using a 2-point non-uniform stencil.
   *
   * Approximates the first derivative f'(x0) using two consecutive points:
   * (x0, y0), (x1, y1). This is a simple forward/backward difference.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   *
   * @return Approximation of f'(x0)
   */
  template <typename Real> inline Real first_derivative_2p( Real x0, Real y0, Real x1, Real y1 )
  {
    return ( y1 - y0 ) / ( x1 - x0 );
  }

  /**
   * @brief First derivative using a 3-point non-uniform stencil.
   *
   * Uses Newton's divided differences for numerical stability.
   * Works for both increasing and decreasing grids.
   * Accurate for polynomials up to degree 2.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   *
   * @return Approximation of f'(x0)
   */
  template <typename Real> inline Real first_derivative_3p( Real x0, Real y0, Real x1, Real y1, Real x2, Real y2 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h02 = x2 - x0;  // = h01 + h12

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;

    // Second divided difference
    const Real f012 = ( f12 - f01 ) / h02;

    // First derivative at x0 using Newton form: f01 + f012 * (x0 - x1)
    return f01 - f012 * h01;
  }

  /**
   * @brief First derivative using a 4-point non-uniform stencil.
   *
   * Uses Newton's divided differences for numerical stability.
   * Works for both increasing and decreasing grids.
   * Accurate for polynomials up to degree 3.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   * @param x3  x-coordinate of fourth point
   * @param y3  y-coordinate of fourth point
   *
   * @return Approximation of f'(x0)
   */
  template <typename Real>
  inline Real first_derivative_4p( Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x3, Real y3 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h23 = x3 - x2;
    const Real h02 = x2 - x0;  // = h01 + h12
    const Real h03 = x3 - x0;  // = h01 + h12 + h23
    const Real h13 = x3 - x1;  // = h12 + h23

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;
    const Real f23 = ( y3 - y2 ) / h23;

    // Second divided differences
    const Real f012 = ( f12 - f01 ) / h02;
    const Real f123 = ( f23 - f12 ) / h13;

    // Third divided difference
    const Real f0123 = ( f123 - f012 ) / h03;

    // First derivative at x0: f01 + f012*(x0-x1) + f0123*(x0-x1)*(x0-x2)
    return f01 - f012 * h01 + f0123 * h01 * h02;
  }

  /**
   * @brief First derivative using a 5-point non-uniform stencil.
   *
   * Uses Newton's divided differences for numerical stability.
   * Works for both increasing and decreasing grids.
   * Accurate for polynomials up to degree 4.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   * @param x3  x-coordinate of fourth point
   * @param y3  y-coordinate of fourth point
   * @param x4  x-coordinate of fifth point
   * @param y4  y-coordinate of fifth point
   *
   * @return Approximation of f'(x0)
   */
  template <typename Real> inline Real first_derivative_5p(
    Real x0,
    Real y0,
    Real x1,
    Real y1,
    Real x2,
    Real y2,
    Real x3,
    Real y3,
    Real x4,
    Real y4 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h23 = x3 - x2;
    const Real h34 = x4 - x3;
    const Real h02 = x2 - x0;  // = h01 + h12
    const Real h03 = x3 - x0;  // = h01 + h12 + h23
    const Real h04 = x4 - x0;  // = h01 + h12 + h23 + h34
    const Real h13 = x3 - x1;  // = h12 + h23
    const Real h14 = x4 - x1;  // = h12 + h23 + h34
    const Real h24 = x4 - x2;  // = h23 + h34

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;
    const Real f23 = ( y3 - y2 ) / h23;
    const Real f34 = ( y4 - y3 ) / h34;

    // Second divided differences
    const Real f012 = ( f12 - f01 ) / h02;
    const Real f123 = ( f23 - f12 ) / h13;
    const Real f234 = ( f34 - f23 ) / h24;

    // Third divided differences
    const Real f0123 = ( f123 - f012 ) / h03;
    const Real f1234 = ( f234 - f123 ) / h14;

    // Fourth divided difference
    const Real f01234 = ( f1234 - f0123 ) / h04;

    // First derivative at x0: f01 + f012*(x0-x1) + f0123*(x0-x1)*(x0-x2) + f01234*(x0-x1)*(x0-x2)*(x0-x3)
    return f01 - f012 * h01 + f0123 * h01 * h02 - f01234 * h01 * h02 * h03;
  }

  /**
   * @brief Second derivative using a 3-point non-uniform stencil via Newton's Divided Differences.
   *
   * Mathematical derivation:
   * P(x) = f[x0] + f[x0,x1](x-x0) + f[x0,x1,x2](x-x0)(x-x1)
   * P''(x) = 2 * f[x0,x1,x2]
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   *
   * @return Approximation of f''(x0)
   */
  template <typename Real> inline Real second_derivative_3p( Real x0, Real y0, Real x1, Real y1, Real x2, Real y2 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h02 = x2 - x0;  // = h01 + h12

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;

    // Second divided difference
    const Real f012 = ( f12 - f01 ) / h02;

    // Second derivative at x0: 2 * f012
    return 2 * f012;
  }

  /**
   * @brief Second derivative using a 4-point non-uniform stencil.
   *
   * Uses Newton's divided differences for numerical stability.
   * Works for both increasing and decreasing grids.
   * Accurate for polynomials up to degree 3.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   * @param x3  x-coordinate of fourth point
   * @param y3  y-coordinate of fourth point
   *
   * @return Approximation of f''(x0)
   */

  template <typename Real>
  inline Real second_derivative_4p( Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x3, Real y3 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h23 = x3 - x2;
    const Real h02 = x2 - x0;
    const Real h03 = x3 - x0;
    const Real h13 = x3 - x1;

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;
    const Real f23 = ( y3 - y2 ) / h23;

    // Second divided differences
    const Real f012 = ( f12 - f01 ) / h02;
    const Real f123 = ( f23 - f12 ) / h13;

    // Third divided difference
    const Real f0123 = ( f123 - f012 ) / h03;

    // Correct formula: P''(x0) = 2*f012 + 2*f0123*(x0-x1 + x0-x2)
    // Since x0-x1 = -h01 and x0-x2 = -h02
    return 2 * ( f012 + f0123 * ( -h01 - h02 ) );
  }

  /**
   * @brief Second derivative using a 5-point non-uniform stencil.
   *
   * Uses Newton's divided differences for numerical stability.
   * Works for both increasing and decreasing grids.
   * Accurate for polynomials up to degree 4.
   *
   * @param x0  x-coordinate of first point
   * @param y0  y-coordinate of first point
   * @param x1  x-coordinate of second point
   * @param y1  y-coordinate of second point
   * @param x2  x-coordinate of third point
   * @param y2  y-coordinate of third point
   * @param x3  x-coordinate of fourth point
   * @param y3  y-coordinate of fourth point
   * @param x4  x-coordinate of fifth point
   * @param y4  y-coordinate of fifth point
   *
   * @return Approximation of f''(x0)
   */

  template <typename Real> inline Real second_derivative_5p(
    Real x0,
    Real y0,
    Real x1,
    Real y1,
    Real x2,
    Real y2,
    Real x3,
    Real y3,
    Real x4,
    Real y4 )
  {
    const Real h01 = x1 - x0;
    const Real h12 = x2 - x1;
    const Real h23 = x3 - x2;
    const Real h34 = x4 - x3;
    const Real h02 = x2 - x0;  // = h01 + h12
    const Real h03 = x3 - x0;  // = h01 + h12 + h23
    const Real h04 = x4 - x0;  // = h01 + h12 + h23 + h34
    const Real h13 = x3 - x1;  // = h12 + h23
    const Real h14 = x4 - x1;  // = h12 + h23 + h34
    const Real h24 = x4 - x2;  // = h23 + h34

    // First divided differences
    const Real f01 = ( y1 - y0 ) / h01;
    const Real f12 = ( y2 - y1 ) / h12;
    const Real f23 = ( y3 - y2 ) / h23;
    const Real f34 = ( y4 - y3 ) / h34;

    // Second divided differences
    const Real f012 = ( f12 - f01 ) / h02;
    const Real f123 = ( f23 - f12 ) / h13;
    const Real f234 = ( f34 - f23 ) / h24;

    // Third divided differences
    const Real f0123 = ( f123 - f012 ) / h03;
    const Real f1234 = ( f234 - f123 ) / h14;

    // Fourth divided difference
    const Real f01234 = ( f1234 - f0123 ) / h04;

    // Second derivative at x0 using Newton polynomial of degree 4:
    // P''(x0) = 2*f012 + 2*f0123*(x0-x1 + x0-x2)
    //          + 2*f01234*[(x0-x1)*(x0-x2) + (x0-x1)*(x0-x3) + (x0-x2)*(x0-x3)]
    // Since x0-x1 = -h01, x0-x2 = -h02, x0-x3 = -h03
    return 2 * ( f012 + f0123 * ( -h01 - h02 ) + f01234 * ( h01 * h02 + h01 * h03 + h02 * h03 ) );
  }

}  // namespace Utils

#endif

//
// eof: Utils_FD.hh
//
