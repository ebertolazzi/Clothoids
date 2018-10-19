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

#ifndef CUBIC_ROOTS_FLOCKE_HH
#define CUBIC_ROOTS_FLOCKE_HH

#include "G2lib.hh"

namespace G2lib {

  using std::abs;
  using std::pow;

  static real_type const third   = 1./3.;
  static real_type const one27th = 1./27.;
  static real_type const two27th = 2./27.;

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/
  void
  solveQuadratic( real_type   a,
                  real_type   b,
                  real_type   c,
                  real_type & r1,
                  real_type & r2,
                  int_type  & nr,
                  int_type  & nc );

  /*\
  ... Calculate the zeros of the cubic A*z^3 + B*z^2 + C*z + D.
  ...
  ... N. FLOCKE, Flash Center for Computational Science, University of Chicago
  ... Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  ... for Physical Applications
  ... ACM Transactions on Mathematical Software, Vol. 41, No. 4, 2015.
  ... DOI: http://dx.doi.org/10.1145/2699468
  \*/

  int_type
  solveCubic( real_type   A,
              real_type   B,
              real_type   C,
              real_type   D,
              real_type & r1,
              real_type & r2,
              real_type & r3,
              int_type  & nr,
              int_type  & nc );

}

#endif

// EOF: CubicRootsFlocke.hh
