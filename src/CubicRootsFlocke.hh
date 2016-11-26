/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

namespace PolynomialRoots {

  typedef double valueType ;
  typedef int    indexType ;

  using std::abs ;
  using std::pow ;

  static valueType const machepsi = std::numeric_limits<valueType>::epsilon() ;
  static valueType const third    = 1./3. ;
  static valueType const one27th  = 1./27. ;
  static valueType const two27th  = 2./27. ;

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/
  void
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2,
                  indexType & nr,
                  indexType & nc ) ;

  /*\
  ... Calculate the zeros of the cubic A*z^3 + B*z^2 + C*z + D.
  ...
  ... N. FLOCKE, Flash Center for Computational Science, University of Chicago
  ... Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  ... for Physical Applications
  ... ACM Transactions on Mathematical Software, Vol. 41, No. 4, 2015.
  ... DOI: http://dx.doi.org/10.1145/2699468
  \*/

  indexType
  solveCubic( valueType   A,
              valueType   B,
              valueType   C,
              valueType   D,
              valueType & r1,
              valueType & r2,
              valueType & r3,
              indexType & nr,
              indexType & nc ) ;

}

// EOF: CubicRootsFlocke.hh
