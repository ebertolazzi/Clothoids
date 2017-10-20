/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex.h"

#include <vector>
#include <sstream>
#include <stdexcept>

#define ASSERT(COND,MSG)                      \
  if ( !(COND) ) {                            \
    std::ostringstream ost ;                  \
    ost << "TriTriOverlap: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;         \
  }

#define arg_p0 prhs[0]
#define arg_p1 prhs[1]
#define arg_p2 prhs[2]
#define arg_q0 prhs[3]
#define arg_q1 prhs[4]
#define arg_q2 prhs[5]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  TriTriOverlap:  Check if two trinagles overlap                      %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    intersect = TriTriOverlap( p0, p1, p2, q0, q1, q2 ) ;             %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%    p0, p1, p2 = coodinates of first triangle                         %\n" \
"%    q0, q1, q2 = coodinates of second triangle                        %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%    intersect = true if the triangle overlap                          %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {

    // Check for proper number of arguments, etc
    if ( nrhs != 6 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    for ( int kk = 0 ; kk < nrhs ; ++kk ) {
      mwSize nDimNum = mxGetNumberOfDimensions(prhs[kk]);
      if ( nDimNum != 2 || mxGetM(prhs[kk])*mxGetN(prhs[kk]) != 2 ) {
        mexErrMsgTxt("Input arguments must vectors of size 2\n");
        return ;
      }
      if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS ) {
        mexErrMsgTxt("Input arguments must be REAL vector\n");
        return ;
      }
    }

    G2lib::valueType const * p0 = mxGetPr(arg_p0) ;
    G2lib::valueType const * p1 = mxGetPr(arg_p1) ;
    G2lib::valueType const * p2 = mxGetPr(arg_p2) ;

    G2lib::valueType const * q0 = mxGetPr(arg_q0) ;
    G2lib::valueType const * q1 = mxGetPr(arg_q1) ;
    G2lib::valueType const * q2 = mxGetPr(arg_q2) ;

    // costruisco bb
    Triangle2D::Triangle2D<G2lib::valueType> T1( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] ) ;
    Triangle2D::Triangle2D<G2lib::valueType> T2( q0[0], q0[1], q1[0], q1[1], q2[0], q2[1] ) ;

    if ( nlhs == 1 ) {
      bool over = T1.overlap(T2) ;
      plhs[0] = mxCreateDoubleScalar( over ? 1.0 : 0.0 ) ;
    } else {
      mexErrMsgTxt("There must be an output argument\n");
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("TriTriOverlap failed\n") ;
  }
}
