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

#define ASSERT(COND,MSG)              \
  if ( !(COND) ) {                    \
    std::ostringstream ost ;          \
    ost << "biarc: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ; \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_x1     prhs[3]
#define arg_y1     prhs[4]
#define arg_theta1 prhs[5]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  biarc:  Compute biarc fitting.                                      %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    [arc1,arc2] = biarc( x0, y0, theta0, x1, y1, theta1 ) ;           %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%  x0, y0 = coordinate of initial point                                %\n" \
"%  theta0 = orientation (angle) at the initial point                   %\n" \
"%  x1, y1 = coordinate of final point                                  %\n" \
"%  theta1 = orientation (angle) at the final point                     %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  arc1, arc2 = rational B-spline of the two arc                       %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

#include <iostream> 

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {
    Clothoid::valueType x0, y0, theta0, x1, y1, theta1, xs, ys, thetas, L0, L1 ;
    Clothoid::valueType knots[7], Poly[4][3] ;

    // Check for proper number of arguments, etc
    bool ok = nrhs == 6;
    if ( !ok ) mexErrMsgTxt("expected 6 input arguments") ;
    ok = nlhs == 2;
    if ( !ok ) mexErrMsgTxt("expected 2 outputs") ;
    if ( ok ) {
      for ( int kk = 0 ; kk < nrhs ; ++kk ) {
        ASSERT( mxGetClassID(prhs[kk]) == mxDOUBLE_CLASS &&
                !mxIsComplex(prhs[kk]),
                "Argument N." << kk+1 <<
                " must be a real double scalar" );
        ASSERT( mxGetM(prhs[kk]) == 1 && mxGetN(prhs[kk]) == 1,
                "Argument N." << kk+1 << " must be a scalar" );
      }

      x0     = mxGetScalar(arg_x0) ;
      y0     = mxGetScalar(arg_y0) ;
      theta0 = mxGetScalar(arg_theta0) ;
      x1     = mxGetScalar(arg_x1) ;
      y1     = mxGetScalar(arg_y1) ;
      theta1 = mxGetScalar(arg_theta1) ;
    }
    if ( ok ) ok = Clothoid::Biarc( x0, y0, theta0, x1, y1, theta1, xs, ys, thetas, L0, L1 ) ;
    if ( ok ) {
      mwSize npt = Clothoid::ArcToNURBS( x0, y0, theta0, xs, ys, knots, Poly ) ;

      mxArray * K = mxCreateDoubleMatrix(1, npt+3, mxREAL);
      mxArray * P = mxCreateDoubleMatrix(3, npt, mxREAL);

      double * pP = mxGetPr(P) ;
      for ( mwSize k = 0 ; k < npt ; ++k ) {
        *pP++ = Poly[k][0] ;
        *pP++ = Poly[k][1] ;
        *pP++ = Poly[k][2] ;
      }
      std::copy( knots, knots+npt+3, mxGetPr(K) ) ;

      char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" } ;

      plhs[0] = mxCreateStructMatrix(1,1,6,fieldnames);
      mxSetFieldByNumber( plhs[0], 0, 0, mxCreateString("rB") );
      mxSetFieldByNumber( plhs[0], 0, 1, mxCreateDoubleScalar(3.0) );
      mxSetFieldByNumber( plhs[0], 0, 2, mxCreateDoubleScalar(2.0) );
      mxSetFieldByNumber( plhs[0], 0, 3, mxCreateDoubleScalar(npt) );
      mxSetFieldByNumber( plhs[0], 0, 4, K );
      mxSetFieldByNumber( plhs[0], 0, 5, P );

      npt = Clothoid::ArcToNURBS( xs, ys, thetas, x1, y1, knots, Poly ) ;
      K   = mxCreateDoubleMatrix(1, npt+3, mxREAL);
      P   = mxCreateDoubleMatrix(3, npt, mxREAL);
      pP = mxGetPr(P) ;
      for ( mwSize k = 0 ; k < npt ; ++k ) {
        *pP++ = Poly[k][0] ;
        *pP++ = Poly[k][1] ;
        *pP++ = Poly[k][2] ;
      }
      std::copy( knots, knots+npt+3, mxGetPr(K) ) ;

      plhs[1] = mxCreateStructMatrix(1,1,6,fieldnames);
      mxSetFieldByNumber( plhs[1], 0, 0, mxCreateString("rB") );
      mxSetFieldByNumber( plhs[1], 0, 1, mxCreateDoubleScalar(3.0) );
      mxSetFieldByNumber( plhs[1], 0, 2, mxCreateDoubleScalar(2.0) );
      mxSetFieldByNumber( plhs[1], 0, 3, mxCreateDoubleScalar(npt) );
      mxSetFieldByNumber( plhs[1], 0, 4, K );
      mxSetFieldByNumber( plhs[1], 0, 5, P );
    }
    if ( !ok ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("biarc failed\n") ;
  }

}
