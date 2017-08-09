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
#define arg_thstar prhs[6]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  biarc:  Compute biarc fitting.                                      %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%   [arc1,arc2,ok] = biarc( x0, y0, theta0, x1, y1, theta1 ) ;         %\n" \
"%   [arc1,arc2,ok] = biarc( x0, y0, theta0, x1, y1, theta1, thstar ) ; %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%  x0, y0 = coordinate of initial point                                %\n" \
"%  theta0 = orientation (angle) at the initial point                   %\n" \
"%  x1, y1 = coordinate of final point                                  %\n" \
"%  theta1 = orientation (angle) at the final point                     %\n" \
"%  thstar = orientation (angle) at the intermediate (optional)         %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  arc1, arc2 = rational B-spline of the two arc                       %\n" \
"%  ok         = false if computation fails                             %\n" \
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
    Clothoid::BiarcData BiData ;
    Clothoid::valueType knots[12], Poly[9][3] ;

    // Check for proper number of arguments, etc
    bool ok = nrhs == 6 || nrhs == 7 ;
    if ( !ok ) mexErrMsgTxt("expected 6 or 7 input arguments") ;
    ok = nlhs == 3;
    if ( !ok ) mexErrMsgTxt("expected 3 outputs") ;
    if ( ok ) {
      for ( int kk = 0 ; kk < nrhs ; ++kk ) {
        ASSERT( mxGetClassID(prhs[kk]) == mxDOUBLE_CLASS &&
                !mxIsComplex(prhs[kk]),
                "Argument N." << kk+1 <<
                " must be a real double scalar" );
        ASSERT( mxGetM(prhs[kk]) == 1 && mxGetN(prhs[kk]) == 1,
                "Argument N." << kk+1 << " must be a scalar" );
      }

      BiData.x0     = mxGetScalar(arg_x0) ;
      BiData.y0     = mxGetScalar(arg_y0) ;
      BiData.theta0 = mxGetScalar(arg_theta0) ;
      BiData.x1     = mxGetScalar(arg_x1) ;
      BiData.y1     = mxGetScalar(arg_y1) ;
      BiData.theta1 = mxGetScalar(arg_theta1) ;
      if ( nrhs == 7 ) BiData.thetas = mxGetScalar(arg_thstar) ;
    }
    if ( ok ) ok = Clothoid::Biarc( BiData, nrhs == 6 ) ;
    if ( ok ) {
      //mwSize npt = Clothoid::ArcToNURBS( BiData.x0, BiData.y0, BiData.theta0,
      //                                   BiData.xs, BiData.ys, knots, Poly ) ;
      mwSize npt = Clothoid::ArcToNURBS( BiData.theta0,
                                         BiData.x0, BiData.y0, BiData.c0, BiData.s0,
                                         BiData.L0, BiData.kappa0,
                                         knots, Poly ) ;
      mxArray * K = mxCreateDoubleMatrix(1, npt+3, mxREAL);
      mxArray * P = mxCreateDoubleMatrix(3, npt, mxREAL);

      double * pP = mxGetPr(P) ;
      for ( mwSize k = 0 ; k < npt ; ++k ) {
        *pP++ = Poly[k][0] ;
        *pP++ = Poly[k][1] ;
        *pP++ = Poly[k][2] ;
      }
      std::copy( knots, knots+npt+3, mxGetPr(K) ) ;

      char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs", "length", "curvature" } ;

      plhs[0] = mxCreateStructMatrix(1,1,8,fieldnames);
      mxSetFieldByNumber( plhs[0], 0, 0, mxCreateString("rB") );
      mxSetFieldByNumber( plhs[0], 0, 1, mxCreateDoubleScalar(3.0) );
      mxSetFieldByNumber( plhs[0], 0, 2, mxCreateDoubleScalar(2.0) );
      mxSetFieldByNumber( plhs[0], 0, 3, mxCreateDoubleScalar(npt) );
      mxSetFieldByNumber( plhs[0], 0, 4, K );
      mxSetFieldByNumber( plhs[0], 0, 5, P );
      mxSetFieldByNumber( plhs[0], 0, 6, mxCreateDoubleScalar(BiData.L0) );
      mxSetFieldByNumber( plhs[0], 0, 7, mxCreateDoubleScalar(BiData.kappa0) );

      //npt = Clothoid::ArcToNURBS( BiData.xs, BiData.ys, BiData.thetas,
      //                            BiData.x1, BiData.y1, knots, Poly ) ;
      npt = Clothoid::ArcToNURBS( BiData.thetas,
                                  BiData.xs, BiData.ys, BiData.cs, BiData.ss,
                                  BiData.L1, BiData.kappa1,
                                  knots, Poly ) ;

      K   = mxCreateDoubleMatrix(1, npt+3, mxREAL);
      P   = mxCreateDoubleMatrix(3, npt, mxREAL);
      pP = mxGetPr(P) ;
      for ( mwSize k = 0 ; k < npt ; ++k ) {
        *pP++ = Poly[k][0] ;
        *pP++ = Poly[k][1] ;
        *pP++ = Poly[k][2] ;
      }
      std::copy( knots, knots+npt+3, mxGetPr(K) ) ;

      plhs[1] = mxCreateStructMatrix(1,1,8,fieldnames);
      mxSetFieldByNumber( plhs[1], 0, 0, mxCreateString("rB") );
      mxSetFieldByNumber( plhs[1], 0, 1, mxCreateDoubleScalar(3.0) );
      mxSetFieldByNumber( plhs[1], 0, 2, mxCreateDoubleScalar(2.0) );
      mxSetFieldByNumber( plhs[1], 0, 3, mxCreateDoubleScalar(npt) );
      mxSetFieldByNumber( plhs[1], 0, 4, K );
      mxSetFieldByNumber( plhs[1], 0, 5, P );
      mxSetFieldByNumber( plhs[1], 0, 6, mxCreateDoubleScalar(BiData.L1) );
      mxSetFieldByNumber( plhs[1], 0, 7, mxCreateDoubleScalar(BiData.kappa1) );
    } else {
      plhs[0] = mxCreateDoubleScalar(0.0);
      plhs[1] = mxCreateDoubleScalar(0.0);
    }
    plhs[2] = mxCreateLogicalScalar(ok);
  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("biarc failed\n") ;
  }

}
