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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <cmath>

#define MEX_ERROR_MESSAGE \
"%=============================================================================%\n" \
"%  buildClothoid:  Compute parameters of the G2 Hermite clothoid fitting      %\n" \
"%                                                                             %\n" \
"%  USAGE: [ S0,S1,SM,SG,iter ] = buildClothoid3arcG2( x0, y0, th0, k0,        %\n" \
"%                                                     x1, y1, th1, k1 ) ;     %\n" \
"%                                                                             %\n" \
"%  USAGE: [ S0,S1,SM,SG,iter ] = buildClothoid3arcG2( x0, y0, th0, k0,        %\n" \
"%                                                     x1, y1, th1, k1,        %\n" \
"%                                                     dmax0, dmax1 ) ;        %\n" \
"%                                                                             %\n" \
"%  On input:                                                                  %\n" \
"%                                                                             %\n" \
"%       x0, y0  = coodinate of initial point                                  %\n" \
"%       theta0  = orientation (angle) of the clothoid at initial point        %\n" \
"%       k0      = initial curvature                                           %\n" \
"%       x1, y1  = coodinate of final point                                    %\n" \
"%       theta1  = orientation (angle) of the clothoid at final point          %\n" \
"%                                                                             %\n" \
"%       dmax0   = rough desidered maximum angle divergence from guess         %\n" \
"%       dmax1   = rough desidered maximum angle variation                     %\n" \
"%                                                                             %\n" \
"%  On output:                                                                 %\n" \
"%                                                                             %\n" \
"%       S0     = initial arc of clothoid                                      %\n" \
"%       SM     = middle arc of clothoid                                       %\n" \
"%       S1     = final arc of clothoid                                        %\n" \
"%       SG     = G1 guess clothoid                                            %\n" \
"%                                                                             %\n" \
"%  Optional Output                                                            %\n" \
"%                                                                             %\n" \
"%       iter   = >0 number of iteration used for the computation              %\n" \
"%                -1 computation failed                                        %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n" \
"%                                                                             %\n" \
"%  Autors: Enrico Bertolazzi and Marco Frego                                  %\n" \
"%          Department of Industrial Engineering                               %\n" \
"%          University of Trento                                               %\n" \
"%          enrico.bertolazzi@unitn.it                                         %\n" \
"%          marco.frego@unitn.it                                               %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n"

#define ASSERT(COND,MSG)                          \
  if ( !(COND) ) {                                \
    std::ostringstream ost ;                      \
    ost << "buildClothoid3arc: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;             \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_kappa0 prhs[3]
#define arg_x1     prhs[4]
#define arg_y1     prhs[5]
#define arg_theta1 prhs[6]
#define arg_kappa1 prhs[7]
#define arg_dmax0  prhs[8]
#define arg_dmax1  prhs[9]

#define out_S0     plhs[0]
#define out_S1     plhs[1]
#define out_SM     plhs[2]
#define out_SG     plhs[3]
#define out_iter   plhs[4]

static
void
save_struct( Clothoid::ClothoidCurve const & curve, mxArray * & plhs ) {
  char const * fieldnames[] = { "x0", "y0", "theta0", "k0", "dk", "L" } ;
  plhs = mxCreateStructMatrix(1,1,6,fieldnames);
  mxSetFieldByNumber( plhs, 0, 0, mxCreateDoubleScalar(curve.getX0()) );
  mxSetFieldByNumber( plhs, 0, 1, mxCreateDoubleScalar(curve.getY0()) );
  mxSetFieldByNumber( plhs, 0, 2, mxCreateDoubleScalar(curve.getTheta0()) );
  mxSetFieldByNumber( plhs, 0, 3, mxCreateDoubleScalar(curve.getKappa()) );
  mxSetFieldByNumber( plhs, 0, 4, mxCreateDoubleScalar(curve.getKappa_D()) );
  mxSetFieldByNumber( plhs, 0, 5, mxCreateDoubleScalar(curve.getSmax()) );
}

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  if ( nrhs == 0 && nlhs == 0 ) {
    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  }

  static Clothoid::G2solve3arc g2solve3arc ;

  Clothoid::valueType x0, y0, th0, k0, x1, y1, th1, k1, Dmax, dmax ;

  try {

    // Check for proper number of arguments, etc
    if ( nrhs != 10 && nrhs != 8 ) {
	    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    } else {
      for ( int kk = 0 ; kk < nrhs ; ++kk )
        if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS ||
             mxGetM(prhs[kk]) != 1 ||
             mxGetN(prhs[kk]) != 1  )
  	      mexErrMsgTxt("Input arguments must be real scalars");
    }

    ASSERT( nlhs >= 4 && nlhs <= 5,
            "wrong number of output arguments\n"
            "expected 4 or 5, found " << nlhs ) ;

    int iter ;
    if ( nrhs == 10 ) {
      x0     = mxGetScalar(arg_x0) ;
      y0     = mxGetScalar(arg_y0) ;
      th0    = mxGetScalar(arg_theta0) ;
      k0     = mxGetScalar(arg_kappa0) ;
      x1     = mxGetScalar(arg_x1) ;
      y1     = mxGetScalar(arg_y1) ;
      th1    = mxGetScalar(arg_theta1) ;
      k1     = mxGetScalar(arg_kappa1) ;
      Dmax   = mxGetScalar(arg_dmax0) ;
      dmax   = mxGetScalar(arg_dmax1) ;

      iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1, Dmax, dmax ) ;

      save_struct( g2solve3arc.getS0(),    out_S0 ) ;
      save_struct( g2solve3arc.getS1(),    out_S1 ) ;
      save_struct( g2solve3arc.getSM(),    out_SM ) ;
      save_struct( g2solve3arc.getGuess(), out_SG ) ;

      if ( nlhs > 4 ) {
        out_iter = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        *static_cast<int *>(mxGetData(out_iter)) = iter ;
      }
    } else {
      x0  = mxGetScalar(arg_x0) ;
      y0  = mxGetScalar(arg_y0) ;
      th0 = mxGetScalar(arg_theta0) ;
      k0  = mxGetScalar(arg_kappa0) ;
      x1  = mxGetScalar(arg_x1) ;
      y1  = mxGetScalar(arg_y1) ;
      th1 = mxGetScalar(arg_theta1) ;
      k1  = mxGetScalar(arg_kappa1) ;

      iter = g2solve3arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;

      save_struct( g2solve3arc.getS0(),    out_S0 ) ;
      save_struct( g2solve3arc.getS1(),    out_S1 ) ;
      save_struct( g2solve3arc.getSM(),    out_SM ) ;
      save_struct( g2solve3arc.getGuess(), out_SG ) ;

      if ( nlhs > 4 ) {
        out_iter = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        *static_cast<int *>(mxGetData(out_iter)) = iter ;
      }
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;
  } catch (...) {
  	mexErrMsgTxt("buildClothoid2arcG2 failed\n") ;
  }
}
