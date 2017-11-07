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

#include <sstream>
#include <stdexcept>

#define MEX_ERROR_MESSAGE \
"%=============================================================================%\n" \
"%  buildClothoid:  Compute parameters of the G1 Hermite clothoid fitting      %\n" \
"%                                                                             %\n" \
"%  USAGE:                                                                     %\n" \
"%    [S,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;             %\n" \
"%    [S,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1, derivative ) ; %\n" \
"%                                                                             %\n" \
"%  On input:                                                                  %\n" \
"%                                                                             %\n" \
"%    x0, y0     = coodinate of initial point                                  %\n" \
"%    theta0     = orientation (angle) of the clothoid at initial point        %\n" \
"%    x1, y1     = coodinate of final point                                    %\n" \
"%    theta1     = orientation (angle) of the clothoid at final point          %\n" \
"%    derivative = if present and true compute additional                      %\n" \
"%                 partial derivative of the solution                          %\n" \
"%                                                                             %\n" \
"%  On output:                                                                 %\n" \
"%                                                                             %\n" \
"%    S.x     = x-coodinate of initial point                                   %\n" \
"%    S.y     = y-coodinate of initial point                                   %\n" \
"%    S.theta = orientation (angle) of the clothoid at initial point           %\n" \
"%    S.k     = curvature at initial point                                     %\n" \
"%    S.dk    = derivative of curvature respect to arclength,                  %\n" \
"%              notice that curvature at final point is k+dk*L                 %\n" \
"%    S.L     = the lenght of the clothoid curve from initial to final point   %\n" \
"%    S.iter  = Newton Iterations used to solve the interpolation problem      %\n" \
"%                                                                             %\n" \
"%  optional output                                                            %\n" \
"%                                                                             %\n" \
"%    S.k_1   = partial derivative of the solution respect to theta0           %\n" \
"%    S.dk_1  = partial derivative of the solution respect to theta0           %\n" \
"%    S.L_1   = partial derivative of the solution respect to theta0           %\n" \
"%                                                                             %\n" \
"%    S.k_2   = partial derivative of the solution respect to theta1           %\n" \
"%    S.dk_2  = partial derivative of the solution respect to theta1           %\n" \
"%    S.L_2   = partial derivative of the solution respect to theta1           %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n" \
"%                                                                             %\n" \
"%  Autors: Enrico Bertolazzi and Marco Frego                                  %\n" \
"%          Department of Industrial Engineering                               %\n" \
"%          University of Trento                                               %\n" \
"%          enrico.bertolazzi@unitn.it                                         %\n" \
"%          m.fregox@gmail.com                                                 %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n"

#define ASSERT(COND,MSG)                      \
  if ( !(COND) ) {                            \
    std::ostringstream ost ;                  \
    ost << "buildClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;         \
  }

#define arg_x0         prhs[0]
#define arg_y0         prhs[1]
#define arg_theta0     prhs[2]
#define arg_x1         prhs[3]
#define arg_y1         prhs[4]
#define arg_theta1     prhs[5]
#define arg_derivative prhs[6]

#define arg_res        plhs[0]
#define arg_iter       plhs[1]

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {

    // Check for proper number of arguments, etc
    if ( nrhs < 6 ) {
	    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    } else if ( mxGetClassID(arg_x0)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_y0)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_theta0) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_x1)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_y1)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_theta1) != mxDOUBLE_CLASS ) {
	    mexErrMsgTxt("Input arguments should be double");
    } else if ( mxIsComplex(arg_x0)     ||
                mxIsComplex(arg_y0)     ||
                mxIsComplex(arg_theta0) ||
                mxIsComplex(arg_x1)     ||
                mxIsComplex(arg_y1)     ||
                mxIsComplex(arg_theta1) ) {
	    mexErrMsgTxt("Input arguments should be real");
    }

    for ( int kk = 0 ; kk < 6 ; ++kk )
      if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	      mexErrMsgTxt("Input arguments must be scalars");

    ASSERT( nlhs == 1 || nlhs == 2,
            "wrong number of output arguments\n"
            "expected 1 or 2 found " << nlhs ) ;
  
    bool deriv = false ;
    if ( nrhs > 6 ) deriv = mxIsLogicalScalarTrue( arg_derivative ) ;
  
    Clothoid::valueType x0     = mxGetScalar(arg_x0) ;
    Clothoid::valueType y0     = mxGetScalar(arg_y0) ;
    Clothoid::valueType theta0 = mxGetScalar(arg_theta0) ;
    Clothoid::valueType x1     = mxGetScalar(arg_x1) ;
    Clothoid::valueType y1     = mxGetScalar(arg_y1) ;
    Clothoid::valueType theta1 = mxGetScalar(arg_theta1) ;

    int iter ;
    Clothoid::valueType k, dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ;
    if ( deriv ) {
      iter = Clothoid::buildClothoid( x0, y0, theta0, x1, y1, theta1,
                                      k, dk, L,
                                      k_1, dk_1, L_1,
                                      k_2, dk_2, L_2 ) ;
      static char const * fieldnames[] = {
        "x0", "y0", "theta0", "k0", "dk", "L",
        "k0_1", "dk_1", "L_1", "k0_2", "dk_2", "L_2"
      } ;
      arg_res = mxCreateStructMatrix(1,1,12,fieldnames);
    } else {
      iter = Clothoid::buildClothoid( x0, y0, theta0, x1, y1, theta1, k, dk, L ) ;
      static char const * fieldnames[] = { "x0", "y0", "theta0", "k0", "dk", "L" } ;
      arg_res = mxCreateStructMatrix(1,1,6,fieldnames);
    }

    mxSetFieldByNumber( arg_res, 0, 0, mxCreateDoubleScalar(x0) );
    mxSetFieldByNumber( arg_res, 0, 1, mxCreateDoubleScalar(y0) );
    mxSetFieldByNumber( arg_res, 0, 2, mxCreateDoubleScalar(theta0) );
    mxSetFieldByNumber( arg_res, 0, 3, mxCreateDoubleScalar(k) );
    mxSetFieldByNumber( arg_res, 0, 4, mxCreateDoubleScalar(dk) );
    mxSetFieldByNumber( arg_res, 0, 5, mxCreateDoubleScalar(L) );

    if ( deriv ) {
      mxSetFieldByNumber( arg_res, 0, 6, mxCreateDoubleScalar(k_1) );
      mxSetFieldByNumber( arg_res, 0, 7, mxCreateDoubleScalar(dk_1) );
      mxSetFieldByNumber( arg_res, 0, 8, mxCreateDoubleScalar(L_1) );

      mxSetFieldByNumber( arg_res, 0, 9,  mxCreateDoubleScalar(k_2) );
      mxSetFieldByNumber( arg_res, 0, 10, mxCreateDoubleScalar(dk_2) );
      mxSetFieldByNumber( arg_res, 0, 11, mxCreateDoubleScalar(L_2) );
    }

    if ( nlhs > 1 ) {
      arg_iter = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL) ;
      *((int*)mxGetData(arg_iter)) = iter ;
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("buildClothoid failed\n") ;
  }
}
