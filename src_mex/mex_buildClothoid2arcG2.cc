/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex_utils.hh"

#include <sstream>
#include <stdexcept>

#define MEX_ERROR_MESSAGE \
"%=============================================================================%\n" \
"%  buildClothoid:  Compute parameters of the G1 Hermite clothoid fitting      %\n" \
"%                                                                             %\n" \
"%  USAGE: [S0,S1,flg] = buildClothoid2arcG2( x0, y0, th0, k0,                 %\n" \
"%                                            x1, y1, th1, k1 ) ;              %\n" \
"%                                                                             %\n" \
"%  On input:                                                                  %\n" \
"%                                                                             %\n" \
"%       x0, y0  = coodinate of initial point                                  %\n" \
"%       theta0  = orientation (angle) of the clothoid at initial point        %\n" \
"%       k0      = initial curvature                                           %\n" \
"%       x1, y1  = coodinate of final point                                    %\n" \
"%       theta1  = orientation (angle) of the clothoid at final point          %\n" \
"%                                                                             %\n" \
"%  On output:                                                                 %\n" \
"%       S0     = initial arc of clothoid                                      %\n" \
"%       S1     = final arc of clothoid                                        %\n" \
"%       flg    = >0 number of iteration used for the computation              %\n" \
"%                -1 computation failed                                        %\n" \
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


namespace G2lib {

  using namespace std;

  #define arg_x0     prhs[0]
  #define arg_y0     prhs[1]
  #define arg_theta0 prhs[2]
  #define arg_kappa0 prhs[3]
  #define arg_x1     prhs[4]
  #define arg_y1     prhs[5]
  #define arg_theta1 prhs[6]
  #define arg_kappa1 prhs[7]

  #define out_S0     plhs[0]
  #define out_S1     plhs[1]
  #define out_flg    plhs[2]

  static
  void
  save_struct( ClothoidCurve const & curve, mxArray * & plhs ) {
    char const * fieldnames[] = { "x0", "y0", "theta0", "k0", "dk", "L" } ;
    plhs = mxCreateStructMatrix(1,1,6,fieldnames);
    mxSetFieldByNumber( plhs, 0, 0, mxCreateDoubleScalar(curve.xBegin()) );
    mxSetFieldByNumber( plhs, 0, 1, mxCreateDoubleScalar(curve.yBegin()) );
    mxSetFieldByNumber( plhs, 0, 2, mxCreateDoubleScalar(curve.thetaBegin()) );
    mxSetFieldByNumber( plhs, 0, 3, mxCreateDoubleScalar(curve.kappaBegin()) );
    mxSetFieldByNumber( plhs, 0, 4, mxCreateDoubleScalar(curve.kappa_D()) );
    mxSetFieldByNumber( plhs, 0, 5, mxCreateDoubleScalar(curve.length()) );
  }

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    if ( nrhs == 0 && nlhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    // Check for proper number of arguments, etc
    if ( nrhs != 8 ) {
	    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    } else if ( mxGetClassID(arg_x0)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_y0)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_theta0) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_kappa0) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_x1)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_y1)     != mxDOUBLE_CLASS ||
                mxGetClassID(arg_theta1) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_kappa1) != mxDOUBLE_CLASS ) {
	    mexErrMsgTxt("Input arguments should be double");
    }

    for ( int kk = 0 ; kk < 8 ; ++kk )
      if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	      mexErrMsgTxt("Input arguments must be scalars");

    MEX_ASSERT( nlhs == 3,
                "wrong number of output arguments\n"
                "expected 3, found " << nlhs ) ;

    G2solve2arc g2solve2arc ;

    valueType x0  = mxGetScalar(arg_x0) ;
    valueType y0  = mxGetScalar(arg_y0) ;
    valueType th0 = mxGetScalar(arg_theta0) ;
    valueType k0  = mxGetScalar(arg_kappa0) ;
    valueType x1  = mxGetScalar(arg_x1) ;
    valueType y1  = mxGetScalar(arg_y1) ;
    valueType th1 = mxGetScalar(arg_theta1) ;
    valueType k1  = mxGetScalar(arg_kappa1) ;

    g2solve2arc.build( x0, y0, th0, k0, x1, y1, th1, k1 ) ;
    int iter = g2solve2arc.solve() ;

    save_struct( g2solve2arc.getS0(), out_S0 ) ;
    save_struct( g2solve2arc.getS1(), out_S1 ) ;

    out_flg = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *static_cast<int *>(mxGetData(out_flg)) = iter ;

  }
}

