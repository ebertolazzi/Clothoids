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

#define ASSERT(COND,MSG)                         \
  if ( !(COND) ) {                               \
    std::ostringstream ost ;                     \
    ost << "pointsOnClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;            \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k0     prhs[3]
#define arg_dk     prhs[4]
#define arg_ss     prhs[5]
#define arg_offs   prhs[6]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  pointsOnClothoid:  Compute points on a clothoid curve.              %\n" \
"%                     Used for plotting purpose.                       %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    XY    = pointsOnClothoid( x0, y0, theta0, k0, dk, ss ) ;          %\n" \
"%    [X,Y] = pointsOnClothoid( x0, y0, theta0, k0, dk, ss ) ;          %\n" \
"%    XY    = pointsOnClothoid( x0, y0, theta0, k0, dk, ss, offs) ;     %\n" \
"%    [X,Y] = pointsOnClothoid( x0, y0, theta0, k0, dk, ss, offs) ;     %\n" \
"%    XY    = pointsOnClothoid( clot, ss ) ;                            %\n" \
"%    [X,Y] = pointsOnClothoid( clot, ss ) ;                            %\n" \
"%    XY    = pointsOnClothoid( clot, ss, offs ) ;                      %\n" \
"%    [X,Y] = pointsOnClothoid( clot, ss, offs ) ;                      %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%    x0, y0  = coodinate of initial point                              %\n" \
"%    theta0  = orientation (angle) of the clothoid at initial point    %\n" \
"%    k0      = curvature at initial point                              %\n" \
"%    dk      = derivative of curvature respect to arclength            %\n" \
"%    ss      = a vector with the curvilinear coordinates where         %\n" \
"%              to compute the clothoid values                          %\n" \
"%    offs    = curve offset                                            %\n" \
"%                                                                      %\n" \
"%    clot    = struct with field `x0`, `y0`, `theta0`, `k0`, `dk`      %\n" \
"%              and a final field                                       %\n" \
"%                `L` the lenght of the clothoid curve                  %\n" \
"%              or in alternative two field                             %\n" \
"%                 `smin` = initial curvilinear coordinate of the curve %\n" \
"%                 `smax` = final curvilinear coordinate of the curve   %\n" \
"%                                                                      %\n" \
"%  On output: (1 argument)                                             %\n" \
"%                                                                      %\n" \
"%    XY = matrix 2 x NPTS whose column are the points of the clothoid  %\n" \
"%                                                                      %\n" \
"%  On output: (2 argument)                                             %\n" \
"%                                                                      %\n" \
"%    X  = matrix 1 x NPTS X coordinate of points of the clothoid       %\n" \
"%    Y  = matrix 1 x NPTS Y coordinate of points of the clothoid       %\n" \
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

    Clothoid::valueType x0, y0, theta0, k0, dk, offs ;
    int                 npts ;
    double              * pSS ;

    // Check for proper number of arguments, etc
    if ( nrhs == 6 || nrhs == 7 ) {
      for ( int kk = 0 ; kk < nrhs ; ++kk ) {
        ASSERT( mxGetClassID(prhs[kk]) == mxDOUBLE_CLASS &&
                !mxIsComplex(prhs[kk]),
                "Argument N." << kk+1 <<
                " must be a real double scalar" );
        if ( kk == 5 ) continue ;
        ASSERT( mxGetM(prhs[kk]) == 1 && mxGetN(prhs[kk]) == 1,
                "Argument N." << kk+1 << " must be a scalar" );
      }

      x0     = mxGetScalar(arg_x0) ;
      y0     = mxGetScalar(arg_y0) ;
      theta0 = mxGetScalar(arg_theta0) ;
      k0     = mxGetScalar(arg_k0) ;
      dk     = mxGetScalar(arg_dk) ;
      offs   = 0 ;
      if ( nrhs == 7 ) offs = mxGetScalar(arg_offs) ;

      ASSERT( mxGetClassID(arg_ss) == mxDOUBLE_CLASS && !mxIsComplex(arg_ss),
              "Argument N.6 must be a real double matrix" );

      pSS    = mxGetPr(arg_ss) ;
      npts   = int(mxGetNumberOfElements(arg_ss)) ;

    } else if ( nrhs == 2 || nrhs == 3 ) {
      ASSERT( mxIsStruct(prhs[0]),
              "First argument must be a struct" ) ;

      ASSERT( mxGetClassID(prhs[1]) == mxDOUBLE_CLASS && !mxIsComplex(prhs[1]),
              "Argument N.2 must be a real double matrix" );

      mxArray * mx_x0     = mxGetField(prhs[0],0,"x0") ;
      mxArray * mx_y0     = mxGetField(prhs[0],0,"y0") ;
      mxArray * mx_theta0 = mxGetField(prhs[0],0,"theta0") ;
      mxArray * mx_k0     = mxGetField(prhs[0],0,"k0") ;
      mxArray * mx_dk     = mxGetField(prhs[0],0,"dk") ;

      ASSERT( mx_x0     != nullptr, "Field `x0` is missing" );
      ASSERT( mx_y0     != nullptr, "Field `y0` is missing" );
      ASSERT( mx_theta0 != nullptr, "Field `theta0` is missing" );
      ASSERT( mx_k0     != nullptr, "Field `k0` is missing" );
      ASSERT( mx_dk     != nullptr, "Field `dk` is missing" );

      ASSERT( mxGetClassID(mx_x0) == mxDOUBLE_CLASS && !mxIsComplex(mx_x0),
              "Field `x0` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_y0) == mxDOUBLE_CLASS && !mxIsComplex(mx_y0),
              "Field `y0` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_theta0) == mxDOUBLE_CLASS && !mxIsComplex(mx_theta0),
              "Field `theta0` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_k0) == mxDOUBLE_CLASS && !mxIsComplex(mx_k0),
              "Field `k0` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_dk) == mxDOUBLE_CLASS && !mxIsComplex(mx_dk),
              "Field `dk` must be a real double scalar" );

      x0     = mxGetScalar(mx_x0) ;
      y0     = mxGetScalar(mx_y0) ;
      theta0 = mxGetScalar(mx_theta0) ;
      k0     = mxGetScalar(mx_k0) ;
      dk     = mxGetScalar(mx_dk) ;

      offs   = 0 ;
      if ( nrhs == 3 ) offs = mxGetScalar(prhs[2]) ;

      pSS  = mxGetPr(prhs[1]) ;
      npts = int(mxGetNumberOfElements(prhs[1])) ;

    } else {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }
    
    Clothoid::ClothoidCurve curve ;
    curve.build( x0, y0, theta0, k0, dk, 1 ) ;

    if ( nlhs == 1 ) {
  	  plhs[0] = mxCreateDoubleMatrix(2, npts, mxREAL);
  	  double * pXY = mxGetPr(plhs[0]);
      for ( int i = 0 ; i < npts ; ++i ) {
        curve.eval( *pSS, offs, pXY[0], pXY[1] ) ;
        pXY += 2 ; ++pSS ;
      }
    } else if ( nlhs == 2 ) {
      plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1, npts, mxREAL);
      double * pX = mxGetPr(plhs[0]);
      double * pY = mxGetPr(plhs[1]);
      for ( int i = 0 ; i < npts ; ++i ) {
        curve.eval( *pSS, offs, *pX, *pY ) ;
        ++pX ; ++pY ; ++pSS ;
      }
    } else {
      mexErrMsgTxt("Output argument must be 1 or 2");
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("pointsOnClothoid failed\n") ;
  }

}
