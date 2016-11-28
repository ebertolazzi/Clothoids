/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex.h"

#include <vector>
#include <sstream>
#include <stdexcept>

#define ASSERT(COND,MSG)                   \
  if ( !(COND) ) {                         \
    std::ostringstream ost ;               \
    ost << "targetClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;      \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_L      prhs[5]

#define out_TV     plhs[0]
#define out_curv2  plhs[1]
#define out_jerk2  plhs[2]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  bbClothoid:  Compute total variation, curvature and jerk            %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    [tv,curv2,jerk2] = targetClothoid( S ) ;                          %\n" \
"%    [tv,curv2,jerk2] = targetClothoid( x0, y0, theta0, k, dk, L ) ;   %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%    x0, y0 = coodinate of initial point                               %\n" \
"%    theta0 = orientation (angle) of the clothoid at initial point     %\n" \
"%    k      = curvature at initial point                               %\n" \
"%    dk     = derivative of curvature respect to arclength             %\n" \
"%    L      = the lenght of the clothoid curve or a vector of length   %\n" \
"%             where to compute the clothoid values                     %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
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

  static Clothoid::ClothoidCurve clot ;

  try {

    Clothoid::valueType x0, y0, theta0, k, dk, L ;

    // Check for proper number of arguments, etc
    if ( nrhs == 6 ) {
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
      k      = mxGetScalar(arg_k) ;
      dk     = mxGetScalar(arg_dk) ;
      L      = mxGetScalar(arg_L) ;

    } else if ( nrhs == 1 ) {
      ASSERT( mxIsStruct(prhs[0]),
              "First argument must be a struct" ) ;

      mxArray * mx_x0     = mxGetField(prhs[0],0,"x") ;
      mxArray * mx_y0     = mxGetField(prhs[0],0,"y") ;
      mxArray * mx_theta0 = mxGetField(prhs[0],0,"theta") ;
      mxArray * mx_k      = mxGetField(prhs[0],0,"k") ;
      mxArray * mx_dk     = mxGetField(prhs[0],0,"dk") ;
      mxArray * mx_L      = mxGetField(prhs[0],0,"L") ;

      ASSERT( mx_x0     != nullptr, "Field `x` is missing" );
      ASSERT( mx_y0     != nullptr, "Field `y` is missing" );
      ASSERT( mx_theta0 != nullptr, "Field `theta` is missing" );
      ASSERT( mx_k      != nullptr, "Field `k` is missing" );
      ASSERT( mx_dk     != nullptr, "Field `dk` is missing" );
      ASSERT( mx_L      != nullptr, "Field `L` is missing" );

      ASSERT( mxGetClassID(mx_x0) == mxDOUBLE_CLASS && !mxIsComplex(mx_x0),
              "Field `x` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_y0) == mxDOUBLE_CLASS && !mxIsComplex(mx_y0),
              "Field `y` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_theta0) == mxDOUBLE_CLASS && !mxIsComplex(mx_theta0),
              "Field `theta` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_k) == mxDOUBLE_CLASS && !mxIsComplex(mx_k),
              "Field `k` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_dk) == mxDOUBLE_CLASS && !mxIsComplex(mx_dk),
              "Field `dk` must be a real double scalar" );
      ASSERT( mxGetClassID(mx_L) == mxDOUBLE_CLASS && !mxIsComplex(mx_L),
              "Field `L` must be a real double scalar" );

      x0     = mxGetScalar(mx_x0) ;
      y0     = mxGetScalar(mx_y0) ;
      theta0 = mxGetScalar(mx_theta0) ;
      k      = mxGetScalar(mx_k) ;
      dk     = mxGetScalar(mx_dk) ;
      L      = mxGetScalar(mx_L) ;

    } else {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }
    
    ASSERT( nlhs == 3,
            "Expected 3 output arguments, found: " << nlhs ) ;

    // costruisco bb
    Clothoid::ClothoidCurve clot( x0, y0, theta0, k, dk, L ) ;

    out_TV    = mxCreateDoubleScalar( clot.thetaTotalVariation() ) ;
    out_curv2 = mxCreateDoubleScalar( clot.integralCurvature2() ) ;
    out_jerk2 = mxCreateDoubleScalar( clot.integralJerk2() ) ;

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("targetClothoid failed\n") ;
  }

}
