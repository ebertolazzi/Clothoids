/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
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

#include <sstream>
#include <stdexcept>

#define MEX_ERROR_MESSAGE \
"%===================================================================%\n" \
"%  Compute Fresnel sine and cosine integrals momenta                %\n" \
"%                                                                   %\n" \
"%  USAGE:                                                           %\n" \
"%    [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;                  %\n" \
"%                                                                   %\n" \
"%  Integrals are defined as:                                        %\n" \
"%                                                                   %\n" \
"%    X_k(a,b,c) = int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt       %\n" \
"%    Y_k(a,b,c) = int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt       %\n" \
"%                                                                   %\n" \
"%  On input:                                                        %\n" \
"%                                                                   %\n" \
"%    nk      = number of momentae to be computed (1 <= nk <= 3)     %\n" \
"%    a, b, c = the parameters of the integrals or vectors of the    %\n" \
"%              same dimensions length(a)=length(b)=length(c)        %\n" \
"%                                                                   %\n" \
"%  On output:                                                       %\n" \
"%                                                                   %\n" \
"%    X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]  %\n" \
"%    Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]  %\n" \
"%                                                                   %\n" \
"%    the dimension of the vectors X are nk x length(a)              %\n" \
"%                                                                   %\n" \
"%===================================================================%\n" \
"%                                                                   %\n" \
"%  Autor: Enrico Bertolazzi                                         %\n" \
"%         Department of Industrial Engineering                      %\n" \
"%         University of Trento                                      %\n" \
"%         enrico.bertolazzi@unitn.it                                %\n" \
"%                                                                   %\n" \
"%===================================================================%\n"


#define arg_nk prhs[0]
#define arg_a  prhs[1]
#define arg_b  prhs[2]
#define arg_c  prhs[3]

#define arg_C  plhs[0]
#define arg_S  plhs[1]

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {

    // Check for proper number of arguments, etc
    if ( nrhs != 4 ) {
	    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    } else if ( mxGetClassID(arg_a) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_b) != mxDOUBLE_CLASS ||
                mxGetClassID(arg_c) != mxDOUBLE_CLASS ) {
	    mexErrMsgTxt("Input arguments should be double");
      return ;
    } else if ( mxIsComplex(arg_nk) ||
                mxIsComplex(arg_a)  ||
                mxIsComplex(arg_b)  ||
                mxIsComplex(arg_c) ) {
	    mexErrMsgTxt("Input arguments should be real");
      return ;
    }

    if ( mxGetM(arg_nk) != 1 || mxGetN(arg_nk) != 1 )
      mexErrMsgTxt("First argument must be a scalar");

    // Output array
    int nk = int(mxGetScalar(arg_nk)) ;
    if ( nk < 1 || nk > 3 )
      mexErrMsgTxt("First argument must be an integer in the range [1:3]");
  
    mwSize na = mxGetNumberOfElements(arg_a) ;
    mwSize nb = mxGetNumberOfElements(arg_b) ;
    mwSize nc = mxGetNumberOfElements(arg_c) ;
  
    if ( na != nb || na != nc )
      mexErrMsgTxt("Second to last arguments must be vectors of the same length");
  
    double const * pa = mxGetPr(arg_a);
    double const * pb = mxGetPr(arg_b);
    double const * pc = mxGetPr(arg_c);

    arg_C = mxCreateDoubleMatrix(nk,na,mxREAL) ;
    arg_S = mxCreateDoubleMatrix(nk,na,mxREAL) ;
    double * pC = mxGetPr(arg_C) ;
    double * pS = mxGetPr(arg_S) ;
  
    for ( mwSize k = 0 ; k < na ; ++k ) {
      Clothoid::GeneralizedFresnelCS( nk, *pa, *pb, *pc, pC, pS ) ;
      ++pa ; ++pb ; ++pc ; pC += nk ; pS += nk ;
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("GeneralizedFresnelCS failed\n") ;
  }
}
