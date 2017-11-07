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
"%======================================================================%\n" \
"% FresnelCS:  Compute Fresnel sine and cosine integrals                %\n" \
"%                                                                      %\n" \
"% USAGE: [FresnelC,FresnelS] = FresnelCS(y) ;                          %\n" \
"%                                                                      %\n" \
"%  Fresnel integral are defined as:                                    %\n" \
"%                                                                      %\n" \
"%  C(y) = int_0^y cos( (pi/2)*t^2 ) dt                                 %\n" \
"%  S(y) = int_0^y sin( (pi/2)*t^2 ) dt                                 %\n" \
"%                                                                      %\n" \
"%  The algorithm is described in:                                      %\n" \
"%    Atlas for computing mathematical functions: an illustrated guide  %\n" \
"%    for practitioners, with programs in C and Mathematica.            %\n" \
"%    William J. Thompson New York : Wiley, c1997.                      %\n" \
"%                                                                      %\n" \
"%  The code is a sligly modification and translation to C language     %\n" \
"%  of original code developed by                                       %\n" \
"%  Venkata Sivakanth Telasula (sivakanth.telasula@gmail.com)           %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%   y = argument of the function C(y) and S(y), it may be a vector     %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  FresnelC = The value(s) of C(y)                                     %\n" \
"%  FresnelS = The value(s) of S(y)                                     %\n" \
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
    if ( nrhs != 1 ) {
	    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    } else if ( mxGetClassID(prhs[0]) != mxDOUBLE_CLASS ) {
	    mexErrMsgTxt("Input argument should be double");
      return ;
    } else if ( mxIsComplex(prhs[0]) ) {
	    mexErrMsgTxt("Input argument should be real");
      return ;
    }

    // Output array
    mwSize         nDimNum = mxGetNumberOfDimensions(prhs[0]);
    mwSize const * pDims   = mxGetDimensions(prhs[0]) ;
    double * pC, *pS ;

    if ( nlhs < 2 ) {
  	  plhs[0] = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxCOMPLEX);
  	  pC      = mxGetPr(plhs[0]);
  	  pS      = mxGetPi(plhs[0]);
    } else {
	    plhs[0] = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
	    plhs[1] = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
	    pC      = mxGetPr(plhs[0]);
	    pS      = mxGetPr(plhs[1]);
    }

    int nElemNum = mxGetNumberOfElements(prhs[0]);
    double * pX = mxGetPr(prhs[0]) ;
    for ( int i = 0 ; i < nElemNum ; ++i ) Clothoid::FresnelCS( *pX++, *pC++, *pS++ ) ;

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("FresnelCS failed\n") ;
  }

}
