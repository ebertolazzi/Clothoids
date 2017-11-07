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

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  intersectClothoid:  Compute intersections betweed clothoids         %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    s1, s2 = intersectClothoid( clot1, clot2 ) ;                      %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%    clot1 = structure with the definition of the first clothoid       %\n" \
"%    clot2 = structure with the definition of the second clothoid      %\n" \
"%                                                                      %\n" \
"%    The structure must contain the field                              %\n" \
"%      x0     = initial x coordinate                                   %\n" \
"%      y0     = initial y coordinate                                   %\n" \
"%      theta0 = orientation (angle) of the clothoid at initial point   %\n" \
"%      k0     = curvature at initial point                             %\n" \
"%      dk     = derivative of curvature respect to arclength           %\n" \
"%                                                                      %\n" \
"%    A final field                                                     %\n" \
"%      L      = the lenght of the clothoid curve                       %\n" \
"%                                                                      %\n" \
"%    or in alternative                                                 %\n" \
"%      smin   = initial curvilinear coordinate of the curve            %\n" \
"%      smax   = final curvilinear coordinate of the curve              %\n" \
"%                                                                      %\n" \
"%    NOTE: (x0,y0,theta0,kappa0) are at curviliear coodinates 0        %\n" \
"%          smax  > smin  and smin may be negative                      %\n" \
"%          specify L is equivalent to pass smin = 0 and smax = L       %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  s1 = curvilinear coordinates of intersections on clot1              %\n" \
"%  s2 = curvilinear coordinates of intersections on clot2              %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

#define arg_clot1 prhs[0]
#define arg_clot2 prhs[1]

#define arg_S1    plhs[0]
#define arg_S2    plhs[1]

#define ASSERT(COND,MSG)                          \
  if ( !(COND) ) {                                \
    std::ostringstream ost ;                      \
    ost << "intersectClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;             \
  }

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {

    // Check for proper number of arguments, etc
    Clothoid::ClothoidCurve cc[2] ;
    Clothoid::valueType     x0, y0, theta0, k, dk, smin, smax ;

    if ( nrhs != 2 ) {
  	  mexErrMsgTxt("expexted 2 input arguments\n") ;
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    for ( Clothoid::indexType ii = 0 ; ii < 2 ; ++ii ) {
      if ( mxGetClassID(prhs[ii]) == mxSTRUCT_CLASS ) {
        mxArray * mxX0     = mxGetField( prhs[ii], 0, "x0"     );
        mxArray * mxY0     = mxGetField( prhs[ii], 0, "y0"     );
        mxArray * mxTheta0 = mxGetField( prhs[ii], 0, "theta0" );
        mxArray * mxK      = mxGetField( prhs[ii], 0, "k0"     );
        mxArray * mxDK     = mxGetField( prhs[ii], 0, "dk"     );
        mxArray * mxL      = mxGetField( prhs[ii], 0, "L"      );
        mxArray * mxSmin   = mxGetField( prhs[ii], 0, "smin"   );
        mxArray * mxSmax   = mxGetField( prhs[ii], 0, "smax"   );
        if ( mxX0 == nullptr ) {
    	    mexErrMsgTxt("missing field `x`\n") ;
          mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
          return ;
        } else {
          x0 = mxGetScalar(mxX0) ;
        }
        if ( mxY0 == nullptr ) {
    	    mexErrMsgTxt("missing field `y`\n") ;
          mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
          return ;
        } else {
          y0 = mxGetScalar(mxY0) ;
        }
        if ( mxTheta0 == nullptr ) {
    	    mexErrMsgTxt("missing field `theta`\n") ;
          mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
          return ;
        } else {
          theta0 = mxGetScalar(mxTheta0) ;
        }
        if ( mxK == nullptr ) {
    	    mexErrMsgTxt("missing field `k`\n") ;
          mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
          return ;
        } else {
          k = mxGetScalar(mxK) ;
        }
        if ( mxDK == nullptr ) {
    	    mexErrMsgTxt("missing field `dk`\n") ;
          mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
          return ;
        } else {
          dk = mxGetScalar(mxDK) ;
        }
        if ( mxL == nullptr ) {
          if ( mxSmin == nullptr || mxSmax == nullptr ) {
      	    mexErrMsgTxt("missing field `L` or `smin` and `smax`\n") ;
            mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
            return ;
          } else {
            smax = mxGetScalar(mxSmax) ;
            smin = mxGetScalar(mxSmin) ;
          }
        } else {
          smax = mxGetScalar(mxL) ;
          smin = 0 ;
        }
        cc[ii].build( x0, y0, theta0, k, dk, smin, smax ) ;
      } else {
  	    mexErrMsgTxt("Arguments expected to be STRUCT\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      }
    }

    std::vector<Clothoid::valueType> s1, s2 ;
    Clothoid::indexType max_iter  = 10 ;
    Clothoid::valueType tolerance = 1e-8 ;
 
    try {
      cc[0].intersect( cc[1], s1, s2, max_iter, tolerance ) ;
    }
    catch (...) {
    	mexErrMsgTxt("Intersection failed\n") ;
    }

    if ( nlhs > 0 ) {
      arg_S1 = mxCreateNumericMatrix(s1.size(),1, mxDOUBLE_CLASS, mxREAL);
      double * pS1 = mxGetPr(arg_S1) ;
      for ( unsigned i = 0 ; i < s1.size() ; ++i ) *pS1++ = s1[i] ;
    }
    if ( nlhs > 1 ) {
      arg_S2 = mxCreateNumericMatrix(s2.size(),1, mxDOUBLE_CLASS, mxREAL);
      double * pS2 = mxGetPr(arg_S2) ;
      for ( unsigned i = 0 ; i < s2.size() ; ++i ) *pS2++ = s2[i] ;
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("intersectClothoid failed\n") ;
  }

}
