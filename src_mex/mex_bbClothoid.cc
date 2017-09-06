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

#define ASSERT(COND,MSG)                   \
  if ( !(COND) ) {                         \
    std::ostringstream ost ;               \
    ost << "bbClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;      \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_L      prhs[5]
#define arg_angle  prhs[6]
#define arg_size   prhs[7]
#define arg_offs   prhs[8]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  bbClothoid:  Compute a series of bounding triangles for             %\n" \
"%               a clothoid curve.                                      %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    TT = bbClothoid( S, angle, size ) ;                               %\n" \
"%    TT = bbClothoid( S, angle, size, offs ) ;                         %\n" \
"%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size ) ;        %\n" \
"%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size, offs ) ;  %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%  x0, y0 = coodinate of initial point                                 %\n" \
"%  theta0 = orientation (angle) of the clothoid at initial point       %\n" \
"%  k      = curvature at initial point                                 %\n" \
"%  dk     = derivative of curvature respect to arclength               %\n" \
"%  L      = the lenght of the clothoid curve or a vector of length     %\n" \
"%           where to compute the clothoid values                       %\n" \
"%  angle  = maximum variation of angle in the bounding box triangle    %\n" \
"%  size   = maximum height of the bounding box triangle                %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  TT = matrix 6 x N whose column are the coordinates of bounding      %\n" \
"%       triangles.                                                     %\n" \
"%       TT(:,i) = [ x0, y0, x1, y1, x2, y2 ].'                         %\n" \
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

    Clothoid::valueType x0, y0, theta0, k, dk, L, max_angle, max_size, offs ;

    // Check for proper number of arguments, etc
    if ( nrhs > 4 ) {
      ASSERT( nrhs == 8 ||  nrhs == 9,
              "expected 8 or 9 arguments, found " << nrhs ) ;
      for ( int kk = 0 ; kk < nrhs ; ++kk ) {
        ASSERT( mxGetClassID(prhs[kk]) == mxDOUBLE_CLASS &&
                !mxIsComplex(prhs[kk]),
                "Argument N." << kk+1 <<
                " must be a real double scalar" );
        ASSERT( mxGetM(prhs[kk]) == 1 && mxGetN(prhs[kk]) == 1,
                "Argument N." << kk+1 << " must be a scalar" );
      }

      x0        = mxGetScalar(arg_x0) ;
      y0        = mxGetScalar(arg_y0) ;
      theta0    = mxGetScalar(arg_theta0) ;
      k         = mxGetScalar(arg_k) ;
      dk        = mxGetScalar(arg_dk) ;
      L         = mxGetScalar(arg_L) ;
      max_angle = mxGetScalar(arg_angle) ;
      max_size  = mxGetScalar(arg_size) ;
      offs      = 0 ;
      if ( nrhs > 8 ) offs = mxGetScalar(arg_offs) ;

    } else if ( nrhs == 3 || nrhs == 4 ) {
      ASSERT( mxIsStruct(prhs[0]),
              "First argument must be a struct" ) ;

      ASSERT( mxGetClassID(prhs[1]) == mxDOUBLE_CLASS && !mxIsComplex(prhs[1]),
              "Argument N.2 must be a real double matrix" );

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

      max_angle = mxGetScalar(prhs[1]) ;
      max_size  = mxGetScalar(prhs[2]) ;
      offs      = 0 ;
      if ( nrhs == 4 ) offs = mxGetScalar(prhs[3]) ;

    } else {
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }
    
    // costruisco bb
    Clothoid::ClothoidCurve clot( x0, y0, theta0, k, dk, L ) ;
    std::vector<Clothoid::ClothoidCurve> c ;
    std::vector<Clothoid::Triangle2D<Clothoid::valueType> > t ;

    clot.bbSplit( max_angle, max_size, offs, c, t ) ;

    if ( nlhs == 1 ) {
      plhs[0] = mxCreateDoubleMatrix(6, t.size(), mxREAL);
      double * pT = mxGetPr(plhs[0]);
      for ( int i = 0 ; i < t.size() ; ++i ) {
        *pT++ = t[i].x1() ; *pT++ = t[i].y1() ;
        *pT++ = t[i].x2() ; *pT++ = t[i].y2() ;
        *pT++ = t[i].x3() ; *pT++ = t[i].y3() ;
      }
    } else {
      mexErrMsgTxt("There must be an output argument\n");
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("buildClothoid failed\n") ;
  }

}
