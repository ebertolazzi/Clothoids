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
"%  ClothoidMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting %\n" \
"%                                                                             %\n" \
"%  USAGE:                                                                     %\n" \
"%    OBJ = ClothoidMexWrapper( 'new', x0, y0, theta0, x1, y1, theta1 ) ;      %\n" \
"%                                                                             %\n" \
"%  On input:                                                                  %\n" \
"%                                                                             %\n" \
"%    x0, y0     = coodinate of initial point                                  %\n" \
"%    theta0     = orientation (angle) of the clothoid at initial point        %\n" \
"%    x1, y1     = coodinate of final point                                    %\n" \
"%    theta1     = orientation (angle) of the clothoid at final point          %\n" \
"%                                                                             %\n" \
"%  On output:                                                                 %\n" \
"%                                                                             %\n" \
"%    OBJ     = pointer to the internal object                                 %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n" \
"%                                                                             %\n" \
"%  Autors: Enrico Bertolazzi and Marco Frego and Paolo Bevilacqua             %\n" \
"%          Department of Industrial Engineering                               %\n" \
"%          Department of XXXXXXXX                                             %\n" \
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

#define arg_in_0 prhs[0]
#define arg_in_1 prhs[1]
#define arg_in_2 prhs[2]
#define arg_in_3 prhs[3]
#define arg_in_4 prhs[4]
#define arg_in_5 prhs[5]
#define arg_in_6 prhs[6]

#define arg_out_0 plhs[0]
#define arg_out_1 plhs[1]

static
double
getScalarValue( mxArray const * arg, char const * msg ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  ASSERT( dims[0] == 1 && dims[1] == 1,
          msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  return mxGetScalar(arg) ;
}

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {


  // the first argument must be a string
  if ( nrhs == 0 ) {
    mexErrorMessage() ;
    return ;
  }

  try {

    ASSERT( mxIsChar(arg_in_0), "First argument must be a string" ) ;
    string cmd = mxArrayToString(arg_in_0) ;

    if ( cmd == "new" ) {

      ASSERT( nrhs == 7, "expected 7 inputs") ;

      Clothoid::valueType x0     = getScalarValue(arg_in_1,"Error in reading x0") ;
      Clothoid::valueType y0     = getScalarValue(arg_in_2,"Error in reading y0") ;
      Clothoid::valueType theta0 = getScalarValue(arg_in_3,"Error in reading theta0") ;
      Clothoid::valueType x1     = getScalarValue(arg_in_4,"Error in reading x1") ;
      Clothoid::valueType y1     = getScalarValue(arg_in_5,"Error in reading y1") ;
      Clothoid::valueType theta1 = getScalarValue(arg_in_6,"Error in reading theta1") ;

      ClothoidCurve * ptr = new ClothoidCurve();
      ptr->build_G1( x0, y0, theta0, x1, y1, theta1 ) ;
      plhs[0] = convertPtr2Mat<ClothoidCurve>(ptr);

    } else if ( cmd == "delete" ) {

      ASSERT( nrhs == 2, "expected 2 inputs") ;
      ASSERT( nlhs == 0, "expected no output") ;
      // Destroy the C++ object
      destroyObject<ClothoidCurve>(prhs[1]);
      // Warn if other commands were ignored
    } else {

    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("clothoid failed\n") ;
  }
}
