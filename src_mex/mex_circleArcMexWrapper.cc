/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Circle.hh"
#include "mex.h"
#include "mex_class_handle.hh"

#include <vector>
#include <sstream>
#include <stdexcept>

#define ASSERT(COND,MSG) \
  if ( !(COND) ) { \
    std::ostringstream ost ; \
    ost << "circleArcMexWrapper: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ; \
  }

#define arg_in_0 prhs[0]
#define arg_in_1 prhs[1]
#define arg_in_2 prhs[2]
#define arg_in_3 prhs[3]
#define arg_in_4 prhs[4]

#define arg_out_0 plhs[0]
#define arg_out_1 plhs[1]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  circle:  Compute cicle arc                                          %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    OBJ = circleArcMexWrapper( 'new', p0, theta0, kur, arclen ) ;     %\n" \
"%    OBJ = circleArcMexWrapper( 'new_by_G1', p0, theta0, p1 ) ;        %\n" \
"%    OBJ = circleArcMexWrapper( 'new_by_3p', p0, p1, p2 ) ;            %\n" \
"%                                                                      %\n" \
"%    circleArcMexWrapper( 'setup', OBJ, p0, theta0, L ) ;              %\n" \
"%    circleArcMexWrapper( 'setup_by_G1', OBJ, p0, theta0, p1  ) ;      %\n" \
"%    circleArcMexWrapper( 'setup_by_3p', OBJ, p0, p1, p2 ) ;           %\n" \
"%    circleArcMexWrapper( 'origin', OBJ, x0, y0 ) ;                    %\n" \
"%    circleArcMexWrapper( 'translate', OBJ, tx, ty ) ;                 %\n" \
"%    circleArcMexWrapper( 'trim', OBJ, smin, smax ) ;                  %\n" \
"%                                                                      %\n" \
"%    circleArcMexWrapper( 'delete', OBJ ) ;                            %\n" \
"%                                                                      %\n" \
"%    P           = circleArcMexWrapper( 'eval', OBJ, s ) ;             %\n" \
"%    [P,theta]   = circleArcMexWrapper( 'eval', OBJ, s ) ;             %\n" \
"%                                                                      %\n" \
"%    dP          = circleArcMexWrapper( 'eval_D', OBJ, s ) ;           %\n" \
"%    [dP,dtheta] = circleArcMexWrapper( 'eval_D', OBJ, s ) ;           %\n" \
"%                                                                      %\n" \
"%    dP          = circleArcMexWrapper( 'eval_DD', OBJ, s ) ;          %\n" \
"%    [dP,dtheta] = circleArcMexWrapper( 'eval_DD', OBJ, s ) ;          %\n" \
"%                                                                      %\n" \
"%    dP          = circleArcMexWrapper( 'eval_DDD', OBJ, s ) ;         %\n" \
"%    [dP,dtheta] = circleArcMexWrapper( 'eval_DDD', OBJ, s ) ;         %\n" \
"%                                                                      %\n" \
"%    C           = circleArcMexWrapper( 'center', OBJ ) ;              %\n" \
"%    nurbs       = circleArcMexWrapper( 'to_nurbs', OBJ ) ;            %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

static
double
getScalarValue( mxArray const * arg, char const msg[] ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  ASSERT( dims[0] == 1 && dims[1] == 1,
          msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  return mxGetScalar(arg) ;
}

static
double const *
getArrayValues( mxArray const * arg, mwSize & size, char const msg[] ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  ASSERT( dims[0] == 1 || dims[1] == 1,
          msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  size = mxGetNumberOfElements(arg);
  return mxGetPr(arg);
}

static
void
setScalarValue( mxArray * & arg, double const & value ) {
  arg = mxCreateDoubleMatrix(1, 1, mxREAL);
  double * pA = mxGetPr(arg);
  *pA = value;
}


extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // the first argument must be a string
  if ( nrhs == 0 ) {
    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  }

  try {

    ASSERT( mxIsChar(arg_in_0), "First argument must be a string" ) ;
    std::string cmd = mxArrayToString(arg_in_0) ;
    mwSize size0, size1, size2 ;


    Circle::CircleArc * ptr ;

    if ( cmd != "new_by_G1" && cmd != "new_by_3p" && cmd != "new_by_3p" )
      ptr = convertMat2Ptr<Circle::CircleArc>(arg_in_1);

    if ( cmd == "new" ) {
      Circle::valueType smin, smax;
      ASSERT(nrhs == 5, "expected 5 inputs");
      ASSERT(nlhs == 1, "expected 1 output");

      Circle::valueType const *   p0 = getArrayValues( arg_in_1, size0, "in new, p0 expected to be a real vector");
      Circle::valueType       theta0 = getScalarValue( arg_in_2, "in new, theta0 expected to be a real scalar");
      Circle::valueType          kur = getScalarValue( arg_in_3, "in new, ray expected to be a real scalar");
      Circle::valueType const * alen = getArrayValues( arg_in_4, size1, "in new, p1 expected to be a real vector");

      ASSERT( size0 == 2, "point p0 expected of size 2, found size0 = " << size0 );

      if ( size1 == 1 ) { smin = 0 ; smax = alen[0] ; }
      else              { smin = alen[0] ; smax = alen[1] ; }

      ptr = new Circle::CircleArc( p0[0], p0[1], theta0, kur, smin, smax );
      plhs[0] = convertPtr2Mat<Circle::CircleArc>(ptr);

    } else if ( cmd == "new_by_G1" ) {

      ASSERT(nrhs == 4, "expected 4 inputs");
      ASSERT(nlhs == 1, "expected 1 output");

      Circle::valueType const * p0 = getArrayValues( arg_in_1, size0, "in new_by_G1, p0 expected to be a real vector");
      Circle::valueType     theta0 = getScalarValue( arg_in_2, "in new_by_G1, theta0 expected to be a real scalar");
      Circle::valueType const * p1 = getArrayValues( arg_in_3, size1, "in new_by_G1, p1 expected to be a real vector");

      ASSERT( size0 == 1 && size1 == 2,
              "point p0 and p1 expected of size 2, found size0 = " << size0 <<
              " size1 = " << size1 );

      ptr = new Circle::CircleArc();
      ptr->build_G1( p0[0], p0[1], theta0, p1[0], p1[1] ) ;
      plhs[0] = convertPtr2Mat<Circle::CircleArc>(ptr);

    } else if ( cmd == "new_by_3p" ) {

      ASSERT(nrhs == 4, "expected 4 inputs");
      ASSERT(nlhs == 1, "expected 1 output");

      Circle::valueType const * p0 = getArrayValues( arg_in_1, size0, "in new_by_3p, p0 expected to be a real vector");
      Circle::valueType const * p1 = getArrayValues( arg_in_2, size1, "in new_by_3p, p1 expected to be a real vector");
      Circle::valueType const * p2 = getArrayValues( arg_in_3, size2, "in new_by_3p, p2 expected to be a real vector");

      ASSERT( size0 == 1 && size1 == 2 && size2 == 2,
              "point p0, p1 and p2 expected of size 2, found size0 = " << size0 <<
              " size1 = " << size1 << " size2 = " << size2 );

      ptr = new Circle::CircleArc();
      ptr->build_3P( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] ) ;
      plhs[0] = convertPtr2Mat<Circle::CircleArc>(ptr);

    } else if ( cmd == "delete" ) {

      ASSERT(nrhs == 2, "expected 2 inputs");
      ASSERT(nlhs == 0, "expected no output");
      // Destroy the C++ object
      destroyObject<Circle::CircleArc>(prhs[1]);

    } else if ( cmd == "setup" ) {
    } else if ( cmd == "setup_by_G1" ) {
    } else if ( cmd == "setup_by_3p" ) {
    } else if ( cmd == "eval" ) {
    } else if ( cmd == "eval_D" ) {
    } else if ( cmd == "eval_DD" ) {
    } else if ( cmd == "eval_DDD" ) {
    } else if ( cmd == "center" ) {
    } else if ( cmd == "origin" ) {
    } else if ( cmd == "translate" ) {
    } else if ( cmd == "trim" ) {
    } else if ( cmd == "to_nurbs" ) {

    } else {
      ASSERT(false, "Unknown command: " << cmd );
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("circleArc failed\n") ;
  }

}
