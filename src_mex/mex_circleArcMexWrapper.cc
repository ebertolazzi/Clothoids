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
#define arg_out_2 plhs[2]
#define arg_out_3 plhs[3]

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

    bool do_new = cmd == "new" ;

    if ( !do_new ) ptr = convertMat2Ptr<Circle::CircleArc>(arg_in_1);
    else           ptr = new Circle::CircleArc();

    if ( do_new || cmd == "build" ) {

      Circle::indexType kk = do_new ? 0 : 1 ;

      ASSERT( nlhs == 1, "expected 1 output" );

      if ( nrhs >= 6+kk && nrhs <= 7+kk ) {

        Circle::valueType x0     = getScalarValue( prhs[1+kk], "CircleArc: x0 expected to be a real scalar" );
        Circle::valueType y0     = getScalarValue( prhs[2+kk], "CircleArc: y0 expected to be a real scalar" );
        Circle::valueType theta0 = getScalarValue( prhs[3+kk], "CircleArc: theta0 expected to be a real scalar" );
        Circle::valueType k0     = getScalarValue( prhs[4+kk], "CircleArc: k0 expected to be a real scalar" );
        Circle::valueType smax   = getScalarValue( prhs[5+kk], "CircleArc: L or s_min expected to be a real scalar" );
        Circle::valueType smin   = 0 ;
        if ( nrhs > 6+kk ) {
          smin = smax ;
          smax = getScalarValue( prhs[6+kk], "CircleArc: s_max expected to be a real scalar" );
        }
        ptr->build( x0, y0, theta0, k0, smin, smax );

      } else if ( nrhs == 3+kk ) {

        Circle::valueType const * p0 = getArrayValues( prhs[1+kk], size0, "CircleArc: p0 expected to be a real vector" );
        Circle::valueType const * p1 = getArrayValues( prhs[2+kk], size1, "CircleArc: p0 expected to be a real vector" );
        Circle::valueType const * p2 = getArrayValues( prhs[3+kk], size2, "CircleArc: p0 expected to be a real vector" );

        ASSERT( size0 == 2 && size2 == 2 && size1 > 0 && size1 <= 2,
                "CircleArc: bad dimension size(p0) = " << size0 <<
                ", size(p1) = " << size1 << ", size(p2) = " << size2 ) ;

        ptr = new Circle::CircleArc();
        if ( size1 == 1 ) ptr->build_G1( p0[0], p0[1], p1[0] /* theta0 */, p2[0], p2[1] ) ;
        else              ptr->build_3P( p0[0], p0[1], p1[0], p1[1],       p2[0], p2[1] ) ;

      } else if ( nrhs == 1 ) {
        // nothing to do
      } else {
        ASSERT(false, "expected " << 4+kk << ", " << 7+kk << " or " << 8+kk << " inputs" );
      }

      plhs[0] = convertPtr2Mat<Circle::CircleArc>(ptr);

    } else if ( cmd == "delete" ) {

      ASSERT(nrhs == 2, "expected 2 inputs");
      ASSERT(nlhs == 0, "expected no output");
      // Destroy the C++ object
      destroyObject<Circle::CircleArc>(prhs[1]);

    } else if ( cmd == "changeOrigin" ) {

      ASSERT(nrhs == 4, "expected 4 inputs");
      ASSERT(nlhs == 0, "expected no output");

      Circle::valueType new_x0 = getScalarValue( arg_in_2, "CircleArc: x0 expected to be a real scalar" );
      Circle::valueType new_y0 = getScalarValue( arg_in_3, "CircleArc: y0 expected to be a real scalar" );

      ptr->changeOrigin( new_x0, new_y0 );

    } else if ( cmd == "translate" ) {

      ASSERT(nrhs == 4, "expected 4 inputs");
      ASSERT(nlhs == 0, "expected no output");

      Circle::valueType tx = getScalarValue( arg_in_2, "CircleArc: tx expected to be a real scalar" );
      Circle::valueType ty = getScalarValue( arg_in_3, "CircleArc: ty expected to be a real scalar" );

      ptr->translate( tx, ty );

    } else if ( cmd == "changeCurvilinearOrigin" ) {

      ASSERT(nrhs == 3, "expected 3 inputs");
      ASSERT(nlhs == 0, "expected no output");

      Circle::valueType new_s = getScalarValue( arg_in_2,  "CircleArc: S expected to be a real scalar" );

      ptr->changeCurvilinearOrigin( new_s );

    } else if ( cmd == "trim" ) {

      ASSERT(nrhs == 4, "expected 4 inputs");
      ASSERT(nlhs == 0, "expected no output");

      Circle::valueType s_begin = getScalarValue( arg_in_2, "CircleArc: s_begin expected to be a real scalar" );
      Circle::valueType s_end   = getScalarValue( arg_in_3, "CircleArc: s_end expected to be a real scalar" );

      ptr->trim( s_begin, s_end );

    } else if ( cmd == "bbTriangle" ) {

      ASSERT(nrhs == 2, "expected 2 inputs");
      ASSERT(nlhs == 4, "expected 4 output");

      arg_out_0 = mxCreateDoubleMatrix(1, 2, mxREAL);
      arg_out_1 = mxCreateDoubleMatrix(1, 2, mxREAL);
      arg_out_2 = mxCreateDoubleMatrix(1, 2, mxREAL);
      arg_out_3 = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);

      bool ok = ptr->bbTriangle( mxGetPr(arg_out_0),
                                 mxGetPr(arg_out_1),
                                 mxGetPr(arg_out_2) );

      *static_cast<mxLogical*>(mxGetData(arg_out_3)) = ok ;

    } else if ( cmd == "to_nurbs" ) {

      Circle::valueType knots[12], Poly[9][3] ;
      Circle::indexType npts = ptr->toNURBS( knots, Poly ); // npt + 2

      static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" } ;
      arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
      mxArray * mx_knots = mxCreateDoubleMatrix(1,npts+2,mxREAL);
      mxArray * mx_Poly  = mxCreateDoubleMatrix(3,npts,mxREAL);

      mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
      mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
      mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
      mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts) );
      mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots );
      mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly );

      double *kb = mxGetPr(mx_knots) ;
      for ( Circle::indexType i = 0 ; i < npts+2 ; ++i ) *kb++ = knots[i] ;

      double *pr = mxGetPr(mx_Poly) ;
      for ( Circle::indexType i = 0 ; i < npts ; ++i ) {
        *pr++ = Poly[i][0] ;
        *pr++ = Poly[i][1] ;
        *pr++ = Poly[i][2] ;
      }
    } else if ( cmd == "length" ) {
      setScalarValue( arg_out_0, ptr->totalLength());
    } else if ( cmd == "eval" ) {
    } else if ( cmd == "eval_D" ) {
    } else if ( cmd == "eval_DD" ) {
    } else if ( cmd == "eval_DDD" ) {

    } else {
      ASSERT(false, "Unknown command: " << cmd );
    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("circleArc failed\n") ;
  }

}
