/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Circle.hh"
#include "mex_utils.hh"

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute cicle arc\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = CircleMexWrapper( 'new', x0, y0, theta0, kur, L ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', x0, y0, theta0, kur, smin, smax ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', p0, p1, p2 ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', p0, theta0, p1 ) ;\n" \
"\n" \
"  CircleMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  CircleMexWrapper( 'build', OBJ, x0, y0, theta0, kur, L ) ;\n" \
"  CircleMexWrapper( 'build', OBJ, x0, y0, theta0, kur, smin, smax );\n" \
"  CircleMexWrapper( 'build', OBJ, p0, p1, p2 ) ;\n" \
"  CircleMexWrapper( 'build', OBJ, p0, theta0, p1 ) ;\n" \
"\n" \
"  CircleMexWrapper( 'changeOrigin', OBJ, x0, y0 ) ;\n" \
"  CircleMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"  CircleMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"  CircleMexWrapper( 'changeCurvilinearOrigin', OBJ, s0 ) ;\n" \
"  CircleMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"  CircleMexWrapper( 'scale', OBJ, scale ) ;\n" \
"  CircleMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  [p0,p1,p2,ok] = CircleMexWrapper( 'bbTriangle', OBJ ) ;\n" \
"\n" \
"  burbs = CircleMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"  res = CircleMexWrapper( 'getX0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getY0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getTheta0', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getKappa', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getSmin', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'getSmax', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'length', OBJ ) ;\n" \
"\n" \
"  [X,Y] = CircleMexWrapper( 'eval', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_D', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_DD', OBJ, s ) ;\n" \
"  [X,Y] = CircleMexWrapper( 'eval_DDD', OBJ, s ) ;\n" \
"\n" \
"  [X,Y] = CircleMexWrapper( 'theta', OBJ, s ) ;\n" \
"\n" \
"  nurbs = CircleMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"%==========================================================================%\n" \
"%                                                                          %\n" \
"%  Autor: Enrico Bertolazzi                                                %\n" \
"%         Department of Industrial Engineering                             %\n" \
"%         University of Trento                                             %\n" \
"%         enrico.bertolazzi@unitn.it                                       %\n" \
"%                                                                          %\n" \
"%==========================================================================%\n"

namespace G2lib {

  using namespace std;

  static
  CircleArc *
  DATA_NEW( mxArray * & mx_id ) {
    CircleArc * ptr = new CircleArc();
    mx_id = convertPtr2Mat<CircleArc>(ptr);
    return ptr ;
  }

  static
  inline
  CircleArc *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<CircleArc>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    CircleArc * ptr = convertMat2Ptr<CircleArc>(mx_id);
    delete ptr ;
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

      MEX_ASSERT( mxIsChar(arg_in_0), "First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;
      mwSize size0, size1, size2 ;

      bool do_new = cmd == "new" ;
      CircleArc * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build" ) {

        indexType kk = do_new ? 0 : 1 ;

        MEX_ASSERT( nlhs == 1, "expected 1 output" );

        if ( nrhs == 6+kk || nrhs == 7+kk ) {

          valueType x0     = getScalarValue( prhs[1+kk], "CircleArc: `x0` expected to be a real scalar" );
          valueType y0     = getScalarValue( prhs[2+kk], "CircleArc: `y0` expected to be a real scalar" );
          valueType theta0 = getScalarValue( prhs[3+kk], "CircleArc: `theta0` expected to be a real scalar" );
          valueType k0     = getScalarValue( prhs[4+kk], "CircleArc: `k0` expected to be a real scalar" );
          valueType smin = 0, smax = 0 ;
          if ( nrhs == 7+kk ) {
            smin = getScalarValue( prhs[5+kk], "CircleArc: `s_min` expected to be a real scalar" );
            smax = getScalarValue( prhs[6+kk], "CircleArc: `s_max` expected to be a real scalar" );
          } else {
            smax = getScalarValue( prhs[6+kk], "CircleArc: `L` expected to be a real scalar" );
          }
          ptr->build( x0, y0, theta0, k0, smin, smax );

        } else if ( nrhs == 3+kk ) {

          valueType const * p0 = getVectorPointer( prhs[1+kk], size0, "CircleArc: `p0` expected to be a real vector" );
          valueType const * p1 = getVectorPointer( prhs[2+kk], size1, "CircleArc: `p1` expected to be a real vector" );
          valueType const * p2 = getVectorPointer( prhs[3+kk], size2, "CircleArc: `p2` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size2 == 2 && size1 > 0 && size1 <= 2,
                      "CircleArc: bad dimension size(p0) = " << size0 <<
                      ", size(p1) = " << size1 << ", size(p2) = " << size2 ) ;

          if ( size1 == 1 ) ptr->build_G1( p0[0], p0[1], p1[0] /* theta0 */, p2[0], p2[1] ) ;
          else              ptr->build_3P( p0[0], p0[1], p1[0], p1[1],       p2[0], p2[1] ) ;

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "expected " << 4+kk << ", " << 7+kk << " or " << 8+kk << " inputs" );
        }

        plhs[0] = convertPtr2Mat<CircleArc>(ptr);

      } else if ( cmd == "delete" ) {

        MEX_ASSERT(nrhs == 2, "expected 2 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 ) ;

      } else if ( cmd == "changeOrigin" ) {

        MEX_ASSERT(nrhs == 4, "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType new_x0 = getScalarValue( arg_in_2, "CircleArc: `x0` expected to be a real scalar" );
        valueType new_y0 = getScalarValue( arg_in_3, "CircleArc: `y0` expected to be a real scalar" );

        ptr->changeOrigin( new_x0, new_y0 );

      } else if ( cmd == "changeCurvilinearOrigin" ) {

        MEX_ASSERT(nrhs == 3, "expected 3 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType new_s = getScalarValue( arg_in_2,  "CircleArc: `S` expected to be a real scalar" );

        ptr->changeCurvilinearOrigin( new_s );

      } else if ( cmd == "translate" ) {

        MEX_ASSERT(nrhs == 4, "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType tx = getScalarValue( arg_in_2, "CircleArc: `tx` expected to be a real scalar" );
        valueType ty = getScalarValue( arg_in_3, "CircleArc: `ty` expected to be a real scalar" );

        ptr->translate( tx, ty );

      } else if ( cmd == "rotate" ) {

        MEX_ASSERT(nrhs == 5, "expected 5 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType angle = getScalarValue( arg_in_2, "CircleArc: `angle` expected to be a real scalar" );
        valueType cx    = getScalarValue( arg_in_3, "CircleArc: `cx` expected to be a real scalar" );
        valueType cy    = getScalarValue( arg_in_4, "CircleArc: `cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );

      } else if ( cmd == "scale" ) {

        MEX_ASSERT(nrhs == 3, "expected 3 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType sc = getScalarValue( arg_in_2, "CircleArc: `scale` expected to be a real scalar" );
        ptr->scale( sc );

      } else if ( cmd == "reverse" ) {

        MEX_ASSERT(nrhs == 2, "expected 2 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");
        ptr->reverse();

      } else if ( cmd == "trim" ) {

        MEX_ASSERT(nrhs == 4, "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, "expected no output");

        valueType s_begin = getScalarValue( arg_in_2, "CircleArc: `s_begin` expected to be a real scalar" );
        valueType s_end   = getScalarValue( arg_in_3, "CircleArc: `s_end` expected to be a real scalar" );

        ptr->trim( s_begin, s_end );

      } else if ( cmd == "theta" ) {

        MEX_ASSERT(nrhs == 3, "expected 3 inputs");
        MEX_ASSERT(nlhs == 1, "expected 1 output");

        valueType s = getScalarValue( arg_in_2, "CircleArc: `s` expected to be a real scalar" );
        setScalarValue( arg_out_0, ptr->theta( s ) ) ;

      } else if ( cmd == "bbTriangle" ) {

        MEX_ASSERT(nrhs == 2, "expected 2 inputs");
        MEX_ASSERT(nlhs == 4, "expected 4 output");

        double * p0 = createMatrixValue( arg_out_0, 1, 2 );
        double * p1 = createMatrixValue( arg_out_1, 1, 2 );
        double * p2 = createMatrixValue( arg_out_2, 1, 2 );

        bool ok = ptr->bbTriangle( p0, p1, p2 );

        arg_out_3 = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
        *static_cast<mxLogical*>(mxGetData(arg_out_3)) = ok ;

      } else if ( cmd == "to_nurbs" ) {

        valueType knots[12], Poly[9][3] ;
        indexType npts = ptr->toNURBS( knots, Poly ); // npt + 2

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
        for ( indexType i = 0 ; i < npts+2 ; ++i ) *kb++ = knots[i] ;

        double *pr = mxGetPr(mx_Poly) ;
        for ( indexType i = 0 ; i < npts ; ++i ) {
          *pr++ = Poly[i][0] ;
          *pr++ = Poly[i][1] ;
          *pr++ = Poly[i][2] ;
        }
      } else {
        if ( nrhs == 3 ) {
          mwSize size ;
          double const * s = getVectorPointer( arg_in_2, size,  "CircleArc: s expected to be a real vector" ) ;
          double *pX = createMatrixValue( arg_out_0, 1,size );
          double *pY = createMatrixValue( arg_out_1, 1,size );
          if ( cmd == "eval" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_D" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_D( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_DD" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_DD( *s, *pX, *pY ) ;
          } else if ( cmd == "eval_DDD" ) {
            for ( mwSize i = 0 ; i < size ; ++i, ++s, ++pX, ++pY )
              ptr->eval_DDD( *s, *pX, *pY ) ;
          } else {
            MEX_ASSERT(false, "Unknown command: " << cmd );
          }
        } else if ( nrhs == 3 ) {
          if      ( cmd == "getX0"     ) setScalarValue( arg_out_0, ptr->getX0());
          else if ( cmd == "getY0"     ) setScalarValue( arg_out_0, ptr->getY0());
          else if ( cmd == "getTheta0" ) setScalarValue( arg_out_0, ptr->getTheta0());
          else if ( cmd == "getKappa"  ) setScalarValue( arg_out_0, ptr->getKappa());
          else if ( cmd == "getSmin"   ) setScalarValue( arg_out_0, ptr->getSmin());
          else if ( cmd == "getSmax"   ) setScalarValue( arg_out_0, ptr->getSmax());
          else if ( cmd == "length"    ) setScalarValue( arg_out_0, ptr->totalLength());
          else {
            MEX_ASSERT(false, "Unknown command: " << cmd );
          }
        } else {
          MEX_ASSERT(false, "Unknown command: " << cmd );
        }
      }

    } catch ( exception const & e ) {
    	mexErrMsgTxt(e.what()) ;
    } catch (...) {
    	mexErrMsgTxt("circleArc failed\n") ;
    }
  }
}
