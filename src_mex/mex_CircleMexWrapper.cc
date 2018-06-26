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
"  OBJ = CircleMexWrapper( 'new' ) ;\n" \
"  OBJ = CircleMexWrapper( 'new', x0, y0, theta0, kur, L ) ;\n" \
"\n" \
"  CircleMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  CircleMexWrapper( 'build', OBJ, x0, y0, theta0, kur, L ) ;\n" \
"  CircleMexWrapper( 'build_G1', OBJ, p0, theta0, p1 ) ;\n" \
"  CircleMexWrapper( 'build_G1', OBJ, x0, y0, theta0, x1, y1 ) ;\n" \
"  CircleMexWrapper( 'build_3P', OBJ, p0, p1, p2 ) ;\n" \
"  CircleMexWrapper( 'build_3P', OBJ, x0, y0, x1, y1, x2, y2 ) ;\n" \
"\n" \
"  CircleMexWrapper( 'changeOrigin', OBJ, x0, y0 ) ;\n" \
"  CircleMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"  CircleMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"  CircleMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"  CircleMexWrapper( 'scale', OBJ, scale ) ;\n" \
"  CircleMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  [d,s] = CircleMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"\n" \
"  [p0,p1,p2,ok] = CircleMexWrapper( 'bbTriangle', OBJ ) ;\n" \
"\n" \
"  burbs = CircleMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"  res = CircleMexWrapper( 'xBegin', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'yBegin', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'thetaBegin', OBJ ) ;\n" \
"  res = CircleMexWrapper( 'kappa', OBJ ) ;\n" \
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
    destroyObject<CircleArc>(mx_id);
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

      MEX_ASSERT( mxIsChar(arg_in_0), "CircleMexWrapper(...): First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;
      mwSize size0, size1, size2 ;

      bool do_new = cmd == "new" ;
      CircleArc * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build" ) {

        int_type kk = do_new ? 0 : 1 ;

        if ( do_new ) {
          MEX_ASSERT( nlhs == 1, "CircleMexWrapper, expected 1 output, nlhs = " << nlhs  );
        } else {
          MEX_ASSERT( nlhs == 0, "CircleMexWrapper, expected no output, nlhs = " << nlhs  );
        }

        if ( nrhs == 6+kk ) {

          #define CMD "CircleMexWrapper(x0,y0,theta0,k0,L): "
          real_type x0     = getScalarValue( prhs[1+kk], CMD "`x0` expected to be a real scalar" );
          real_type y0     = getScalarValue( prhs[2+kk], CMD "`y0` expected to be a real scalar" );
          real_type theta0 = getScalarValue( prhs[3+kk], CMD "`theta0` expected to be a real scalar" );
          real_type k0     = getScalarValue( prhs[4+kk], CMD "`k0` expected to be a real scalar" );
          real_type L      = getScalarValue( prhs[5+kk], CMD "`L` expected to be a real scalar" );

          ptr->build( x0, y0, theta0, k0, L );
          #undef CMD

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "CircleArc, expected 1 or " << 6+kk << " inputs, nrhs = " << nrhs );
        }

        plhs[0] = convertPtr2Mat<CircleArc>(ptr);

      } else if ( cmd == "build_G1" ) {

        MEX_ASSERT( nlhs == 0 || nlhs ==1, "CircleMexWrapper, expected 1 or no output, nlhs = " << nlhs );

        real_type x0(0), y0(0), x1(0), y1(0), theta0(0);
        if ( nrhs == 5 ) {

          #define CMD "CircleMexWrapper('build_G1',OBJ,p0,theta0,p1): "
          real_type const * p0 = getVectorPointer( arg_in_2, size0, CMD "`p0` expected to be a real vector" );
          theta0 = getScalarValue( arg_in_3, CMD "`theta0` expected to be a real vector" );
          real_type const * p1 = getVectorPointer( arg_in_4, size1, CMD "`p1` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size1 == 2,
                      CMD "bad dimension size(p0) = " << size0 << ", size(p1) = " << size1 ) ;

          #undef CMD

          x0 = p0[0] ; y0 = p0[1] ;
          x1 = p1[0] ; y1 = p1[1] ;

        } else if ( nrhs == 7 ) {

          #define CMD "CircleMexWrapper('build_G1',OBJ,x0,x1,theta0,x1,y1): "
          x0     = getScalarValue( arg_in_2,CMD "`x0` expected to be a scalar value" );
          y0     = getScalarValue( arg_in_3,CMD "`y0` expected to be a scalar value" );
          theta0 = getScalarValue( arg_in_4,CMD "`theta0` expected to be a scalar value" );
          x1     = getScalarValue( arg_in_5,CMD "`x1` expected to be a scalar value" );
          y1     = getScalarValue( arg_in_6,CMD "`y1` expected to be a scalar value" );
          #undef CMD
        } else {
          MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = " << nrhs );
        }

        bool ok = ptr->build_G1( x0, y0, theta0, x1, y1 ) ;
        if ( nlhs == 1 ) {
          arg_out_0 = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
          *static_cast<mxLogical*>(mxGetData(arg_out_0)) = ok ;
        }
        #undef CMD

      } else if ( cmd == "build_3P" ) {

        MEX_ASSERT( nlhs == 0 || nlhs ==1, "CircleMexWrapper, expected 1 or no output, nlhs = " << nlhs  );

        real_type x0(0), y0(0), x1(0), y1(0), x2(0), y2(0);
        if ( nrhs == 5 ) {

          #define CMD "CircleMexWrapper('build_G1',OBJ,p0,p1,p2): "
          real_type const * p0 = getVectorPointer( arg_in_2, size0, CMD "`p0` expected to be a real vector" );
          real_type const * p1 = getVectorPointer( arg_in_3, size1, CMD "`p1` expected to be a real vector" );
          real_type const * p2 = getVectorPointer( arg_in_4, size2, CMD "`p2` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size1 == 2 && size2 == 2,
                      CMD "bad dimension size(p0) = " << size0 <<
                      ", size(p1) = " << size1 << ", size(p2) = " << size2 ) ;

          #undef CMD

          x0 = p0[0] ; y0 = p0[1] ;
          x1 = p1[0] ; y1 = p1[1] ;
          x2 = p2[0] ; y2 = p2[1] ;

        } else if ( nrhs == 8 ) {

          #define CMD "CircleMexWrapper('build_G1',OBJ,x0,x1,x1,y1,x2,y2): "
          x0     = getScalarValue( arg_in_2,CMD "`x0` expected to be a scalar value" );
          y0     = getScalarValue( arg_in_3,CMD "`y0` expected to be a scalar value" );
          x1     = getScalarValue( arg_in_4,CMD "`x1` expected to be a scalar value" );
          y1     = getScalarValue( arg_in_5,CMD "`y1` expected to be a scalar value" );
          x2     = getScalarValue( arg_in_6,CMD "`x2` expected to be a scalar value" );
          y2     = getScalarValue( arg_in_7,CMD "`y2` expected to be a scalar value" );
          #undef CMD
        } else {
          MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = " << nrhs );
        }

        bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 ) ;
        if ( nlhs == 1 ) setScalarBool(arg_out_0,ok);

      } else if ( cmd == "delete" ) {

        #define CMD "CircleMexWrapper('delete',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 ) ;
        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "CircleMexWrapper('changeOrigin',OBJ,x0,y0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type new_x0 = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
        real_type new_y0 = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );

        ptr->changeOrigin( new_x0, new_y0 );
        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "CircleMexWrapper('translate',OBJ,tx,ty): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type tx = getScalarValue( arg_in_2, CMD "`tx` expected to be a real scalar" );
        real_type ty = getScalarValue( arg_in_3, CMD "`ty` expected to be a real scalar" );

        ptr->translate( tx, ty );
        #undef CMD

      } else if ( cmd == "changeCurvilinearOrigin" ) {

        #define CMD "CircleMexWrapper('changeCurvilinearOrigin',OBJ,s0,L): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );

        real_type s0 = getScalarValue(arg_in_2,CMD "Error in reading s0") ;
        real_type L  = getScalarValue(arg_in_3,CMD "Error in reading L") ;
        ptr->changeCurvilinearOrigin(s0,L);
        #undef CMD

      } else if ( cmd == "rotate" ) {
        #define CMD "CircleMexWrapper('rotate',OBJ,angle,cx,cy): "

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
        real_type cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
        real_type cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );
        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "CircleMexWrapper('scale',OBJ,scale): "
        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type sc = getScalarValue( arg_in_2, CMD "`scale` expected to be a real scalar" );
        ptr->scale( sc );
        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "CircleMexWrapper('reverse',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        ptr->reverse();
        #undef CMD

      } else if ( cmd == "trim" ) {

        #define CMD "CircleMexWrapper('trim',OBJ,s_begin,s_end): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type s_begin = getScalarValue( arg_in_2, CMD "`s_begin` expected to be a real scalar" );
        real_type s_end   = getScalarValue( arg_in_3, CMD "`s_end` expected to be a real scalar" );

        ptr->trim( s_begin, s_end );
        #undef CMD

      } else if ( cmd == "theta" ) {

        #define CMD "CircleMexWrapper('theta',OBJ,s): "
        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

        real_type s = getScalarValue( arg_in_2, CMD "`s` expected to be a real scalar" );
        setScalarValue( arg_out_0, ptr->theta( s ) ) ;
        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "CircleMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output, nlhs = " << nlhs );
          mwSize nrx, ncx, nry, ncy;
          real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" ) ;
          real_type const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" ) ;
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          real_type * dst = createMatrixValue( arg_out_0, nrx, ncx ) ;

          mwSize size = nrx*ncx ;
          if ( nlhs > 1 ) {
            real_type * s = createMatrixValue( arg_out_1, nrx, ncx ) ;
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distance( *x++, *y++, *s++ ) ;
          } else {
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distance( *x++, *y++ ) ;
          }
        }
        #undef CMD

      } else if ( cmd == "bbTriangle" ) {

        #define CMD "CircleMexWrapper('bbTriangle',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 4, CMD "expected 4 output, nlhs = " << nlhs );

        double * p0 = createMatrixValue( arg_out_0, 1, 2 );
        double * p1 = createMatrixValue( arg_out_1, 1, 2 );
        double * p2 = createMatrixValue( arg_out_2, 1, 2 );

        bool ok = ptr->bbTriangle( p0, p1, p2 );
        setScalarBool(arg_out_3,ok);

        #undef CMD

      } else if ( cmd == "to_nurbs" ) {

        #define CMD "CircleMexWrapper('to_nurbs',OBJ): "

        int_type npts = ptr->toNURBS( nullptr, nullptr, true );

        mxArray * mx_knots, * mx_Poly ;
        double * knots = createMatrixValue( mx_knots, 1, npts+3 );
        double * poly  = createMatrixValue( mx_Poly,  3, npts );

        ptr->toNURBS( knots, poly, false );

        static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" } ;
        arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);

        mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
        mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
        mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
        mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts) );
        mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots );
        mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly );

        #undef CMD
      } else {

        #define CMD "CircleMexWrapper('eval*'',OBJ,s): "

        if ( nrhs == 3 ) {
          mwSize size ;
          double const * s = getVectorPointer( arg_in_2, size, CMD "`s` expected to be a real vector" ) ;
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
        } else if ( nrhs == 2 ) {
          if      ( cmd == "xBegin"      ) setScalarValue( arg_out_0, ptr->xBegin());
          else if ( cmd == "xEnd"        ) setScalarValue( arg_out_0, ptr->xEnd());
          else if ( cmd == "yBegin"      ) setScalarValue( arg_out_0, ptr->yBegin());
          else if ( cmd == "yEnd"        ) setScalarValue( arg_out_0, ptr->yEnd());
          else if ( cmd == "thetaBegin"  ) setScalarValue( arg_out_0, ptr->thetaBegin());
          else if ( cmd == "thetaEnd"    ) setScalarValue( arg_out_0, ptr->thetaEnd());
          else if ( cmd == "kappa"       ) setScalarValue( arg_out_0, ptr->kappa());
          else if ( cmd == "length"      ) setScalarValue( arg_out_0, ptr->length());
          else {
            MEX_ASSERT(false, "CircleMexWrapper unknown command: " << cmd );
          }
        } else {
          MEX_ASSERT(false, "CircleMexWrapper unknown command: " << cmd );
        }

        #undef CMD

      }

    } catch ( exception const & e ) {
    	mexErrMsgTxt(e.what()) ;
    } catch (...) {
    	mexErrMsgTxt("CircleMexWrapper failed\n") ;
    }
  }
}
