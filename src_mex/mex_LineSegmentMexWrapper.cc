/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Line.hh"
#include "mex_utils.hh"

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute cicle arc\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = LineSegmentMexWrapper( 'new', x0, y0, theta0, L ) ;\n" \
"  OBJ = LineSegmentMexWrapper( 'new', p0, p1 ) ;\n" \
"\n" \
"  LineSegmentMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  LineSegmentMexWrapper( 'build', OBJ, x0, y0, theta0, L ) ;\n" \
"  LineSegmentMexWrapper( 'build', OBJ, p0, p1 ) ;\n" \
"  LineSegmentMexWrapper( 'build', OBJ, p0, theta0, L ) ;\n" \
"\n" \
"  LineSegmentMexWrapper( 'changeOrigin', OBJ, x0, y0 ) ;\n" \
"  LineSegmentMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"  LineSegmentMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"  LineSegmentMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"  LineSegmentMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  burbs = LineSegmentMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"  res = LineSegmentMexWrapper( 'getX0', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'getY0', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'getTheta0', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'getSmin', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'getSmax', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'length', OBJ ) ;\n" \
"\n" \
"  [X,Y] = LineSegmentMexWrapper( 'eval', OBJ, s ) ;\n" \
"  [X,Y] = LineSegmentMexWrapper( 'eval_D', OBJ, s ) ;\n" \
"  [X,Y] = LineSegmentMexWrapper( 'eval_DD', OBJ, s ) ;\n" \
"  [X,Y] = LineSegmentMexWrapper( 'eval_DDD', OBJ, s ) ;\n" \
"\n" \
"  [d,s] = LineSegmentMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"\n" \
"  nurbs = LineSegmentMexWrapper( 'to_nurbs', OBJ ) ;\n" \
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

  using namespace std ;

  static
  LineSegment *
  DATA_NEW( mxArray * & mx_id ) {
    LineSegment * ptr = new LineSegment();
    mx_id = convertPtr2Mat<LineSegment>(ptr);
    return ptr ;
  }

  static
  inline
  LineSegment *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<LineSegment>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<LineSegment>(mx_id);
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
      mwSize size0, size1 ;

      bool do_new = cmd == "new" ;

      LineSegment * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build" ) {

        indexType kk = do_new ? 0 : 1 ;

        #define BUILD "LineSegmentMexWrapper('build',OBJ,x0,y0,theta0,L): "
        MEX_ASSERT( nlhs == 1, BUILD "expected 1 output" );

        if ( nrhs == 5+kk ) {

          valueType x0     = getScalarValue( prhs[1+kk], BUILD "`x0` expected to be a real scalar" );
          valueType y0     = getScalarValue( prhs[2+kk], BUILD "`y0` expected to be a real scalar" );
          valueType theta0 = getScalarValue( prhs[3+kk], BUILD "`theta0` expected to be a real scalar" );
          valueType L      = getScalarValue( prhs[4+kk], BUILD "`L` expected to be a real scalar" );

          ptr->build( x0, y0, theta0, L );

        } else if ( nrhs == 3+kk ) {

          #undef BUILD
          #define BUILD "LineSegmentMexWrapper('build',OBJ,p0,p1): "
          valueType const * p0 = getVectorPointer( prhs[1+kk], size0, BUILD "`p0` expected to be a real vector" );
          valueType const * p1 = getVectorPointer( prhs[2+kk], size1, BUILD "`p1` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size1 == 2,
                      BUILD "bad dimension size(p0) = " << size0 << ", size(p1) = " << size1 ) ;
          #undef BUILD

          ptr->build_2P( p0[0], p0[1], p1[0], p1[1] ) ;

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "nrhs = " << nrhs << " expected " << 2+kk << " or " << 5+kk << " inputs" );
        }

        plhs[0] = convertPtr2Mat<LineSegment>(ptr);

      } else if ( cmd == "delete" ) {

        #define CMD "LineSegmentMexWrapper('delete',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 ) ;
        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "LineSegmentMexWrapper('changeOrigin',OBJ,x0,y0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");

        valueType new_x0 = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
        valueType new_y0 = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );

        ptr->changeOrigin( new_x0, new_y0 );
        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "LineSegmentMexWrapper('translate',OBJ,t0,t0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");

        valueType tx = getScalarValue( arg_in_2, CMD "`tx` expected to be a real scalar" );
        valueType ty = getScalarValue( arg_in_3, CMD "`ty` expected to be a real scalar" );

        ptr->translate( tx, ty );
        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "LineSegmentMexWrapper('rotate',OBJ,angle,cx,cy): "
        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");

        valueType angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
        valueType cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
        valueType cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );
        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "LineSegmentMexWrapper('reverse',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");
        ptr->reverse();
        #undef CMD

      } else if ( cmd == "trim" ) {

        #define CMD "LineSegmentMexWrapper('trim',OBJ,s_begin,s_end): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");

        valueType s_begin = getScalarValue( arg_in_2, CMD "`s_begin` expected to be a real scalar" );
        valueType s_end   = getScalarValue( arg_in_3, CMD "`s_end` expected to be a real scalar" );

        ptr->trim( s_begin, s_end );
        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "LineSegmentMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input");
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output");
          mwSize nrx, ncx, nry, ncy;
          valueType const * x = getMatrixPointer( arg_in_2, nrx, ncx, "`x` expected to be a real vector/matrix" ) ;
          valueType const * y = getMatrixPointer( arg_in_3, nry, ncy, "`y` expected to be a real vector/matrix" ) ;
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          valueType * dst = createMatrixValue( arg_out_0, nrx, ncx ) ;

          mwSize size = nrx*ncx ;
          if ( nlhs > 1 ) {
            valueType * s = createMatrixValue( arg_out_1, nrx, ncx ) ;
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distance( *x++, *y++, *s++ ) ;
          } else {
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distance( *x++, *y++ ) ;
          }
        }

        #undef CMD

      } else if ( cmd == "to_nurbs" ) {

        valueType knots[12], Poly[9][3] ;
        indexType npts = ptr->toNURBS( knots, Poly ); // npt + 2

        static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" } ;
        arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
        mxArray * mx_knots = mxCreateDoubleMatrix(1,npts+2,mxREAL);
        mxArray * mx_Poly  = mxCreateDoubleMatrix(3,npts,mxREAL);

        mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
        mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(2) );
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
          double const * s = getVectorPointer( arg_in_2, size, "Line: s expected to be a real vector" ) ;
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

          if      ( cmd == "getX0"     ) setScalarValue( arg_out_0, ptr->getX0());
          else if ( cmd == "getY0"     ) setScalarValue( arg_out_0, ptr->getY0());
          else if ( cmd == "getTheta0" ) setScalarValue( arg_out_0, ptr->getTheta0());
          else if ( cmd == "length"    ) setScalarValue( arg_out_0, ptr->getL());
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
    	mexErrMsgTxt("Line failed\n") ;
    }

  }

}

