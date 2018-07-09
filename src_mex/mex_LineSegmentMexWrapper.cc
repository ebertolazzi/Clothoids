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
"  LineSegmentMexWrapper( 'copy', OBJ, OBJ1 ) ;\n" \
"  [p1,p2] = LineSegmentMexWrapper( 'points', OBJ ) ;\n" \
"\n" \
"  LineSegmentMexWrapper( 'changeOrigin', OBJ, x0, y0 ) ;\n" \
"  LineSegmentMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"  LineSegmentMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"  LineSegmentMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"  LineSegmentMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  burbs = LineSegmentMexWrapper( 'to_nurbs', OBJ ) ;\n" \
"\n" \
"  res = LineSegmentMexWrapper( 'xBegin', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'yBegin', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'theta',  OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'sMin', OBJ ) ;\n" \
"  res = LineSegmentMexWrapper( 'sMmax', OBJ ) ;\n" \
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

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "expected 1 output, nlhs = " << nlhs  );

        if ( nrhs == 5 ) {

          #define CMD "LineSegmentMexWrapper('new',x0,y0,theta0,L): "

          real_type x0     = getScalarValue( arg_in_1, CMD "`x0` expected to be a real scalar" );
          real_type y0     = getScalarValue( arg_in_2, CMD "`y0` expected to be a real scalar" );
          real_type theta0 = getScalarValue( arg_in_3, CMD "`theta0` expected to be a real scalar" );
          real_type L      = getScalarValue( arg_in_4, CMD "`L` expected to be a real scalar" );

          ptr->build( x0, y0, theta0, L );

          #undef CMD

        } else if ( nrhs == 3 ) {

          #define CMD "LineSegmentMexWrapper('new',OBJ,p0,p1): "
          real_type const * p0 = getVectorPointer( arg_in_1, size0,
                                 CMD "`p0` expected to be a real vector" );
          real_type const * p1 = getVectorPointer( arg_in_2, size1,
                                 CMD "`p1` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size1 == 2,
                      CMD "bad dimension size(p0) = " << size0 << ", size(p1) = " << size1 ) ;
          #undef CMD

          ptr->build_2P( p0[0], p0[1], p1[0], p1[1] ) ;

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "nrhs = " << nrhs << " expected 3 or 5 inputs, nrhs = " << nrhs  );
        }

        plhs[0] = convertPtr2Mat<LineSegment>(ptr);

      } else if ( cmd == "build" ) {

        MEX_ASSERT( nlhs == 0, "expected no output, nlhs = " << nlhs  );

        if ( nrhs == 6 ) {

          #define CMD "LineSegmentMexWrapper('build',OBJ,x0,y0,theta0,L): "

          real_type x0     = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
          real_type y0     = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );
          real_type theta0 = getScalarValue( arg_in_4, CMD "`theta0` expected to be a real scalar" );
          real_type L      = getScalarValue( arg_in_5, CMD "`L` expected to be a real scalar" );

          ptr->build( x0, y0, theta0, L );

          #undef CMD

        } else if ( nrhs == 4 ) {

          #define CMD "LineSegmentMexWrapper('build',OBJ,p0,p1): "
          real_type const * p0 = getVectorPointer( arg_in_2, size0,
                                 CMD "`p0` expected to be a real vector" );
          real_type const * p1 = getVectorPointer( arg_in_3, size1,
                                 CMD "`p1` expected to be a real vector" );

          MEX_ASSERT( size0 == 2 && size1 == 2,
                      CMD "bad dimension size(p0) = " << size0 << ", size(p1) = " << size1 ) ;
          #undef CMD

          ptr->build_2P( p0[0], p0[1], p1[0], p1[1] ) ;

        } else if ( nrhs == 1 ) {
          // nothing to do
        } else {
          MEX_ASSERT(false, "nrhs = " << nrhs << " expected 4 or 6 inputs, nrhs = " << nrhs  );
        }

      } else if ( cmd == "delete" ) {

        #define CMD "LineSegmentMexWrapper('delete',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 ) ;
        #undef CMD

      } else if ( cmd == "copy" ) {

        #define CMD "LineSegmentMexWrapper('copy',OBJ,OBJ1): "
        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        LineSegment const * LS = convertMat2Ptr<LineSegment>(arg_in_2);
        ptr->copy(*LS) ;

        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "LineSegmentMexWrapper('changeOrigin',OBJ,x0,y0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type new_x0 = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
        real_type new_y0 = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );

        ptr->changeOrigin( new_x0, new_y0 );
        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "LineSegmentMexWrapper('translate',OBJ,t0,t0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type tx = getScalarValue( arg_in_2, CMD "`tx` expected to be a real scalar" );
        real_type ty = getScalarValue( arg_in_3, CMD "`ty` expected to be a real scalar" );

        ptr->translate( tx, ty );
        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "LineSegmentMexWrapper('rotate',OBJ,angle,cx,cy): "
        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
        real_type cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
        real_type cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );
        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "LineSegmentMexWrapper('reverse',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        ptr->reverse();
        #undef CMD

      } else if ( cmd == "trim" ) {

        #define CMD "LineSegmentMexWrapper('trim',OBJ,s_begin,s_end): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type s_begin = getScalarValue( arg_in_2, CMD "`s_begin` expected to be a real scalar" );
        real_type s_end   = getScalarValue( arg_in_3, CMD "`s_end` expected to be a real scalar" );

        ptr->trim( s_begin, s_end );
        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "LineSegmentMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output, nlhs = " << nlhs );
          mwSize nrx, ncx, nry, ncy;
          real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx, "`x` expected to be a real vector/matrix" ) ;
          real_type const * y = getMatrixPointer( arg_in_3, nry, ncy, "`y` expected to be a real vector/matrix" ) ;
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

      } else if ( cmd == "to_nurbs" ) {

        real_type knots[12], Poly[9][3] ;
        int_type  npts = ptr->toNURBS( knots, Poly ); // npt + 2

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
        for ( int_type i = 0 ; i < npts+2 ; ++i ) *kb++ = knots[i] ;

        double *pr = mxGetPr(mx_Poly) ;
        for ( int_type i = 0 ; i < npts ; ++i ) {
          *pr++ = Poly[i][0] ;
          *pr++ = Poly[i][1] ;
          *pr++ = Poly[i][2] ;
        }

      } else if ( cmd == "points" ) {

        #define CMD "LineSegmentMexWrapper('points',OBJ): "
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

        real_type * p1 = createMatrixValue( arg_out_0, 2, 1 ) ;
        real_type * p2 = createMatrixValue( arg_out_1, 2, 1 ) ;

        ptr->p1p2(p1,p2);

        #undef CMD

      } else if ( cmd == "info" ) {

        #define CMD "LineSegmentMexWrapper('info',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs ) ;
        MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = " << nlhs ) ;
        
        ptr->info(cout) ;

        #undef CMD

      } else {

        if ( nrhs == 3 ) {

          #define CMD "LineSegmentMexWrapper('eval*',OBJ,s): "

          mwSize size ;
          real_type const * s = getVectorPointer( arg_in_2, size, CMD "`s` expected to be a real vector" ) ;
          if ( nlhs == 1 ) {
            real_type *pXY = createMatrixValue( arg_out_0, 2,size );
            if ( cmd == "eval" ) {
              for ( mwSize i = 0 ; i < size ; ++i, ++s, pXY += 2 )
                ptr->eval( *s, pXY[0], pXY[1] ) ;
            } else if ( cmd == "eval_D" ) {
              for ( mwSize i = 0 ; i < size ; ++i, ++s, pXY += 2 )
                ptr->eval_D( *s, pXY[0], pXY[1] ) ;
            } else if ( cmd == "eval_DD" ) {
              for ( mwSize i = 0 ; i < size ; ++i, ++s, pXY += 2 )
                ptr->eval_DD( *s, pXY[0], pXY[1] ) ;
            } else if ( cmd == "eval_DDD" ) {
              for ( mwSize i = 0 ; i < size ; ++i, ++s, pXY += 2 )
                ptr->eval_DDD( *s, pXY[0], pXY[1] ) ;
            } else {
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else if ( nlhs == 2 ) {
            real_type *pX = createMatrixValue( arg_out_0, 1,size );
            real_type *pY = createMatrixValue( arg_out_1, 1,size );
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
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else {
            MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs ) ;
          }

          #undef CMD

        } else if ( nrhs == 2 ) {

          if      ( cmd == "xBegin" ) setScalarValue( arg_out_0, ptr->xBegin());
          else if ( cmd == "yBegin" ) setScalarValue( arg_out_0, ptr->yBegin());
          else if ( cmd == "theta"  ) setScalarValue( arg_out_0, ptr->theta());
          else if ( cmd == "length" ) setScalarValue( arg_out_0, ptr->length());
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

