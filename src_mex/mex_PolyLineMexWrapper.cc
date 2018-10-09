/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "PolyLine.hh"
#include "Line.hh"
#include "Circle.hh"
#include "Biarc.hh"
#include "Clothoid.hh"

#include "mex_utils.hh"

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute polyline\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = PolyLineMexWrapper( 'new', x0, y0, theta0, L );\n" \
"  OBJ = PolyLineMexWrapper( 'new', p0, p1 );\n" \
"\n" \
"  PolyLineMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  PolyLineMexWrapper( 'build', OBJ, x0, y0, theta0, L );\n" \
"  PolyLineMexWrapper( 'build', OBJ, p0, p1 );\n" \
"  PolyLineMexWrapper( 'build', OBJ, p0, theta0, L );\n" \
"  PolyLineMexWrapper( 'copy', OBJ, OBJ1 );\n" \
"  [p1,p2] = PolyLineMexWrapper( 'points', OBJ );\n" \
"\n" \
"  PolyLineMexWrapper( 'translate', OBJ, tx, ty );\n" \
"  PolyLineMexWrapper( 'rotate', OBJ, angle, cx, cy );\n" \
"  PolyLineMexWrapper( 'reverse', OBJ );\n" \
"\n" \
"  burbs = PolyLineMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
"  res = PolyLineMexWrapper( 'xBegin', OBJ );\n" \
"  res = PolyLineMexWrapper( 'yBegin', OBJ );\n" \
"  res = PolyLineMexWrapper( 'xEnd', OBJ );\n" \
"  res = PolyLineMexWrapper( 'yEnd', OBJ );\n" \
"  res = PolyLineMexWrapper( 'length', OBJ );\n" \
"\n" \
"  [X,Y] = PolyLineMexWrapper( 'eval', OBJ, s [,t] );\n" \
"  [X,Y] = PolyLineMexWrapper( 'eval_D', OBJ, s [,t] );\n" \
"  [X,Y] = PolyLineMexWrapper( 'eval_DD', OBJ, s [,t] );\n" \
"  [X,Y] = PolyLineMexWrapper( 'eval_DDD', OBJ, s [,t] );\n" \
"\n" \
"  [s0,s1] = PolyLineMexWrapper( 'intersect', OBJ, OBJ1 );\n" \
"  [d,s]   = PolyLineMexWrapper( 'distance', OBJ, x, y );\n" \
"\n" \
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
  PolyLine *
  DATA_NEW( mxArray * & mx_id ) {
    PolyLine * ptr = new PolyLine();
    mx_id = convertPtr2Mat<PolyLine>(ptr);
    return ptr;
  }

  static
  inline
  PolyLine *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<PolyLine>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<PolyLine>(mx_id);
  }

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);
      mwSize size0, size1;

      bool do_new = cmd == "new";

      PolyLine * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 1, "expected 1 output, nrhs = " << nrhs );
        plhs[0] = convertPtr2Mat<PolyLine>(ptr);

      } else if ( cmd == "build" ) {

        MEX_ASSERT( nlhs == 0, "expected no output, nlhs = " << nlhs );
        #define CMD "PolyLineMexWrapper('build',OBJ,x,y): "

        real_type const * x = getVectorPointer( arg_in_2, size0,
                              CMD "`x` expected to be a real vector" );
        real_type const * y = getVectorPointer( arg_in_3, size1,
                              CMD "`y` expected to be a real vector" );

        MEX_ASSERT( size0 == size1,
                    CMD "expected size(x) = " << size0 <<
                    " = size(y) = " << size1 );

        ptr->build( x, y, size0 );

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "PolyLineMexWrapper('delete',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
        // Destroy the C++ object
        DATA_DELETE( arg_in_1 );
        #undef CMD

      } else if ( cmd == "copy" ) {

        #define CMD "PolyLineMexWrapper('copy',OBJ,OBJ1): "
        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        PolyLine const * LS = convertMat2Ptr<PolyLine>(arg_in_2);
        ptr->copy(*LS);

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "PolyLineMexWrapper('translate',OBJ,t0,t0): "
        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type tx = getScalarValue( arg_in_2, CMD "`tx` expected to be a real scalar" );
        real_type ty = getScalarValue( arg_in_3, CMD "`ty` expected to be a real scalar" );

        ptr->translate( tx, ty );
        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "PolyLineMexWrapper('rotate',OBJ,angle,cx,cy): "
        MEX_ASSERT( nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        real_type angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
        real_type cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
        real_type cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

        ptr->rotate( angle, cx, cy );
        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "PolyLineMexWrapper('reverse',OBJ): "
        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        ptr->reverse();
        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "PolyLineMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output, nlhs = " << nlhs );
          mwSize nrx, ncx, nry, ncy;
          real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx,
                                CMD "`x` expected to be a real vector/matrix" );
          real_type const * y = getMatrixPointer( arg_in_3, nry, ncy,
                                CMD "`y` expected to be a real vector/matrix" );
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          real_type * dst = createMatrixValue( arg_out_0, nrx, ncx );

          mwSize size = nrx*ncx;
          if ( nlhs > 1 ) {
            real_type * s = createMatrixValue( arg_out_1, nrx, ncx );
            for ( mwSize i = 0; i < size; ++i )
              *dst++ = ptr->distance( *x++, *y++, *s++ );
          } else {
            for ( mwSize i = 0; i < size; ++i )
              *dst++ = ptr->distance( *x++, *y++ );
          }
        }

        #undef CMD

      } else if ( cmd == "polygon" ) {

        #define CMD "PolyLineMexWrapper('polygon',OBJ): "
        MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

        real_type * x = createMatrixValue( arg_out_0, ptr->numPoints(), 1 );
        real_type * y = createMatrixValue( arg_out_1, ptr->numPoints(), 1 );

        ptr->polygon( x, y );

        #undef CMD

      } else if ( cmd == "approx" ) {

        #define CMD "PolyLineMexWrapper('approx',OBJ,OBJ1,tol,type): "
        MEX_ASSERT( nrhs == 5, CMD "expected 5 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

        real_type tol = getScalarValue( arg_in_3,
                                        CMD "`tol` expected to be a real scalar" );

        MEX_ASSERT( mxIsChar(arg_in_4), CMD "'type' argument must be a string" );
        string kind = mxArrayToString( arg_in_4 );

        if ( kind == "LineSegment" ) {
          ptr->build( *convertMat2Ptr<LineSegment>(arg_in_2) );
        } else if ( kind == "CircleArc" ) {
          ptr->build( *convertMat2Ptr<CircleArc>(arg_in_2), tol );
        } else if ( kind == "BiArc" ) {
          ptr->build( *convertMat2Ptr<Biarc>(arg_in_2), tol );
        } else if ( kind == "ClothoidCurve" ) {
          ptr->build( *convertMat2Ptr<ClothoidCurve>(arg_in_2), tol );
        } else if ( kind == "ClothoidList" ) {
          ptr->build( *convertMat2Ptr<ClothoidList>(arg_in_2), tol );
        } else {
          MEX_ASSERT( false, CMD "'type' = '" << kind << "' unknown type" );
        }

        #undef CMD

      } else if ( cmd == "intersect" ) {

        #define CMD "PolyLineMexWrapper('intersect',OBJ,OBJ1,type): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

        MEX_ASSERT( mxIsChar(arg_in_3), CMD "'type' argument must be a string" );
        string kind = mxArrayToString( arg_in_3 );

        MEX_ASSERT( kind == "PolyLine",
                    CMD " can intersect only with PolyLine objects" );

        PolyLine const * PL = convertMat2Ptr<PolyLine>(arg_in_2);

        if ( nlhs == 2 ) {
          std::vector<real_type> s1, s2;
          ptr->intersect( *PL, s1, s2 );
          real_type * ps1 = createMatrixValue( arg_out_0, s1.size(), 1 );
          real_type * ps2 = createMatrixValue( arg_out_1, s2.size(), 1 );
          std::copy( s1.begin(), s1.end(), ps1 );
          std::copy( s2.begin(), s2.end(), ps2 );
        } else if ( nlhs == 1 ) {
          bool ok = ptr->intersect( *PL );
          setScalarBool( arg_out_0, ok );
        } else {
          MEX_ASSERT( false, CMD "expected 1 or 2 output, nlhs = " << nlhs );
        }

        #undef CMD

      #define LOOPXY1 \
      real_type *pXY = createMatrixValue( arg_out_0, 2, npts ); \
      for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )

      #define LOOPXY2 \
      real_type *pX = createMatrixValue( arg_out_0, 1, npts ); \
      real_type *pY = createMatrixValue( arg_out_1, 1, npts ); \
      for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

      } else if ( cmd == "eval" ) {

        #define CMD "PolyLineMexWrapper('eval',OBJ,s): "
        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );

        mwSize npts;
        real_type const * s = getVectorPointer( arg_in_2, npts,
                                                CMD "`s` expected to be a real vector" );
        if ( nlhs == 1 ) {
          LOOPXY1 ptr->eval( *s, pXY[0], pXY[1] );
        } else if ( nlhs == 2 ) {
          LOOPXY2 ptr->eval( *s, *pX, *pY );
        } else {
          MEX_ASSERT( false, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
        }
        #undef CMD

      } else if ( cmd == "eval_D" ) {

        #define CMD "PolyLineMexWrapper('eval_D',OBJ,s): "
        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );

        mwSize npts;
        real_type const * s = getVectorPointer( arg_in_2, npts,
                                                CMD "`s` expected to be a real vector" );

        if ( nlhs == 1 ) {
          LOOPXY1 ptr->eval_D( *s, pXY[0], pXY[1] );
        } else if ( nlhs == 2 ) {
          LOOPXY2 ptr->eval_D( *s, *pX, *pY );
        } else {
          MEX_ASSERT( false, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
        }
        #undef CMD

      } else if ( cmd == "eval_DD" ) {

        #define CMD "PolyLineMexWrapper('eval_DD',OBJ,s): "
        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );

        mwSize npts;
        real_type const * s = getVectorPointer( arg_in_2, npts,
                                                CMD "`s` expected to be a real vector" );

        if ( nlhs == 1 ) {
          LOOPXY1 ptr->eval_DD( *s, pXY[0], pXY[1] );
        } else if ( nlhs == 2 ) {
          LOOPXY2 ptr->eval_DD( *s, *pX, *pY );
        } else {
          MEX_ASSERT( false, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
        }
        #undef CMD

      } else if ( cmd == "eval_DDD" ) {

        #define CMD "PolyLineMexWrapper('eval_DDD',OBJ,s): "
        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );

        mwSize npts;
        real_type const * s = getVectorPointer( arg_in_2, npts,
                                                CMD "`s` expected to be a real vector" );

        if ( nlhs == 1 ) {
          LOOPXY1 ptr->eval_DDD( *s, pXY[0], pXY[1] );
        } else if ( nlhs == 2 ) {
          LOOPXY2 ptr->eval_DDD( *s, *pX, *pY );
        } else {
          MEX_ASSERT( false, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
        }
        #undef CMD

      } else if ( cmd == "xBegin" ) {
        setScalarValue( arg_out_0, ptr->xBegin());
      } else if ( cmd == "yBegin" ) {
        setScalarValue( arg_out_0, ptr->yBegin());
      } else if ( cmd == "xEnd" ) {
        setScalarValue( arg_out_0, ptr->xEnd());
      } else if ( cmd == "yEnd" ) {
        setScalarValue( arg_out_0, ptr->yEnd());
      } else if ( cmd == "length" ) {
        setScalarValue( arg_out_0, ptr->length());
      } else {
        MEX_ASSERT( false, "Unknown command: " << cmd );
      }

    } catch ( exception const & e ) {
    	mexErrMsgTxt(e.what());
    } catch (...) {
    	mexErrMsgTxt("Line failed\n");
    }

  }

}
