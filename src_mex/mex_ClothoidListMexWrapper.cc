/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "Triangle2D.hh"
#include "mex_utils.hh"

#include <fstream>

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"ClothoidListMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = ClothoidListMexWrapper( 'new' );\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    ClothoidListMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  - Build:\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, CLOT );\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, kappa0, dkappa, L );\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, x0, y0, theta0, kappa0, dkappa, L );\n" \
"    ClothoidListMexWrapper( 'push_back_G1', OBJ, x1, y1, theta1 );\n" \
"    ClothoidListMexWrapper( 'push_back_G1', OBJ, x0, y0, theta0, x1, y1, theta1 );\n" \
"    ClothoidListMexWrapper( 'copy', OBJ, OBJ1 );\n" \
"\n" \
"  - Eval:\n" \
"    [x,y,theta,kappa] = ClothoidListMexWrapper( 'evaluate', OBJ, ss );\n" \
"    s = ClothoidListMexWrapper( 'sBegin', OBJ );\n" \
"    s = ClothoidListMexWrapper( 'sEnd',   OBJ );\n" \
"    x = ClothoidListMexWrapper( 'xBegin', OBJ );\n" \
"    x = ClothoidListMexWrapper( 'xEnd',   OBJ );\n" \
"    y = ClothoidListMexWrapper( 'yBegin', OBJ );\n" \
"    y = ClothoidListMexWrapper( 'yEnd',   OBJ );\n" \
"    theta = ClothoidListMexWrapper( 'thetaBegin', OBJ );\n" \
"    theta = ClothoidListMexWrapper( 'thetaEnd', OBJ );\n" \
"    kappa = ClothoidListMexWrapper( 'kappaBegin', OBJ );\n" \
"    kappa = ClothoidListMexWrapper( 'kappaEnd', OBJ );\n" \
"\n" \
"    [x,y]           = ClothoidListMexWrapper( 'eval', OBJ, s[, t] );\n" \
"    [x_D,y_D]       = ClothoidListMexWrapper( 'eval_D', OBJ, s[, t] );\n" \
"    [x_DD,y_DD]     = ClothoidListMexWrapper( 'eval_DD', OBJ, s[, t] );\n" \
"    [x_DDD,y_DDD]   = ClothoidListMexWrapper( 'eval_DDD', OBJ, s[, t] );\n" \
"    [s,theta,kappa] = ClothoidListMexWrapper( 'getSTK', OBJ );\n" \
"    [x,y]           = ClothoidListMexWrapper( 'getXY', OBJ );\n" \
"\n" \
"  - Transform:\n" \
"    ClothoidListMexWrapper( 'changeOrigin', OBJ, newX0, newY0 );\n" \
"    ClothoidListMexWrapper( 'rotate', OBJ, angle, cx, cy );\n" \
"    ClothoidListMexWrapper( 'translate', OBJ, tx, ty );\n" \
"    ClothoidListMexWrapper( 'scale', OBJ, scaling );\n" \
"    ClothoidListMexWrapper( 'reverse', OBJ );\n" \
"\n" \
"  - Distance:\n" \
"    [X,Y,s,dst] = ClothoidListMexWrapper( 'closestPoint', OBJ, x, y );\n" \
"    [dst,s]     = ClothoidListMexWrapper( 'distance', OBJ, x, y );\n" \
"  - Intersection:\n" \
"    [s1,s2] = ClothoidCurveMexWrapper( 'intersect_line', OBJ, OBJ2 );%\n" \
"    [s1,s2] = ClothoidCurveMexWrapper( 'intersect_circle', OBJ, OBJ2 );%\n" \
"    [s1,s2] = ClothoidCurveMexWrapper( 'intersect_clothoid', OBJ, OBJ2 );%\n" \
"    [s1,s2] = ClothoidCurveMexWrapper( 'intersect_clothoid_list', OBJ, OBJ2 );%\n" \
"\n" \
"  - Bounding Box:\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size );%\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size, offs );%\n" \
"\n" \
"  - G2 spline:\n" \
"    ok = ClothoidListMexWrapper( 'build_3arcG2' || 'build_2arcG2' || 'build_3arcCLC, ...\n" \
"                                 OBJ, ...\n" \
"                                 x0, y0, theta0, kappa0, ...\n" \
"                                 x1, y1, theta1, kappa1 );%\n" \
"    ok = ClothoidListMexWrapper( 'build_3arcG2fixed', OBJ, ...\n" \
"                                 s0, x0, y0, theta0, kappa0, ...\n" \
"                                 s1, x1, y1, theta1, kappa1 );%\n" \
"    ok = ClothoidListMexWrapper( 'build_G1', OBJ, x, y [,theta] );%\n" \
"    [theta,ok] = ClothoidListMexWrapper( 'build_theta', OBJ, x, y );%\n" \
"    dtheta = ClothoidListMexWrapper( 'deltaTheta', OBJ );%\n" \
"    dkappa = ClothoidListMexWrapper( 'deltaKappa', OBJ );%\n" \
"\n" \
"=====================================================================================\n" \
"\n" \
"Autor: Enrico Bertolazzi\n" \
"  Department of Industrial Engineering\n" \
"  University of Trento\n" \
"  enrico.bertolazzi@unitn.it\n" \
"\n" \
"=====================================================================================\n"

namespace G2lib {

  using namespace std;

  static
  ClothoidList *
  DATA_NEW( mxArray * & mx_id ) {
    ClothoidList * ptr = new ClothoidList();
    mx_id = convertPtr2Mat<ClothoidList>(ptr);
    return ptr;
  }

  static
  inline
  ClothoidList *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<ClothoidList>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<ClothoidList>(mx_id);
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

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidListMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      ClothoidList * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "ClothoidListMexWrapper, expected 1 output, nlhs=" << nlhs );

      } else if ( cmd == "push_back" ) {

        #define CMD "ClothoidListMexWrapper('push_back',OBJ,CLOT|[kappa0,dkappa,L]|[x0,y0,theta0,kappa0,dkappa,L]): "

        if ( nrhs == 8 ) {
          real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
          real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
          real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
          real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading kappa0" );
          real_type dkappa = getScalarValue( arg_in_6, CMD "Error in reading dkappa" );
          real_type L      = getScalarValue( arg_in_7, CMD "Error in reading L" );
          ptr->push_back( x0, y0, theta0, kappa0, dkappa, L );
        } else if ( nrhs == 5 ) {
          real_type kappa0 = getScalarValue( arg_in_2, CMD "Error in reading kappa0" );
          real_type dkappa = getScalarValue( arg_in_3, CMD "Error in reading dkappa" );
          real_type L      = getScalarValue( arg_in_4, CMD "Error in reading L" );
          ptr->push_back( kappa0, dkappa, L );
        } else if ( nrhs == 3 ) {
          ClothoidCurve * cc = convertMat2Ptr<ClothoidCurve>(arg_in_2);
          ptr->push_back( *cc );
        } else {
          MEX_ASSERT( false, CMD "expected 3, 5 or 8 inputs nrhs = " << nrhs );
        }

        #undef CMD

      } else if ( cmd == "push_back_G1" ) {

        #define CMD "ClothoidListMexWrapper('push_back_G1',OBJ,[x0,y0,theta0,x1,y1,theta1]|[CLOT]): "

        if ( nrhs == 8 ) {
          real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
          real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
          real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
          real_type x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" );
          real_type y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" );
          real_type theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" );
          ptr->push_back_G1( x0, y0, theta0, x1, y1, theta1 );
        } else if ( nrhs == 5 ) {
          real_type x1     = getScalarValue( arg_in_2, CMD "Error in reading x1" );
          real_type y1     = getScalarValue( arg_in_3, CMD "Error in reading y1" );
          real_type theta1 = getScalarValue( arg_in_4, CMD "Error in reading theta1" );
          ptr->push_back_G1( x1, y1, theta1 );
        } else {
          MEX_ASSERT( false, CMD "expected 5 or 8 inputs, nrhs = " << nrhs );
        }

        #undef CMD

      } else if ( cmd == "reserve" ) {

        #define CMD "ClothoidListMexWrapper('reserve',OBJ,N): "

        MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0 , CMD "expected no outputs, nlhs = " << nlhs );

        int64_t N = getInt( arg_in_2, CMD "Error in reading N" );
        ptr->reserve( N );

        #undef CMD

      } else if ( cmd == "evaluate" ) {

        #define CMD "ClothoidListMexWrapper('evaluate',OBJ,s): "

        MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs, nrhs = " << nrhs );

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, CMD "Error in reading s" );

        if ( nlhs == 4 ) {
          double * xVals     = createMatrixValue( arg_out_0, size, 1 );
          double * yVals     = createMatrixValue( arg_out_1, size, 1 );
          double * thetaVals = createMatrixValue( arg_out_2, size, 1 );
          double * kappaVals = createMatrixValue( arg_out_3, size, 1 );
          for ( mwSize i=0; i < size; ++i )
            ptr->eval( sVals[i], thetaVals[i], kappaVals[i], xVals[i], yVals[i] );
        } else if ( nlhs == 3 ) {
          double * xVals     = createMatrixValue( arg_out_0, size, 1 );
          double * yVals     = createMatrixValue( arg_out_1, size, 1 );
          double * thetaVals = createMatrixValue( arg_out_2, size, 1 );
          double   dummy;
          for ( mwSize i=0; i < size; ++i )
            ptr->eval( sVals[i], thetaVals[i], dummy, xVals[i], yVals[i] );
        } else if ( nlhs == 2 ) {
          double * xVals = createMatrixValue( arg_out_0, size, 1 );
          double * yVals = createMatrixValue( arg_out_1, size, 1 );
          for ( mwSize i=0; i < size; ++i )
            ptr->eval( sVals[i], xVals[i], yVals[i] );
        } else {
          MEX_ASSERT( false, CMD "expected 2, 3 or 4 outputs, nrhs = " << nrhs );
        }

        #undef CMD

      } else if ( cmd == "eval"    || cmd == "eval_D" ||
                  cmd == "eval_DD" || cmd == "eval_DDD" ) {

        if ( nrhs == 4 ) {

          #define CMD "ClothoidListMexWrapper('eval*',OBJ,s,t): "

          mwSize size, sizet;
          real_type const * s = getVectorPointer( arg_in_2, size,
                                CMD "`s` expected to be a real vector" );
          real_type const * t = getVectorPointer( arg_in_3, sizet,
                                CMD "`t` expected to be a real vector" );

          MEX_ASSERT( size == sizet || size == 1 || sizet ==1,
                      CMD " size(s) = " << size <<
                      " must be equal to size(t) = " << sizet <<
                      " or size(s|t) == 1" );

          mwSize incs = size  == 1 ? 0 : 1;
          mwSize inct = sizet == 1 ? 0 : 1;
          mwSize npts = max(size,sizet);

          #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, pXY += 2 )
          #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, s += incs, t += inct, ++pX, ++pY )

          if ( nlhs == 1 ) {
            real_type *pXY = createMatrixValue( arg_out_0, 2,size );
            if ( cmd == "eval" ) {
              LOOPXY1 ptr->eval( *s, *t, pXY[0], pXY[1] );
            } else if ( cmd == "eval_D" ) {
              LOOPXY1 ptr->eval_D( *s, *t, pXY[0], pXY[1] );
            } else if ( cmd == "eval_DD" ) {
              LOOPXY1 ptr->eval_DD( *s, *t, pXY[0], pXY[1] );
            } else if ( cmd == "eval_DDD" ) {
              LOOPXY1 ptr->eval_DDD( *s, *t, pXY[0], pXY[1] );
            } else {
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else if ( nlhs == 2 ) {
            real_type *pX = createMatrixValue( arg_out_0, 1,size );
            real_type *pY = createMatrixValue( arg_out_1, 1,size );
            if ( cmd == "eval" ) {
              LOOPXY2 ptr->eval( *s, *t, *pX, *pY );
            } else if ( cmd == "eval_D" ) {
              LOOPXY2 ptr->eval_D( *s, *t, *pX, *pY );
            } else if ( cmd == "eval_DD" ) {
              LOOPXY2 ptr->eval_DD( *s, *t, *pX, *pY );
            } else if ( cmd == "eval_DDD" ) {
              LOOPXY2 ptr->eval_DDD( *s, *t, *pX, *pY );
            } else {
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else {
            MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
          }

          #undef CMD
          #undef LOOPXY1
          #undef LOOPXY2

        } else if ( nrhs == 3 ) {

          #define CMD "ClothoidListMexWrapper('eval*',OBJ,s): "

          mwSize npts;
          real_type const * s = getVectorPointer( arg_in_2, npts, CMD "`s` expected to be a real vector" );

          #define LOOPXY1 for ( mwSize i = 0; i < npts; ++i, ++s, pXY += 2 )
          #define LOOPXY2 for ( mwSize i = 0; i < npts; ++i, ++s, ++pX, ++pY )

          if ( nlhs == 1 ) {
            real_type *pXY = createMatrixValue( arg_out_0, 2, npts );
            if ( cmd == "eval" ) {
              LOOPXY1 ptr->eval( *s, pXY[0], pXY[1] );
            } else if ( cmd == "eval_D" ) {
              LOOPXY1 ptr->eval_D( *s, pXY[0], pXY[1] );
            } else if ( cmd == "eval_DD" ) {
              LOOPXY1 ptr->eval_DD( *s, pXY[0], pXY[1] );
            } else if ( cmd == "eval_DDD" ) {
              LOOPXY1 ptr->eval_DDD( *s, pXY[0], pXY[1] );
            } else {
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else if ( nlhs == 2 ) {
            real_type *pX = createMatrixValue( arg_out_0, 1, npts );
            real_type *pY = createMatrixValue( arg_out_1, 1, npts );
            if ( cmd == "eval" ) {
              LOOPXY2 ptr->eval( *s, *pX, *pY );
            } else if ( cmd == "eval_D" ) {
              LOOPXY2 ptr->eval_D( *s, *pX, *pY );
            } else if ( cmd == "eval_DD" ) {
              LOOPXY2 ptr->eval_DD( *s, *pX, *pY );
            } else if ( cmd == "eval_DDD" ) {
              LOOPXY2 ptr->eval_DDD( *s, *pX, *pY );
            } else {
              MEX_ASSERT(false, CMD "Unknown command: " << cmd );
            }
          } else {
            MEX_ASSERT( nlhs == 0, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
          }
          #undef CMD
          #undef LOOPXY1
          #undef LOOPXY2

        } else {
          MEX_ASSERT(false, "ClothoidListMexWrapper('eval*',OBJ,...) bad number of arguments, nrhs = " << nrhs );
        }

      } else if ( cmd == "getSTK" ) {

        #define CMD "ClothoidListMexWrapper('getSTK',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = " << nlhs );

        int_type n = ptr->numSegment();

        double * s     = createMatrixValue( arg_out_0, 1, n+1 );
        double * theta = createMatrixValue( arg_out_1, 1, n+1 );
        double * kappa = createMatrixValue( arg_out_2, 1, n+1 );

        ptr->getSTK( s, theta, kappa );

        #undef CMD

      } else if ( cmd == "getXY" ) {

        #define CMD "ClothoidListMexWrapper('getXY',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );

        int_type n = ptr->numSegment();
        double * x = createMatrixValue( arg_out_0, 1, n+1 );
        double * y = createMatrixValue( arg_out_1, 1, n+1 );

        ptr->getXY( x, y );

        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "ClothoidListMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output, nlhs = " << nlhs );
          mwSize nrx, ncx, nry, ncy;
          real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" );
          real_type const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" );
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

      } else if ( cmd == "sBegin"     || cmd == "sEnd"     ||
                  cmd == "xBegin"     || cmd == "xEnd"     ||
                  cmd == "yBegin"     || cmd == "yEnd"     ||
                  cmd == "thetaBegin" || cmd == "thetaEnd" ||
                  cmd == "kappaBegin" || cmd == "kappaEnd" ||
                  cmd == "length" ) {


        #define CMD "ClothoidListMexWrapper('...',OBJ[,n]): "
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        if ( nrhs == 3 ) {
          int64_t n = getInt( arg_in_2, CMD "Error in reading n" );
          MEX_ASSERT( n > 0 && n <= ptr->numSegment(),
                      CMD "n =  " << n << " must be >= 1 and <= " << ptr->numSegment() );
          --n;
          if      ( cmd == "sBegin"     ) setScalarValue(arg_out_0, ptr->sBegin(n));
          else if ( cmd == "sEnd"       ) setScalarValue(arg_out_0, ptr->sEnd(n));
          else if ( cmd == "xBegin"     ) setScalarValue(arg_out_0, ptr->xBegin(n));
          else if ( cmd == "xEnd"       ) setScalarValue(arg_out_0, ptr->xEnd(n));
          else if ( cmd == "yBegin"     ) setScalarValue(arg_out_0, ptr->yBegin(n));
          else if ( cmd == "yEnd"       ) setScalarValue(arg_out_0, ptr->yEnd(n));
          else if ( cmd == "thetaBegin" ) setScalarValue(arg_out_0, ptr->thetaBegin(n));
          else if ( cmd == "thetaEnd"   ) setScalarValue(arg_out_0, ptr->thetaEnd(n));
          else if ( cmd == "kappaBegin" ) setScalarValue(arg_out_0, ptr->kappaBegin(n));
          else if ( cmd == "kappaEnd"   ) setScalarValue(arg_out_0, ptr->kappaEnd(n));
          else if ( cmd == "length"     ) setScalarValue(arg_out_0, ptr->length(n));
        } else {
          MEX_ASSERT(nrhs == 2, CMD "expected 2 or 3 inputs, nrhs = " << nrhs );
          if      ( cmd == "sBegin"     ) setScalarValue(arg_out_0, ptr->sBegin());
          else if ( cmd == "sEnd"       ) setScalarValue(arg_out_0, ptr->sEnd());
          else if ( cmd == "xBegin"     ) setScalarValue(arg_out_0, ptr->xBegin());
          else if ( cmd == "xEnd"       ) setScalarValue(arg_out_0, ptr->xEnd());
          else if ( cmd == "yBegin"     ) setScalarValue(arg_out_0, ptr->yBegin());
          else if ( cmd == "yEnd"       ) setScalarValue(arg_out_0, ptr->yEnd());
          else if ( cmd == "thetaBegin" ) setScalarValue(arg_out_0, ptr->thetaBegin());
          else if ( cmd == "thetaEnd"   ) setScalarValue(arg_out_0, ptr->thetaEnd());
          else if ( cmd == "kappaBegin" ) setScalarValue(arg_out_0, ptr->kappaBegin());
          else if ( cmd == "kappaEnd"   ) setScalarValue(arg_out_0, ptr->kappaEnd());
          else if ( cmd == "length"     ) setScalarValue(arg_out_0, ptr->totalLength());
        }

        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "ClothoidListMexWrapper('rotate',OBJ,angle,cx,cy): "

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );

        real_type angle = getScalarValue( arg_in_2, CMD "Error in reading angle" );
        real_type cx    = getScalarValue( arg_in_3, CMD "Error in reading cx" );
        real_type cy    = getScalarValue( arg_in_4, CMD "Error in reading cy" );
        ptr->rotate( angle, cx, cy );

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "ClothoidListMexWrapper('translate',OBJ,tx,ty): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        real_type tx = getScalarValue( arg_in_2, CMD "Error in reading tx" );
        real_type ty = getScalarValue( arg_in_3, CMD "Error in reading ty" );
        ptr->translate( tx, ty );

        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "ClothoidListMexWrapper('changeOrigin',OBJ,newX0,newY0): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        real_type newX0 = getScalarValue( arg_in_2, CMD "Error in reading newX0" );
        real_type newY0 = getScalarValue( arg_in_3, CMD "Error in reading newY0" );
        ptr->changeOrigin( newX0, newY0 );

        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "ClothoidListMexWrapper('scale',OBJ,s): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        real_type s = getScalarValue( arg_in_2, CMD "Error in reading s" );
        ptr->scale( s );

        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "ClothoidListMexWrapper('reverse',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        ptr->reverse();

        #undef CMD

      } else if ( cmd == "closestPoint" ) {

        #define CMD "ClothoidListMexWrapper('closestPoint',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 4, CMD "expected 4 outputs, nlhs = " << nlhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output, nlhs = " << nlhs );
          mwSize nrx, ncx, nry, ncy;
          real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" );
          real_type const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" );
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          real_type * X   = createMatrixValue( arg_out_0, nrx, ncx );
          real_type * Y   = createMatrixValue( arg_out_1, nrx, ncx );
          real_type * S   = createMatrixValue( arg_out_2, nrx, ncx );
          real_type * dst = createMatrixValue( arg_out_3, nrx, ncx );

          mwSize size = nrx*ncx;
          for ( mwSize i = 0; i < size; ++i )
            *dst++ = ptr->closestPoint( *x++, *y++, *X++, *Y++, *S++ );
        }
        #undef CMD

      } else if ( cmd == "intersect_line"   ||
                  cmd == "intersect_circle" ||
                  cmd == "intersect_clothoid" ||
                  cmd == "intersect_clothoid_list" ) {

        #define CMD "ClothoidListMexWrapper('intersect_*',OBJ,OBJ2): "
        MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );

        std::vector<real_type> s1, s2;
        int_type               max_iter  = 10;
        real_type              tolerance = 1e-8;

        if ( cmd == "intersect_line" ) {
          LineSegment const * ptr1 = convertMat2Ptr<LineSegment>(arg_in_2);
          ptr->intersect( *ptr1, s1, s2, max_iter, tolerance );
        } else if (  cmd == "intersect_circle" ) {
          CircleArc const * ptr1 = convertMat2Ptr<CircleArc>(arg_in_2);
          ptr->intersect( *ptr1, s1, s2, max_iter, tolerance );
        } else if (  cmd == "intersect_clothoid" ) {
          ClothoidCurve const * ptr1 = convertMat2Ptr<ClothoidCurve>(arg_in_2);
          ptr->intersect( *ptr1, s1, s2, max_iter, tolerance );
        } else if ( cmd == "intersect_clothoid_list" ) {
          ClothoidList const * ptr1 = convertMat2Ptr<ClothoidList>(arg_in_2);
          ptr->intersect( *ptr1, s1, s2, max_iter, tolerance );
        }

        real_type * S1 = createMatrixValue( arg_out_0, s1.size(), 1 );
        real_type * S2 = createMatrixValue( arg_out_1, s2.size(), 1 );

        std::copy( s1.begin(), s1.end(), S1 );
        std::copy( s2.begin(), s2.end(), S2 );

        #undef CMD

      } else if ( cmd == "bbox" ) {

        #define CMD "ClothoidListMexWrapper('bbox', OBJ, max_angle, max_size [,offs]): "

        MEX_ASSERT(nrhs == 4 || nrhs == 5, CMD "expected 4 or 5 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlsh = " << nlhs );

        real_type max_angle = getScalarValue( arg_in_2, CMD "Error in reading max_angle" );
        real_type max_size  = getScalarValue( arg_in_3, CMD "Error in reading max_size" );
        real_type offs      = 0;
        if ( nrhs == 5 ) offs = getScalarValue( arg_in_4, CMD "Error in reading offs" );

        vector<ClothoidCurve::bbData> bb;
        ptr->bbSplit( max_angle, max_size, offs, bb );

        plhs[0] = mxCreateDoubleMatrix(6, bb.size(), mxREAL);
        double * pT = mxGetPr(plhs[0]);
        for ( int i = 0; i < bb.size(); ++i ) {
          Triangle2D const & t = bb[i].t;
          *pT++ = t.x1(); *pT++ = t.y1();
          *pT++ = t.x2(); *pT++ = t.y2();
          *pT++ = t.x3(); *pT++ = t.y3();
        }
        #undef CMD

      } else if ( cmd == "get" ) {

        #define CMD "ClothoidListMexWrapper('get', OBJ, n): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 6, CMD "expected 6 output, nlhs = " << nlhs );

        int64_t n = getInt( arg_in_2, CMD "Error in reading n" );

        MEX_ASSERT( n > 0 && n <= ptr->numSegment(),
                    CMD "n = " << n << " must be >= 1 and <= " << ptr->numSegment() );

        --n;
        ClothoidCurve const & c = ptr->get(n);

        setScalarValue(arg_out_0, c.xBegin());
        setScalarValue(arg_out_1, c.yBegin());
        setScalarValue(arg_out_2, c.thetaBegin());
        setScalarValue(arg_out_3, c.kappaBegin());
        setScalarValue(arg_out_4, c.dkappa());
        setScalarValue(arg_out_5, c.length());

        #undef CMD

      } else if ( cmd == "numSegment" ) {

        #define CMD "ClothoidListMexWrapper('numSegment', OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

        setScalarInt( arg_out_0, ptr->numSegment() );

        #undef CMD

      } else if ( cmd == "build_3arcG2" || cmd == "build_2arcG2" || cmd == "build_CLC" ) {

        #define CMD "ClothoidListMexWrapper('build_[2|3]arcG2', OBJ, ...): "

        MEX_ASSERT(nrhs == 10, CMD "expected 10 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

        real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading `x0`" );
        real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading `y0`" );
        real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading `theta0`" );
        real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading `kappa0`" );
        real_type x1     = getScalarValue( arg_in_6, CMD "Error in reading `x1`" );
        real_type y1     = getScalarValue( arg_in_7, CMD "Error in reading `y1`" );
        real_type theta1 = getScalarValue( arg_in_8, CMD "Error in reading `theta1`" );
        real_type kappa1 = getScalarValue( arg_in_9, CMD "Error in reading `kappa1`" );

        int iter;
        if ( cmd == "build_3arcG2" ) {
          static G2solve3arc g2sol;
          iter = g2sol.build( x0, y0, theta0, kappa0,
                              x1, y1, theta1, kappa1 );
          if ( iter >= 0 ) {
            ptr->init();
            ptr->reserve(3);
            ptr->push_back(g2sol.getS0());
            ptr->push_back(g2sol.getSM());
            ptr->push_back(g2sol.getS1());
          }
        } else if ( cmd == "build_2arcG2" ) {

          static G2solve2arc g2sol;
          iter = g2sol.build( x0, y0, theta0, kappa0,
                              x1, y1, theta1, kappa1 );
          if ( iter >= 0 ) {
            ptr->init();
            ptr->reserve(2);
            ptr->push_back(g2sol.getS0());
            ptr->push_back(g2sol.getS1());
          }
        } else {
          static G2solveCLC g2sol;
          iter = g2sol.build( x0, y0, theta0, kappa0,
                              x1, y1, theta1, kappa1 );
          if ( iter >= 0 ) {
            ptr->init();
            ptr->reserve(3);
            ptr->push_back(g2sol.getS0());
            ptr->push_back(g2sol.getSM());
            ptr->push_back(g2sol.getS1());
          }
        }
        setScalarInt( arg_out_0, iter );

        #undef CMD

      } else if ( cmd == "build_3arcG2fixed" ) {

        #define CMD "ClothoidListMexWrapper('build_3arcG2fixed', OBJ, ...): "

        MEX_ASSERT(nrhs == 12, CMD "expected 12 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

        real_type s0     = getScalarValue( arg_in_2,  CMD "Error in reading `s0`" );
        real_type x0     = getScalarValue( arg_in_3,  CMD "Error in reading `x0`" );
        real_type y0     = getScalarValue( arg_in_4,  CMD "Error in reading `y0`" );
        real_type theta0 = getScalarValue( arg_in_5,  CMD "Error in reading `theta0`" );
        real_type kappa0 = getScalarValue( arg_in_6,  CMD "Error in reading `kappa0`" );
        real_type s1     = getScalarValue( arg_in_7,  CMD "Error in reading `s1`" );
        real_type x1     = getScalarValue( arg_in_8,  CMD "Error in reading `x1`" );
        real_type y1     = getScalarValue( arg_in_9,  CMD "Error in reading `y1`" );
        real_type theta1 = getScalarValue( arg_in_10, CMD "Error in reading `theta1`" );
        real_type kappa1 = getScalarValue( arg_in_11, CMD "Error in reading `kappa1`" );

        static G2solve3arc g2sol;
        int iter = g2sol.build_fixed_length( s0, x0, y0, theta0, kappa0,
                                             s1, x1, y1, theta1, kappa1 );
        if ( iter >= 0 ) {
          ptr->init();
          ptr->reserve(3);
          ptr->push_back(g2sol.getS0());
          ptr->push_back(g2sol.getSM());
          ptr->push_back(g2sol.getS1());
        }

        setScalarInt( arg_out_0, iter );

        #undef CMD

      } else if ( cmd == "build_G1" ) {

        #define CMD "ClothoidListMexWrapper('build_G1', OBJ, x, y [, theta]): "

        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs);

        bool ok = true;

        if ( nrhs == 4 ) {
          mwSize nx, ny;
          real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
          real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

          MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " != length(y) = " << ny );

          ok = ptr->build_G1( nx, x, y );

        } else if ( nrhs == 5 ) {

          mwSize nx, ny, nt;
          real_type const * x     = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
          real_type const * y     = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );
          real_type const * theta = getVectorPointer( arg_in_4, nt, CMD "Error in reading theta" );

          MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " != length(y) = " << ny );
          MEX_ASSERT( nx == nt, CMD "length(theta) = " << nt << " != length(x) = length(y) = " << ny );

          ok = ptr->build_G1( nx, x, y, theta );

        } else {
          MEX_ASSERT( false, CMD "expected 4 or 5 input, nrhs = " << nrhs);
        }

        setScalarBool( arg_out_0, ok );

        #undef CMD

      } else if ( cmd == "build_theta" ) {

        #define CMD "ClothoidListMexWrapper('build_theta', OBJ, x, y): "

        MEX_ASSERT(nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs);
        MEX_ASSERT(nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs);

        bool ok = true;

        mwSize nx, ny;
        real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
        real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

        MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " != length(y) = " << ny );

        real_type * theta = createMatrixValue( arg_out_0, nx, 1 );

        ok = ptr->build_theta( nx, x, y, theta );

        setScalarBool( arg_out_1, ok );

        #undef CMD

      } else if ( cmd == "deltaTheta" ) {

        #define CMD "ClothoidListMexWrapper('deltaTheta', OBJ): "

        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs);
        MEX_ASSERT(nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs);

        int_type nseg = ptr->numSegment();

        real_type * dtheta = createMatrixValue( arg_out_0, nseg, 1 );
        ptr->getDeltaTheta( dtheta );

        #undef CMD

      } else if ( cmd == "deltaKappa" ) {

        #define CMD "ClothoidListMexWrapper('deltaKappa', OBJ): "

        MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs);
        MEX_ASSERT(nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs);

        int_type nseg = ptr->numSegment();

        real_type * dkappa = createMatrixValue( arg_out_0, nseg, 1 );
        ptr->getDeltaKappa( dkappa );

        #undef CMD

      } else if ( cmd == "findST" ) {

        #define CMD "ClothoidListMexWrapper('findST',OBJ,x,y[,nseg]): "
        MEX_ASSERT( nrhs == 4 || nrhs == 5,
                    CMD "expected 4 or 5 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 3,
                    CMD "expected 3 output, nlhs = " << nlhs );
        mwSize nrx, ncx, nry, ncy;
        real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx,
                              CMD "`x` expected to be a real vector/matrix" );
        real_type const * y = getMatrixPointer( arg_in_3, nry, ncy,
                              CMD "`y` expected to be a real vector/matrix" );
        MEX_ASSERT( nrx == nry && ncx == ncy,
                    CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                    nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

        real_type * s   = createMatrixValue( arg_out_0, nrx, ncx );
        real_type * t   = createMatrixValue( arg_out_1, nrx, ncx );
        real_type * idx = createMatrixValue( arg_out_2, nrx, ncx );

        mwSize size = nrx*ncx;
        if ( nrhs == 4 ) {
          for ( mwSize i = 0; i < size; ++i )
            int_type nseg = ptr->findST( *x++, *y++, *s++, *t++ );
            *idx++ = nseg >= 0 ? nseg+1 : nseg;
        } else {
          real_type const * idx_guess = getMatrixPointer( arg_in_3, nrx, ncx,
                                        CMD "`nseg` expected to be a real vector/matrix" );
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` and `nseg` expected to be of the same size, found size(x) =  size(y) = " <<
                      nrx << " x " << nry << " size(nseg) = " << nry << " x " << ncy );
          for ( mwSize i = 0; i < size; ++i ) {
            int_type nseg = int_type(*idx_guess++);
            bool ok = ptr->findST( nseg-1, *x++, *y++, *s++, *t++ );
            *idx = ok ? nseg+1 : -(1+nseg);
            ++idx; ++idx_guess;
          }
        }

        #undef CMD
      } else if ( cmd == "export_table" || cmd == "export_ruby" ) {

        #define CMD "ClothoidListMexWrapper('export_[table|ruby]', OBJ, filename ): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nlhs = " << nlhs);
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nrhs = " << nrhs);

        MEX_ASSERT( mxIsChar(arg_in_2), CMD "filename must be a string" );
        string filename = mxArrayToString(arg_in_2);

        std::ofstream file(filename.c_str());

        MEX_ASSERT( file.good(), CMD " cannot open file: `" << filename << "`" );

        if ( cmd == "export_table" ) ptr->export_table(file);
        else                         ptr->export_ruby(file);

        file.close();

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "ClothoidListMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs);
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs);

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        #undef CMD

      } else if ( cmd == "copy" ) {

        #define CMD "ClothoidListMexWrapper('copy',OBJ,OBJ1): "
        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        ClothoidList const * CL = convertMat2Ptr<ClothoidList>(arg_in_2);
        ptr->copy(*CL);

        #undef CMD

        // Warn if other commands were ignored
      } else {
        MEX_ASSERT(false, "Unknown command: " << cmd );
      }

    } catch ( std::exception const & e ) {
    	mexErrMsgTxt(e.what());

    } catch (...) {
  	  mexErrMsgTxt("clothoid failed\n");
    }
  }
}
