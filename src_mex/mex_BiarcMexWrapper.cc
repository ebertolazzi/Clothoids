/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Biarc.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"BiarcMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = BiarcMexWrapper( 'new' );\n" \
"\n" \
"    On output:\n" \
"       OBJ   = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    BiarcMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  - Build:\n" \
"    BiarcMexWrapper( 'build', OBJ, x0, y0, theta0, x1, y1, theta1 );\n" \
"    BiarcMexWrapper( 'build_3P', OBJ, x0, y0, x1, y1, x2, y2 );\n" \
"    [arc0,arc1] = BiarcMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
"  - Eval:\n" \
"    [x,y,theta,kappa] = BiarcMexWrapper( 'evaluate', OBJ, ss );\n" \
"    [x0,y0,theta0,kappa0,L0,x1,y1,theta1,kappa1,L1] = BiarcMexWrapper( 'getPars', OBJ );\n" \
"\n" \
"    [x,y]         = BiarcMexWrapper( 'eval', OBJ, s[, t] );\n" \
"    [x_D,y_D]     = BiarcMexWrapper( 'eval_D', OBJ, s[, t] );\n" \
"    [x_DD,y_DD]   = BiarcMexWrapper( 'eval_DD', OBJ, s[, t] );\n" \
"    [x_DDD,y_DDD] = BiarcMexWrapper( 'eval_DDD', OBJ, s[, t] );\n" \
"\n" \
"  res = CircleMexWrapper( 'xBegin0', OBJ );\n" \
"  res = CircleMexWrapper( 'xEnd0', OBJ );\n" \
"  res = CircleMexWrapper( 'yBegin0', OBJ );\n" \
"  res = CircleMexWrapper( 'yEnd0', OBJ );\n" \
"  res = CircleMexWrapper( 'thetaBegin0', OBJ );\n" \
"  res = CircleMexWrapper( 'thetaEnd0', OBJ );\n" \
"  res = CircleMexWrapper( 'kappa0', OBJ );\n" \
"  res = CircleMexWrapper( 'length0', OBJ );\n" \
"  res = CircleMexWrapper( 'xBegin1', OBJ );\n" \
"  res = CircleMexWrapper( 'xEnd1', OBJ );\n" \
"  res = CircleMexWrapper( 'yBegin1', OBJ );\n" \
"  res = CircleMexWrapper( 'yEnd1', OBJ );\n" \
"  res = CircleMexWrapper( 'thetaBegin1', OBJ );\n" \
"  res = CircleMexWrapper( 'thetaEnd1', OBJ );\n" \
"  res = CircleMexWrapper( 'kappa1', OBJ );\n" \
"  res = CircleMexWrapper( 'length1', OBJ );\n" \
"\n" \
"  - Transform:\n" \
"    BiarcMexWrapper( 'changeOrigin', OBJ, newX0, newY0 );\n" \
"    BiarcMexWrapper( 'rotate', OBJ, angle, cx, cy );\n" \
"    BiarcMexWrapper( 'translate', OBJ, tx, ty );\n" \
"    BiarcMexWrapper( 'scale', OBJ, scaling );\n" \
"    BiarcMexWrapper( 'reverse', OBJ );\n" \
"  - Distance:\n" \
"    [X,Y,s,dst] = BiarcMexWrapper( 'closestPoint', OBJ, x, y );\n" \
"    [dst,s]     = BiarcMexWrapper( 'distance', OBJ, x, y );\n" \
"    [X,Y,s,dst] = BiarcMexWrapper( 'closestPointBySample', OBJ, x, y, ds );\n" \
"    [s,t]       = BiarcMexWrapper( 'findST', OBJ, x, y );\n" \
"\n" \
"=====================================================================================\n" \
"\n" \
"Autors: Enrico Bertolazzi^(1), Marco Frego^(2)\n" \
"  (1) Department of Industrial Engineering\n" \
"  (2) Department of Information Engineering and Computer Science\n" \
"  University of Trento\n" \
"  enrico.bertolazzi@unitn.it\n" \
"  m.fregox@gmail.com\n" \
"\n" \
"=====================================================================================\n"

namespace G2lib {

  using namespace std;

  static
  Biarc *
  DATA_NEW( mxArray * & mx_id ) {
    Biarc * ptr = new Biarc();
    mx_id = convertPtr2Mat<Biarc>(ptr);
    return ptr;
  }

  static
  inline
  Biarc *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<Biarc>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<Biarc>(mx_id);
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

      MEX_ASSERT( mxIsChar(arg_in_0), "BiarcMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      Biarc * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "BiarcMexWrapper, expected 1 output, nlhs = " << nlhs );
        MEX_ASSERT( nrhs == 1, "BiarcMexWrapper, expected 1 input, nrhs = " << nrhs );

      } else if ( cmd == "build" ) {

        #define CMD "BiarcMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1): "

        MEX_ASSERT( nrhs == 8 , CMD "expected 8 inputs, nrhs = " << nrhs );

        real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
        real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
        real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
        real_type x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" );
        real_type y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" );
        real_type theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" );

        bool ok = ptr->build( x0, y0, theta0, x1, y1, theta1 );

        // returns the status of the interpolation
        setScalarBool( arg_out_0, ok );

        #undef CMD

      } else if ( cmd == "build_3P" ) {

        if ( nrhs == 8 ) {
          #define CMD "BiarcMexWrapper('build_3P',OBJ,x0,y0,x1,y1,x2,y2): "

          real_type x0 = getScalarValue( arg_in_2, CMD "Error in reading x0" );
          real_type y0 = getScalarValue( arg_in_3, CMD "Error in reading y0" );
          real_type x1 = getScalarValue( arg_in_4, CMD "Error in reading x1" );
          real_type y1 = getScalarValue( arg_in_5, CMD "Error in reading y1" );
          real_type x2 = getScalarValue( arg_in_6, CMD "Error in reading x2" );
          real_type y2 = getScalarValue( arg_in_7, CMD "Error in reading y2" );

          bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 );

          // returns the status of the interpolation
          setScalarBool(arg_out_0,ok);

          #undef CMD
        } else if ( nrhs == 5 ) {
          #define CMD "BiarcMexWrapper('build_3P',OBJ,p0,p1,p2): "

          mwSize n;
          real_type const * p0 = getVectorPointer( arg_in_2, n, CMD "Error in reading p0" );
          MEX_ASSERT( n == 2 , CMD "Error in reading length(p0) == " << n << " expect length(p0) == 2" );
          real_type const * p1 = getVectorPointer( arg_in_3, n, CMD "Error in reading p1" );
          MEX_ASSERT( n == 2 , CMD "Error in reading length(p1) == " << n << " expect length(p1) == 2" );
          real_type const * p2 = getVectorPointer( arg_in_4, n, CMD "Error in reading p2" );
          MEX_ASSERT( n == 2 , CMD "Error in reading length(p2) == " << n << " expect length(p2) == 2" );

          bool ok = ptr->build_3P( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );

          // returns the status of the interpolation
          setScalarBool(arg_out_0,ok);

          #undef CMD
        } else {
          MEX_ASSERT( false, "BiarcMexWrapper('build_3P',OBJ,...) expected 5 or 8 arguments");
        }

      } else if ( cmd == "evaluate" ) {

        #define CMD "BiarcMexWrapper('evaluate',OBJ,s): "

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
        } else if ( nlhs == 2 ) {
          double * xVals = createMatrixValue( arg_out_0, size, 1 );
          double * yVals = createMatrixValue( arg_out_1, size, 1 );
          for ( mwSize i=0; i < size; ++i )
            ptr->eval( sVals[i], xVals[i], yVals[i] );
        } else {
          MEX_ASSERT( false, CMD "expected 2 or 4 outputs, nlhs = " << nlhs );
        }

        #undef CMD

      } else if ( cmd == "eval"    || cmd == "eval_D" ||
                  cmd == "eval_DD" || cmd == "eval_DDD" ) {

        if ( nrhs == 4 ) {

          #define CMD "BiarcMexWrapper('eval*',OBJ,s,t): "

          mwSize size, sizet;
          real_type const * s = getVectorPointer( arg_in_2, size,  CMD "`s` expected to be a real vector" );
          real_type const * t = getVectorPointer( arg_in_3, sizet, CMD "`t` expected to be a real vector" );

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
          #undef LOOPXY1
          #undef LOOPXY2

          #undef CMD

        } else if ( nrhs == 3 ) {

          #define CMD "BiarcMexWrapper('eval*',OBJ,s): "

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
          #undef LOOPXY1
          #undef LOOPXY2

        } else {
          MEX_ASSERT(false, CMD " bad number of arguments nrhs = " << nrhs );
        }

        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "BiarcMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 outputs, nlhs = " << nlhs );
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

      } else if ( cmd == "closestPoint" ) {

        #define CMD "BiarcMexWrapper('closestPoint',OBJ,x,y): "
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

      } else if ( cmd == "findST" ) {

        #define CMD "BiarcMexWrapper('findST',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );
        mwSize nrx, ncx, nry, ncy;
        real_type const * x = getMatrixPointer( arg_in_2, nrx, ncx,
                              CMD "`x` expected to be a real vector/matrix" );
        real_type const * y = getMatrixPointer( arg_in_3, nry, ncy,
                              CMD "`y` expected to be a real vector/matrix" );
        MEX_ASSERT( nrx == nry && ncx == ncy,
                    CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                    nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

        real_type * s = createMatrixValue( arg_out_0, nrx, ncx );
        real_type * t = createMatrixValue( arg_out_1, nrx, ncx );

        mwSize size = nrx*ncx;
        for ( mwSize i = 0; i < size; ++i )
          ptr->findST( *x++, *y++, *s++, *t++ );

        #undef CMD

      } else if ( cmd == "getPars" ) {

        MEX_ASSERT(nrhs == 2, "BiarcMexWrapper('getPars',OBJ): expected 2 inputs, nrhs = " << nrhs );
        if ( nlhs > 0 ) setScalarValue(arg_out_0,ptr->xBegin0());
        if ( nlhs > 1 ) setScalarValue(arg_out_1,ptr->yBegin0());
        if ( nlhs > 2 ) setScalarValue(arg_out_2,ptr->thetaBegin0());
        if ( nlhs > 3 ) setScalarValue(arg_out_3,ptr->kappa0());
        if ( nlhs > 4 ) setScalarValue(arg_out_4,ptr->length0());

        if ( nlhs > 5 ) setScalarValue(arg_out_5,ptr->xBegin1());
        if ( nlhs > 6 ) setScalarValue(arg_out_6,ptr->yBegin1());
        if ( nlhs > 7 ) setScalarValue(arg_out_7,ptr->thetaBegin1());
        if ( nlhs > 8 ) setScalarValue(arg_out_8,ptr->kappa1());
        if ( nlhs > 9 ) setScalarValue(arg_out_9,ptr->length1());

      } else if ( cmd == "xBegin0" ) {

        #define CMD "BiarcMexWrapper('xBegin0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->xBegin0());

        #undef CMD

      } else if ( cmd == "xEnd0" ) {

        #define CMD "BiarcMexWrapper('xEnd0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->xEnd0());

        #undef CMD

      } else if ( cmd == "xBegin1" ) {

        #define CMD "BiarcMexWrapper('xBegin1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->xBegin1());

        #undef CMD

      } else if ( cmd == "xEnd1" ) {

        #define CMD "BiarcMexWrapper('xEnd1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->xEnd1());

        #undef CMD

      } else if ( cmd == "yBegin0" ) {

        #define CMD "BiarcMexWrapper('yBegin0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->yBegin0());

        #undef CMD

      } else if ( cmd == "yEnd0" ) {

        #define CMD "BiarcMexWrapper('yEnd0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->yEnd0());

        #undef CMD

      } else if ( cmd == "yBegin1" ) {

        #define CMD "BiarcMexWrapper('yBegin1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->yBegin1());

        #undef CMD

      } else if ( cmd == "yEnd1" ) {

        #define CMD "BiarcMexWrapper('yEnd1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->yEnd1());

        #undef CMD

      } else if ( cmd == "thetaBegin0" ) {

        #define CMD "BiarcMexWrapper('thetaBegin0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->thetaBegin0());

        #undef CMD

      } else if ( cmd == "thetaEnd0" ) {

        #define CMD "BiarcMexWrapper('thetaEnd0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->thetaEnd0());

        #undef CMD

      } else if ( cmd == "thetaBegin1" ) {

        #define CMD "BiarcMexWrapper('thetaBegin1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->thetaBegin1());

        #undef CMD

      } else if ( cmd == "thetaEnd1" ) {

        #define CMD "BiarcMexWrapper('thetaEnd1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->thetaEnd1());

        #undef CMD

      } else if ( cmd == "kappa0" ) {

        #define CMD "BiarcMexWrapper('kappa0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->kappa0());

        #undef CMD

      } else if ( cmd == "kappa1" ) {

        #define CMD "BiarcMexWrapper('kappa1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->kappa1());

        #undef CMD

      } else if ( cmd == "length0" ) {

        #define CMD "BiarcMexWrapper('length0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->length0());

        #undef CMD

      } else if ( cmd == "length1" ) {

        #define CMD "BiarcMexWrapper('length1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->length1());

        #undef CMD

      } else if ( cmd == "length" ) {

        #define CMD "BiarcMexWrapper('length',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );
        setScalarValue(arg_out_0, ptr->length());

        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "BiarcMexWrapper('rotate',OBJ,angle,cx,cy): "

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs, nrhs = " << nrhs );

        real_type angle = getScalarValue( arg_in_2, CMD "Error in reading angle" );
        real_type cx    = getScalarValue( arg_in_3, CMD "Error in reading cx" );
        real_type cy    = getScalarValue( arg_in_4, CMD "Error in reading cy" );
        ptr->rotate(angle, cx, cy);

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "BiarcMexWrapper('translate',OBJ,tx,ty): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        real_type tx = getScalarValue( arg_in_2, CMD "Error in reading tx" );
        real_type ty = getScalarValue( arg_in_3, CMD "Error in reading ty" );
        ptr->translate(tx, ty);

        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "BiarcMexWrapper('changeOrigin',OBJ,newX0,newY0): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        real_type newX0 = getScalarValue( arg_in_2, CMD "Error in reading newX0" );
        real_type newY0 = getScalarValue( arg_in_3, CMD "Error in reading newY0" );
        ptr->changeOrigin(newX0, newY0);

        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "BiarcMexWrapper('scale',OBJ,s): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        real_type s = getScalarValue( arg_in_2, CMD "Error in reading s" );
        ptr->scale(s);

        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "BiarcMexWrapper('reverse',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        ptr->reverse();

        #undef CMD

      } else if ( cmd == "to_nurbs" ) {

        #define CMD "BiarcMexWrapper('to_nurbs',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );

        CircleArc const & C0 = ptr->getC0();
        CircleArc const & C1 = ptr->getC1();

        int_type npts0, nknots0, npts1, nknots1;
        C0.paramNURBS( nknots0, npts0 );
        C1.paramNURBS( nknots1, npts1 );

        mxArray * mx_knots0, * mx_Poly0, * mx_knots1, * mx_Poly1;

        double * knots0 = createMatrixValue( mx_knots0, 1, nknots0 );
        double * poly0  = createMatrixValue( mx_Poly0,  3, npts0   );
        double * knots1 = createMatrixValue( mx_knots1, 1, nknots1 );
        double * poly1  = createMatrixValue( mx_Poly1,  3, npts1   );

        C0.toNURBS( knots0, poly0 );
        C1.toNURBS( knots1, poly1 );

        static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" };
        arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
        arg_out_1 = mxCreateStructMatrix(1,1,6,fieldnames);

        mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
        mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
        mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
        mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts0) );
        mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots0 );
        mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly0 );

        mxSetFieldByNumber( arg_out_1, 0, 0, mxCreateString("rB") );
        mxSetFieldByNumber( arg_out_1, 0, 1, mxCreateDoubleScalar(3) );
        mxSetFieldByNumber( arg_out_1, 0, 2, mxCreateDoubleScalar(2) );
        mxSetFieldByNumber( arg_out_1, 0, 3, mxCreateDoubleScalar(npts1) );
        mxSetFieldByNumber( arg_out_1, 0, 4, mx_knots1 );
        mxSetFieldByNumber( arg_out_1, 0, 5, mx_Poly1 );

        #undef CMD

      } else if ( cmd == "info" ) {

        #define CMD "BiarcMexWrapper('info',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = " << nlhs );

        ptr->info(cout);

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "BiarcMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT(nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        #undef CMD

        // Warn if other commands were ignored
      } else {
        MEX_ASSERT(false, "Unknown command: " << cmd );
      }

    } catch ( std::exception const & e ) {
    	mexErrMsgTxt(e.what());

    } catch (...) {
  	  mexErrMsgTxt("Biarc failed\n");
    }
  }
}
