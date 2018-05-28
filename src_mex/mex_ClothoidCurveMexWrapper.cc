/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"ClothoidCurveMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = ClothoidCurveMexWrapper( 'new' ) ;\n" \
"    OBJ = ClothoidCurveMexWrapper( 'new', x0, y0, theta0, k0, dk, L ) ;\n" \
"\n" \
"    On input:\n" \
"      x0, y0 = coordinate of initial point\n" \
"      theta0 = orientation (angle) of the clothoid at initial point" \
"      k0     = curvature of the clothoid at initial point\n" \
"      dk     = derivative of curvature respect to arclength\n" \
"      L      = length of the clothoid curve from initial to final point\n" \
"\n" \
"     On output:\n" \
"       OBJ   = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    ClothoidCurveMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  - Build:\n" \
"    ClothoidCurveMexWrapper( 'build', OBJ, x0, y0, theta0, k0, dk, L ) ;\n" \
"    ClothoidCurveMexWrapper( 'build_G1', OBJ, x0, y0, theta0, x1, y1, theta1 ) ;\n" \
"    res = ClothoidCurveMexWrapper( 'build_forward', OBJ,x0,y0,theta0,k0,x1,y1 ) ;\n" \
"\n" \
"  - Eval:\n" \
"    [x,y,theta,kappa] = ClothoidCurveMexWrapper( 'evaluate', OBJ, ss ) ;\n" \
"    [x0,y0,theta0,k0,dk,smin,smax] = ClothoidCurveMexWrapper( 'getPars', OBJ ) ;\n" \
"\n" \
"    [x,y]         = ClothoidCurveMexWrapper( 'eval', OBJ, ss, offs ) ;\n" \
"    [x_D,y_D]     = ClothoidCurveMexWrapper( 'eval_D', OBJ, ss, offs ) ;\n" \
"    [x_DD,y_DD]   = ClothoidCurveMexWrapper( 'eval_DD', OBJ, ss, offs ) ;\n" \
"    [x_DDD,y_DDD] = ClothoidCurveMexWrapper( 'eval_DDD', OBJ, ss, offs ) ;\n" \
"\n" \
"  - Transform:\n" \
"    ClothoidCurveMexWrapper( 'trim', OBJ, smin, smax ) ;\n" \
"    ClothoidCurveMexWrapper( 'changeOrigin', OBJ, newX0, newY0 ) ;\n" \
"    ClothoidCurveMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"    ClothoidCurveMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"    ClothoidCurveMexWrapper( 'scale', OBJ, scaling ) ;\n" \
"    ClothoidCurveMexWrapper( 'reverse', OBJ ) ;\n" \
"    [xp, yp, xm, ym] = ClothoidCurveMexWrapper( 'infinity', OBJ ) ;\n" \
"\n" \
"  - Distance:\n" \
"    [X,Y,s,dst] = ClothoidCurveMexWrapper( 'closestPoint', OBJ, x, y ) ;\n" \
"    [dst,s] = ClothoidCurveMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"    [X,Y,s,dst] = ClothoidCurveMexWrapper( 'closestPointBySample', OBJ, x, y, ds ) ;\n" \
"    [dst,s] = ClothoidCurveMexWrapper( 'distanceBySample', OBJ, x, y, ds ) ;\n" \
"\n" \
"=====================================================================================\n" \
"\n" \
"Autors: Enrico Bertolazzi^(1), Marco Frego^(2), Paolo Bevilacqua^(2)\n" \
"  (1) Department of Industrial Engineering\n" \
"  (2) Department of Information Engineering and Computer Science\n" \
"  University of Trento\n" \
"  enrico.bertolazzi@unitn.it\n" \
"  m.fregox@gmail.com\n" \
"  paolo.bevilacqua@unitn.it\n" \
"\n" \
"=====================================================================================\n"

namespace G2lib {

  using namespace std;

  static
  ClothoidCurve *
  DATA_NEW( mxArray * & mx_id ) {
    ClothoidCurve * ptr = new ClothoidCurve();
    mx_id = convertPtr2Mat<ClothoidCurve>(ptr);
    return ptr ;
  }

  static
  inline
  ClothoidCurve *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<ClothoidCurve>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<ClothoidCurve>(mx_id);
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

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidCurveMexWrapper: First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;

      bool do_new = cmd == "new" ;
      ClothoidCurve * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build"  ) {

        if ( do_new ) {
          MEX_ASSERT( nlhs == 1, "ClothoidCurveMexWrapper, expected 1 output" );
        } else {
          MEX_ASSERT( nlhs == 0, "ClothoidCurveMexWrapper, expected no output" );
        }

        indexType kk = do_new ? 0 : 1 ;

        valueType x0(0), y0(0), theta0(0), k0(0), dk(0), L(0);

        if ( do_new && nrhs == 1 ) {
          // nothing to do
        } else if ( nrhs == 1+kk ) {

          #define CMD "ClothoidCurveMexWrapper('build',OBJ,struct): "

          MEX_ASSERT( mxIsStruct(arg_in_1), CMD << "Argument n. " << 1+kk << " must be a struct" ) ;

          mxArray * mx_x0     = mxGetField( prhs[1+kk], 0, "x0" ) ;
          mxArray * mx_y0     = mxGetField( prhs[1+kk], 0, "y0" ) ;
          mxArray * mx_theta0 = mxGetField( prhs[1+kk], 0, "theta0" ) ;
          mxArray * mx_k0     = mxGetField( prhs[1+kk], 0, "k0" ) ;
          mxArray * mx_dk     = mxGetField( prhs[1+kk], 0, "dk" ) ;

          MEX_ASSERT( mx_x0     != nullptr, CMD "Field `x0` is missing" );
          MEX_ASSERT( mx_y0     != nullptr, CMD "Field `y0` is missing" );
          MEX_ASSERT( mx_theta0 != nullptr, CMD "Field `theta0` is missing" );
          MEX_ASSERT( mx_k0     != nullptr, CMD "Field `k0` is missing" );
          MEX_ASSERT( mx_dk     != nullptr, CMD "Field `dk` is missing" );

          x0     = getScalarValue( mx_x0,     CMD "Field `x0` must be a real double scalar" ) ;
          y0     = getScalarValue( mx_y0,     CMD "Field `y0` must be a real double scalar" ) ;
          theta0 = getScalarValue( mx_theta0, CMD "Field `theta0` must be a real double scalar" ) ;
          k0     = getScalarValue( mx_k0,     CMD "Field `k0` must be a real double scalar" ) ;
          dk     = getScalarValue( mx_dk,     CMD "Field `dk` must be a real double scalar" ) ;

          #undef CMD

        } else if ( nrhs == 7+kk ) {

          #define CMD "ClothoidCurveMexWrapper('build',OBJ,x0,y0,theta0,k0,dk,L): "

          x0     = getScalarValue(prhs[1+kk],CMD "Error in reading x0") ;
          y0     = getScalarValue(prhs[2+kk],CMD "Error in reading y0") ;
          theta0 = getScalarValue(prhs[3+kk],CMD "Error in reading theta0") ;
          k0     = getScalarValue(prhs[4+kk],CMD "Error in reading k0") ;
          dk     = getScalarValue(prhs[5+kk],CMD "Error in reading dk") ;
          L      = getScalarValue(prhs[6+kk],CMD "Error in reading L") ;
        } else {
          MEX_ASSERT( false, CMD  << "expected " << kk+1 << " or " << kk+7 << " inputs") ;
        }

        ptr->build( x0, y0, theta0, k0, dk, L );

        #undef CMD

      } else if ( cmd == "build_G1" ) {

        #define CMD "ClothoidCurveMexWrapper('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): "

        MEX_ASSERT( nrhs == 8 , CMD "expected 8 inputs") ;

        valueType x0(0), y0(0), theta0(0), x1(0), y1(0), theta1(0);

        x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" ) ;
        y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" ) ;
        theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" ) ;
        x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" ) ;
        y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" ) ;
        theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" ) ;

        ptr->build_G1(x0, y0, theta0, x1, y1, theta1);

        #undef CMD

      } else if ( cmd == "build_forward" ) {

        #define CMD "ClothoidCurveMexWrapper('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): "

        MEX_ASSERT( nrhs == 8 , CMD "expected 8 inputs") ;

        valueType x0(0), y0(0), theta0(0), kappa0(0), x1(0), y1(0);

        x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" ) ;
        y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" ) ;
        theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" ) ;
        kappa0 = getScalarValue( arg_in_5, CMD "Error in reading kappa0" ) ;
        x1     = getScalarValue( arg_in_6, CMD "Error in reading x1" ) ;
        y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" ) ;

        bool res = ptr->build_forward(x0, y0, theta0, kappa0, x1, y1);

        // returns the status of the interpolation
        mwSize dims[2] = {1,1} ;
        arg_out_0 = mxCreateLogicalArray(2, dims);
        ((bool*)mxGetPr(arg_out_0))[0] = res;
        //mexPrintf("Result: %d\n", res);

        #undef CMD

      } else if ( cmd == "evaluate" ) {

        #define CMD "ClothoidCurveMexWrapper('evaluate',OBJ,s): "

        MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, CMD "Error in reading s" );

        if ( nlhs == 4 ) {
          double * xVals     = createMatrixValue( arg_out_0, size, 1 );
          double * yVals     = createMatrixValue( arg_out_1, size, 1 );
          double * thetaVals = createMatrixValue( arg_out_2, size, 1 );
          double * kappaVals = createMatrixValue( arg_out_3, size, 1 );
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval( sVals[i], thetaVals[i], kappaVals[i], xVals[i], yVals[i] );
        } else if ( nlhs == 2 ) {
          double * xVals = createMatrixValue( arg_out_0, size, 1 );
          double * yVals = createMatrixValue( arg_out_1, size, 1 );
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval( sVals[i], xVals[i], yVals[i] );
        } else {
          MEX_ASSERT( false, CMD "expected 2 or 4 outputs") ;
        }

        #undef CMD

      } else if ( cmd == "eval"    || cmd == "eval_D" ||
                  cmd == "eval_DD" || cmd == "eval_DDD" ) {

        #define CMD "ClothoidCurveMexWrapper('eval*',OBJ,s[,offs]): "

        MEX_ASSERT( nrhs == 3 || nrhs == 4, CMD "expected 3 or 4 inputs") ;
        MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, CMD "Error in reading s" );

        double offs = 0 ;
        if ( nrhs == 4 ) offs = getScalarValue( arg_in_3, CMD "Error in reading offs" ) ;

        double * xVals = createMatrixValue( arg_out_0, size, 1 );
        double * yVals = createMatrixValue( arg_out_1, size, 1 );
        if ( cmd == "eval" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval( sVals[i], offs, xVals[i], yVals[i] );
        } else if ( cmd == "eval_D" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_D( sVals[i], offs, xVals[i], yVals[i] );
        } else if ( cmd == "eval_DD" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_DD( sVals[i], offs, xVals[i], yVals[i] );
        } else {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_DDD( sVals[i], offs, xVals[i], yVals[i] );
        }

        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "ClothoidCurveMexWrapper('distance',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input");
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output");
          mwSize nrx, ncx, nry, ncy;
          valueType const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" ) ;
          valueType const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" ) ;
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

      } else if ( cmd == "distanceBySample" ) {

        #define CMD "ClothoidCurveMexWrapper('distance',OBJ,x,y,ds): "
        MEX_ASSERT( nrhs == 5, CMD "expected 5 input");
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output");
          mwSize nrx, ncx, nry, ncy;
          valueType const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" ) ;
          valueType const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" ) ;
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );
          valueType ds = getScalarValue( arg_in_4, CMD "`ds` expected to be a real scalar" );
          MEX_ASSERT( ds > 0, CMD "`ds` = " << ds << " must be a positive number" );

          valueType * dst = createMatrixValue( arg_out_0, nrx, ncx ) ;

          mwSize size = nrx*ncx ;
          if ( nlhs > 1 ) {
            valueType * s = createMatrixValue( arg_out_1, nrx, ncx ) ;
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distanceBySample( ds, *x++, *y++, *s++ ) ;
          } else {
            for ( mwSize i = 0 ; i < size ; ++i )
              *dst++ = ptr->distanceBySample( ds, *x++, *y++ ) ;
          }
        }
        #undef CMD

      } else if ( cmd == "closestPoint" ) {

        #define CMD "ClothoidCurveMexWrapper('closestPoint',OBJ,x,y): "
        MEX_ASSERT( nrhs == 4, CMD "expected 4 input");
        MEX_ASSERT( nlhs == 4, CMD "expected 4 outputs") ;
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output");
          mwSize nrx, ncx, nry, ncy;
          valueType const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" ) ;
          valueType const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" ) ;
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          valueType * X   = createMatrixValue( arg_out_0, nrx, ncx ) ;
          valueType * Y   = createMatrixValue( arg_out_1, nrx, ncx ) ;
          valueType * S   = createMatrixValue( arg_out_2, nrx, ncx ) ;
          valueType * dst = createMatrixValue( arg_out_3, nrx, ncx ) ;

          mwSize size = nrx*ncx ;
          for ( mwSize i = 0 ; i < size ; ++i )
            *dst++ = ptr->closestPoint( *x++, *y++, *X++, *Y++, *S++ ) ;
        }
        #undef CMD

      } else if ( cmd == "closestPointBySample" ) {

        #define CMD "ClothoidCurveMexWrapper('closestPointBySample',OBJ,x,y,ds): "
        MEX_ASSERT( nrhs == 5, CMD "expected 5 input");
        MEX_ASSERT( nlhs == 4, CMD "expected 4 outputs") ;
        if ( nlhs > 0 ) {
          MEX_ASSERT(nlhs <= 2, CMD "expected 1 or 2 output");
          mwSize nrx, ncx, nry, ncy;
          valueType const * x = getMatrixPointer( arg_in_2, nrx, ncx, CMD "`x` expected to be a real vector/matrix" ) ;
          valueType const * y = getMatrixPointer( arg_in_3, nry, ncy, CMD "`y` expected to be a real vector/matrix" ) ;
          MEX_ASSERT( nrx == nry && ncx == ncy,
                      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
                      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy );

          valueType ds = getScalarValue( arg_in_4, CMD "`ds` expected to be a real scalar" );
          MEX_ASSERT( ds > 0, CMD "`ds` = " << ds << " must be a positive number" );

          valueType * X   = createMatrixValue( arg_out_0, nrx, ncx ) ;
          valueType * Y   = createMatrixValue( arg_out_1, nrx, ncx ) ;
          valueType * S   = createMatrixValue( arg_out_2, nrx, ncx ) ;
          valueType * dst = createMatrixValue( arg_out_3, nrx, ncx ) ;

          mwSize size = nrx*ncx ;
          for ( mwSize i = 0 ; i < size ; ++i )
            *dst++ = ptr->closestPointBySample( ds, *x++, *y++, *X++, *Y++, *S++ ) ;
        }
        #undef CMD

      } else if ( cmd == "getPars" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurveMexWrapper('getPars',OBJ): expected 2 inputs");
        if ( nlhs > 0 ) setScalarValue(arg_out_0,ptr->getX0()) ;
        if ( nlhs > 1 ) setScalarValue(arg_out_1,ptr->getY0()) ;
        if ( nlhs > 2 ) setScalarValue(arg_out_2,ptr->getTheta0()) ;
        if ( nlhs > 3 ) setScalarValue(arg_out_3,ptr->getKappa()) ;
        if ( nlhs > 4 ) setScalarValue(arg_out_4,ptr->getKappa_D()) ;
        if ( nlhs > 5 ) setScalarValue(arg_out_5,ptr->getL()) ;

      } else if ( cmd == "getX0" ) {

        #define CMD "ClothoidCurveMexWrapper('getX0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getX0());

        #undef CMD

      } else if ( cmd == "getY0" ) {

        #define CMD "ClothoidCurveMexWrapper('getY0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getY0());

        #undef CMD

      } else if ( cmd == "getTheta0" ) {

        #define CMD "ClothoidCurveMexWrapper('getTheta0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getTheta0());

        #undef CMD

      } else if ( cmd == "getKappa0" ) {

        #define CMD "ClothoidCurveMexWrapper('getKappa0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa());

        #undef CMD

      } else if ( cmd == "getKappa_D" ) {

        #define CMD "ClothoidCurveMexWrapper('getKappa_D',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa_D());

        #undef CMD

      } else if ( cmd == "length" ) {

        #define CMD "ClothoidCurveMexWrapper('length',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getL());

        #undef CMD

      } else if ( cmd == "trim" ) {

        #define CMD "ClothoidCurveMexWrapper('trim',OBJ): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");

        valueType smin = getScalarValue(arg_in_2,CMD "Error in reading smin") ;
        valueType smax = getScalarValue(arg_in_3,CMD "Error in reading smax") ;

        ptr->trim(smin, smax);

        #undef CMD

      } else if ( cmd == "changeCurvilinearOrigin" ) {

        #define CMD "ClothoidCurveMexWrapper('changeCurvilinearOrigin',OBJ,s0,L): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");

        valueType s0 = getScalarValue(arg_in_2,CMD "Error in reading s0") ;
        valueType L  = getScalarValue(arg_in_3,CMD "Error in reading L") ;
        ptr->changeCurvilinearOrigin(s0,L);

        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "ClothoidCurveMexWrapper('rotate',OBJ,angle,cx,cy): "

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs");

        valueType angle = getScalarValue( arg_in_2, CMD "Error in reading angle" ) ;
        valueType cx    = getScalarValue( arg_in_3, CMD "Error in reading cx" ) ;
        valueType cy    = getScalarValue( arg_in_4, CMD "Error in reading cy" ) ;
        ptr->rotate(angle, cx, cy);

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "ClothoidCurveMexWrapper('translate',OBJ,tx,ty): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType tx = getScalarValue( arg_in_2, CMD "Error in reading tx" ) ;
        valueType ty = getScalarValue( arg_in_3, CMD "Error in reading ty" ) ;
        ptr->translate(tx, ty);

        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "ClothoidCurveMexWrapper('changeOrigin',OBJ,newX0,newY0): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType newX0 = getScalarValue( arg_in_2, CMD "Error in reading newX0" ) ;
        valueType newY0 = getScalarValue( arg_in_3, CMD "Error in reading newY0" ) ;
        ptr->changeOrigin(newX0, newY0);

        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "ClothoidCurveMexWrapper('scale',OBJ,s): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs");
        valueType s = getScalarValue( arg_in_2, CMD "Error in reading s" ) ;
        ptr->scale(s);

        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "ClothoidCurveMexWrapper('reverse',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        ptr->reverse();

        #undef CMD

      } else if ( cmd == "infinity" ) {

        #define CMD "ClothoidCurveMexWrapper('infinity',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 4, CMD "expected 4 output");

        valueType xp, yp, xm, ym ;
        ptr->Pinfinity( xp, yp, true );
        ptr->Pinfinity( xm, ym, false );

        setScalarValue(arg_out_0, xp ) ;
        setScalarValue(arg_out_1, yp );
        setScalarValue(arg_out_2, xm );
        setScalarValue(arg_out_3, ym );

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "ClothoidCurveMexWrapper('delete',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 0, CMD "expected no output");

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

        #undef CMD

        // Warn if other commands were ignored
      } else {
        MEX_ASSERT(false, "Unknown command: " << cmd );
      }

    } catch ( std::exception const & e ) {
    	mexErrMsgTxt(e.what()) ;

    } catch (...) {
  	  mexErrMsgTxt("clothoid failed\n") ;
    }
  }
}
