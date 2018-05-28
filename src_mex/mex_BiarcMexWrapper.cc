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
"    OBJ = BiarcMexWrapper( 'new' ) ;\n" \
"\n" \
"    On output:\n" \
"       OBJ   = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    BiarcMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  - Build:\n" \
"    BiarcMexWrapper( 'build', OBJ, x0, y0, theta0, x1, y1, theta1 ) ;\n" \
"\n" \
"  - Eval:\n" \
"    [x,y,theta,kappa] = BiarcMexWrapper( 'evaluate', OBJ, ss ) ;\n" \
"    [x0,y0,theta0,kappa0,L0,x1,y1,theta1,kappa1,L1] = BiarcMexWrapper( 'getPars', OBJ ) ;\n" \
"\n" \
"    [x,y]         = BiarcMexWrapper( 'eval', OBJ, ss ) ;\n" \
"    [x_D,y_D]     = BiarcMexWrapper( 'eval_D', OBJ, ss ) ;\n" \
"    [x_DD,y_DD]   = BiarcMexWrapper( 'eval_DD', OBJ, ss ) ;\n" \
"    [x_DDD,y_DDD] = BiarcMexWrapper( 'eval_DDD', OBJ, ss ) ;\n" \
"\n" \
"  - Transform:\n" \
"    BiarcMexWrapper( 'changeOrigin', OBJ, newX0, newY0 ) ;\n" \
"    BiarcMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"    BiarcMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"    BiarcMexWrapper( 'scale', OBJ, scaling ) ;\n" \
"    BiarcMexWrapper( 'reverse', OBJ ) ;\n" \
"  - Distance:\n" \
"    [X,Y,s,dst] = BiarcMexWrapper( 'closestPoint', OBJ, x, y ) ;\n" \
"    [dst,s] = BiarcMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"    [X,Y,s,dst] = BiarcMexWrapper( 'closestPointBySample', OBJ, x, y, ds ) ;\n" \
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
    return ptr ;
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
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "BiarcMexWrapper: First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;

      bool do_new = cmd == "new" ;
      Biarc * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "BiarcMexWrapper, expected 1 output" );

      } else if ( cmd == "build" ) {

        #define CMD "BiarcMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1): "

        MEX_ASSERT( nrhs == 8 , CMD "expected 8 inputs") ;

        valueType x0(0), y0(0), theta0(0), x1(0), y1(0), theta1(0);

        x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" ) ;
        y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" ) ;
        theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" ) ;
        x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" ) ;
        y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" ) ;
        theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" ) ;

        ptr->build( x0, y0, theta0, x1, y1, theta1 );

        #undef CMD

      } else if ( cmd == "evaluate" ) {

        #define CMD "BiarcMexWrapper('evaluate',OBJ,s): "

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

        #define CMD "BiarcMexWrapper('eval*',OBJ,s): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs") ;
        MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, CMD "Error in reading s" );

        double * xVals = createMatrixValue( arg_out_0, size, 1 );
        double * yVals = createMatrixValue( arg_out_1, size, 1 );
        if ( cmd == "eval" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval( sVals[i], xVals[i], yVals[i] );
        } else if ( cmd == "eval_D" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_D( sVals[i], xVals[i], yVals[i] );
        } else if ( cmd == "eval_DD" ) {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_DD( sVals[i], xVals[i], yVals[i] );
        } else {
          for ( mwSize i=0; i < size ; ++i )
            ptr->eval_DDD( sVals[i], xVals[i], yVals[i] );
        }

        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "BiarcMexWrapper('distance',OBJ,x,y): "
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

      } else if ( cmd == "closestPoint" ) {

        #define CMD "BiarcMexWrapper('closestPoint',OBJ,x,y): "
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

      } else if ( cmd == "getPars" ) {

        MEX_ASSERT(nrhs == 2, "BiarcMexWrapper('getPars',OBJ): expected 2 inputs");
        if ( nlhs > 0 ) setScalarValue(arg_out_0,ptr->getX0()) ;
        if ( nlhs > 1 ) setScalarValue(arg_out_1,ptr->getY0()) ;
        if ( nlhs > 2 ) setScalarValue(arg_out_2,ptr->getTheta0()) ;
        if ( nlhs > 3 ) setScalarValue(arg_out_3,ptr->getKappa0()) ;
        if ( nlhs > 4 ) setScalarValue(arg_out_4,ptr->getL0()) ;

        if ( nlhs > 5 ) setScalarValue(arg_out_5,ptr->getX1()) ;
        if ( nlhs > 6 ) setScalarValue(arg_out_6,ptr->getY1()) ;
        if ( nlhs > 7 ) setScalarValue(arg_out_7,ptr->getTheta1()) ;
        if ( nlhs > 8 ) setScalarValue(arg_out_8,ptr->getKappa1()) ;
        if ( nlhs > 9 ) setScalarValue(arg_out_9,ptr->getL1()) ;

      } else if ( cmd == "getX0" ) {

        #define CMD "BiarcMexWrapper('getX0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getX0());

        #undef CMD

      } else if ( cmd == "getX1" ) {

        #define CMD "BiarcMexWrapper('getX1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getX1());

        #undef CMD

      } else if ( cmd == "getY0" ) {

        #define CMD "BiarcMexWrapper('getY0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getY0());

        #undef CMD

      } else if ( cmd == "getY1" ) {

        #define CMD "BiarcMexWrapper('getY1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getY1());

        #undef CMD

      } else if ( cmd == "getTheta0" ) {

        #define CMD "BiarcMexWrapper('getTheta0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getTheta0());

        #undef CMD

      } else if ( cmd == "getTheta1" ) {

        #define CMD "BiarcMexWrapper('getTheta1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getTheta1());

        #undef CMD

      } else if ( cmd == "getKappa0" ) {

        #define CMD "BiarcMexWrapper('getKappa0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa0());

        #undef CMD

      } else if ( cmd == "getKappa1" ) {

        #define CMD "BiarcMexWrapper('getKappa1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa1());

        #undef CMD

      } else if ( cmd == "getL0" ) {

        #define CMD "BiarcMexWrapper('getL0',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getL0());

        #undef CMD

      } else if ( cmd == "getL1" ) {

        #define CMD "BiarcMexWrapper('getL1',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getL1());

        #undef CMD

      } else if ( cmd == "length" ) {

        #define CMD "BiarcMexWrapper('length',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getL());

        #undef CMD

      } else if ( cmd == "rotate" ) {

        #define CMD "BiarcMexWrapper('rotate',OBJ,angle,cx,cy): "

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs");

        valueType angle = getScalarValue( arg_in_2, CMD "Error in reading angle" ) ;
        valueType cx    = getScalarValue( arg_in_3, CMD "Error in reading cx" ) ;
        valueType cy    = getScalarValue( arg_in_4, CMD "Error in reading cy" ) ;
        ptr->rotate(angle, cx, cy);

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "BiarcMexWrapper('translate',OBJ,tx,ty): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType tx = getScalarValue( arg_in_2, CMD "Error in reading tx" ) ;
        valueType ty = getScalarValue( arg_in_3, CMD "Error in reading ty" ) ;
        ptr->translate(tx, ty);

        #undef CMD

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "BiarcMexWrapper('changeOrigin',OBJ,newX0,newY0): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType newX0 = getScalarValue( arg_in_2, CMD "Error in reading newX0" ) ;
        valueType newY0 = getScalarValue( arg_in_3, CMD "Error in reading newY0" ) ;
        ptr->changeOrigin(newX0, newY0);

        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "BiarcMexWrapper('scale',OBJ,s): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs");
        valueType s = getScalarValue( arg_in_2, CMD "Error in reading s" ) ;
        ptr->scale(s);

        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "BiarcMexWrapper('reverse',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        ptr->reverse();

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "BiarcMexWrapper('delete',OBJ): "

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
  	  mexErrMsgTxt("Biarc failed\n") ;
    }
  }
}
