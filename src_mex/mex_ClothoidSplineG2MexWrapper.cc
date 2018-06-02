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

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"ClothoidSplineG2MexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = ClothoidListMexWrapper( 'new' ) ;\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    ClothoidListMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  - Build:\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, CLOT ) ;\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, x1, y1, theta1 ) ;\n" \
"\n" \
"  - Eval:\n" \
"    [x,y,theta,kappa] = ClothoidListMexWrapper( 'evaluate', OBJ, ss ) ;\n" \
"    x = ClothoidListMexWrapper( 'xBegin', OBJ ) ;\n" \
"    x = ClothoidListMexWrapper( 'xEnd', OBJ ) ;\n" \
"    y = ClothoidListMexWrapper( 'yBegin', OBJ ) ;\n" \
"    y = ClothoidListMexWrapper( 'yEnd', OBJ ) ;\n" \
"    theta = ClothoidListMexWrapper( 'thetaBegin', OBJ ) ;\n" \
"    theta = ClothoidListMexWrapper( 'thetaEnd', OBJ ) ;\n" \
"    kappa = ClothoidListMexWrapper( 'kappaBegin', OBJ ) ;\n" \
"    kappa = ClothoidListMexWrapper( 'kappaEnd', OBJ ) ;\n" \
"\n" \
"    [x,y]         = ClothoidListMexWrapper( 'eval', OBJ, ss, offs ) ;\n" \
"    [x_D,y_D]     = ClothoidListMexWrapper( 'eval_D', OBJ, ss, offs ) ;\n" \
"    [x_DD,y_DD]   = ClothoidListMexWrapper( 'eval_DD', OBJ, ss, offs ) ;\n" \
"    [x_DDD,y_DDD] = ClothoidListMexWrapper( 'eval_DDD', OBJ, ss, offs ) ;\n" \
"\n" \
"  - Transform:\n" \
"    ClothoidListMexWrapper( 'changeOrigin', OBJ, newX0, newY0 ) ;\n" \
"    ClothoidListMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"    ClothoidListMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"    ClothoidListMexWrapper( 'scale', OBJ, scaling ) ;\n" \
"    ClothoidListMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"  - Distance:\n" \
"    [X,Y,s,dst] = ClothoidListMexWrapper( 'closestPoint', OBJ, x, y ) ;\n" \
"    [dst,s]     = ClothoidListMexWrapper( 'distance', OBJ, x, y ) ;\n" \
"\n" \
"  - Bounding Box:\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size ) ;%\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size, offs ) ;%\n" \
"\n" \
"  - G2 spline:\n" \
"    ok = ClothoidListMexWrapper( 'build3arcG2' || 'build2arcG2' || 'build3arcCLC, ...\n" \
"                                 OBJ, ...\n" \
"                                 x0, y0, theta0, kappa0, ...\n" \
"                                 x1, y1, theta1, kappa1 ) ;%\n" \
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
  ClothoidSplineG2 *
  DATA_NEW( mxArray * & mx_id ) {
    ClothoidSplineG2 * ptr = new ClothoidSplineG2();
    mx_id = convertPtr2Mat<ClothoidSplineG2>(ptr);
    return ptr ;
  }

  static
  inline
  ClothoidSplineG2 *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<ClothoidSplineG2>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<ClothoidSplineG2>(mx_id);
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

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidListMexWrapper: First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;

      bool do_new = cmd == "new" ;
      ClothoidSplineG2 * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "ClothoidSplineG2MexWrapper, expected 1 output" );

      } else if ( cmd == "setup" ) {

        #define CMD "ClothoidSplineG2MexWrapper('setup',OBJ,x,y): "

        MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs") ;
        MEX_ASSERT( nlhs == 0, CMD "expected no output") ;

        mwSize nx, ny;
        valueType const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
        valueType const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

        MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " must be equal to size(y) = " << ny ) ;

        ptr->setup( x, y, nx ) ;

        #undef CMD

      } else if ( cmd == "target" ) {

        #define CMD "ClothoidSplineG2MexWrapper('target',OBJ,target[,theta0,theta1]): "

        MEX_ASSERT( nrhs >= 3, CMD "expected at least 3 inputs") ;
        MEX_ASSERT( nlhs == 0, CMD "expected no output") ;

        MEX_ASSERT( mxIsChar(arg_in_2), CMD "Third argument must be a string" ) ;
        string obj = mxArrayToString(arg_in_2) ;

        if ( obj == "P1" ) {
          MEX_ASSERT( nrhs == 5, CMD "expected at 5 inputs") ;
          valueType theta0 = getScalarValue( arg_in_3, CMD "Error in reading theta0" ) ;
          valueType theta1 = getScalarValue( arg_in_4, CMD "Error in reading theta1" ) ;
          ptr->setP1( theta0, theta1 );
        } else {
          MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs") ;
          if      ( obj == "P2" ) ptr->setP2();
          else if ( obj == "P3" ) ptr->setP3();
          else if ( obj == "P4" ) ptr->setP4();
          else if ( obj == "P5" ) ptr->setP5();
          else if ( obj == "P6" ) ptr->setP6();
          else if ( obj == "P7" ) ptr->setP7();
          else if ( obj == "P8" ) ptr->setP8();
          else if ( obj == "P9" ) ptr->setP9();
          else {
            MEX_ASSERT( false, CMD "Unknown target " << obj ) ;
          }
        }
        #undef CMD

      } else if ( cmd == "guess" ) {

        #define CMD "ClothoidSplineG2MexWrapper('guess',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs") ;
        MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs") ;

        indexType N = ptr->numPnts() ;

        valueType * theta_guess = createMatrixValue( arg_out_0, N, 1 );
        valueType * theta_min   = createMatrixValue( arg_out_1, N, 1 );
        valueType * theta_max   = createMatrixValue( arg_out_2, N, 1 );

        ptr->guess( theta_guess, theta_min, theta_max );

        #undef CMD

      } else if ( cmd == "objective" ) {

        #define CMD "ClothoidSplineG2MexWrapper('objective',OBJ,theta): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs") ;
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs") ;

        mwSize ntheta;
        valueType const * theta = getVectorPointer( arg_in_2, ntheta, CMD "Error in reading theta" );
        MEX_ASSERT( ntheta == ptr->numPnts(),
                    CMD "length(theta) = " << ntheta << " must be " << ptr->numPnts() ) ;
        valueType f ;
        bool ok = ptr->objective( theta, f );
        if ( ok ) setScalarValue( arg_out_0, f );

        #undef CMD
#if 0
      } else if ( cmd == "evaluate" ) {

        #define CMD "ClothoidListMexWrapper('evaluate',OBJ,s): "

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

        #define CMD "ClothoidListMexWrapper('eval*',OBJ,s[,offs]): "

        MEX_ASSERT( nrhs == 3 || nrhs == 4, CMD "expected 3 or 4 inputs") ;
        MEX_ASSERT( nlhs == 1 || nlhs == 2, CMD "expected 1 or 2 outputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, CMD "Error in reading s" );

        double offs = 0 ;
        if ( nrhs == 4 ) offs = getScalarValue( arg_in_3, CMD "Error in reading offs" ) ;

        if ( nlhs == 1 ) {
          double * xyVals = createMatrixValue( arg_out_0, 2, size );
          if ( cmd == "eval" ) {
            for ( mwSize i=0; i < size ; ++i )
              ptr->eval( sVals[i], offs, xyVals[2*i], xyVals[2*i+1] );
          } else if ( cmd == "eval_D" ) {
            for ( mwSize i=0; i < size ; ++i )
              ptr->eval_D( sVals[i], offs, xyVals[2*i], xyVals[2*i+1] );
          } else if ( cmd == "eval_DD" ) {
            for ( mwSize i=0; i < size ; ++i )
              ptr->eval_DD( sVals[i], offs, xyVals[2*i], xyVals[2*i+1] );
          } else {
            for ( mwSize i=0; i < size ; ++i )
              ptr->eval_DDD( sVals[i], offs, xyVals[2*i], xyVals[2*i+1] );
          }
        } else {
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
        }

        #undef CMD

      } else if ( cmd == "distance" ) {

        #define CMD "ClothoidListMexWrapper('distance',OBJ,x,y): "
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

        #define CMD "ClothoidListMexWrapper('closestPoint',OBJ,x,y): "
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

      } else if ( cmd == "xBegin"     || cmd == "xEnd"     ||
                  cmd == "yBegin"     || cmd == "yEnd"     ||
                  cmd == "thetaBegin" || cmd == "thetaEnd" ||
                  cmd == "kappaBegin" || cmd == "kappaEnd" ||
                  cmd == "length" ) {


        #define CMD "ClothoidListMexWrapper('...',OBJ[,n]): "
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");

        if ( nrhs == 3 ) {
          int64_t n = getInt( arg_in_2, CMD "Error in reading n" );
          MEX_ASSERT( n > 0 && n <= ptr->numSegment(),
                      CMD "n =  " << n << " must be >= 1 and <= " << ptr->numSegment() );
          ClothoidCurve const & c = ptr->get(n-1);
          if      ( cmd == "xBegin"     ) setScalarValue(arg_out_0, c.xBegin());
          else if ( cmd == "xEnd"       ) setScalarValue(arg_out_0, c.xEnd());
          else if ( cmd == "yBegin"     ) setScalarValue(arg_out_0, c.yBegin());
          else if ( cmd == "yEnd"       ) setScalarValue(arg_out_0, c.yEnd());
          else if ( cmd == "thetaBegin" ) setScalarValue(arg_out_0, c.thetaBegin());
          else if ( cmd == "thetaEnd"   ) setScalarValue(arg_out_0, c.thetaEnd());
          else if ( cmd == "kappaBegin" ) setScalarValue(arg_out_0, c.kappaBegin());
          else if ( cmd == "kappaEnd"   ) setScalarValue(arg_out_0, c.kappaEnd());
          else if ( cmd == "length"     ) setScalarValue(arg_out_0, c.length());
        } else {
          MEX_ASSERT(nrhs == 2, CMD "expected 2 or 3 inputs");
          if      ( cmd == "xBegin"     ) setScalarValue(arg_out_0, ptr->xBegin());
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

        MEX_ASSERT(nrhs == 5, CMD "expected 5 inputs");

        valueType angle = getScalarValue( arg_in_2, CMD "Error in reading angle" ) ;
        valueType cx    = getScalarValue( arg_in_3, CMD "Error in reading cx" ) ;
        valueType cy    = getScalarValue( arg_in_4, CMD "Error in reading cy" ) ;
        ptr->rotate(angle, cx, cy);

        #undef CMD

      } else if ( cmd == "translate" ) {

        #define CMD "ClothoidListMexWrapper('translate',OBJ,tx,ty): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType tx = getScalarValue( arg_in_2, CMD "Error in reading tx" ) ;
        valueType ty = getScalarValue( arg_in_3, CMD "Error in reading ty" ) ;
        ptr->translate(tx, ty);

        #undef CMD
#endif

      } else if ( cmd == "delete" ) {

        #define CMD "ClothoidListMexWrapper('delete',OBJ): "

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
