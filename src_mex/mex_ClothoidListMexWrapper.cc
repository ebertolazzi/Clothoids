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
"ClothoidListMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
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
"    x = ClothoidListMexWrapper( 'Xbegin', OBJ ) ;\n" \
"    x = ClothoidListMexWrapper( 'Xend', OBJ ) ;\n" \
"    y = ClothoidListMexWrapper( 'Ybegin', OBJ ) ;\n" \
"    y = ClothoidListMexWrapper( 'Yend', OBJ ) ;\n" \
"    theta = ClothoidListMexWrapper( 'ThetaBegin', OBJ ) ;\n" \
"    theta = ClothoidListMexWrapper( 'ThetaEnd', OBJ ) ;\n" \
"    kappa = ClothoidListMexWrapper( 'KappaBegin', OBJ ) ;\n" \
"    kappa = ClothoidListMexWrapper( 'KappaEnd', OBJ ) ;\n" \
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
"%   TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size ) ;%\n" \
"%   TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size, offs ) ;%\n" \
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
    return ptr ;
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
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidListMexWrapper: First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;

      bool do_new = cmd == "new" ;
      ClothoidList * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "ClothoidListMexWrapper, expected 1 output" );

      } else if ( cmd == "push_back" ) {

        #define CMD "ClothoidListMexWrapper('push_back',OBJ,[x1,y1,theta1]|[CLOT]): "

        MEX_ASSERT( nrhs == 3 || nrhs == 5 , CMD "expected 3 or 5 inputs") ;

        if ( nrhs == 5 ) {
          valueType x1     = getScalarValue( arg_in_2, CMD "Error in reading x1" ) ;
          valueType y1     = getScalarValue( arg_in_3, CMD "Error in reading y1" ) ;
          valueType theta1 = getScalarValue( arg_in_4, CMD "Error in reading theta1" ) ;
          ptr->push_back( x1, y1, theta1 );
        } else {
          ClothoidCurve * cc = convertMat2Ptr<ClothoidCurve>(arg_in_2);
          ptr->push_back( *cc );
        }

        #undef CMD

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

      } else if ( cmd == "xBegin" ) {

        #define CMD "ClothoidListMexWrapper('xBegin',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->xBegin());

        #undef CMD

      } else if ( cmd == "xEnd" ) {

        #define CMD "ClothoidListMexWrapper('xEnd',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->xEnd());

        #undef CMD

      } else if ( cmd == "yBegin" ) {

        #define CMD "ClothoidListMexWrapper('yBegin',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->yBegin());

        #undef CMD

      } else if ( cmd == "yEnd" ) {

        #define CMD "ClothoidListMexWrapper('yEnd',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->yEnd());

        #undef CMD

      } else if ( cmd == "thetaBegin" ) {

        #define CMD "ClothoidListMexWrapper('thetaBegin',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->thetaBegin());

        #undef CMD

      } else if ( cmd == "thetaEnd" ) {

        #define CMD "ClothoidListMexWrapper('thetaEnd',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->thetaEnd());

        #undef CMD

      } else if ( cmd == "kappaBegin" ) {

        #define CMD "ClothoidListMexWrapper('kappaBegin',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->kappaBegin());

        #undef CMD

      } else if ( cmd == "kappaEnd" ) {

        #define CMD "ClothoidListMexWrapper('kappaEnd',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        setScalarValue(arg_out_0, ptr->kappaEnd());

        #undef CMD

      } else if ( cmd == "length" ) {

        #define CMD "ClothoidListMexWrapper('length',OBJ[,k]): "

        MEX_ASSERT(nrhs == 2 || nrhs == 3, CMD "expected 2 or 3 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 outputs");
        if ( nrhs == 2 ) {
          setScalarValue(arg_out_0, ptr->totalLength());
        } else {
          int64_t idx = getInt( arg_in_2, CMD "Error in reading k" );
          ClothoidCurve const & c = ptr->get(idx);
          setScalarValue(arg_out_0, c.length());
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

      } else if ( cmd == "changeOrigin" ) {

        #define CMD "ClothoidListMexWrapper('changeOrigin',OBJ,newX0,newY0): "

        MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs");
        valueType newX0 = getScalarValue( arg_in_2, CMD "Error in reading newX0" ) ;
        valueType newY0 = getScalarValue( arg_in_3, CMD "Error in reading newY0" ) ;
        ptr->changeOrigin(newX0, newY0);

        #undef CMD

      } else if ( cmd == "scale" ) {

        #define CMD "ClothoidListMexWrapper('scale',OBJ,s): "

        MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs");
        valueType s = getScalarValue( arg_in_2, CMD "Error in reading s" ) ;
        ptr->scale(s);

        #undef CMD

      } else if ( cmd == "reverse" ) {

        #define CMD "ClothoidListMexWrapper('reverse',OBJ): "

        MEX_ASSERT(nrhs == 2, CMD "expected 2 inputs");
        ptr->reverse();

        #undef CMD

      } else if ( cmd == "bbox" ) {

        #define CMD "ClothoidListMexWrapper('bbox', OBJ, max_angle, max_size [,offs]): "

        MEX_ASSERT(nrhs == 4 || nrhs == 5, CMD "expected 4 or 5 inputs");
        MEX_ASSERT(nlhs == 1, CMD "expected 1 output");

        valueType max_angle = getScalarValue( arg_in_2, CMD "Error in reading max_angle" ) ;
        valueType max_size  = getScalarValue( arg_in_3, CMD "Error in reading max_size" ) ;
        valueType offs      = 0 ;
        if ( nrhs == 5 ) offs = getScalarValue( arg_in_4, CMD "Error in reading offs" ) ;

        vector<ClothoidCurve::bbData> bb ;
        ptr->bbSplit( max_angle, max_size, offs, bb ) ;

        plhs[0] = mxCreateDoubleMatrix(6, bb.size(), mxREAL);
        double * pT = mxGetPr(plhs[0]);
        for ( int i = 0 ; i < bb.size() ; ++i ) {
          T2D const & t = bb[i].t ;
          *pT++ = t.x1() ; *pT++ = t.y1() ;
          *pT++ = t.x2() ; *pT++ = t.y2() ;
          *pT++ = t.x3() ; *pT++ = t.y3() ;
        }
        #undef CMD

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
