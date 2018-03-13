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
"    ClothoidCurveMexWrapper( 'changeOrigin', OBJ, s0, L ) ;\n" \
"    ClothoidCurveMexWrapper( 'moveOrigin', OBJ, newX0, newY0 ) ;\n" \
"    ClothoidCurveMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;\n" \
"    ClothoidCurveMexWrapper( 'translate', OBJ, tx, ty ) ;\n" \
"    ClothoidCurveMexWrapper( 'scale', OBJ, scaling ) ;\n" \
"    ClothoidCurveMexWrapper( 'reverse', OBJ ) ;\n" \
"\n" \
"=====================================================================================\n" \
"\n" \
"Autors: Enrico Bertolazzi^(1), Marco Frego^(2) and Paolo Bevilacqua^(2)\n" \
"  (1) Department of Industrial Engineering\n" \
"  (2)Department of Information Engineering and Computer Science\n" \
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

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidCurve: First argument must be a string" ) ;
      string cmd = mxArrayToString(arg_in_0) ;

      bool do_new = cmd == "new" ;
      ClothoidCurve * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new || cmd == "build"  ) {

        if ( do_new ) {
          MEX_ASSERT( nlhs == 1, "ClothoidCurve, expected 1 output" );
        } else {
          MEX_ASSERT( nlhs == 0, "ClothoidCurve, expected no output" );
        }

        indexType kk = do_new ? 0 : 1 ;

        valueType x0(0), y0(0), theta0(0), k0(0), dk(0), L(0);

        if ( do_new && nrhs == 1 ) {
          // nothing to do
        } else if ( nrhs == 1+kk ) {

          MEX_ASSERT( mxIsStruct(arg_in_1), "Second argument must be a struct" ) ;

          mxArray * mx_x0     = mxGetField( prhs[1+kk], 0, "x0" ) ;
          mxArray * mx_y0     = mxGetField( prhs[1+kk], 0, "y0" ) ;
          mxArray * mx_theta0 = mxGetField( prhs[1+kk], 0, "theta0" ) ;
          mxArray * mx_k0     = mxGetField( prhs[1+kk], 0, "k0" ) ;
          mxArray * mx_dk     = mxGetField( prhs[1+kk], 0, "dk" ) ;

          MEX_ASSERT( mx_x0     != nullptr, "ClothoidCurve('build',OBJ,struct): Field `x0` is missing" );
          MEX_ASSERT( mx_y0     != nullptr, "ClothoidCurve('build',OBJ,struct): Field `y0` is missing" );
          MEX_ASSERT( mx_theta0 != nullptr, "ClothoidCurve('build',OBJ,struct): Field `theta0` is missing" );
          MEX_ASSERT( mx_k0     != nullptr, "ClothoidCurve('build',OBJ,struct): Field `k0` is missing" );
          MEX_ASSERT( mx_dk     != nullptr, "ClothoidCurve('build',OBJ,struct): Field `dk` is missing" );

          x0     = getScalarValue( mx_x0,     "ClothoidCurve('build',OBJ,struct): Field `x0` must be a real double scalar" ) ;
          y0     = getScalarValue( mx_y0,     "ClothoidCurve('build',OBJ,struct): Field `y0` must be a real double scalar" ) ;
          theta0 = getScalarValue( mx_theta0, "ClothoidCurve('build',OBJ,struct): Field `theta0` must be a real double scalar" ) ;
          k0     = getScalarValue( mx_k0,     "ClothoidCurve('build',OBJ,struct): Field `k0` must be a real double scalar" ) ;
          dk     = getScalarValue( mx_dk,     "ClothoidCurve('build',OBJ,struct): Field `dk` must be a real double scalar" ) ;

        } else if ( nrhs == 7+kk ) {

          x0     = getScalarValue(prhs[1+kk],"Error in reading x0") ;
          y0     = getScalarValue(prhs[2+kk],"Error in reading y0") ;
          theta0 = getScalarValue(prhs[3+kk],"Error in reading theta0") ;
          k0     = getScalarValue(prhs[4+kk],"Error in reading k0") ;
          dk     = getScalarValue(prhs[5+kk],"Error in reading dk") ;
          L      = getScalarValue(prhs[6+kk],"Error in reading L") ;
        } else {
          MEX_ASSERT( false, "expected " << kk+1 << " or " << kk+7 << " inputs") ;
        }

        ptr->build( x0, y0, theta0, k0, dk, L );

      } else if ( cmd == "build_G1" ) {

        MEX_ASSERT( nrhs == 8 , "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1):  expected 8 inputs") ;

        valueType x0(0), y0(0), theta0(0), x1(0), y1(0), theta1(0);

        x0     = getScalarValue( arg_in_2, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading x0" ) ;
        y0     = getScalarValue( arg_in_3, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading y0" ) ;
        theta0 = getScalarValue( arg_in_4, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading theta0" ) ;
        x1     = getScalarValue( arg_in_5, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading x1" ) ;
        y1     = getScalarValue( arg_in_6, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading y1" ) ;
        theta1 = getScalarValue( arg_in_7, "ClothoidCurve('build_G1',OBJ,x0,y0,theta0,x1,y1,theta1): Error in reading theta1" ) ;

        ptr->build_G1(x0, y0, theta0, x1, y1, theta1);

      } else if ( cmd == "build_forward" ) {

        MEX_ASSERT( nrhs == 8 , "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): expected 8 inputs") ;

        valueType x0(0), y0(0), theta0(0), kappa0(0), x1(0), y1(0);

        x0     = getScalarValue( arg_in_2, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading x0" ) ;
        y0     = getScalarValue( arg_in_3, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading y0" ) ;
        theta0 = getScalarValue( arg_in_4, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading theta0" ) ;
        kappa0 = getScalarValue( arg_in_5, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading kappa0" ) ;
        x1     = getScalarValue( arg_in_6, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading x1" ) ;
        y1     = getScalarValue( arg_in_6, "ClothoidCurve('build_forward',OBJ,x0,y0,theta0,kappa0,x1,y1): Error in reading y1" ) ;

        bool res = ptr->build_forward(x0, y0, theta0, kappa0, x1, y1);

        // returns the status of the interpolation
        mwSize dims[2] = {1,1} ;
        arg_out_0 = mxCreateLogicalArray(2, dims);
        ((bool*)mxGetPr(arg_out_0))[0] = res;
        //mexPrintf("Result: %d\n", res);

      } else if ( cmd == "evaluate" ) {

        MEX_ASSERT( nrhs == 3 , "ClothoidCurve('evaluate',OBJ,s): expected 3 inputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, "Error in reading s" );

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
          MEX_ASSERT( false, "ClothoidCurve('evaluate',OBJ,s): expected 2 or 4 outputs") ;
        }

      } else if ( cmd == "eval"    || cmd == "eval_D" ||
                  cmd == "eval_DD" || cmd == "eval_DDD" ) {

        MEX_ASSERT( nrhs == 3 || nrhs == 4, "ClothoidCurve('evaluate',OBJ,s[,offs]): expected 3 or 4 inputs") ;
        MEX_ASSERT( nlhs == 2, "ClothoidCurve('evaluate',OBJ,s): expected 2 outputs") ;

        mwSize size;
        double const * sVals = getVectorPointer( arg_in_2, size, "ClothoidCurve('eval*',OBJ,s[,offs]): Error in reading s" );

        double offs = 0 ;
        if ( nrhs == 4 ) offs = getScalarValue( arg_in_3, "ClothoidCurve('eval*',OBJ,s[,offs]): Error in reading offs" ) ;

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

      } else if ( cmd == "closestPoint" ) {

        MEX_ASSERT( nrhs == 5 , "ClothoidCurve('closestPoint',OBJ,x,y,ds): expected 5 inputs") ;
        MEX_ASSERT( nlhs == 4 , "ClothoidCurve('closestPoint',OBJ,x,y,ds): expected 4 outputs") ;

        valueType x, y, ds, X, Y, S ;

        x  = getScalarValue(arg_in_2,"ClothoidCurve('closestPoint',OBJ,x,y,ds): Error in reading x") ;
        y  = getScalarValue(arg_in_3,"ClothoidCurve('closestPoint',OBJ,x,y,ds): Error in reading y") ;
        ds = getScalarValue(arg_in_4,"ClothoidCurve('closestPoint',OBJ,x,y,ds): Error in reading ds") ;

        valueType DST = ptr->closestPoint( x, y, ds, X, Y, S );

        if ( nlhs > 0 ) setScalarValue(arg_out_0,S) ;
        if ( nlhs > 1 ) setScalarValue(arg_out_1,X) ;
        if ( nlhs > 2 ) setScalarValue(arg_out_2,Y) ;
        if ( nlhs > 3 ) setScalarValue(arg_out_3,DST) ;

      } else if ( cmd == "getPars" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getPars',OBJ): expected 2 inputs");
        if ( nlhs > 0 ) setScalarValue(arg_out_0,ptr->getX0()) ;
        if ( nlhs > 1 ) setScalarValue(arg_out_1,ptr->getY0()) ;
        if ( nlhs > 2 ) setScalarValue(arg_out_2,ptr->getTheta0()) ;
        if ( nlhs > 3 ) setScalarValue(arg_out_3,ptr->getKappa()) ;
        if ( nlhs > 4 ) setScalarValue(arg_out_4,ptr->getKappa_D()) ;
        if ( nlhs > 5 ) setScalarValue(arg_out_5,ptr->getL()) ;

      } else if ( cmd == "getX0" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getX0',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('getX0',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getX0());

      } else if ( cmd == "getY0" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getY0',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('getY0',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getY0());

      } else if ( cmd == "getTheta0" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getTheta0',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('getTheta0',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getTheta0());

      } else if ( cmd == "getKappa0" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getKappa0',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('getKappa0',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa());

      } else if ( cmd == "getKappa_D" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('getKappa_D',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('getKappa_D',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getKappa_D());

      } else if ( cmd == "length" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('length',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 1, "ClothoidCurve('length',OBJ): expected 1 outputs");
        setScalarValue(arg_out_0, ptr->getL());

      } else if ( cmd == "trim" ) {

        MEX_ASSERT(nrhs == 4, "ClothoidCurve('trim',OBJ,smin,smax): expected 4 inputs");

        valueType smin = getScalarValue(arg_in_2,"ClothoidCurve('trim',OBJ,smin,smax): Error in reading smin") ;
        valueType smax = getScalarValue(arg_in_3,"ClothoidCurve('trim',OBJ,smin,smax): Error in reading smax") ;

        ptr->trim(smin, smax);

      } else if ( cmd == "changeCurvilinearOrigin" ) {

        MEX_ASSERT(nrhs == 4, "ClothoidCurve('changeCurvilinearOrigin',OBJ,s0,L): expected 4 inputs");

        valueType s0 = getScalarValue(arg_in_2,"ClothoidCurve('changeCurvilinearOrigin',OBJ,s0,L): Error in reading s0") ;
        valueType L  = getScalarValue(arg_in_3,"ClothoidCurve('changeCurvilinearOrigin',OBJ,s0,L): Error in reading L") ;
        ptr->changeCurvilinearOrigin(s0,L);

      } else if ( cmd == "rotate" ) {

        MEX_ASSERT(nrhs == 5, "ClothoidCurve('rotate',OBJ,angle,cx,cy): expected 5 inputs");

        valueType angle = getScalarValue( arg_in_2, "ClothoidCurve('rotate',OBJ,angle,cx,cy): Error in reading angle" ) ;
        valueType cx    = getScalarValue( arg_in_3, "ClothoidCurve('rotate',OBJ,angle,cx,cy): Error in reading cx" ) ;
        valueType cy    = getScalarValue( arg_in_4, "ClothoidCurve('rotate',OBJ,angle,cx,cy): Error in reading cy" ) ;
        ptr->rotate(angle, cx, cy);

      } else if ( cmd == "translate" ) {

        MEX_ASSERT(nrhs == 4, "expected 4 inputs");
        valueType tx = getScalarValue( arg_in_2, "ClothoidCurve('translate',OBJ,tx,ty): Error in reading tx" ) ;
        valueType ty = getScalarValue( arg_in_3, "ClothoidCurve('translate',OBJ,tx,ty): Error in reading ty" ) ;
        ptr->translate(tx, ty);

      } else if ( cmd == "moveOrigin" ) {

        MEX_ASSERT(nrhs == 4, "ClothoidCurve('moveOrigin',OBJ,newX0,newY0): expected 4 inputs");
        valueType newX0 = getScalarValue( arg_in_2, "ClothoidCurve('moveOrigin',OBJ,newX0,newY0): Error in reading newX0" ) ;
        valueType newY0 = getScalarValue( arg_in_3, "ClothoidCurve('moveOrigin',OBJ,newX0,newY0): Error in reading newY0" ) ;
        ptr->moveOrigin(newX0, newY0);

      } else if ( cmd == "scale" ) {

        MEX_ASSERT(nrhs == 3, "ClothoidCurve('scale',OBJ,s): expected 3 inputs");
        valueType s = getScalarValue( arg_in_2, "ClothoidCurve('scale',OBJ,s): Error in reading s" ) ;
        ptr->scale(s);

      } else if ( cmd == "reverse" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('reverse',OBJ): expected 2 inputs");
        ptr->reverse();

      } else if ( cmd == "delete" ) {

        MEX_ASSERT(nrhs == 2, "ClothoidCurve('delete',OBJ): expected 2 inputs");
        MEX_ASSERT(nlhs == 0, "ClothoidCurve('delete',OBJ): expected no output");

        // Destroy the C++ object
        DATA_DELETE(arg_in_1);

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
