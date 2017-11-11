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
"    OBJ = ClothoidCurveMexWrapper( 'new', x0, y0, theta0, k0, dk, smin, smax ) ;\n" \
"\n" \
"    On input:\n" \
"      x0, y0 = coordinate of initial point\n" \
"      theta0 = orientation (angle) of the clothoid at initial point" \
"      k0     = curvature of the clothoid at initial point\n" \
"      dk     = derivative of curvature respect to arclength\n" \
"      L      = length of the clothoid curve from initial to final point\n" \
"      smin   = initial curvilinear coordinate of the curve\n" \
"      smax   = final curvilinear coordinate of the curve\n" \
"\n" \
"     On output:\n" \
"       OBJ   = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    ClothoidCurveMexWrapper( 'delete', OBJ ) ;\n" \
"\n" \
"  - Build:\n" \
"    ClothoidCurveMexWrapper( 'build', OBJ, x0, y0, theta0, k0, dk, L ) ;\n" \
"    ClothoidCurveMexWrapper( 'build', OBJ, x0, y0, theta0, k0, dk, smin, smax ) ;\n" \
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
"    ClothoidCurveMexWrapper( 'changeOrigin', OBJ, s0 ) ;\n" \
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

static
Clothoid::ClothoidCurve *
DATA_NEW( mxArray * & mx_id ) {
  Clothoid::ClothoidCurve * ptr = new Clothoid::ClothoidCurve();
  mx_id = convertPtr2Mat<Clothoid::ClothoidCurve>(ptr);
  return ptr ;
}

static
inline
Clothoid::ClothoidCurve *
DATA_GET( mxArray const * & mx_id ) {
  return convertMat2Ptr<Clothoid::ClothoidCurve>(mx_id);
}

static
void
DATA_DELETE( mxArray const * & mx_id ) {
  Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(mx_id);
  delete ptr ;
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
    std::string cmd = mxArrayToString(arg_in_0) ;

    bool do_new = cmd == "new" ;
    Clothoid::ClothoidCurve * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

    if ( do_new || cmd == "build"  ) {

      Clothoid::valueType x0(0), y0(0), theta0(0), k0(0), dk(0), smin(0), smax(0);

      if ( nrhs == 2 ) {

        MEX_ASSERT( mxIsStruct(arg_in_1), "Second argument must be a struct" ) ;

        mxArray * mx_x0     = mxGetField( arg_in_1, 0, "x0" ) ;
        mxArray * mx_y0     = mxGetField( arg_in_1, 0, "y0" ) ;
        mxArray * mx_theta0 = mxGetField( arg_in_1, 0, "theta0" ) ;
        mxArray * mx_k0     = mxGetField( arg_in_1, 0, "k0" ) ;
        mxArray * mx_dk     = mxGetField( arg_in_1, 0, "dk" ) ;

        MEX_ASSERT( mx_x0     != nullptr, "Field `x0` is missing" );
        MEX_ASSERT( mx_y0     != nullptr, "Field `y0` is missing" );
        MEX_ASSERT( mx_theta0 != nullptr, "Field `theta0` is missing" );
        MEX_ASSERT( mx_k0     != nullptr, "Field `k0` is missing" );
        MEX_ASSERT( mx_dk     != nullptr, "Field `dk` is missing" );

        x0     = getScalarValue( mx_x0,     "Field `x0` must be a real double scalar" ) ;
        y0     = getScalarValue( mx_y0,     "Field `y0` must be a real double scalar" ) ;
        theta0 = getScalarValue( mx_theta0, "Field `theta0` must be a real double scalar" ) ;
        k0     = getScalarValue( mx_k0,     "Field `k0` must be a real double scalar" ) ;
        dk     = getScalarValue( mx_dk,     "Field `dk` must be a real double scalar" ) ;
      } else if ( 7 <= nrhs && nrhs <= 8 ) {
        x0     = getScalarValue(arg_in_1,"Error in reading x0") ;
        y0     = getScalarValue(arg_in_2,"Error in reading y0") ;
        theta0 = getScalarValue(arg_in_3,"Error in reading theta0") ;
        k0     = getScalarValue(arg_in_4,"Error in reading k0") ;
        dk     = getScalarValue(arg_in_5,"Error in reading dk") ;
        if ( nrhs == 8 ) {
          smin = getScalarValue(arg_in_6,"Error in reading smin") ;
          smax = getScalarValue(arg_in_7,"Error in reading smax") ;
        } else {  // the curve starts at curvilinear coordinate 0 and ends at L
          smin = 0;
          smax = getScalarValue(arg_in_6,"Error in reading L") ;
        }
      } else {
        MEX_ASSERT( false, "expected 2, 7 or 8 inputs") ;
      }

      ptr->build( x0, y0, theta0, k0, dk, smin, smax );

    } else if ( cmd == "build_G1" ) {

      MEX_ASSERT( nrhs == 8 , "expected 8 inputs") ;

      Clothoid::valueType x0(0), y0(0), theta0(0), x1(0), y1(0), theta1(0);

      x0     = getScalarValue( arg_in_2, "Error in reading x0" ) ;
      y0     = getScalarValue( arg_in_3, "Error in reading y0" ) ;
      theta0 = getScalarValue( arg_in_4, "Error in reading theta0" ) ;
      x1     = getScalarValue( arg_in_5, "Error in reading x1" ) ;
      y1     = getScalarValue( arg_in_6, "Error in reading y1" ) ;
      theta1 = getScalarValue( arg_in_7, "Error in reading theta1" ) ;

      ptr->build_G1(x0, y0, theta0, x1, y1, theta1);

    } else if ( cmd == "build_forward" ) {

      MEX_ASSERT( nrhs == 8 , "expected 8 inputs") ;

      Clothoid::valueType x0(0), y0(0), theta0(0), kappa0(0), x1(0), y1(0);

      x0     = getScalarValue( arg_in_2, "Error in reading x0" ) ;
      y0     = getScalarValue( arg_in_3, "Error in reading y0" ) ;
      theta0 = getScalarValue( arg_in_4, "Error in reading theta0" ) ;
      kappa0 = getScalarValue( arg_in_5, "Error in reading kappa0" ) ;
      x1     = getScalarValue( arg_in_6, "Error in reading x1" ) ;
      y1     = getScalarValue( arg_in_6, "Error in reading y1" ) ;

      bool res = ptr->build_forward(x0, y0, theta0, kappa0, x1, y1);

      // returns the status of the interpolation
      mwSize dims[2] = {1,1} ;
      arg_out_0 = mxCreateLogicalArray(2, dims);
      ((bool*)mxGetPr(arg_out_0))[0] = res;
      //mexPrintf("Result: %d\n", res);

    } else if ( cmd == "evaluate" ) {

      MEX_ASSERT( nrhs == 3 , "expected 3 inputs") ;

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
        MEX_ASSERT( false, "expected 2 or 4 outputs") ;
      }

    } else if ( cmd == "eval"    || cmd == "eval_D" ||
                cmd == "eval_DD" || cmd == "eval_DDD" ) {

      MEX_ASSERT( nrhs == 3 || nrhs == 4, "expected 3 or 4 inputs") ;
      MEX_ASSERT( nlhs == 2, "expected 2 outputs") ;

      mwSize size;
      double const * sVals = getVectorPointer( arg_in_2, size, "Error in reading s" );

      double offs = 0 ;
      if ( nrhs == 4 ) offs = getScalarValue( arg_in_3, "Error in reading offs" ) ;

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

      MEX_ASSERT( nrhs == 5 , "expected 5 inputs") ;
      MEX_ASSERT( nlhs == 4 , "expected 4 outputs") ;

      Clothoid::valueType x, y, ds, X, Y, S ;

      x  = getScalarValue(arg_in_2,"Error in reading x") ;
      y  = getScalarValue(arg_in_3,"Error in reading y") ;
      ds = getScalarValue(arg_in_4,"Error in reading ds") ;

      Clothoid::valueType DST = ptr->closestPoint( x, y, ds, X, Y, S );

      if ( nlhs > 0 ) setScalarValue(arg_out_0,S) ;
      if ( nlhs > 1 ) setScalarValue(arg_out_1,X) ;
      if ( nlhs > 2 ) setScalarValue(arg_out_2,Y) ;
      if ( nlhs > 3 ) setScalarValue(arg_out_3,DST) ;

    } else if ( cmd == "getPars" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      if ( nlhs > 0 ) setScalarValue(arg_out_0,ptr->getX0()) ;
      if ( nlhs > 1 ) setScalarValue(arg_out_1,ptr->getY0()) ;
      if ( nlhs > 2 ) setScalarValue(arg_out_2,ptr->getTheta0()) ;
      if ( nlhs > 3 ) setScalarValue(arg_out_3,ptr->getKappa()) ;
      if ( nlhs > 4 ) setScalarValue(arg_out_3,ptr->getKappa_D()) ;
      if ( nlhs > 5 ) setScalarValue(arg_out_3,ptr->getSmin()) ;
      if ( nlhs > 6 ) setScalarValue(arg_out_3,ptr->getSmax()) ;

    } else if ( cmd == "getX0" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getX0());

    } else if ( cmd == "getY0" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getY0());

    } else if ( cmd == "getTheta0" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getTheta0());

    } else if ( cmd == "getKappa0" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getKappa());

    } else if ( cmd == "getKappa_D" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getKappa_D());

    } else if ( cmd == "getSmin" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getSmin());

    } else if ( cmd == "getSmax" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getSmax());

    } else if ( cmd == "length" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 1, "expected 1 outputs");
      setScalarValue(arg_out_0, ptr->getSmax()-ptr->getSmin());

    } else if ( cmd == "trim" ) {

      MEX_ASSERT(nrhs == 4, "expected 4 inputs");

      Clothoid::valueType smin, smax;

      smin = getScalarValue(arg_in_2,"Error in reading smin") ;
      smax = getScalarValue(arg_in_3,"Error in reading smax") ;

      ptr->trim(smin, smax);

    } else if ( cmd == "changeOrigin" ) {

      MEX_ASSERT(nrhs == 3, "expected 3 inputs");

      Clothoid::valueType s0 = getScalarValue(arg_in_2,"Error in reading s0") ;
      ptr->change_origin(s0);

    } else if ( cmd == "rotate" ) {

      MEX_ASSERT(nrhs == 5, "expected 5 inputs");

      Clothoid::valueType angle = getScalarValue( arg_in_2, "Error in reading angle" ) ;
      Clothoid::valueType cx    = getScalarValue( arg_in_3, "Error in reading cx" ) ;
      Clothoid::valueType cy    = getScalarValue( arg_in_4, "Error in reading cy" ) ;
      ptr->rotate(angle, cx, cy);

    } else if ( cmd == "translate" ) {

      MEX_ASSERT(nrhs == 4, "expected 4 inputs");
      Clothoid::valueType tx = getScalarValue( arg_in_2, "Error in reading tx" ) ;
      Clothoid::valueType ty = getScalarValue( arg_in_3, "Error in reading ty" ) ;
      ptr->translate(tx, ty);

    } else if ( cmd == "moveOrigin" ) {

      MEX_ASSERT(nrhs == 4, "expected 4 inputs");
      Clothoid::valueType newX0 = getScalarValue( arg_in_2, "Error in reading newX0" ) ;
      Clothoid::valueType newY0 = getScalarValue( arg_in_3, "Error in reading newY0" ) ;
      ptr->moveOrigin(newX0, newY0);

    } else if ( cmd == "scale" ) {

      MEX_ASSERT(nrhs == 3, "expected 3 inputs");
      Clothoid::valueType s = getScalarValue( arg_in_2, "Error in reading s" ) ;
      ptr->scale(s);

    } else if ( cmd == "reverse" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      ptr->reverse();

    } else if ( cmd == "delete" ) {

      MEX_ASSERT(nrhs == 2, "expected 2 inputs");
      MEX_ASSERT(nlhs == 0, "expected no output");

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
