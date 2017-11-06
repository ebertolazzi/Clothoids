/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex_class_handle.hh"
#include "mex.h"

#include <sstream>
#include <stdexcept>

#define MEX_ERROR_MESSAGE \
"%=====================================================================================%\n" \
"%  ClothoidMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting         %\n" \
"%                                                                                     %\n" \
"%  USAGE:                                                                             %\n" \
"%  - Constructors:                                                                    %\n" \
"%      OBJ = ClothoidMexWrapper( 'new' ) ;                                            %\n" \
"%      OBJ = ClothoidMexWrapper( 'new', x0, y0, theta0, k0, dk, L ) ;                 %\n" \
"%      OBJ = ClothoidMexWrapper( 'new', x0, y0, theta0, k0, dk, smin, smax ) ;        %\n" \
"%                                                                                     %\n" \
"%      On input:                                                                      %\n" \
"%        x0, y0     = coordinate of initial point                                     %\n" \
"%        theta0     = orientation (angle) of the clothoid at initial point            %\n" \
"%        k0         = curvature of the clothoid at initial point                      %\n" \
"%        dk         = derivative of curvature respect to arclength                    %\n" \
"%        L          = length of the clothoid curve from initial to final point        %\n" \
"%        smin       = initial curvilinear coordinate of the curve                     %\n" \
"%        smax       = final curvilinear coordinate of the curve                       %\n" \
"%                                                                                     %\n" \
"%      On output:                                                                     %\n" \
"%        OBJ     = pointer to the internal object                                     %\n" \
"%                                                                                     %\n" \
"%  - Destructor:                                                                      %\n" \
"%      ClothoidMexWrapper( 'delete', OBJ ) ;                                          %\n" \
"%                                                                                     %\n" \
"%  - Build                                                                            %\n" \
"%      ClothoidMexWrapper( 'build', OBJ, x0, y0, theta0, k0, dk, L ) ;                %\n" \
"%      ClothoidMexWrapper( 'build', OBJ, x0, y0, theta0, k0, dk, smin, smax ) ;       %\n" \
"%      ClothoidMexWrapper( 'build_G1', OBJ, x0, y0, theta0, x1, y1, theta1 ) ;        %\n" \
"%      res = ClothoidMexWrapper( 'build_forward', OBJ, x0, y0, theta0, k0, x1, y1 ) ; %\n" \
"%                                                                                     %\n" \
"%  - Eval                                                                             %\n" \
"%      [x,y] = ClothoidMexWrapper( 'eval', OBJ, ss ) ;                                %\n" \
"%      [x,y,theta,kappa] = ClothoidMexWrapper( 'eval', OBJ, ss ) ;                    %\n" \
"%      [x0,y0,theta0,k0,dk,smin,smax] = ClothoidMexWrapper( 'getPars', OBJ ) ;        %\n" \
"%                                                                                     %\n" \
"%  - Transform                                                                        %\n" \
"%      ClothoidMexWrapper( 'trim', OBJ, smin, smax ) ;                                %\n" \
"%      ClothoidMexWrapper( 'changeOrigin', OBJ, s0 ) ;                                %\n" \
"%      ClothoidMexWrapper( 'rotate', OBJ, angle, cx, cy ) ;                           %\n" \
"%      ClothoidMexWrapper( 'translate', OBJ, tx, ty ) ;                               %\n" \
"%      ClothoidMexWrapper( 'moveOrigin', OBJ, newX0, newY0 ) ;                        %\n" \
"%      ClothoidMexWrapper( 'scale', OBJ, scaling ) ;                                  %\n" \
"%      ClothoidMexWrapper( 'reverse', OBJ ) ;                                         %\n" \
"%                                                                                     %\n" \
"%=====================================================================================%\n" \
"%                                                                                     %\n" \
"%  Autors: Enrico Bertolazzi and Marco Frego and Paolo Bevilacqua                     %\n" \
"%          Department of Industrial Engineering                                       %\n" \
"%          Department of Information Engineering and Computer Science                 %\n" \
"%          University of Trento                                                       %\n" \
"%          enrico.bertolazzi@unitn.it                                                 %\n" \
"%          m.fregox@gmail.com                                                         %\n" \
"%          paolo.bevilacqua@unitn.it                                                  %\n" \
"%                                                                                     %\n" \
"%=====================================================================================%\n"

#define ASSERT(COND,MSG)                           \
  if ( !(COND) ) {                                 \
    std::ostringstream ost ;                       \
    ost << "ClothoidMexWrapper: " << MSG << '\n' ; \
    ost << MEX_ERROR_MESSAGE ;                     \
    mexErrMsgTxt(ost.str().c_str()) ;              \
  }

#define arg_in_0 prhs[0]
#define arg_in_1 prhs[1]
#define arg_in_2 prhs[2]
#define arg_in_3 prhs[3]
#define arg_in_4 prhs[4]
#define arg_in_5 prhs[5]
#define arg_in_6 prhs[6]
#define arg_in_7 prhs[7]
#define arg_in_8 prhs[8]

#define arg_out_0 plhs[0]
#define arg_out_1 plhs[1]
#define arg_out_2 plhs[2]
#define arg_out_3 plhs[3]
#define arg_out_4 plhs[4]
#define arg_out_5 plhs[5]
#define arg_out_6 plhs[6]
#define arg_out_7 plhs[7]

static
double
getScalarValue( mxArray const * arg, char const * msg ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  ASSERT( dims[0] == 1 && dims[1] == 1,
          msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  return mxGetScalar(arg) ;
}

static
double *
getArrayValues( mxArray const * arg, int & size, char const * msg) {
    mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
    ASSERT( number_of_dimensions == 2, msg ) ;
    mwSize const * dims = mxGetDimensions(arg) ;
    ASSERT( dims[0] == 1 || dims[1] == 1,
            msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
    size = mxGetNumberOfElements(arg);
    return mxGetPr(arg);
}

static
void
setScalarValue( mxArray * & arg, double const & value ) {
    arg = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * pA = mxGetPr(arg);
    *pA = value;
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

    ASSERT( mxIsChar(arg_in_0), "First argument must be a string" ) ;
    std::string cmd = mxArrayToString(arg_in_0) ;


    if ( cmd == "new" ) {

      ASSERT( nrhs == 1 || nrhs == 7 || nrhs == 8, "expected 1, 7 or 8 inputs") ;

      Clothoid::valueType x0, y0, theta0, k0, dk, smin, smax;

      if (nrhs>=7) {
        x0     = getScalarValue(arg_in_1,"Error in reading x0") ;
        y0     = getScalarValue(arg_in_2,"Error in reading y0") ;
        theta0 = getScalarValue(arg_in_3,"Error in reading theta0") ;
        k0     = getScalarValue(arg_in_4,"Error in reading k0") ;
        dk     = getScalarValue(arg_in_5,"Error in reading dk") ;
        if (nrhs==8) {
          smin = getScalarValue(arg_in_6,"Error in reading smin") ;
          smax = getScalarValue(arg_in_7,"Error in reading smax") ;
        }
        else {  // the curve starts at curvilinear coordinate 0 and ends at L
          smin = 0;
          smax = getScalarValue(arg_in_6,"Error in reading L") ;
        }
      }

      Clothoid::ClothoidCurve * ptr = new Clothoid::ClothoidCurve(x0, y0, theta0, k0, dk, smin, smax);
      plhs[0] = convertPtr2Mat<Clothoid::ClothoidCurve>(ptr);

    } else if ( cmd == "build" ) {

        ASSERT( nrhs == 8 || nrhs == 9, "expected 8 or 9 inputs") ;

        Clothoid::valueType x0, y0, theta0, k0, dk, smin, smax;

        x0     = getScalarValue(arg_in_2,"Error in reading x0") ;
        y0     = getScalarValue(arg_in_3,"Error in reading y0") ;
        theta0 = getScalarValue(arg_in_4,"Error in reading theta0") ;
        k0     = getScalarValue(arg_in_5,"Error in reading k0") ;
        dk     = getScalarValue(arg_in_6,"Error in reading dk") ;
        if (nrhs==9) {
            smin = getScalarValue(arg_in_7,"Error in reading smin") ;
            smax = getScalarValue(arg_in_8,"Error in reading smax") ;
        }
        else {  // the curve starts at curvilinear coordinate 0 and ends at L
            smin = 0;
            smax = getScalarValue(arg_in_7,"Error in reading L") ;
        }

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->build(x0, y0, theta0, k0, dk, smin, smax);

    } else if ( cmd == "build_G1" ) {

        ASSERT( nrhs == 8 , "expected 8 inputs") ;

        Clothoid::valueType x0, y0, theta0, x1, y1, theta1;

        x0     = getScalarValue(arg_in_2,"Error in reading x0") ;
        y0     = getScalarValue(arg_in_3,"Error in reading y0") ;
        theta0 = getScalarValue(arg_in_4,"Error in reading theta0") ;
        x1     = getScalarValue(arg_in_5,"Error in reading x1") ;
        y1     = getScalarValue(arg_in_6,"Error in reading y1") ;
        theta1 = getScalarValue(arg_in_7,"Error in reading theta1") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->build_G1(x0, y0, theta0, x1, y1, theta1);

    } else if ( cmd == "build_forward" ) {

        ASSERT( nrhs == 8 , "expected 8 inputs") ;

        Clothoid::valueType x0, y0, theta0, kappa0, x1, y1;

        x0     = getScalarValue(arg_in_2,"Error in reading x0") ;
        y0     = getScalarValue(arg_in_3,"Error in reading y0") ;
        theta0 = getScalarValue(arg_in_4,"Error in reading theta0") ;
        kappa0 = getScalarValue(arg_in_5,"Error in reading kappa0") ;
        x1     = getScalarValue(arg_in_6,"Error in reading x1") ;
        y1     = getScalarValue(arg_in_6,"Error in reading y1") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        bool res = ptr->build_forward(x0, y0, theta0, kappa0, x1, y1);

        // returns the status of the interpolation
        mwSize dims[2] = {1,1} ;
        arg_out_0 = mxCreateLogicalArray(2, dims);
        ((bool*)mxGetPr(arg_out_0))[0] = res;
        mexPrintf("Result: %d\n", res);

    } else if ( cmd == "eval" ) {

        ASSERT( nrhs == 3 , "expected 3 inputs") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);

        int size;
        double * sVals = getArrayValues(arg_in_2,size,"Error in reading s");

        bool fullData = nlhs==4; // output also theta(s) and kappa(s)

        double *xVals, *yVals, *thetaVals, *kappaVals;

        arg_out_0 = mxCreateDoubleMatrix(size, 1, mxREAL);
        xVals = mxGetPr(arg_out_0);
        arg_out_1 = mxCreateDoubleMatrix(size, 1, mxREAL);
        yVals = mxGetPr(arg_out_1);
        if (fullData) {
            arg_out_2 = mxCreateDoubleMatrix(size, 1, mxREAL);
            thetaVals = mxGetPr(arg_out_2);
            arg_out_3 = mxCreateDoubleMatrix(size, 1, mxREAL);
            kappaVals = mxGetPr(arg_out_3);
        }

        for (int i=0; i<size; ++i) {
            if (fullData)
                ptr->eval(sVals[i], thetaVals[i], kappaVals[i], xVals[i], yVals[i]);
            else
                ptr->eval(sVals[i], xVals[i], yVals[i]);
        }

    } else if ( cmd == "getPars" ) {

        ASSERT(nrhs == 2, "expected 2 inputs");
        ASSERT(nlhs == 7, "expected 7 outputs");

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        setScalarValue(arg_out_0, ptr->getX0());      // x0
        setScalarValue(arg_out_1, ptr->getY0());      // y0
        setScalarValue(arg_out_2, ptr->getTheta0());  // theta0
        setScalarValue(arg_out_3, ptr->getKappa());   // kappa0
        setScalarValue(arg_out_4, ptr->getKappa_D()); // kappa_D
        setScalarValue(arg_out_5, ptr->getSmin());    // smin
        setScalarValue(arg_out_6, ptr->getSmax());    // smax

    } else if ( cmd == "trim" ) {

        ASSERT(nrhs == 4, "expected 4 inputs");

        Clothoid::valueType smin, smax;

        smin = getScalarValue(arg_in_2,"Error in reading smin") ;
        smax = getScalarValue(arg_in_3,"Error in reading smax") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->trim(smin, smax);

    } else if ( cmd == "changeOrigin" ) {

        ASSERT(nrhs == 3, "expected 3 inputs");

        Clothoid::valueType s0;

        s0 = getScalarValue(arg_in_2,"Error in reading s0") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->change_origin(s0);

    } else if ( cmd == "rotate" ) {

        ASSERT(nrhs == 5, "expected 5 inputs");

        Clothoid::valueType angle, cx, cy;

        angle = getScalarValue(arg_in_2,"Error in reading angle") ;
        cx    = getScalarValue(arg_in_3,"Error in reading cx") ;
        cy    = getScalarValue(arg_in_4,"Error in reading cy") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->rotate(angle, cx, cy);

    } else if ( cmd == "translate" ) {

        ASSERT(nrhs == 4, "expected 4 inputs");

        Clothoid::valueType tx, ty;

        tx = getScalarValue(arg_in_2,"Error in reading tx") ;
        ty = getScalarValue(arg_in_3,"Error in reading ty") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->translate(tx, ty);

    } else if ( cmd == "moveOrigin" ) {

        ASSERT(nrhs == 4, "expected 4 inputs");

        Clothoid::valueType newX0, newY0;

        newX0 = getScalarValue(arg_in_2,"Error in reading newX0") ;
        newY0 = getScalarValue(arg_in_3,"Error in reading newY0") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->moveOrigin(newX0, newY0);

    } else if ( cmd == "scale" ) {

        ASSERT(nrhs == 3, "expected 3 inputs");

        Clothoid::valueType s;

        s = getScalarValue(arg_in_2,"Error in reading s") ;

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->scale(s);

    } else if ( cmd == "reverse" ) {

        ASSERT(nrhs == 2, "expected 2 inputs");

        Clothoid::ClothoidCurve * ptr = convertMat2Ptr<Clothoid::ClothoidCurve>(arg_in_1);
        ptr->reverse();

    } else if ( cmd == "delete" ) {

        ASSERT(nrhs == 2, "expected 2 inputs");
        ASSERT(nlhs == 0, "expected no output");

        // Destroy the C++ object
        destroyObject<Clothoid::ClothoidCurve>(prhs[1]);

        // Warn if other commands were ignored
    } else {

    }

  } catch ( std::exception const & e ) {
  	mexErrMsgTxt(e.what()) ;

  } catch (...) {
  	mexErrMsgTxt("clothoid failed\n") ;
  }
}
