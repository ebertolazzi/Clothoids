/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "ClothoidList.hh"
#include "Triangle2D.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"ClothoidSplineG2MexWrapper:  Compute parameters of the G2 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = ClothoidSplineG2MexWrapper( 'new' );\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
"n" \
"  - Destructor:\n" \
"    ClothoidSplineG2MexWrapper( 'delete', OBJ );\n" \
"\n" \
"  - Build:\n" \
"    ClothoidSplineG2MexWrapper('build',OBJ,x,y);\n" \
"    ClothoidSplineG2MexWrapper('target',OBJ,target[,theta0,theta1]);\n" \
"\n" \
"  - Eval:\n" \
"    theta_guess = ClothoidSplineG2MexWrapper('guess',OBJ);\n" \
"    obj         = ClothoidSplineG2MexWrapper('objective',OBJ,theta);\n" \
"    g           = ClothoidSplineG2MexWrapper('gradient',OBJ,theta);\n" \
"    c           = ClothoidSplineG2MexWrapper('constraints',OBJ,theta);\n" \
"    jac         = ClothoidSplineG2MexWrapper('jacobian',OBJ,theta);\n" \
"    jac_pattern = ClothoidSplineG2MexWrapper('jacobian_pattern',OBJ);\n" \
"    [n,nc]      = ClothoidSplineG2MexWrapper('dims',OBJ);\n" \
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
    return ptr;
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
  mexFunction(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {

      MEX_ASSERT( mxIsChar(arg_in_0), "ClothoidListMexWrapper: First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);

      bool do_new = cmd == "new";
      ClothoidSplineG2 * ptr = do_new ? DATA_NEW(arg_out_0) : DATA_GET(arg_in_1);

      if ( do_new ) {

        MEX_ASSERT( nlhs == 1, "ClothoidSplineG2MexWrapper, expected 1 output, nlhs = " << nlhs  );

      } else if ( cmd == "build" ) {

        #define CMD "ClothoidSplineG2MexWrapper('build',OBJ,x,y): "

        MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        mwSize nx, ny;
        real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
        real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

        MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " must be equal to size(y) = " << ny );

        ptr->build( x, y, nx );

        #undef CMD

      } else if ( cmd == "target" ) {

        #define CMD "ClothoidSplineG2MexWrapper('target',OBJ,target[,theta0,theta1]): "

        MEX_ASSERT( nrhs >= 3, CMD "expected at least 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );

        MEX_ASSERT( mxIsChar(arg_in_2), CMD "Third argument must be a string" );
        string obj = mxArrayToString(arg_in_2);

        if ( obj == "P1" ) {
          MEX_ASSERT( nrhs == 5, CMD "expected at 5 inputs, nrhs = " << nrhs );
          real_type theta0 = getScalarValue( arg_in_3, CMD "Error in reading theta0" );
          real_type theta1 = getScalarValue( arg_in_4, CMD "Error in reading theta1" );
          ptr->setP1( theta0, theta1 );
        } else {
          MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
          if      ( obj == "P2" ) ptr->setP2();
          else if ( obj == "P3" ) ptr->setP3();
          else if ( obj == "P4" ) ptr->setP4();
          else if ( obj == "P5" ) ptr->setP5();
          else if ( obj == "P6" ) ptr->setP6();
          else if ( obj == "P7" ) ptr->setP7();
          else if ( obj == "P8" ) ptr->setP8();
          else if ( obj == "P9" ) ptr->setP9();
          else {
            MEX_ASSERT( false, CMD "Unknown target " << obj );
          }
        }
        #undef CMD

      } else if ( cmd == "guess" ) {

        #define CMD "ClothoidSplineG2MexWrapper('guess',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = " << nlhs );

        int_type N = ptr->numPnts();

        real_type * theta_guess = createMatrixValue( arg_out_0, N, 1 );
        real_type * theta_min   = createMatrixValue( arg_out_1, N, 1 );
        real_type * theta_max   = createMatrixValue( arg_out_2, N, 1 );

        ptr->guess( theta_guess, theta_min, theta_max );

        #undef CMD

      } else if ( cmd == "objective" ) {

        #define CMD "ClothoidSplineG2MexWrapper('objective',OBJ,theta): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        mwSize ntheta;
        real_type const * theta = getVectorPointer( arg_in_2, ntheta, CMD "Error in reading theta" );
        MEX_ASSERT( ntheta == ptr->numPnts(),
                    CMD "length(theta) = " << ntheta << " must be " << ptr->numPnts() );
        real_type f;
        bool ok = ptr->objective( theta, f );
        if ( ok ) setScalarValue( arg_out_0, f );

        #undef CMD

      } else if ( cmd == "gradient" ) {

        #define CMD "ClothoidSplineG2MexWrapper('gradient',OBJ,theta): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        mwSize ntheta;
        real_type const * theta = getVectorPointer( arg_in_2, ntheta, CMD "Error in reading theta" );
        MEX_ASSERT(
          ntheta == ptr->numPnts(),
          CMD "length(theta) = " << ntheta << " must be " << ptr->numPnts()
        );
        double * g = createMatrixValue( arg_out_0, ntheta, 1 );
        bool ok = ptr->gradient( theta, g );
        MEX_ASSERT( ok, CMD "bad gradient computation");

        #undef CMD

      } else if ( cmd == "constraints" ) {

        #define CMD "ClothoidSplineG2MexWrapper('gradient',OBJ,theta): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        mwSize ntheta;
        real_type const * theta = getVectorPointer( arg_in_2, ntheta, CMD "Error in reading theta" );
        MEX_ASSERT(
          ntheta == ptr->numPnts(),
          CMD "length(theta) = " << ntheta << " must be " << ptr->numPnts()
        );
        double * c = createMatrixValue( arg_out_0, ptr->numConstraints(), 1 );
        bool ok = ptr->constraints( theta, c );
        MEX_ASSERT( ok, CMD "bad constraints computation");

        #undef CMD

      } else if ( cmd == "jacobian" ) {

        #define CMD "ClothoidSplineG2MexWrapper('jacobian',OBJ,theta): "

        MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        mwSize ntheta;
        real_type const * theta = getVectorPointer( arg_in_2, ntheta, CMD "Error in reading theta" );
        MEX_ASSERT(
          ntheta == ptr->numPnts(),
          CMD "length(theta) = " << ntheta << " must be " << ptr->numPnts()
        );

        int_type n   = ptr->numConstraints();
        int_type m   = ptr->numTheta();
        int_type nnz = ptr->jacobian_nnz();

        mxArray *args[5];

        real_type * I = createMatrixValue( args[0], 1, nnz );
        real_type * J = createMatrixValue( args[1], 1, nnz );
        real_type * V = createMatrixValue( args[2], 1, nnz );
        setScalarValue( args[3], n );
        setScalarValue( args[4], m );

        ptr->jacobian_pattern_matlab( I, J );
        ptr->jacobian( theta, V );

        int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );
        MEX_ASSERT( ok == 0, CMD "failed the call sparse(...)" );

        #undef CMD

      } else if ( cmd == "jacobian_pattern" ) {

        #define CMD "ClothoidSplineG2MexWrapper('jacobian_pattern',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = " << nlhs );

        int_type n   = ptr->numConstraints();
        int_type m   = ptr->numTheta();
        int_type nnz = ptr->jacobian_nnz();

        mxArray *args[5];

        real_type * I = createMatrixValue( args[0], 1, nnz );
        real_type * J = createMatrixValue( args[1], 1, nnz );
        real_type * V = createMatrixValue( args[2], 1, nnz );
        setScalarValue( args[3], n );
        setScalarValue( args[4], m );

        ptr->jacobian_pattern_matlab( I, J );
        std::fill( V, V+nnz, 1 );

        int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );
        MEX_ASSERT( ok == 0, CMD "failed the call sparse(...)" );

        #undef CMD

      } else if ( cmd == "dims" ) {

        #define CMD "ClothoidSplineG2MexWrapper('dims',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );

        int_type m = ptr->numTheta();
        int_type n = ptr->numConstraints();

        setScalarInt( arg_out_0, m );
        setScalarInt( arg_out_1, n );

        #undef CMD

      } else if ( cmd == "info" ) {

        #define CMD "ClothoidSplineG2MexWrapper('info',OBJ): "

        MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
        MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = " << nlhs );

        ptr->info(cout);

        #undef CMD

      } else if ( cmd == "delete" ) {

        #define CMD "ClothoidListMexWrapper('delete',OBJ): "

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
  	  mexErrMsgTxt("clothoid failed\n");
    }
  }
}
