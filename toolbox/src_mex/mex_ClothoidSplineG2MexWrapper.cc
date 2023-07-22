/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

#ifdef _MSC_VER
  #pragma comment(lib, "IPHLPAPI.lib")
  #pragma comment(lib, "ws2_32.lib")
  #pragma comment(lib, "Shlwapi.lib")
  #pragma comment(lib, "Advapi32.lib")
  #pragma comment(lib, "Shell32.lib")
  #pragma comment(lib, "kernel32.lib")
#endif

#include "Clothoids.hh"
#include "Utils_mex.hh"
#include "mex_info.hxx"

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

#include <unordered_map>

namespace G2lib {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *[]
  ) {
    #define CMD "ClothoidSplineG2MexWrapper('new'): "
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nlhs = {}\n", nrhs );
    arg_out_0 = Utils::mex_convert_ptr_to_mx<ClothoidSplineG2>(new ClothoidSplineG2());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('delete',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    // Destroy the C++ object
    Utils::mex_destroy_object<ClothoidSplineG2>(arg_in_1);

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidSplineG2MexWrapper('build',OBJ,x,y): "

    UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    mwSize nx, ny;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD "Error in reading x" );
    real_type const * y = Utils::mex_vector_pointer( arg_in_3, ny, CMD "Error in reading y" );

    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} must be equal to size(y) = {}\n", nx, ny );

    ptr->build( x, y, nx );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_target(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('target',OBJ,target[,theta0,theta1]): "

    UTILS_MEX_ASSERT( nrhs >= 3, CMD "expected at least 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), CMD "Third argument must be a string" );
    string obj = mxArrayToString(arg_in_2);

    if ( obj == "P1" ) {
      UTILS_MEX_ASSERT( nrhs == 5, CMD "expected at 5 inputs, nrhs = {}\n", nrhs );
      real_type theta0 = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading theta0" );
      real_type theta1 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta1" );
      ptr->setP1( theta0, theta1 );
    } else {
      UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
      if      ( obj == "P2" ) ptr->setP2();
      else if ( obj == "P3" ) ptr->setP3();
      else if ( obj == "P4" ) ptr->setP4();
      else if ( obj == "P5" ) ptr->setP5();
      else if ( obj == "P6" ) ptr->setP6();
      else if ( obj == "P7" ) ptr->setP7();
      else if ( obj == "P8" ) ptr->setP8();
      else if ( obj == "P9" ) ptr->setP9();
      else {
        UTILS_MEX_ASSERT( false, CMD "Unknown target {}\n", obj );
      }
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_guess(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('guess',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    integer N = ptr->numPnts();

    real_type * theta_guess = Utils::mex_create_matrix_value( arg_out_0, N, 1 );
    real_type * theta_min   = Utils::mex_create_matrix_value( arg_out_1, N, 1 );
    real_type * theta_max   = Utils::mex_create_matrix_value( arg_out_2, N, 1 );

    ptr->guess( theta_guess, theta_min, theta_max );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_objective(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('objective',OBJ,theta): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    mwSize ntheta;
    real_type const * theta = Utils::mex_vector_pointer( arg_in_2, ntheta, CMD "Error in reading theta" );
    UTILS_MEX_ASSERT(
      ntheta == static_cast<mwSize>(ptr->numPnts()),
      CMD "length(theta) = {} must be {}\n", ntheta, ptr->numPnts()
    );
    real_type f;
    bool ok = ptr->objective( theta, f );
    if ( ok ) Utils::mex_set_scalar_value( arg_out_0, f );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_gradient(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('gradient',OBJ,theta): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    mwSize ntheta;
    real_type const * theta = Utils::mex_vector_pointer( arg_in_2, ntheta, CMD "Error in reading theta" );
    UTILS_MEX_ASSERT(
      ntheta == static_cast<mwSize>(ptr->numPnts()),
      CMD "length(theta) = {} must be {}\n", ntheta, ptr->numPnts()
    );
    double * g = Utils::mex_create_matrix_value( arg_out_0, ntheta, 1 );
    bool ok = ptr->gradient( theta, g );
    UTILS_MEX_ASSERT0( ok, CMD "bad gradient computation");

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_constraints(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('gradient',OBJ,theta): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    mwSize ntheta;
    real_type const * theta = Utils::mex_vector_pointer( arg_in_2, ntheta, CMD "Error in reading theta" );
    UTILS_MEX_ASSERT(
      ntheta == static_cast<mwSize>(ptr->numPnts()),
      CMD "length(theta) = {} must be {}\n", ntheta, ptr->numPnts()
    );
    double * c = Utils::mex_create_matrix_value( arg_out_0, ptr->numConstraints(), 1 );
    bool ok = ptr->constraints( theta, c );
    UTILS_MEX_ASSERT0( ok, CMD "bad constraints computation");

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_jacobian(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('jacobian',OBJ,theta): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    mwSize ntheta;
    real_type const * theta = Utils::mex_vector_pointer( arg_in_2, ntheta, CMD "Error in reading theta" );
    UTILS_MEX_ASSERT(
      ntheta == static_cast<mwSize>(ptr->numPnts()),
      CMD "length(theta) = {} must be {}\n", ntheta, ptr->numPnts()
    );

    integer n   = ptr->numConstraints();
    integer m   = ptr->numTheta();
    integer nnz = ptr->jacobian_nnz();

    mxArray *args[5];

    real_type * I = Utils::mex_create_matrix_value( args[0], 1, nnz );
    real_type * J = Utils::mex_create_matrix_value( args[1], 1, nnz );
    real_type * V = Utils::mex_create_matrix_value( args[2], 1, nnz );
    Utils::mex_set_scalar_value( args[3], n );
    Utils::mex_set_scalar_value( args[4], m );

    ptr->jacobian_pattern_matlab( I, J );
    ptr->jacobian( theta, V );

    int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );
    UTILS_MEX_ASSERT0( ok == 0, CMD "failed the call sparse(...)" );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_jacobian_pattern(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('jacobian_pattern',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    integer n   = ptr->numConstraints();
    integer m   = ptr->numTheta();
    integer nnz = ptr->jacobian_nnz();

    mxArray *args[5];

    real_type * I = Utils::mex_create_matrix_value( args[0], 1, nnz );
    real_type * J = Utils::mex_create_matrix_value( args[1], 1, nnz );
    real_type * V = Utils::mex_create_matrix_value( args[2], 1, nnz );
    Utils::mex_set_scalar_value( args[3], n );
    Utils::mex_set_scalar_value( args[4], m );

    ptr->jacobian_pattern_matlab( I, J );
    std::fill( V, V+nnz, 1 );

    int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );
    UTILS_MEX_ASSERT0( ok == 0, CMD "failed the call sparse(...)" );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_dims(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('dims',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    integer m = ptr->numTheta();
    integer n = ptr->numConstraints();

    Utils::mex_set_scalar_int32( arg_out_0, m );
    Utils::mex_set_scalar_int32( arg_out_1, n );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_info(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidSplineG2MexWrapper('info',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO outputs, nlhs = {}\n", nlhs );

    ClothoidSplineG2 * ptr = Utils::mex_convert_mx_to_ptr<ClothoidSplineG2>(arg_in_1);

    ptr->info( std::cout );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"delete",do_delete},
    {"build",do_build},
    {"target",do_target},
    {"guess",do_guess},
    {"objective",do_objective},
    {"gradient",do_gradient},
    {"constraints",do_constraints},
    {"jacobian",do_jacobian},
    {"jacobian_pattern",do_jacobian_pattern},
    {"dims",do_dims},
    {"info",do_info},
  };

  extern "C"
  void
  mexFunction(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    char cmd[256];

    // the first argument must be a string
    if ( nrhs == 0 ) { mexErrMsgTxt(MEX_ERROR_MESSAGE); return; }

    try {
      UTILS_MEX_ASSERT0(
        mxIsChar(arg_in_0),
        "ClothoidListMexWrapper: First argument must be a string"
      );
      mxGetString( arg_in_0, cmd, 256 );
      cmd_to_fun.at(cmd)( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "Clothoid Error: {}", e.what() ).c_str() );
    } catch (...) {
  	  mexErrMsgTxt("Clothoid failed\n");
    }
  }
}
