/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2023
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
"DubinsMexWrapper:  Compute solution of Dubins problem\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = DubinsMexWrapper( 'new' );\n" \
"\n" \
"    On output:\n" \
"    OBJ = pointer to the internal object\n" \
"\n" \
"  - Build:\n" \
"    [arc0,arc1,arc2] = DubinsMexWrapper( 'build', OBJ, x0, y0, theta0, x1, y1, theta1, k_max );\n" \
"\n" \
"  - Eval:\n" \
"    [x0,y0,theta0,kappa0,...\n" \
"     L0,x1,y1,theta1,kappa1,...\n" \
"     L1,x2,y2,theta2,kappa2,L2] = DubinsMexWrapper( 'get_pars', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("DubinsMexWrapper") \
MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  /*\
   |  ____    _  _____  _
   | |  _ \  / \|_   _|/ \
   | | | | |/ _ \ | | / _ \
   | | |_| / ___ \| |/ ___ \
   | |____/_/   \_\_/_/   \_\
   |
  \*/

  /*\
   *                      _____                 _   _
   *  _ __ ___   _____  _|  ___|   _ _ __   ___| |_(_) ___  _ __
   * | '_ ` _ \ / _ \ \/ / |_ | | | | '_ \ / __| __| |/ _ \| '_ \
   * | | | | | |  __/>  <|  _|| |_| | | | | (__| |_| | (_) | | | |
   * |_| |_| |_|\___/_/\_\_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
   *
  \*/

  #define CMD_BASE "DubinsMexWrapper"

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *[] ) {

    #define CMD CMD_BASE "('new'): "
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    #undef CMD

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Dubins>(new Dubins());
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD CMD_BASE "('delete',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );
    // Destroy the C++ object
    Utils::mex_destroy_object<Dubins>( arg_in_1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_a_copy(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD CMD_BASE "('make_a_copy',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins const * ptr = Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1);
    arg_out_0 = Utils::mex_convert_ptr_to_mx<Dubins>(new Dubins( *ptr ));
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1,k_max): "

    UTILS_MEX_ASSERT( nrhs == 9, CMD "expected 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr = Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1);

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading x1" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading y1" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading theta1" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading k_max" );

    bool ok = ptr->build( x0, y0, theta0, x1, y1, theta1, k_max );

    // returns the status of the interpolation
    Utils::mex_set_scalar_bool( arg_out_0, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_pars( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Dubins * ptr = Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1);

    #define CMD "DubinsMexWrapper('get_pars',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    CircleArc const & C0 = ptr->C0();
    CircleArc const & C1 = ptr->C1();
    CircleArc const & C2 = ptr->C2();

    static char const * fieldnames[] = {
      "x0", "y0", "theta0", "kappa0", "L0",
      "x1", "y1", "theta1", "kappa1", "L1",
      "x2", "y2", "theta2", "kappa2", "L2"
    };

    arg_out_0 = mxCreateStructMatrix(1,1,15,fieldnames);

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateDoubleScalar(C0.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(C0.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(C0.theta_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(C0.kappa_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mxCreateDoubleScalar(C0.length()) );

    mxSetFieldByNumber( arg_out_0, 0, 5, mxCreateDoubleScalar(C1.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 6, mxCreateDoubleScalar(C1.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 7, mxCreateDoubleScalar(C1.theta_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 8, mxCreateDoubleScalar(C1.kappa_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 9, mxCreateDoubleScalar(C1.length()) );

    mxSetFieldByNumber( arg_out_0, 0, 10, mxCreateDoubleScalar(C2.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 11, mxCreateDoubleScalar(C2.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 12, mxCreateDoubleScalar(C2.theta_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 13, mxCreateDoubleScalar(C2.kappa_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 14, mxCreateDoubleScalar(C2.length()) );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"delete",do_delete},
    {"make_a_copy",do_make_a_copy},
    {"build",do_build},
    {"get_pars",do_get_pars}
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  extern "C"
  void
  mexFunction( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    char cmd[256];

    // the first argument must be a string
    if ( nrhs == 0 ) { mexErrMsgTxt(MEX_ERROR_MESSAGE); return; }

    try {
      UTILS_MEX_ASSERT0( mxIsChar(arg_in_0), "First argument must be a string" );
      mxGetString( arg_in_0, cmd, 256 );
      cmd_to_fun.at(cmd)( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "Dubins Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Dubins failed\n");
    }

  }
}
