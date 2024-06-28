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
#include "Clothoids_fmt.hh"
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

  #define CMD_BASE "DubinsMexWrapper"
  #define G2LIB_CLASS Dubins
  #include "mex_common.hxx"
  //#undef CMD_BASE
  #undef G2LIB_CLASS

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

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Dubins>(new Dubins("dubins"));
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build( int nlhs, mxArray       *plhs[],
            int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1,k_max): "

    UTILS_MEX_ASSERT( nrhs == 9, CMD "expected 9 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

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

    #define CMD "DubinsMexWrapper('get_pars',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    CircleArc const & C0{ ptr->C0() };
    CircleArc const & C1{ ptr->C1() };
    CircleArc const & C2{ ptr->C2() };

    static char const * fieldnames[] = {
      "x0", "y0", "theta0",
      "x1", "y1", "theta1",
      "x2", "y2", "theta2",
      "x3", "y3", "theta3",
      "kappa1", "kappa2", "kappa3",
      "L1", "L2", "L3", "L123",
      "dubins_type"
    };

    arg_out_0 = mxCreateStructMatrix(1,1,20,fieldnames);

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateDoubleScalar(C0.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(C0.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(C0.theta_begin()) );

    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(C1.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mxCreateDoubleScalar(C1.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 5, mxCreateDoubleScalar(C1.theta_begin()) );

    mxSetFieldByNumber( arg_out_0, 0, 6, mxCreateDoubleScalar(C2.x_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 7, mxCreateDoubleScalar(C2.y_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 8, mxCreateDoubleScalar(C2.theta_begin()) );

    mxSetFieldByNumber( arg_out_0, 0, 9, mxCreateDoubleScalar(C2.x_end()) );
    mxSetFieldByNumber( arg_out_0, 0, 10, mxCreateDoubleScalar(C2.y_end()) );
    mxSetFieldByNumber( arg_out_0, 0, 11, mxCreateDoubleScalar(C2.theta_end()) );

    mxSetFieldByNumber( arg_out_0, 0, 12, mxCreateDoubleScalar(C0.kappa_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 13, mxCreateDoubleScalar(C1.kappa_begin()) );
    mxSetFieldByNumber( arg_out_0, 0, 14, mxCreateDoubleScalar(C2.kappa_begin()) );

    mxSetFieldByNumber( arg_out_0, 0, 15, mxCreateDoubleScalar(C0.length()) );
    mxSetFieldByNumber( arg_out_0, 0, 16, mxCreateDoubleScalar(C1.length()) );
    mxSetFieldByNumber( arg_out_0, 0, 17, mxCreateDoubleScalar(C2.length()) );

    mxSetFieldByNumber( arg_out_0, 0, 18, mxCreateDoubleScalar(C0.length()+C1.length()+C2.length()) );

    mxSetFieldByNumber( arg_out_0, 0, 19, mxCreateDoubleScalar(ptr->icode()) );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_length( int nlhs, mxArray       *plhs[],
                 int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('get_length',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2,               CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1 || nlhs == 3 , CMD "expected 1 or 3 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->length() );
    if ( nlhs == 3 ) {
      real_type grad[2];
      ptr->length_grad( grad );
      Utils::mex_set_scalar_value( arg_out_1, grad[0] );
      Utils::mex_set_scalar_value( arg_out_2, grad[1] );
    }

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_curve_type( int nlhs, mxArray       *plhs[],
                 int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('curve_type',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    Utils::mex_set_scalar_value( arg_out_0, static_cast<int>(ptr->solution_type()) );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_curve_type_string( int nlhs, mxArray       *plhs[],
                        int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('curve_type_string',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    string s1{ ptr->solution_type_string() };
    string s2{ ptr->solution_type_string_short() };

    plhs[0] = mxCreateString( s1.c_str() );
    plhs[1] = mxCreateString( s2.c_str() );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_range_angles_begin( int nlhs, mxArray       *plhs[],
                             int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('get_range_angles_begin',OBJ,x0,y0,x1,y1,theta1,k_max): "

    UTILS_MEX_ASSERT( nrhs == 8, CMD "expected 8 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading x1" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading y1" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading theta1" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading k_max" );

    real_type angles[12];
    integer npts = ptr->get_range_angles_begin( x0, y0, x1, y1, theta1, k_max, angles );

    double * Angles = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
    std::copy_n( angles, npts, Angles );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_range_angles_end( int nlhs, mxArray       *plhs[],
                           int nrhs, mxArray const *prhs[] ) {

    #define CMD "DubinsMexWrapper('get_range_angles_init',OBJ,x0,y0,theta0,x1,y1,k_max): "

    UTILS_MEX_ASSERT( nrhs == 8, CMD "expected 8 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins * ptr{ Utils::mex_convert_mx_to_ptr<Dubins>(arg_in_1) };

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading x1" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading y1" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading k_max" );

    real_type angles[12];
    integer npts = ptr->get_range_angles_end( x0, y0, theta0, x1, y1, k_max, angles );

    double * Angles = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
    std::copy_n( angles, npts, Angles );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"build",do_build},
    {"get_pars",do_get_pars},
    {"get_length",do_get_length},
    {"curve_type",do_curve_type},
    {"curve_type_string",do_curve_type_string},
    {"get_range_angles_begin",do_get_range_angles_begin},
    {"get_range_angles_end",do_get_range_angles_end},
    CMD_MAP_FUN
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
