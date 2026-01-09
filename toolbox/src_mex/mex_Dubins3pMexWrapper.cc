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
#pragma comment( lib, "IPHLPAPI.lib" )
#pragma comment( lib, "ws2_32.lib" )
#pragma comment( lib, "Shlwapi.lib" )
#pragma comment( lib, "Advapi32.lib" )
#pragma comment( lib, "Shell32.lib" )
#pragma comment( lib, "kernel32.lib" )
#endif

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"
#include "Utils_string.hh"
#include "Utils_mex.hh"
#include "mex_info.hxx"

#define MEX_ERROR_MESSAGE                                                                                       \
  "=====================================================================================\n"                     \
  "Dubins3pMexWrapper:  Compute solution of Dubins3p problem\n"                                                 \
  "\n"                                                                                                          \
  "USAGE:\n"                                                                                                    \
  "  - Constructors:\n"                                                                                         \
  "    OBJ = Dubins3pMexWrapper( 'new' );\n"                                                                    \
  "\n"                                                                                                          \
  "    On output:\n"                                                                                            \
  "    OBJ = pointer to the internal object\n"                                                                  \
  "\n"                                                                                                          \
  "  - Build:\n"                                                                                                \
  "    [arc0,arc1,arc2] = Dubins3pMexWrapper( 'build', OBJ, x0, y0, theta0, xm, ym, x1, y1, theta1, k_max );\n" \
  "\n" MEX_INFO_MESSAGE( "Dubins3pMexWrapper" ) MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib
{
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

#define CMD_BASE "Dubins3pMexWrapper"
#define G2LIB_CLASS Dubins3p
#include "mex_common.hxx"
// #undef CMD_BASE
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

  static void do_new( int nlhs, mxArray * plhs[], int nrhs, mxArray const *[] )
  {
#define CMD CMD_BASE "('new'): "
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
#undef CMD

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Dubins3p>( new Dubins3p( "Dubins3p" ) );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_build( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('build',OBJ,x0,y0,theta0,xm,ym,x1,y1,theta1,k_max,method): "

    UTILS_MEX_ASSERT( nrhs == 12, CMD "expected 12 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
    real_type xm     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading xm" );
    real_type ym     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading ym" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading x1" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading y1" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading theta1" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_10, CMD "Error in reading k_max" );

    char method_str[256];
    UTILS_MEX_ASSERT0( mxIsChar( arg_in_11 ), CMD "last argument must be a string" );
    mxGetString( arg_in_11, method_str, 256 );

    Dubins3pBuildType method{ string_to_Dubins3pBuildType( method_str ) };
    bool              ok = ptr->build( x0, y0, theta0, xm, ym, x1, y1, theta1, k_max, method );

    // returns the status of the interpolation
    Utils::mex_set_scalar_bool( arg_out_0, ok );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_get_pars( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('get_pars',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    CircleArc const & C0{ ptr->C0() };
    CircleArc const & C1{ ptr->C1() };
    CircleArc const & C2{ ptr->C2() };
    CircleArc const & C3{ ptr->C3() };
    CircleArc const & C4{ ptr->C4() };
    CircleArc const & C5{ ptr->C5() };

    static char const * fieldnames[] = { "x0",          "y0",
                                         "theta0",      "x1",
                                         "y1",          "theta1",
                                         "x2",          "y2",
                                         "theta2",      "x3",
                                         "y3",          "theta3",
                                         "x4",          "y4",
                                         "theta4",      "x5",
                                         "y5",          "theta5",
                                         "x6",          "y6",
                                         "theta6",      "kappa1",
                                         "kappa2",      "kappa3",
                                         "kappa4",      "kappa5",
                                         "kappa6",      "L1",
                                         "L2",          "L3",
                                         "L4",          "L5",
                                         "L6",          "L123",
                                         "L456",        "L123456",
                                         "tolerance",   "max_evaluation",
                                         "evaluation",  "dubins1_type",
                                         "dubins2_type" };

    arg_out_0 = mxCreateStructMatrix( 1, 1, 41, fieldnames );

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateDoubleScalar( C0.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar( C0.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar( C0.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar( C1.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mxCreateDoubleScalar( C1.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 5, mxCreateDoubleScalar( C1.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 6, mxCreateDoubleScalar( C2.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 7, mxCreateDoubleScalar( C2.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 8, mxCreateDoubleScalar( C2.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 9, mxCreateDoubleScalar( C3.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 10, mxCreateDoubleScalar( C3.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 11, mxCreateDoubleScalar( C3.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 12, mxCreateDoubleScalar( C4.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 13, mxCreateDoubleScalar( C4.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 14, mxCreateDoubleScalar( C4.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 15, mxCreateDoubleScalar( C5.x_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 16, mxCreateDoubleScalar( C5.y_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 17, mxCreateDoubleScalar( C5.theta_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 18, mxCreateDoubleScalar( C5.x_end() ) );
    mxSetFieldByNumber( arg_out_0, 0, 19, mxCreateDoubleScalar( C5.y_end() ) );
    mxSetFieldByNumber( arg_out_0, 0, 20, mxCreateDoubleScalar( C5.theta_end() ) );

    mxSetFieldByNumber( arg_out_0, 0, 21, mxCreateDoubleScalar( C0.kappa_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 22, mxCreateDoubleScalar( C1.kappa_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 23, mxCreateDoubleScalar( C2.kappa_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 24, mxCreateDoubleScalar( C3.kappa_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 25, mxCreateDoubleScalar( C4.kappa_begin() ) );
    mxSetFieldByNumber( arg_out_0, 0, 26, mxCreateDoubleScalar( C5.kappa_begin() ) );

    mxSetFieldByNumber( arg_out_0, 0, 27, mxCreateDoubleScalar( C0.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 28, mxCreateDoubleScalar( C1.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 29, mxCreateDoubleScalar( C2.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 30, mxCreateDoubleScalar( C3.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 31, mxCreateDoubleScalar( C4.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 32, mxCreateDoubleScalar( C5.length() ) );

    mxSetFieldByNumber( arg_out_0, 0, 33, mxCreateDoubleScalar( C0.length() + C1.length() + C2.length() ) );
    mxSetFieldByNumber( arg_out_0, 0, 34, mxCreateDoubleScalar( C3.length() + C4.length() + C5.length() ) );
    mxSetFieldByNumber(
      arg_out_0,
      0,
      35,
      mxCreateDoubleScalar( C0.length() + C1.length() + C2.length() + C3.length() + C4.length() + C5.length() ) );

    mxSetFieldByNumber( arg_out_0, 0, 36, mxCreateDoubleScalar( ptr->tolerance() ) );
    mxSetFieldByNumber( arg_out_0, 0, 37, mxCreateDoubleScalar( ptr->max_num_evaluation() ) );
    mxSetFieldByNumber( arg_out_0, 0, 38, mxCreateDoubleScalar( ptr->num_evaluation() ) );

    mxSetFieldByNumber( arg_out_0, 0, 39, mxCreateDoubleScalar( ptr->icode0() ) );
    mxSetFieldByNumber( arg_out_0, 0, 40, mxCreateDoubleScalar( ptr->icode1() ) );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_curve_type( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('curve_type',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, static_cast<int>( ptr->solution_type0() ) );
    Utils::mex_set_scalar_value( arg_out_1, static_cast<int>( ptr->solution_type1() ) );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_curve_type_string( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('curve_type_string',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    string s1{ ptr->solution_type_string() };
    string s2{ ptr->solution_type_string_short() };

    plhs[0] = mxCreateString( s1.c_str() );
    plhs[1] = mxCreateString( s2.c_str() );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_num_evaluation( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('num_evaluation',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    Utils::mex_set_scalar_value( arg_out_0, ptr->num_evaluation() );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_set_tolerance( int nlhs, mxArray *[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('set_tolerance',OBJ,tol): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    real_type tol = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading tol" );
    ptr->set_tolerance( tol );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_set_sample_angle( int nlhs, mxArray *[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('set_sample_angle',OBJ,ang): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    real_type ang = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading ang" );
    ptr->set_sample_angle( ang );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_set_sample_points( int nlhs, mxArray *[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('set_sample_points',OBJ,ang): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    Utils::int64_t npts = Utils::mex_get_int64( arg_in_2, CMD "Error in reading npts" );
    ptr->set_sample_points( npts );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_set_max_evaluation( int nlhs, mxArray *[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "Dubins3pMexWrapper('set_max_evaluation',OBJ,tol): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    Utils::int64_t max_eval = Utils::mex_get_int64( arg_in_2, CMD "Error in reading tol" );
    ptr->set_max_evaluation( max_eval );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_get_range_angles( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "DubinsMexWrapper('get_range_angles',OBJ,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max): "

    UTILS_MEX_ASSERT( nrhs == 11, CMD "expected 8 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    real_type xi     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading xi" );
    real_type yi     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading yi" );
    real_type thetai = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading thetai" );
    real_type xm     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading xm" );
    real_type ym     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading ym" );
    real_type xf     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading xf" );
    real_type yf     = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading yf" );
    real_type thetaf = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading thetaf" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_10, CMD "Error in reading k_max" );

    real_type angles[12];
    integer   npts = ptr->get_range_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, angles );

    double * Angles = Utils::mex_create_matrix_value( arg_out_0, 1, npts );
    std::copy_n( angles, npts, Angles );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static void do_get_sample_angles( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
#define CMD "DubinsMexWrapper('get_sample_angles',OBJ,xi,yi,thetai,xm,ym,xf,yf,thetaf,k_max,tolerance): "

    UTILS_MEX_ASSERT( nrhs == 12, CMD "expected 12 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Dubins3p * ptr{ Utils::mex_convert_mx_to_ptr<Dubins3p>( arg_in_1 ) };

    real_type xi     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading xi" );
    real_type yi     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading yi" );
    real_type thetai = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading thetai" );
    real_type xm     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading xm" );
    real_type ym     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading ym" );
    real_type xf     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading xf" );
    real_type yf     = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading yf" );
    real_type thetaf = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading thetaf" );
    real_type k_max  = Utils::mex_get_scalar_value( arg_in_10, CMD "Error in reading k_max" );
    real_type tol    = Utils::mex_get_scalar_value( arg_in_11, CMD "Error in reading tolerance" );

    vector<real_type> angles;
    ptr->get_sample_angles( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, tol, angles );

    double * Angles = Utils::mex_create_matrix_value( arg_out_0, 1, angles.size() );
    std::copy( angles.begin(), angles.end(), Angles );

#undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void ( *DO_CMD )( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] );

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string, DO_CMD> cmd_to_fun{ { "new", do_new },
                                                   { "build", do_build },
                                                   { "get_pars", do_get_pars },
                                                   { "curve_type", do_curve_type },
                                                   { "curve_type_string", do_curve_type_string },
                                                   { "num_evaluation", do_num_evaluation },
                                                   { "set_tolerance", do_set_tolerance },
                                                   { "set_sample_angle", do_set_sample_angle },
                                                   { "set_sample_points", do_set_sample_points },
                                                   { "set_max_evaluation", do_set_max_evaluation },
                                                   { "get_range_angles", do_get_range_angles },
                                                   { "get_sample_angles", do_get_sample_angles },
                                                   CMD_MAP_FUN };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  extern "C" void mexFunction( int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[] )
  {
    char cmd[256];

    // the first argument must be a string
    if ( nrhs == 0 )
    {
      mexErrMsgTxt( MEX_ERROR_MESSAGE );
      return;
    }

    try
    {
      UTILS_MEX_ASSERT0( mxIsChar( arg_in_0 ), "First argument must be a string" );
      mxGetString( arg_in_0, cmd, 256 );
      cmd_to_fun.at( cmd )( nlhs, plhs, nrhs, prhs );
    }
    catch ( std::exception const & e )
    {
      mexErrMsgTxt( fmt::format( "Dubins3p Error: {}", e.what() ).c_str() );
    }
    catch ( ... )
    {
      mexErrMsgTxt( "Dubins3p failed\n" );
    }
  }
}  // namespace G2lib
