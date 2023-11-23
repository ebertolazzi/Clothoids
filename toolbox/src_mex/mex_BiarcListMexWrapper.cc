/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2019
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

#include <fstream>

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"BiarcListMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = BiarcListMexWrapper( 'new' );\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
"\n" \
"  - Build:\n" \
"    BiarcListMexWrapper( 'push_back', OBJ, CLOT );\n" \
"    BiarcListMexWrapper( 'push_back', OBJ, kappa0, dkappa, L );\n" \
"    BiarcListMexWrapper( 'push_back', OBJ, x0, y0, theta0, kappa0, dkappa, L );\n" \
"    BiarcListMexWrapper( 'push_back_G1', OBJ, x1, y1, theta1 );\n" \
"    BiarcListMexWrapper( 'push_back_G1', OBJ, x0, y0, theta0, x1, y1, theta1 );\n" \
"    BiarcListMexWrapper( 'copy', OBJ, OBJ1 );\n" \
"\n" \
"    [s,theta,kappa] = BiarcListMexWrapper( 'get_STK', OBJ );\n" \
"    [x,y]           = BiarcListMexWrapper( 'get_XY', OBJ );\n" \
"\n" \
"  - Bounding Box:\n" \
"    TT = BiarcListMexWrapper( 'bbox', OBJ, max_angle, max_size );%\n" \
"    TT = BiarcListMexWrapper( 'bbox', OBJ, max_angle, max_size, offs );%\n" \
"\n" \
"    ok = BiarcListMexWrapper( 'build_G1', OBJ, x, y [,theta] );%\n" \
"    [theta,ok] = BiarcListMexWrapper( 'build_theta', OBJ, x, y );%\n" \
"\n" \
MEX_INFO_MESSAGE("BiarcListMexWrapper") \
MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define CMD_BASE "BiarcListMexWrapper"
  #define G2LIB_CLASS BiarcList
  #include "mex_common.hxx"
  #undef CMD_BASE
  #undef G2LIB_CLASS

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

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *[]
  ) {

    #define CMD "BiarcListMexWrapper('new'): "
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = {}\n", nrhs );

    arg_out_0 = Utils::mex_convert_ptr_to_mx<BiarcList>(new BiarcList());

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back_G1(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('push_back_G1',OBJ,[x0,y0,theta0,x1,y1,theta1]|[x1,y1,theta1]): "

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    if ( nrhs == 8 ) {
      real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
      real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
      real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
      real_type x1     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading x1" );
      real_type y1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading y1" );
      real_type theta1 = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading theta1" );
      ptr->push_back_G1( x0, y0, theta0, x1, y1, theta1 );
    } else if ( nrhs == 5 ) {
      real_type x1     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x1" );
      real_type y1     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y1" );
      real_type theta1 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta1" );
      ptr->push_back_G1( x1, y1, theta1 );
    } else {
      UTILS_MEX_ASSERT( false, CMD "expected 5 or 8 inputs, nrhs = {|\n", nrhs );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_reserve(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('reserve',OBJ,N): "

    UTILS_MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0 , CMD "expected no outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    int64_t N = Utils::mex_get_int64( arg_in_2, CMD "Error in reading N" );
    ptr->reserve( N );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_STK(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('get_STK',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    integer n = ptr->num_segments();

    real_type * s     = Utils::mex_create_matrix_value( arg_out_0, 1, n+1 );
    real_type * theta = Utils::mex_create_matrix_value( arg_out_1, 1, n+1 );
    real_type * kappa = Utils::mex_create_matrix_value( arg_out_2, 1, n+1 );

    ptr->get_STK( s, theta, kappa );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get_XY(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('get_XY',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    integer     n = ptr->num_segments();
    real_type * x = Utils::mex_create_matrix_value( arg_out_0, 1, n+1 );
    real_type * y = Utils::mex_create_matrix_value( arg_out_1, 1, n+1 );

    ptr->get_XY( x, y );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_G1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('build_G1', OBJ, x, y [, theta]): "

    UTILS_MEX_ASSERT(nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs);

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    bool ok = true;

    if ( nrhs == 4 ) {
      mwSize nx, ny;
      real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y = Utils::mex_vector_pointer( arg_in_3, ny, CMD "Error in reading y" );

      UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );

      ok = ptr->build_G1( nx, x, y );

    } else if ( nrhs == 5 ) {

      mwSize nx, ny, nt;
      real_type const * x     = Utils::mex_vector_pointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y     = Utils::mex_vector_pointer( arg_in_3, ny, CMD "Error in reading y" );
      real_type const * theta = Utils::mex_vector_pointer( arg_in_4, nt, CMD "Error in reading theta" );

      UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );
      UTILS_MEX_ASSERT( nx == nt, CMD "length(theta) = {} != length(x) = length(y) = {}\n", nt, ny );

      ok = ptr->build_G1( nx, x, y, theta );

    } else {
      UTILS_MEX_ASSERT( false, CMD "expected 4 or 5 input, nrhs = {}\n", nrhs );
    }

    Utils::mex_set_scalar_bool( arg_out_0, ok );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_theta(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('build_theta',OBJ,x,y): "

    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    mwSize nx, ny;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD "Error in reading x" );
    real_type const * y = Utils::mex_vector_pointer( arg_in_3, ny, CMD "Error in reading y" );

    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );

    real_type * theta = Utils::mex_create_matrix_value( arg_out_0, nx, 1 );

    bool ok = build_guess_theta( nx, x, y, theta );

    Utils::mex_set_scalar_bool( arg_out_1, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('get',OBJ,n): "

    UTILS_MEX_ASSERT(nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT(nlhs == 6, CMD "expected 6 output, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    int64_t n = Utils::mex_get_int64( arg_in_2, CMD "Error in reading n" );

    UTILS_MEX_ASSERT(
      n > 0 && n <= ptr->num_segments(),
      CMD "n = {} must be >= 1 and <= {}\n", n, ptr->num_segments()
    );

    Biarc const & c = ptr->get(n-1);

    Utils::mex_set_scalar_value(arg_out_0, c.x_begin());
    Utils::mex_set_scalar_value(arg_out_1, c.y_begin());
    Utils::mex_set_scalar_value(arg_out_2, c.theta_begin());
    Utils::mex_set_scalar_value(arg_out_3, c.x_end());
    Utils::mex_set_scalar_value(arg_out_4, c.y_end());
    Utils::mex_set_scalar_value(arg_out_5, c.theta_end());

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_num_segments(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "BiarcListMexWrapper('num_segments', OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    Utils::mex_set_scalar_int32( arg_out_0, ptr->num_segments() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_findST1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('findST1',OBJ,x,y): "
    UTILS_MEX_ASSERT(
      nrhs == 4,
      CMD "expected 4 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 3,
      CMD "expected 3 output, nlhs = {}\n", nlhs
    );

    BiarcList * ptr = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_1);

    mwSize nrx, ncx, nry, ncy;
    real_type const * x;
    real_type const * y;

    x = Utils::mex_matrix_pointer(
      arg_in_2, nrx, ncx,
      CMD "`x` expected to be a real vector/matrix"
    );

    y = Utils::mex_matrix_pointer(
      arg_in_3, nry, ncy,
      CMD "`y` expected to be a real vector/matrix"
    );

    UTILS_MEX_ASSERT(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found\n"
      "size(x) = {} x {} size(y) = {} x {}\n",
      nrx, ncx, nry, ncy
    );

    real_type * s   = Utils::mex_create_matrix_value( arg_out_0, nrx, ncx );
    real_type * t   = Utils::mex_create_matrix_value( arg_out_1, nrx, ncx );
    real_type * idx = Utils::mex_create_matrix_value( arg_out_2, nrx, ncx );

    mwSize size = nrx*ncx;
    for ( mwSize i = 0; i < size; ++i ) {
      integer nseg = ptr->findST1( *x++, *y++, *s++, *t++ );
      *idx++ = nseg >= 0 ? nseg+1 : nseg;
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"push_back_G1",do_push_back_G1},
    {"reserve",do_reserve},
    {"get_STK",do_get_STK},
    {"get_XY",do_get_XY},
    {"build_G1",do_build_G1},
    {"build_theta",do_build_theta},
    {"get",do_get},
    {"num_segments",do_num_segments},
    {"findST1",do_findST1},
    CMD_MAP_FUN
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
      UTILS_MEX_ASSERT0( mxIsChar(arg_in_0), "First argument must be a string" );
      mxGetString( arg_in_0, cmd, 256 );
      cmd_to_fun.at(cmd)( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "BiarcList Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt( "BiarcList failed" );
    }
  }
}
