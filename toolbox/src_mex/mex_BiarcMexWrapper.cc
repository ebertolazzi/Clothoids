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
"BiarcMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = BiarcMexWrapper( 'new' );\n" \
"\n" \
"    On output:\n" \
"    OBJ = pointer to the internal object\n" \
"\n" \
"  - Build:\n" \
"    BiarcMexWrapper( 'build_G1', OBJ, x0, y0, theta0, x1, y1, theta1 );\n" \
"    BiarcMexWrapper( 'build_3P', OBJ, x0, y0, x1, y1, x2, y2 );\n" \
"    [arc0,arc1] = BiarcMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
"  - Eval:\n" \
"    [x0,y0,theta0,kappa0,L0,x1,y1,theta1,kappa1,L1] = BiarcMexWrapper( 'get_pars', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("BiarcMexWrapper") \
MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define CMD_BASE "BiarcMexWrapper"
  #define G2LIB_CLASS Biarc
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

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new( int nlhs, mxArray       *plhs[],
          int nrhs, mxArray const *[] ) {

    #define CMD "BiarcMexWrapper('new'): "
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    #undef CMD

    arg_out_0 = Utils::mex_convert_ptr_to_mx<Biarc>(new Biarc());
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_G1( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define CMD "BiarcMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1): "

    UTILS_MEX_ASSERT( nrhs == 8, CMD "expected 8 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading x1" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading y1" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading theta1" );

    bool ok = ptr->build( x0, y0, theta0, x1, y1, theta1 );

    // returns the status of the interpolation
    Utils::mex_set_scalar_bool( arg_out_0, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_3P( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('build_3P',OBJ,...): "
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    #undef CMD

    if ( nrhs == 8 ) {
      #define CMD "BiarcMexWrapper('build_3P',OBJ,x0,y0,x1,y1,x2,y2): "

      real_type x0 = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
      real_type y0 = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
      real_type x1 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading x1" );
      real_type y1 = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading y1" );
      real_type x2 = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading x2" );
      real_type y2 = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading y2" );

      bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 );

      // returns the status of the interpolation
      Utils::mex_set_scalar_bool(arg_out_0,ok);

      #undef CMD
    } else if ( nrhs == 5 ) {
      #define CMD "BiarcMexWrapper('build_3P',OBJ,p0,p1,p2): "
      real_type const * p0;
      real_type const * p1;
      real_type const * p2;

      mwSize n;
      p0 = Utils::mex_vector_pointer( arg_in_2, n, CMD "Error in reading p0" );
      UTILS_MEX_ASSERT(
        n == 2,
        CMD "Error in reading length(p0) == {} expect length(p0) == 2\n", n
      );
      p1 = Utils::mex_vector_pointer( arg_in_3, n, CMD "Error in reading p1" );
      UTILS_MEX_ASSERT(
        n == 2,
        CMD "Error in reading length(p1) == {} expect length(p1) == 2\n", n
      );
      p2 = Utils::mex_vector_pointer( arg_in_4, n, CMD "Error in reading p2" );
      UTILS_MEX_ASSERT(
        n == 2,
        CMD "Error in reading length(p2) == {} expect length(p2) == 2\n", n
      );

      bool ok = ptr->build_3P( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );

      // returns the status of the interpolation
      Utils::mex_set_scalar_bool(arg_out_0,ok);

      #undef CMD
    } else {
      UTILS_MEX_ASSERT0(
        false, "BiarcMexWrapper('build_3P',OBJ,...) expected 5 or 8 arguments\n"
      );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('x_middle',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->x_middle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('y_middle',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->y_middle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_theta_middle( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('theta_middle',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->theta_middle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_s_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('s_middle',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->length0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_kappa0( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('kappa0',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->kappa0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_kappa1( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('kappa1',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->kappa1());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_length0( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('length0',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->length0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_length1( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('length1',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    Utils::mex_set_scalar_value( arg_out_0, ptr->length1());;
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_to_nurbs( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_1);

    #define CMD "BiarcMexWrapper('to_nurbs',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

    CircleArc const & C0 = ptr->C0();
    CircleArc const & C1 = ptr->C1();

    integer npts0, nknots0, npts1, nknots1;
    C0.paramNURBS( nknots0, npts0 );
    C1.paramNURBS( nknots1, npts1 );

    mxArray * mx_knots0, * mx_Poly0, * mx_knots1, * mx_Poly1;

    double * knots0 = Utils::mex_create_matrix_value( mx_knots0, 1, nknots0 );
    double * poly0  = Utils::mex_create_matrix_value( mx_Poly0,  3, npts0   );
    double * knots1 = Utils::mex_create_matrix_value( mx_knots1, 1, nknots1 );
    double * poly1  = Utils::mex_create_matrix_value( mx_Poly1,  3, npts1   );

    C0.toNURBS( knots0, reinterpret_cast<real_type (*)[3]>(poly0) );
    C1.toNURBS( knots1, reinterpret_cast<real_type (*)[3]>(poly1) );

    static char const * fieldnames[] = {
      "form",
      "order",
      "dim",
      "number",
      "knots",
      "coefs"
    };

    arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
    arg_out_1 = mxCreateStructMatrix(1,1,6,fieldnames);

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts0) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots0 );
    mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly0 );

    mxSetFieldByNumber( arg_out_1, 0, 0, mxCreateString("rB") );
    mxSetFieldByNumber( arg_out_1, 0, 1, mxCreateDoubleScalar(3) );
    mxSetFieldByNumber( arg_out_1, 0, 2, mxCreateDoubleScalar(2) );
    mxSetFieldByNumber( arg_out_1, 0, 3, mxCreateDoubleScalar(npts1) );
    mxSetFieldByNumber( arg_out_1, 0, 4, mx_knots1 );
    mxSetFieldByNumber( arg_out_1, 0, 5, mx_Poly1 );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    //{"build",do_build},
    {"build_3P",do_build_3P},
    {"build_G1",do_build_G1},
    {"x_middle",do_x_middle},
    {"y_middle",do_y_middle},
    {"theta_middle",do_theta_middle},
    {"s_middle",do_s_middle},
    {"kappa0",do_kappa0},
    {"kappa1",do_kappa1},
    {"length0",do_length0},
    {"length1",do_length1},
    {"to_nurbs",do_to_nurbs},
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
      mexErrMsgTxt( fmt::format( "Biarc Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Biarc failed\n");
    }

  }
}
