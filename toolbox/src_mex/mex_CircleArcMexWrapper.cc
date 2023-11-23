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

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute cicle arc\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = CircleArcMexWrapper( 'new' );\n" \
"  OBJ = CircleArcMexWrapper( 'new', x0, y0, theta0, kur, L );\n" \
"\n" \
"  CircleArcMexWrapper( 'build', OBJ, x0, y0, theta0, kur, L );\n" \
"  CircleArcMexWrapper( 'build_G1', OBJ, p0, theta0, p1 );\n" \
"  CircleArcMexWrapper( 'build_G1', OBJ, x0, y0, theta0, x1, y1 );\n" \
"  CircleArcMexWrapper( 'build_3P', OBJ, p0, p1, p2 );\n" \
"  CircleArcMexWrapper( 'build_3P', OBJ, x0, y0, x1, y1, x2, y2 );\n" \
"\n" \
"  [p0,p1,p2] = CircleArcMexWrapper( 'bbTriangles', OBJ );\n" \
"\n" \
"  burbs = CircleArcMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("CircleArcMexWrapper") \
MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define CMD_BASE "CircleArcMexWrapper"
  #define G2LIB_CLASS CircleArc
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

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('new',...): "
    UTILS_MEX_ASSERT(
      nrhs == 1 || nrhs == 6,
      CMD "expected 1 or 6 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 1,
      CMD "expected 1 output, nlhs = {}\n", nlhs
    );
    #undef CMD

    CircleArc * ptr = new CircleArc();
    arg_out_0 = Utils::mex_convert_ptr_to_mx<CircleArc>(ptr);

    if ( nrhs == 6 ) {
      #define CMD "CircleArcMexWrapper('new',x0,y0,theta0,k0,L): "
      real_type x0, y0, theta0, k0, L;
      x0     = Utils::mex_get_scalar_value( arg_in_1, CMD "`x0` expected to be a real scalar" );
      y0     = Utils::mex_get_scalar_value( arg_in_2, CMD "`y0` expected to be a real scalar" );
      theta0 = Utils::mex_get_scalar_value( arg_in_3, CMD "`theta0` expected to be a real scalar" );
      k0     = Utils::mex_get_scalar_value( arg_in_4, CMD "`k0` expected to be a real scalar" );
      L      = Utils::mex_get_scalar_value( arg_in_5, CMD "`L` expected to be a real scalar" );

      ptr->build( x0, y0, theta0, k0, L );
      #undef CMD
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_arc(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('build_arc',OBJ,x0,y0,theta0,k0,L): "
    UTILS_MEX_ASSERT(
      nrhs == 7, CMD "expected 7 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs
    );

    CircleArc * ptr = Utils::mex_convert_mx_to_ptr<CircleArc>( arg_in_1 );

    real_type x0, y0, theta0, k0, L;
    x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "`x0` expected to be a real scalar" );
    y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "`y0` expected to be a real scalar" );
    theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "`theta0` expected to be a real scalar" );
    k0     = Utils::mex_get_scalar_value( arg_in_5, CMD "`k0` expected to be a real scalar" );
    L      = Utils::mex_get_scalar_value( arg_in_6, CMD "`L` expected to be a real scalar" );

    ptr->build( x0, y0, theta0, k0, L );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_3P(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('build_3P',...): "
    UTILS_MEX_ASSERT(
      nlhs == 0 || nlhs ==1,
      CMD " expected 1 or no output, nlhs = {}\n", nlhs
    );
    #undef CMD

    CircleArc * ptr = Utils::mex_convert_mx_to_ptr<CircleArc>( arg_in_1 );

    real_type x0(0), y0(0), x1(0), y1(0), x2(0), y2(0);
    if ( nrhs == 5 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,p0,p1,p2): "
      real_type const * p0;
      real_type const * p1;
      real_type const * p2;
      mwSize size0, size1, size2;

      p0 = Utils::mex_vector_pointer(
        arg_in_2, size0, CMD "`p0` expected to be a real vector"
      );
      p1 = Utils::mex_vector_pointer(
        arg_in_3, size1, CMD "`p1` expected to be a real vector"
      );
      p2 = Utils::mex_vector_pointer(
        arg_in_4, size2, CMD "`p2` expected to be a real vector"
      );

      UTILS_MEX_ASSERT(
        size0 == 2 && size1 == 2 && size2 == 2,
        CMD "bad dimension size(p0) = {}, size(p1) = {}, size(p2) = {}\n",
        size0, size1, size2
      );
      #undef CMD

      x0 = p0[0]; y0 = p0[1];
      x1 = p1[0]; y1 = p1[1];
      x2 = p2[0]; y2 = p2[1];

    } else if ( nrhs == 8 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,x0,x1,x1,y1,x2,y2): "
      x0 = Utils::mex_get_scalar_value( arg_in_2, CMD "`x0` expected to be a scalar value" );
      y0 = Utils::mex_get_scalar_value( arg_in_3, CMD "`y0` expected to be a scalar value" );
      x1 = Utils::mex_get_scalar_value( arg_in_4, CMD "`x1` expected to be a scalar value" );
      y1 = Utils::mex_get_scalar_value( arg_in_5, CMD "`y1` expected to be a scalar value" );
      x2 = Utils::mex_get_scalar_value( arg_in_6, CMD "`x2` expected to be a scalar value" );
      y2 = Utils::mex_get_scalar_value( arg_in_7, CMD "`y2` expected to be a scalar value" );
      #undef CMD
    } else {
      UTILS_MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = {}\n", nrhs );
    }

    bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 );
    if ( nlhs == 1 ) Utils::mex_set_scalar_bool( arg_out_0, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_G1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('build_G1',...): "
    UTILS_MEX_ASSERT(
      nlhs == 0 || nlhs ==1,
      CMD "expected 1 or no output, nlhs = {}\n", nlhs
    );
    #undef CMD

    CircleArc * ptr = Utils::mex_convert_mx_to_ptr<CircleArc>( arg_in_1 );

    real_type x0(0), y0(0), x1(0), y1(0), theta0(0);
    if ( nrhs == 5 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,p0,theta0,p1): "
      real_type const * p0;
      real_type const * p1;
      mwSize size0, size1;
      p0 = Utils::mex_vector_pointer(
        arg_in_2, size0, CMD "`p0` expected to be a real vector"
      );
      p1 = Utils::mex_vector_pointer(
        arg_in_4, size1, CMD "`p1` expected to be a real vector"
      );
      theta0 = Utils::mex_get_scalar_value(
        arg_in_3, CMD "`theta0` expected to be a real vector"
      );

      UTILS_MEX_ASSERT(
        size0 == 2 && size1 == 2,
        CMD "bad dimension size(p0) = {}, size(p1) = {}\n", size0, size1
      );
      #undef CMD

      x0 = p0[0]; y0 = p0[1];
      x1 = p1[0]; y1 = p1[1];

    } else if ( nrhs == 7 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,x0,x1,theta0,x1,y1): "
      x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "`x0` expected to be a scalar value" );
      y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "`y0` expected to be a scalar value" );
      theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "`theta0` expected to be a scalar value" );
      x1     = Utils::mex_get_scalar_value( arg_in_5, CMD "`x1` expected to be a scalar value" );
      y1     = Utils::mex_get_scalar_value( arg_in_6, CMD "`y1` expected to be a scalar value" );
      #undef CMD
    } else {
      UTILS_MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = {}\n", nrhs );
    }

    bool ok = ptr->build_G1( x0, y0, theta0, x1, y1 );
    if ( nlhs == 1 ) Utils::mex_set_scalar_bool( arg_out_0, ok );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_change_curvilinear_origin(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('change_curvilinear_origin',OBJ,s0,L): "
    UTILS_MEX_ASSERT(nlhs == 0, CMD "expected NOoutputs, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );

    CircleArc * ptr = Utils::mex_convert_mx_to_ptr<CircleArc>( arg_in_1 );

    real_type s0 = Utils::mex_get_scalar_value(arg_in_2,CMD "Error in reading s0");
    real_type L  = Utils::mex_get_scalar_value(arg_in_3,CMD "Error in reading L");
    ptr->change_curvilinear_origin(s0,L);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_to_nurbs(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    CircleArc * ptr = Utils::mex_convert_mx_to_ptr<CircleArc>( arg_in_1 );

    #define CMD "CircleArcMexWrapper('to_nurbs',OBJ): "
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    integer npts, nknots;
    ptr->paramNURBS( nknots, npts );

    mxArray * mx_knots, * mx_Poly;
    double * knots = Utils::mex_create_matrix_value( mx_knots, 1, nknots );
    double * poly  = Utils::mex_create_matrix_value( mx_Poly,  3, npts );

    ptr->toNURBS( knots, reinterpret_cast<real_type (*)[3]>(poly) );

    static char const * fieldnames[] = { "form", "order", "dim", "number", "knots", "coefs" };
    arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(3) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots );
    mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"build_arc",do_build_arc},
    {"build_3P",do_build_3P},
    {"build_G1",do_build_G1},
    {"change_curvilinear_origin",do_change_curvilinear_origin},
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
      mexErrMsgTxt( fmt::format( "CircleArc Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("CircleArc failed\n");
    }

  }
}
