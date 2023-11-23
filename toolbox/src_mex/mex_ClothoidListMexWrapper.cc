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

#include <fstream>

#define MEX_ERROR_MESSAGE \
"=====================================================================================\n" \
"ClothoidListMexWrapper:  Compute parameters of the G1 Hermite clothoid fitting\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    OBJ = ClothoidListMexWrapper( 'new' );\n" \
"\n" \
"  On output:\n" \
"    OBJ = pointer to the internal object\n" \
"\n" \
"  - Build:\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, CLOT );\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, kappa0, dkappa, L );\n" \
"    ClothoidListMexWrapper( 'push_back', OBJ, x0, y0, theta0, kappa0, dkappa, L );\n" \
"    ClothoidListMexWrapper( 'push_back_G1', OBJ, x1, y1, theta1 );\n" \
"    ClothoidListMexWrapper( 'push_back_G1', OBJ, x0, y0, theta0, x1, y1, theta1 );\n" \
"    ClothoidListMexWrapper( 'copy', OBJ, OBJ1 );\n" \
"\n" \
"    [s,theta,kappa] = ClothoidListMexWrapper( 'get_STK', OBJ );\n" \
"    [x,y]           = ClothoidListMexWrapper( 'get_XY', OBJ );\n" \
"\n" \
"  - Bounding Box:\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size );%\n" \
"    TT = ClothoidListMexWrapper( 'bbox', OBJ, max_angle, max_size, offs );%\n" \
"\n" \
"  - G2 spline:\n" \
"    ok = ClothoidListMexWrapper( 'build_3arcG2' || 'build_2arcG2' || 'build_3arcCLC, ...\n" \
"                                 OBJ, ...\n" \
"                                 x0, y0, theta0, kappa0, ...\n" \
"                                 x1, y1, theta1, kappa1 );%\n" \
"    ok = ClothoidListMexWrapper( 'build_3arcG2fixed', OBJ, ...\n" \
"                                 s0, x0, y0, theta0, kappa0, ...\n" \
"                                 s1, x1, y1, theta1, kappa1 );%\n" \
"    ok = ClothoidListMexWrapper( 'build_G1', OBJ, x, y [,theta] );%\n" \
"    ok = ClothoidListMexWrapper( 'build_raw', OBJ, x, y, abscissa, theta, kappa );%\n" \
"    [theta,ok] = ClothoidListMexWrapper( 'build_theta', OBJ, x, y );%\n" \
"    dtheta = ClothoidListMexWrapper( 'deltaTheta', OBJ );%\n" \
"    dkappa = ClothoidListMexWrapper( 'deltaKappa', OBJ );%\n" \
"\n" \
MEX_INFO_MESSAGE("ClothoidListMexWrapper") \
MEX_INFO_MESSAGE_END

#include <unordered_map>

namespace G2lib {

  #define CMD_BASE "ClothoidListMexWrapper"
  #define G2LIB_CLASS ClothoidList
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

    #define CMD "ClothoidListMexWrapper('new'): "
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = {}\n", nrhs );

    //ClothoidList * ptr =
    arg_out_0 = Utils::mex_convert_ptr_to_mx<ClothoidList>(new ClothoidList());

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('push_back',OBJ,[type,OBJIN]|[kappa0,dkappa,L]|[x0,y0,theta0,kappa0,dkappa,L]): "

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    if ( nrhs == 8 ) {
      real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading x0" );
      real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading y0" );
      real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading theta0" );
      real_type kappa0 = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading kappa0" );
      real_type dkappa = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading dkappa" );
      real_type L      = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading L" );
      ptr->push_back( x0, y0, theta0, kappa0, dkappa, L );
    } else if ( nrhs == 5 ) {
      real_type kappa0 = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading kappa0" );
      real_type dkappa = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading dkappa" );
      real_type L      = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading L" );
      ptr->push_back( kappa0, dkappa, L );
    } else if ( nrhs == 4 ) {
      UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), "Third argument must be a string" );
      string type = mxArrayToString(arg_in_2);
      if ( type == "LineSegment" ) {
        LineSegment * cc = Utils::mex_convert_mx_to_ptr<LineSegment>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "BiArc" ) {
        Biarc * cc = Utils::mex_convert_mx_to_ptr<Biarc>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "BiarcList" ) {
        BiarcList * cc = Utils::mex_convert_mx_to_ptr<BiarcList>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "CircleArc" ) {
        CircleArc * cc = Utils::mex_convert_mx_to_ptr<CircleArc>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "ClothoidCurve" ) {
        ClothoidCurve const * cc = Utils::mex_convert_mx_to_ptr<ClothoidCurve>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "ClothoidList" ) {
        ClothoidList const * cc = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_3);
        ptr->push_back( *cc );
      } else if ( type == "PolyLine" ) {
        PolyLine const * cc = Utils::mex_convert_mx_to_ptr<PolyLine>(arg_in_3);
        ptr->push_back( *cc );
      } else {
        UTILS_MEX_ASSERT( false, CMD "unknown type = {}\n", type );
      }
    } else {
      UTILS_MEX_ASSERT( false, CMD "expected 4, 5 or 8 inputs nrhs = {}\n", nrhs );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back_G1(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('push_back_G1',OBJ,[x0,y0,theta0,x1,y1,theta1]|[x1,y1,theta1]): "

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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
      UTILS_MEX_ASSERT( false, CMD "expected 5 or 8 inputs, nrhs = {}\n", nrhs );
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

    #define CMD "ClothoidListMexWrapper('reserve',OBJ,N): "

    UTILS_MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0 , CMD "expected no outputs, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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

    #define CMD "ClothoidListMexWrapper('get_STK',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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

    #define CMD "ClothoidListMexWrapper('get_XY',OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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

    #define CMD "ClothoidListMexWrapper('build_G1', OBJ, x, y [, theta]): "

    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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
  do_build_raw(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('build_raw', OBJ, x, y, ascissa, theta, kappa): "

    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 7, CMD "expected 7 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    bool ok = true;

    mwSize nx, ny;
    real_type const * x = Utils::mex_vector_pointer( arg_in_2, nx, CMD "Error in reading x" );
    real_type const * y = Utils::mex_vector_pointer( arg_in_3, ny, CMD "Error in reading y" );
    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );

    real_type const * a = Utils::mex_vector_pointer( arg_in_4, ny, CMD "Error in reading abscissa" );
    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(abscissa) = {}\n", nx, ny );

    real_type const * t = Utils::mex_vector_pointer( arg_in_5, ny, CMD "Error in reading theta" );
    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(theta) = {}\n", nx, ny );

    real_type const * k = Utils::mex_vector_pointer( arg_in_6, ny, CMD "Error in reading kappa" );
    UTILS_MEX_ASSERT( nx == ny, CMD "length(x) = {} != length(kappa) = {}\n", nx, ny );

    ok = ptr->build_raw( nx, x, y, a, t, k );

    Utils::mex_set_scalar_bool( arg_out_0, ok );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_3arcG2(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('build_3arcG2',OBJ,x0,y0,theta0,kappa0,x1,y1,theta1,kappa1): "

    UTILS_MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading `kappa1`" );

    G2solve3arc g2sol;
    integer iter = g2sol.build( x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    Utils::mex_set_scalar_int32( arg_out_0, iter );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_2arcG2(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('build_2arcG2',OBJ,x0,y0,theta0,kappa0,x1,y1,theta1,kappa1): "

    UTILS_MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading `kappa1`" );

    G2solve2arc g2sol;
    integer iter = g2sol.build( x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(2);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getS1());
    }

    Utils::mex_set_scalar_int32( arg_out_0, iter );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_CLC(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('build_CLC',OBJ,x0,y0,theta0,kappa0,x1,y1,theta1,kappa1): "

    UTILS_MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = Utils::mex_get_scalar_value( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = Utils::mex_get_scalar_value( arg_in_9, CMD "Error in reading `kappa1`" );

    G2solveCLC g2sol;
    integer iter = g2sol.build( x0, y0, theta0, kappa0, x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    Utils::mex_set_scalar_int32( arg_out_0, iter );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_3arcG2fixed(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('build_3arcG2fixed',OBJ,s0,x0,y0,theta0,kappa0,s1,x1,y1,theta1,kappa1): "

    UTILS_MEX_ASSERT( nrhs == 12, CMD "expected 12 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type s0     = Utils::mex_get_scalar_value( arg_in_2,  CMD "Error in reading `s0`" );
    real_type x0     = Utils::mex_get_scalar_value( arg_in_3,  CMD "Error in reading `x0`" );
    real_type y0     = Utils::mex_get_scalar_value( arg_in_4,  CMD "Error in reading `y0`" );
    real_type theta0 = Utils::mex_get_scalar_value( arg_in_5,  CMD "Error in reading `theta0`" );
    real_type kappa0 = Utils::mex_get_scalar_value( arg_in_6,  CMD "Error in reading `kappa0`" );
    real_type s1     = Utils::mex_get_scalar_value( arg_in_7,  CMD "Error in reading `s1`" );
    real_type x1     = Utils::mex_get_scalar_value( arg_in_8,  CMD "Error in reading `x1`" );
    real_type y1     = Utils::mex_get_scalar_value( arg_in_9,  CMD "Error in reading `y1`" );
    real_type theta1 = Utils::mex_get_scalar_value( arg_in_10, CMD "Error in reading `theta1`" );
    real_type kappa1 = Utils::mex_get_scalar_value( arg_in_11, CMD "Error in reading `kappa1`" );

    G2solve3arc g2sol;
    int iter = g2sol.build_fixed_length( s0, x0, y0, theta0, kappa0,
                                         s1, x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    Utils::mex_set_scalar_int32( arg_out_0, iter );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('build',OBJ,x0,y0,theta0,s,kappa): "

    UTILS_MEX_ASSERT( nrhs == 7, CMD "expected 7 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type x0, y0, theta0;
    real_type const * s;
    real_type const * kappa;
    mwSize ns, nk;

    x0     = Utils::mex_get_scalar_value( arg_in_2, CMD "Error in reading `x0`" );
    y0     = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading `y0`" );
    theta0 = Utils::mex_get_scalar_value( arg_in_4, CMD "Error in reading `theta0`" );
    s      = Utils::mex_vector_pointer( arg_in_5, ns, CMD "Error in reading `s1`" );
    kappa  = Utils::mex_vector_pointer( arg_in_6, nk, CMD "Error in reading `x1`" );

    UTILS_MEX_ASSERT( ns == nk, CMD "length(s) = {} != length(kappa) = {}\n", ns, nk );

    bool ok = ptr->build( x0, y0, theta0, ns, s, kappa );

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

    #define CMD "ClothoidListMexWrapper('build_theta',OBJ,x,y): "

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
  do_make_closed( int nlhs, mxArray       *[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('make_closed',OBJ): "

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);
    ptr->make_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_open( int nlhs, mxArray       *[],
                int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('make_open',OBJ): "

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);
    ptr->make_open();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_closed( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('is_closed',OBJ): "

    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);
    Utils::mex_set_scalar_bool( arg_out_0, ptr->is_closed() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_get(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('get',OBJ,n): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 6, CMD "expected 6 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    int64_t n = Utils::mex_get_int64( arg_in_2, CMD "Error in reading n" );

    UTILS_MEX_ASSERT(
      n > 0 && n <= ptr->num_segments(),
      CMD "n = {} must be >= 1 and <= {}\n", n, ptr->num_segments()
    );

    ClothoidCurve const & c = ptr->get(n-1);

    Utils::mex_set_scalar_value( arg_out_0, c.x_begin()     );
    Utils::mex_set_scalar_value( arg_out_1, c.y_begin()     );
    Utils::mex_set_scalar_value( arg_out_2, c.theta_begin() );
    Utils::mex_set_scalar_value( arg_out_3, c.kappa_begin() );
    Utils::mex_set_scalar_value( arg_out_4, c.dkappa()      );
    Utils::mex_set_scalar_value( arg_out_5, c.length()      );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_num_segments(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('num_segments', OBJ): "

    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    Utils::mex_set_scalar_int32( arg_out_0, ptr->num_segments() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_deltaTheta(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('deltaTheta',OBJ): "

    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    integer nseg = ptr->num_segments();

    real_type * dtheta = Utils::mex_create_matrix_value( arg_out_0, nseg, 1 );
    ptr->get_delta_theta( dtheta );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_deltaKappa(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('deltaKappa',OBJ): "

    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    integer nseg = ptr->num_segments();

    real_type * dkappa = Utils::mex_create_matrix_value( arg_out_0, nseg, 1 );
    ptr->get_delta_kappa( dkappa );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_export_table(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('export_table',OBJ,filename ): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nrhs = {}\n", nrhs );

    UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), CMD "filename must be a string" );
    string filename = mxArrayToString(arg_in_2);

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    std::ofstream file(filename.c_str());
    UTILS_MEX_ASSERT( file.good(), CMD " cannot open file: `{}'\n", filename );
    ptr->export_table(file);
    file.close();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_export_ruby(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {

   #define CMD "ClothoidListMexWrapper('export_ruby',OBJ,filename ): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nlhs = {}\n", nlhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected no output, nrhs = {}\n", nrhs );

    UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), CMD "filename must be a string" );
    string filename = mxArrayToString(arg_in_2);

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    std::ofstream file(filename.c_str());
    UTILS_MEX_ASSERT( file.good(), CMD " cannot open file: `{}'\n", filename );
    ptr->export_ruby(file);
    file.close();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_findST1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('findST1',OBJ,x,y): "
    UTILS_MEX_ASSERT(
      nrhs == 4,
      CMD "expected 4 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 3,
      CMD "expected 3 output, nlhs = {}\n", nlhs
    );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

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

  static
  void
  do_closest_segment(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('closest_segment',OBJ,qx,qy): "

    UTILS_MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type qx = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`qx` expected to be a real scalar"
    );
    real_type qy = Utils::mex_get_scalar_value(
      arg_in_3, CMD "`qx` expected to be a real scalar"
    );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->closest_segment( qx, qy ) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_closest_point_in_range(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('closest_point_in_range',OBJ,qx,qy,[ibegin,iend],['ISO'/'SAE']): "

    UTILS_MEX_ASSERT(
      4 <= nrhs && nrhs <= 7,
      CMD "expected 4, 5, 6 or 7 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 1 || nlhs == 7,
      CMD "expected 1 or 7 outputs, nlhs = {}\n", nlhs
    );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type qx = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`qx` expected to be a real scalar"
    );
    real_type qy = Utils::mex_get_scalar_value(
      arg_in_3, CMD "`qx` expected to be a real scalar"
    );
    int64_t icurve_begin = 0;
    int64_t icurve_end   = ptr->num_segments()-1;

    if ( nrhs >= 6 ) {
      icurve_begin = Utils::mex_get_int64(
        arg_in_4, CMD "`ibegin` expected to be an integer"
      );
      icurve_end = Utils::mex_get_int64(
        arg_in_5, CMD "`iend` expected to be an integer"
      );
    }

    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string" );
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string" );

    integer iflag, icurve;
    real_type x, y, s, t, dst;
    if ( ISO ) {
      iflag = ptr->closest_point_in_range_ISO(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
    } else {
      iflag = ptr->closest_point_in_range_SAE(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
    }

    if ( nlhs == 1 ) {
      static char const * fieldnames[] = {
        "icurve", "x", "y", "s", "t", "iflag", "dst"
      };

      mxArray * mx_icurve; Utils::mex_set_scalar_int32( mx_icurve, icurve );
      mxArray * mx_iflag;  Utils::mex_set_scalar_int32( mx_iflag,  iflag  );

      arg_out_0 = mxCreateStructMatrix(1,1,7,fieldnames);
      mxSetFieldByNumber( arg_out_0, 0, 0, mx_icurve );
      mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(x) );
      mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(y) );
      mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(s) );
      mxSetFieldByNumber( arg_out_0, 0, 4, mxCreateDoubleScalar(t) );
      mxSetFieldByNumber( arg_out_0, 0, 5, mx_iflag );
      mxSetFieldByNumber( arg_out_0, 0, 6, mxCreateDoubleScalar(dst) );

    } else {
      Utils::mex_set_scalar_int32  ( arg_out_0, icurve );
      Utils::mex_set_scalar_value( arg_out_1, x      );
      Utils::mex_set_scalar_value( arg_out_2, y      );
      Utils::mex_set_scalar_value( arg_out_3, s      );
      Utils::mex_set_scalar_value( arg_out_4, t      );
      Utils::mex_set_scalar_int32  ( arg_out_5, iflag  );
      Utils::mex_set_scalar_value( arg_out_6, dst    );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_closest_point_in_s_range(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('closest_point_in_s_range',OBJ,qx,qy,s_begin,s_end,['ISO'/'SAE']): "

    UTILS_MEX_ASSERT(
      4 <= nrhs && nrhs <= 7,
      CMD "expected 4, 5, 6 or 7 inputs, nrhs = {}\n", nrhs
    );
    UTILS_MEX_ASSERT(
      nlhs == 1 || nlhs == 7,
      CMD "expected 1 or 7 output, nlhs = {}\n", nlhs
    );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type qx = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`qx` expected to be a real scalar"
    );
    real_type qy = Utils::mex_get_scalar_value(
      arg_in_3, CMD "`qx` expected to be a real scalar"
    );
    int64_t s_curve_begin = 0;
    int64_t s_curve_end   = ptr->length();

    if ( nrhs >= 6 ) {
      s_curve_begin = Utils::mex_get_scalar_value(
        arg_in_4, CMD "`s_begin` expected to be a scalar"
      );
      s_curve_end = Utils::mex_get_scalar_value(
        arg_in_5, CMD "`s_end` expected to be a scalar"
      );
    }

    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string" );
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string" );

    integer iflag, icurve;
    real_type x, y, s, t, dst;
    if ( ISO ) {
      iflag = ptr->closest_point_in_s_range_ISO(
        qx, qy, s_curve_begin, s_curve_end, x, y, s, t, dst, icurve
      );
    } else {
      iflag = ptr->closest_point_in_s_range_SAE(
        qx, qy, s_curve_begin, s_curve_end, x, y, s, t, dst, icurve
      );
    }

    if ( nlhs == 1 ) {
      static char const * fieldnames[] = {
        "icurve", "x", "y", "s", "t", "iflag", "dst"
      };

      mxArray * mx_icurve; Utils::mex_set_scalar_int32( mx_icurve, icurve );
      mxArray * mx_iflag;  Utils::mex_set_scalar_int32( mx_iflag,  iflag  );

      arg_out_0 = mxCreateStructMatrix(1,1,7,fieldnames);
      mxSetFieldByNumber( arg_out_0, 0, 0, mx_icurve );
      mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(x) );
      mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(y) );
      mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(s) );
      mxSetFieldByNumber( arg_out_0, 0, 4, mxCreateDoubleScalar(t) );
      mxSetFieldByNumber( arg_out_0, 0, 5, mx_iflag );
      mxSetFieldByNumber( arg_out_0, 0, 6, mxCreateDoubleScalar(dst) );

    } else {
      Utils::mex_set_scalar_int32  ( arg_out_0, icurve );
      Utils::mex_set_scalar_value( arg_out_1, x      );
      Utils::mex_set_scalar_value( arg_out_2, y      );
      Utils::mex_set_scalar_value( arg_out_3, s      );
      Utils::mex_set_scalar_value( arg_out_4, t      );
      Utils::mex_set_scalar_int32  ( arg_out_5, iflag  );
      Utils::mex_set_scalar_value( arg_out_6, dst    );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_s_to_index(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('s_to_index',OBJ,s): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    real_type s = Utils::mex_get_scalar_value(
      arg_in_2, CMD "`s` expected to be a real scalar"
    );

    Utils::mex_set_scalar_int32( arg_out_0, ptr->find_at_s( s ) );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_load(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('load',OBJ,filename,[epsi]): "
    UTILS_MEX_ASSERT(
      nrhs == 3 || nrhs == 4,
      CMD "expected 3 or 4 inputs, nrhs = {}\n", nrhs
    );

    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), "Third argument must be a string" );

    std::string fname = mxArrayToString(arg_in_2);
    std::ifstream file(fname.c_str());

    if ( nrhs == 4 ) {
      real_type epsi = Utils::mex_get_scalar_value( arg_in_3, CMD "Error in reading epsi" );
      ptr->load( file, epsi );
    } else {
      ptr->load( file );
    }
    file.close();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_save(
    int nlhs, mxArray       *[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('save',OBJ,filename): "

    UTILS_MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    UTILS_MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = {}\n", nlhs );

    ClothoidList * ptr = Utils::mex_convert_mx_to_ptr<ClothoidList>(arg_in_1);

    UTILS_MEX_ASSERT0( mxIsChar(arg_in_2), "Third argument must be a string" );

    std::string fname = mxArrayToString(arg_in_2);
    std::ofstream file(fname.c_str());
    ptr->save( file );
    file.close();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun{
    {"new",do_new},
    {"push_back",do_push_back},
    {"push_back_G1",do_push_back_G1},
    {"reserve",do_reserve},
    {"get_STK",do_get_STK},
    {"get_XY",do_get_XY},
    {"build_G1",do_build_G1},
    {"build_raw",do_build_raw},
    {"build_3arcG2",do_build_3arcG2},
    {"build_2arcG2",do_build_2arcG2},
    {"build_CLC",do_build_CLC},
    {"build_3arcG2fixed",do_build_3arcG2fixed},
    {"build",do_build},
    {"build_theta",do_build_theta},
    {"make_closed",do_make_closed},
    {"make_open",do_make_open},
    {"is_closed",do_is_closed},
    {"get",do_get},
    {"num_segments",do_num_segments},
    {"deltaTheta",do_deltaTheta},
    {"deltaKappa",do_deltaKappa},
    {"export_table",do_export_table},
    {"export_ruby",do_export_ruby},
    {"findST1",do_findST1},
    {"closest_segment",do_closest_segment},
    {"closest_point_in_range",do_closest_point_in_range},
    {"closest_point_in_s_range",do_closest_point_in_s_range},
    {"s_to_index",do_s_to_index},
    {"save",do_save},
    {"load",do_load},
    CMD_MAP_FUN
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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
      mexErrMsgTxt( fmt::format( "ClothoidList Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("ClothoidList failed\n");
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

}
