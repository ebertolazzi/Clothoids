/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2019
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoids.hh"
#include "mex_utils.hh"
#include "mex_info.hxx"

#include <fstream>

#include "mex_Workaround.hxx"

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
"    [s,theta,kappa] = BiarcListMexWrapper( 'getSTK', OBJ );\n" \
"    [x,y]           = BiarcListMexWrapper( 'getXY', OBJ );\n" \
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

  static
  BiarcList *
  DATA_NEW( mxArray * & mx_id ) {
    BiarcList * ptr = new BiarcList();
    mx_id = convertPtr2Mat<BiarcList>(ptr);
    return ptr;
  }

  static
  inline
  BiarcList *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<BiarcList>(mx_id);
  }

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
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('new'): "
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 1, CMD "expected 1 input, nlhs = {}\n", nrhs );

    //BiarcList * ptr =
    DATA_NEW(arg_out_0);

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back_G1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('push_back_G1',OBJ,[x0,y0,theta0,x1,y1,theta1]|[x1,y1,theta1]): "

    MEX_ASSERT2( nlhs == 0, CMD "expected NO output, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    if ( nrhs == 8 ) {
      real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
      real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
      real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
      real_type x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" );
      real_type y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" );
      real_type theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" );
      ptr->push_back_G1( x0, y0, theta0, x1, y1, theta1 );
    } else if ( nrhs == 5 ) {
      real_type x1     = getScalarValue( arg_in_2, CMD "Error in reading x1" );
      real_type y1     = getScalarValue( arg_in_3, CMD "Error in reading y1" );
      real_type theta1 = getScalarValue( arg_in_4, CMD "Error in reading theta1" );
      ptr->push_back_G1( x1, y1, theta1 );
    } else {
      MEX_ASSERT2( false, CMD "expected 5 or 8 inputs, nrhs = {|\n", nrhs );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_reserve(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('reserve',OBJ,N): "

    MEX_ASSERT2( nrhs == 3 , CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0 , CMD "expected no outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    int64_t N = getInt( arg_in_2, CMD "Error in reading N" );
    ptr->reserve( N );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_getSTK(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('getSTK',OBJ): "

    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 3, CMD "expected 3 outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    int_type n = ptr->numSegments();

    real_type * s     = createMatrixValue( arg_out_0, 1, n+1 );
    real_type * theta = createMatrixValue( arg_out_1, 1, n+1 );
    real_type * kappa = createMatrixValue( arg_out_2, 1, n+1 );

    ptr->getSTK( s, theta, kappa );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_getXY(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "BiarcListMexWrapper('getXY',OBJ): "

    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 2, CMD "expected 2 outputs, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    int_type    n = ptr->numSegments();
    real_type * x = createMatrixValue( arg_out_0, 1, n+1 );
    real_type * y = createMatrixValue( arg_out_1, 1, n+1 );

    ptr->getXY( x, y );

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

    MEX_ASSERT2(nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs);

    BiarcList * ptr = DATA_GET(arg_in_1);

    bool ok = true;

    if ( nrhs == 4 ) {
      mwSize nx, ny;
      real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

      MEX_ASSERT2( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );

      ok = ptr->build_G1( nx, x, y );

    } else if ( nrhs == 5 ) {

      mwSize nx, ny, nt;
      real_type const * x     = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y     = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );
      real_type const * theta = getVectorPointer( arg_in_4, nt, CMD "Error in reading theta" );

      MEX_ASSERT2( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );
      MEX_ASSERT2( nx == nt, CMD "length(theta) = {} != length(x) = length(y) = {}\n", nt, ny );

      ok = ptr->build_G1( nx, x, y, theta );

    } else {
      MEX_ASSERT2( false, CMD "expected 4 or 5 input, nrhs = {}\n", nrhs );
    }

    setScalarBool( arg_out_0, ok );
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

    MEX_ASSERT2( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );

    mwSize nx, ny;
    real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
    real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

    MEX_ASSERT2( nx == ny, CMD "length(x) = {} != length(y) = {}\n", nx, ny );

    real_type * theta = createMatrixValue( arg_out_0, nx, 1 );

    bool ok = build_guess_theta( nx, x, y, theta );

    setScalarBool( arg_out_1, ok );

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

    MEX_ASSERT2(nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2(nlhs == 6, CMD "expected 6 output, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    int64_t n = getInt( arg_in_2, CMD "Error in reading n" );

    MEX_ASSERT2(
      n > 0 && n <= ptr->numSegments(),
      CMD "n = {} must be >= 1 and <= {}\n", n, ptr->numSegments()
    );

    Biarc const & c = ptr->get(n-1);

    setScalarValue(arg_out_0, c.xBegin());
    setScalarValue(arg_out_1, c.yBegin());
    setScalarValue(arg_out_2, c.thetaBegin());
    setScalarValue(arg_out_3, c.xEnd());
    setScalarValue(arg_out_4, c.yEnd());
    setScalarValue(arg_out_5, c.thetaEnd());

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_numSegments(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "BiarcListMexWrapper('numSegments', OBJ): "

    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    BiarcList * ptr = DATA_GET(arg_in_1);

    setScalarInt( arg_out_0, ptr->numSegments() );

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
    MEX_ASSERT2(
      nrhs == 4,
      CMD "expected 4 inputs, nrhs = {}\n", nrhs
    );
    MEX_ASSERT2(
      nlhs == 3,
      CMD "expected 3 output, nlhs = {}\n", nlhs
    );

    BiarcList * ptr = DATA_GET(arg_in_1);

    mwSize nrx, ncx, nry, ncy;
    real_type const * x;
    real_type const * y;

    x = getMatrixPointer(
      arg_in_2, nrx, ncx,
      CMD "`x` expected to be a real vector/matrix"
    );

    y = getMatrixPointer(
      arg_in_3, nry, ncy,
      CMD "`y` expected to be a real vector/matrix"
    );

    MEX_ASSERT2(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found\n"
      "size(x) = {} x {} size(y) = {} x {}\n",
      nrx, ncx, nry, ncy
    );

    real_type * s   = createMatrixValue( arg_out_0, nrx, ncx );
    real_type * t   = createMatrixValue( arg_out_1, nrx, ncx );
    real_type * idx = createMatrixValue( arg_out_2, nrx, ncx );

    mwSize size = nrx*ncx;
    for ( mwSize i = 0; i < size; ++i ) {
      int_type nseg = ptr->findST1( *x++, *y++, *s++, *t++ );
      *idx++ = nseg >= 0 ? nseg+1 : nseg;
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"push_back_G1",do_push_back_G1},
    {"reserve",do_reserve},
    {"getSTK",do_getSTK},
    {"getXY",do_getXY},
    {"build_G1",do_build_G1},
    {"build_theta",do_build_theta},
    {"get",do_get},
    {"numSegments",do_numSegments},
    {"findST1",do_findST1},
    CMD_MAP_FUN
  };

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
      MEX_ASSERT( mxIsChar(arg_in_0), "First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);
      DO_CMD pfun = cmd_to_fun.at(cmd);
      pfun( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "BiarcList Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt( "BiarcList failed" );
    }
  }
}
