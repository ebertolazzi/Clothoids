/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoids.hh"
#include "mex_utils.hh"
#include "mex_info.hxx"

#include "mex_Workaround.hxx"

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
"    [x0,y0,theta0,kappa0,L0,x1,y1,theta1,kappa1,L1] = BiarcMexWrapper( 'getPars', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("BiarcMexWrapper") \
MEX_INFO_MESSAGE_END

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

  static
  Biarc *
  DATA_NEW( mxArray * & mx_id ) {
    Biarc * ptr = new Biarc();
    mx_id = convertPtr2Mat<Biarc>(ptr);
    return ptr;
  }

  static
  inline
  Biarc *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<Biarc>(mx_id);
  }

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
          int nrhs, mxArray const *prhs[] ) {

    #define CMD "BiarcMexWrapper('new'): "
    MEX_ASSERT2( nrhs == 1, CMD "expected 1 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    #undef CMD

    DATA_NEW( arg_out_0 );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_G1( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    #define CMD "BiarcMexWrapper('build',OBJ,x0,y0,theta0,x1,y1,theta1): "

    MEX_ASSERT2( nrhs == 8, CMD "expected 8 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Biarc * ptr = DATA_GET(arg_in_1);

    real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
    real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
    real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
    real_type x1     = getScalarValue( arg_in_5, CMD "Error in reading x1" );
    real_type y1     = getScalarValue( arg_in_6, CMD "Error in reading y1" );
    real_type theta1 = getScalarValue( arg_in_7, CMD "Error in reading theta1" );

    bool ok = ptr->build( x0, y0, theta0, x1, y1, theta1 );

    // returns the status of the interpolation
    setScalarBool( arg_out_0, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build_3P( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('build_3P',OBJ,...): "
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    #undef CMD

    if ( nrhs == 8 ) {
      #define CMD "BiarcMexWrapper('build_3P',OBJ,x0,y0,x1,y1,x2,y2): "

      real_type x0 = getScalarValue( arg_in_2, CMD "Error in reading x0" );
      real_type y0 = getScalarValue( arg_in_3, CMD "Error in reading y0" );
      real_type x1 = getScalarValue( arg_in_4, CMD "Error in reading x1" );
      real_type y1 = getScalarValue( arg_in_5, CMD "Error in reading y1" );
      real_type x2 = getScalarValue( arg_in_6, CMD "Error in reading x2" );
      real_type y2 = getScalarValue( arg_in_7, CMD "Error in reading y2" );

      bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 );

      // returns the status of the interpolation
      setScalarBool(arg_out_0,ok);

      #undef CMD
    } else if ( nrhs == 5 ) {
      #define CMD "BiarcMexWrapper('build_3P',OBJ,p0,p1,p2): "
      real_type const * p0;
      real_type const * p1;
      real_type const * p2;

      mwSize n;
      p0 = getVectorPointer( arg_in_2, n, CMD "Error in reading p0" );
      MEX_ASSERT2(
        n == 2,
        CMD "Error in reading length(p0) == {} expect length(p0) == 2\n", n
      );
      p1 = getVectorPointer( arg_in_3, n, CMD "Error in reading p1" );
      MEX_ASSERT2(
        n == 2,
        CMD "Error in reading length(p1) == {} expect length(p1) == 2\n", n
      );
      p2 = getVectorPointer( arg_in_4, n, CMD "Error in reading p2" );
      MEX_ASSERT2(
        n == 2,
        CMD "Error in reading length(p2) == {} expect length(p2) == 2\n", n
      );

      bool ok = ptr->build_3P( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );

      // returns the status of the interpolation
      setScalarBool(arg_out_0,ok);

      #undef CMD
    } else {
      MEX_ASSERT(
        false, "BiarcMexWrapper('build_3P',OBJ,...) expected 5 or 8 arguments\n"
      );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_x_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('xMiddle',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->xMiddle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_y_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('yMiddle',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->yMiddle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_theta_middle( int nlhs, mxArray       *plhs[],
                   int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('thetaMiddle',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->thetaMiddle());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_s_middle( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('sMiddle',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->length0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_kappa0( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('kappa0',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->kappa0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_kappa1( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('kappa1',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->kappa1());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_length0( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('length0',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->length0());
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_length1( int nlhs, mxArray       *plhs[],
              int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('length1',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );
    setScalarValue( arg_out_0, ptr->length1());;
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_to_nurbs( int nlhs, mxArray       *plhs[],
               int nrhs, mxArray const *prhs[] ) {

    Biarc * ptr = DATA_GET(arg_in_1);

    #define CMD "BiarcMexWrapper('to_nurbs',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 2, CMD "expected 2 output, nlhs = {}\n", nlhs );

    CircleArc const & C0 = ptr->C0();
    CircleArc const & C1 = ptr->C1();

    int_type npts0, nknots0, npts1, nknots1;
    C0.paramNURBS( nknots0, npts0 );
    C1.paramNURBS( nknots1, npts1 );

    mxArray * mx_knots0, * mx_Poly0, * mx_knots1, * mx_Poly1;

    double * knots0 = createMatrixValue( mx_knots0, 1, nknots0 );
    double * poly0  = createMatrixValue( mx_Poly0,  3, npts0   );
    double * knots1 = createMatrixValue( mx_knots1, 1, nknots1 );
    double * poly1  = createMatrixValue( mx_Poly1,  3, npts1   );

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

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    //{"build",do_build},
    {"build_3P",do_build_3P},
    {"build_G1",do_build_G1},
    {"xMiddle",do_x_middle},
    {"yMiddle",do_y_middle},
    {"thetaMiddle",do_theta_middle},
    {"sMiddle",do_s_middle},
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
      mexErrMsgTxt( fmt::format( "Biarc Error: {}", e.what() ).c_str() );
    } catch (...) {
      mexErrMsgTxt("Biarc failed\n");
    }

  }
}
