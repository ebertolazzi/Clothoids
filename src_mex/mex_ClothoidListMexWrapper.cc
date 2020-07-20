/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Line.hh"
#include "Circle.hh"
#include "Biarc.hh"
#include "PolyLine.hh"
#include "ClothoidList.hh"
#include "Triangle2D.hh"

#include "mex_utils.hh"
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
"    [s,theta,kappa] = ClothoidListMexWrapper( 'getSTK', OBJ );\n" \
"    [x,y]           = ClothoidListMexWrapper( 'getXY', OBJ );\n" \
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
"    [theta,ok] = ClothoidListMexWrapper( 'build_theta', OBJ, x, y );%\n" \
"    dtheta = ClothoidListMexWrapper( 'deltaTheta', OBJ );%\n" \
"    dkappa = ClothoidListMexWrapper( 'deltaKappa', OBJ );%\n" \
"\n" \
MEX_INFO_MESSAGE("ClothoidListMexWrapper") \
MEX_INFO_MESSAGE_END



namespace G2lib {

  using namespace std;

  /*\
   |  ____    _  _____  _
   | |  _ \  / \|_   _|/ \
   | | | | |/ _ \ | | / _ \
   | | |_| / ___ \| |/ ___ \
   | |____/_/   \_\_/_/   \_\
   |
  \*/

  static
  ClothoidList *
  DATA_NEW( mxArray * & mx_id ) {
    ClothoidList * ptr = new ClothoidList();
    mx_id = convertPtr2Mat<ClothoidList>(ptr);
    return ptr;
  }

  static
  inline
  ClothoidList *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<ClothoidList>(mx_id);
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

    #define CMD "ClothoidListMexWrapper('new'): "
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nlhs = " << nrhs );

    //ClothoidList * ptr =
    DATA_NEW(arg_out_0);

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('push_back',OBJ,CLOT|[kappa0,dkappa,L]|[x0,y0,theta0,kappa0,dkappa,L]): "

    MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    if ( nrhs == 8 ) {
      real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading x0" );
      real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading y0" );
      real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading theta0" );
      real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading kappa0" );
      real_type dkappa = getScalarValue( arg_in_6, CMD "Error in reading dkappa" );
      real_type L      = getScalarValue( arg_in_7, CMD "Error in reading L" );
      ptr->push_back( x0, y0, theta0, kappa0, dkappa, L );
    } else if ( nrhs == 5 ) {
      real_type kappa0 = getScalarValue( arg_in_2, CMD "Error in reading kappa0" );
      real_type dkappa = getScalarValue( arg_in_3, CMD "Error in reading dkappa" );
      real_type L      = getScalarValue( arg_in_4, CMD "Error in reading L" );
      ptr->push_back( kappa0, dkappa, L );
    } else if ( nrhs == 3 ) {
      ClothoidCurve * cc = convertMat2Ptr<ClothoidCurve>(arg_in_2);
      ptr->push_back( *cc );
    } else {
      MEX_ASSERT( false, CMD "expected 3, 5 or 8 inputs nrhs = " << nrhs );
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_push_back_G1(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('push_back_G1',OBJ,[x0,y0,theta0,x1,y1,theta1]|[CLOT]): "

    MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

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
      MEX_ASSERT( false, CMD "expected 5 or 8 inputs, nrhs = " << nrhs );
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

    #define CMD "ClothoidListMexWrapper('reserve',OBJ,N): "

    MEX_ASSERT( nrhs == 3 , CMD "expected 3 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 0 , CMD "expected no outputs, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

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

    #define CMD "ClothoidListMexWrapper('getSTK',OBJ): "

    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 3, CMD "expected 3 outputs, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    int_type n = ptr->numSegment();

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

    #define CMD "ClothoidListMexWrapper('getXY',OBJ): "

    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 outputs, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    int_type    n = ptr->numSegment();
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

    #define CMD "ClothoidListMexWrapper('build_G1', OBJ, x, y [, theta]): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    bool ok = true;

    if ( nrhs == 4 ) {
      mwSize nx, ny;
      real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

      MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " != length(y) = " << ny );

      ok = ptr->build_G1( nx, x, y );

    } else if ( nrhs == 5 ) {

      mwSize nx, ny, nt;
      real_type const * x     = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
      real_type const * y     = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );
      real_type const * theta = getVectorPointer( arg_in_4, nt, CMD "Error in reading theta" );

      MEX_ASSERT(
        nx == ny,
        CMD "length(x) = " << nx << " != length(y) = " << ny
      );
      MEX_ASSERT(
        nx == nt,
        CMD "length(theta) = " << nt << " != length(x) = length(y) = " << ny
      );

      ok = ptr->build_G1( nx, x, y, theta );

    } else {
      MEX_ASSERT( false, CMD "expected 4 or 5 input, nrhs = " << nrhs );
    }

    setScalarBool( arg_out_0, ok );
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

    MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = getScalarValue( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = getScalarValue( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = getScalarValue( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = getScalarValue( arg_in_9, CMD "Error in reading `kappa1`" );

    static G2solve3arc g2sol;
    int_type iter = g2sol.build( x0, y0, theta0, kappa0,
                                 x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    setScalarInt( arg_out_0, iter );

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

    MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = getScalarValue( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = getScalarValue( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = getScalarValue( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = getScalarValue( arg_in_9, CMD "Error in reading `kappa1`" );

    static G2solve2arc g2sol;
    int_type iter = g2sol.build( x0, y0, theta0, kappa0,
                                 x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(2);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getS1());
    }

    setScalarInt( arg_out_0, iter );

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

    MEX_ASSERT( nrhs == 10, CMD "expected 10 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type x0     = getScalarValue( arg_in_2, CMD "Error in reading `x0`" );
    real_type y0     = getScalarValue( arg_in_3, CMD "Error in reading `y0`" );
    real_type theta0 = getScalarValue( arg_in_4, CMD "Error in reading `theta0`" );
    real_type kappa0 = getScalarValue( arg_in_5, CMD "Error in reading `kappa0`" );
    real_type x1     = getScalarValue( arg_in_6, CMD "Error in reading `x1`" );
    real_type y1     = getScalarValue( arg_in_7, CMD "Error in reading `y1`" );
    real_type theta1 = getScalarValue( arg_in_8, CMD "Error in reading `theta1`" );
    real_type kappa1 = getScalarValue( arg_in_9, CMD "Error in reading `kappa1`" );

    static G2solveCLC g2sol;
    int_type iter = g2sol.build( x0, y0, theta0, kappa0,
                                 x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    setScalarInt( arg_out_0, iter );

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

    MEX_ASSERT( nrhs == 12, CMD "expected 12 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1,  CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type s0     = getScalarValue( arg_in_2,  CMD "Error in reading `s0`" );
    real_type x0     = getScalarValue( arg_in_3,  CMD "Error in reading `x0`" );
    real_type y0     = getScalarValue( arg_in_4,  CMD "Error in reading `y0`" );
    real_type theta0 = getScalarValue( arg_in_5,  CMD "Error in reading `theta0`" );
    real_type kappa0 = getScalarValue( arg_in_6,  CMD "Error in reading `kappa0`" );
    real_type s1     = getScalarValue( arg_in_7,  CMD "Error in reading `s1`" );
    real_type x1     = getScalarValue( arg_in_8,  CMD "Error in reading `x1`" );
    real_type y1     = getScalarValue( arg_in_9,  CMD "Error in reading `y1`" );
    real_type theta1 = getScalarValue( arg_in_10, CMD "Error in reading `theta1`" );
    real_type kappa1 = getScalarValue( arg_in_11, CMD "Error in reading `kappa1`" );

    static G2solve3arc g2sol;
    int iter = g2sol.build_fixed_length( s0, x0, y0, theta0, kappa0,
                                         s1, x1, y1, theta1, kappa1 );
    if ( iter >= 0 ) {
      ptr->init();
      ptr->reserve(3);
      ptr->push_back(g2sol.getS0());
      ptr->push_back(g2sol.getSM());
      ptr->push_back(g2sol.getS1());
    }

    setScalarInt( arg_out_0, iter );

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

    MEX_ASSERT( nrhs == 7, CMD "expected 7 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type x0, y0, theta0;
    real_type const * s;
    real_type const * kappa;
    mwSize ns, nk;

    x0     = getScalarValue( arg_in_2, CMD "Error in reading `x0`" );
    y0     = getScalarValue( arg_in_3, CMD "Error in reading `y0`" );
    theta0 = getScalarValue( arg_in_4, CMD "Error in reading `theta0`" );
    s      = getVectorPointer( arg_in_5, ns, CMD "Error in reading `s1`" );
    kappa  = getVectorPointer( arg_in_6, nk, CMD "Error in reading `x1`" );

    MEX_ASSERT( ns == nk, CMD "length(s) = " << ns << " != length(kappa) = " << nk );

    bool ok = ptr->build( x0, y0, theta0, ns, s, kappa );

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

    #define CMD "ClothoidListMexWrapper('build_theta',OBJ,x,y): "

    MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, CMD "expected 4 input, nrhs = " << nrhs );

    mwSize nx, ny;
    real_type const * x = getVectorPointer( arg_in_2, nx, CMD "Error in reading x" );
    real_type const * y = getVectorPointer( arg_in_3, ny, CMD "Error in reading y" );

    MEX_ASSERT( nx == ny, CMD "length(x) = " << nx << " != length(y) = " << ny );

    real_type * theta = createMatrixValue( arg_out_0, nx, 1 );

    bool ok = build_guess_theta( nx, x, y, theta );

    setScalarBool( arg_out_1, ok );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_closed( int nlhs, mxArray       *plhs[],
                  int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('make_closed',OBJ): "

    MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    ClothoidList * ptr = DATA_GET( arg_in_1 );
    ptr->make_closed();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_make_open( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('make_open',OBJ): "

    MEX_ASSERT( nlhs == 0, CMD "expected 0 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    ClothoidList * ptr = DATA_GET( arg_in_1 );
    ptr->make_open();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_is_closed( int nlhs, mxArray       *plhs[],
                int nrhs, mxArray const *prhs[] ) {

    #define CMD "ClothoidListMexWrapper('is_closed',OBJ): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    ClothoidList * ptr = DATA_GET( arg_in_1 );
    setScalarBool( arg_out_0, ptr->is_closed() );

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

    MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 6, CMD "expected 6 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    int64_t n = getInt( arg_in_2, CMD "Error in reading n" );

    MEX_ASSERT(
      n > 0 && n <= ptr->numSegment(),
      CMD "n = " << n << " must be >= 1 and <= " << ptr->numSegment()
    );

    ClothoidCurve const & c = ptr->get(n-1);

    setScalarValue( arg_out_0, c.xBegin()     );
    setScalarValue( arg_out_1, c.yBegin()     );
    setScalarValue( arg_out_2, c.thetaBegin() );
    setScalarValue( arg_out_3, c.kappaBegin() );
    setScalarValue( arg_out_4, c.dkappa()     );
    setScalarValue( arg_out_5, c.length()     );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_numSegment(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "ClothoidListMexWrapper('numSegment', OBJ): "

    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    setScalarInt( arg_out_0, ptr->numSegment() );

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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    int_type nseg = ptr->numSegment();

    real_type * dtheta = createMatrixValue( arg_out_0, nseg, 1 );
    ptr->getDeltaTheta( dtheta );

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

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    int_type nseg = ptr->numSegment();

    real_type * dkappa = createMatrixValue( arg_out_0, nseg, 1 );
    ptr->getDeltaKappa( dkappa );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_export_table(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('export_table',OBJ,filename ): "

    MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nlhs = " << nlhs );
    MEX_ASSERT( nlhs == 0, CMD "expected no output, nrhs = " << nrhs );

    MEX_ASSERT( mxIsChar(arg_in_2), CMD "filename must be a string" );
    string filename = mxArrayToString(arg_in_2);

    ClothoidList * ptr = DATA_GET(arg_in_1);

    std::ofstream file(filename.c_str());
    MEX_ASSERT( file.good(), CMD " cannot open file: `" << filename << "`" );
    ptr->export_table(file);
    file.close();

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_export_ruby(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

   #define CMD "ClothoidListMexWrapper('export_ruby',OBJ,filename ): "

    MEX_ASSERT( nrhs == 3, CMD "expected 3 inputs, nlhs = " << nlhs );
    MEX_ASSERT( nlhs == 0, CMD "expected no output, nrhs = " << nrhs );

    MEX_ASSERT( mxIsChar(arg_in_2), CMD "filename must be a string" );
    string filename = mxArrayToString(arg_in_2);

    ClothoidList * ptr = DATA_GET(arg_in_1);

    std::ofstream file(filename.c_str());
    MEX_ASSERT( file.good(), CMD " cannot open file: `" << filename << "`" );
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

    #define CMD "ClothoidListMexWrapper('findST1',OBJ,x,y[,ibegin,iend]): "
    MEX_ASSERT(
      nrhs == 4 || nrhs == 6,
      CMD "expected 4 or 6 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 3,
      CMD "expected 3 output, nlhs = " << nlhs
    );

    ClothoidList * ptr = DATA_GET(arg_in_1);

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

    MEX_ASSERT(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found size(x) = " <<
      nrx << " x " << nry << " size(y) = " << nry << " x " << ncy
    );

    int64_t ibegin = 0;
    int64_t iend   = ptr->numSegment()-1;
    if ( nrhs == 6 ) {
      ibegin = getInt(
        arg_in_4, CMD "`ibegin` expected to be a scalar integer"
      );
      iend = getInt(
        arg_in_5, CMD "`iend` expected to be a scalar integer"
      );
    }

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

  #define CMD_BASE "ClothoidListMexWrapper"
  #define G2LIB_CLASS ClothoidList
  #include "mex_common.hxx"
  #undef CMD_BASE
  #undef G2LIB_CLASS

  static
  void
  do_bbTriangles(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('bbTriangles',OBJ[,max_angle,max_size,offs,'ISO'/'SAE']): "

    MEX_ASSERT(
      nrhs >= 2 && nrhs <= 6,
      CMD "expected 2 up to 6 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 3,
      CMD "expected 3 output, nlhs = " << nlhs
    );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type max_angle = m_pi/18;
    real_type max_size  = 1e100;
    real_type offs      = 0;
    if ( nrhs >= 3 )
      max_angle = getScalarValue(
        arg_in_2, CMD "`max_angle` expected to be a real scalar"
      );
    if ( nrhs >= 4 )
      max_size = getScalarValue(
        arg_in_3, CMD "`max_size` expected to be a real scalar"
      );
    if ( nrhs >= 5 )
      offs = getScalarValue(
        arg_in_4, CMD "`offs` expected to be a real scalar"
      );

    bool ISO = true;
    if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_5, CMD " last argument must be a string");

    std::vector<Triangle2D> tvec;
    if ( nrhs >= 5 ) {
      if ( ISO ) ptr->bbTriangles_ISO( offs, tvec, max_angle, max_size );
      else       ptr->bbTriangles_SAE( offs, tvec, max_angle, max_size );
    } else {
      ptr->bbTriangles( tvec, max_angle, max_size );
    }

    mwSize nt = tvec.size();

    real_type * p0 = createMatrixValue( arg_out_0, 2, nt );
    real_type * p1 = createMatrixValue( arg_out_1, 2, nt );
    real_type * p2 = createMatrixValue( arg_out_2, 2, nt );

    for ( mwSize i = 0; i < nt; ++i ) {
      Triangle2D const & t = tvec[i];
      *p0++ = t.x1();
      *p0++ = t.y1();
      *p1++ = t.x2();
      *p1++ = t.y2();
      *p2++ = t.x3();
      *p2++ = t.y3();
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_aabb_true(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('aabb_true',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
    G2lib::yesAABBtree();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_aabb_false(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('aabb_false',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );
    G2lib::noAABBtree();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_closestSegment(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('closestSegment',OBJ,qx,qy): "

    MEX_ASSERT( nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type qx = getScalarValue(
      arg_in_2, CMD "`qx` expected to be a real scalar"
    );
    real_type qy = getScalarValue(
      arg_in_3, CMD "`qx` expected to be a real scalar"
    );

    setScalarInt( arg_out_0, ptr->closestSegment( qx, qy ) );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_closestPointInRange(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "ClothoidListMexWrapper('closestPointInRange',OBJ,qx,qy,[ibegin,iend],['ISO'/'SAE']): "

    MEX_ASSERT(
      4 <= nrhs && nrhs <= 7,
      CMD "expected 4, 5, 6 or 7 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 7,
      CMD "expected 1 output, nlhs = " << nlhs
    );

    ClothoidList * ptr = DATA_GET(arg_in_1);

    real_type qx = getScalarValue(
      arg_in_2, CMD "`qx` expected to be a real scalar"
    );
    real_type qy = getScalarValue(
      arg_in_3, CMD "`qx` expected to be a real scalar"
    );
    int64_t icurve_begin = 0;
    int64_t icurve_end   = ptr->numSegment()-1;

    if ( nrhs >= 6 ) {
      icurve_begin = getInt(
        arg_in_4, CMD "`ibegin` expected to be an integer"
      );
      icurve_end = getInt(
        arg_in_5, CMD "`iend` expected to be an integer"
      );
    }

    bool ISO = true;
    if ( nrhs == 7 ) ISO = do_is_ISO( arg_in_6, CMD " last argument must be a string" );
    if ( nrhs == 5 ) ISO = do_is_ISO( arg_in_4, CMD " last argument must be a string" );

    int_type  iflag, icurve;
    real_type x, y, s, t, dst;
    if ( ISO ) {
      iflag = ptr->closestPointInRange_ISO(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
    } else {
      iflag = ptr->closestPointInRange_SAE(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
    }

    setScalarInt( arg_out_0, icurve );
    setScalarValue( arg_out_1, x );
    setScalarValue( arg_out_2, y );
    setScalarValue( arg_out_3, s );
    setScalarValue( arg_out_4, t );
    setScalarInt( arg_out_5, iflag );
    setScalarValue( arg_out_6, dst );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_NEW,
    CMD_PUSH_BACK,
    CMD_PUSH_BACK_G1,
    CMD_RESERVE,
    CMD_GET_STK,
    CMD_GET_XY,
    CMD_BUILD_G1,
    CMD_3_ARC_G2,
    CMD_2_ARC_G2,
    CMD_3_ARC_CLC,
    CMD_3_ARC_G2_FIXED,
    CMD_BUILD,
    CMD_BUILD_THETA,
    CMD_MAKE_CLOSED,
    CMD_MAKE_OPEN,
    CMD_IS_CLOSED,
    CMD_GET,
    CMD_NUM_SEGMENT,
    CMD_DELTA_THETA,
    CMD_DELTA_KAPPA,
    CMD_EXPORT_TABLE,
    CMD_EXPORT_RUBY,
    CMD_FIND_ST1,
    CMD_BB_TRIANGLES,
    CMD_AABB_TRUE,
    CMD_AABB_FALSE,
    CMD_CLOSEST_SEGMENT,
    CMD_CLOSEST_POINT_IN_RANGE,
    CMD_VIRTUAL_LIST
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"new",CMD_NEW},
    {"push_back",CMD_PUSH_BACK},
    {"push_back_G1",CMD_PUSH_BACK_G1},
    {"reserve",CMD_RESERVE},
    {"getSTK",CMD_GET_STK},
    {"getXY",CMD_GET_XY},
    {"build_G1",CMD_BUILD_G1},
    {"build_3arcG2",CMD_3_ARC_G2},
    {"build_2arcG2",CMD_2_ARC_G2},
    {"build_CLC",CMD_3_ARC_CLC},
    {"build_3arcG2fixed",CMD_3_ARC_G2_FIXED},
    {"build",CMD_BUILD},
    {"build_theta",CMD_BUILD_THETA},
    {"make_closed",CMD_MAKE_CLOSED},
    {"make_open",CMD_MAKE_OPEN},
    {"is_closed",CMD_IS_CLOSED},
    {"get",CMD_GET},
    {"numSegment",CMD_NUM_SEGMENT},
    {"deltaTheta",CMD_DELTA_THETA},
    {"deltaKappa",CMD_DELTA_KAPPA},
    {"export_table",CMD_EXPORT_TABLE},
    {"export_ruby",CMD_EXPORT_RUBY},
    {"findST1",CMD_FIND_ST1},
    {"bbTriangles",CMD_BB_TRIANGLES},
    {"aabb_true",CMD_AABB_TRUE},
    {"aabb_false",CMD_AABB_FALSE},
    {"closestSegment",CMD_CLOSEST_SEGMENT},
    {"closestPointInRange",CMD_CLOSEST_POINT_IN_RANGE},
    CMD_MAP_LIST
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

      switch ( cmd_to_idx.at(cmd) ) {
      case CMD_NEW:
        do_new( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_PUSH_BACK:
        do_push_back( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_PUSH_BACK_G1:
        do_push_back_G1( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_RESERVE:
        do_reserve( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_GET_STK:
        do_getSTK( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_GET_XY:
        do_getXY( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD_G1:
        do_build_G1( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_3_ARC_G2:
        do_build_3arcG2( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_2_ARC_G2:
        do_build_2arcG2( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_3_ARC_CLC:
        do_build_CLC( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_3_ARC_G2_FIXED:
        do_build_3arcG2fixed( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD:
        do_build( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD_THETA:
        do_build_theta( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_MAKE_CLOSED:
        do_make_closed( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_MAKE_OPEN:
        do_make_open( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_IS_CLOSED:
        do_is_closed( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_GET:
        do_get( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_NUM_SEGMENT:
        do_numSegment( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_DELTA_THETA:
        do_deltaTheta( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_DELTA_KAPPA:
        do_deltaKappa( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EXPORT_TABLE:
        do_export_table( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_EXPORT_RUBY:
        do_export_ruby( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_FIND_ST1:
        do_findST1( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BB_TRIANGLES:
        do_bbTriangles( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_AABB_TRUE:
        do_aabb_true( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_AABB_FALSE:
        do_aabb_false( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_CLOSEST_SEGMENT:
        do_closestSegment( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_CLOSEST_POINT_IN_RANGE:
        do_closestPointInRange( nlhs, plhs, nrhs, prhs );
        break;
      CMD_CASE_LIST;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("ClothoidList failed\n");
    }
  }
}
