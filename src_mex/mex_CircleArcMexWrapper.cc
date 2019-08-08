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
#include "Clothoid.hh"
#include "ClothoidList.hh"

#include "mex_utils.hh"
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
  CircleArc *
  DATA_NEW( mxArray * & mx_id ) {
    CircleArc * ptr = new CircleArc();
    mx_id = convertPtr2Mat<CircleArc>(ptr);
    return ptr;
  }

  static
  inline
  CircleArc *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<CircleArc>(mx_id);
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

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('new',...): "
    MEX_ASSERT(
      nrhs == 1 || nrhs == 6,
      CMD "expected 1 or 6 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 1,
      CMD "expected 1 output, nlhs = " << nlhs
    );
    #undef CMD

    CircleArc * ptr = DATA_NEW( arg_out_0 );

    if ( nrhs == 6 ) {
      #define CMD "CircleArcMexWrapper('new',x0,y0,theta0,k0,L): "
      real_type x0, y0, theta0, k0, L;
      x0     = getScalarValue( arg_in_1, CMD "`x0` expected to be a real scalar" );
      y0     = getScalarValue( arg_in_2, CMD "`y0` expected to be a real scalar" );
      theta0 = getScalarValue( arg_in_3, CMD "`theta0` expected to be a real scalar" );
      k0     = getScalarValue( arg_in_4, CMD "`k0` expected to be a real scalar" );
      L      = getScalarValue( arg_in_5, CMD "`L` expected to be a real scalar" );

      ptr->build( x0, y0, theta0, k0, L );
      #undef CMD
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('build',OBJ,x0,y0,theta0,k0,L): "
    MEX_ASSERT(
      nrhs == 7, CMD "expected 7 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 0, CMD "expected NO output, nlhs = " << nlhs
    );

    CircleArc * ptr = DATA_GET( arg_in_1 );

    real_type x0, y0, theta0, k0, L;
    x0     = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
    y0     = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );
    theta0 = getScalarValue( arg_in_4, CMD "`theta0` expected to be a real scalar" );
    k0     = getScalarValue( arg_in_5, CMD "`k0` expected to be a real scalar" );
    L      = getScalarValue( arg_in_6, CMD "`L` expected to be a real scalar" );

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
    MEX_ASSERT(
      nlhs == 0 || nlhs ==1,
      CMD " expected 1 or no output, nlhs = " << nlhs
    );
    #undef CMD

    CircleArc * ptr = DATA_GET(arg_in_1);

    real_type x0(0), y0(0), x1(0), y1(0), x2(0), y2(0);
    if ( nrhs == 5 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,p0,p1,p2): "
      real_type const * p0;
      real_type const * p1;
      real_type const * p2;
      mwSize size0, size1, size2;

      p0 = getVectorPointer(
        arg_in_2, size0, CMD "`p0` expected to be a real vector"
      );
      p1 = getVectorPointer(
        arg_in_3, size1, CMD "`p1` expected to be a real vector"
      );
      p2 = getVectorPointer(
        arg_in_4, size2, CMD "`p2` expected to be a real vector"
      );

      MEX_ASSERT(
        size0 == 2 && size1 == 2 && size2 == 2,
        CMD "bad dimension size(p0) = " << size0 <<
        ", size(p1) = " << size1 << ", size(p2) = " << size2
      );
      #undef CMD

      x0 = p0[0]; y0 = p0[1];
      x1 = p1[0]; y1 = p1[1];
      x2 = p2[0]; y2 = p2[1];

    } else if ( nrhs == 8 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,x0,x1,x1,y1,x2,y2): "
      x0 = getScalarValue( arg_in_2, CMD "`x0` expected to be a scalar value" );
      y0 = getScalarValue( arg_in_3, CMD "`y0` expected to be a scalar value" );
      x1 = getScalarValue( arg_in_4, CMD "`x1` expected to be a scalar value" );
      y1 = getScalarValue( arg_in_5, CMD "`y1` expected to be a scalar value" );
      x2 = getScalarValue( arg_in_6, CMD "`x2` expected to be a scalar value" );
      y2 = getScalarValue( arg_in_7, CMD "`y2` expected to be a scalar value" );
      #undef CMD
    } else {
      MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = " << nrhs );
    }

    bool ok = ptr->build_3P( x0, y0, x1, y1, x2, y2 );
    if ( nlhs == 1 ) setScalarBool( arg_out_0, ok );

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
    MEX_ASSERT(
      nlhs == 0 || nlhs ==1,
      CMD "expected 1 or no output, nlhs = " << nlhs
    );
    #undef CMD

    CircleArc * ptr = DATA_GET(arg_in_1);

    real_type x0(0), y0(0), x1(0), y1(0), theta0(0);
    if ( nrhs == 5 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,p0,theta0,p1): "
      real_type const * p0;
      real_type const * p1;
      mwSize size0, size1;
      p0 = getVectorPointer(
        arg_in_2, size0, CMD "`p0` expected to be a real vector"
      );
      p1 = getVectorPointer(
        arg_in_4, size1, CMD "`p1` expected to be a real vector"
      );
      theta0 = getScalarValue(
        arg_in_3, CMD "`theta0` expected to be a real vector"
      );

      MEX_ASSERT(
        size0 == 2 && size1 == 2,
        CMD "bad dimension size(p0) = " << size0 <<
        ", size(p1) = " << size1
      );
      #undef CMD

      x0 = p0[0]; y0 = p0[1];
      x1 = p1[0]; y1 = p1[1];

    } else if ( nrhs == 7 ) {

      #define CMD "CircleArcMexWrapper('build_G1',OBJ,x0,x1,theta0,x1,y1): "
      x0     = getScalarValue( arg_in_2, CMD "`x0` expected to be a scalar value" );
      y0     = getScalarValue( arg_in_3, CMD "`y0` expected to be a scalar value" );
      theta0 = getScalarValue( arg_in_4, CMD "`theta0` expected to be a scalar value" );
      x1     = getScalarValue( arg_in_5, CMD "`x1` expected to be a scalar value" );
      y1     = getScalarValue( arg_in_6, CMD "`y1` expected to be a scalar value" );
      #undef CMD
    } else {
      MEX_ASSERT(false, "CircleArc, expected 5 or 7 inputs, nrhs = " << nrhs );
    }

    bool ok = ptr->build_G1( x0, y0, theta0, x1, y1 );
    if ( nlhs == 1 ) setScalarBool( arg_out_0, ok );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_change_curvilinear_origin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "CircleArcMexWrapper('changeCurvilinearOrigin',OBJ,s0,L): "
    MEX_ASSERT(nrhs == 4, CMD "expected 4 inputs, nrhs = " << nrhs );

    CircleArc * ptr = DATA_GET(arg_in_1);

    real_type s0 = getScalarValue(arg_in_2,CMD "Error in reading s0");
    real_type L  = getScalarValue(arg_in_3,CMD "Error in reading L");
    ptr->changeCurvilinearOrigin(s0,L);
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_to_nurbs(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    CircleArc * ptr = DATA_GET(arg_in_1);

    #define CMD "CircleArcMexWrapper('to_nurbs',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    int_type npts, nknots;
    ptr->paramNURBS( nknots, npts );

    mxArray * mx_knots, * mx_Poly;
    double * knots = createMatrixValue( mx_knots, 1, nknots );
    double * poly  = createMatrixValue( mx_Poly,  3, npts );

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

  #define CMD_BASE "CircleArcMexWrapper"
  #define G2LIB_CLASS CircleArc
  #include "mex_common.hxx"
  #undef CMD_BASE
  #undef G2LIB_CLASS

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_bbTriangles(
    int nlhs, mxArray       * plhs[],
    int nrhs, mxArray const * prhs[]
  ) {

    CircleArc * ptr = DATA_GET(arg_in_1);

    #define CMD "CircleArcMexWrapper('bbTriangles',OBJ[,max_angle,max_size,offs,['ISO'/'SAE']): "

    MEX_ASSERT(
      nrhs >= 2 && nrhs <= 6,
      CMD "expected 2 up to 6 inputs, nrhs = " << nrhs
    );
    MEX_ASSERT(
      nlhs == 3,
      CMD "expected 3 output, nlhs = " << nlhs
    );

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

    std::vector<Triangle2D> tvec;
    if ( nrhs >= 5 ) {
      bool ISO = true;
      if ( nrhs == 6 ) ISO = do_is_ISO( arg_in_7, CMD " last argument must be a string");
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

  typedef enum {
    CMD_NEW,
    CMD_BUILD,
    CMD_BUILD_3P,
    CMD_BUILD_G1,
    CMD_CHANGE_CURVILINEAR_ORIGIN,
    CMD_TO_NURBS,
    CMD_BB_TRIANGLES,
    CMD_VIRTUAL_LIST
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"new",CMD_NEW},
    {"build",CMD_BUILD},
    {"build_3P",CMD_BUILD_3P},
    {"build_G1",CMD_BUILD_G1},
    {"changeCurvilinearOrigin",CMD_CHANGE_CURVILINEAR_ORIGIN},
    {"to_nurbs",CMD_TO_NURBS},
    {"bbTriangles",CMD_BB_TRIANGLES},
    CMD_MAP_LIST
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

      switch ( cmd_to_idx.at(cmd) ) {
      case CMD_NEW:
        do_new( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD:
        do_build( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD_3P:
        do_build_3P( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BUILD_G1:
        do_build_G1( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_CHANGE_CURVILINEAR_ORIGIN:
        do_change_curvilinear_origin( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_TO_NURBS:
        do_to_nurbs( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_BB_TRIANGLES:
        do_bbTriangles( nlhs, plhs, nrhs, prhs );
        break;
      CMD_CASE_LIST;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("CircleArc failed\n");
    }

  }
}
