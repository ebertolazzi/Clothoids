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

#include <vector>

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute polyline\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = PolyLineMexWrapper( 'new', x0, y0, theta0, L );\n" \
"  OBJ = PolyLineMexWrapper( 'new', p0, p1 );\n" \
"\n" \
"  PolyLineMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  PolyLineMexWrapper( 'build', OBJ, x0, y0, theta0, L );\n" \
"  PolyLineMexWrapper( 'build', OBJ, p0, p1 );\n" \
"  PolyLineMexWrapper( 'build', OBJ, p0, theta0, L );\n" \
"  PolyLineMexWrapper( 'copy', OBJ, OBJ1 );\n" \
"  [p1,p2] = PolyLineMexWrapper( 'points', OBJ );\n" \
"\n" \
"  burbs = PolyLineMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("PolyLineMexWrapper") \
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
  PolyLine *
  DATA_NEW( mxArray * & mx_id ) {
    PolyLine * ptr = new PolyLine();
    mx_id = convertPtr2Mat<PolyLine>(ptr);
    return ptr;
  }

  static
  inline
  PolyLine *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<PolyLine>(mx_id);
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
    #define CMD "PolyLineMexWrapper('new'): "
    MEX_ASSERT( nlhs == 1, "expected 1 input, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 1, "expected 1 output, nrhs = " << nrhs );
    DATA_NEW( arg_out_0 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "PolyLineMexWrapper('build',OBJ,x,y): "
    MEX_ASSERT( nlhs == 0, "expected no output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 4, "expected 4 input, nrhs = " << nrhs );

    mwSize size0, size1;

    real_type const * x = getVectorPointer(
      arg_in_2, size0,
      CMD "`x` expected to be a real vector"
    );
    real_type const * y = getVectorPointer(
      arg_in_3, size1,
      CMD "`y` expected to be a real vector"
    );

    MEX_ASSERT(
      size0 == size1,
      CMD "expected size(x) = " << size0 <<
      " = size(y) = " << size1
    );

    PolyLine * ptr = DATA_GET( arg_in_1 );
    ptr->build( x, y, size0 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_polygon(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "PolyLineMexWrapper('polygon',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

    PolyLine * ptr = DATA_GET( arg_in_1 );
    real_type * x = createMatrixValue( arg_out_0, ptr->numPoints(), 1 );
    real_type * y = createMatrixValue( arg_out_1, ptr->numPoints(), 1 );
    ptr->polygon( x, y );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_approx(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "PolyLineMexWrapper('approx',OBJ,OBJ1,tol,type): "
    MEX_ASSERT( nrhs == 5, CMD "expected 5 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 0, CMD "expected NO output, nlhs = " << nlhs );

    real_type tol = getScalarValue(
      arg_in_3, CMD "`tol` expected to be a real scalar"
    );

    MEX_ASSERT( mxIsChar(arg_in_4), CMD "'type' argument must be a string" );
    string kind = mxArrayToString( arg_in_4 );

    PolyLine * ptr = DATA_GET( arg_in_1 );
    if ( kind == "LineSegment" ) {
      ptr->build( *convertMat2Ptr<LineSegment>(arg_in_2) );
    } else if ( kind == "CircleArc" ) {
      ptr->build( *convertMat2Ptr<CircleArc>(arg_in_2), tol );
    } else if ( kind == "BiArc" ) {
      ptr->build( *convertMat2Ptr<Biarc>(arg_in_2), tol );
    } else if ( kind == "ClothoidCurve" ) {
      ptr->build( *convertMat2Ptr<ClothoidCurve>(arg_in_2), tol );
    } else if ( kind == "ClothoidList" ) {
      ptr->build( *convertMat2Ptr<ClothoidList>(arg_in_2), tol );
    } else {
      MEX_ASSERT( false, CMD "'type' = '" << kind << "' unsupported" );
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define CMD_BASE "PolyLineMexWrapper"
  #define G2LIB_CLASS PolyLine
  #include "mex_common.hxx"
  #undef CMD_BASE
  #undef G2LIB_CLASS

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_NEW,
    CMD_BUILD,
    CMD_POLYGON,
    CMD_APPROX,
    CMD_VIRTUAL_LIST
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"new",CMD_NEW},
    {"build",CMD_BUILD},
    {"polygon",CMD_POLYGON},
    {"approx",CMD_APPROX},
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
      case CMD_BUILD:
        do_build( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_POLYGON:
        do_polygon( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_APPROX:
        do_approx( nlhs, plhs, nrhs, prhs );
        break;
      CMD_CASE_LIST;
      }

    } catch ( exception const & e ) {
      mexErrMsgTxt(e.what());
    } catch (...) {
      mexErrMsgTxt("ClothoidCurve failed\n");
    }
  }
}
