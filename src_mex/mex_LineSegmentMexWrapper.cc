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
"  OBJ = LineSegmentMexWrapper( 'new', x0, y0, theta0, L );\n" \
"  OBJ = LineSegmentMexWrapper( 'new', p0, p1 );\n" \
"\n" \
"  LineSegmentMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  LineSegmentMexWrapper( 'build', OBJ, x0, y0, theta0, L );\n" \
"  LineSegmentMexWrapper( 'build', OBJ, p0, p1 );\n" \
"  LineSegmentMexWrapper( 'build', OBJ, p0, theta0, L );\n" \
"  [p1,p2] = LineSegmentMexWrapper( 'points', OBJ );\n" \
"\n" \
"  nurbs = LineSegmentMexWrapper( 'to_nurbs', OBJ );\n" \
"\n" \
MEX_INFO_MESSAGE("LineSegmentMexWrapper") \
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
  LineSegment *
  DATA_NEW( mxArray * & mx_id ) {
    LineSegment * ptr = new LineSegment();
    mx_id = convertPtr2Mat<LineSegment>(ptr);
    return ptr;
  }

  static
  inline
  LineSegment *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<LineSegment>(mx_id);
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

    #define CMD "LineSegmentMexWrapper('new',...): "
    MEX_ASSERT(
      nlhs == 1,
      CMD "expected 1 output, nlhs = " << nlhs
    );
    MEX_ASSERT(
      nrhs == 1 || nrhs == 3 || nrhs == 5,
      CMD "expected 1, 3 or 5 inputs, nrhs = " << nrhs
    );
    #undef CMD

    LineSegment * ptr = DATA_NEW( arg_out_0 );

    if ( nrhs == 5 ) {

      real_type x0, y0, th0, L;
      #define CMD "LineSegmentMexWrapper('new',x0,y0,theta0,L): "

      x0 = getScalarValue(
        arg_in_1, CMD "`x0` expected to be a real scalar"
      );
      y0 = getScalarValue(
        arg_in_2, CMD "`y0` expected to be a real scalar"
      );
      th0 = getScalarValue(
        arg_in_3, CMD "`theta0` expected to be a real scalar"
      );
      L = getScalarValue(
        arg_in_4, CMD "`L` expected to be a real scalar"
      );
      ptr->build( x0, y0, th0, L );

      #undef CMD

    } else if ( nrhs == 3 ) {

      real_type const * p0;
      real_type const * p1;
      #define CMD "LineSegmentMexWrapper('new',OBJ,p0,p1): "

      mwSize size0, size1;
      p0 = getVectorPointer(
        arg_in_1, size0,
        CMD "`p0` expected to be a real vector"
      );
      p1 = getVectorPointer(
        arg_in_2, size1,
        CMD "`p1` expected to be a real vector"
      );

      MEX_ASSERT(
        size0 == 2 && size1 == 2,
        CMD "bad dimension size(p0) = " << size0 <<
        ", size(p1) = " << size1
      );
      #undef CMD

      ptr->build_2P( p0[0], p0[1], p1[0], p1[1] );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "LineSegmentMexWrapper('build',OBJ,...): "
    MEX_ASSERT(
      nlhs == 0,
      CMD "expected NO output, nlhs = " << nlhs
    );
    MEX_ASSERT(
      nrhs == 4 || nrhs == 6,
      CMD "expected 4 or 6 inputs, nrhs = " << nrhs
    );
    #undef CMD

    LineSegment * ptr = DATA_GET(arg_in_1);

    if ( nrhs == 6 ) {
      real_type x0, y0, th0, L;

      #define CMD "LineSegmentMexWrapper('build',OBJ,x0,y0,theta0,L): "
      x0 = getScalarValue(
        arg_in_2, CMD "`x0` expected to be a real scalar"
      );
      y0 = getScalarValue(
        arg_in_3, CMD "`y0` expected to be a real scalar"
      );
      th0 = getScalarValue(
        arg_in_4, CMD "`theta0` expected to be a real scalar"
      );
      L = getScalarValue(
        arg_in_5, CMD "`L` expected to be a real scalar"
      );

      ptr->build( x0, y0, th0, L );

      #undef CMD

    } else if ( nrhs == 4 ) {
      real_type const * p0;
      real_type const * p1;

      #define CMD "LineSegmentMexWrapper('build',OBJ,p0,p1): "

      mwSize size0, size1;
      p0 = getVectorPointer(
        arg_in_2, size0,
        CMD "`p0` expected to be a real vector"
      );
      p1 = getVectorPointer(
        arg_in_3, size1,
        CMD "`p1` expected to be a real vector"
      );

      MEX_ASSERT(
        size0 == 2 && size1 == 2,
        CMD "bad dimension size(p0) = " << size0 <<
        ", size(p1) = " << size1
      );
      #undef CMD

      ptr->build_2P( p0[0], p0[1], p1[0], p1[1] );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_to_nurbs(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    LineSegment * ptr = DATA_GET(arg_in_1);

    #define CMD "LineSegmentMexWrapper('to_nurbs',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 inputs, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );

    int_type npts, nknots;
    ptr->paramNURBS( nknots, npts );

    real_type knots[12], Poly[9][3];

    ptr->toNURBS( knots, Poly ); // npt + 2

    static char const * fieldnames[] = {
      "form",
      "order",
      "dim",
      "number",
      "knots",
      "coefs"
    };

    arg_out_0 = mxCreateStructMatrix(1,1,6,fieldnames);
    mxArray * mx_knots = mxCreateDoubleMatrix(1,nknots,mxREAL);
    mxArray * mx_Poly  = mxCreateDoubleMatrix(3,npts,mxREAL);

    mxSetFieldByNumber( arg_out_0, 0, 0, mxCreateString("rB") );
    mxSetFieldByNumber( arg_out_0, 0, 1, mxCreateDoubleScalar(2) );
    mxSetFieldByNumber( arg_out_0, 0, 2, mxCreateDoubleScalar(2) );
    mxSetFieldByNumber( arg_out_0, 0, 3, mxCreateDoubleScalar(npts) );
    mxSetFieldByNumber( arg_out_0, 0, 4, mx_knots );
    mxSetFieldByNumber( arg_out_0, 0, 5, mx_Poly );

    double *kb = mxGetPr(mx_knots);
    for ( int_type i = 0; i < nknots; ++i ) *kb++ = knots[i];

    double *pr = mxGetPr(mx_Poly);
    for ( int_type i = 0; i < npts; ++i ) {
      *pr++ = Poly[i][0];
      *pr++ = Poly[i][1];
      *pr++ = Poly[i][2];
    }

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_points(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    LineSegment * ptr = DATA_GET(arg_in_1);

    #define CMD "LineSegmentMexWrapper('points',OBJ): "
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
    MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );

    real_type * p1 = createMatrixValue( arg_out_0, 2, 1 );
    real_type * p2 = createMatrixValue( arg_out_1, 2, 1 );

    ptr->p1p2(p1,p2);

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  #define CMD_BASE "LineSegmentMexWrapper"
  #define G2LIB_CLASS LineSegment
  #include "mex_common.hxx"
  #undef CMD_BASE
  #undef G2LIB_CLASS

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef enum {
    CMD_NEW,
    CMD_BUILD,
    CMD_TO_NURBS,
    CMD_POINTS,
    CMD_VIRTUAL_LIST
  } CMD_LIST;

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static map<string,unsigned> cmd_to_idx = {
    {"new",CMD_NEW},
    {"build",CMD_BUILD},
    {"to_nurbs",CMD_TO_NURBS},
    {"points",CMD_POINTS},
    CMD_MAP_LIST
  };

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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
      case CMD_TO_NURBS:
        do_to_nurbs( nlhs, plhs, nrhs, prhs );
        break;
      case CMD_POINTS:
        do_points( nlhs, plhs, nrhs, prhs );
        break;
      CMD_CASE_LIST;
      }

    } catch ( exception const & e ) {
    	mexErrMsgTxt(e.what());
    } catch (...) {
    	mexErrMsgTxt("Line failed\n");
    }

  }

}
