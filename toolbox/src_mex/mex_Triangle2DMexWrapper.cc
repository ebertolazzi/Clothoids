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

#include <vector>

#include "mex_Workaround.hxx"

#define MEX_ERROR_MESSAGE \
"==========================================================================\n" \
"Compute Triangle2D.hh\n" \
"\n" \
"USAGE:\n" \
"\n" \
"  OBJ = Triangle2DMexWrapper( 'new' );\n" \
"  OBJ = Triangle2DMexWrapper( 'new', x0, y0, x1, y1, x2, y2 );\n" \
"  OBJ = Triangle2DMexWrapper( 'new', p0, p1, p2 );\n" \
"\n" \
"  Triangle2DMexWrapper( 'delete', OBJ );\n" \
"\n" \
"  Triangle2DMexWrapper( 'build', OBJ, x0, y0, x1, y1, x2, y2 );\n" \
"  Triangle2DMexWrapper( 'build', OBJ, p0, p1, p2 );\n" \
"  [p1,p2,p3] = Triangle2DMexWrapper( 'points', OBJ );\n" \
"\n" \
"  Triangle2DMexWrapper( 'translate', OBJ, tx, ty );\n" \
"  Triangle2DMexWrapper( 'rotate', OBJ, angle, cx, cy );\n" \
"  Triangle2DMexWrapper( 'scale', OBJ, scale );\n" \
"\n" \
"  [dmin,dmax] = Triangle2DMexWrapper( 'distance', OBJ, x, y );\n" \
"  [icode] = Triangle2DMexWrapper( 'isInside', OBJ, x, y );\n" \
"\n" \
MEX_INFO_MESSAGE_END

namespace G2lib {

  using namespace std;

  static
  Triangle2D *
  DATA_NEW( mxArray * & mx_id ) {
    Triangle2D * ptr = new Triangle2D();
    mx_id = convertPtr2Mat<Triangle2D>(ptr);
    return ptr;
  }

  static
  inline
  Triangle2D *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<Triangle2D>(mx_id);
  }

  static
  void
  DATA_DELETE( mxArray const * & mx_id ) {
    destroyObject<Triangle2D>(mx_id);
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    MEX_ASSERT2( nrhs >= 2, "Triangle2DMexWrapper, expected at least 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, "Triangle2DMexWrapper, expected 1 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_NEW(arg_out_0);

    if ( nrhs == 7 ) {
      #define CMD "Triangle2DMexWrapper('new',x0,y0,x1,y1,x2,y2): "
      real_type x0 = getScalarValue( arg_in_1, CMD "`x0` expected to be a real scalar" );
      real_type y0 = getScalarValue( arg_in_2, CMD "`y0` expected to be a real scalar" );
      real_type x1 = getScalarValue( arg_in_3, CMD "`x1` expected to be a real scalar" );
      real_type y1 = getScalarValue( arg_in_4, CMD "`y1` expected to be a real scalar" );
      real_type x2 = getScalarValue( arg_in_5, CMD "`x2` expected to be a real scalar" );
      real_type y2 = getScalarValue( arg_in_6, CMD "`y2` expected to be a real scalar" );
      ptr->build( x0, y0, x1, y1, x2, y2, 0, 0, 0 );
      #undef CMD

   } else if ( nrhs == 4 ) {

      #define CMD "Triangle2DMexWrapper('new',p0,p1,p2): "
      mwSize sz0, sz1, sz2;
      real_type const * p0 = getVectorPointer( arg_in_1, sz0, CMD "`p0` expected to be a real vector" );
      real_type const * p1 = getVectorPointer( arg_in_2, sz1, CMD "`p1` expected to be a real vector" );
      real_type const * p2 = getVectorPointer( arg_in_3, sz2, CMD "`p2` expected to be a real vector" );
      MEX_ASSERT( sz0 == 2, "Triangle2D, expected a vector of 2 elements for `p0`" );
      MEX_ASSERT( sz1 == 2, "Triangle2D, expected a vector of 2 elements for `p1`" );
      MEX_ASSERT( sz2 == 2, "Triangle2D, expected a vector of 2 elements for `p2`" );
      ptr->build( p0, p1, p2, 0, 0, 0  );
      #undef CMD

    } else {
     MEX_ASSERT2( false, "Triangle2D, expected 4 or 7 inputs, nrhs = {}\n", nrhs );
    }

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_build(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    MEX_ASSERT2( nrhs >= 2, "Triangle2DMexWrapper, expected at least 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, "Triangle2DMexWrapper, expected no output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    if ( nrhs == 8 ) {

      #define CMD "Triangle2DMexWrapper('build',OBJ,x0,y0,x1,y1,x2,y2): "
      real_type x0 = getScalarValue( arg_in_2, CMD "`x0` expected to be a real scalar" );
      real_type y0 = getScalarValue( arg_in_3, CMD "`y0` expected to be a real scalar" );
      real_type x1 = getScalarValue( arg_in_4, CMD "`x1` expected to be a real scalar" );
      real_type y1 = getScalarValue( arg_in_5, CMD "`y1` expected to be a real scalar" );
      real_type x2 = getScalarValue( arg_in_6, CMD "`x2` expected to be a real scalar" );
      real_type y2 = getScalarValue( arg_in_7, CMD "`y2` expected to be a real scalar" );
      ptr->build( x0, y0, x1, y1, x2, y2, 0, 0, 0  );
      #undef CMD

    } else if ( nrhs == 5 ) {

      #define CMD "Triangle2DMexWrapper('build',OBJ,p0,p1,p2): "
      mwSize sz0, sz1, sz2;
      real_type const * p0 = getVectorPointer( arg_in_2, sz0, CMD "`p0` expected to be a real vector" );
      real_type const * p1 = getVectorPointer( arg_in_3, sz1, CMD "`p1` expected to be a real vector" );
      real_type const * p2 = getVectorPointer( arg_in_4, sz2, CMD "`p2` expected to be a real vector" );
      MEX_ASSERT( sz0 == 2, "Triangle2D, expected a vector of 2 elements for `p0`" );
      MEX_ASSERT( sz1 == 2, "Triangle2D, expected a vector of 2 elements for `p1`" );
      MEX_ASSERT( sz2 == 2, "Triangle2D, expected a vector of 2 elements for `p2`" );
      ptr->build( p0, p1, p2, 0, 0, 0  );
      #undef CMD
    } else {
      MEX_ASSERT2(false, "Triangle2D, expected 5 or 8 inputs, nrhs = {}\n", nrhs );
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('delete',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );
    // Destroy the C++ object
    DATA_DELETE( arg_in_1 );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_translate(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "Triangle2DMexWrapper('translate',OBJ,tx,ty): "
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    real_type tx = getScalarValue( arg_in_2, CMD "`tx` expected to be a real scalar" );
    real_type ty = getScalarValue( arg_in_3, CMD "`ty` expected to be a real scalar" );

    ptr->translate( tx, ty );
    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_rotate(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('rotate',OBJ,angle,cx,cy): "

    MEX_ASSERT2( nrhs == 5, CMD "expected 5 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    real_type angle = getScalarValue( arg_in_2, CMD "`angle` expected to be a real scalar" );
    real_type cx    = getScalarValue( arg_in_3, CMD "`cx` expected to be a real scalar" );
    real_type cy    = getScalarValue( arg_in_4, CMD "`cy` expected to be a real scalar" );

    ptr->rotate( angle, cx, cy );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_scale(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('scale',OBJ,scale): "
    MEX_ASSERT2( nrhs == 3, CMD "expected 3 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected no output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    real_type sc = getScalarValue( arg_in_2, CMD "`scale` expected to be a real scalar" );
    ptr->scale( sc );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_distanceMin(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('distanceMin',OBJ,x,y): "
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    mwSize nrx, ncx, nry, ncy;
    real_type const * x = getMatrixPointer(
      arg_in_2, nrx, ncx,
      CMD "`x` expected to be a real vector/matrix"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, nry, ncy,
      CMD "`y` expected to be a real vectormatrix"
    );
    MEX_ASSERT2(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found\n"
      "size(x) = {} x {} size(y) = {} x {}\n",
      nrx, ncx, nry, ncy
    );

    real_type * dst = createMatrixValue( arg_out_0, nrx, ncx );

    mwSize size = nrx*ncx;
    for ( mwSize i = 0; i < size; ++i )
      *dst++ = ptr->distMin( *x++, *y++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_distanceMax(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('distanceMax',OBJ,x,y): "
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    mwSize nrx, ncx, nry, ncy;
    real_type const * x = getMatrixPointer(
      arg_in_2, nrx, ncx,
      CMD "`x` expected to be a real vector/matrix"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, nry, ncy,
      CMD "`y` expected to be a real vectormatrix"
    );
    MEX_ASSERT2(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found\n"
      "size(x) = {} x {} size(y) = {} x {}\n",
      nrx, ncx, nry, ncy
    );

    real_type * dst = createMatrixValue( arg_out_0, nrx, ncx );

    mwSize size = nrx*ncx;
    for ( mwSize i = 0; i < size; ++i )
      *dst++ = ptr->distMax( *x++, *y++ );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_isInside(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('isInside',OBJ,x,y): "
    MEX_ASSERT2( nrhs == 4, CMD "expected 4 input, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    mwSize nrx, ncx, nry, ncy;
    real_type const * x = getMatrixPointer(
      arg_in_2, nrx, ncx,
      CMD "`x` expected to be a real vector/matrix"
    );
    real_type const * y = getMatrixPointer(
      arg_in_3, nry, ncy,
      CMD "`y` expected to be a real vectormatrix"
    );
    MEX_ASSERT2(
      nrx == nry && ncx == ncy,
      CMD "`x` and `y` expected to be of the same size, found\n"
      "size(x) = {} x {} size(y) = {} x {}",
      nrx, ncx, nry, ncy
    );

    int64_t * icode = createMatrixInt64( arg_out_0, nrx, ncx );

    mwSize size = nrx*ncx;
    for ( mwSize i = 0; i < size; ++i )
      *icode++ = int64_t(ptr->isInside( *x++, *y++ ));

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_overlap(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('overlap',OBJ,OBJ2): "
    MEX_ASSERT2( nrhs == 3, CMD "expected 3 input, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 1, CMD "expected 1 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr  = DATA_GET(arg_in_1);
    Triangle2D * ptr2 = DATA_GET(arg_in_2);

    setScalarBool( arg_out_0, ptr->overlap( * ptr2 ) ) ;

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_points(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('points',OBJ): "
    MEX_ASSERT2( nrhs == 2, CMD "expected 2 input, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 3, CMD "expected 3 output, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    real_type * p1 = createMatrixValue( arg_out_0, 2, 1 );
    real_type * p2 = createMatrixValue( arg_out_1, 2, 1 );
    real_type * p3 = createMatrixValue( arg_out_2, 2, 1 );

    std::copy( ptr->P1(), ptr->P1()+2, p1 );
    std::copy( ptr->P2(), ptr->P2()+2, p2 );
    std::copy( ptr->P3(), ptr->P3()+2, p3 );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_info(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "Triangle2DMexWrapper('info',OBJ): "

    MEX_ASSERT2( nrhs == 2, CMD "expected 2 inputs, nrhs = {}\n", nrhs );
    MEX_ASSERT2( nlhs == 0, CMD "expected NO outputs, nlhs = {}\n", nlhs );

    Triangle2D * ptr = DATA_GET(arg_in_1);

    ptr->info( std::cout );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"build",do_build},
    {"delete",do_delete},
    {"translate", do_translate},
    {"rotate", do_rotate},
    {"scale", do_scale},
    {"distanceMin", do_distanceMin},
    {"distanceMax", do_distanceMax},
    {"isInside", do_isInside},
    {"overlap", do_overlap},
    {"points", do_points},
    {"info", do_info}
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
      MEX_ASSERT( mxIsChar(arg_in_0), "Triangle2DMexWrapper(...): First argument must be a string" );
      string cmd = mxArrayToString(arg_in_0);
      DO_CMD pfun = cmd_to_fun.at(cmd);
      pfun( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "Triangle2D Error: {}", e.what() ).c_str() );
    } catch (...) {
    	mexErrMsgTxt("Triangle2D failed\n");
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

}
