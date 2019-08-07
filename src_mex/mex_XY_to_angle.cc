/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "G2lib.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"======================================================================\n" \
"\n" \
"XY_to_angle:\n" \
"\n" \
"USAGE:\n" \
"  [theta,theta_min,theta_max,omega,d] = XY_to_angle(x,y);\n" \
"\n" \
"\n" \
"======================================================================\n" \
"\n" \
"Autor: Enrico Bertolazzi\n" \
"       Department of Industrial Engineering\n" \
"       University of Trento\n" \
"       enrico.bertolazzi@unitn.it\n" \
"\n" \
"======================================================================\n"

namespace G2lib {

  extern "C"
  void
  mexFunction(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    if ( nrhs == 0 && nlhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {

      if ( nrhs == 2 ) {

        mwSize nx, ny;
        real_type const * x = getVectorPointer( arg_in_0, nx, "XY_to_angle: argument `x` expected to be a real scalar/vector");
        real_type const * y = getVectorPointer( arg_in_1, ny, "XY_to_angle: argument `y` expected to be a real scalar/vector");

        MEX_ASSERT(
          nx == ny,
          "XY_to_angle( x, y ), size(x) [" << nx << "] != size(y) [" << ny << "]"
        );

        real_type * theta     = createMatrixValue( arg_out_0, 1, nx );
        real_type * theta_min = createMatrixValue( arg_out_1, 1, nx );
        real_type * theta_max = createMatrixValue( arg_out_2, 1, nx );
        real_type * omega     = createMatrixValue( arg_out_3, 1, nx );
        real_type * len       = createMatrixValue( arg_out_4, 1, nx );
        xy_to_guess_angle( nx, x, y, theta, theta_min, theta_max, omega, len );

      } else {
        MEX_ASSERT( false, "XY_to_angle: expected 2 arguments" );
      }

    } catch ( std::exception const & e ) {
    	mexErrMsgTxt(e.what());

    } catch (...) {
  	  mexErrMsgTxt("XY_to_angle failed\n");
    }

  }

}
