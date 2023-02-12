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

    if ( nrhs == 0 && nlhs == 0 ) { mexErrMsgTxt(MEX_ERROR_MESSAGE); return; }

    try {

      if ( nrhs == 2 ) {

        mwSize nx, ny;
        real_type const * x = Utils::mex_vector_pointer( arg_in_0, nx, "XY_to_angle: argument `x` expected to be a real scalar/vector");
        real_type const * y = Utils::mex_vector_pointer( arg_in_1, ny, "XY_to_angle: argument `y` expected to be a real scalar/vector");

        UTILS_MEX_ASSERT(
          nx == ny,
          "XY_to_angle( x, y ), size(x) [{}] != size(y) [{}]\n", nx, ny
        );

        real_type * theta     = Utils::mex_create_matrix_value( arg_out_0, 1, nx );
        real_type * theta_min = Utils::mex_create_matrix_value( arg_out_1, 1, nx );
        real_type * theta_max = Utils::mex_create_matrix_value( arg_out_2, 1, nx );
        real_type * omega     = Utils::mex_create_matrix_value( arg_out_3, 1, nx );
        real_type * len       = Utils::mex_create_matrix_value( arg_out_4, 1, nx );
        xy_to_guess_angle( nx, x, y, theta, theta_min, theta_max, omega, len );

      } else {
        UTILS_MEX_ASSERT0( false, "XY_to_angle: expected 2 arguments" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "XY_to_angle Error: {}", e.what() ).c_str() );
    } catch (...) {
  	  mexErrMsgTxt("XY_to_angle failed\n");
    }

  }

}
