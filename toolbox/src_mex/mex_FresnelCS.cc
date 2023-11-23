/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoids.hh"
#include "Utils_mex.hh"

#ifdef _MSC_VER
  #pragma comment(lib, "IPHLPAPI.lib")
  #pragma comment(lib, "ws2_32.lib")
  #pragma comment(lib, "Shlwapi.lib")
  #pragma comment(lib, "Advapi32.lib")
  #pragma comment(lib, "Shell32.lib")
  #pragma comment(lib, "kernel32.lib")
#endif

#define MEX_ERROR_MESSAGE \
"======================================================================\n" \
"\n" \
"FresnelCS: Compute Fresnel sine and cosine integrals\n" \
"           or Fresnel sine and cosine integrals momenta\n" \
"\n" \
"USAGE:\n" \
"  [C,S] = FresnelCS(y);\n" \
"  [X,Y] = FresnelCS( nk, a, b, c );\n" \
"\n" \
"Fresnel integral are defined as:\n" \
"\n" \
"  C(y) = int_0^y cos( (pi/2)*t^2 ) dt\n" \
"  S(y) = int_0^y sin( (pi/2)*t^2 ) dt\n" \
"\n" \
"  X_k(a,b,c) = int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt\n" \
"  Y_k(a,b,c) = int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt\n" \
"\n" \
"The algorithm for the computation of C(y) and S(y) is described in:\n" \
"  Atlas for computing mathematical functions: an illustrated guide\n" \
"  for practitioners, with programs in C and Mathematica.\n" \
"  William J. Thompson New York : Wiley, 1997.\n" \
"\n" \
"The code is a sligly modification and translation to C language\n" \
"of original code developed by\n" \
"  Venkata Sivakanth Telasula (sivakanth.telasula@gmail.com)\n" \
"\n" \
"On input:\n" \
"\n" \
"  y = argument of the function C(y) and S(y), it may be a vector\n" \
"\n" \
"On output:\n" \
"\n" \
"  C = The value(s) of C(y)\n" \
"  S = The value(s) of S(y)\n" \
"  X = The value(s) of X_k(a,b,c)\n" \
"  Y = The value(s) of Y_k(a,b,c)\n" \
"\n" \
"======================================================================\n" \
"\n" \
"Autor: Enrico Bertolazzi\n" \
"       Department of Industrial Engineering\n" \
"       University of Trento\n" \
"       enrico.bertolazzi@unitn.it\n" \
"\n" \
"======================================================================\n"

#include <unordered_map>

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

      if ( nrhs == 1 ) {

        mwSize size0;
        double const * y = Utils::mex_vector_pointer(
          arg_in_0, size0,
          "FresnelCS: argument `y` expected to be a real scalar/vector"
        );
        double * C = Utils::mex_create_matrix_value( arg_out_0, 1, size0 );
        double * S = Utils::mex_create_matrix_value( arg_out_1, 1, size0 );
        for ( mwSize i = 0; i < size0; ++i ) FresnelCS( *y++, *C++, *S++ );

      } else if ( nrhs == 4 ) {

        integer nk = Utils::mex_get_int64(
          arg_in_0, "FresnelCS: argument `nk` expected to be and integer"
        );
        UTILS_MEX_ASSERT(
          nk >= 1 && nk <= 3,
          "FresnelCS: argument `nk` = {} must be in [1,2,3]\n", nk
        );
        mwSize na, nb, nc;
        double const * a = Utils::mex_vector_pointer(
          arg_in_1, na, "FresnelCS: argument `a` expected to be a real scalar/vector"
        );
        double const * b = Utils::mex_vector_pointer(
          arg_in_2, nb, "FresnelCS: argument `b` expected to be a real scalar/vector"
        );
        double const * c = Utils::mex_vector_pointer(
          arg_in_3, nc, "FresnelCS: argument `c` expected to be a real scalar/vector"
        );
        UTILS_MEX_ASSERT(
          na == nb && nb == nc,
          "FresnelCS: Second to last arguments must be vectors of the same length\n"
          "found length(a) = {}, length(b) = {}, length(c) = {}\n",
          na, nb, nc
        );
        double * X = Utils::mex_create_matrix_value( arg_out_0, nk, na );
        double * Y = Utils::mex_create_matrix_value( arg_out_1, nk, na );
        for ( mwSize k = 0; k < na; ++k ) {
          GeneralizedFresnelCS( nk, *a, *b, *c, X, Y );
          ++a; ++b; ++c; X += nk; Y += nk;
        }
      } else {
        UTILS_MEX_ASSERT0( false, "FresnelCS: expected 1 or 4 airguments" );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt( fmt::format( "FresnelCS Error: {}", e.what() ).c_str() );
    } catch (...) {
  	  mexErrMsgTxt("FresnelCS failed\n");
    }

  }

}
