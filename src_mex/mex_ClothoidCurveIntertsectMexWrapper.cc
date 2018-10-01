/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex_utils.hh"

#include <vector>

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  intersectClothoid:  Compute intersections betweed clothoids         %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    s1, s2 = ClothoidCurveIntertsectMexWrapper( clot1, clot2 );       %\n" \
"%    s1, s2 = ClothoidCurveIntertsectMexWrapper( clot1, clot2,         %\n" \
"%                                                offs1, offs2 );       %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%    clot1 = object pointer of the first clothoid                      %\n" \
"%    clot2 = object pointer of the second clothoid                     %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%   s1 = curvilinear coordinates of intersections on clot1             %\n" \
"%   s2 = curvilinear coordinates of intersections on clot2             %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

#define arg_clot1 prhs[0]
#define arg_clot2 prhs[1]
#define arg_offs1 prhs[2]
#define arg_offs2 prhs[3]

#define arg_S1    plhs[0]
#define arg_S2    plhs[1]

namespace G2lib {

  using namespace std;

  static
  inline
  ClothoidCurve *
  DATA_GET( mxArray const * & mx_id ) {
    return convertMat2Ptr<ClothoidCurve>(mx_id);
  }

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

      MEX_ASSERT( nrhs == 2 || nrhs == 4, "ClothoidCurveMexWrapper, expected 2 or 4 input, nrhs = " << nrhs  );

      ClothoidCurve * c1 = DATA_GET(arg_clot1);
      ClothoidCurve * c2 = DATA_GET(arg_clot2);

      std::vector<real_type> s1, s2;
      int_type  max_iter  = 10;
      real_type tolerance = 1e-8;

      try {
        c1->intersect( *c2, s1, s2, max_iter, tolerance );
      }
      catch (...) {
        mexErrMsgTxt("Intersection failed\n");
      }

      if ( nlhs > 0 ) {
        arg_S1 = mxCreateNumericMatrix( s1.size(),1, mxDOUBLE_CLASS, mxREAL );
        double * pS1 = mxGetPr(arg_S1);
        for ( unsigned i = 0; i < s1.size(); ++i ) *pS1++ = s1[i];
      }
      if ( nlhs > 1 ) {
        arg_S2 = mxCreateNumericMatrix( s2.size(),1, mxDOUBLE_CLASS, mxREAL );
        double * pS2 = mxGetPr(arg_S2);
        for ( unsigned i = 0; i < s2.size(); ++i ) *pS2++ = s2[i];
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("clothoid intersection failed\n");
    }
  }
}
