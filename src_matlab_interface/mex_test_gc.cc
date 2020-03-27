/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  See file license.txt
\****************************************************************************/

#include "GenericContainer.hh"
#include "GenericContainerMatlabInterface.hh"
#include "mex.h"

#include <sstream>

namespace GenericContainerNamespace {

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  try {

    // Check for proper number of arguments, etc
    /*
    if ( nrhs != 1 ) {
	    mexErrMsgTxt(
"%======================================================================%\n"
"%                                                                      %\n"
"%  Autor: Enrico Bertolazzi                                            %\n"
"%         Department of Industrial Engineering                         %\n"
"%         University of Trento                                         %\n"
"%         enrico.bertolazzi@unitn.it                                   %\n"
"%                                                                      %\n"
"%======================================================================%\n" );
    }
    */

    GenericContainer gc;
    mxArray_to_GenericContainer( prhs[0], gc );
    GenericContainer_to_mxArray( gc, plhs[0] );
    //mexPrint(gc);
  }
  catch ( std::exception & exc ) {
    mexPrintf("Error: %s\n", exc.what() );
  }
  catch (...) {
    mexPrintf("Unknown erroe\n");
  }
}

}
