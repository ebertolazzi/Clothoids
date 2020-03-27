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
    GenericContainer gc;
    mxArray_to_GenericContainer( prhs[0], gc );
    mexPrint(gc);
  }
  catch ( std::exception & exc ) {
    mexPrintf("Error: %s\n", exc.what() );
  }
  catch (...) {
    mexPrintf("Unknown erroe\n");
  }
}

}
