/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
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
"%======================================================================%\n" ) ;
    }
    */
  
    GenericContainer gc ;
    mxArray_to_GenericContainer( prhs[0], gc ) ;
    GenericContainer_to_mxArray( gc, plhs[0] ) ;
    //mexPrint(gc) ;
  }
  catch ( std::exception & exc ) {
    mexPrintf("Error: %s\n", exc.what() ) ;
  }
  catch (...) {
    mexPrintf("Unknown erroe\n") ;
  }
}

}
