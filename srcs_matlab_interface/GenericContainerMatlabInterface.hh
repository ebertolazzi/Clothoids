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

#ifndef GENERIC_CONTAINER_HH
  #include "GenericContainer.hh"
#endif
#include "mex.h"

namespace GenericContainerNamespace {

  void mxArray_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) ;
  void GenericContainer_to_mxArray( GenericContainer const & gc, mxArray * & mx ) ;
  void mexPrint( GenericContainer const & gc ) ;
  void mxSparse_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) ;

}
