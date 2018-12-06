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

  void mxArray_to_GenericContainer( mxArray const * mx, GenericContainer & gc );
  void mxSparse_to_GenericContainer( mxArray const * mx, GenericContainer & gc );

  void GenericContainer_to_mxArray( GenericContainer const & gc, mxArray * & mx );

  void mexPrint( GenericContainer const & gc );

  // ===========================================================================

  void to_mxArray( bool             const & val, mxArray * & mx );
  void to_mxArray( int_type         const & val, mxArray * & mx );
  void to_mxArray( long_type        const & val, mxArray * & mx );
  void to_mxArray( real_type        const & val, mxArray * & mx );
  void to_mxArray( complex_type     const & val, mxArray * & mx );
  void to_mxArray( string_type      const & val, mxArray * & mx );
  void to_mxArray( vec_bool_type    const & val, mxArray * & mx );
  void to_mxArray( vec_int_type     const & val, mxArray * & mx );
  void to_mxArray( vec_long_type    const & val, mxArray * & mx );
  void to_mxArray( vec_real_type    const & val, mxArray * & mx );
  void to_mxArray( vec_complex_type const & val, mxArray * & mx );
  void to_mxArray( vec_string_type  const & val, mxArray * & mx );
  void to_mxArray( mat_int_type     const & val, mxArray * & mx );
  void to_mxArray( mat_long_type    const & val, mxArray * & mx );
  void to_mxArray( mat_real_type    const & val, mxArray * & mx );
  void to_mxArray( mat_complex_type const & val, mxArray * & mx );

}
