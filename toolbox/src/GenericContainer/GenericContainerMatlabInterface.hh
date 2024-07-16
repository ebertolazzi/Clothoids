/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  See file license.txt
\****************************************************************************/

#ifndef GENERIC_CONTAINER_HH
  #include "GenericContainer.hh"
#endif
#include "mex.h"

namespace GC_namespace {

  //!
  //! Convert an `mxArray` to a `GenericContainer`
  //!
  void mxArray_to_GenericContainer( mxArray const * mx, GenericContainer & gc );

  //!
  //! Convert an `mxArray` containing a sparse matrix to a `GenericContainer`.
  //! The `GenericContainer` will constain a map with the sparse matrix in
  //! compressed column format (https://en.wikipedia.org/wiki/Sparse_matrix)
  //!
  //! - "jc" pointer to the columns
  //! - "ir" index of the rows
  //! - "values" vector of real numbers with nonzeros of the matrix
  //!
  void mxSparse_to_GenericContainer( mxArray const * mx, GenericContainer & gc );

  //!
  //! Convert `GenericContainer` to an `mxArray`
  //!
  void GenericContainer_to_mxArray( GenericContainer const & gc, mxArray * & mx );

  //!
  //! Print the contents of `GenericContainer` to teh `MATLAB` console
  //!
  void mexPrint( GenericContainer const & gc );

  // ===========================================================================

  //!
  //! Convert a boolean to a `mxArray`
  //!
  void to_mxArray( bool const & val, mxArray * & mx );

  //!
  //! Convert an integer to a `mxArray`
  //!
  void to_mxArray( int_type const & val, mxArray * & mx );

  //!
  //! Convert a long integer to a `mxArray`
  //!
  void to_mxArray( long_type const & val, mxArray * & mx );

  //!
  //! Convert a real number to a `mxArray`
  //!
  void to_mxArray( real_type const & val, mxArray * & mx );

  //!
  //! Convert a complex number to a `mxArray`
  //!
  void to_mxArray( complex_type const & val, mxArray * & mx );

  //!
  //! Convert a string integer to a `mxArray`
  //!
  void to_mxArray( string_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of boolean to a `mxArray`
  //!
  void to_mxArray( vec_bool_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of integer to a `mxArray`
  //!
  void to_mxArray( vec_int_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of long integer to a `mxArray`
  //!
  void to_mxArray( vec_long_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of real number to a `mxArray`
  //!
  void to_mxArray( vec_real_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of complex number to a `mxArray`
  //!
  void to_mxArray( vec_complex_type const & val, mxArray * & mx );

  //!
  //! Convert a vector of string to a `mxArray`.
  //! The data will be a cell array.
  //!
  void to_mxArray( vec_string_type const & val, mxArray * & mx );

  //!
  //! Convert a matrix of integer to a `mxArray`
  //!
  void to_mxArray( mat_int_type const & val, mxArray * & mx );

  //!
  //! Convert a matrix of long integer to a `mxArray`
  //!
  void to_mxArray( mat_long_type const & val, mxArray * & mx );

  //!
  //! Convert a matrix of real number to a `mxArray`
  //!
  void to_mxArray( mat_real_type const & val, mxArray * & mx );

  //!
  //! Convert a matrix of complex number to a `mxArray`
  //!
  void to_mxArray( mat_complex_type const & val, mxArray * & mx );

}
