//
//  main.cpp
//  genconjson
//
//  Created by Nicola Dal Bianco on 08/11/17.
//  Copyright © 2017 Nicola Dal Bianco. All rights reserved.
//

#ifndef GC_JSON_HH
#define GC_JSON_HH

#include "GenericContainer.hh"

#define GC_JSON_IM_UNIT      "imaginary_unit"
#define GC_JSON_MAT_ORDER    "matrix_ordering"
#define GC_JSON_PRETTY       "prettify"
#define GC_JSON_INDENT_CHAR  "indent_char"
#define GC_JSON_INDENT_NUM   "indent_num"

namespace GC_namespace {

  enum GCJsonMatrixOrder {
    column_major = 0,   // a matrix is written in Json as vector of vectors
                        // using the column major ordering
    row_major    = 1,   // a matrix is written in Json as vector of vectors
                        // using the row major ordering
    none         = 2    // (option valid only for decoding a json to a
                        // GenericContainer): vector of vectors are not
                        // converted back to matrices
  };

  /*!
   * Converts the content of a GenericContainer (argument gc) to a json string.
   * Matrices are  written as vectors of vectors.
   * gc_options is used to set some string formatting options.
   * gc_options is expcted to contain the following entries:
   *
   * - gc_options[GC_JSON_IM_UNIT]->string (default "i") :
   *   character to use as imaginary unit for complex number representation
   *
   * - gc_options[GC_JSON_MAT_ORDER]->GCJsonMatrixOrder (int_type expected) :
   *   specifies the representation of matrices in the json.
   *   See description of the enum GCJsonMatrixOrder
   *
   * - gc_options[GC_JSON_PRETTY]->bool (default false):
   *   if true the json string comes in prettified-human readable format.
   *
   * - gc_options[GC_JSON_INDENT_CHAR]->string of length 1 (default ' '):
   *   indentation character to use for the prettified output.
   *   Allowed characters are: ' ', '\\t', '\\r', '\\n'
   *
   * - gc_options[GC_JSON_INDENT_NUM]->int (default 4):
   *   number of indentation characters to use for the prettified output.
   */
  std::string
  genericContainerToJsonString(
    GenericContainer const & gc,
    GenericContainer const & gc_options = GenericContainer()
  );

  /*!
   *  Converts a json string to a GenericContainer.
   *  It throws an exception if the parsing of the json string fails.
   *  Strings representing pointers are not converted to pointers but remain strings.
   *  gc_options is used as above, but the 'prettified' and 'indentation'
   *  options does not bother here.
   */
  void
  jsonStringToGenericContainer(
    std::string      const & json,
    GenericContainer       & gc_output,
    GenericContainer const & gc_options = GenericContainer()
  );

  void
  real_to_stream( real_type number, ostream_type & out );

  void
  complex_to_stream(
    complex_type const & number,
    ostream_type       & out,
    std::string  const & im_unit
  );

}

#endif
