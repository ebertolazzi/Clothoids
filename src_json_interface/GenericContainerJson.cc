//
//  main.cpp
//  genconjson
//
//  Created by Nicola Dal Bianco on 08/11/17.
//  Copyright © 2017 Nicola Dal Bianco. All rights reserved.
//

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wcovered-switch-default"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wimplicit-fallthrough"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wcovered-switch-default"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include "GenericContainerJson.hh"
#include "GenericContainerJsonHandler.hh"

#include <iostream>
#include <string>
#include <iomanip>

#ifdef USE_SYSTEM_JSON
  #include <rapidjson/prettywriter.h>
  #include <rapidjson/writer.h>
  #include <rapidjson/reader.h>
#else
  #include "rapidjson/prettywriter.h"
  #include "rapidjson/writer.h"
  #include "rapidjson/reader.h"
#endif


using namespace GC;
using namespace rapidjson;
using namespace std;

/*
 The encoding of a GenericContainer to a Json is much simpler than the reverse operation, and it is
 done basically by the C-like function 'gc_to_writer'.
 The implementation of the header GenericContainerJson.hh then follows.
 */

namespace GenericContainerNamespace {

  static
  inline
  bool isZero ( real_type x )
  { return FP_ZERO == fpclassify( x ); }

  void
  real_to_stream ( real_type number, ostream_type & out ) {
    out << std::setprecision( std::numeric_limits<real_type>::digits10 + 1 )
        << std::scientific
        << number;
  }

  void
  complex_to_stream ( complex_type const & number,
                      ostream_type       & out,
                      string       const & im_unit )
  {
    if ( number == complex_type ( 0, 0 ) ) {
      out << "0";
      return;
    }
    if ( !isZero(number.real()) ) {
      real_to_stream ( number.real(), out );
    }
    if ( isZero(number.imag()) ) {
      return;
    }
    if ( !isZero(number.real()) ) {
      out << "+";
    }
    real_to_stream ( number.imag(), out );
    out << im_unit;
  }

  // Writer and PrettyWriter do not share any virtual method, so we need to use a template
  template <typename W>
  static
  void
  gc_to_writer( GenericContainer const & gc,
                W                      & writer,
                string const           & im_unit,
                GCJsonMatrixOrder        mat_order ) {

    switch ( gc.get_type() ) {
    case GC_NOTYPE: {
      string null = "null";
      writer.String ( null.c_str() );
      break;
    }
    case GC_POINTER: {
      stringstream ss;
      ss << gc.get_pointer<void *>();
      writer.String ( ss.str().c_str(), SizeType(ss.str().length()) );
      break;
    }
    case GC_BOOL: {
      writer.Bool ( gc.get_bool() );
      break;
    }
    case GC_INTEGER: {
      writer.Int ( gc.get_int() );
      break;
    }
    case GC_LONG: {
      writer.Int64 ( gc.get_long() );
      break;
    }
    case GC_REAL: {
      writer.Double ( gc.get_real() );
      break;
    }
    case GC_STRING: {
      writer.String ( gc.get_string().c_str(), SizeType(gc.get_string().length()) );
      break;
    }
    case GC_COMPLEX: {
      stringstream ss;
      complex_to_stream ( gc.get_complex(), ss, im_unit );
      writer.String ( ss.str().c_str(), SizeType(ss.str().length()) );
      break;
    }
    //vector type
    case GC_VEC_POINTER: {
      writer.StartArray();
      for ( void * pointer : gc.get_vec_pointer() ) {
        stringstream ss;
        ss << pointer;
        writer.String ( ss.str().c_str(), SizeType(ss.str().length()) );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_BOOL: {
      writer.StartArray();
      for ( bool value : gc.get_vec_bool() ) {
        writer.Bool ( value );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_INTEGER: {
      writer.StartArray();
      for ( int_type value : gc.get_vec_int() ) {
        writer.Int ( value );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_LONG: {
      writer.StartArray();
      for ( long_type value : gc.get_vec_long() ) {
        writer.Int64 ( value );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_REAL: {
      writer.StartArray();
      for ( real_type value : gc.get_vec_real() ) {
        writer.Double ( value );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_STRING: {
      writer.StartArray();
      for ( string const & value : gc.get_vec_string() ) {
        writer.String ( value.c_str(), SizeType(value.length()) );
      }
      writer.EndArray();
      break;
    }
    case GC_VECTOR: {
      writer.StartArray();
      for ( GenericContainer const & value : gc.get_vector() ) {
        gc_to_writer ( value, writer, im_unit, mat_order );
      }
      writer.EndArray();
      break;
    }
    case GC_VEC_COMPLEX: {
      writer.StartArray();
      for ( complex_type value : gc.get_vec_complex() ) {
        stringstream ss;
        complex_to_stream ( value, ss, im_unit );
        writer.String ( ss.str().c_str(), SizeType(ss.str().length()) );
      }
      writer.EndArray();
      break;
    }
    // map type
    case GC_MAP: {
      writer.StartObject();
      map_type const & m = gc.get_map();
      for ( map_type::const_iterator it = m.begin(); it != m.end(); ++it ) {
        writer.Key ( it->first.c_str() );
        gc_to_writer ( it->second, writer, im_unit, mat_order );
      }
      writer.EndObject();
      break;
    }
    // matrix type
    case GC_MAT_INTEGER: {
      writer.StartArray();
      mat_int_type mat = gc.get_mat_int();
      if ( mat_order == row_major ) { // if row_major
        for ( unsigned int i_row = 0; i_row < mat.numRows(); i_row++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_int();
          vec_int_type & vec = tmp_gc.get_vec_int();
          mat.getRow ( i_row, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      } else { // else: column majot (default)
        for ( unsigned int i_col = 0; i_col < mat.numCols(); i_col++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_int();
          vec_int_type & vec = tmp_gc.get_vec_int();
          mat.getColumn ( i_col, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      }
      writer.EndArray();
      break;
    }
    case GC_MAT_LONG: {
      writer.StartArray();
      mat_long_type mat = gc.get_mat_long();
      if ( mat_order == row_major ) { // if row_major
        for ( unsigned int i_row = 0; i_row < mat.numRows(); i_row++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_long();
          vec_long_type & vec = tmp_gc.get_vec_long();
          mat.getRow ( i_row, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      } else { // else: column majot (default)
        for ( unsigned int i_col = 0; i_col < mat.numCols(); i_col++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_long();
          vec_long_type & vec = tmp_gc.get_vec_long();
          mat.getColumn ( i_col, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      }
      writer.EndArray();
      break;
    }
    case GC_MAT_REAL: {
      writer.StartArray();
      mat_real_type mat = gc.get_mat_real();
      if ( mat_order == row_major ) { // if row_major
        for ( unsigned int i_row = 0; i_row < mat.numRows(); i_row++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_real();
          vec_real_type & vec = tmp_gc.get_vec_real();
          mat.getRow ( i_row, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      } else { // else: column majot (default)
        for ( unsigned int i_col = 0; i_col < mat.numCols(); i_col++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_real();
          vec_real_type & vec = tmp_gc.get_vec_real();
          mat.getColumn ( i_col, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      }
      writer.EndArray();
      break;
    }
    case GC_MAT_COMPLEX: {
      writer.StartArray();
      mat_complex_type mat = gc.get_mat_complex();
      if ( mat_order == row_major ) { // if row_major
        for ( unsigned int i_row = 0; i_row < mat.numRows(); i_row++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_complex();
          vec_complex_type & vec = tmp_gc.get_vec_complex();
          mat.getRow ( i_row, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      } else { // else: column majot (default)
        for ( unsigned int i_col = 0; i_col < mat.numCols(); i_col++ ) {
          GenericContainer tmp_gc;
          tmp_gc.set_vec_complex();
          vec_complex_type & vec = tmp_gc.get_vec_complex();
          mat.getColumn ( i_col, vec );
          gc_to_writer ( tmp_gc, writer, im_unit, mat_order );
        }
      }
      writer.EndArray();
      break;
    }
    default: {
      string null = "null";
      writer.String ( null.c_str() );
    }
    } // end swicth
  }

  std::string
  genericContainerToJsonString(
    GenericContainer const & gc,
    GenericContainer const & gc_options
  ) {
    // check if prettify
    bool pretty = false;
    gc_options.get_if_exists ( GC_JSON_PRETTY, pretty );

    // check immaginary unit
    string im_unit = "i";
    gc_options.get_if_exists ( GC_JSON_IM_UNIT, im_unit );

    // check matrix ordering
    int_type mat_order = GCJsonMatrixOrder::column_major;
    gc_options.get_if_exists ( GC_JSON_MAT_ORDER, mat_order );

    StringBuffer buffer;

    // prettify if pretty
    if ( pretty ) {
      PrettyWriter<StringBuffer> writer ( buffer );
      writer.SetIndent ( ' ', 0 );

      int num_char = 4;
      char ch = ' ';

      gc_options.get_if_exists ( GC_JSON_INDENT_NUM, num_char );
      if ( num_char < 0 ) {
        num_char = 0;
      }

      string ch_str;
      if ( gc_options.get_if_exists ( GC_JSON_INDENT_CHAR, ch_str ) ) {
        ch = * ( ch_str.c_str() );
      }

      // now parse
      writer.SetIndent ( ch, unsigned(num_char) );
      gc_to_writer ( gc, writer, im_unit, GCJsonMatrixOrder(mat_order) );
      return string ( buffer.GetString(), buffer.GetSize() );
    }

    // Otherwise use standard writer
    Writer<StringBuffer> writer ( buffer );

    // now parse
    gc_to_writer ( gc, writer, im_unit, GCJsonMatrixOrder(mat_order) );
    return string ( buffer.GetString(), buffer.GetSize() );
  }

  void
  jsonStringToGenericContainer(
    std::string      const & json,
    GenericContainer       & gc_output,
    GenericContainer const & gc_options_in
  ) {
    // clean output
    gc_output.clear();

    //set the defaults options
    GenericContainer gc_options;

    gc_options[GC_JSON_MAT_ORDER].set_int ( GCJsonMatrixOrder::column_major );
    gc_options_in.get_if_exists ( GC_JSON_MAT_ORDER, gc_options[GC_JSON_MAT_ORDER].get_int() );

    gc_options[GC_JSON_IM_UNIT].set_string ( "i" );
    gc_options_in.get_if_exists ( GC_JSON_IM_UNIT, gc_options[GC_JSON_IM_UNIT].get_string() );

    // create the handler and the requested objects
    GenericContainerJsonHandler gc_handler ( gc_output, gc_options );
    Reader reader;
    StringStream ss ( json.c_str() );

    // now parse
    reader.Parse ( ss, gc_handler );
  }

}
