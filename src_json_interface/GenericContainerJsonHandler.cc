//
//  main.cpp
//  genconjson
//
//  Created by Nicola Dal Bianco on 08/11/17.
//  Copyright © 2017 Nicola Dal Bianco. All rights reserved.
//

#include "GenericContainerJsonHandler.hh"
#include <iostream>
#include <iomanip>
#include <regex>
#include <string>
#include <algorithm>

#ifdef USE_SYSTEM_JSON
  #include <rapidjson/reader.h>
#else
  #include "rapidjson/reader.h"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wimplicit-fallthrough"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wcovered-switch-default"
#endif

#define ASSERT_STACK_SIZE \
  if (_gc_stack.size()==0) \
    throw runtime_error("Stack size is null");

using namespace GC;
using namespace rapidjson;
using namespace std;

/*\
  Json objects are quite general objects that here we want to convert
  into a GenericContainer.
  Due to the etherogeneus nature of Json, some 'awful' C-like functions
  are required for parsing: these functions are written below.
  After such functions, the implementation of the GenericContainerJsonHandler class follows.
\*/

static
inline
bool isZero ( real_type x )
{ return FP_ZERO == fpclassify( x ); }

static
inline
bool isInteger32 ( real_type x )
{ return isZero ( x - static_cast<int_type> ( floor ( x ) ) ); }

static
inline
bool isInteger64 ( real_type x )
{ return isZero ( x - static_cast<long_type> ( floor ( x ) ) ); }

static
vector<string>
splitString( string const & str, string const & sep )
{
  char * cstr = const_cast<char *> ( str.c_str() );
  char * current;
  vector<string> arr;
  current = strtok ( cstr, sep.c_str() );
  while ( current != nullptr ) {
    arr.push_back ( current );
    current = strtok ( nullptr, sep.c_str() );
  }
  return arr;
}

static
bool
isStringARealNumber( string const & str_input, double & number ) {
  if ( str_input.find_first_not_of ( "0123456789.e-E" ) != std::string::npos ) {
    return false;
  }
  try {
    double real = stod ( str_input );
    number = real;
    return true;
  } catch ( ... ) {}
  return false;
}

static
bool
isStringAnImaginaryNumber(
  string const & str_input,
  double       & imaginary_part,
  string const & im_unit
) {
  if ( str_input.length() == 0 ) {
    return false;
  }
  if ( str_input.find ( im_unit ) == std::string::npos ) {
    return false;
  }
  if ( str_input.compare ( im_unit ) == 0 ) {
    imaginary_part = 1;
    return true;
  }
  string str = str_input.substr ( 0, str_input.length() - im_unit.length() ); // remove imaginary unit (if present)
  return isStringARealNumber ( str, imaginary_part );
}

static
bool
isStringAComplexNumber(
  string const & str_input,
  complex_type & out_number,
  string const & im_unit
) {
  if ( str_input.find ( "+" ) == std::string::npos ) {
    // in this case, the + is not found, thus it is either a pure real or immaginary number
    double num;
    if ( isStringARealNumber ( str_input, num ) ) {
      out_number = num;
      return true;
    }
    if ( isStringAnImaginaryNumber ( str_input, num, im_unit ) ) {
      out_number = complex_type ( 0, num );
      return true;
    }
    return false;
  }
  // in this case, the + was found
  vector<string> real_imag = splitString ( str_input, "+" );
  if ( real_imag.size() != 2 ) {
    return false;
  }
  double real, imag;
  if ( ( isStringARealNumber ( real_imag[0], real ) ) && ( isStringAnImaginaryNumber ( real_imag[1], imag, im_unit ) ) ) {
    out_number = complex_type ( real, imag );
    return true;
  }
  return false;
}

static
TypeAllowed
findMinimumCommonTypeAllowed( TypeAllowed type1, TypeAllowed type2 ) {

  if ( type1 == type2 ) return type1;

  switch ( type1 ) {
  //the following types are not compatible with any other type different from them
  case GC_NOTYPE:
  case GC_BOOL:
  case GC_VEC_BOOL:
  case GC_MAP:
  case GC_POINTER:
  case GC_VEC_POINTER:
  case GC_VEC_STRING:
  case GC_VECTOR:
    return GC_NOTYPE;

  case GC_STRING: {
    if ( ( type2 == GC_COMPLEX ) ||
         ( type2 == GC_REAL ) ||
         ( type2 == GC_LONG ) ||
         ( type2 == GC_INTEGER ) ) {
      return GC_STRING;
    }
    break;
  }

  //numbers
  case GC_INTEGER: {
    if ( ( type2 == GC_LONG ) ||
         ( type2 == GC_REAL ) ||
         ( type2 == GC_COMPLEX ) ) {
      return type2;
    }
    if ( type2 == GC_STRING ) {
      return GC_STRING;
    }
    break;
  }
  case GC_LONG: {
    if ( type2 == GC_INTEGER ) {
      return GC_LONG;
    }
    if ( ( type2 == GC_REAL ) ||
         ( type2 == GC_COMPLEX ) ||
         ( type2 == GC_STRING ) ) {
      return type2;
    }
    break;
  }
  case GC_REAL: {
    if ( ( type2 == GC_INTEGER ) ||
         ( type2 == GC_LONG ) ) {
      return GC_REAL;
    }
    if ( ( type2 == GC_COMPLEX ) ||
         ( type2 == GC_STRING ) ) {
      return type2;
    }
    break;
  }
  case GC_COMPLEX: {
    if ( ( type2 == GC_INTEGER ) ||
         ( type2 == GC_LONG ) ||
         ( type2 == GC_REAL ) ) {
      return GC_COMPLEX;
    }
    if ( type2 == GC_STRING ) {
      return GC_STRING;
    }
    break;
  }

  //vector type
  case GC_VEC_INTEGER: {
    if ( ( type2 == GC_VEC_LONG ) ||
         ( type2 == GC_VEC_REAL ) ||
         ( type2 == GC_VEC_COMPLEX ) ) {
      return type2;
    }
    break;
  }
  case GC_VEC_LONG: {
    if ( type2 == GC_VEC_INTEGER ) {
      return GC_VEC_LONG;
    }
    if ( ( type2 == GC_VEC_REAL ) ||
         ( type2 == GC_VEC_COMPLEX ) ) {
      return type2;
    }
    break;
  }
  case GC_VEC_REAL: {
    if ( ( type2 == GC_VEC_INTEGER ) ||
         ( type2 == GC_VEC_LONG ) ) {
      return GC_VEC_REAL;
    }
    if ( type2 == GC_VEC_COMPLEX ) {
      return type2;
    }
    break;
  }
  case GC_VEC_COMPLEX: {
    if ( ( type2 == GC_VEC_INTEGER ) ||
         ( type2 == GC_VEC_LONG ) ||
         ( type2 == GC_VEC_REAL ) ) {
      return GC_VEC_COMPLEX;
    }
    break;
  }

  //matrix type
  case GC_MAT_INTEGER: {
    if ( ( type2 == GC_MAT_LONG ) ||
         ( type2 == GC_MAT_REAL ) ||
         ( type2 == GC_MAT_COMPLEX ) ) {
      return type2;
    }
    break;
  }
  case GC_MAT_LONG: {
    if ( type2 == GC_MAT_INTEGER ) {
      return GC_MAT_LONG;
    }
    if ( ( type2 == GC_MAT_REAL ) ||
         ( type2 == GC_MAT_COMPLEX ) ) {
      return type2;
    }
    break;
  }
  case GC_MAT_REAL: {
    if ( ( type2 == GC_MAT_INTEGER ) ||
         ( type2 == GC_MAT_LONG ) ) {
      return GC_MAT_REAL;
    }
    if ( type2 == GC_MAT_COMPLEX ) {
      return type2;
    }
    break;
  }
  case GC_MAT_COMPLEX: {
    if ( ( type2 == GC_MAT_INTEGER ) ||
         ( type2 == GC_MAT_LONG ) ||
         ( type2 == GC_MAT_REAL ) ) {
      return GC_MAT_COMPLEX;
    }
    break;
  }

  }
  return GC_NOTYPE;
}

static
void
convertGenericVectorToBoolVector ( GenericContainer & root ) {
  vector<bool> new_vec;
  for ( GenericContainer const & gc : root.get_vector() )
    new_vec.push_back ( gc.get_bool() );
  root.set_vec_bool ( new_vec );
}

static
void
convertGenericVectorToStringVector(
  GenericContainer & root,
  string     const & im_unit
) {
  vector<string> new_vec;

  for ( GenericContainer const & gc : root.get_vector() ) {
    if ( gc.get_type() == GC_STRING ) {
      new_vec.push_back ( gc.get_string() );
    } else if ( gc.get_type() == GC_INTEGER ) {
      stringstream ss;
      ss << gc.get_int();
      new_vec.push_back ( ss.str() );
    } else if ( gc.get_type() == GC_LONG ) {
      stringstream ss;
      ss << gc.get_long();
      new_vec.push_back ( ss.str() );
    } else if ( gc.get_type() == GC_REAL ) {
      stringstream ss;
      real_to_stream ( gc.get_real(), ss );
      new_vec.push_back ( ss.str() );
    } else if ( gc.get_type() == GC_COMPLEX ) {
      stringstream ss;
      complex_to_stream ( gc.get_complex(), ss, im_unit );
      new_vec.push_back ( ss.str() );
    }
  }
  root.set_vec_string ( new_vec );
}

static
void
convertGenericVectorToIntVector( GenericContainer & root ) {
  vector<int_type> new_vec;
  for ( GenericContainer const & gc : root.get_vector() ) {
    new_vec.push_back ( gc.get_int() );
  }
  root.set_vec_int ( new_vec );
}

static
void
convertGenericVectorToLongVector( GenericContainer & root ) {
  vector<long_type> new_vec;
  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_LONG ) {
      gc.promote_to_long();
    }
    new_vec.push_back ( gc.get_long() );
  }
  root.set_vec_long ( new_vec );
}

static
void
convertGenericVectorToRealVector( GenericContainer & root ) {
  vector<real_type> new_vec;

  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_REAL ) {
      gc.promote_to_real();
    }
    new_vec.push_back ( gc.get_real() );
  }
  root.set_vec_real ( new_vec );
}

static
void
convertGenericVectorToComplexVector( GenericContainer & root ) {
  vector<complex_type> new_vec;

  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_COMPLEX ) {
      gc.promote_to_complex();
    }
    new_vec.push_back ( gc.get_complex() );
  }
  root.set_vec_complex ( new_vec );
}

static
void
convertGenericVectorToIntMat(
  GenericContainer & root,
  GCJsonMatrixOrder  order
) {

  if ( root.get_vector().size() == 0 ) return;

  //first check if they have the same size
  unsigned num_rows = unsigned(root.get_vector() [0].get_vec_int().size());
  for ( GenericContainer const & gc : root.get_vector() ) {
    if ( gc.get_vec_int().size() != num_rows )return;
  }

  unsigned num_cols = unsigned(root.get_vector().size());
  unsigned num_mat_rows, num_mat_cols;
  if ( order == column_major ) {
    num_mat_cols = num_cols;
    num_mat_rows = num_rows;
  } else {
    num_mat_cols = num_rows;
    num_mat_rows = num_cols;
  }
  mat_int_type mat ( num_mat_rows, num_mat_cols );

  for ( unsigned i_col = 0; i_col < num_cols; i_col++ ) {
    vec_int_type vec = root.get_vector() [i_col].get_vec_int();
    for ( unsigned i_row = 0; i_row < num_rows; i_row++ ) {
      if ( order == column_major ) {
        mat ( i_row, i_col ) = vec[i_row];
      } else {
        mat ( i_col, i_row ) = vec[i_row];
      }
    }
  }
  root.set_mat_int ( mat );
}

static
void
convertGenericVectorToLongMat(
  GenericContainer & root,
  GCJsonMatrixOrder  order
) {

  if ( root.get_vector().size() == 0 ) return;

  //first check if they have the same size
  bool started = false;
  unsigned num_rows = 0;
  for ( GenericContainer & gc : root.get_vector() ) {
    switch ( gc.get_type() ) {
    case GC_VEC_INTEGER: {
      if ( started ) {
        if ( gc.get_vec_int().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_int().size());
        started = true;
      }
      break;
    }
    case GC_VEC_LONG: {
      if ( started ) {
        if ( gc.get_vec_long().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_long().size());
        started = true;
      }
      break;
    }
    default:
      return;
    }
  }

  // now promote
  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_VEC_LONG ) gc.promote_to_vec_long();
  }

  unsigned num_cols = unsigned(root.get_vector().size());
  unsigned num_mat_rows, num_mat_cols;
  if ( order == column_major ) {
    num_mat_cols = num_cols;
    num_mat_rows = num_rows;
  } else {
    num_mat_cols = num_rows;
    num_mat_rows = num_cols;
  }
  mat_long_type mat ( num_mat_rows, num_mat_cols );

  for ( unsigned i_col = 0; i_col < num_cols; i_col++ ) {
    vec_long_type vec = root.get_vector() [i_col].get_vec_long();
    for ( unsigned i_row = 0; i_row < num_rows; i_row++ ) {
      if ( order == column_major ) mat( i_row, i_col ) = vec[i_row];
      else                         mat( i_col, i_row ) = vec[i_row];
    }
  }
  root.set_mat_long ( mat );
}

static
void
convertGenericVectorToRealMat(
  GenericContainer & root,
  GCJsonMatrixOrder  order
) {

  if ( root.get_vector().size() == 0 ) return;

  //first check if they have the same size
  bool started = false;
  unsigned num_rows = 0;
  for ( GenericContainer & gc : root.get_vector() ) {
    switch ( gc.get_type() ) {
    case GC_VEC_INTEGER: {
      if ( started ) {
        if ( gc.get_vec_int().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_int().size());
        started = true;
      }
      break;
    }
    case GC_VEC_LONG: {
      if ( started ) {
        if ( gc.get_vec_long().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_long().size());
        started = true;
      }
      break;
    }
    case GC_VEC_REAL: {
      if ( started ) {
        if ( gc.get_vec_real().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_real().size());
        started = true;
      }
      break;
    }
    default:
      return;
    }
  }

  // now promote
  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_VEC_REAL ) {
      gc.promote_to_vec_real();
    }
  }

  unsigned num_cols = unsigned(root.get_vector().size());
  unsigned num_mat_rows, num_mat_cols;
  if ( order == column_major ) {
    num_mat_cols = num_cols;
    num_mat_rows = num_rows;
  } else {
    num_mat_cols = num_rows;
    num_mat_rows = num_cols;
  }
  mat_real_type mat ( num_mat_rows, num_mat_cols );

  for ( unsigned i_col = 0; i_col < num_cols; i_col++ ) {
    vec_real_type vec = root.get_vector() [i_col].get_vec_real();
    for ( unsigned i_row = 0; i_row < num_rows; i_row++ ) {
      if ( order == column_major ) {
        mat ( i_row, i_col ) = vec[i_row];
      } else {
        mat ( i_col, i_row ) = vec[i_row];
      }
    }
  }
  root.set_mat_real ( mat );
}

static
void
convertGenericVectorToComplexMat(
  GenericContainer & root,
  GCJsonMatrixOrder  order
) {

  if ( root.get_vector().size() == 0 ) return;

  //first check if they have the same size
  bool started = false;
  unsigned num_rows = 0;
  for ( GenericContainer & gc : root.get_vector() ) {
    switch ( gc.get_type() ) {
    case GC_VEC_INTEGER: {
      if ( started ) {
        if ( gc.get_vec_int().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_int().size());
        started = true;
      }
      break;
    }
    case GC_VEC_LONG: {
      if ( started ) {
        if ( gc.get_vec_long().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_long().size());
        started = true;
      }
      break;
    }
    case GC_VEC_REAL: {
      if ( started ) {
        if ( gc.get_vec_real().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_real().size());
        started = true;
      }
      break;
    }
    case GC_VEC_COMPLEX: {
      if ( started ) {
        if ( gc.get_vec_complex().size() != num_rows ) {
          return;
        }
      } else {
        num_rows = unsigned(root.get_vector() [0].get_vec_complex().size());
        started = true;
      }
      break;
    }
    default:
      return;
    }
  }

  // now promote
  for ( GenericContainer & gc : root.get_vector() ) {
    if ( gc.get_type() != GC_VEC_COMPLEX ) gc.promote_to_vec_complex();
  }

  unsigned num_cols = unsigned(root.get_vector().size());
  unsigned num_mat_rows, num_mat_cols;
  if ( order == column_major ) {
    num_mat_cols = num_cols;
    num_mat_rows = num_rows;
  } else {
    num_mat_cols = num_rows;
    num_mat_rows = num_cols;
  }
  mat_complex_type mat ( num_mat_rows, num_mat_cols );

  for ( unsigned i_col = 0; i_col < num_cols; i_col++ ) {
    vec_complex_type vec = root.get_vector() [i_col].get_vec_complex();
    for ( unsigned i_row = 0; i_row < num_rows; i_row++ ) {
      if ( order == column_major ) {
        mat ( i_row, i_col ) = vec[i_row];
      } else {
        mat ( i_col, i_row ) = vec[i_row];
      }
    }
  }
  root.set_mat_complex ( mat );
}

// GenericContainerJsonHandler Class Implementation

GenericContainerJsonHandler::GenericContainerJsonHandler(
  GenericContainer       & gc_output,
  GenericContainer const & gc_options
) {
  int_type mat_order = GCJsonMatrixOrder::column_major;
  gc_options.get_if_exists ( GC_JSON_MAT_ORDER, mat_order );
  _mat_order = GCJsonMatrixOrder(mat_order);

  _im_unit = "i";
  gc_options.get_if_exists ( GC_JSON_IM_UNIT, _im_unit );

  _gc_stack = {stack_entry ( &gc_output, false ) };
}

GenericContainer *
GenericContainerJsonHandler::getCurrentGCPointer() const {
  return _gc_stack.back().first;
}

void
GenericContainerJsonHandler::setCurrentGCPointerArrayType( bool is_array ) {
  _gc_stack.back().second = is_array;
}

bool
GenericContainerJsonHandler::isCurrentGCPointerArrayType() const {
  return _gc_stack.back().second;
}

bool
GenericContainerJsonHandler::Null() {
  ASSERT_STACK_SIZE
  // do nothing
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Bool( bool b ) {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_bool( b );
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Int( int i ) {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_int ( i );
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Uint ( unsigned int i ) {
  ASSERT_STACK_SIZE
  if ( isInteger32 ( i ) ) {
    getCurrentGCPointer()->set_int ( static_cast<int> ( i ) );
  } else {
    getCurrentGCPointer()->set_long ( static_cast<long_type> ( i ) );
  }
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Int64 ( int64_t i ) {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_long ( i );
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Uint64 ( uint64_t i ) {
  ASSERT_STACK_SIZE
  if ( isInteger64 ( i ) ) {
    getCurrentGCPointer()->set_long ( static_cast<long_type> ( i ) );
  } else {
    getCurrentGCPointer()->set_real ( i );
  }
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::Double ( double d ) {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_real ( d );
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::String(
  const char * str,
  SizeType     length,
  bool         /* copy */
) {
  ASSERT_STACK_SIZE
  string stringa ( str, length );
  string stringa_lower = stringa.substr ( 0, 4 );
  transform ( stringa_lower.begin(),
              stringa_lower.end(),
              stringa_lower.begin(),
              [] ( unsigned char c )->unsigned char
                   { return static_cast<unsigned char>(std::tolower( c )); } );

  if ( ( stringa_lower.length() != stringa.length() ) || ( stringa_lower.compare ( "null" ) != 0 ) ) {
    // check if it is a complex number
    complex_type complex;
    if ( isStringAComplexNumber ( stringa, complex, _im_unit ) ) {
      if ( isZero( complex.imag() ) ) {
        real_type re = complex.real();
        if ( isInteger32 ( re ) ) {
          getCurrentGCPointer()->set_int ( static_cast<int_type> ( re ) );
        } else {
          getCurrentGCPointer()->set_real ( re );
        }
      } else {
        getCurrentGCPointer()->set_complex ( complex );
      }
    } else { //in this case it is a string
      getCurrentGCPointer()->set_string ( stringa );
    }
  } else {
    // do not insert nothing
  }
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::StartObject() {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_map();
  return true;
}

bool
GenericContainerJsonHandler::Key(
  const char * str,
  SizeType     length,
  bool         /* copy */
) {
  string key ( str, length );
  _gc_stack.push_back ( stack_entry ( & ( getCurrentGCPointer()->operator[] ( key ) ), false ) );
  return true;
}

bool
GenericContainerJsonHandler::EndObject( SizeType /* member_count */ ) {
  ASSERT_STACK_SIZE
  advanceCurrentGCPointer();
  return true;
}

bool
GenericContainerJsonHandler::StartArray() {
  ASSERT_STACK_SIZE
  getCurrentGCPointer()->set_vector();
  setCurrentGCPointerArrayType( true );
  vector_type & vec = getCurrentGCPointer()->get_vector();
  vec.push_back( GenericContainer() );
  _gc_stack.push_back ( stack_entry ( & ( vec.back() ), false ) );
  return true;
}

bool
GenericContainerJsonHandler::EndArray( SizeType /* member_count */ ) {
  ASSERT_STACK_SIZE
  // go one level up to the stack: point to the real vector being constructed
  _gc_stack.pop_back();

  //remove last empty element
  vector_type & vec = getCurrentGCPointer()->get_vector();
  vec.pop_back();

  //check if the array is homogeneous or etherogeneus
  finalizeArrayProcess();

  //advance pointer
  advanceCurrentGCPointer();
  return true;
}

void
GenericContainerJsonHandler::advanceCurrentGCPointer() {
  ASSERT_STACK_SIZE
  _gc_stack.pop_back();

  if ( ( _gc_stack.size() > 0 ) && ( isCurrentGCPointerArrayType() ) ) {
    vector_type & vec = getCurrentGCPointer()->get_vector();
    vec.push_back ( GenericContainer() );
    _gc_stack.push_back ( stack_entry ( & ( vec.back() ), false ) );
  }
}

void
GenericContainerJsonHandler::finalizeArrayProcess() {
  vector_type & vec = getCurrentGCPointer()->get_vector();

  if ( vec.size() == 0 ) return;

  TypeAllowed minimum_common_type = vec[0].get_type();
  for ( GenericContainer const  & gc : vec ) {
    minimum_common_type = findMinimumCommonTypeAllowed ( minimum_common_type, gc.get_type() );
  }

  if ( minimum_common_type == GC_NOTYPE ) return;

  if ( minimum_common_type == GC_BOOL ) {
    convertGenericVectorToBoolVector( *getCurrentGCPointer() );
    return;
  } else if ( minimum_common_type == GC_STRING ) {
    convertGenericVectorToStringVector( *getCurrentGCPointer(), _im_unit );
    return;
  } else if ( minimum_common_type == GC_INTEGER ) {
    convertGenericVectorToIntVector( *getCurrentGCPointer() );
    return;
  } else if ( minimum_common_type == GC_LONG ) {
    convertGenericVectorToLongVector( *getCurrentGCPointer() );
    return;
  } else if ( minimum_common_type == GC_REAL ) {
    convertGenericVectorToRealVector( *getCurrentGCPointer() );
    return;
  } else if ( minimum_common_type == GC_COMPLEX ) {
    convertGenericVectorToComplexVector( *getCurrentGCPointer() );
    return;
  }

  if ( _mat_order == GCJsonMatrixOrder::none ) return;

  if ( minimum_common_type == GC_VEC_INTEGER ) {
    convertGenericVectorToIntMat( *getCurrentGCPointer(), _mat_order );
  } else if ( minimum_common_type == GC_VEC_LONG ) {
    convertGenericVectorToLongMat( *getCurrentGCPointer(), _mat_order );
  } else if ( minimum_common_type == GC_VEC_REAL ) {
    convertGenericVectorToRealMat( *getCurrentGCPointer(), _mat_order );
  } else if ( minimum_common_type == GC_VEC_COMPLEX ) {
    convertGenericVectorToComplexMat( *getCurrentGCPointer(), _mat_order );
  }

}
