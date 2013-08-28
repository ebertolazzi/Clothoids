/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: GenericContainerCinterface.cc
///

/*
 \file GenericContainerCinterface.cc
 This file contains the sources for the C interface to `GenericContainer`
 */

#include "GenericContainer.hh"
#include "GenericContainerCinterface.h"

#include <vector>
#include <map>
#include <string>
#include <deque>

using namespace std ;

#define EXTERN_C extern "C"

static std::map<std::string,GenericContainer> gc_stored ;
std::deque<GenericContainer*>                 head ;

static
int
check( int data_type ) {
  if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD ;
  if ( data_type == head.back() -> get_type() ) {
    if ( head.back() == nullptr ) return GENERIC_CONTAINER_NO_DATA ;
    return GENERIC_CONTAINER_OK ;
  } else {
    return GENERIC_CONTAINER_BAD_TYPE ;
  }
}

static
int
check_no_data( int data_type ) {
  if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD ;
  if ( GenericContainer::NOTYPE == head.back() -> get_type() ||
       data_type == head.back() -> get_type() ) return GENERIC_CONTAINER_OK ;
  return GENERIC_CONTAINER_NOT_EMPTY ;
}

EXTERN_C
int
GC_select( char const name[] ) {
  head.clear() ;
  head.push_back(&gc_stored[name]) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_delete( char const name[] ) {
  gc_stored[name].clear() ;
  head.clear() ;
  return GENERIC_CONTAINER_OK ;
}

int
GC_pop_head() {
  if ( head.empty() ) {
    return GENERIC_CONTAINER_NO_DATA ;
  } else {
    head.pop_back() ;
  }
  return GENERIC_CONTAINER_OK ;
}

int
GC_print() {
  if ( head.empty() ) {
    return GENERIC_CONTAINER_NO_DATA ;
  } else {
    head.front() -> print(cout) ;
    return GENERIC_CONTAINER_OK ;
  }
}

void *
GC_mem_ptr( char const name[] ) {
  // check if exists ?
  return (void*) & gc_stored[name] ;
}

int
GC_set_bool( bool const a ) {
  int ok = check_no_data( GenericContainer::BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_bool(a) ;
  return ok ;
}

int
GC_set_int( int const a ) {
  int ok = check_no_data( GenericContainer::INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_int(a) ;
  return ok ;
}

int
GC_set_real( double const a ) {
  int ok = check_no_data( GenericContainer::REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_real(a) ;
  return ok ;
}

int
GC_set_string( char const a[] ) {
  int ok = check_no_data( GenericContainer::STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_string(a) ;
  return ok ;
}

int
GC_set_vector_of_int( int const a[], int nelem ) {
  int ok = check_no_data( GenericContainer::VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_int( nelem ) ;
    GenericContainer::vec_int_type & v = head.back() -> get_vec_int() ;
    std::copy( a, a + nelem, v.begin() ) ;
  }
  return ok ;
}

int
GC_set_empty_vector_of_int() {
  int ok = check_no_data( GenericContainer::VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_int() ;
  return ok ;
}

int
GC_push_int( int const a ) {
  int ok = check( GenericContainer::VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GenericContainer::vec_int_type & v = head.back() -> get_vec_int() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GenericContainer::VECTOR )) == GENERIC_CONTAINER_OK ) {
    GenericContainer::vector_type & v = head.back() -> get_vector() ;
    GenericContainer tmp ;
    tmp . set_int( a ) ;
    v.push_back( tmp ) ;
  }
  return ok ;  
}

int
GC_set_vector_of_real( double const a[], int nelem ) {
  int ok = check_no_data( GenericContainer::VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_real( nelem ) ;
    GenericContainer::vec_real_type & v = head.back() -> get_vec_real() ;
    std::copy( a, a + nelem, v.begin() ) ;
  }
  return ok ;
}

int
GC_set_empty_vector_of_real() {
  int ok = check_no_data( GenericContainer::VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_real() ;
  return ok ;
}

int
GC_push_real( double const a ) {
  int ok = check( GenericContainer::VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GenericContainer::vec_real_type & v = head.back() -> get_vec_real() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GenericContainer::VECTOR )) == GENERIC_CONTAINER_OK ) {
    GenericContainer::vector_type & v = head.back() -> get_vector() ;
    v.push_back( a ) ;
  }
  return ok ;
}

int
GC_set_vector_of_string( char const *a[], int nelem ) {
  int ok = check_no_data( GenericContainer::VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_string( nelem ) ;
    GenericContainer::vec_string_type & v = head.back() -> get_vec_string() ;
    for ( unsigned i = 0 ; i < nelem ; ++i ) v[i] = a[i] ;
  }
  return ok ;
}

int
GC_set_empty_vector_of_string() {
  int ok = check_no_data( GenericContainer::VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_string() ;
  return ok ;
}

int
GC_push_string( char const a[] ) {
  int ok = check( GenericContainer::VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GenericContainer::vec_string_type & v = head.back() -> get_vec_string() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GenericContainer::VECTOR )) == GENERIC_CONTAINER_OK ) {
    GenericContainer::vector_type & v = head.back() -> get_vector() ;
    v.push_back( a ) ;
  }
  return ok ;
}

int
GC_set_vector( int nelem ) {
  int ok = check_no_data( GenericContainer::VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vector( nelem ) ;
  return ok ;
}

int
GC_set_empty_vector() {
  int ok = check_no_data( GenericContainer::VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vector() ;
  return ok ;
}

int
GC_set_vector_position( int pos ) {
  int ok = check( GenericContainer::VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GenericContainer * gc = &(*head.back())[pos] ;
    head.push_back(gc) ;
  }
  return ok ;
}

int
GC_set_map() {
  int ok = check_no_data( GenericContainer::MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_map() ;
  return ok ;
}

/*! \brief
 *  Set the position of insertion point is at the `pos` element
 *  of the actual map.
 */
int
GC_set_map_position( char const pos[] ) {
  int ok = check( GenericContainer::MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GenericContainer * gc = &(*head.back())[pos] ;
    head.push_back(gc) ;
  }
  return ok ;
}

///
/// eof: GenericContainerCinterface.cc
///
