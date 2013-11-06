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

//
// file: GenericContainerCinterface.cc
//

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

typedef std::map<std::string,GC::GenericContainer> MAP_GC ;

static MAP_GC                            gc_stored ;
static std::deque<GC::GenericContainer*> head ;

static GC::map_type         * ptr_map ;
static GC::map_type::iterator map_iterator ;

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
  if ( GC::GC_NOTYPE == head.back() -> get_type() ||
       data_type == head.back() -> get_type() ) return GENERIC_CONTAINER_OK ;
  return GENERIC_CONTAINER_NOT_EMPTY ;
}

EXTERN_C
int
GC_new( char const id[] ) {
  // ckeck if exists
  MAP_GC::iterator it = gc_stored.find(id) ;
  head.clear() ;
  if ( it == gc_stored.end() ) head.push_back(&gc_stored[id]) ;
  else                         head.push_back(&it->second) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_select( char const id[] ) {
  head.clear() ;
  MAP_GC::iterator it = gc_stored.find(id) ;
  if ( it != gc_stored.end() ) {
    head.push_back(&it->second) ;
  } else {
    head.push_back(&gc_stored[id]) ;
  }
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_delete( char const id[] ) {
  head.clear() ;
  MAP_GC::iterator it = gc_stored.find(id) ;
  if ( it != gc_stored.end() ) {
    it->second.clear() ;
    gc_stored.erase(it) ;
  }
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_pop_head() {
  if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA ;
  head.pop_back() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_reset_head() {
  if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA ;
  while ( head.size() > 1 ) head.pop_back() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_print() {
  if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA ;
  head.back() -> print(cout) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_get_type() {
  if ( head.empty() ) return -1 ;
  return head.back() -> get_type() ;
}

EXTERN_C
char const *
GC_get_type_name() {
  static char const empty[] = "" ;
  if ( head.empty() ) return empty ;
  return head.back() -> get_type_name() ;
}

EXTERN_C
void *
GC_mem_ptr( char const name[] ) {
  // check if exists ?
  return (void*) & gc_stored[name] ;
}

EXTERN_C
int
GC_set_bool( int const a ) {
  int ok = check_no_data( GC::GC_BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_bool(a) ;
  return ok ;
}

EXTERN_C
int
GC_get_bool() {
  int ok = check_no_data( GC::GC_BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) return head.back() -> get_bool() ? 1 : 0 ;
  return 0 ;
}

EXTERN_C
int
GC_set_int( int const a ) {
  int ok = check_no_data( GC::GC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_int(a) ;
  return ok ;
}

EXTERN_C
int
GC_get_int() {
  int ok = check_no_data( GC::GC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) return head.back() -> get_int() ;
  return 0 ;
}

EXTERN_C
int
GC_set_real( double const a ) {
  int ok = check_no_data( GC::GC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_real(a) ;
  return ok ;
}

EXTERN_C
double
GC_get_real( ) {
  int ok = check_no_data( GC::GC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) return head.back() -> get_real() ;
  return 0 ;
}

EXTERN_C
int
GC_set_string( char const a[] ) {
  int ok = check_no_data( GC::GC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_string(a) ;
  return ok ;
}

EXTERN_C
char const *
GC_get_string( ) {
  int ok = check_no_data( GC::GC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) return (head.back() -> get_string() . c_str()) ;
  return nullptr ;
}

EXTERN_C
int
GC_set_vector_of_bool( int const a[], int nelem ) {
  int ok = check_no_data( GC::GC_VEC_BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_bool( nelem ) ;
    GC::vec_bool_type & v = head.back() -> get_vec_bool() ;
    for ( int i = 0 ; i < nelem ; ++i ) v[i] = a[i] != 0 ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_empty_vector_of_bool() {
  int ok = check_no_data( GC::GC_VEC_BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_bool() ;
  return ok ;
}

EXTERN_C
int
GC_push_bool( int const a ) {
  int ok = check( GC::GC_VEC_BOOL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GC::vec_bool_type & v = head.back() -> get_vec_bool() ;
    v.push_back( a != 0 ) ;
  } else if ( (ok = check( GC::GC_VECTOR )) == GENERIC_CONTAINER_OK ) {
    GC::vector_type & v = head.back() -> get_vector() ;
    GC::GenericContainer tmp ;
    tmp . set_bool( a != 0 ) ;
    v.push_back( tmp ) ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_vector_of_int( int const a[], int nelem ) {
  int ok = check_no_data( GC::GC_VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_int( nelem ) ;
    GC::vec_int_type & v = head.back() -> get_vec_int() ;
    std::copy( a, a + nelem, v.begin() ) ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_empty_vector_of_int() {
  int ok = check_no_data( GC::GC_VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_int() ;
  return ok ;
}

EXTERN_C
int
GC_push_int( int const a ) {
  int ok = check( GC::GC_VEC_INT ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GC::vec_int_type & v = head.back() -> get_vec_int() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GC::GC_VECTOR )) == GENERIC_CONTAINER_OK ) {
    GC::vector_type & v = head.back() -> get_vector() ;
    GC::GenericContainer tmp ;
    tmp . set_int( a ) ;
    v.push_back( tmp ) ;
  }
  return ok ;  
}

EXTERN_C
int
GC_set_vector_of_real( double const a[], int nelem ) {
  int ok = check_no_data( GC::GC_VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_real( nelem ) ;
    GC::vec_real_type & v = head.back() -> get_vec_real() ;
    std::copy( a, a + nelem, v.begin() ) ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_empty_vector_of_real() {
  int ok = check_no_data( GC::GC_VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_real() ;
  return ok ;
}

EXTERN_C
int
GC_push_real( double const a ) {
  int ok = check( GC::GC_VEC_REAL ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GC::vec_real_type & v = head.back() -> get_vec_real() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GC::GC_VECTOR )) == GENERIC_CONTAINER_OK ) {
    GC::vector_type & v = head.back() -> get_vector() ;
    v.push_back( a ) ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_vector_of_string( char const *a[], unsigned nelem ) {
  int ok = check_no_data( GC::GC_VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    head.back() -> set_vec_string( nelem ) ;
    GC::vec_string_type & v = head.back() -> get_vec_string() ;
    for ( unsigned i = 0 ; i < nelem ; ++i )
      v[i] = a[i] ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_empty_vector_of_string() {
  int ok = check_no_data( GC::GC_VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vec_string() ;
  return ok ;
}

EXTERN_C
int
GC_push_string( char const a[] ) {
  int ok = check( GC::GC_VEC_STRING ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    GC::vec_string_type & v = head.back() -> get_vec_string() ;
    v.push_back( a ) ;
  } else if ( (ok = check( GC::GC_VECTOR )) == GENERIC_CONTAINER_OK ) {
    GC::vector_type & v = head.back() -> get_vector() ;
    v.push_back( a ) ;
  }
  return ok ;
}

EXTERN_C
int
GC_set_vector( int nelem ) {
  int ok = check_no_data( GC::GC_VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vector( nelem ) ;
  return ok ;
}

EXTERN_C
int
GC_set_empty_vector() {
  int ok = check_no_data( GC::GC_VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_vector() ;
  return ok ;
}

EXTERN_C
int
GC_get_vector_size() {
  if ( head.empty() ) return 0 ;
  if ( head.back() == nullptr ) return 0 ;

  switch ( head.back() -> get_type() ) {
  case GC::GC_VEC_BOOL:
    return int(head.back() -> get_vec_bool().size()) ;
  case GC::GC_VEC_INT:
    return int(head.back() -> get_vec_int().size()) ;
  case GC::GC_VEC_REAL:
    return int(head.back() -> get_vec_real().size()) ;
  case GC::GC_VEC_STRING:
    return int(head.back() -> get_vec_string().size()) ;
  case GC::GC_VECTOR:
    return int(head.back() -> get_vector().size()) ;
  default:
    return 0 ;
  }
}

EXTERN_C
int
GC_set_vector_position( int pos ) {
  if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD ;
  if ( head.back() == nullptr ) return GENERIC_CONTAINER_NO_DATA ;
  int ok = check( GC::GC_VECTOR ) ;
  if ( ok == GENERIC_CONTAINER_OK  ) {
    GC::GenericContainer * gc = &(*head.back())[pos] ;
    head.push_back(gc) ;
  }
  return ok ;
}

EXTERN_C
int
GC_get_bool_at_pos( int pos ) {
  if ( ! head.empty() && 
       head.back() != nullptr &&
       head.back() -> get_type() == GC::GC_VEC_BOOL ) {
    return head.back()->get_bool(pos) ? 1 : 0 ;
  }
  return 0 ;
}

EXTERN_C
int
GC_get_int_at_pos( int pos ) {
  if ( ! head.empty() && 
       head.back() != nullptr &&
       head.back() -> get_type() == GC::GC_VEC_INT ) {
    return head.back()->get_int(pos) ;
  }
  return 0 ;
}

EXTERN_C
double
GC_get_real_at_pos( int pos ) {
  if ( ! head.empty() && 
       head.back() != nullptr &&
       head.back() -> get_type() == GC::GC_VEC_REAL ) {
    return head.back()->get_real(pos) ;
  }
  return 0 ;
}

EXTERN_C
char const *
GC_get_string_at_pos( int pos ) {
  if ( ! head.empty() && 
       head.back() != nullptr &&
       head.back() -> get_type() == GC::GC_VEC_STRING ) {
    return head.back()->get_string(pos).c_str() ;
  }
  return 0 ;
}

EXTERN_C
int
GC_set_map() {
  int ok = check_no_data( GC::GC_MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK ) head.back() -> set_map() ;
  return ok ;
}

/*! \brief
 *  Set the position of insertion point is at the `pos` element
 *  of the actual map.
 */
EXTERN_C
int
GC_set_map_position( char const pos[] ) {
  int ok = check( GC::GC_MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    std::string tmp(pos) ; // temporary to make happy clang compiler
    GC::GenericContainer * gc = &((*head.back())[tmp]) ;
    head.push_back(gc) ;
  }
  return ok ;
}

EXTERN_C
int
GC_get_map() {
  int ok = check_no_data( GC::GC_MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK ) {
    ptr_map = &(head.back() -> get_map()) ;
    map_iterator = ptr_map->begin() ;
  }
  return ok ;
}

EXTERN_C
char const *
GC_get_key() {
  int ok = check_no_data( GC::GC_MAP ) ;
  if ( ok == GENERIC_CONTAINER_OK && map_iterator != ptr_map->end() ) {
    return map_iterator++ -> first.c_str() ;
  }
  return nullptr ;
}

//
// eof: GenericContainerCinterface.cc
//
