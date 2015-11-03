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

#define EXTERN_C extern "C" GENERIC_CONTAINER_API_DLL 

std::map<std::string,GC::GenericContainerExplorer> gc_explorer ;

GC::GenericContainerExplorer * gc_active = nullptr ;

EXTERN_C
int
GC_new( char const id[] ) {
  // ckeck if exists
  gc_active = &gc_explorer[id] ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_select( char const id[] ) {
  // ckeck if exists
  gc_active = &gc_explorer[id] ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
void *
GC_mem_ptr( char const id[] ) {
  GC_select( id ) ;
  if ( gc_active == nullptr ) return nullptr ;
  return gc_active->mem_ptr() ;
}

EXTERN_C
int
GC_delete( char const id[] ) {
  gc_explorer.erase(id) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_fill_for_test( char const id[] ) {
  // ckeck if exists
  gc_active = &gc_explorer[id] ;

  GC::GenericContainer & gc = *(GC::GenericContainer*)gc_active->mem_ptr() ;

  GC::vector_type & v = gc.set_vector() ;
  v.resize(11) ;
  v[0] = 1 ;
  v[1].set_vec_real() ;
  v[2].set_map() ;
  v[3].set_vec_string() ;
  v[4] = 1.3 ;
  v[5] = "pippo" ;
  v[6].set_map() ;
  v[7].set_vector() ;
  v[10] = true ;
  GC::vec_real_type & vv = v[1].get_vec_real() ;
  vv.resize(10) ;
  vv[2] = 123 ;
  GC::map_type & mm = v[2].get_map() ;
  mm["pippo"]    = 13 ;
  mm["pluto"]    = 1  ;
  mm["paperino"] = 3  ;
  GC::GenericContainer & gmm = v[2] ; // access element 2 as GenericContainer
  gmm["aaa"]     = "stringa1"  ; // is the same as mm["aaa"] = "stringa"
  gmm["bbb"]     = "stringa2"  ; // is the same as mm["aaa"] = "stringa"
  GC::vec_string_type & vs = v[3].get_vec_string() ;
  vs.push_back("string1");
  vs.push_back("string2");
  vs.push_back("string3");
  vs.push_back("string4");
  GC::map_type & m = v[6].get_map() ;
  m["aaa"]    = 123 ;
  m["bbb"]    = 3.4 ;
  m["vector"].set_vec_int() ;
  GC::vec_int_type & vi = m["vector"].get_vec_int() ;
  vi.push_back(12) ;
  vi.push_back(10) ;
  vi.push_back(1) ;

  GC::vector_type & vg = v[7].get_vector() ;
  vg.resize(4) ;
  vg[0] = 123 ;
  vg[1] = 3.14 ;
  vg[2] = "nonna papera" ;
  GC::vec_complex_type & vg1 = vg[3].set_vec_complex() ;
  vg1.push_back(1) ;
  vg1.push_back(2) ;
  vg1.push_back(3) ;
  vg1.push_back(std::complex<double>(-12,2)) ;

  GC::GenericContainer & gm = v[8] ;
  gm.set_mat_real(2,2) ;
  gm.get_real(1,1) = 2 ;
  gm.get_real(0,1) = 3 ;

  GC::GenericContainer & gm1 = v[9] ;
  gm1.set_mat_complex(2,2) ;
  gm1.get_complex(1,1) = std::complex<double>(2,2) ;
  gm1.get_complex(0,1) = std::complex<double>(1,-1) ;

  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_pop_head() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  return gc_active->pop() ;
}

EXTERN_C
int
GC_reset_head() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  return gc_active->reset() ;
}

EXTERN_C
int
GC_print() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> print(cout) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_get_type() {
  if ( gc_active == nullptr ) return -1 ;
  return gc_active -> top() -> get_type() ;
}

EXTERN_C
char const *
GC_get_type_name() {
  static char const empty[] = "" ;
  if ( gc_active == nullptr ) return empty ;
  return gc_active -> top() -> get_type_name() ;
}

//..............................................................................

EXTERN_C
int
GC_set_bool( int const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_bool(a?true:false) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_int( int_type const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_int(a) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_real( real_type const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_real(a) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_complex( complex_type const * a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_complex(a->real,a->imag) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_complex2( real_type const re, real_type const im ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_complex(re,im) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_string( char const a[] ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_string(a) ;
  return GENERIC_CONTAINER_OK ;
}

//..............................................................................

EXTERN_C
int
GC_get_bool() {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_bool() ? 1 : 0 ;
}

EXTERN_C
int_type
GC_get_int() {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_int() ;
}

EXTERN_C
long_type
GC_get_long() {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_long() ;
}

EXTERN_C
real_type
GC_get_real( ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_real() ;
}

EXTERN_C
complex_type
GC_get_complex( ) {
  complex_type tmp ;
  tmp.real = tmp.imag = 0 ;
  if ( gc_active != nullptr ) tmp = *((complex_type*)&gc_active -> top() -> get_complex()) ;
  return tmp ;
}

EXTERN_C
real_type
GC_get_complex_re( ) {
  if ( gc_active != nullptr ) return gc_active -> top() -> get_complex() . real() ;
  return 0 ;
}

EXTERN_C
real_type
GC_get_complex_im( ) {
  if ( gc_active != nullptr ) return gc_active -> top() -> get_complex() . imag() ;
  return 0 ;
}

EXTERN_C
char const *
GC_get_string( ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_string() . c_str() ;
}

//..............................................................................

EXTERN_C
int
GC_push_bool( int const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_bool( a != 0 ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_push_int( int_type const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_int( a ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_push_real( real_type const a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_real( a ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_push_complex( complex_type const * a ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_complex( a->real, a->imag ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_push_complex2( real_type const re, real_type const im ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_complex( re, im ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_push_string( char const a[] ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> push_string( a ) ;
  return GENERIC_CONTAINER_OK ;
}

//..............................................................................

EXTERN_C
int
GC_get_bool_at_pos( int pos ) {
  if ( gc_active == nullptr ) return 0 ;
  bool val = gc_active -> top() -> get_bool( unsigned(pos) ) ;
  return val ? 1 : 0 ;
}

EXTERN_C
int_type
GC_get_int_at_pos( int pos ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_int( unsigned(pos) ) ;
}

EXTERN_C
real_type
GC_get_real_at_pos( int pos ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_number( unsigned(pos) ) ;
}

EXTERN_C
complex_type
GC_get_complex_at_pos( int pos ) {
  complex_type tmp ;
  tmp.real = tmp.imag = 0 ;
  if ( gc_active != nullptr ) gc_active -> top() -> get_complex_number( unsigned(pos), tmp.real, tmp.imag ) ;
  return tmp ;
}

EXTERN_C
real_type
GC_get_complex_real_at_pos( int pos ) {
  real_type re = 0, im = 0 ;
  if ( gc_active != nullptr ) gc_active -> top() -> get_complex_number( unsigned(pos), re, im ) ;
  return re ;
}

EXTERN_C
real_type
GC_get_complex_imag_at_pos( int pos ) {
  real_type re = 0, im = 0 ;
  if ( gc_active != nullptr ) gc_active -> top() -> get_complex_number( unsigned(pos), re, im ) ;
  return im ;
}

EXTERN_C
char const *
GC_get_string_at_pos( int pos ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_string( unsigned(pos) ).c_str() ;
}

//..............................................................................

EXTERN_C
real_type
GC_get_real_at_coor( int i, int j ) {
  if ( gc_active == nullptr ) return 0 ;
  return gc_active -> top() -> get_real( unsigned(i), unsigned(j) ) ;
}

EXTERN_C
complex_type
GC_get_complex_at_coor( int i, int j ) {
  complex_type tmp ;
  tmp.real = tmp.imag = 0 ;
  if ( gc_active != nullptr ) {
    GC::complex_type res = gc_active -> top() -> get_complex( unsigned(i), unsigned(j) ) ;
    tmp.real = res.real() ;
    tmp.imag = res.imag() ;
  }
  return tmp ;
}

EXTERN_C
real_type
GC_get_complex_real_at_coor( int i, int j ) {
  if ( gc_active != nullptr ) {
    GC::complex_type res = gc_active -> top() -> get_complex( unsigned(i), unsigned(j) ) ;
    return res.real() ;
  }
  return 0 ;
}

EXTERN_C
real_type
GC_get_complex_imag_at_coor( int i, int j ) {
  if ( gc_active != nullptr ) {
    GC::complex_type res = gc_active -> top() -> get_complex( unsigned(i), unsigned(j) ) ;
    return res.imag() ;
  }
  return 0 ;
}

//..............................................................................

EXTERN_C
int
GC_set_empty_vector_of_bool() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vec_bool() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_empty_vector_of_int() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vec_int() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_empty_vector_of_real() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vec_real() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_empty_vector_of_complex() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vec_complex() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_empty_vector_of_string() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vec_string() ;
  return GENERIC_CONTAINER_OK ;
}

//..............................................................................

EXTERN_C
int
GC_set_vector_of_bool( int const a[], int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  GC::vec_bool_type & v = gc_active -> top() -> set_vec_bool( unsigned(nelem) ) ;
  for ( unsigned i = 0 ; i < unsigned(nelem) ; ++i ) v[i] = a[i] != 0 ;
  return GENERIC_CONTAINER_OK ;
}

//..............................................................................

EXTERN_C
int
GC_set_vector_of_int( int_type const a[], int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  GC::vec_int_type & v = gc_active -> top() -> set_vec_int( unsigned(nelem) ) ;
  std::copy( a, a + nelem, v.begin() ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_vector_of_real( real_type const a[], int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  GC::vec_real_type & v = gc_active -> top() -> set_vec_real( unsigned(nelem) ) ;
  std::copy( a, a + nelem, v.begin() ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_vector_of_complex( real_type const re[], real_type const im[], int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  GC::vec_complex_type & v = gc_active -> top() -> set_vec_complex( unsigned(nelem) ) ;
  for ( unsigned i = 0 ; i < unsigned(nelem) ; ++i ) v[i] = GC::complex_type(re[i],im[i]) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_vector_of_string( char const *a[], int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  GC::vec_string_type & v = gc_active -> top() -> set_vec_string( unsigned(nelem) ) ;
  for ( unsigned i = 0 ; i < unsigned(nelem) ; ++i ) v[i] = a[i] ;
  return GENERIC_CONTAINER_OK ;
}

//..............................................................................

EXTERN_C
int
GC_set_vector( int nelem ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vector( unsigned(nelem) ) ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_set_empty_vector() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active -> top() -> set_vector() ;
  return GENERIC_CONTAINER_OK ;
}

EXTERN_C
int
GC_get_vector_size() {
  if ( gc_active == nullptr ) return 0 ;
  return int(gc_active -> top() -> get_num_elements()) ;
}

EXTERN_C
int
GC_get_matrix_num_rows() {
  if ( gc_active == nullptr ) return 0 ;
  return int(gc_active -> top() -> get_numRows()) ;
}

EXTERN_C
int
GC_get_matrix_num_cols() {
  if ( gc_active == nullptr ) return 0 ;
  return int(gc_active -> top() -> get_numCols()) ;
}

EXTERN_C
int
GC_push_vector_position( int pos ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  return gc_active -> push_vector_position( unsigned(pos) ) ;
}

//..............................................................................

EXTERN_C
int
GC_set_map() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  gc_active->top()->set_map() ;
  return GENERIC_CONTAINER_OK ;
}

/*! \brief
 *  Set the position of insertion point is at the `pos` element
 *  of the actual map.
 */
EXTERN_C
int
GC_push_map_position( char const pos[] ) {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  return gc_active -> push_map_position( pos ) ;
}

EXTERN_C
int
GC_init_map_key() {
  if ( gc_active == nullptr ) return GENERIC_CONTAINER_BAD_HEAD ;
  int ok = gc_active -> init_map_key() ;
  return ok ;
}

EXTERN_C
char const *
GC_get_next_key() {
  if ( gc_active == nullptr ) return nullptr ;
  return gc_active -> next_map_key() ;
}

//
// eof: GenericContainerCinterface.cc
//
