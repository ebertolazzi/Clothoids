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

/*
// file: GenericContainerCinterface.h
*/

/*!
  \addtogroup Generic Container C interface
 */

/* @{ */

#ifndef GENERIC_CONTAINER_C_INTERFACE_H
#define GENERIC_CONTAINER_C_INTERFACE_H

#ifndef GENERIC_CONTAINER_API_DLL
  #if defined(_WIN32) || defined(_WIN64)
    #ifdef GENERIC_CONTAINER_EXPORT
      #define GENERIC_CONTAINER_API_DLL __declspec(dllexport)
    #elif defined(GENERIC_CONTAINER_IMPORT)
     #define GENERIC_CONTAINER_API_DLL __declspec(dllimport)
    #else
      #define GENERIC_CONTAINER_API_DLL
    #endif
  #else
    #define GENERIC_CONTAINER_API_DLL
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum {
  GENERIC_CONTAINER_OK = 0,
  GENERIC_CONTAINER_BAD_TYPE,
  GENERIC_CONTAINER_NO_DATA,
  GENERIC_CONTAINER_NOT_EMPTY,
  GENERIC_CONTAINER_BAD_HEAD
} ;

typedef int    int_type ;
typedef long   long_type ;
typedef double real_type ;

typedef struct {
  real_type real ;
  real_type imag ;
} complex_type ;

/*! Create a new `GenericContainer` object 'id' */
int GENERIC_CONTAINER_API_DLL GC_new( char const id[] ) ;

/*! Select an old `GenericContainer` object 'id' */
int GENERIC_CONTAINER_API_DLL GC_select( char const id[] ) ;

/*! Delete the `GenericContainer` object 'id' */
int GENERIC_CONTAINER_API_DLL GC_delete( char const id[] ) ;

/*! Fill the `GenericContainer` object 'id' with data for test purposes */
int GENERIC_CONTAINER_API_DLL GC_fill_for_test( char const id[] ) ;

/*! Move `head` up to a level */
int GENERIC_CONTAINER_API_DLL GC_pop_head() ;

/*! Move `head` to the first level */
int GENERIC_CONTAINER_API_DLL GC_reset_head() ;

/*! Print the actual `GenericContainer` */
int GENERIC_CONTAINER_API_DLL GC_print() ;

/*! Get type of actual pointed element of `GenericContainer` */
int GENERIC_CONTAINER_API_DLL GC_get_type() ;
  
/*! Get type of actual pointed element of `GenericContainer` */
char const * GENERIC_CONTAINER_API_DLL GC_get_type_name() ;
  
/*! Get pointer to the internal `GenericContainer` object 'id' */
void * GENERIC_CONTAINER_API_DLL GC_mem_ptr( char const id[] ) ;

// -----------------------------------------------------------------------------

/*! Set actual pointed element of `GenericContainer` to `bool` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_bool( int const a ) ;

/*! Set actual pointed element of `GenericContainer` to `int` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_int( int_type const a ) ;

/*! Set actual pointed element of `GenericContainer` to `real_type` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_real( real_type const a ) ;

/*! Set actual pointed element of `GenericContainer` to `complex` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_complex( complex_type const * a ) ;

/*! Set actual pointed element of `GenericContainer` to `complex` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_complex2( real_type const re, real_type const im ) ;

/*! Set actual pointed element of `GenericContainer` to `string` with value `a` */
int GENERIC_CONTAINER_API_DLL GC_set_string( char const a[] ) ;

// -----------------------------------------------------------------------------

/*! Get actual pointed element of `GenericContainer` of type `bool` */
int GENERIC_CONTAINER_API_DLL GC_get_bool() ;

/*! Get actual pointed element of `GenericContainer` of type `int` */
int_type GENERIC_CONTAINER_API_DLL GC_get_int() ;

/*! Get actual pointed element of `GenericContainer` of type `int` */
long_type GENERIC_CONTAINER_API_DLL GC_get_long() ;

/*! Get actual pointed element of `GenericContainer` of type `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_real() ;

/*! Get actual pointed element of `GenericContainer` of type `complex` */
complex_type GENERIC_CONTAINER_API_DLL GC_get_complex() ;

/*! Get actual pointed element of `GenericContainer` of type `complex` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_re() ;

/*! Get actual pointed element of `GenericContainer` of type `complex` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_im() ;

/*! Get actual pointed element of `GenericContainer` of type `string` */
char const * GENERIC_CONTAINER_API_DLL GC_get_string() ;

// -----------------------------------------------------------------------------

/*! Push `a` to a vector of `bool` */
int GENERIC_CONTAINER_API_DLL GC_push_bool( int const a ) ;

/*! Push `a` to a vector of `integer` */
int GENERIC_CONTAINER_API_DLL GC_push_int( int_type const a ) ;

/*! Push `a` to a vector of `real_type` */
int GENERIC_CONTAINER_API_DLL GC_push_real( real_type const a ) ;

/*! Push `a` to a vector of `complex`  */
int GENERIC_CONTAINER_API_DLL GC_push_complex( complex_type const * a ) ;
  
/*! Push `a` to a vector of `complex`  */
int GENERIC_CONTAINER_API_DLL GC_push_complex2( real_type const re, real_type const im ) ;

/*! Push `a` to a vector of string */
int GENERIC_CONTAINER_API_DLL GC_push_string( char const a[] ) ;

// -----------------------------------------------------------------------------

/*! Get boolean at position `pos` of a vector of bool  */
int GENERIC_CONTAINER_API_DLL GC_get_bool_at_pos( int pos ) ;

/*! Get `int_type` at position `pos` of a vector of `int_type` */
int_type GENERIC_CONTAINER_API_DLL GC_get_int_at_pos( int pos ) ;

/*! Get `real_type` at position `pos` of a vector of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_real_at_pos( int pos ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
complex_type GENERIC_CONTAINER_API_DLL GC_get_complex_at_pos( int pos ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_real_at_pos( int pos ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_imag_at_pos( int pos ) ;

/*! Get string at position `pos` of a vector of string */
char const * GENERIC_CONTAINER_API_DLL GC_get_string_at_pos( int pos ) ;

// -----------------------------------------------------------------------------

/*! Get `real_type` at position `i,j` of a matrix of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_real_at_coor( int i, int j ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
complex_type GENERIC_CONTAINER_API_DLL GC_get_complex_at_coor( int i, int j ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_real_at_coor( int i, int j ) ;

/*! Get `complex` at position `pos` of a vector of `real_type` */
real_type GENERIC_CONTAINER_API_DLL GC_get_complex_imag_at_coor( int i, int j ) ;

// -----------------------------------------------------------------------------

/*! Set actual pointed element of `GenericContainer` to a vector of `bool` of size `0`. */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector_of_bool() ;

/*! Set actual pointed element of `GenericContainer` to a vector of `integer` of size `0`. */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector_of_int() ;

/*! Set actual pointed element of `GenericContainer` to a vector of `real_type` of size `0`. */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector_of_real() ;

/*! Set actual pointed element of `GenericContainer` to a vector of `complex` of size `0`. */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector_of_complex() ;
  
/*! Set actual pointed element of `GenericContainer` to a vector of string of size `0`. */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector_of_string() ;

// -----------------------------------------------------------------------------

/*! \brief
 *  Set actual pointed element of `GenericContainer` to 
 *  a vector of bool of size `nelem` and copy vector 
 *  of int `a` to the element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector_of_bool( int const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to 
 *  a vector of integer of size `nelem` and copy vector 
 *  of int `a` to the element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector_of_int( int_type const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `real_type` of size `nelem` and copy vector
 *  of real_type `a` to the element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector_of_real( real_type const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `complex` of size `nelem` and copy vector `a`
 *  of `real_type` to the element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector_of_complex( real_type const re[], real_type const im[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of strings of size `nelem` and copy vector
 *  of strings `a` to the element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector_of_string( char const *a[], int nelem ) ;

// -----------------------------------------------------------------------------

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `GenericContainer` of size `nelem`.
 *  The position of insertion point is at the first element.
 */
int GENERIC_CONTAINER_API_DLL GC_set_vector( int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `GenericContainer` of size `0`.
 */
int GENERIC_CONTAINER_API_DLL GC_set_empty_vector() ;

/*! \brief
 *  Get the size of the actual vector.
 */
int GENERIC_CONTAINER_API_DLL GC_get_vector_size() ;

/*! \brief
 *  Get the size of the actual matrix.
 */
int GENERIC_CONTAINER_API_DLL GC_get_matrix_num_rows() ;

/*! \brief
 *  Get the size of the actual matrix.
 */
int GENERIC_CONTAINER_API_DLL GC_get_matrix_num_cols() ;

/*! \brief
 *  Set the position of insertion point is at the `pos` element 
 *  of the actual generic vector.
 */
int GENERIC_CONTAINER_API_DLL GC_push_vector_position( int pos ) ;

// -----------------------------------------------------------------------------

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a map of `GenericContainer`.
 */
int GENERIC_CONTAINER_API_DLL GC_set_map() ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a map of `GenericContainer`.
 */
int GENERIC_CONTAINER_API_DLL GC_init_map_key() ;

/*! \brief
 *  Return key of the actual element of a map 
 */
char const * GENERIC_CONTAINER_API_DLL GC_get_next_key() ;

/*! \brief
 *  Set the position of insertion point is at the `pos` element
 *  of the actual map.
 */
int GENERIC_CONTAINER_API_DLL GC_push_map_position( char const pos[] ) ;

#ifdef __cplusplus
}
#endif

#endif

/* @} */

/*
// eof: GenericContainerCinterface.hh
*/
