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

/*! Select or create a new `GenericContainer` object 'id' */
int GC_select( char const id[] ) ;

/*! Create the `GenericContainer` object 'id' */
int GC_new( char const id[] ) ;

/*! Delete the `GenericContainer` object 'id' */
int GC_delete( char const id[] ) ;

/*! Move `head` up to a level */
int GC_pop_head() ;

/*! Move `head` to the first level */
int GC_reset_head() ;

/*! Print the actual `GenericContainer` */
int GC_print() ;

/*! Get type of actual pointed element of `GenericContainer` */
int GC_get_type() ;
  
/*! Get type of actual pointed element of `GenericContainer` */
char const * GC_get_type_name() ;
  
/*! Get pointer to the internal `GenericContainer` object 'id' */
void * GC_mem_ptr( char const id[] ) ;

/*! Set actual pointed element of `GenericContainer` to `bool` with value `a` */
int GC_set_bool( int const a ) ;

/*! Get actual pointed element of `GenericContainer` of type `bool` */
int GC_get_bool() ;

/*! Set actual pointed element of `GenericContainer` to `int` with value `a` */
int GC_set_int( int const a ) ;

/*! Get actual pointed element of `GenericContainer` of type `int` */
int GC_get_int() ;

/*! Set actual pointed element of `GenericContainer` to `double` with value `a` */
int GC_set_real( double const a ) ;

/*! Get actual pointed element of `GenericContainer` of type `double` */
double GC_get_real() ;

/*! Set actual pointed element of `GenericContainer` to `string` with value `a` */
int GC_set_string( char const a[] ) ;

/*! Get actual pointed element of `GenericContainer` of type `string` */
char const * GC_get_string() ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to 
 *  a vector of bool of size `nelem` and copy vector 
 *  of int `a` to the element.
 */
int GC_set_vector_of_bool( int const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `bool` of size `0`.
 */
int GC_set_empty_vector_of_bool() ;

/*! \brief
 *  Push `a` to a vector of `bool`
 */
int GC_push_bool( int const a ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to 
 *  a vector of integer of size `nelem` and copy vector 
 *  of int `a` to the element.
 */
int GC_set_vector_of_int( int const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `integer` of size `0`.
 */
int GC_set_empty_vector_of_int() ;

/*! \brief
 *  Push `a` to a vector of `integer`
 */
int GC_push_int( int const a ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `double` of size `nelem` and copy vector
 *  of double `a` to the element.
 */
int GC_set_vector_of_real( double const a[], int nelem ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `double` of size `0`.
 */
int GC_set_empty_vector_of_real() ;

/*! \brief
 *  Push `a` to a vector of integer
 */
int GC_push_real( double const a ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of strings of size `nelem` and copy vector
 *  of strings `a` to the element.
 */
int GC_set_vector_of_string( char const *a[], unsigned nelem ) ;
  
/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of string of size `0`.
 */
int GC_set_empty_vector_of_string() ;
  
/*! \brief
 *  Push `a` to a vector of string
 */
int GC_push_string( char const a[] ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `GenericContainer` of size `nelem`.
 *  The position of insertion point is at the first element.
 */
int GC_set_vector( int nelem ) ;
  
/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a vector of `GenericContainer` of size `0`.
 */
int GC_set_empty_vector() ;

/*! \brief
 *  Get the size of the actual vector.
 */
int GC_get_vector_size() ;

/*! \brief
 *  Set the position of insertion point is at the `pos` element 
 *  of the actual generic vector.
 */
int GC_set_vector_position( int pos ) ;


/*! \brief
 *  Get boolean at position `pos` of a vector of bool
 */
int GC_get_bool_at_pos( int pos ) ;

/*! \brief
 *  Get integer at position `pos` of a vector of integer
 */
int GC_get_int_at_pos( int pos ) ;

/*! \brief
 *  Get double at position `pos` of a vector of double
 */
double GC_get_real_at_pos( int pos ) ;

/*! \brief
 *  Get string at position `pos` of a vector of string
 */
char const * GC_get_string_at_pos( int pos ) ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a map of `GenericContainer`.
 */
int GC_set_map() ;

/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a map of `GenericContainer`.
 */
int GC_get_map() ;

/*! \brief
 *  Return key of the actual element of a map 
 */
char const * GC_get_key() ;

/*! \brief
 *  Set the position of insertion point is at the `pos` element
 *  of the actual map.
 */
int GC_set_map_position( char const pos[] ) ;

#ifdef __cplusplus
}
#endif

#endif

/* @} */

/*
// eof: GenericContainerCinterface.hh
*/
