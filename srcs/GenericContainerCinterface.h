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

/*! Select or create a new `GenericContainer` object 'name' */
int GC_select( char const name[] ) ;

/*! Delete the `GenericContainer` object 'name' */
int GC_delete( char const name[] ) ;

/*! Move `head` up to a level */
int GC_pop_head() ;

/*! Print the actual `GenericContainer` */
int GC_print() ;
  
/*! Get pointer to the internal `GenericContainer` object 'name' */
void * GC_mem_ptr( char const name[] ) ;
  
/*! Set actual pointed element of `GenericContainer` to `bool` with value `a` */
int GC_set_bool( int const a ) ;

/*! Set actual pointed element of `GenericContainer` to `int` with value `a` */
int GC_set_int( int const a ) ;

/*! Set actual pointed element of `GenericContainer` to `double` with value `a` */
int GC_set_real( double const a ) ;

/*! Set actual pointed element of `GenericContainer` to `string` with value `a` */
int GC_set_string( char const a[] ) ;

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
 *  Set the position of insertion point is at the `pos` element 
 *  of the actual generic vector.
 */
int GC_set_vector_position( int pos ) ;
  
/*! \brief
 *  Set actual pointed element of `GenericContainer` to
 *  a map of `GenericContainer`.
 */
int GC_set_map() ;
  
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
