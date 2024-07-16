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

//!
//! \addtogroup Generic Container C interface
//!

/* @{ */

#ifndef GENERIC_CONTAINER_C_INTERFACE_H
#define GENERIC_CONTAINER_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
enum {
  GENERIC_CONTAINER_OK = 0,
  GENERIC_CONTAINER_BAD_TYPE,
  GENERIC_CONTAINER_NO_DATA,
  GENERIC_CONTAINER_NOT_EMPTY,
  GENERIC_CONTAINER_BAD_HEAD
};
#endif

typedef struct {
  double real;
  double imag;
} c_complex_type;

//!
//! Create a new `GenericContainer` object 'id'
//!
//! \param[in] id string `id` of the new `GenericContainer`
//! \return error code, 0 = OK
//!
int GC_new( char const id[] );

//!
//! Select an old `GenericContainer` object 'id'
//!
//! \param[in] id string `id` of the `GenericContainer`
//! \return error code, 0 = OK
//!
int GC_select( char const id[] );

//!
//! Delete the `GenericContainer` object 'id'
//!
//! \param[in] id string `id` of the `GenericContainer`
//! \return error code, 0 = OK
//!
int GC_delete( char const id[] );

//!
//! Fill the `GenericContainer` object 'id' with data for test purposes
//!
//! \param[in] id string `id` of the `GenericContainer`
//! \return error code, 0 = OK
//!
int GC_fill_for_test( char const id[] );

//!
//! Move `head` up to a level
//!
//! \return error code, 0 = OK
//!
int GC_pop_head();

//!
//! Move `head` to the first level
//!
//! \return error code, 0 = OK
//!
int GC_reset_head();

//!
//! Print the actual `GenericContainer`
//!
//! \return error code, 0 = OK
//!
int GC_dump();

//!
//! Print the actual `GenericContainer`
//!
//! \return error code, 0 = OK
//!
int GC_print_content_types();

//!
//! Get type of actual pointed element of `GenericContainer`
//!
//! \return error code, 0 = OK
//!
int GC_get_type();

//!
//! Get type of actual pointed element of `GenericContainer`
//!
//! \return id name of the `GenericContainer`
//!
char const * GC_get_type_name();

//!
//! Get pointer to the internal `GenericContainer` object 'id'
//!
//! \param[in] id string `id` of the `GenericContainer`
//! \return pointer to the internal object
//!
void * GC_mem_ptr( char const id[] );

// -----------------------------------------------------------------------------

//!
//! Set actual pointed element of `GenericContainer` to `bool` with value `a`
//!
//! \param[in] a boolean to be stored
//! \return error code, 0 = OK
//!
int GC_set_bool( int a );

//!
//! Set actual pointed element of `GenericContainer` to `int` with value `a`
//!
//! \param[in] a integer to be stored
//! \return error code, 0 = OK
//!
int GC_set_int( int a );

//!
//! Set actual pointed element of `GenericContainer` to `real_type` with value `a`
//!
//! \param[in] a double to be stored
//! \return error code, 0 = OK
//!
int GC_set_real( double a );

//!
//! Set actual pointed element of `GenericContainer` to `complex` with value `a`
//!
//! \param[in] a complex number to be stored
//! \return error code, 0 = OK
//!
int GC_set_complex( c_complex_type const * a );

//!
//! Set actual pointed element of `GenericContainer` to `complex` with value `a`
//!
//! \param[in] re complex number to be stored, real part
//! \param[in] im complex number to be stored, imaginary part
//! \return error code, 0 = OK
//!
int GC_set_complex2( double re, double im );

//!
//! Set actual pointed element of `GenericContainer` to `string` with value `a`
//!
//! \param[in] a string to be stored
//! \return error code, 0 = OK
//!
int GC_set_string( char const a[] );

// -----------------------------------------------------------------------------

//!
//! Get actual pointed element of `GenericContainer` of type `bool`
//!
//! \return the value
//!
int GC_get_bool();

//!
//! Get actual pointed element of `GenericContainer` of type `int`
//!
//! \return the value
//!
int GC_get_int();

//!
//! Get actual pointed element of `GenericContainer` of type `int`
//!
//! \return the value
//!
long GC_get_long();

//!
//! Get actual pointed element of `GenericContainer` of type `real_type`
//!
//! \return the value
//!
double GC_get_real();

//!
//! Get actual pointed element of `GenericContainer` of type `complex`
//!
//! \return the value
//!
c_complex_type GC_get_complex();

//!
//! Get actual pointed element of `GenericContainer` of type `complex`
//!
//! \return the value
//!
double GC_get_complex_re();

//!
//! Get actual pointed element of `GenericContainer` of type `complex`
//!
//! \return the value
//!
double GC_get_complex_im();

//!
//! Get actual pointed element of `GenericContainer` of type `string`
//!
//! \return the value
//!
char const * GC_get_string();

// -----------------------------------------------------------------------------

//!
//! Push `a` to a vector of `bool`
//!
//! \param[in] a the value to be stored
//! \return error code, 0 = OK
//!
int GC_push_bool( int a );

//!
//! Push `a` to a vector of `integer`
//!
//! \param[in] a the value to be stored
//! \return error code, 0 = OK
//!
int GC_push_int( int a );

//!
//! Push `a` to a vector of `real_type`
//!
//! \param[in] a the value to be stored
//! \return error code, 0 = OK
//!
int GC_push_real( double a );

//!
//! Push `a` to a vector of `complex`
//!
//! \param[in] a the value to be stored
//! \return error code, 0 = OK
//!
int GC_push_complex( c_complex_type const * a );

//!
//! Push `a` to a vector of `complex`
//!
//! \param[in] re the value to be stored, real part
//! \param[in] im the value to be stored, imaginary part
//! \return error code, 0 = OK
//!
int GC_push_complex2( double re, double im );

//!
//! Push `a` to a vector of string
//!
//! \param[in] a the value to be stored
//! \return error code, 0 = OK
//!
int GC_push_string( char const a[] );

// -----------------------------------------------------------------------------

//!
//! Get boolean at position `pos` of a vector of bool
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
int GC_get_bool_at_pos( int pos );

//!
//! Get `int_type` at position `pos` of a vector of `int_type`
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
int GC_get_int_at_pos( int pos );

//!
//! Get `real_type` at position `pos` of a vector of `real_type`
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
double GC_get_real_at_pos( int pos );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
c_complex_type GC_get_complex_at_pos( int pos );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
double GC_get_complex_real_at_pos( int pos );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
double GC_get_complex_imag_at_pos( int pos );

//!
//! Get string at position `pos` of a vector of string
//!
//! \param[in] pos position of the value to be extracted
//! \return the stored value
//!
char const * GC_get_string_at_pos( int pos );

// -----------------------------------------------------------------------------

//!
//! Get `real_type` at position `i,j` of a matrix of `real_type`
//!
//! \param[in] i row
//! \param[in] j column
//! \return the stored value
//!
double GC_get_real_at_coor( int i, int j );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] i row
//! \param[in] j column
//! \return the stored value
//!
c_complex_type GC_get_complex_at_coor( int i, int j );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] i row
//! \param[in] j column
//! \return the stored value
//!
double GC_get_complex_real_at_coor( int i, int j );

//!
//! Get `complex` at position `pos` of a vector of `real_type`
//!
//! \param[in] i row
//! \param[in] j column
//! \return the stored value
//!
double GC_get_complex_imag_at_coor( int i, int j );

// -----------------------------------------------------------------------------

//!
//! Set actual pointed element of `GenericContainer` to a vector of `bool` of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector_of_bool();

//!
//! Set actual pointed element of `GenericContainer` to a vector of `integer` of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector_of_int();

//!
//! Set actual pointed element of `GenericContainer` to a vector of `real_type` of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector_of_real();

//!
//! Set actual pointed element of `GenericContainer` to a vector of `complex` of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector_of_complex();

//!
//! Set actual pointed element of `GenericContainer` to a vector of string of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector_of_string();

// -----------------------------------------------------------------------------

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of bool of size `nelem` and copy vector
//! of int `a` to the element.
//! \param a     vector of bool (stored as integer) to be stored
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector_of_bool( int const * a, int nelem );

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of integer of size `nelem` and copy vector
//! of int `a` to the element.
//!
//! \param a     vector of integer to be stored
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector_of_int( int const * a, int nelem );

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of `real_type` of size `nelem` and copy vector
//! of real_type `a` to the element.
//!
//! \param a     vector of double to be stored
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector_of_real( double const * a, int nelem );

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of `complex` of size `nelem` and copy vector `a`
//! of `real_type` to the element.
//!
//! \param re    vector of real part to be stored
//! \param im    vector of imaginary to be stored
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector_of_complex( double const * re, double const * im, int nelem );

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of strings of size `nelem` and copy vector
//! of strings `a` to the element.
//!
//! \param a     vector of string to be stored
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector_of_string( char const ** a, int nelem );

// -----------------------------------------------------------------------------

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of `GenericContainer` of size `nelem`.
//! The position of insertion point is at the first element.
//!
//! \param nelem number of element of the vector
//! \return error code, 0 = OK
//!
int GC_set_vector( int nelem );

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a vector of `GenericContainer` of size `0`.
//!
//! \return error code, 0 = OK
//!
int GC_set_empty_vector();

//!
//! \brief
//! Get the size of the actual vector.
//!
//! \return error code, 0 = OK
//!
int GC_get_vector_size();

//!
//! \brief
//! Get the size of the actual matrix.
//!
//! \return error code, 0 = OK
//!
int GC_get_matrix_num_rows();

//!
//! \brief
//! Get the size of the actual matrix.
//!
//! \return error code, 0 = OK
//!
int GC_get_matrix_num_cols();

//!
//! \brief
//! Set the position of insertion point is at the `pos` element
//! of the actual generic vector.
//!
//! \param[in] pos insertion position
//! \return error code, 0 = OK
//!
int GC_push_vector_position( int pos );

// -----------------------------------------------------------------------------

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a map of `GenericContainer`.
//!
//! \return error code, 0 = OK
//!
int GC_set_map();

//!
//! \brief
//! Set actual pointed element of `GenericContainer` to
//! a map of `GenericContainer`.
//!
//! \return error code, 0 = OK
//!
int GC_init_map_key();

//!
//! \brief
//! Return key of the actual element of a map
//!
//! \return error code, 0 = OK
//!
char const * GC_get_next_key( void );

//!
//! \brief
//! Set the position of insertion point is at the `pos` element
//! of the actual map.
//!
//! \param[in] pos key of the map
//! \return error code, 0 = OK
//!
int GC_push_map_position( char const pos[] );

#ifdef __cplusplus
}
#endif

#endif

/* @} */

/*
// eof: GenericContainerCinterface.hh
*/
