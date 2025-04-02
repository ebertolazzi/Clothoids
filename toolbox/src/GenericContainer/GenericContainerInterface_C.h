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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once

/*
// file: GenericContainerCinterface.h
*/

//!
//! \addtogroup GC
//!
//! @{

#ifndef GENERIC_CONTAINER_INTERFACE_C_H
#define GENERIC_CONTAINER_INTERFACE_C_H

#ifdef __cplusplus
extern "C" {
#endif

//! \enum GC_ErrorCodes
//! \brief Error codes returned by the GenericContainer functions.
enum {
  GENERIC_CONTAINER_OK = 0,
  GENERIC_CONTAINER_BAD_TYPE,
  GENERIC_CONTAINER_NO_DATA,
  GENERIC_CONTAINER_NOT_EMPTY,
  GENERIC_CONTAINER_BAD_HEAD
};

//! \struct c_complex_type
//! \brief A structure representing a complex number, with real and imaginary components.
typedef struct {
  double real;  //!< Real part of the complex number
  double imag;  //!< Imaginary part of the complex number
} c_complex_type;

//!
//! \brief Create a new GenericContainer object.
//!
//! \param[in] id A string identifying the new GenericContainer object.
//! \return Error code: 0 = OK.
//!
int GC_new( char const id[] );

//!
//! \brief Select an existing GenericContainer object.
//!
//! \param[in] id A string identifying the GenericContainer to select.
//! \return Error code: 0 = OK.
//!
int GC_select( char const id[] );

//!
//! \brief Delete a GenericContainer object.
//!
//! \param[in] id A string identifying the GenericContainer to delete.
//! \return Error code: 0 = OK.
//!
int GC_delete( char const id[] );

//!
//! \brief Fill a GenericContainer with test data.
//!
//! \param[in] id A string identifying the GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_fill_for_test( char const id[] );

//!
//! \brief Move the head pointer up one level in the GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_pop_head();

//!
//! \brief Reset the head pointer to the first level in the GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_reset_head();

//!
//! \brief Dump the contents of the current GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_dump();

//!
//! \brief Print the content types of the current GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_print_content_types();

//!
//! \brief Retrieve the type of the current element pointed by the head in the GenericContainer.
//! \return Error code: 0 = OK.
//!
int GC_get_type();

//!
//! \brief Retrieve the name of the type of the current element in the GenericContainer.
//! \return The type name as a string.
//!
char const * GC_get_type_name();

//!
//! \brief Get a pointer to the internal memory of the GenericContainer object.
//!
//! \param[in] id A string identifying the GenericContainer.
//! \return Pointer to the internal memory of the object.
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
//! \brief Set the current element of the `GenericContainer` to an integer value.
//!
//! This function assigns the integer value `a` to the currently selected element
//! of the `GenericContainer`. If the element is not of the correct type,
//! an error code will be returned.
//!
//! \param[in] a Integer value to be stored.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_int( int a );

//!
//! \brief Set the current element of the `GenericContainer` to a real number.
//!
//! This function assigns the double-precision floating-point value `a` to the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type, an error code will be returned.
//!
//! \param[in] a Double value (real number) to be stored.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_real( double a );

//!
//! \brief Set the current element of the `GenericContainer` to a complex number.
//!
//! This function stores a complex number represented by the structure
//! `c_complex_type` in the currently selected element of the `GenericContainer`.
//! If the element is not of the correct type, an error code will be returned.
//!
//! \param[in] a Pointer to a `c_complex_type` structure containing the real and
//!              imaginary parts of the complex number.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_complex( c_complex_type const * a );

//!
//! \brief Set the current element of the `GenericContainer` to a complex number using real and imaginary parts.
//!
//! This function stores a complex number directly using its real (`re`) and
//! imaginary (`im`) parts in the currently selected element of the
//! `GenericContainer`. If the element is not of the correct type, an error code
//! will be returned.
//!
//! \param[in] re Real part of the complex number.
//! \param[in] im Imaginary part of the complex number.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_complex2( double re, double im );

//!
//! \brief Set the current element of the `GenericContainer` to a string.
//!
//! This function stores a string `a` in the currently selected element of the
//! `GenericContainer`. If the element is not of the correct type, an error code
//! will be returned.
//!
//! \param[in] a C-string to be stored (null-terminated).
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_string( char const a[] );

// -----------------------------------------------------------------------------

//!
//! \brief Retrieve the boolean value from the current element of the `GenericContainer`.
//!
//! This function returns the boolean value (stored as an integer) from the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type, an error code or unexpected value may be returned.
//!
//! \return The boolean value (as integer, 0 = false, 1 = true).
//!
int GC_get_bool();

//!
//! \brief Retrieve the integer value from the current element of the `GenericContainer`.
//!
//! This function returns the integer value from the currently selected element
//! of the `GenericContainer`. If the element is not of the correct type, an
//! error code or unexpected value may be returned.
//!
//! \return The integer value.
//!
int GC_get_int();

//!
//! \brief Retrieve the long integer value from the current element of the `GenericContainer`.
//!
//! This function returns the long integer value from the currently selected
//! element of the `GenericContainer`. If the element is not of the correct type,
//! an error code or unexpected value may be returned.
//!
//! \return The long integer value.
//!
long GC_get_long();

//!
//! \brief Retrieve the real number value from the current element of the `GenericContainer`.
//!
//! This function returns the floating-point value from the currently selected
//! element of the `GenericContainer`. If the element is not of the correct type,
//! an error code or unexpected value may be returned.
//!
//! \return The floating-point value.
//!
double GC_get_real();

//!
//! \brief Retrieve the complex number from the current element of the `GenericContainer`.
//!
//! This function returns the complex number (as a `c_complex_type` structure)
//! from the currently selected element of the `GenericContainer`. If the element
//! is not of the correct type, an error code or unexpected value may be returned.
//!
//! \return The complex number value as `c_complex_type`.
//!
c_complex_type GC_get_complex();

//!
//! \brief Retrieve the real part of a complex number from the current element of the `GenericContainer`.
//!
//! This function returns the real part of a complex number stored in the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type, an error code or unexpected value may be returned.
//!
//! \return The real part of the complex number.
//!
double GC_get_complex_re();

//!
//! \brief Retrieve the imaginary part of a complex number from the current element of the `GenericContainer`.
//!
//! This function returns the imaginary part of a complex number stored in the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type, an error code or unexpected value may be returned.
//!
//! \return The imaginary part of the complex number.
//!
double GC_get_complex_im();

//!
//! \brief Retrieve the string from the current element of the `GenericContainer`.
//!
//! This function returns the C-string stored in the currently selected element
//! of the `GenericContainer`. If the element is not of the correct type, an
//! error code or unexpected value may be returned.
//!
//! \return The C-string (null-terminated).
//!
char const * GC_get_string();

// -----------------------------------------------------------------------------

//!
//! \brief Append a boolean value to the currently selected vector of booleans in the `GenericContainer`.
//!
//! This function pushes a boolean value (stored as an integer) to the end of
//! the vector of booleans in the currently selected element of the
//! `GenericContainer`. If the element is not of the correct type or is not a
//! vector, an error code will be returned.
//!
//! \param[in] a Boolean value (stored as an integer) to be appended.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_bool( int a );

//!
//! \brief Append an integer value to the currently selected vector of integers in the `GenericContainer`.
//!
//! This function pushes an integer value `a` to the end of the vector of
//! integers in the currently selected element of the `GenericContainer`. If the
//! element is not of the correct type or is not a vector, an error code will be
//! returned.
//!
//! \param[in] a Integer value to be appended.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_int( int a );

//!
//! \brief Append a floating-point value to the currently selected vector of real numbers in the `GenericContainer`.
//!
//! This function pushes a real number `a` (double precision) to the end of the
//! vector of real numbers in the currently selected element of the
//! `GenericContainer`. If the element is not of the correct type or is not a
//! vector, an error code will be returned.
//!
//! \param[in] a Real number (double precision) to be appended.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_real( double a );

//!
//! \brief Append a complex number to the currently selected vector of complex numbers in the `GenericContainer`.
//!
//! This function pushes a complex number `a` to the end of the vector of complex
//! numbers in the currently selected element of the `GenericContainer`. If the
//! element is not of the correct type or is not a vector, an error code will be
//! returned.
//!
//! \param[in] a Pointer to a `c_complex_type` structure containing the complex number.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_complex( c_complex_type const * a );

//!
//! \brief Append a complex number to the currently selected vector using real and imaginary parts.
//!
//! This function pushes a complex number, specified by its real (`re`) and
//! imaginary (`im`) parts, to the end of the vector of complex numbers in the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type or is not a vector, an error code will be returned.
//!
//! \param[in] re Real part of the complex number.
//! \param[in] im Imaginary part of the complex number.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_complex2( double re, double im );

//!
//! \brief Append a string to the currently selected vector of strings in the `GenericContainer`.
//!
//! This function pushes a string `a` to the end of the vector of strings in the
//! currently selected element of the `GenericContainer`. If the element is not
//! of the correct type or is not a vector, an error code will be returned.
//!
//! \param[in] a C-string (null-terminated) to be appended.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_string( char const a[] );

// -----------------------------------------------------------------------------

//!
//! \brief Retrieve a boolean value at a specified position from a vector of booleans in the `GenericContainer`.
//!
//! This function retrieves the boolean value (stored as an integer) at position
//! `pos` from a vector of booleans in the currently selected element of the
//! `GenericContainer`. If the element is not a vector or if `pos` is out of
//! range, an error code or unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The boolean value at the specified position.
//!
int GC_get_bool_at_pos( int pos );

//!
//! \brief Retrieve an integer value at a specified position from a vector of integers in the `GenericContainer`.
//!
//! This function retrieves the integer value at position `pos` from a vector of
//! integers in the currently selected element of the `GenericContainer`. If the
//! element is not a vector or if `pos` is out of range, an error code or
//! unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The integer value at the specified position.
//!
int GC_get_int_at_pos( int pos );

//!
//! \brief Retrieve a real number at a specified position from a vector of real numbers in the `GenericContainer`.
//!
//! This function retrieves the floating-point value at position `pos` from a
//! vector of real numbers in the currently selected element of the
//! `GenericContainer`. If the element is not a vector or if `pos` is out of
//! range, an error code or unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The floating-point value at the specified position.
//!
double GC_get_real_at_pos( int pos );

//!
//! \brief Retrieve a complex number at a specified position from a vector of complex numbers in the `GenericContainer`.
//!
//! This function retrieves the complex number (as a `c_complex_type` structure)
//! at position `pos` from a vector of complex numbers in the currently selected
//! element of the `GenericContainer`. If the element is not a vector or if `pos`
//! is out of range, an error code or unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The complex number at the specified position.
//!
c_complex_type GC_get_complex_at_pos( int pos );

//!
//! \brief Retrieve the real part of a complex number at a specified position from a vector of complex numbers.
//!
//! This function retrieves the real part of a complex number at position `pos`
//! from a vector of complex numbers in the currently selected element of the
//! `GenericContainer`. If the element is not a vector or if `pos` is out of
//! range, an error code or unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The real part of the complex number at the specified position.

//!
double GC_get_complex_real_at_pos( int pos );

//!
//! \brief Retrieve the imaginary part of a complex number at a specified position from a vector of complex numbers.
//!
//! This function retrieves the imaginary part of a complex number at position
//! `pos` from a vector of complex numbers in the currently selected element of
//! the `GenericContainer`. If the element is not a vector or if `pos` is out of
//! range, an error code or unexpected value may be returned.
//!
//! \param[in] pos Index position of the value to be retrieved.
//! \return The imaginary part of the complex number at the specified position.
//!
double GC_get_complex_imag_at_pos( int pos );

//!
//! \brief Retrieve a string at a specified position from a vector of strings in the `GenericContainer`.
//!
//! This function retrieves the C-string stored at position `pos` from a vector
//! of strings in the currently selected element of the `GenericContainer`. If
//! the element is not a vector or if `pos` is out of range, an error code or
//! unexpected value may be returned.
//!
//! \param[in] pos Index position of the string to be retrieved.
//! \return The C-string (null-terminated) at the specified position.
//!
char const * GC_get_string_at_pos( int pos );

// -----------------------------------------------------------------------------

//!
//! \brief Retrieve a real number at the specified row and column from a matrix of real numbers in the `GenericContainer`.
//!
//! This function retrieves a floating-point value at the matrix coordinates
//! (row `i`, column `j`) from a matrix of real numbers in the currently selected
//! element of the `GenericContainer`. If the element is not a matrix or if the
//! coordinates are out of range, an error code or unexpected value may be
//! returned.
//!
//! \param[in] i Row index of the value to be retrieved.
//! \param[in] j Column index of the value to be retrieved.
//! \return The floating-point value at the specified coordinates.
//!
double GC_get_real_at_coor( int i, int j );

//!
//! \brief Retrieve a complex number at the specified row and column from a matrix of complex numbers.
//!
//! This function retrieves the complex number (as a `c_complex_type` structure)
//! at the matrix coordinates (row `i`, column `j`) from a matrix of complex
//! numbers in the currently selected element of the `GenericContainer`. If the
//! element is not a matrix or if the coordinates are out of range, an error code
//! or unexpected value may be returned.
//!
//! \param[in] i Row index of the value to be retrieved.
//! \param[in] j Column index of the value to be retrieved.
//! \return The complex number at the specified coordinates.
//!
c_complex_type GC_get_complex_at_coor( int i, int j );

//!
//! \brief Retrieve the real part of a complex number at the specified row and column from a matrix.
//!
//! This function retrieves the real part of a complex number at the matrix
//! coordinates (row `i`, column `j`) from a matrix of complex numbers in the
//! currently selected element of the `GenericContainer`. If the element is not
//! a matrix or if the coordinates are out of range, an error code or unexpected
//! value may be returned.
//!
//! \param[in] i Row index of the value to be retrieved.
//! \param[in] j Column index of the value to be retrieved.
//! \return The real part of the complex number at the specified coordinates.
//!
double GC_get_complex_real_at_coor( int i, int j );

//!
//! \brief Retrieve the imaginary part of a complex number at the specified row and column from a matrix.
//!
//! This function retrieves the imaginary part of a complex number at the matrix
//! coordinates (row `i`, column `j`) from a matrix of complex numbers in the
//! currently selected element of the `GenericContainer`. If the element is not
//! a matrix or if the coordinates are out of range, an error code or unexpected
//! value may be returned.
//!
//! \param[in] i Row index of the value to be retrieved.
//! \param[in] j Column index of the value to be retrieved.
//! \return The imaginary part of the complex number at the specified coordinates.
//!
double GC_get_complex_imag_at_coor( int i, int j );

// -----------------------------------------------------------------------------

//!
//! \brief Initialize the current element of the `GenericContainer` to an empty vector of booleans.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty vector of booleans (size 0).
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector_of_bool();

//!
//! \brief Initialize the current element of the `GenericContainer` to an empty vector of integers.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty vector of integers (size 0).
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector_of_int();

//!
//! \brief Initialize the current element of the `GenericContainer` to an empty vector of real numbers.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty vector of real numbers (size 0).
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector_of_real();

//!
//! \brief Initialize the current element of the `GenericContainer` to an empty vector of complex numbers.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty vector of complex numbers (size 0).
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector_of_complex();

//!
//! \brief Initialize the current element of the `GenericContainer` to an empty vector of strings.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty vector of strings (size 0).
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector_of_string();

// -----------------------------------------------------------------------------

//!
//! \brief Set the current element of the `GenericContainer` to a vector of booleans, and copy the provided values.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to a vector of booleans with `nelem` elements, copying the boolean values
//! from the provided array `a`.
//!
//! \param[in] a     Array of boolean values (stored as integers) to be copied.
//! \param[in] nelem Number of elements in the array `a`.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector_of_bool( int const * a, int nelem );

//!
//! \brief Set the current element of the `GenericContainer` to a vector of integers, and copy the provided values.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to a vector of integers with `nelem` elements, copying the integer values
//! from the provided array `a`.
//!
//! \param[in] a     Array of integer values to be copied.
//! \param[in] nelem Number of elements in the array `a`.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector_of_int( int const * a, int nelem );

//!
//! \brief Set the current element of the `GenericContainer` to a vector of real numbers, and copy the provided values.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to a vector of real numbers (`double`) with `nelem` elements, copying the
//! values from the provided array `a`.
//!
//! \param[in] a     Array of floating-point values to be copied.
//! \param[in] nelem Number of elements in the array `a`.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector_of_real( double const * a, int nelem );

//!
//! \brief Set the current element of the `GenericContainer` to a vector of complex numbers, and copy the provided values.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to a vector of complex numbers with `nelem` elements, copying the real and
//! imaginary parts of each complex number from the provided arrays `re` and `im`.
//!
//! \param[in] re    Array of real parts of complex numbers to be copied.
//! \param[in] im    Array of imaginary parts of complex numbers to be copied.
//! \param[in] nelem Number of elements in the arrays `re` and `im`.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector_of_complex( double const * re, double const * im, int nelem );

//!
//! \brief Set the current element of the `GenericContainer` to a vector of strings, and copy the provided values.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to a vector of strings with `nelem` elements, copying the strings from the
//! provided array `a`.
//!
//! \param[in] a     Array of C-strings (null-terminated) to be copied.
//! \param[in] nelem Number of elements in the array `a`.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector_of_string( char const ** a, int nelem );

// -----------------------------------------------------------------------------

//!
//! \brief Set the current element of the `GenericContainer` to a vector of `GenericContainer` elements with the specified size.
//!
//! This function initializes the currently selected element of the
//! `GenericContainer` to a vector of `GenericContainer` elements, with `nelem`
//! elements. The insertion point is set to the first element.
//!
//! \param[in] nelem Number of `GenericContainer` elements to be created.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_vector( int nelem );

//!
//! \brief Set the current element of the `GenericContainer` to an empty vector of `GenericContainer` elements.
//!
//! This function initializes the currently selected element of the
//! `GenericContainer` to an empty vector of `GenericContainer` elements.
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_empty_vector();

//!
//! \brief Retrieve the number of elements in the current vector of the `GenericContainer`.
//!
//! This function returns the number of elements in the currently selected vector
//! of the `GenericContainer`. If the element is not a vector, an error code or
//! unexpected value may be returned.
//!
//! \return The number of elements in the current vector, or an error code.
//!
int GC_get_vector_size();

//!
//! \brief Retrieve the number of rows in the current matrix of the `GenericContainer`.
//!
//! This function returns the number of rows in the currently selected matrix of
//! the `GenericContainer`. If the element is not a matrix, an error code or
//! unexpected value may be returned.
//!
//! \return The number of rows in the current matrix, or an error code.
//!
int GC_get_matrix_num_rows();

//!
//! \brief Retrieve the number of columns in the current matrix of the `GenericContainer`.
//!
//! This function returns the number of columns in the currently selected matrix
//! of the `GenericContainer`. If the element is not a matrix, an error code or
//! unexpected value may be returned.
//!
//! \return The number of columns in the current matrix, or an error code.
//!
int GC_get_matrix_num_cols();

//!
//! \brief Set the insertion point in the current vector of the `GenericContainer`.
//!
//! This function sets the insertion point to the element at position `pos` in
//! the currently selected vector of the `GenericContainer`. If the element is
//! not a vector or if `pos` is out of range, an error code will be returned.
//!
//! \param[in] pos Index position where the insertion point is set.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_vector_position( int pos );

// -----------------------------------------------------------------------------

//!
//! \brief Initialize the current element of the `GenericContainer` to a map.
//!
//! This function sets the currently selected element of the `GenericContainer`
//! to an empty map. If the element is not of the correct type, an error code
//! will be returned.
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_set_map();

//!
//! \brief Initialize the current element as a map key in the `GenericContainer`.
//!
//! This function sets the current element of the `GenericContainer` as a key
//! for a map. If the element is not of the correct type, an error code will be
//! returned.
//!
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_init_map_key();

//!
//! \brief Retrieve the key of the next element in a map from the `GenericContainer`.
//!
//! This function returns the key of the next element in a map within the
//! `GenericContainer`. If the element is not a map or there are no more keys,
//! an error code or unexpected value may be returned.
//!
//! \return The key as a C-string (null-terminated), or an error code.
//!
char const * GC_get_next_key( void );

//!
//! \brief Set the insertion point to the specified key in the current map of the `GenericContainer`.
//!
//! This function sets the insertion point to the map element associated with
//! the specified key `pos`. If the element is not a map or the key does not
//! exist, an error code will be returned.
//!
//! \param[in] pos C-string representing the key.
//! \return Error code: 0 = OK, other values indicate failure.
//!
int GC_push_map_position( char const pos[] );

#ifdef __cplusplus
}
#endif

#endif

//!
//! @}
//!

/*
// eof: GenericContainerCinterface.hh
*/
