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

#ifndef GENERIC_CONTAINER_FFI_H
#define GENERIC_CONTAINER_FFI_H

#include <stdint.h>

//!
//! \addtogroup GC
//!
//! @{
//!
//! A minimal, ABI-stable, pure-C wrapper around `GenericContainer`, meant as
//! a base for writing FFI bindings in other languages (Python `ctypes`/`cffi`,
//! Ruby `FFI`, Julia `ccall`, Rust `bindgen`/`extern "C"`, Go `cgo`, ...).
//! It is a single header/implementation pair, exporting only C-compatible
//! types (opaque handle, fixed-width integers, `double`, `char*`) -- no C++
//! type ever crosses this boundary, and no C++ exception ever escapes it.
//!
//! Design:
//!  - Lifecycle: gc_new() / gc_free().
//!  - A handful of scalar get/set functions cover the common case where the
//!    whole container is just one value.
//!  - Every function that can fail clears the thread-local last-error
//!    message on entry and sets it if it catches a C++ exception; check
//!    gc_last_error() right after a call if its return value alone (0,
//!    0.0, NULL, an empty string) cannot distinguish success from failure.
//!  - Every function accepts a NULL handle safely (sets the last error,
//!    returns a sentinel) instead of crashing the host process.
//!

#ifdef __cplusplus
extern "C" {
#endif

//! Opaque handle to a `GenericContainer` instance.
typedef void * gc_handle_t;

// -----------------------------------------------------------------------------
//  Lifecycle
// -----------------------------------------------------------------------------

//! Create a new, empty `GenericContainer`. Returns NULL on allocation failure.
gc_handle_t gc_new( void );

//! Destroy a `GenericContainer` created by gc_new(). NULL is a safe no-op.
void gc_free( gc_handle_t h );

// -----------------------------------------------------------------------------
//  Scalar convenience accessors
// -----------------------------------------------------------------------------

void gc_set_bool( gc_handle_t h, int v );
void gc_set_int32( gc_handle_t h, int32_t v );
void gc_set_int64( gc_handle_t h, int64_t v );
void gc_set_real( gc_handle_t h, double v );
void gc_set_string( gc_handle_t h, char const * v );

//! \return The stored value, or 0/0.0 on type mismatch or NULL handle.
int gc_get_bool( gc_handle_t h );
//! \return The stored value, or 0 on type mismatch or NULL handle.
int32_t gc_get_int32( gc_handle_t h );
//! \return The stored value, or 0 on type mismatch or NULL handle.
int64_t gc_get_int64( gc_handle_t h );
//! \return The stored value, or 0.0 on type mismatch or NULL handle.
double gc_get_real( gc_handle_t h );
//!
//! \return The stored string, borrowed -- valid until the next call made on
//!         `h`, do not free it. Empty string on type mismatch or NULL handle.
//!
char const * gc_get_string( gc_handle_t h );

// -----------------------------------------------------------------------------
//  Introspection
// -----------------------------------------------------------------------------

//! \return The numeric GC_type of the value stored in `h` (see GC_namespace::GC_type), or -1 on a NULL handle.
int gc_type( gc_handle_t h );

//! \return The name of the GC_type of the value stored in `h` (e.g. "int_type", "map_type"), or "" on a NULL handle.
char const * gc_type_name( gc_handle_t h );

// -----------------------------------------------------------------------------
//  Error handling
// -----------------------------------------------------------------------------

//!
//! \return The last error message recorded on the calling thread by any
//!         gc_* function, or an empty string if the most recent call
//!         succeeded. Never NULL.
//!
char const * gc_last_error( void );

#ifdef __cplusplus
}
#endif

#endif

//!
//! @}
//!

/*
// eof: GenericContainer_ffi.h
*/
