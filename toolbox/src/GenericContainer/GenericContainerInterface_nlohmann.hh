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

#ifndef GENERIC_CONTAINER_INTERFACE_NLOHMANN_HH
#define GENERIC_CONTAINER_INTERFACE_NLOHMANN_HH

#include "GenericContainer.hh"

#include <nlohmann/json.hpp>

#include <cstdint>
#include <limits>
#include <sstream>
#include <utility>

//!
//! Opt-in `nlohmann::json` support for `GenericContainer`.
//!
//! Including this header (and only this header) enables conversion between
//! `GC_namespace::GenericContainer` and `nlohmann::json` via the standard
//! `adl_serializer` customization point:
//!
//! \code
//! GC_namespace::GenericContainer gc; gc["a"] = 42;
//! nlohmann::json j = gc;                                            // to_json
//! GC_namespace::GenericContainer gc2 = j.get<GC_namespace::GenericContainer>(); // from_json
//! \endcode
//!
//! `nlohmann::json` is header-only, and so is this adapter: `GenericContainer`
//! itself has no hard dependency on it, only translation units that include
//! this header need `nlohmann/json.hpp` on their include path.
//!
namespace nlohmann
{

  //!
  //! \addtogroup JSON
  //!
  //! @{

  template <>
  struct adl_serializer<GC_namespace::GenericContainer>
  {
    using GenericContainer = GC_namespace::GenericContainer;
    using GC_type          = GC_namespace::GC_type;

    static void
    to_json( json & j, GenericContainer const & gc )
    {
      switch ( gc.get_type() )
      {
        case GC_type::NOTYPE: j = nullptr; break;
        case GC_type::BOOL: j = gc.get_bool(); break;
        case GC_type::INTEGER: j = gc.get_int(); break;
        case GC_type::LONG: j = gc.get_long(); break;
        case GC_type::REAL: j = gc.get_real(); break;
        case GC_type::COMPLEX: j = GC_namespace::to_string( gc.get_complex() ); break;
        case GC_type::STRING: j = gc.get_string(); break;
        case GC_type::POINTER: j = pointer_to_hex( gc.get_pvoid() ); break;
        case GC_type::VEC_BOOL: j = gc.get_vec_bool(); break;
        case GC_type::VEC_INTEGER: j = gc.get_vec_int(); break;
        case GC_type::VEC_LONG: j = gc.get_vec_long(); break;
        case GC_type::VEC_REAL: j = gc.get_vec_real(); break;
        case GC_type::VEC_STRING: j = gc.get_vec_string(); break;
        case GC_type::VEC_COMPLEX:
        {
          json arr = json::array();
          for ( auto const & vi : gc.get_vec_complex() ) arr.push_back( GC_namespace::to_string( vi ) );
          j = std::move( arr );
          break;
        }
        case GC_type::VEC_POINTER:
        {
          json arr = json::array();
          for ( auto const & vi : gc.get_vec_pointer() ) arr.push_back( pointer_to_hex( vi ) );
          j = std::move( arr );
          break;
        }
        case GC_type::VECTOR:
        {
          json arr = json::array();
          for ( auto const & vi : gc.get_vector() )
          {
            json jj;
            to_json( jj, vi );
            arr.push_back( std::move( jj ) );
          }
          j = std::move( arr );
          break;
        }
        case GC_type::MAP:
        {
          json obj = json::object();
          for ( auto const & [key, val] : gc.get_map() )
          {
            json jj;
            to_json( jj, val );
            obj[key] = std::move( jj );
          }
          j = std::move( obj );
          break;
        }
        case GC_type::MAT_INTEGER: j = matrix_to_json( gc.get_mat_int() ); break;
        case GC_type::MAT_LONG: j = matrix_to_json( gc.get_mat_long() ); break;
        case GC_type::MAT_REAL: j = matrix_to_json( gc.get_mat_real() ); break;
        case GC_type::MAT_COMPLEX:
        {
          auto const & M{ gc.get_mat_complex() };
          json         arr = json::array();
          std::size_t const NR{ M.num_rows() };
          std::size_t const NC{ M.num_cols() };
          for ( std::size_t jc{ 0 }; jc < NC; ++jc )
          {
            json col = json::array();
            for ( std::size_t ir{ 0 }; ir < NR; ++ir ) col.push_back( GC_namespace::to_string( M( ir, jc ) ) );
            arr.push_back( std::move( col ) );
          }
          j = std::move( arr );
          break;
        }
      }
    }

    static void
    from_json( json const & j, GenericContainer & gc )
    {
      gc.clear();
      from_json_value( j, gc );
      // Collapse homogeneous vectors of scalars (and, when possible, vectors of
      // matching-length vectors) into vec_int_type/vec_real_type/... and matrices,
      // mirroring GenericContainer::from_json()'s own JSON path. Without this, an
      // array like [1,2,3] would decode into a generic vector<GenericContainer> of
      // three INTEGER-typed elements, forcing callers such as copyto_vec_real() to
      // promote element-by-element instead of getting the fast, expected
      // vec_real_type/vec_int_type representation directly.
      gc.collapse();
    }

  private:
    template <typename Mat>
    static json
    matrix_to_json( Mat const & M )
    {
      json           arr = json::array();
      std::size_t const NR{ M.num_rows() };
      std::size_t const NC{ M.num_cols() };
      for ( std::size_t jc{ 0 }; jc < NC; ++jc )
      {
        json col = json::array();
        for ( std::size_t ir{ 0 }; ir < NR; ++ir ) col.push_back( M( ir, jc ) );
        arr.push_back( std::move( col ) );
      }
      return arr;
    }

    static std::string
    pointer_to_hex( void const * p )
    {
      std::ostringstream ss;
      ss << std::hex << std::showbase << reinterpret_cast<std::uintptr_t>( p );
      return ss.str();
    }

    static void
    from_json_value( json const & j, GenericContainer & gc )
    {
      if ( j.is_null() ) { gc.clear(); }
      else if ( j.is_boolean() ) { gc = j.get<bool>(); }
      else if ( j.is_number_integer() )
      {
        if ( j.is_number_unsigned() )
        {
          auto const v{ j.get<std::uint64_t>() };
          if ( v <= static_cast<std::uint64_t>( std::numeric_limits<GenericContainer::uint_type>::max() ) )
            gc = static_cast<GenericContainer::uint_type>( v );
          else
            gc = static_cast<GenericContainer::ulong_type>( v );
        }
        else
        {
          auto const v{ j.get<std::int64_t>() };
          if ( v >= std::numeric_limits<GenericContainer::int_type>::min() &&
               v <= std::numeric_limits<GenericContainer::int_type>::max() )
            gc = static_cast<GenericContainer::int_type>( v );
          else
            gc = static_cast<GenericContainer::long_type>( v );
        }
      }
      else if ( j.is_number_float() ) { gc = j.get<GenericContainer::real_type>(); }
      else if ( j.is_string() ) { gc = j.get<std::string>(); }
      else if ( j.is_array() )
      {
        auto & V{ gc.set_vector( j.size() ) };
        for ( std::size_t i{ 0 }; i < j.size(); ++i ) from_json_value( j[i], V[i] );
      }
      else if ( j.is_object() )
      {
        auto & M{ gc.set_map() };
        for ( auto it{ j.begin() }; it != j.end(); ++it ) from_json_value( it.value(), M[it.key()] );
      }
    }
  };

  //!
  //! @}
  //!

}  // namespace nlohmann

#endif

//
// eof: GenericContainerInterface_nlohmann.hh
//
