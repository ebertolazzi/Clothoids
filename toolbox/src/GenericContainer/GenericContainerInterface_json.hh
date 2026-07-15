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

#ifndef GENERIC_CONTAINER_INTERFACE_JSON_HH
#define GENERIC_CONTAINER_INTERFACE_JSON_HH

#include "GenericContainer.hh"

#include <fstream>

namespace GC_namespace
{

  using std::ifstream;

  inline bool file_JSON_to_GC( string_view file_name, GenericContainer & gc )
  {
    ifstream stream( file_name.data() );
    gc.clear();
    return gc.from_json( stream );
  }

  inline bool JSON_to_GC( istream_type & stream, GenericContainer & gc )
  {
    gc.clear();
    return gc.from_json( stream );
  }

  inline bool JSON_to_GC( string const & data, GenericContainer & gc )
  {
    gc.clear();
    return gc.from_json( data );
  }

  inline bool JSON_to_GC( vec_string_type const & chunks, GenericContainer & gc )
  {
    string_type json;
    for ( auto const & chunk : chunks ) json += chunk;
    gc.clear();
    return gc.from_json( json );
  }

  inline void GC_to_JSON( GenericContainer const & gc, std::string & res )
  { res = gc.to_json(); }

  inline void GC_to_JSON( GenericContainer const & gc, ostream_type & stream )
  { gc.to_json( stream ); }

  inline void GC_to_JSON( GenericContainer const & gc, vec_string_type & chunks )
  { chunks.assign( 1, gc.to_json() ); }

}  // namespace GC_namespace

#endif
