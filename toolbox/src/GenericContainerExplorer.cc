/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2021                                                      |
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
// file: GenericContainerExplorer.cc
//

#include "GenericContainer/GenericContainer.hh"
#include <cstring>

namespace GC_namespace {

  GenericContainer *
  GenericContainerExplorer::top() {
    GC_ASSERT(
      head.size() > 0,
      "GenericContainerExplorer::top() empty stack!"
    )
    GC_ASSERT(
      head.back() != nullptr,
      "GenericContainerExplorer::top() bad top pointer!"
    )
    return head.back();
  }

  GenericContainer const *
  GenericContainerExplorer::top() const {
    GC_ASSERT(
      head.size() > 0,
      "GenericContainerExplorer::top() empty stack!"
    )
    GC_ASSERT(
      head.back() != nullptr,
      "GenericContainerExplorer::top() bad top pointer!"
    )
    return head.back();
  }

  int
  GenericContainerExplorer::check( GC_type data_type ) const {
    if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD;
    if ( data_type == head.back()->get_type() ) {
      if ( head.back() == nullptr ) return GENERIC_CONTAINER_NO_DATA;
      return GENERIC_CONTAINER_OK;
    } else {
      return GENERIC_CONTAINER_BAD_TYPE;
    }
  }

  int
  GenericContainerExplorer::check_no_data( GC_type data_type ) const {
    if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD;
    if ( GC_type::NOTYPE == head.back()->get_type() ||
         data_type       == head.back()->get_type() ) return GENERIC_CONTAINER_OK;
    return GENERIC_CONTAINER_NOT_EMPTY;
  }

  int
  GenericContainerExplorer::pop() {
    if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA;
    head.pop_back();
    return GENERIC_CONTAINER_OK;
  }

  int
  GenericContainerExplorer::push( GenericContainer * gc ) {
    head.push_back( gc );
    return GENERIC_CONTAINER_OK;
  }

  int
  GenericContainerExplorer::push_vector_position( unsigned pos ) {
    int ok = check( GC_type::VECTOR );
    if ( ok == GENERIC_CONTAINER_OK  ) {
      GenericContainer * gc = &((*head.back())[pos]);
      head.push_back( gc );
    }
    return ok;
  }

  int
  GenericContainerExplorer::push_map_position( char const pos[] ) {
    int ok = check( GC_type::MAP );
    if ( ok == GENERIC_CONTAINER_OK  ) {
      GenericContainer * gc = &((*head.back())[pos]);
      head.push_back( gc );
    }
    return ok;
  }

  int
  GenericContainerExplorer::init_map_key() {
    int ok = check( GC_type::MAP );
    if ( ok == GENERIC_CONTAINER_OK ) {
      ptr_map = &head.back()->get_map();
      map_iterator = ptr_map->begin();
    }
    return ok;
  }

  char const *
  GenericContainerExplorer::next_map_key() {
    if ( map_iterator != ptr_map->end() )
      return map_iterator++->first.c_str();
    else
      return nullptr;
  }

  int
  GenericContainerExplorer::reset() {
    if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA;
    while ( head.size() > 1 ) head.pop_back();
    return GENERIC_CONTAINER_OK;
  }

}

//
// eof: GenericContainerExplorer.cc
//
