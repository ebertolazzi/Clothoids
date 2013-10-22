/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2000                                                      |
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

#ifndef GENERIC_CONTAINER_LUA_INTERFACE_HH
#define GENERIC_CONTAINER_LUA_INTERFACE_HH

#include "GenericContainer.hh"

class LuaInterpreter {
  /* lua_State * */ void * void_L ; //!< interpreter status
  void fromLua( /* mrb_value */ void * p_v, GenericContainer & gc, std::string const & trace ) ;
  void toLua( GenericContainer const & gc, /* mrb_value */ void * p_v, std::string const & trace ) ;
  void global_to_GC( GenericContainer & gc ) ;
public:
  LuaInterpreter() ;
  ~LuaInterpreter() ;
  void dump( std::ostream & stream ) ;
  void load( char const fname[] ) ;
  void parseString( char const cmd[] ) ;
  void toLua( GenericContainer const & gc, char const [] ) ;
  void global_to_GC( char const var[], GenericContainer & gc ) ;
} ;

#endif
