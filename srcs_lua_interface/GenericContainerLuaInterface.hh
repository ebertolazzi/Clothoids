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
public:
  LuaInterpreter() ;
  ~LuaInterpreter() ;
  void dump( std::ostream & stream ) ;
  void execute( char const cmd[] ) ;
  void do_file( char const fname[], bool check_syntax_only = false ) ;
  void GC_to_global( GenericContainer const & gc, char const [] ) ; // not yet implemented
  void global_to_GC( char const var[], GenericContainer & gc ) ;
} ;

#endif
