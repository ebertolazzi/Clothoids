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

#ifndef GENERIC_CONTAINER_LUA_INTERFACE_HH
#define GENERIC_CONTAINER_LUA_INTERFACE_HH

#include "GenericContainer.hh"

namespace GenericContainerNamespace {

  class GENERIC_CONTAINER_API_DLL LuaInterpreter {
    /* lua_State * */ void * void_L ; //!< interpreter status
  public:
    LuaInterpreter() ;
    ~LuaInterpreter() ;
    void dump( std::basic_ostream<char> & stream ) ;
    void execute( char const cmd[] ) ;
    void call( GenericContainer const & args, GenericContainer & res ) ;
    void do_file( char const fname[] ) ;
    void GC_to_global( GenericContainer const & gc, char const [] ) ; // not yet implemented
    void global_to_GC( char const var[], GenericContainer & gc ) ;
    int  interactive( int argc,
                      char const * argv[],
                      char const * messages[],
                      char const * prompt ) ; // launch interpret mode
  } ;
}

#endif
