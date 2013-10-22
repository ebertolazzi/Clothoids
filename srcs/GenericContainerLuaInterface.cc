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
 |      E.Bertolazzi, F.Biral, P.Bosetti                                    |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "GenericContainerLuaInterface.hh"

#ifndef ASSERT
  #include <sstream>
  #include <stdexcept>
  #define ASSERT(COND,MSG)                               \
    if ( !(COND) ) {                                     \
      std::ostringstream ost ;                           \
      ost << "in GenericContainer: " << MSG << '\n' ;    \
      GenericContainer::exception( ost.str().c_str() ) ; \
    }
#endif

#ifdef DEBUG
  #ifndef ASSERT_DEBUG
    #define ASSERT_DEBUG(COND,MSG) ASSERT(COND,MSG)
  #endif
#else
  #ifndef ASSERT_DEBUG
    #define ASSERT_DEBUG(COND,MSG)
  #endif
#endif

#include <lua.hpp>

using namespace std ;

typedef double valueType ;
typedef int    indexType ;

static
void
lua_to_GC( lua_State        * L,
           GenericContainer & gc,
           string const     & indent ) ;

static
void
push_vec_element( lua_State        * L,
                  GenericContainer & gc,
                  string const     & indent ) {
  // assegna il valore
  int idx  = lua_tonumber(L, -2)-1 ; // index start from 1 in LUA
  int type = lua_type(L, -1) ;
  switch( type ) {
    case LUA_TBOOLEAN:
      {
        gc.get_bool(idx) ;
        GenericContainer::vec_bool_type & bv = gc.get_vec_bool() ;
        bv[idx] = lua_toboolean(L, -1) ;
      }
      break ;
    case LUA_TNUMBER:
      {
        valueType val = lua_tonumber(L, -1) ;
        if ( gc.get_type() == GenericContainer::GC_VEC_REAL ) {
          gc.get_real(idx) = val ;
        } else if ( indexType(val) == val ) {
          gc.get_int(idx)  = indexType(val) ;
        } else {
          gc.get_real(idx) = val ;
        }
      }
      break ;
    case LUA_TSTRING:
      gc.get_string(idx) = lua_tostring(L, -1) ;
      break ;
    case LUA_TTABLE:
      lua_to_GC( L, gc[idx], indent+"  " ) ;
      break ;
  }
}

static
void
push_hash_element( lua_State        * L,
                   GenericContainer & gc,
                   string const     & indent ) {
  // assegna il valore
  string key  = lua_tostring(L, -2) ;
  int    type = lua_type(L, -1) ;
  switch( type ) {
    case LUA_TBOOLEAN:
      gc[key].set_bool(lua_toboolean(L, -1)) ;
      break ;
    case LUA_TNUMBER:
      {
        valueType  val = lua_tonumber(L, -1) ;
        indexType ival = indexType(val) ;
        if ( ival == val ) gc[key].set_int(ival) ;
        else               gc[key].set_real(val) ;
      }
      break ;
    case LUA_TSTRING:
      gc[key].set_string(lua_tostring(L, -1)) ;
      break ;
    case LUA_TTABLE:
      lua_to_GC( L, gc[key], indent+"  " ) ;
      break ;
  }
}

static
void
lua_to_GC( lua_State        * L,
           GenericContainer & gc,
           string const     & indent ) {
  // ---
  lua_pushnil(L) ; // first key
  // ---
  while ( lua_next(L,-2) != 0 ) {
    switch( lua_type(L, -2) ) {
      case LUA_TNUMBER:
        push_vec_element( L, gc, indent ) ;
        break ;
      case LUA_TSTRING:
        push_hash_element( L, gc, indent ) ;
        break ;
    }
    lua_pop(L, 1);  // removes `value'; keeps `key' for next iteration
  }
}

LuaInterpreter::LuaInterpreter() {
  lua_State *& L = *((lua_State**)&void_L) ;
  L = luaL_newstate() ; // opens Lua
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::openSession() sessionLuaState invalid!" ) ;
  luaL_openlibs(L) ;
}

LuaInterpreter::~LuaInterpreter() {
  lua_State *& L = *((lua_State**)&void_L) ;
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::~LuaInterpreter() lua_State invalid!" ) ;
  lua_close(L);
}
  
void
LuaInterpreter::load( char const filename[] ) {
  lua_State *& L = *((lua_State**)&void_L) ;
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::~LuaInterpreter() lua_State invalid!" ) ;
  if ( luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)  ) {
    ASSERT( lua_isnil(L,-1),
            "In LuaInterpreter::load('" << filename << "')\n" << lua_tostring(L, -1) ) ;
  }
}
  
void
LuaInterpreter::parseString( char const cmd[] ) {
  lua_State *& L = *((lua_State**)&void_L) ;
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::~LuaInterpreter() lua_State invalid!" ) ;
  if ( luaL_loadbuffer(L, cmd, strlen(cmd), "line") || lua_pcall(L, 0, 0, 0) )
    ASSERT( false, "In Lua::parseString('" << cmd << "')\ncannot run the command" ) ;
}

void
LuaInterpreter::global_to_GC( char const global_var[], GenericContainer & gc ) {
  lua_State *& L = *((lua_State**)&void_L) ;
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::global_to_GC(...) lua_State invalid!" ) ;
    
  lua_getglobal( L, global_var ) ;
  ASSERT( !lua_isnil(L,-1),
          "LuaInterpreter::global_to_GC(...) cannot find global variable: '" << global_var << "'" ) ;
  ASSERT( lua_istable(L,-1),
          "LuaInterpreter::global_to_GC(...) global variable '" << global_var << "' is not a table" ) ;
    
  gc.clear() ;
  global_to_GC( gc ) ;
  lua_settop(L, 0);
}

void
LuaInterpreter::global_to_GC( GenericContainer & gc ) {
  lua_State *& L = *((lua_State**)&void_L) ;
  ASSERT_DEBUG( L != NULL, "LuaInterpreter::global_to_GC(...) sessionLuaState invalid!" ) ;
  lua_to_GC( L, gc, "" ) ;
}

void
LuaInterpreter::toLua( GenericContainer const & gc, char const var[] ) {
  // creo variabile
}
