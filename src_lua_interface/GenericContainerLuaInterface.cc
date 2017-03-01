/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           \|                   \|                        |
 |                                                                          |
 |      E.Bertolazzi, F.Biral, P.Bosetti                                    |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "GenericContainerLuaInterface.hh"

#ifndef GC_ASSERT
  #include <sstream>
  #include <stdexcept>
  #define GC_ASSERT(COND,MSG)                            \
    if ( !(COND) ) {                                     \
      std::ostringstream ost ;                           \
      ost << "in GenericContainer: " << MSG << '\n' ;    \
      GenericContainer::exception( ost.str().c_str() ) ; \
    }
#endif

#ifdef DEBUG
  #ifndef GC_ASSERT_DEBUG
    #define GC_ASSERT_DEBUG(COND,MSG) GC_ASSERT(COND,MSG)
  #endif
#else
  #ifndef GC_ASSERT_DEBUG
    #define GC_ASSERT_DEBUG(COND,MSG)
  #endif
#endif

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lua.hpp>

// load string.h for strlen
#include <string.h>

using namespace std ;

namespace GenericContainerNamespace {

  static
  inline
  bool isZero( real_type x )
  { return FP_ZERO == fpclassify(x) ; }

  static
  inline
  bool isInteger( real_type x )
  { return isZero( x-static_cast<long>(floor(x)) ) ; }

  /*
  //   _                _           ____  ____
  //  | |_   _  __ _   | |_ ___    / ___|/ ___|
  //  | | | | |/ _` |  | __/ _ \  | |  _| |
  //  | | |_| | (_| |  | || (_) | | |_| | |___
  //  |_|\__,_|\__,_|___\__\___/___\____|\____|
  //               |_____|    |_____|
  */

  static
  void
  lua_table_to_GC( lua_State * L, GenericContainer & gc ) ;

  static
  void
  push_vec_element( lua_State * L, GenericContainer & gc ) {
    // assegna il valore
    // index start from 1 in LUA
    unsigned    idx  = unsigned(lua_tointeger(L, -2)-1) ;
    lua_Integer type = lua_type(L, -1) ;
    switch( type ) {
    case LUA_TBOOLEAN:
      {
        gc.get_bool_at(idx) ;
        vec_bool_type & bv = gc.get_vec_bool() ;
        bv[idx] = lua_toboolean(L, -1) ? true : false ;
      }
      break ;
    case LUA_TNUMBER:
      {
        real_type val = lua_tonumber(L, -1) ;
        if ( gc.get_type() == GC_VEC_REAL || gc.get_type() == GC_MAT_REAL ) {
          gc.get_real_at(idx) = val ;
        } else if ( isInteger(val) ) {
          gc.get_long_at(idx) = long_type(val) ;
        } else {
          gc.get_real_at(idx) = val ;
        }
      }
      break ;
    case LUA_TSTRING:
      gc.get_string_at(idx) = lua_tostring(L, -1) ;
      break ;
    case LUA_TTABLE:
      lua_table_to_GC( L, gc[idx] ) ;
      break ;
    }
  }

  // -----------------------------------------------------------------------------

  static
  void
  push_hash_element( lua_State * L, GenericContainer & gc ) {
    // assegna il valore
    string key  = lua_tostring(L, -2) ;
    int    type = lua_type(L, -1) ;
    switch( type ) {
      case LUA_TBOOLEAN:
        gc[key].set_bool( lua_toboolean(L, -1) ? true : false ) ;
        break ;
      case LUA_TNUMBER:
        {
          real_type val = lua_tonumber(L, -1) ;
          if ( isInteger(val) ) gc[key].set_long(long_type(val)) ;
          else                  gc[key].set_real(val) ;
        }
        break ;
      case LUA_TSTRING:
        gc[key].set_string(lua_tostring(L, -1)) ;
        break ;
      case LUA_TTABLE:
        lua_table_to_GC( L, gc[key] ) ;
        break ;
    }
  }

  // -----------------------------------------------------------------------------

  static
  void
  lua_table_to_GC( lua_State * L, GenericContainer & gc ) {
    // ---
    lua_pushnil(L) ; // first key
    // ---
    while ( lua_next(L,-2) != 0 ) {
      switch( lua_type(L, -2) ) {
        case LUA_TNUMBER:
          push_vec_element( L, gc ) ;
          break ;
        case LUA_TSTRING:
          push_hash_element( L, gc ) ;
          break ;
      }
      lua_pop(L, 1);  // removes `value'; keeps `key' for next iteration
    }
  }

  static
  void
  lua_to_GC( lua_State * L, GenericContainer & gc ) {
    gc.clear() ;
    switch( lua_type(L, -1) ) {
      case LUA_TBOOLEAN:
        gc.set_bool(lua_toboolean(L, -1) ? true : false ) ;
        break ;
      case LUA_TNUMBER:
        {
        real_type val = lua_tonumber(L, -1) ;
        if ( isInteger(val) ) gc.set_long(long_type(val)) ;
        else                  gc.set_real(val) ;
        }
        break ;
      case LUA_TSTRING:
        gc.set_string(lua_tostring(L, -1)) ;
        break ;
      case LUA_TTABLE:
        lua_table_to_GC( L, gc ) ;
        //global_to_GC( gc ) ;
        lua_settop(L, 0);
        break ;
    }
  }

  /*
  //    ____  ____    _           _
  //   / ___|/ ___|  | |_ ___    | |_   _  __ _
  //  | |  _| |      | __/ _ \   | | | | |/ _` |
  //  | |_| | |___   | || (_) |  | | |_| | (_| |
  //   \____|\____|___\__\___/___|_|\__,_|\__,_|
  //             |_____|    |_____|
  */

  static
  void
  GC_to_lua( lua_State * L, GenericContainer const & gc ) {
    // inizializzazione
    switch ( gc.get_type() ) {
    case GC_NOTYPE:
      lua_pushnil(L) ;
      break;
    case GC_POINTER:
      lua_pushnil(L) ;
      break;
    case GC_BOOL:
      lua_pushboolean(L, gc.get_bool() ? 1 : 0 ) ;
      break;
    case GC_INTEGER:
      lua_pushnumber(L,gc.get_int()) ;
      break;
    case GC_LONG:
      lua_pushnumber(L,gc.get_long()) ;
      break;
    case GC_REAL:
      lua_pushnumber(L,gc.get_real()) ;
      break;
    case GC_STRING:
      lua_pushstring(L, gc.get_string().c_str());
      break;
    case GC_VEC_POINTER:
      lua_pushnil(L) ;
      break;
    case GC_VEC_BOOL:
      { vec_bool_type const & vb = gc.get_vec_bool() ;
        lua_createtable(L, int(vb.size()), 0);
        for ( unsigned i = 0 ; i < vb.size() ; ++i ) {
          lua_pushboolean(L, vb[i] ? 1 : 0 ) ;
          lua_rawseti (L, -2, int(i+1));
        }
      }
      lua_pushnil(L) ;
      break;
    case GC_VEC_INTEGER:
      { vec_int_type const & vi = gc.get_vec_int() ;
        lua_createtable(L, int(vi.size()), 0);
        for ( unsigned i = 0 ; i < vi.size() ; ++i ) {
          lua_pushnumber(L, vi[i]);
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type const & vi = gc.get_vec_long() ;
        lua_createtable(L, int(vi.size()), 0);
        for ( unsigned i = 0 ; i < vi.size() ; ++i ) {
          lua_pushnumber(L, vi[i]);
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type const & vr = gc.get_vec_real() ;
        lua_createtable(L, int(vr.size()), 0);
        for ( unsigned i = 0 ; i < vr.size() ; ++i ) {
          lua_pushnumber(L, vr[i]);
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC_VEC_STRING:
      { vec_string_type const & vs = gc.get_vec_string() ;
        lua_createtable(L, int(vs.size()), 0);
        for ( unsigned i = 0 ; i < vs.size() ; ++i ) {
          lua_pushstring(L, vs[i].c_str());
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC_VECTOR:
      { vector_type const & v = gc.get_vector() ;
        lua_createtable(L, int(v.size()), 0);
        for ( unsigned i = 0 ; i < v.size() ; ++i ) {
          GC_to_lua( L, v[i] ) ;
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC_MAP:
      { map_type const & m = gc.get_map() ;
        lua_createtable(L, int(m.size()), 0);
        for ( map_type::const_iterator it = m.begin() ; it != m.end() ; ++it ) {
          lua_pushstring(L, it->first.c_str());
          GC_to_lua( L, it->second ) ;
          lua_settable(L, -3);
        }
      }
      break;
    case GC_COMPLEX:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    //default:
      lua_pushnil(L) ;
      break;
    }
  }

  /*
  //   _                ___       _                           _
  //  | |   _   _  __ _|_ _|_ __ | |_ ___ _ __ _ __  _ __ ___| |_ ___ _ __
  //  | |  | | | |/ _` || || '_ \| __/ _ \ '__| '_ \| '__/ _ \ __/ _ \ '__|
  //  | |__| |_| | (_| || || | | | ||  __/ |  | |_) | | |  __/ ||  __/ |
  //  |_____\__,_|\__,_|___|_| |_|\__\___|_|  | .__/|_|  \___|\__\___|_|
  //                                          |_|
  */

  LuaInterpreter::LuaInterpreter() {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    L = luaL_newstate() ; // opens Lua
    GC_ASSERT_DEBUG( L != NULL, "LuaInterpreter::LuaInterpreter() lua_State invalid!" ) ;
    luaL_openlibs(L) ;
  }

  // -----------------------------------------------------------------------------

  LuaInterpreter::~LuaInterpreter() {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    GC_ASSERT_DEBUG( L != NULL, "LuaInterpreter::~LuaInterpreter() lua_State invalid!" ) ;
    lua_close(L);
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::do_file( char const filename[] ) {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    GC_ASSERT_DEBUG( L != NULL, "LuaInterpreter::do_file('" << filename <<
                                "')\nlua_State invalid!" ) ;
    if ( luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)  ) {
      GC_ASSERT( lua_isnil(L,-1),
                 "In LuaInterpreter::do_file('" << filename << "')\n" << lua_tostring(L, -1) ) ;
    }
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::execute( char const cmd[] ) {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    GC_ASSERT( L != NULL &&
               !luaL_loadbuffer(L, cmd, strlen(cmd), "line") &&
               !lua_pcall(L, 0, 0, 0),
               "In LuaInterpreter::execute('" << cmd <<
               "')\ncannot run the command or lua_State invalid" ) ;
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::call( GenericContainer const & fun_io, GenericContainer & res ) {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;

    // args must be of type MAP
    string_type      const & fname = fun_io("function").get_string() ;
    GenericContainer const & args  = fun_io("args") ;

    // push functions and arguments
    lua_getglobal( L, fname.c_str() );  // function to be called
    GC_to_lua( L, args ) ;

    /* do the call (1 arguments, 1 result) */
    GC_ASSERT( lua_pcall(L, 1, 1, 0) == 0,
               "GenericContainer: error running function `" << fname << "'\n");

    /* retrieve result */
    lua_to_GC( L, res ) ;
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::global_to_GC( char const global_var[], GenericContainer & gc ) {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    GC_ASSERT_DEBUG( L != NULL, "LuaInterpreter::global_to_GC(...) lua_State invalid!" ) ;
      
    lua_getglobal( L, global_var ) ;
    GC_ASSERT( !lua_isnil(L,-1),
               "LuaInterpreter::global_to_GC(...) cannot find global variable: '" << global_var << "'" ) ;

    gc.clear() ;
    switch( lua_type(L, -1) ) {
    case LUA_TBOOLEAN:
      gc.set_bool(lua_toboolean(L, -1) ? true : false ) ;
      break ;
    case LUA_TNUMBER:
      {
        real_type val = lua_tonumber(L, -1) ;
        if ( isInteger(val) ) gc.set_long(long_type(val)) ;
        else                  gc.set_real(val) ;
      }
      break ;
    case LUA_TSTRING:
      gc.set_string(lua_tostring(L, -1)) ;
      break ;
    case LUA_TTABLE:
      lua_table_to_GC( L, gc ) ;
      //global_to_GC( gc ) ;
      lua_settop(L, 0);
      break ;
    default:
      GC_ASSERT( false,
                 "LuaInterpreter::global_to_GC(...) global variable '" <<
                 global_var << "' cannot be converted!" ) ;

    }

  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::GC_to_global( GenericContainer const & gc, char const global_var[] ) {
    lua_State *& L = *(reinterpret_cast<lua_State**>(&void_L)) ;
    GC_to_lua( L, gc ) ;
    lua_setglobal( L, global_var ) ;
  }

  // -----------------------------------------------------------------------------
  extern "C"
  int
  pmain ( lua_State *L ) ;

  extern
  int
  report (lua_State *L, int status);

}

//
// EOF: GenericContainerLuaInterface.cc
//
