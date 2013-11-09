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

// load string.h for strlen
#include <string.h>

using namespace std ;

namespace GC {

  typedef double valueType ;
  typedef int    indexType ;

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
  lua_to_GC( lua_State        * L,
             GenericContainer & gc,
             string const     & indent ) ;

  static
  void
  push_vec_element( lua_State        * L,
                    GenericContainer & gc,
                    string const     & indent ) {
    // assegna il valore
    // index start from 1 in LUA
    unsigned    idx  = unsigned(lua_tointeger(L, -2)-1) ;
    lua_Integer type = lua_type(L, -1) ;
    switch( type ) {
    case LUA_TBOOLEAN:
      {
        gc.get_bool(idx) ;
        vec_bool_type & bv = gc.get_vec_bool() ;
        bv[idx] = lua_toboolean(L, -1) ;
      }
      break ;
    case LUA_TNUMBER:
      {
        valueType val = lua_tonumber(L, -1) ;
        if ( gc.get_type() == GC::GC_VEC_REAL ) {
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

  // -----------------------------------------------------------------------------

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

  // -----------------------------------------------------------------------------

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
    case GC::GC_NOTYPE:
      lua_pushnil(L) ;
      break;
    case GC::GC_POINTER:
      lua_pushnil(L) ;
      break;
    case GC::GC_BOOL:
      lua_pushboolean(L, gc.get_bool() ? 1 : 0 ) ;
      break;
    case GC::GC_INT:
      lua_pushnumber(L,gc.get_int()) ;
      break;
    case GC::GC_REAL:
      lua_pushnumber(L,gc.get_real()) ;
      break;
    case GC::GC_STRING:
      lua_pushstring(L, gc.get_string().c_str());
      break;
    case GC::GC_VEC_POINTER:
      lua_pushnil(L) ;
      break;
    case GC::GC_VEC_BOOL:
      { GC::vec_bool_type const & vb = gc.get_vec_bool() ;
        lua_createtable(L, int(vb.size()), 0);
        for ( unsigned i = 0 ; i < vb.size() ; ++i ) {
          lua_pushboolean(L, vb[i] ? 1 : 0 ) ;
          lua_rawseti (L, -2, int(i+1));
        }
      }
      lua_pushnil(L) ;
      break;
    case GC::GC_VEC_INT:
      { GC::vec_int_type const & vi = gc.get_vec_int() ;
        lua_createtable(L, int(vi.size()), 0);
        for ( unsigned i = 0 ; i < vi.size() ; ++i ) {
          lua_pushnumber(L, vi[i]);
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC::GC_VEC_REAL:
      { GC::vec_real_type const & vr = gc.get_vec_real() ;
        lua_createtable(L, int(vr.size()), 0);
        for ( unsigned i = 0 ; i < vr.size() ; ++i ) {
          lua_pushnumber(L, vr[i]);
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC::GC_VEC_STRING:
      { GC::vec_string_type const & vs = gc.get_vec_string() ;
        lua_createtable(L, int(vs.size()), 0);
        for ( unsigned i = 0 ; i < vs.size() ; ++i ) {
          lua_pushstring(L, vs[i].c_str());
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC::GC_VECTOR:
      { GC::vector_type const & v = gc.get_vector() ;
        lua_createtable(L, int(v.size()), 0);
        for ( unsigned i = 0 ; i < v.size() ; ++i ) {
          GC_to_lua( L, v[i] ) ;
          lua_rawseti (L, -2, int(i+1));
        }
      }
      break;
    case GC::GC_MAP:
      { GC::map_type const & m = gc.get_map() ;
        lua_createtable(L, int(m.size()), 0);
        for ( GC::map_type::const_iterator it = m.begin() ; it != m.end() ; ++it ) {
          lua_pushstring(L, it->first.c_str());
          GC_to_lua( L, it->second ) ;
          lua_settable(L, -3);
        }
      }
      break;
    default:
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
    lua_State *& L = *((lua_State**)&void_L) ;
    L = luaL_newstate() ; // opens Lua
    ASSERT_DEBUG( L != NULL, "LuaInterpreter::LuaInterpreter() lua_State invalid!" ) ;
    luaL_openlibs(L) ;
  }

  // -----------------------------------------------------------------------------

  LuaInterpreter::~LuaInterpreter() {
    lua_State *& L = *((lua_State**)&void_L) ;
    ASSERT_DEBUG( L != NULL, "LuaInterpreter::~LuaInterpreter() lua_State invalid!" ) ;
    lua_close(L);
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::do_file( char const filename[], bool check_syntax_only ) {
    lua_State *& L = *((lua_State**)&void_L) ;
    ASSERT_DEBUG( L != NULL, "LuaInterpreter::do_file('" << filename <<
                             "')\nlua_State invalid!" ) ;
    if ( luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)  ) {
      ASSERT( lua_isnil(L,-1),
              "In LuaInterpreter::do_file('" << filename << "')\n" << lua_tostring(L, -1) ) ;
    }
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::execute( char const cmd[] ) {
    lua_State *& L = *((lua_State**)&void_L) ;
    ASSERT( L != NULL &&
            !luaL_loadbuffer(L, cmd, strlen(cmd), "line") &&
            !lua_pcall(L, 0, 0, 0),
            "In LuaInterpreter::execute('" << cmd <<
            "')\ncannot run the command or lua_State invalid" ) ;
  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::global_to_GC( char const global_var[], GenericContainer & gc ) {
    lua_State *& L = *((lua_State**)&void_L) ;
    ASSERT_DEBUG( L != NULL, "LuaInterpreter::global_to_GC(...) lua_State invalid!" ) ;
      
    lua_getglobal( L, global_var ) ;
    ASSERT( !lua_isnil(L,-1),
            "LuaInterpreter::global_to_GC(...) cannot find global variable: '" << global_var << "'" ) ;

    gc.clear() ;
    switch( lua_type(L, -1) ) {
    case LUA_TBOOLEAN:
      gc.set_bool(lua_toboolean(L, -1)) ;
      break ;
    case LUA_TNUMBER:
      {
        valueType  val = lua_tonumber(L, -1) ;
        indexType ival = indexType(val) ;
        if ( ival == val ) gc.set_int(ival) ;
        else               gc.set_real(val) ;
      }
      break ;
    case LUA_TSTRING:
      gc.set_string(lua_tostring(L, -1)) ;
      break ;
    case LUA_TTABLE:
      lua_to_GC( L, gc, "" ) ;
      //global_to_GC( gc ) ;
      lua_settop(L, 0);
      break ;
    default:
      ASSERT( false,
              "LuaInterpreter::global_to_GC(...) global variable '" <<
              global_var << "' cannot be converted!" ) ;

    }

  }

  // -----------------------------------------------------------------------------

  void
  LuaInterpreter::GC_to_global( GenericContainer const & gc, char const global_var[] ) {
    lua_State *& L = *((lua_State**)&void_L) ;
    GC_to_lua( L, gc ) ;
    lua_setglobal( L, global_var ) ;
  }
}
//
// EOF: GenericContainerLuaInterface.cc
//
