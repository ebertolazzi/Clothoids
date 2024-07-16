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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "GenericContainer.hh"

namespace GC_namespace {

  //!
  //! Convert a global variable in `Lua` to a `GenericContainer`
  //!
  //! \param void_L pointer to a Lua interpreter
  //! \param[in]    global_var string containing the name of the Lua global variable to be converted
  //! \param[out]   gc resulting `GenericContainer`
  //!
  void
  Lua_global_to_GC(
    void             * void_L,
    char const       * global_var,
    GenericContainer & gc
  );

  //!
  //! Convert a `GenericContainer` to a global variable in the `Lua` interpreter
  //!
  //! \param void_L pointer to a Lua interpreter
  //! \param[in]    gc input `GenericContainer`
  //! \param[in]    global_var string containing the name of the Lua global variable storing the result of conversion
  //!
  void
  Lua_GC_to_global(
    void                   * void_L,
    GenericContainer const & gc,
    char const             * global_var
  );

  //!
  //! C++ class imnplementing a simple Lua interpreter that can be used
  //! to read and interpret data file.
  //!
  class LuaInterpreter {
    //!
    //! Interpreter status
    //!
    /* lua_State * */ void * void_L;
  public:
    LuaInterpreter();
    ~LuaInterpreter();

    //!
    //! Interpret the string `cmd` as a Lua statement
    //!
    void execute( char const cmd[] );

    //!
    //! Execute a function in Lua with arguments passed by the `GenericContainer`
    //!
    //! - `arguments` must contain the field
    //! - "function" of type string with the name of the Lua function to be called
    //! - "args"     a generic container storing the arguments of the function
    //!
    //! The result of computation is returned in `res`.
    //!
    void call( GenericContainer const & arguments, GenericContainer & res );

    //!
    //! Interpret the file `fname` as a Lua script.
    //!
    void do_file( char const fname[] );

    //!
    //! Get `gc` and store in the Lua global variable `global_var`.
    //!
    void
    GC_to_global( GenericContainer const & gc, char const global_var[] ) {
      Lua_GC_to_global( void_L, gc, global_var );
    }

    //!
    //! Get a Lua global variable `global_var` and store in `gc`.
    //!
    void
    global_to_GC( char const var[], GenericContainer & gc ) {
      Lua_global_to_GC( void_L, var, gc );
    }

    //!
    //! Execute the Lua interpreter interatively.
    //!
    int
    interactive(
      int argc,
      char const ** argv,
      char const ** messages,
      char const *  prompt
    ); // launch interpret mode
  };
}

#endif
