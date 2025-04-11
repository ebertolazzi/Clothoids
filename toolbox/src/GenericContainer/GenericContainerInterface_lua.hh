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

#ifndef GENERIC_CONTAINER_INTERFACE_LUA_HH
#define GENERIC_CONTAINER_INTERFACE_LUA_HH

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include "GenericContainer.hh"

namespace GC_namespace {

  //!
  //! \addtogroup LUA
  //!
  //! @{
  //!

  //!
  //! \brief Convert a Lua global variable to a `GenericContainer`.
  //!
  //! This function retrieves the value of a Lua global variable and converts it into
  //! a `GenericContainer` object. This allows easy access to Lua data from C++.
  //!
  //! @param[in]  void_L     Pointer to the Lua interpreter (lua_State).
  //! @param[in]  global_var Name of the Lua global variable to convert.
  //! @param[out] gc         The resulting `GenericContainer`.
  //!
  void
  Lua_global_to_GC(
    void             * void_L,
    char const       * global_var,
    GenericContainer & gc
  );

  //!
  //! \brief Convert a `GenericContainer` to a Lua global variable.
  //!
  //! This function converts the contents of a `GenericContainer` and stores it in a
  //! Lua global variable. It allows passing data from C++ to Lua.
  //!
  //! \param[in] void_L     Pointer to the Lua interpreter (lua_State).
  //! \param[in] gc         The `GenericContainer` to be converted.
  //! \param[in] global_var Name of the Lua global variable where the result will be stored.
  //!
  void
  Lua_GC_to_global(
    void                   * void_L,
    GenericContainer const & gc,
    char const             * global_var
  );

  //!
  //! \brief A class implementing a simple Lua interpreter.
  //!
  //! The `LuaInterpreter` class provides an interface for loading and executing Lua scripts,
  //! as well as interacting with Lua global variables using the `GenericContainer`.
  //!
  class LuaInterpreter {
    //! Lua interpreter status (lua_State pointer)
    /* lua_State * */ void * void_L;
  public:
    //! Constructor that initializes the Lua interpreter.
    LuaInterpreter();

    //! Destructor that closes the Lua interpreter.
    ~LuaInterpreter();

    //!
    //! \brief Execute a Lua command string.
    //!
    //! This function allows you to run arbitrary Lua code passed as a string.
    //!
    //! \param[in] cmd Lua command to execute.
    //!
    void execute( char const cmd[] );

    //!
    //! \brief Call a Lua function with arguments from a `GenericContainer`.
    //!
    //! This function invokes a Lua function and passes the arguments stored in a `GenericContainer`.
    //! The result of the function is returned in another `GenericContainer`.
    //!
    //! \param[in]  arguments `GenericContainer` containing the Lua function name and arguments.
    //! \param[out] res       `GenericContainer` where the result of the Lua function will be stored.
    //!
    void call( GenericContainer const & arguments, GenericContainer & res );

    //!
    //! \brief Load and execute a Lua script file.
    //!
    //! This function reads and executes a Lua script from a file.
    //!
    //! \param[in] fname Name of the Lua script file.
    //!
    void do_file( char const fname[] );

    //!
    //! \brief Store a `GenericContainer` as a Lua global variable.
    //!
    //! This function converts a `GenericContainer` into a Lua global variable.
    //!
    //! \param[in] gc         The `GenericContainer` to convert.
    //! \param[in] global_var Name of the Lua global variable.
    //!
    void
    GC_to_global( GenericContainer const & gc, char const global_var[] ) {
      Lua_GC_to_global( void_L, gc, global_var );
    }

    //!
    //! \brief Convert a Lua global variable into a `GenericContainer`.
    //!
    //! This function retrieves a Lua global variable and converts it into a `GenericContainer`.
    //!
    //! \param[in]  var Name of the Lua global variable.
    //! \param[out] gc  `GenericContainer` where the Lua variable will be stored.
    //!
    void
    global_to_GC( char const var[], GenericContainer & gc ) {
      Lua_global_to_GC( void_L, var, gc );
    }

    //!
    //! \brief Launch an interactive Lua interpreter mode.
    //!
    //! This function starts an interactive Lua session, allowing users to enter and execute Lua commands.
    //!
    //! \param[in] argc Number of arguments.
    //! \param[in] argv List of command-line arguments.
    //! \param[in] messages List of messages to display during the session.
    //! \param[in] prompt Prompt string to display for the user.
    //! \return Status code of the interpreter session.
    //!
    int
    interactive(
      int argc,
      char const ** argv,
      char const ** messages,
      char const *  prompt
    ); // launch interpret mode
  };

  //!
  //! @}
  //!

}

#endif
