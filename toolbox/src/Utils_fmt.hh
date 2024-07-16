/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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

///
/// file: Utils_fmt.hh
///

#pragma once

#ifndef UTILS_FMT_dot_HH
#define UTILS_FMT_dot_HH

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Utils/fmt/printf.h"
#include "Utils/fmt/chrono.h"
#include "Utils/fmt/ostream.h"
#include "Utils/fmt/color.h"
#include "Utils/fmt/std.h"
#endif

///
/// file: Trace.hxx
///

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <string.h>
#ifndef __FILENAME__
  #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr("/" __FILE__, '/') + 1 : __FILE__)
#endif

#ifndef UTILS_ERROR0
  #define UTILS_ERROR0(MSG) \
  throw Utils::Runtime_Error( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT0
  #define UTILS_ASSERT0(COND,MSG) if ( !(COND) ) UTILS_ERROR0( MSG )
#endif

#ifndef UTILS_WARNING0
  #define UTILS_WARNING0(COND,MSG) if ( !(COND) ) std::cerr << MSG
#endif

#ifndef UTILS_ERROR
  #define UTILS_ERROR(...) \
  throw Utils::Runtime_Error( fmt::format(__VA_ARGS__), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT
  #define UTILS_ASSERT(COND,...) if ( !(COND) ) UTILS_ERROR( __VA_ARGS__ )
#endif

#ifndef UTILS_WARNING
  #define UTILS_WARNING(COND,...) if ( !(COND) ) fmt::print( __VA_ARGS__ )
#endif

#ifdef UTILS_NO_DEBUG
  #ifndef UTILS_ASSERT0_DEBUG
    #define UTILS_ASSERT0_DEBUG(COND,MSG)
  #endif
  #ifndef UTILS_ASSERT_DEBUG
    #define UTILS_ASSERT_DEBUG(COND,...)
  #endif
#else
  #ifndef UTILS_ASSERT0_DEBUG
    #define UTILS_ASSERT0_DEBUG(COND,MSG) UTILS_ASSERT0(COND,MSG)
  #endif
  #ifndef UTILS_ASSERT_DEBUG
    #define UTILS_ASSERT_DEBUG(COND,...) UTILS_ASSERT(COND,__VA_ARGS__)
  #endif
#endif

#endif

namespace Utils {

  using std::runtime_error;

  class Runtime_Error : public runtime_error {
  public:
    explicit
    Runtime_Error(
      std::string const & reason,
      char const *        file,
      int                 line
    )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    explicit
    Runtime_Error(
      char const * reason,
      char const * file,
      int          line
    )
    : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    { }

    char const * what() const noexcept override;
  };

}

#endif

///
/// EOF: Utils_fmt.hh
///
