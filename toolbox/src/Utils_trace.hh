/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
/// file: Utils_trace.hh
///

#pragma once

#ifndef UTILS_TRACE_HH
#define UTILS_TRACE_HH

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_ERROR_TRACE0
  #define UTILS_ERROR_TRACE0(MSG) \
  throw Utils::Runtime_TraceError( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT_TRACE0
  #define UTILS_ASSERT_TRACE0(COND,MSG) if ( !(COND) ) UTILS_ERROR_TRACE0( MSG )
#endif

#ifndef UTILS_ERROR_TRACE
  #define UTILS_ERROR_TRACE(...) \
  throw Utils::Runtime_TraceError( fmt::format(__VA_ARGS__), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT_TRACE
  #define UTILS_ASSERT_TRACE(COND,...) if ( !(COND) ) UTILS_ERROR_TRACE( __VA_ARGS__ )
#endif

#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace Utils {

  using std::string;
  using std::runtime_error;

  void
  print_trace(
    int                 line,
    char const *        file,
    std::string const & msg,
    ostream_type      & stream
  );

  inline
  void
  printTrace(
    int                 line,
    char const *        file,
    std::string const & msg,
    ostream_type      & stream
  ) {
    print_trace( line, file, msg, stream );
  }

  class Runtime_TraceError : public runtime_error {
  private:
    std::string
    grab_backtrace(
      std::string const & reason,
      char const *        file,
      int                 line
    ) const;

  public:
    explicit
    Runtime_TraceError(
      std::string const & reason,
      char const *        file,
      int                 line
    )
    : std::runtime_error( grab_backtrace( reason, file, line ) )
    { }

    char const * what() const noexcept override;
  };

}

#endif

///
/// eof: Utils_trace.hh
///
