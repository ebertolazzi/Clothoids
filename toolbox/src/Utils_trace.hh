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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_trace.hh
//

#pragma once

#ifndef UTILS_TRACE_HH
#define UTILS_TRACE_HH

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_ERROR_TRACE0
#define UTILS_ERROR_TRACE0( MSG ) throw Utils::Runtime_TraceError( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT_TRACE0
#define UTILS_ASSERT_TRACE0( COND, MSG )                                                                               \
  if ( !( COND ) ) UTILS_ERROR_TRACE0( MSG )
#endif

#ifndef UTILS_ERROR_TRACE
#define UTILS_ERROR_TRACE( ... ) throw Utils::Runtime_TraceError( fmt::format( __VA_ARGS__ ), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT_TRACE
#define UTILS_ASSERT_TRACE( COND, ... )                                                                                \
  if ( !( COND ) ) UTILS_ERROR_TRACE( __VA_ARGS__ )
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

namespace Utils
{

  using std::runtime_error;
  using std::string;

  /*!
   * \addtogroup OS
   * @{
   */

  //!
  //! \brief Prints a formatted trace message to the specified stream.
  //!
  //! This function outputs a trace message that includes the line number, file
  //! name, and a custom message to help in debugging or logging operations.
  //!
  //! \param line   The line number in the source file where the trace is
  //! generated.
  //! \param file   The name of the source file where the trace is generated.
  //! \param msg    A custom message that provides additional information for
  //! the trace.
  //! \param stream The output stream where the trace will be printed (e.g.,
  //! std::cout, std::cerr).
  //!
  //! \note This function is useful for logging and tracking the execution flow,
  //! especially in
  //!       debugging scenarios.
  //!/
  void print_trace( int line, string_view file, string_view msg, ostream_type & stream );

  //!
  //! \deprecated Use print_trace() instead.
  //!
  inline void
  printTrace( int line, string_view file, string_view reason, ostream_type & stream )
  {
    print_trace( line, file, reason, stream );
  }

  //!
  //! \class Runtime_TraceError
  //! \brief A custom exception class that captures and stores a backtrace on
  //! error.
  //!
  //! `Runtime_TraceError` is a subclass of `std::runtime_error` that captures
  //! the backtrace information (file, line, and reason) when an error occurs,
  //! making it easier to debug.
  //!
  //! The exception message includes the backtrace details, which can be
  //! accessed via the `what()` method.
  //!
  class Runtime_TraceError : public runtime_error
  {
  private:
    //!
    //! \brief Captures and formats a backtrace.
    //!
    //! This function generates a string containing backtrace information, which
    //! includes the reason for the error, the file where it occurred, and the
    //! line number. This formatted string is used as the exception message.
    //!
    //! \param reason The reason or description of the error.
    //! \param file   The name of the source file where the error occurred.
    //! \param line   The line number where the error occurred.
    //!
    //! \return A formatted string containing the backtrace information.
    //!
    static string grab_backtrace( string_view reason, string_view file, int line );

  public:
    //!
    //! \brief Constructs a `Runtime_TraceError` with a backtrace.
    //!
    //! This constructor initializes the exception with a backtrace that
    //! includes the error reason, the file name, and the line number where the
    //! error occurred.
    //!
    //! \param reason A description of the error or exception reason.
    //! \param file   The name of the source file where the error occurred.
    //! \param line   The line number where the error occurred.
    //!
    explicit Runtime_TraceError( string_view reason, string_view file, int line )
      : std::runtime_error( grab_backtrace( reason, file, line ) )
    {
    }

    char const * what() const noexcept override;
  };

  /*! @} */

}  // namespace Utils

#endif

//
// eof: Utils_trace.hh
//
