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
// file: Utils_trace.hh (header-only version)
//
// \brief Header-only implementation of trace and stack trace utilities
// \author Enrico Bertolazzi
// \date 2017
//

#pragma once

#ifndef UTILS_TRACE_HH
#define UTILS_TRACE_HH

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/**
 * \def UTILS_ERROR_TRACE0(MSG)
 * \brief Macro to throw a Runtime_TraceError with a simple message
 * \param MSG The error message string
 */
#ifndef UTILS_ERROR_TRACE0
#define UTILS_ERROR_TRACE0( MSG ) throw Utils::Runtime_TraceError( MSG, __FILENAME__, __LINE__ )
#endif

/**
 * \def UTILS_ASSERT_TRACE0(COND, MSG)
 * \brief Macro to assert a condition and throw Runtime_TraceError if false
 * \param COND The condition to check
 * \param MSG The error message to display if condition is false
 */
#ifndef UTILS_ASSERT_TRACE0
#define UTILS_ASSERT_TRACE0( COND, MSG ) \
  if ( !( COND ) ) UTILS_ERROR_TRACE0( MSG )
#endif

/**
 * \def UTILS_ERROR_TRACE(...)
 * \brief Macro to throw a Runtime_TraceError with formatted message
 * \param ... Format string and arguments (compatible with fmt::format)
 */
#ifndef UTILS_ERROR_TRACE
#define UTILS_ERROR_TRACE( ... ) throw Utils::Runtime_TraceError( fmt::format( __VA_ARGS__ ), __FILENAME__, __LINE__ )
#endif

/**
 * \def UTILS_ASSERT_TRACE(COND, ...)
 * \brief Macro to assert a condition and throw Runtime_TraceError with formatted message if false
 * \param COND The condition to check
 * \param ... Format string and arguments for the error message
 */
#ifndef UTILS_ASSERT_TRACE
#define UTILS_ASSERT_TRACE( COND, ... ) \
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
  //! especially in debugging scenarios.
  //!
  inline void print_trace( int line, string_view file, string_view msg, ostream_type & stream );

  //!
  //! \deprecated Use print_trace() instead.
  //! \brief Legacy function for printing trace information
  //! \param line The line number where the trace is generated
  //! \param file The source file name
  //! \param reason The reason or message for the trace
  //! \param stream The output stream
  //!
  inline void printTrace( int line, string_view file, string_view reason, ostream_type & stream )
  {
    print_trace( line, file, reason, stream );
  }

  //!
  //! \class Runtime_TraceError
  //! \brief A custom exception class that captures and stores a backtrace on error.
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

    /**
     * \brief Returns the exception message
     * \return C-string containing the exception message with backtrace
     */
    const char * what() const noexcept override { return runtime_error::what(); }
  };

  /*! @} */

}  // namespace Utils

// ============================================================================
// IMPLEMENTATION
// ============================================================================

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils_fmt.hh"

#ifdef UTILS_OS_WINDOWS
#include <windows.h>
#else
#include <cstdio>
#include <cstdlib>
#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <vector>
#ifdef UTILS_OS_OSX
#include <mach-o/dyld.h>
#endif
#endif

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

#ifdef UTILS_OS_OSX
// #define UNW_LOCAL_ONLY
#include <libunwind.h>
#endif

namespace Utils
{

  using std::dec;
  using std::hex;
  using std::ostringstream;
  using std::runtime_error;
  using std::string;

#ifdef UTILS_OS_WINDOWS

  /**
   * \brief Windows implementation of print_trace
   * \details Prints error information and captures stack backtrace on Windows
   * \param line The line number where the error occurred
   * \param file The source file name
   * \param msg The error message
   * \param stream The output stream
   */
  inline void print_trace( int const line, string_view const file, string_view const msg, ostream_type & stream )
  {
    fmt::print(
      stream,
      "---------------------------------------------------------\n"
      "file: {}:{}\n{}\n"
      "---------------------------------------------------------\n",
      file,
      line,
      msg );
#ifndef __MINGW32__
    ULONG const  framesToSkip    = 0;
    ULONG const  framesToCapture = 64;
    void *       backTrace[framesToCapture]{};
    ULONG        backTraceHash = 0;
    USHORT const nFrame        = CaptureStackBackTrace( framesToSkip, framesToCapture, backTrace, &backTraceHash );
    for ( USHORT iFrame = 0; iFrame < nFrame; ++iFrame ) fmt::print( stream, "[{}] = {}\n", iFrame, backTrace[iFrame] );
    fmt::print( stream, "backTraceHash = {0:x}\n", backTraceHash );
#endif
  }

  /**
   * \brief Windows implementation of grab_backtrace
   * \details Creates a formatted error message with file and line information
   * \param reason The reason for the error
   * \param file The source file name
   * \param line The line number
   * \return Formatted error message string
   */
  inline string Runtime_TraceError::grab_backtrace( string_view const reason, string_view const file, int const line )
  {
    return fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line );
  }

#else

#ifdef UTILS_OS_OSX
  /**
   * \brief macOS-specific function to convert address to source line information
   * \details Uses the 'atos' tool to convert memory addresses to source file and line numbers
   * \param binary_path Path to the executable binary
   * \param addr Memory address to resolve
   * \return String containing function name, source file, and line number
   */
  static std::string addr_to_line( std::string const & binary_path, void * addr )
  {
    string cmd{ fmt::format( "atos -o {} {}", binary_path, addr ) };

    FILE * pipe = popen( cmd.data(), "r" );
    if ( !pipe ) return "??";

    char        buffer[256];
    std::string result;
    if ( fgets( buffer, sizeof( buffer ), pipe ) )
    {
      result = buffer;
      result.erase( result.find_last_not_of( " \n\r\t" ) + 1 );  // trim
    }

    pclose( pipe );
    return result;
  }
#else
  /**
   * \brief Demangles C++ symbol names on non-macOS systems
   * \details Uses the C++ ABI demangling function to convert mangled names to human-readable form
   * \param mangled_name The mangled symbol name
   * \return Demangled symbol name
   */
  static string demang( string_view const mangled_name )
  {
    if ( mangled_name.data() == nullptr ) return string{ "" };
    int    status = 0;
    string retval{ mangled_name };
    char * name = abi::__cxa_demangle( mangled_name.data(), nullptr, nullptr, &status );
    if ( status == 0 )
    {
      retval = name;
      // extract only name
      if ( char const * p{ strchr( name, '(' ) }; p != nullptr ) retval = retval.substr( 0, p - name );
    }
    if ( name != nullptr ) free( name );
    return retval;
  }
#endif

  /**
   * \brief Non-Windows implementation of print_trace
   * \details Prints detailed stack trace information including process IDs and function names
   * \param line The line number where the error occurred
   * \param file The source file name
   * \param reason The error message
   * \param stream The output stream
   */
  inline void print_trace( int const line, string_view const file, string_view const reason, ostream_type & stream )
  {
    fmt::print(
      stream,
      "\n{}\nOn File:{}:{}\nprocess ID:{}, parent process "
      "ID:{}\nstack trace:\n",
      reason,
      std::filesystem::path( file ).parent_path().string(),
      line,
      getpid(),
      getppid() );

#ifdef UTILS_OS_OSX
    int const max_frames{ 64 };
    void *    callstack[max_frames];
    int       frames{ backtrace( callstack, max_frames ) };

    Dl_info info;
    char    exe_path[1024];
    ssize_t len = readlink( "/proc/self/exe", exe_path, sizeof( exe_path ) - 1 );
    if ( len == -1 )
    {
      // fallback for macOS
      if ( _NSGetExecutablePath( exe_path, (uint32_t *) &len ) != 0 )
      {
        std::cerr << "Cannot get executable path" << std::endl;
        return;
      }
    }

    for ( int i = 1; i < frames; ++i )
    {
      if ( dladdr( callstack[i], &info ) )
      {
        const char * symname{ info.dli_sname };
        int          status{ 0 };
        char *       demangled{ abi::__cxa_demangle( symname, nullptr, nullptr, &status ) };

        std::string symbol    = ( status == 0 && demangled ) ? demangled : symname;
        std::string file_line = addr_to_line( info.dli_fname, callstack[i] );

        stream << "#" << i << " " << symbol << " at " << file_line << '\n';

        free( demangled );
      }
    }
#else
    // record stack trace upto 128 frames
    void * callstack[128] = {};

    // collect stack frames
    int frames = backtrace( callstack, 128 );

    // get the human-readable symbols (mangled)
    char ** strs = backtrace_symbols( callstack, frames );

    for ( int i = 1; i < frames; ++i )
    {
      Dl_info dlinfo;
      if ( !dladdr( callstack[i], &dlinfo ) ) continue;
      fmt::print( stream, "{:2} {}\n", i, demang( dlinfo.dli_sname ) );
    }
    free( strs );
#endif
  }

  /**
   * \brief Non-Windows implementation of grab_backtrace
   * \details Captures and formats stack trace information for exception messages
   * \param reason The reason for the error
   * \param file The source file name
   * \param line The line number
   * \return Formatted stack trace string
   */
  inline string Runtime_TraceError::grab_backtrace( string_view const reason, string_view const file, int const line )
  {
    ostringstream ost;
    print_trace( line, file, reason, ost );
    return ost.str();
  }
#endif

}  // namespace Utils

#endif  // DOXYGEN_SHOULD_SKIP_THIS

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#endif  // UTILS_TRACE_HH

//
// eof: Utils_trace.hh
//
