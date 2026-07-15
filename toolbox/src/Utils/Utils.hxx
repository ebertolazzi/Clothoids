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
// file: Utils.hxx
//

// ============================================================================
// Compiler and System Detection
// ============================================================================

// select computer architecture
#if defined( __APPLE__ ) && defined( __MACH__ )
// osx architecture
#define UTILS_OS_OSX 1
#if defined( __i386__ )
#define UTILS_ARCH32 1
#elif defined( __x86_64__ )
#define UTILS_ARCH64 1
#endif
#elif defined( __unix__ )
// linux architecture
#define UTILS_OS_LINUX 1
#if defined( __i386__ )
#define UTILS_ARCH32 1
#elif defined( __x86_64__ )
#define UTILS_ARCH64 1
#endif
#elif defined( _WIN32 ) || defined( WIN32 ) || defined( _WIN64 ) || defined( WIN64 )
// windows architecture
#define UTILS_OS_WINDOWS 1
// mingw subsystem
#if defined( __MINGW64__ )
#define UTILS_OS_MINGW 1
#define UTILS_ARCH64 1
#elif defined( __MINGW32__ )
#define UTILS_OS_MINGW 1
#define UTILS_ARCH32 1
#else
#if defined( _M_X64 ) || defined( _M_AMD64 ) || defined( _WIN64 ) || defined( WIN64 )
#define UTILS_ARCH64 1
#else
#define UTILS_ARCH32 1
#endif
#endif
// Windows headers: include order matters.
// winsock2.h must come before windows.h.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Winsock2.h>
#include <Ws2tcpip.h>
#include <Windows.h>
#include <Iphlpapi.h>
#include <iptypes.h>
// --------------------
#include <stdio.h>
#include <tchar.h>
#else
#error "unsupported OS!"
#endif

// check if compiler is C++17
#if ( defined( _MSC_VER ) && _MSC_VER >= 1910 ) || ( defined( __cplusplus ) && __cplusplus >= 201703L )

#if ( defined( __cplusplus ) && __cplusplus <= 201103L )
#define UTILS_DEFAULT \
  {                   \
  }
#else
#define UTILS_DEFAULT = default
#endif
#else
#error "Utils library must be compiled using C++ >= C++17"
#endif

// ============================================================================
// Standard Library Headers
// ============================================================================

#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>


// disable mingw-std-threads for mingw on MATLAB
#if ( defined( __MINGW32__ ) || defined( __MINGW64__ ) ) && !defined( MATLAB_MEX_FILE )
#include <_mingw.h>
#if defined( __MINGW64_VERSION_MAJOR )
#if __MINGW64_VERSION_MAJOR < 8
#define UTILS_USE_MINGW_PORTABLE_THREADS
#endif
#endif
#if defined( __MINGW32_VERSION_MAJOR )
#if __MINGW32_VERSION_MAJOR < 8
#define UTILS_USE_MINGW_PORTABLE_THREADS
#endif
#endif
#endif

#ifdef UTILS_USE_MINGW_PORTABLE_THREADS
#include "mingw-std-threads/mingw.condition_variable.h"
#include "mingw-std-threads/mingw.future.h"
#include "mingw-std-threads/mingw.invoke.h"
#include "mingw-std-threads/mingw.mutex.h"
#include "mingw-std-threads/mingw.shared_mutex.h"
#include "mingw-std-threads/mingw.thread.h"
#else
#include <atomic>
#include <condition_variable>
#include <future>
#include <mutex>
#include <shared_mutex>
#include <thread>
#endif

#ifdef _MSC_VER
// Workaround for visual studio
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef ERROR
#undef ERROR
#endif
#endif

namespace Utils
{
  using string       = std::string;
  using string_view  = std::string_view;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;
  using std::runtime_error;

  //!
  //! \brief Custom runtime error class for handling runtime exceptions.
  //!
  //! This class extends the standard `std::runtime_error` to include additional
  //! context information, specifically the file name and line number where the
  //! error occurred. It provides constructors that accept a reason for the
  //! error and formats the error message accordingly.
  //!
  //! **Usage**
  //!
  //! \code
  //! try {
  //!     throw Runtime_Error("An error occurred", __FILE__, __LINE__);
  //! } catch (const Runtime_Error& e) {
  //!     std::cerr << e.what();
  //! }
  //! \endcode

  class Runtime_Error : public runtime_error
  {
  public:
    //!
    //! \brief Constructs a Runtime_Error instance with a given reason.
    //!
    //! This constructor initializes the error with a specified reason,
    //! the file where the error occurred, and the line number. It formats
    //! the error message accordingly.
    //!
    //! \param reason A string that describes the reason for the error.
    //! \param file The name of the file where the error occurred.
    //! \param line The line number in the file where the error occurred.
    //!
    explicit Runtime_Error( string_view reason, string_view file, int line )
      : std::runtime_error( std::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
    {
    }

    //!
    //! \brief Returns a C-style string describing the error.
    //!
    //! This method overrides the `what()` method from `std::runtime_error`
    //! to provide a more detailed error message, including the reason for
    //! the error, the file name, and the line number.
    //!
    //! \return A C-style string representing the error message.
    //!
    char const * what() const noexcept override { return runtime_error::what(); }
  };

  /**
   * @brief Error handling utilities with source location and formatted messages.
   *
   * This module provides a set of utilities for error handling, assertions,
   * warnings, and debug checks with automatic source location tracking and
   * C++20 format string support.
   *
   * @note Requires C++20 for std::source_location and std::format.
   * @see Utils::Runtime_Error
   */

  /**
   * @brief Throws a runtime error with the given message and source location.
   *
   * This function constructs a Utils::Runtime_Error exception using the provided
   * error message and captures the source location where the error was triggered.
   *
   * @param msg The error message to be included in the exception.
   * @param loc The source location (file name, line number) where the error occurred.
   *            Defaults to the caller's location via std::source_location::current().
   *
   * @throws Utils::Runtime_Error Always throws with the provided message and location.
   *
   * @example
   * @code
   * if (data.empty()) {
   *   Error("Data container is empty");
   * }
   * @endcode
   */
  inline void Error( std::string_view msg, std::source_location loc = std::source_location::current() )
  { throw Utils::Runtime_Error( std::string{ msg }, loc.file_name(), loc.line() ); }

  /**
   * @brief Asserts a condition and throws an exception if it evaluates to false.
   *
   * This function checks the provided condition and throws a Utils::Runtime_Error
   * with the specified error message and source location if the condition is false.
   *
   * @param ok The condition to check. If false, an exception is thrown.
   * @param msg The error message to include in the exception.
   * @param loc The source location where the assertion was triggered.
   *            Defaults to the caller's location.
   *
   * @throws Utils::Runtime_Error If the condition is false.
   *
   * @example
   * @code
   * Assert(ptr != nullptr, "Pointer must not be null");
   * Assert(value > 0, "Value must be positive");
   * @endcode
   */
  inline void Assert( bool ok, std::string_view msg, std::source_location loc = std::source_location::current() )
  {
    if ( !ok ) throw Utils::Runtime_Error( std::string{ msg }, loc.file_name(), loc.line() );
  }

  /**
   * @brief Displays a warning message if the condition is false.
   *
   * This function checks the provided condition and outputs a warning to std::cout
   * with the error message and source location if the condition is false.
   * Unlike Assert, this function does not throw an exception.
   *
   * @param ok The condition to check. If false, a warning is displayed.
   * @param msg The warning message to display.
   * @param loc The source location where the warning was triggered.
   *            Defaults to the caller's location.
   *
   * @example
   * @code
   * Warning(config.is_valid(), "Configuration has invalid parameters");
   * Warning(file_exists, "File not found, continuing anyway");
   * @endcode
   */
  inline void Warning( bool ok, std::string_view msg, std::source_location loc = std::source_location::current() )
  {
    if ( !ok ) std::cout << std::format( "{}\nfile: {}, line: {}\n", msg, loc.file_name(), loc.line() );
  }

  /**
   * @brief Helper struct for format strings with source location.
   *
   * This struct combines a C++20 format string with its source location,
   * enabling formatted error messages that automatically capture the call site.
   *
   * @tparam Args The types of the arguments to be formatted.
   *
   * @note This struct is designed to be used with the format-enabled versions
   *       of Error, Check, Assert, Warning, and Debug.
   *
   * @see Format_With_Location::Format_With_Location
   */
  template <class... Args> struct Format_With_Location
  {
    std::format_string<Args...> fmt;  ///< The format string with placeholders for arguments.
    std::source_location        loc;  ///< The source location where the format string was defined.

    /**
     * @brief Constructs a Format_With_Location object with compile-time format validation.
     *
     * This constructor is consteval (compile-time evaluated) to ensure that
     * the format string is valid at compile time.
     *
     * @tparam S The type of the format string (deduced).
     * @param s The format string containing placeholders for arguments.
     * @param l The source location of the caller. Defaults to current location.
     *
     * @note The format string is validated at compile time, preventing runtime
     *       format errors.
     */
    template <class S>
    consteval Format_With_Location( S && s, std::source_location l = std::source_location::current() )
      : fmt( std::forward<S>( s ) ), loc( l )
    {
    }
  };

  /**
   * @brief Throws a runtime error with a formatted message and source location.
   *
   * This function uses C++20 formatting to construct the error message and
   * captures the source location where the error was triggered.
   *
   * @tparam Args The types of the arguments to format (deduced).
   * @param f The format string with source location.
   * @param args The arguments to format into the message.
   *
   * @throws Utils::Runtime_Error Always throws with the formatted message and location.
   *
   * @example
   * @code
   * int id = 42;
   * std::string name = "item";
   * Error("Cannot process {} with ID {}", name, id);
   * @endcode
   */
  template <class... Args> inline void Error( Format_With_Location<std::type_identity_t<Args>...> f, Args &&... args )
  {
    throw Utils::Runtime_Error( std::format( f.fmt, std::forward<Args>( args )... ), f.loc.file_name(), f.loc.line() );
  }

  /**
   * @brief Checks a condition and throws if false, with a formatted error message.
   *
   * This function evaluates a condition and throws a formatted error message
   * if the condition evaluates to false. The error message is constructed using
   * C++20 format string syntax.
   *
   * @tparam Args The types of the arguments to format (deduced).
   * @param ok The condition to check. If false, an exception is thrown.
   * @param f The format string with source location.
   * @param args The arguments to format into the error message.
   *
   * @throws Utils::Runtime_Error If the condition is false.
   *
   * @example
   * @code
   * Check(value > 0, "Value {} must be positive", value);
   * Check(index < size, "Index {} out of bounds (size={})", index, size);
   * @endcode
   */
  template <class... Args>
  inline void Check( bool ok, Format_With_Location<std::type_identity_t<Args>...> f, Args &&... args )
  {
    if ( !ok )
      throw Utils::Runtime_Error(
        std::format( f.fmt, std::forward<Args>( args )... ),
        f.loc.file_name(),
        f.loc.line() );
  }

  /**
   * @brief Alias of Check for compatibility with standard assertion naming.
   *
   * This function is functionally identical to Check() and is provided as
   * a convenience alias to maintain compatibility with code that expects
   * an Assert() function with formatted messages.
   *
   * @tparam Args The types of the arguments to format (deduced).
   * @param ok The condition to check. If false, an exception is thrown.
   * @param f The format string with source location.
   * @param args The arguments to format into the error message.
   *
   * @throws Utils::Runtime_Error If the condition is false.
   *
   * @see Check
   *
   * @example
   * @code
   * Assert(ptr != nullptr, "Pointer {} is null", ptr_name);
   * Assert(value > 0, "Invalid value: {}", value);
   * @endcode
   */
  template <class... Args>
  inline void Assert( bool ok, Format_With_Location<std::type_identity_t<Args>...> f, Args &&... args )
  {
    if ( !ok )
      throw Utils::Runtime_Error(
        std::format( f.fmt, std::forward<Args>( args )... ),
        f.loc.file_name(),
        f.loc.line() );
  }

  /**
   * @brief Displays a formatted warning message if the condition is false.
   *
   * This function checks the provided condition and outputs a formatted warning
   * to std::cout if the condition is false. Unlike Assert/Check, this function
   * does not throw an exception.
   *
   * @tparam Args The types of the arguments to format (deduced).
   * @param ok The condition to check. If false, a warning is displayed.
   * @param f The format string with source location.
   * @param args The arguments to format into the warning message.
   *
   * @example
   * @code
   * Warning(retry_count == 0, "Retry attempt {} of {}", retry_count, max_retries);
   * Warning(config.is_valid(), "Invalid configuration: {}", config.error_message());
   * @endcode
   */
  template <class... Args>
  inline void Warning( bool ok, Format_With_Location<std::type_identity_t<Args>...> f, Args &&... args )
  {
    if ( !ok )
      std::cout << std::format(
        "{}\nfile: {}, line: {}\n",
        std::format( f.fmt, std::forward<Args>( args )... ),
        f.loc.file_name(),
        f.loc.line() );
  }

/**
 * @brief Debug assertion that can be conditionally compiled out.
 *
 * This function behaves like Assert() when NDEBUG is not defined,
 * and becomes a no-op when NDEBUG is defined. This allows debug
 * assertions to be removed from release builds for performance reasons.
 *
 * @note When NDEBUG is defined, this function does nothing and the
 *       condition and arguments are not evaluated (due to template
 *       instantiation behavior).
 *
 * @tparam Args The types of the arguments to format (deduced).
 * @param ok The condition to check in debug builds.
 * @param f The format string with source location.
 * @param args The arguments to format into the debug assertion message.
 *
 * @see Assert
 * @see NDEBUG
 *
 * @example
 * @code
 * // This will only be checked in debug builds
 * Debug(vector.size() > 0, "Vector is empty after operation");
 *
 * // More complex debug check with formatting
 * Debug(value > 0, "Invalid value {} at iteration {}", value, iter);
 * @endcode
 */
#ifdef NDEBUG
  // When NDEBUG is defined, Debug() becomes a no-op
  template <class... Args>
  inline void Debug( bool /*ok*/, Format_With_Location<std::type_identity_t<Args>...> /*f*/, Args &&... /*args*/ )
  {
    // Empty implementation - all debug checks are removed in release builds
  }
  inline void Debug( bool ok, std::string_view msg, std::source_location loc = std::source_location::current() )
  {
    if ( !ok ) throw Utils::Runtime_Error( std::string{ msg }, loc.file_name(), loc.line() );
  }

#else
  // When NDEBUG is NOT defined, Debug() behaves like Assert()
  template <class... Args>
  inline void Debug( bool ok, Format_With_Location<std::type_identity_t<Args>...> f, Args &&... args )
  {
    if ( !ok )
      throw Utils::Runtime_Error(
        std::format( f.fmt, std::forward<Args>( args )... ),
        f.loc.file_name(),
        f.loc.line() );
  }
  inline void Debug( bool ok, std::string_view msg, std::source_location loc = std::source_location::current() )
  {
    if ( !ok ) throw Utils::Runtime_Error( std::string{ msg }, loc.file_name(), loc.line() );
  }
#endif
}  // namespace Utils


// OBSOLETE MASCROS WILL BE REMOVED IN FUTURE RELEASE --- BEGIN

#ifndef __FILENAME__
#define __FILENAME__ ( strrchr( __FILE__, '/' ) ? strrchr( "/" __FILE__, '/' ) + 1 : __FILE__ )
#endif

#ifndef UTILS_ERROR0
#define UTILS_ERROR0( MSG ) throw Utils::Runtime_Error( MSG, __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT0
#define UTILS_ASSERT0( COND, MSG ) \
  if ( !( COND ) ) UTILS_ERROR0( MSG )
#endif

#ifndef UTILS_WARNING0
#define UTILS_WARNING0( COND, MSG ) \
  if ( !( COND ) ) std::cerr << MSG
#endif

#ifndef UTILS_ERROR
#define UTILS_ERROR( ... ) throw Utils::Runtime_Error( std::format( __VA_ARGS__ ), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT
#define UTILS_ASSERT( COND, ... ) \
  if ( !( COND ) ) Utils::Error( __VA_ARGS__ )
#endif

#ifndef UTILS_WARNING
#define UTILS_WARNING( COND, ... ) \
  if ( !( COND ) ) fmt::print( __VA_ARGS__ )
#endif

#ifdef UTILS_NO_DEBUG
#ifndef UTILS_ASSERT0_DEBUG
#define UTILS_ASSERT0_DEBUG( COND, MSG )
#endif
#ifndef UTILS_ASSERT_DEBUG
#define UTILS_ASSERT_DEBUG( COND, ... )
#endif
#else
#ifndef UTILS_ASSERT0_DEBUG
#define UTILS_ASSERT0_DEBUG( COND, MSG ) UTILS_ASSERT0( COND, MSG )
#endif
#ifndef UTILS_ASSERT_DEBUG
#define UTILS_ASSERT_DEBUG( COND, ... ) Utils::Check( COND, __VA_ARGS__ )
#endif
#endif

// OBSOLETE MASCROS WILL BE REMOVED IN FUTURE RELEASE --- END

#include "Malloc.hxx"
#include "Numbers.hxx"

// order must be preserved
#include "ThreadPoolBase.hxx"
#include "ThreadUtils.hxx"
#include "ThreadPool0.hxx"
#include "ThreadPool1.hxx"
#include "ThreadPool2.hxx"
#include "ThreadPool3.hxx"
#include "ThreadPool4.hxx"
#include "ThreadPool5.hxx"
#include "ThreadPoolEigen.hxx"
// -----------------------

// ============================================================================
// Utility Functions
// ============================================================================

namespace Utils
{

  // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
  static inline unsigned iLog2( uint32_t N )
  {
    static unsigned const tab32[32] = { 0, 9,  1,  10, 13, 21, 2,  29, 11, 14, 16, 18, 22, 25, 3, 30,
                                        8, 12, 20, 28, 15, 17, 24, 7,  19, 27, 23, 6,  26, 5,  4, 31 };
    N |= N >> 1;
    N |= N >> 2;
    N |= N >> 4;
    N |= N >> 8;
    N |= N >> 16;
    return tab32[uint32_t( N * 0x07C4ACDD ) >> 27];
  }

  static inline unsigned iLog2( uint64_t N )
  {
    static unsigned const tab64[64] = { 63, 0,  58, 1,  59, 47, 53, 2,  60, 39, 48, 27, 54, 33, 42, 3,
                                        61, 51, 37, 40, 49, 18, 28, 20, 55, 30, 34, 11, 43, 14, 22, 4,
                                        62, 57, 46, 52, 38, 26, 32, 41, 50, 36, 17, 19, 29, 10, 13, 21,
                                        56, 45, 25, 31, 35, 16, 9,  12, 44, 24, 15, 8,  23, 7,  6,  5 };
    N |= N >> 1;
    N |= N >> 2;
    N |= N >> 4;
    N |= N >> 8;
    N |= N >> 16;
    N |= N >> 32;
    return tab64[uint64_t( ( N - ( N >> 1 ) ) * 0x07EDD5E59A4E28C2 ) >> 58];
  }

}  // namespace Utils

//
// eof: Utils.hxx
//
