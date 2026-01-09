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
// windows headers, order matters!
#include <Iphlpapi.h>
#include <Windows.h>
#include <Winsock2.h>
#include <Ws2tcpip.h>
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
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// disable mingw-std-threads for mingw on MATLAB
#if defined( __MINGW32__ ) || defined( __MINGW64__ ) && !defined( MATLAB_MEX_FILE )
#include <_mingw.h>
#if defined( __MINGW64_VERSION_MAJOR )
#if __MINGW64_VERSION_MAJOR < 8
#define UTILS_USE_MINGW_PORTABLE_THREADS
#endif
#endif
#if defined( __MINGW32_VERSION_MAJOR )
#if __MINGW32_VERSION_MAJORs < 8
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
#endif

namespace Utils
{
  using string       = std::string;
  using string_view  = std::string_view;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;
}  // namespace Utils

#include "Utils/3rd/spdlog/spdlog.h"
#include "Utils/3rd/spdlog/fmt/bundled/std.h"
#include "Utils/3rd/spdlog/fmt/bundled/chrono.h"
#include "Utils/3rd/spdlog/fmt/bundled/color.h"
#include "Utils/3rd/spdlog/fmt/bundled/ostream.h"
#include "Utils/3rd/spdlog/fmt/bundled/printf.h"

namespace Utils
{
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
      : std::runtime_error( fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line ) )
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
}  // namespace Utils

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
#define UTILS_ERROR( ... ) throw Utils::Runtime_Error( fmt::format( __VA_ARGS__ ), __FILENAME__, __LINE__ )
#endif

#ifndef UTILS_ASSERT
#define UTILS_ASSERT( COND, ... ) \
  if ( !( COND ) ) UTILS_ERROR( __VA_ARGS__ )
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
#define UTILS_ASSERT_DEBUG( COND, ... ) UTILS_ASSERT( COND, __VA_ARGS__ )
#endif
#endif

#include "Malloc.hxx"
#include "Console.hxx"
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
