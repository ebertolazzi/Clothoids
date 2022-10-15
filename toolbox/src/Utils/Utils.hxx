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
/// file: Utils.hxx
///

// select computer architecture
#if defined(__APPLE__) && defined(__MACH__)
  // osx architecture
  #define UTILS_OS_OSX 1
  #if defined(__i386__)
    #define UTILS_ARCH32 1
  #elif defined(__x86_64__)
    #define UTILS_ARCH64 1
  #endif
#elif defined(__unix__)
  // linux architecture
  #define UTILS_OS_LINUX 1
  #if defined(__i386__)
    #define UTILS_ARCH32 1
  #elif defined(__x86_64__)
    #define UTILS_ARCH64 1
  #endif
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define UTILS_OS_WINDOWS 1
  // mingw subsystem
  #if defined(__MINGW64__)
    #define UTILS_OS_MINGW 1
    #define UTILS_ARCH64   1
  #elif defined(__MINGW32__)
    #define UTILS_OS_MINGW 1
    #define UTILS_ARCH32   1
  #else
    #if defined(_M_X64) || defined(_M_AMD64) || defined(_WIN64) || defined(WIN64)
      #define UTILS_ARCH64 1
    #else
      #define UTILS_ARCH32 1
    #endif
  #endif
  // windows headers, order matters!
  #include <Winsock2.h>
  #include <Windows.h>
  #include <Ws2tcpip.h>
  #include <iptypes.h>
  #include <Iphlpapi.h>
  // --------------------
  #include <tchar.h>
  #include <stdio.h>
#else
  #error "unsupported OS!"
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
#else
  #error "Lapack Wrapper must be compiled using C++ >= C++11"
#endif

#define UTILS_PURE_VIRTUAL = 0
#ifndef UTILS_OS_WINDOWS
  #define UTILS_OVERRIDE  override
  #define UTILS_CONSTEXPR constexpr
#else
  #define UTILS_OVERRIDE
  #define UTILS_CONSTEXPR
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "fmt/printf.h"
#include "fmt/chrono.h"
#include "fmt/ostream.h"
#endif

// STL
#include <cassert>
#include <iterator>
#include <utility>	    // For std::move(), std::forward()
#include <algorithm>
#include <type_traits>  // For std::remove_reference()
#include <functional>		// For std::bind()

#include <string>
#include <vector>
#include <map>
#include <limits>

// I/O
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

// C/C++
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <memory>

// disable mingw-std-threads for mingw on MATLAB
#if defined(__MINGW32__) || defined(__MINGW64__) && !defined(MATLAB_MEX_FILE )
  #include <_mingw.h>
  #if defined(__MINGW64_VERSION_MAJOR)
    #if  __MINGW64_VERSION_MAJOR < 8
      #define UTILS_USE_MINGW_PORTABLE_THREADS
    #endif
  #endif
  #if defined(__MINGW32_VERSION_MAJOR)
    #if  __MINGW32_VERSION_MAJORs < 8
      #define UTILS_USE_MINGW_PORTABLE_THREADS
    #endif
  #endif
#endif

#ifdef UTILS_USE_MINGW_PORTABLE_THREADS
  #include "mingw-std-threads/mingw.future.h"
  #include "mingw-std-threads/mingw.mutex.h"
  #include "mingw-std-threads/mingw.invoke.h"
  #include "mingw-std-threads/mingw.shared_mutex.h"
  #include "mingw-std-threads/mingw.thread.h"
  #include "mingw-std-threads/mingw.condition_variable.h"
#else
  //#include <future>
  #include <mutex>
  //#include <shared_mutex>
  #include <thread>
  #include <condition_variable>
  #include <atomic>
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

#include "rang.hxx"
#include "Trace.hxx"
#include "Console.hxx"
#include "Malloc.hxx"
#include "Numbers.hxx"
#include "TicToc.hxx"
#include "Quaternion.hxx"
#include "Table.hxx"
#include "Token.hxx"

// order must be preserved
#include "ThreadUtils.hxx"
#include "ThreadPoolBase.hxx"
#include "ThreadPool0.hxx"
#include "ThreadPool1.hxx"
#include "ThreadPool2.hxx"
#include "ThreadPool3.hxx"
#include "ThreadPool4.hxx"
#include "ThreadPool5.hxx"
// -----------------------

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::string;
  #endif

  string basename( char const filename[] );

  void   get_MAC_address( std::map<string,string> & addr );
  string get_host_name();
  void   get_IP_address( std::vector<string> & addr );
  string get_date();
  string get_day_time();
  string get_day_time_and_date();
  string get_log_date_time();
  string get_user_name();
  string get_home_directory();
  string get_executable_path_name();
  bool   check_if_file_exists( char const fname[] );
  bool   check_if_dir_exists( char const dirname[] );
  bool   make_directory( char const dirname[], unsigned mode = 0777 );


  template <typename T_int, typename T_real>
  void
  search_interval(
    T_int        npts,
    T_real const X[],
    T_real     & x,
    T_int      & lastInterval,
    bool         closed,
    bool         can_extend
  );

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  extern template void search_interval(
    int32_t     npts,
    float const X[],
    float     & x,
    int32_t   & lastInterval,
    bool        closed,
    bool        can_extend
  );

  extern template void search_interval(
    int32_t      npts,
    double const X[],
    double     & x,
    int32_t    & lastInterval,
    bool         closed,
    bool         can_extend
  );

  extern template void search_interval(
    int64_t     npts,
    float const X[],
    float     & x,
    int64_t   & lastInterval,
    bool        closed,
    bool        can_extend
  );

  extern template void search_interval(
    int64_t      npts,
    double const X[],
    double     & x,
    int64_t    & lastInterval,
    bool         closed,
    bool         can_extend
  );
  #endif

  template <typename T_int, typename T_real>
  inline
  void
  searchInterval(
    T_int        npts,
    T_real const X[],
    T_real     & x,
    T_int      & lastInterval,
    bool         closed,
    bool         can_extend
  ) {
    search_interval( npts, X, x, lastInterval, closed, can_extend );
  }

  static
  inline
  void
  to_upper( string & str ) {
    for ( auto & c: str ) c = char(toupper(int(c)));
  }

  static
  inline
  void
  to_lower( string & str ) {
    for ( auto & c: str ) c = char(tolower(int(c)));
  }

  static
  inline
  bool
  is_lower( string const & s ) {
    return std::all_of( s.begin(), s.end(), islower );
  }

  static
  inline
  bool
  is_upper( string const & s ) {
    return std::all_of( s.begin(), s.end(), isupper );
  }

  static
  inline
  bool
  is_alpha( string const & s ) {
    return std::all_of( s.begin(), s.end(), isalpha );
  }

  static
  inline
  bool
  is_alphanum( string const & s ) {
    return std::all_of( s.begin(), s.end(), isalnum );
  }

  static
  inline
  bool
  is_digits( string const & s ) {
    return std::all_of( s.begin(), s.end(), isdigit );
  }

  static
  inline
  bool
  is_xdigits( string const & s ) {
    return std::all_of( s.begin(), s.end(), isxdigit );
  }

  // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
  static
  inline
  unsigned
  iLog2( uint32_t N ) {
    static unsigned const tab32[32] = {
       0,  9,  1, 10, 13, 21,  2, 29,
      11, 14, 16, 18, 22, 25,  3, 30,
       8, 12, 20, 28, 15, 17, 24,  7,
      19, 27, 23,  6, 26,  5,  4, 31
    };
    N |= N >> 1;
    N |= N >> 2;
    N |= N >> 4;
    N |= N >> 8;
    N |= N >> 16;
    return tab32[uint32_t(N*0x07C4ACDD)>>27];
  }

  static
  inline
  unsigned
  iLog2( uint64_t N ) {
    static unsigned const tab64[64] = {
      63,  0, 58,  1, 59, 47, 53,  2,
      60, 39, 48, 27, 54, 33, 42,  3,
      61, 51, 37, 40, 49, 18, 28, 20,
      55, 30, 34, 11, 43, 14, 22,  4,
      62, 57, 46, 52, 38, 26, 32, 41,
      50, 36, 17, 19, 29, 10, 13, 21,
      56, 45, 25, 31, 35, 16,  9, 12,
      44, 24, 15,  8, 23,  7,  6,  5
    };
    N |= N >> 1;
    N |= N >> 2;
    N |= N >> 4;
    N |= N >> 8;
    N |= N >> 16;
    N |= N >> 32;
    return tab64[uint64_t((N - (N>>1))*0x07EDD5E59A4E28C2) >> 58];
  }
}

///
/// eof: Utils.hxx
///
