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
  #if defined(_M_X64) || defined(_M_AMD64) || defined(_WIN64) || defined(WIN64)
    #define UTILS_ARCH64 1
  #else
    #define UTILS_ARCH32 1
  #endif
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

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <stdexcept>

// disable mingw-std-threads for mingw on MATLAB
#if (defined(MINGW) || defined(__MINGW32__)) && !defined(MATLAB_MEX_FILE )
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

#include "rang.hxx"
#include "Trace.hxx"
#include "Console.hxx"
#include "Malloc.hxx"
#include "Numbers.hxx"
#include "TicToc.hxx"
#include "ThreadPool.hxx"
// not used for the moment
//#include "quickpool.hxx"
#include "Quaternion.hxx"
#include "Table.hxx"
#include "Token.hxx"

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::string;
  #endif

  string basename( char const * const filename );

  template <typename T_int, typename T_real>
  void
  search_interval(
    T_int                npts,
    T_real const * const X,
    T_real             & x,
    T_int              & lastInterval,
    bool                 closed,
    bool                 can_extend
  );

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  extern template void search_interval(
    int32_t             npts,
    float const * const X,
    float             & x,
    int32_t           & lastInterval,
    bool                closed,
    bool                can_extend
  );

  extern template void search_interval(
    int32_t              npts,
    double const * const X,
    double             & x,
    int32_t            & lastInterval,
    bool                 closed,
    bool                 can_extend
  );

  extern template void search_interval(
    int64_t             npts,
    float const * const X,
    float             & x,
    int64_t           & lastInterval,
    bool                closed,
    bool                can_extend
  );

  extern template void search_interval(
    int64_t              npts,
    double const * const X,
    double             & x,
    int64_t            & lastInterval,
    bool                 closed,
    bool                 can_extend
  );
  #endif

  template <typename T_int, typename T_real>
  inline
  void
  searchInterval(
    T_int                npts,
    T_real const * const X,
    T_real             & x,
    T_int              & lastInterval,
    bool                 closed,
    bool                 can_extend
  ) {
    search_interval( npts, X, x, lastInterval, closed, can_extend );
  }

  static
  inline
  void
  to_upper( std::string & s ) {
    std::transform( s.begin(), s.end(), s.begin(), toupper );
  }

  static
  inline
  void
  to_lower( std::string & s ) {
    std::transform( s.begin(), s.end(), s.begin(), tolower );
  }

  static
  inline
  bool
  is_lower( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), islower );
  }

  static
  inline
  bool
  is_upper( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), isupper );
  }

  static
  inline
  bool
  is_alpha( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), isalpha );
  }

  static
  inline
  bool
  is_alphanum( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), isalnum );
  }

  static
  inline
  bool
  is_digits( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), isdigit );
  }

  static
  inline
  bool
  is_xdigits( std::string const & s ) {
    return std::all_of( s.begin(), s.end(), isxdigit );
  }

}

///
/// eof: Utils.hxx
///
