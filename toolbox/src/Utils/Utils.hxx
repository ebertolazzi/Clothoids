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

// check if compiler is C++11
#if ( defined( _MSC_VER ) && _MSC_VER >= 1800 ) || ( defined( __cplusplus ) && __cplusplus > 199711L )

#if ( defined( __cplusplus ) && __cplusplus <= 201103L )
#define UTILS_DEFAULT                                                                                                  \
  {                                                                                                                    \
  }
#else
#define UTILS_DEFAULT = default
#endif
#else
#error "Lapack Wrapper must be compiled using C++ >= C++11"
#endif

#define UTILS_PURE_VIRTUAL = 0
#ifndef UTILS_OS_WINDOWS
#define UTILS_OVERRIDE override
#define UTILS_CONSTEXPR constexpr
#else
#define UTILS_OVERRIDE
#define UTILS_CONSTEXPR
#endif

// STL
#include <algorithm>
#include <cassert>
#include <cctype>
#include <functional>  // For std::bind()
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>  // For std::remove_reference()
#include <utility>      // For std::move(), std::forward()
#include <vector>

// I/O
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

// C/C++
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>

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
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using string       = std::string;
  using string_view  = std::string_view;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;
#endif
}  // namespace Utils

#include "Malloc.hxx"
#include "rang.hxx"
#include "Console.hxx"
#include "Numbers.hxx"
#include "Quaternion.hxx"
#include "Table.hxx"
#include "TicToc.hxx"
#include "Token.hxx"

// order must be preserved
#include "ThreadPoolBase.hxx"
#include "ThreadUtils.hxx"
#include "ThreadPool0.hxx"
#include "ThreadPool1.hxx"
#include "ThreadPool2.hxx"
#include "ThreadPool3.hxx"
#include "ThreadPool4.hxx"
#include "ThreadPool5.hxx"
// -----------------------

namespace Utils
{

  /*\
  :|:  _
  :|: | |__  __ _ ___ ___ _ _  __ _ _ __  ___
  :|: | '_ \/ _` (_-</ -_) ' \/ _` | '  \/ -_)
  :|: |_.__/\__,_/__/\___|_||_\__,_|_|_|_\___|
  :|:
  \*/

  inline string
  get_basename( string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).parent_path().string();
  }

  inline string
  get_filename( string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).filename().string();
  }

  inline string
  get_extension( string_view path )
  {
    namespace fs = std::filesystem;
    return fs::path( path ).extension().string();
  }

  /*!
   * \addtogroup OS
   * @{
   */

  //! Fetches the value of an environment variable.
  /*!
   * This function retrieves the value of the environment variable with the name
   * specified by `ename` and stores it in the reference `res`.
   *
   * \param ename Name of the environment variable to retrieve.
   * \param res Reference to a string where the result will be stored.
   * \return True if the environment variable was found and its value retrieved,
   *         false otherwise.
   */
  bool get_environment( string_view ename, string & res );

  //! Sets the value of an environment variable.
  /*!
   * This function sets the value of the environment variable specified by
   * `ename` to `newval`. If the variable already exists, it will be overwritten
   * if the `overwrite` flag is true.
   *
   * \param ename Name of the environment variable to set.
   * \param newval The new value to set for the environment variable.
   * \param overwrite Flag to indicate whether the environment variable should
   * be overwritten if it already exists.
   */
  void set_environment( string_view ename, string_view newval, bool overwrite );

  //! Retrieves the MAC addresses of network interfaces.
  /*!
   * This function retrieves the MAC addresses for all available network
   * interfaces on the system and stores them in the provided map, with
   * interface names as the keys and MAC addresses as the values.
   *
   * \param addr A reference to a map where the MAC addresses will be stored.
   */
  void get_MAC_address( std::map<string, string> & addr );

  //! Retrieves the hostname of the system.
  /*!
   * \return A string containing the system's hostname.
   */
  string get_host_name();

  //! Fetches the IP addresses of the system.
  /*!
   * This function retrieves all IP addresses associated with the current
   * system's network interfaces and stores them in the provided vector `addr`.
   *
   * \param addr A reference to a vector where the IP addresses will be stored.
   */
  void get_IP_address( std::vector<string> & addr );

  //! Retrieves the username of the current user.
  /*!
   * \return A string containing the username of the current user.
   */
  string get_user_name();

  //! Retrieves the home directory of the current user.
  /*!
   * \return A string containing the home directory path of the current user.
   */
  string get_home_directory();

  //! Retrieves the full path to the current executable.
  /*!
   * \return A string containing the full path to the executable.
   */
  string get_executable_path_name();

  //! Checks if a file exists.
  /*!
   * \param fname The path to the file to check.
   * \return True if the file exists and is a regular file, false otherwise.
   */
  bool check_if_file_exists( string_view fname );

  //! Checks if a directory exists.
  /*!
   * \param dirname The path to the directory to check.
   * \return True if the directory exists and is valid, false otherwise.
   */
  bool check_if_dir_exists( string_view dirname );

  //! Creates a directory if it does not exist.
  /*!
   * This function creates a directory with the specified mode if it does not
   * already exist.
   *
   * \param dirname The path to the directory to create.
   * \param mode The permissions mode to set for the new directory.
   * \return True if the directory was created or already exists, false
   * otherwise.
   */
  bool make_directory( string_view dirname, unsigned mode = 0777 );

  string get_date();
  string get_day_time();
  string get_day_time_and_date();
  string get_log_date_time();

  /*! @} */  // End of OS group

  template <typename T_int, typename T_real>
  void search_interval( T_int npts, T_real const X[], T_real & x, T_int & last_interval, bool closed, bool can_extend );

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  extern template void search_interval( int32_t     npts,
                                        float const X[],
                                        float &     x,
                                        int32_t &   last_interval,
                                        bool        closed,
                                        bool        can_extend );

  extern template void search_interval( int32_t      npts,
                                        double const X[],
                                        double &     x,
                                        int32_t &    last_interval,
                                        bool         closed,
                                        bool         can_extend );

  extern template void search_interval( int64_t     npts,
                                        float const X[],
                                        float &     x,
                                        int64_t &   last_interval,
                                        bool        closed,
                                        bool        can_extend );

  extern template void search_interval( int64_t      npts,
                                        double const X[],
                                        double &     x,
                                        int64_t &    last_interval,
                                        bool         closed,
                                        bool         can_extend );
#endif

  template <typename T_int, typename T_real>
  inline void
  searchInterval( T_int npts, T_real const X[], T_real & x, T_int & last_interval, bool closed, bool can_extend )
  {
    search_interval( npts, X, x, last_interval, closed, can_extend );
  }

  static inline void
  to_upper( string & str )
  {
    std::transform( str.begin(), str.end(), str.begin(), []( unsigned char c ) -> unsigned char
                    { return static_cast<unsigned char>( std::toupper( c ) ); } );
  }
  static inline void
  to_lower( string & str )
  {
    std::transform( str.begin(), str.end(), str.begin(), []( unsigned char c ) -> unsigned char
                    { return static_cast<unsigned char>( std::tolower( c ) ); } );
  }
  static inline bool
  is_lower( string_view s )
  {
    return std::all_of( s.begin(), s.end(), islower );
  }
  static inline bool
  is_upper( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isupper );
  }
  static inline bool
  is_alpha( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isalpha );
  }
  static inline bool
  is_alphanum( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isalnum );
  }
  static inline bool
  is_digits( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isdigit );
  }
  static inline bool
  is_xdigits( string_view s )
  {
    return std::all_of( s.begin(), s.end(), isxdigit );
  }

  // https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
  static inline unsigned
  iLog2( uint32_t N )
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

  static inline unsigned
  iLog2( uint64_t N )
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

  string progress_bar( double progress, int width );
  string progress_bar2( double const progress, int const width, string_view const msg );

  void progress_bar( ostream_type &, double progress, int width, string_view msg );
  void progress_bar2( ostream_type &, double progress, int width, string_view msg );

}  // namespace Utils

//
// eof: Utils.hxx
//
