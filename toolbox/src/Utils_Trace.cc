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
// file: Trace.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include "Utils_trace.hh"

#include "Utils.hh"
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

  const char *
  Runtime_Error::what() const noexcept
  {
    return runtime_error::what();
  }

  const char *
  Runtime_TraceError::what() const noexcept
  {
    return runtime_error::what();
  }

#ifdef UTILS_OS_WINDOWS

  void
  print_trace( int const line, string_view const file, string_view const msg, ostream_type & stream )
  {
    fmt::print( stream,
                "---------------------------------------------------------\n"
                "file: {}:{}\n{}\n"
                "---------------------------------------------------------\n",
                file, line, msg );
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

  string
  Runtime_TraceError::grab_backtrace( string_view const reason, string_view const file, int const line )
  {
    return fmt::format( "\n{}\nOn File:{}:{}\n", reason, file, line );
  }

#else

#ifdef UTILS_OS_OSX
  static std::string
  addr_to_line( std::string const & binary_path, void * addr )
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
  static string
  demang( string_view const mangled_name )
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

  //! print a trace stack used in debug
  void
  print_trace( int const line, string_view const file, string_view const reason, ostream_type & stream )
  {
    fmt::print( stream,
                "\n{}\nOn File:{}:{}\nprocess ID:{}, parent process "
                "ID:{}\nstack trace:\n",
                reason, Utils::get_basename( file ), line, getpid(), getppid() );

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

  string
  Runtime_TraceError::grab_backtrace( string_view const reason, string_view const file, int const line )
  {
    ostringstream ost;
    print_trace( line, file, reason, ost );
    return ost.str();
  }
#endif

}  // namespace Utils

#endif

//
// eof: Trace.cc
//
