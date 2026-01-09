/*\

  Taken from: https://github.com/agauniyal/rang/blob/master/include/rang.hpp



\*/

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>

// System headers
#if defined( __APPLE__ ) || defined( __linux__ )
#include <unistd.h>
#elif defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#include <memory>
#include <windows.h>

// Only defined in windows 10 onwards, redefining in lower windows since it
// doesn't gets used in lower versions
// https://docs.microsoft.com/en-us/windows/console/getconsolemode
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#endif

namespace rang
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::ostream;
#endif

  /* For better compability with most of terminals do not use any style settings
   * except of reset, bold and reversed.
   * Note that on Windows terminals bold style is same as fgB color.
   */
  enum class style
  {
    reset     = 0,
    bold      = 1,
    dim       = 2,
    italic    = 3,
    underline = 4,
    blink     = 5,
    rblink    = 6,
    reversed  = 7,
    conceal   = 8,
    crossed   = 9
  };

  enum class fg
  {
    black   = 30,
    red     = 31,
    green   = 32,
    yellow  = 33,
    blue    = 34,
    magenta = 35,
    cyan    = 36,
    gray    = 37,
    reset   = 39
  };

  enum class bg
  {
    black   = 40,
    red     = 41,
    green   = 42,
    yellow  = 43,
    blue    = 44,
    magenta = 45,
    cyan    = 46,
    gray    = 47,
    reset   = 49
  };

  enum class fgB
  {
    black   = 90,
    red     = 91,
    green   = 92,
    yellow  = 93,
    blue    = 94,
    magenta = 95,
    cyan    = 96,
    gray    = 97
  };

  enum class bgB
  {
    black   = 100,
    red     = 101,
    green   = 102,
    yellow  = 103,
    blue    = 104,
    magenta = 105,
    cyan    = 106,
    gray    = 107
  };

  enum class control
  {             // Behaviour of rang function calls
    Off   = 0,  // toggle off rang style/color calls
    Auto  = 1,  // (Default) autodect terminal and colorize if needed
    Force = 2   // force ansi color output to non terminal streams
  };
  // Use rang::setControlMode to set rang control mode

  enum class winTerm
  {              // Windows Terminal Mode
    Auto   = 0,  // (Default) automatically detects wheter Ansi or Native API
    Ansi   = 1,  // Force use Ansi API
    Native = 2   // Force use Native API
  };
  // Use rang::setWinTermMode to explicitly set terminal API for Windows
  // Calling rang::setWinTermMode have no effect on other OS

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  namespace rang_implementation
  {

    using std::any_of;
    using std::atomic;
    using std::begin;
    using std::cerr;
    using std::clog;
    using std::cout;
    using std::enable_if;
    using std::end;
    using std::getenv;
    using std::is_same;
    using std::nothrow;
    using std::ostream;
    using std::streambuf;
    using std::string;
    using std::strstr;
    using std::unique_ptr;
    using std::wstring;

    inline std::atomic<control> & controlMode() noexcept
    {
      static std::atomic<control> value( control::Auto );
      return value;
    }

    inline std::atomic<winTerm> & winTermMode() noexcept
    {
      static std::atomic<winTerm> termMode( winTerm::Auto );
      return termMode;
    }

    template <typename T> using enableStd = typename enable_if<
      is_same<T, rang::style>::value || is_same<T, rang::fg>::value || is_same<T, rang::bg>::value ||
        is_same<T, rang::fgB>::value || is_same<T, rang::bgB>::value,
      ostream &>::type;

    inline bool supportsColor() noexcept
    {
#if defined( __APPLE__ ) || defined( __linux__ )
      static const bool result = []
      {
        char const * terms[] = { "ansi",  "color", "console", "cygwin", "gnome",  "konsole", "kterm",
                                 "linux", "msys",  "putty",   "rxvt",   "screen", "vt100",   "xterm" };
        char const * env_p   = getenv( "TERM" );
        if ( env_p == nullptr ) return false;
        return any_of(
          begin( terms ),
          end( terms ),
          [&]( const char * term ) { return strstr( env_p, term ) != nullptr; } );
      }();
#elif defined( _WIN32 ) || defined( _WIN64 )
      // All windows versions support colors through native console methods
      static constexpr bool result = true;
#else
      static constexpr bool result = false;
#endif
      return result;
    }

#ifdef _WIN32

    inline bool isMsysPty( int fd ) noexcept
    {
// disable for mingw on MATLAB
#if ( defined( MINGW ) || defined( __MINGW32__ ) ) && defined( MATLAB_MEX_FILE )
      return false;
#else
      // Dynamic load for binary compability with old Windows
      auto const ptrGetFileInformationByHandleEx = reinterpret_cast<decltype( &GetFileInformationByHandleEx )>(
        GetProcAddress( GetModuleHandle( TEXT( "kernel32.dll" ) ), "GetFileInformationByHandleEx" ) );

      if ( !ptrGetFileInformationByHandleEx ) return false;

      HANDLE h = reinterpret_cast<HANDLE>( _get_osfhandle( fd ) );
      if ( h == INVALID_HANDLE_VALUE ) return false;

      // Check that it's a pipe:
      if ( GetFileType( h ) != FILE_TYPE_PIPE ) return false;

      // POD type is binary compatible with FILE_NAME_INFO from WinBase.h
      // It have the same alignment and used to avoid UB in caller code
      struct MY_FILE_NAME_INFO
      {
        DWORD FileNameLength;
        WCHAR FileName[MAX_PATH];
      };

      auto pNameInfo = unique_ptr<MY_FILE_NAME_INFO>( new ( nothrow ) MY_FILE_NAME_INFO() );
      if ( !pNameInfo ) return false;

      // Check pipe name is template of
      // {"cygwin-","msys-"}XXXXXXXXXXXXXXX-ptyX-XX
      if ( !ptrGetFileInformationByHandleEx( h, FileNameInfo, pNameInfo.get(), sizeof( MY_FILE_NAME_INFO ) ) )
        return false;

      wstring name( pNameInfo->FileName, pNameInfo->FileNameLength / sizeof( WCHAR ) );

      if (
        ( name.find( L"msys-" ) == wstring::npos && name.find( L"cygwin-" ) == wstring::npos ) ||
        name.find( L"-pty" ) == wstring::npos )
        return false;

      return true;
#endif
    }

#endif

    inline bool isTerminal( streambuf const * osbuf ) noexcept
    {
#if defined( __APPLE__ ) || defined( __linux__ )
      if ( osbuf == cout.rdbuf() )
      {
        static const bool cout_term = isatty( fileno( stdout ) ) != 0;
        return cout_term;
      }
      if ( osbuf == cerr.rdbuf() || osbuf == clog.rdbuf() )
      {
        static const bool cerr_term = isatty( fileno( stderr ) ) != 0;
        return cerr_term;
      }
#elif defined( _WIN32 ) || defined( _WIN64 )
      if ( osbuf == cout.rdbuf() )
      {
        static const bool cout_term = ( _isatty( _fileno( stdout ) ) || isMsysPty( _fileno( stdout ) ) );
        return cout_term;
      }
      if ( osbuf == cerr.rdbuf() || osbuf == clog.rdbuf() )
      {
        static const bool cerr_term = ( _isatty( _fileno( stderr ) ) || isMsysPty( _fileno( stderr ) ) );
        return cerr_term;
      }
#endif
      return false;
    }

#ifdef _WIN32

    struct SGR
    {                     // Select Graphic Rendition parameters for Windows console
      BYTE    fgColor;    // foreground color (0-15) lower 3 rgb bits + intense bit
      BYTE    bgColor;    // background color (0-15) lower 3 rgb bits + intense bit
      BYTE    bold;       // emulated as FOREGROUND_INTENSITY bit
      BYTE    underline;  // emulated as BACKGROUND_INTENSITY bit
      BOOLEAN inverse;    // swap foreground/bold & background/underline
      BOOLEAN conceal;    // set foreground/bold to background/underline
    };

    enum class AttrColor : BYTE
    {  // Color attributes for console screen buffer
      black   = 0,
      red     = 4,
      green   = 2,
      yellow  = 6,
      blue    = 1,
      magenta = 5,
      cyan    = 3,
      gray    = 7
    };

    inline HANDLE getConsoleHandle( streambuf const * osbuf ) noexcept
    {
      if ( osbuf == cout.rdbuf() )
      {
        static const HANDLE hStdout = GetStdHandle( STD_OUTPUT_HANDLE );
        return hStdout;
      }
      else if ( osbuf == cerr.rdbuf() || osbuf == clog.rdbuf() )
      {
        static const HANDLE hStderr = GetStdHandle( STD_ERROR_HANDLE );
        return hStderr;
      }
      return INVALID_HANDLE_VALUE;
    }

    inline bool setWinTermAnsiColors( streambuf const * osbuf ) noexcept
    {
      HANDLE h = getConsoleHandle( osbuf );
      if ( h == INVALID_HANDLE_VALUE ) return false;
      DWORD dwMode = 0;
      if ( !GetConsoleMode( h, &dwMode ) ) return false;
      dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
      if ( !SetConsoleMode( h, dwMode ) ) return false;
      return true;
    }

    inline bool supportsAnsi( streambuf const * osbuf ) noexcept
    {
      if ( osbuf == cout.rdbuf() )
      {
        static const bool cout_ansi = ( isMsysPty( _fileno( stdout ) ) || setWinTermAnsiColors( osbuf ) );
        return cout_ansi;
      }
      else if ( osbuf == cerr.rdbuf() || osbuf == clog.rdbuf() )
      {
        static const bool cerr_ansi = ( isMsysPty( _fileno( stderr ) ) || setWinTermAnsiColors( osbuf ) );
        return cerr_ansi;
      }
      return false;
    }

    inline SGR const & defaultState() noexcept
    {
      static const SGR defaultSgr = []() -> SGR
      {
        CONSOLE_SCREEN_BUFFER_INFO info;
        WORD                       attrib = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;
        if (
          GetConsoleScreenBufferInfo( GetStdHandle( STD_OUTPUT_HANDLE ), &info ) ||
          GetConsoleScreenBufferInfo( GetStdHandle( STD_ERROR_HANDLE ), &info ) )
        {
          attrib = info.wAttributes;
        }
        SGR sgr     = { 0, 0, 0, 0, FALSE, FALSE };
        sgr.fgColor = attrib & 0x0F;
        sgr.bgColor = ( attrib & 0xF0 ) >> 4;
        return sgr;
      }();
      return defaultSgr;
    }

    inline BYTE ansi2attr( BYTE rgb ) noexcept
    {
      static const AttrColor rev[8] = { AttrColor::black, AttrColor::red,     AttrColor::green, AttrColor::yellow,
                                        AttrColor::blue,  AttrColor::magenta, AttrColor::cyan,  AttrColor::gray };
      return static_cast<BYTE>( rev[rgb] );
    }

    inline void setWinSGR( rang::bg col, SGR & state ) noexcept
    {
      if ( col != rang::bg::reset ) { state.bgColor = ansi2attr( static_cast<BYTE>( col ) - 40 ); }
      else
      {
        state.bgColor = defaultState().bgColor;
      }
    }

    inline void setWinSGR( rang::fg col, SGR & state ) noexcept
    {
      if ( col != rang::fg::reset ) { state.fgColor = ansi2attr( static_cast<BYTE>( col ) - 30 ); }
      else
      {
        state.fgColor = defaultState().fgColor;
      }
    }

    inline void setWinSGR( rang::bgB col, SGR & state ) noexcept
    {
      state.bgColor = ( BACKGROUND_INTENSITY >> 4 ) | ansi2attr( static_cast<BYTE>( col ) - 100 );
    }

    inline void setWinSGR( rang::fgB col, SGR & state ) noexcept
    {
      state.fgColor = FOREGROUND_INTENSITY | ansi2attr( static_cast<BYTE>( col ) - 90 );
    }

    inline void setWinSGR( rang::style style, SGR & state ) noexcept
    {
      switch ( style )
      {
        case rang::style::reset: state = defaultState(); break;
        case rang::style::bold: state.bold = FOREGROUND_INTENSITY; break;
        case rang::style::underline:
        case rang::style::blink: state.underline = BACKGROUND_INTENSITY; break;
        case rang::style::reversed: state.inverse = TRUE; break;
        case rang::style::conceal: state.conceal = TRUE; break;
        default: break;
      }
    }

    inline SGR & current_state() noexcept
    {
      static SGR state = defaultState();
      return state;
    }

    inline WORD SGR2Attr( const SGR & state ) noexcept
    {
      WORD attrib = 0;
      if ( state.conceal )
      {
        if ( state.inverse )
        {
          attrib = ( state.fgColor << 4 ) | state.fgColor;
          if ( state.bold ) attrib |= FOREGROUND_INTENSITY | BACKGROUND_INTENSITY;
        }
        else
        {
          attrib = ( state.bgColor << 4 ) | state.bgColor;
          if ( state.underline ) attrib |= FOREGROUND_INTENSITY | BACKGROUND_INTENSITY;
        }
      }
      else if ( state.inverse )
      {
        attrib = ( state.fgColor << 4 ) | state.bgColor;
        if ( state.bold ) attrib |= BACKGROUND_INTENSITY;
        if ( state.underline ) attrib |= FOREGROUND_INTENSITY;
      }
      else
      {
        attrib = state.fgColor | ( state.bgColor << 4 ) | state.bold | state.underline;
      }
      return attrib;
    }

    template <typename T> inline void setWinColorNative_tmpl( ostream & os, T const value )
    {
      HANDLE const h = getConsoleHandle( os.rdbuf() );
      if ( h != INVALID_HANDLE_VALUE )
      {
        setWinSGR( value, current_state() );
        // Out all buffered text to console with previous settings:
        os.flush();
        SetConsoleTextAttribute( h, SGR2Attr( current_state() ) );
      }
    }

#define DO_INSTANTIATION( A )                            \
  inline void setWinColorNative( ostream & os, A value ) \
  {                                                      \
    setWinColorNative_tmpl<A>( os, value );              \
  }

    DO_INSTANTIATION( enum rang::style );
    DO_INSTANTIATION( enum rang::bg );
    DO_INSTANTIATION( enum rang::fg );
    DO_INSTANTIATION( enum rang::bgB );
    DO_INSTANTIATION( enum rang::fgB );

#undef DO_INSTANTIATION

#endif

#ifdef _WIN32
    template <typename T> inline void setWinColorAnsi( ostream & os, T const value )
    {
      os << "\033[" << static_cast<int>( value ) << "m";
    }

    template <typename T> inline enableStd<T> setColor( ostream & os, T const value )
    {
      if ( winTermMode() == winTerm::Auto )
      {
        if ( supportsAnsi( os.rdbuf() ) ) { setWinColorAnsi( os, value ); }
        else
        {
          setWinColorNative( os, value );
        }
      }
      else if ( winTermMode() == winTerm::Ansi ) { setWinColorAnsi( os, value ); }
      else
      {
        setWinColorNative( os, value );
      }
      return os;
    }
#else
    template <typename T> inline enableStd<T> setColor( ostream & os, T const value )
    {
      return os << "\033[" << static_cast<int>( value ) << "m";
    }
#endif

  }  // namespace rang_implementation

  template <typename T> inline rang_implementation::enableStd<T> operator<<( ostream & os, T const value )
  {
    control const option = rang_implementation::controlMode();
    switch ( option )
    {
      case control::Auto:
        return rang_implementation::supportsColor() && rang_implementation::isTerminal( os.rdbuf() )
                 ? rang_implementation::setColor( os, value )
                 : os;
      case control::Force: return rang_implementation::setColor( os, value );
      default: return os;
    }
  }

  inline void setWinTermMode( rang::winTerm const value ) noexcept
  {
    rang_implementation::winTermMode() = value;
  }

  inline void setControlMode( control const value ) noexcept
  {
    rang_implementation::controlMode() = value;
  }
#endif

}  // namespace rang

//
// eof: rang.hh
//
