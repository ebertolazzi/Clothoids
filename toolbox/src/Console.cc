/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
// file: Console.cc
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#if defined( __llvm__ ) || defined( __clang__ )
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

namespace Utils
{

  using std::lock_guard;
  using std::mutex;
  using std::string;

  void
  Console::change_level( int new_level )
  {
    UTILS_ASSERT( new_level >= -1 && new_level <= 4,
                  "Console::change_level( new_level = {})\n"
                  "new_level must be in the range [-1,4]\n",
                  new_level );
    m_level = new_level;
  }

  void
  Console::black( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::black << msg << rang::fg::reset;
  }

  void
  Console::red( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::red << msg << rang::fg::reset;
  }

  void
  Console::green( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::green << msg << rang::fg::reset;
  }

  void
  Console::yellow( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::yellow << msg << rang::fg::reset;
  }

  void
  Console::blue( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::blue << msg << rang::fg::reset;
  }

  void
  Console::magenta( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::magenta << msg << rang::fg::reset;
  }

  void
  Console::cyan( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::cyan << msg << rang::fg::reset;
  }

  void
  Console::gray( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level ) ( *m_stream ) << rang::fg::gray << msg << rang::fg::reset;
  }

  void
  Console::black_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::black << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::red_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::red << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::green_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::green << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::yellow_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::yellow << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::blue_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::blue << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::magenta_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::magenta << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::cyan_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::cyan << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  void
  Console::gray_reversed( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::fg::gray << rang::style::reversed << msg << rang::style::reset << rang::fg::reset;
  }

  Console::Console( ostream_type * stream, int const level ) : m_stream( stream ), m_level( level )
  {
    m_message_style.s = rang::style::reset;
    m_message_style.f = rang::fg::reset;
    m_message_style.b = rang::bg::reset;

    m_warning_style.s = rang::style::reset;
    m_warning_style.f = rang::fg::yellow;
    m_warning_style.b = rang::bg::reset;

    m_error_style.s = rang::style::italic;
    m_error_style.f = rang::fg::red;
    m_error_style.b = rang::bg::reset;

    m_fatal_style.s = rang::style::underline;
    m_fatal_style.f = rang::fg::red;
    m_fatal_style.b = rang::bg::reset;
  }

  void
  Console::semaphore( unsigned const ryg, string_view const msg, int const msg_level ) const
  {
    lock_guard const          lock_access( m_message_mutex );
    static constexpr rang::fg ryg_color[3]{ rang::fg::red, rang::fg::yellow, rang::fg::green };
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::style::reset << rang::bg::reset << ryg_color[ryg % 3] << msg << rang::fg::reset;
  }

  void
  Console::colors( unsigned const c, string_view const msg, int const msg_level ) const
  {
    lock_guard const          lock_access( m_message_mutex );
    static constexpr rang::fg rvg_color[5]{ rang::fg::red, rang::fg::magenta, rang::fg::yellow, rang::fg::cyan,
                                            rang::fg::green };
    if ( msg_level <= m_level )
      ( *m_stream ) << rang::style::reset << rang::bg::reset << rvg_color[c % 5] << msg << rang::fg::reset;
  }

  void
  Console::message( string_view const msg, int const msg_level ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( msg_level <= m_level )
      ( *m_stream ) << m_message_style.s << m_message_style.f << m_message_style.b << msg << rang::style::reset
                    << rang::fg::reset << rang::bg::reset;
  }

  void
  Console::warning( string_view const msg ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( m_level >= 2 )
      ( *m_stream ) << m_warning_style.s << m_warning_style.f << m_warning_style.b << msg << rang::style::reset
                    << rang::fg::reset << rang::bg::reset;
  }

  void
  Console::error( string_view const msg ) const
  {
    lock_guard const lock_access( m_message_mutex );
    if ( m_level >= 1 )
      ( *m_stream ) << m_error_style.s << m_error_style.f << m_error_style.b << msg << rang::style::reset
                    << rang::fg::reset << rang::bg::reset;
  }

  void
  Console::fatal( string_view const msg ) const
  {
    lock_guard const lock_access( m_message_mutex );
    ( *m_stream ) << m_fatal_style.s << m_fatal_style.f << m_fatal_style.b << msg << rang::style::reset
                  << rang::fg::reset << rang::bg::reset;
  }

}  // namespace Utils

#endif

//
// eof: Console.cc
//
