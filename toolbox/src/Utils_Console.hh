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

/**
 * \file Utils_Console.hh
 * \brief Console utility for formatted output with different message levels.
 */

#ifndef UTILS_CONSOLE_HXX
#define UTILS_CONSOLE_HXX

#include "Utils_fmt.hh"

#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace Utils
{

  //! Get the base name of a file.
  /*!
   * \param filename C-style string representing the file name.
   * \return Base name of the file.
   */
  [[nodiscard]] inline constexpr string basename( string_view filename )
  {
    size_t pos = filename.find_last_of( "/\\" );
    return ( pos == string::npos ) ? string( filename ) : string( filename.substr( pos + 1 ) );
  }

  //!
  //! \brief Class to handle console output with different styles and levels.
  //!
  //! \anchor console_level
  //! \note `change_level` sets the threshold in the range [-1,4]. A message is
  //!       printed when `msg_level <= m_level`. Level 0 is the highest-priority
  //!       regular level, while level -1 suppresses every message. The default
  //!       threshold is 4, which prints all valid message levels.
  //!       `fatal` messages are always printed regardless of their level,
  //!       with the sole exception of `m_level == -1` (all messages suppressed).
  //!
  class Console
  {
  public:
    //! Legacy wrapper retained for source compatibility.
    //! \deprecated prefer `fmt::text_style` directly.
    class Console_style
    {
    public:
      fmt::text_style ts;
    };
    
    using integer = int;

  private:
    mutable std::mutex m_message_mutex;  //!< Mutex protecting m_stream/m_level and all output operations

    ostream_type * m_stream        = &std::cout; //!< Output stream pointer
    integer        m_level         = 4;          //!< Message level threshold
    bool           m_color_enabled = true;       //!< If false, styling is stripped and plain text is written instead

    //! Text style used for a given message category.
    fmt::text_style m_message_style;  //!< Message style
    fmt::text_style m_warning_style;  //!< Warning style
    fmt::text_style m_error_style;    //!< Error style
    fmt::text_style m_fatal_style;    //!< Fatal style

    [[nodiscard]] static ostream_type * validate_stream( ostream_type * stream )
    {
      Utils::Warning( stream != nullptr, "Utils::Console: output stream must not be null, revert to std::cout" );
      return stream != nullptr ? stream : &std::cout;
    }

    [[nodiscard]] static integer validate_level( integer level )
    {
      Utils::Assert( level >= -1 && level <= 4, "Utils::Console: level={} must be in the range [-1,4]", level );
      return level;
    }

    //! Writes styled text while ensuring that ANSI reset codes precede every
    //! newline. This prevents a background color from leaking into the next
    //! terminal line. The caller must hold `m_message_mutex`.
    void write_styled_unlocked( fmt::text_style const & ts, string_view const msg ) const
    {
      if ( !m_color_enabled )
      {
        *m_stream << msg;
        return;
      }

      size_t line_begin = 0;
      while ( line_begin < msg.size() )
      {
        size_t const newline = msg.find( '\n', line_begin );
        size_t const line_end = newline == string_view::npos ? msg.size() : newline;
        string_view const line = msg.substr( line_begin, line_end - line_begin );

        // fmt appends ESC[0m after the styled fragment. By excluding '\n'
        // from the fragment, the reset is emitted before moving to next line.
        if ( !line.empty() ) { *m_stream << fmt::format( ts, "{}", line ); }

        if ( newline == string_view::npos ) { break; }
        *m_stream << '\n';
        line_begin = newline + 1;
      }
    }

    //! Common implementation shared by every "leveled" print method
    //! (message/semaphore/colors/warning/error and all named-color helpers).
    //! Locks the mutex, checks the level, and writes the styled message.
    void print_styled( fmt::text_style const & ts, string_view const msg, integer const msg_level ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level >= 0 && msg_level <= m_level )
      {
        write_styled_unlocked( ts, msg );
      }
    }

    template <typename... Args>
    void print_styled(
      fmt::text_style const &       ts,
      integer const                 msg_level,
      fmt::format_string<Args...>   fmt_msg,
      Args &&...                    args
    ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level >= 0 && msg_level <= m_level )
      {
        string const msg = fmt::format( fmt_msg, std::forward<Args>( args )... );
        write_styled_unlocked( ts, msg );
      }
    }

    //! Common implementation for `fatal`: always printed, the only exception
    //! being `m_level == -1` (which suppresses every message, fatal included).
    void print_always( fmt::text_style const & ts, string_view const msg ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level > -1 )
      {
        write_styled_unlocked( ts, msg );
      }
    }

    template <typename... Args>
    void print_always(
      fmt::text_style const &       ts,
      fmt::format_string<Args...>   fmt_msg,
      Args &&...                    args
    ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level >= 0 )
      {
        string const msg = fmt::format( fmt_msg, std::forward<Args>( args )... );
        write_styled_unlocked( ts, msg );
      }
    }

  public:
    //! Deleted default constructor.
    Console() = delete;

    // Rule of five: this class owns a std::mutex (non copyable, non movable),
    // so copy/move are made explicitly unavailable instead of relying on the
    // implicit deletion caused by the mutex member.
    Console( Console const & )             = delete;  //!< Deleted copy constructor.
    Console & operator=( Console const & ) = delete;   //!< Deleted copy assignment.
    Console( Console && )                  = delete;   //!< Deleted move constructor.
    Console & operator=( Console && )      = delete;   //!< Deleted move assignment.

    //! Constructor with stream and level parameters.
    /*!
     * \param stream Pointer to output stream (default is std::cout).
     * \param level Minimum message \ref console_level "level" to output
     * (default is 4).
     * \throws std::invalid_argument if `stream` is null.
     * \throws std::out_of_range if `level` is outside [-1,4].
     */
    explicit Console( ostream_type * stream = &std::cout, integer level = 4 )
    : m_stream( validate_stream( stream ) )
    , m_level( validate_level( level ) )
    {
      // Initialize default styles using fmt
      m_message_style = fmt::text_style();
      m_warning_style = fmt::fg( fmt::color::yellow );
      m_error_style   = fmt::emphasis::italic | fmt::fg( fmt::color::red );
      m_fatal_style   = fmt::emphasis::underline | fmt::fg( fmt::color::red );
    }

    //! Destructor.
    ~Console() = default;

    //! Change the message level.
    /*!
     * \param new_level New level for message output.
     *
     * The admissible levels run from -1 up to 4.
     * - Level -1 means that all messages are suppressed.
     * - Level 4 means that all messages are printed.
     * \throws std::out_of_range if `new_level` is outside [-1,4].
     */
    void change_level( integer new_level )
    {
      new_level = validate_level( new_level );
      std::lock_guard lock_access( m_message_mutex );
      m_level = new_level;
    }

    //! Change the message \ref console_level "level".
    /*!
     * \param new_level New level for message output.
     * \deprecated
     */
    void changeLevel( integer new_level ) { this->change_level( new_level ); }

    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream; must not be null.
     * \throws std::invalid_argument if `new_stream` is null.
     */
    void change_stream( ostream_type * new_stream )
    {
      ostream_type * checked_stream = validate_stream( new_stream );
      std::lock_guard lock_access( m_message_mutex );
      m_stream = checked_stream;
    }
    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream.
     * \deprecated
     */
    void changeStream( ostream_type * new_stream ) { this->change_stream( new_stream ); }

    //! Get the current message \ref console_level "level".
    /*!
     * \return Current message level.
     */
    [[nodiscard]] integer get_level() const
    {
      std::lock_guard lock_access( m_message_mutex );
      return m_level;
    }
    //! Get the current message level.
    /*!
     * \return Current message level.
     * \deprecated
     */
    [[nodiscard]] integer getLevel() const { return this->get_level(); }

    //! Get the current output stream.
    /*!
     * \return A snapshot of the current output-stream pointer.
     * \note The caller owns the stream and must keep it alive. Another thread
     * may change the configured stream immediately after this method returns.
     */
    [[nodiscard]] ostream_type * get_stream() const
    {
      std::lock_guard lock_access( m_message_mutex );
      return m_stream;
    }
    //! Get the current output stream.
    /*!
     * \return Pointer to the current output stream.
     * \deprecated
     */
    [[nodiscard]] ostream_type * getStream() const { return this->get_stream(); }

    //! Flush the output stream.
    void flush() const
    {
      std::lock_guard lock_access( m_message_mutex );
      m_stream->flush();
    }

    //! Output a message at a specified \ref console_level "level".
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 4).
     */
    void message( string_view const msg, integer const msg_level = 4 ) const { print_styled( m_message_style, msg, msg_level ); }

    //! Output a formatted message at a specified \ref console_level "level".
    /*!
     * \param msg_level The \ref console_level "level of the message".
     * \param fmt_msg   The fmt format string.
     * \param args      Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void message( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( m_message_style, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a semaphore message.
    /*!
     * \param ryg Semaphore indicator (0=red, 1=yellow, 2=green; wraps modulo 3).
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void semaphore( integer const ryg, string_view const msg, integer const msg_level = 0 ) const
    {
      static constexpr fmt::color ryg_color[3]{ fmt::color::red, fmt::color::yellow, fmt::color::green };
      print_styled( fmt::fg( ryg_color[ryg % 3] ), msg, msg_level );
    }

    //! Output a formatted semaphore message at a specified \ref console_level "level".
    /*!
     * \param msg_level The \ref console_level "level of the message".
     * \param ryg       Semaphore indicator (red, yellow, green).
     * \param fmt_msg   The fmt format string.
     * \param args      Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void semaphore( integer const msg_level, integer const ryg, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      static constexpr fmt::color ryg_color[3]{ fmt::color::red, fmt::color::yellow, fmt::color::green };
      print_styled( fmt::fg( ryg_color[ryg % 3] ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message with specified colors.
    /*!
     * \param c Color code (0=red,1=magenta,2=yellow,3=cyan,4=green; wraps modulo 5).
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void colors( integer const c, string_view const msg, integer const msg_level = 0 ) const
    {
      static constexpr fmt::color rvg_color[5]{
        fmt::color::red, fmt::color::magenta, fmt::color::yellow, fmt::color::cyan, fmt::color::green
      };
      print_styled( fmt::fg( rvg_color[c % 5] ), msg, msg_level );
    }

    //! Output a formatted message with specified colors at a specified \ref console_level "level".
    /*!
     * \param msg_level The \ref console_level "level of the message".
     * \param c         Color code.
     * \param fmt_msg   The fmt format string.
     * \param args      Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void colors( integer const msg_level, integer const c, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      static constexpr fmt::color rvg_color[5]{
        fmt::color::red, fmt::color::magenta, fmt::color::yellow, fmt::color::cyan, fmt::color::green
      };
      print_styled( fmt::fg( rvg_color[c % 5] ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a warning message. Printed whenever `m_level >= 2`.
    /*!
     * \param msg The warning message to output.
     */
    void warning( string_view const msg ) const { print_styled( m_warning_style, msg, 2 ); }

    //! Output a formatted warning message.
    /*!
     * \param fmt_msg The fmt format string.
     * \param args    Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void warning( fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( m_warning_style, 2, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output an error message. Printed whenever `m_level >= 1`.
    /*!
     * \param msg The error message to output.
     */
    void error( string_view const msg ) const { print_styled( m_error_style, msg, 1 ); }

    //! Output a formatted error message.
    /*!
     * \param fmt_msg The fmt format string.
     * \param args    Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void error( fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( m_error_style, 1, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a fatal message. Always printed, unless `m_level == -1`.
    /*!
     * \param msg The fatal message to output.
     */
    void fatal( string_view const msg ) const { print_always( m_fatal_style, msg ); }

    //! Output a formatted fatal message. Always printed, unless `m_level == -1`.
    /*!
     * \param fmt_msg The fmt format string.
     * \param args    Arguments passed to `fmt::format`.
     */
    template <typename... Args>
    void fatal( fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_always( m_fatal_style, fmt_msg, std::forward<Args>( args )... );
    }

    // ------------------------------------------------------------------
    // Named-color helpers.
    // Every XXX/XXX_reversed pair below is a thin wrapper around
    // print_styled(): the locking/level-check logic lives in a single
    // place, so it can no longer drift out of sync between colors
    // (this is what previously caused `fatal()` to forget the level
    // check that every other method performed).
    // ------------------------------------------------------------------

    //! Output a message in black color.
    void black( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::black ), msg, msg_level );
    }
    //! Output a formatted message in black color at a specified \ref console_level "level".
    template <typename... Args>
    void black( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::black ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in black reversed color.
    void black_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::black ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in black reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void black_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::black ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in red color.
    void red( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::red ), msg, msg_level );
    }
    //! Output a formatted message in red color at a specified \ref console_level "level".
    template <typename... Args>
    void red( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::red ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in red reversed color.
    void red_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::red ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in red reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void red_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::red ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in green color.
    void green( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::green ), msg, msg_level );
    }
    //! Output a formatted message in green color at a specified \ref console_level "level".
    template <typename... Args>
    void green( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::green ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in green reversed color.
    void green_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::green ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in green reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void green_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::green ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in yellow color.
    void yellow( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::yellow ), msg, msg_level );
    }
    //! Output a formatted message in yellow color at a specified \ref console_level "level".
    template <typename... Args>
    void yellow( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::yellow ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in yellow reversed color.
    void yellow_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::yellow ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in yellow reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void yellow_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::yellow ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in blue color.
    void blue( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::blue ), msg, msg_level );
    }
    //! Output a formatted message in blue color at a specified \ref console_level "level".
    template <typename... Args>
    void blue( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::blue ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in blue reversed color.
    void blue_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::blue ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in blue reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void blue_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::blue ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in magenta color.
    void magenta( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::magenta ), msg, msg_level );
    }
    //! Output a formatted message in magenta color at a specified \ref console_level "level".
    template <typename... Args>
    void magenta( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::magenta ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in magenta reversed color.
    void magenta_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::magenta ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in magenta reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void magenta_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::magenta ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in cyan color.
    void cyan( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::cyan ), msg, msg_level );
    }
    //! Output a formatted message in cyan color at a specified \ref console_level "level".
    template <typename... Args>
    void cyan( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::cyan ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in cyan reversed color.
    void cyan_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::cyan ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in cyan reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void cyan_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::cyan ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in gray color.
    void gray( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::gray ), msg, msg_level );
    }
    //! Output a formatted message in gray color at a specified \ref console_level "level".
    template <typename... Args>
    void gray( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::gray ), msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    //! Output a message in gray reversed color.
    void gray_reversed( string_view const msg, integer const msg_level = 0 ) const
    {
      print_styled( fmt::fg( fmt::color::gray ) | fmt::emphasis::reverse, msg, msg_level );
    }
    //! Output a formatted message in gray reversed color at a specified \ref console_level "level".
    template <typename... Args>
    void gray_reversed( integer const msg_level, fmt::format_string<Args...> fmt_msg, Args &&... args ) const
    {
      print_styled( fmt::fg( fmt::color::gray ) | fmt::emphasis::reverse, msg_level, fmt_msg, std::forward<Args>( args )... );
    }

    // ------------------------------------------------------------------
    // Style setters
    // ------------------------------------------------------------------

    //! Sets the message style.
    //! \param s The text style (e.g., bold, underline).
    //! \param f The foreground color of the text.
    //! \param b The background color of the text.
    void set_message_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      std::lock_guard lock_access( m_message_mutex );
      m_message_style = s | fmt::fg( f ) | fmt::bg( b );
    }
    //! \deprecated use `set_message_style`
    void setMessageStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_message_style( s, f, b ); }

    //! Sets the warning style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_warning_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      std::lock_guard lock_access( m_message_mutex );
      m_warning_style = s | fmt::fg( f ) | fmt::bg( b );
    }
    //! \deprecated use `set_warning_style`
    void setWarningStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_warning_style( s, f, b ); }

    //! Sets the error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_error_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      std::lock_guard lock_access( m_message_mutex );
      m_error_style = s | fmt::fg( f ) | fmt::bg( b );
    }
    //! \deprecated use `set_error_style`
    void setErrorStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_error_style( s, f, b ); }

    //! Sets the fatal error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_fatal_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      std::lock_guard lock_access( m_message_mutex );
      m_fatal_style = s | fmt::fg( f ) | fmt::bg( b );
    }
    //! \deprecated use `set_fatal_style`
    void setFatalStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_fatal_style( s, f, b ); }

    //! Enables or disables ANSI styling explicitly.
    void set_color_enabled( bool enabled )
    {
      std::lock_guard lock_access( m_message_mutex );
      m_color_enabled = enabled;
    }

    //! Reports whether ANSI styling is enabled.
    [[nodiscard]] bool color_enabled() const
    {
      std::lock_guard lock_access( m_message_mutex );
      return m_color_enabled;
    }

    //! Disables coloring unconditionally.
    /*!
     * Forces every subsequent message to be written as plain text, with no
     * ANSI escape codes at all -- even if the destination stream is a
     * color-capable terminal. This is a real, stateful switch (unlike a
     * "best effort" detection): once `set_off()` is called, styling is
     * stripped by `print_styled`/`print_always` regardless of what the
     * terminal could actually display.
     */
    void set_off()
    {
      this->set_color_enabled( false );
    }

    //! Disables coloring unconditionally.
    //! \deprecated use `set_off`
    void setOff() { this->set_off(); }

    //! Restores coloring (the default state set by the constructor).
    /*!
     * Re-enables emission of ANSI escape codes for every subsequent message.
     * \note this class does not perform terminal-capability detection
     * (e.g. `isatty`): "auto" here means "styling is emitted normally",
     * i.e. it cancels a previous `set_off()`. If the destination stream is
     * later redirected to a file/pipe, ANSI codes will still be written to
     * it; use `set_off()` explicitly when writing to a non-terminal sink.
     */
    void set_auto()
    {
      this->set_color_enabled( true );
    }

    //! Restores coloring (the default state set by the constructor).
    //! \deprecated use `set_auto`
    void setAuto() { this->set_auto(); }
  };

}  // namespace Utils

#endif  // UTILS_CONSOLE_HXX
