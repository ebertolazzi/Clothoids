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
 * \file Console.hxx
 * \brief Console utility for formatted output with different message levels.
 */

#ifndef UTILS_CONSOLE_HXX
#define UTILS_CONSOLE_HXX

namespace Utils
{

  using istream_type = std::basic_istream<char>;  //!< Type for input stream
  using ostream_type = std::basic_ostream<char>;  //!< Type for output stream
  using string       = std::string;               //!< Type for string
  using string_view  = std::string_view;          //!< Type for string view

  //! Get the base name of a file.
  /*!
   * \param filename C-style string representing the file name.
   * \return Base name of the file.
   */
  inline string basename( string_view filename )
  {
    size_t pos = filename.find_last_of( "/\\" );
    return ( pos == string::npos ) ? string( filename ) : string( filename.substr( pos + 1 ) );
  }

  //!
  //! \brief Class to handle console output with different styles and levels.
  //!
  //! \anchor console_level
  //! \note messages are printed only if the level priority is greather
  //!       than internal level stored in internal variable `m_level`.
  //!       The method `change_level` set `m_level` in the range from
  //!       -1 to 4. Each print message has a level, if level is less or
  //!       equal to `m_level` the message is printed.
  //!       Message with level=0 are high priority message that are
  //!       always printed unless `m_level` is set to -1 which suppress all
  //!       message. The default for `m_level` is 4 which means all messages are
  //!       printed.
  //!
  class Console
  {
    mutable std::mutex m_message_mutex;  //!< Mutex for critical section

  public:
    //! Structure representing console style.
    class Console_style
    {
    public:
      fmt::text_style ts;  //!< Text style using fmt::text_style
    };

  private:
    ostream_type * m_stream{ nullptr };  //!< Output stream pointer
    int            m_level{ 4 };         //!< Message level threshold

    Console_style m_message_style;  //!< Message style
    Console_style m_warning_style;  //!< Warning style
    Console_style m_error_style;    //!< Error style
    Console_style m_fatal_style;    //!< Fatal style

  public:
    //! Deleted default constructor.
    Console() = delete;

    //! Deleted copy constructor.
    Console( Console const & ) = delete;

    //! Constructor with stream and level parameters.
    /*!
     * \param stream Pointer to output stream (default is std::cout).
     * \param level Minimum message \ref console_level "level" to output
     * (default is 4).
     */
    explicit Console( ostream_type * stream = &std::cout, int level = 4 ) : m_stream( stream ), m_level( level )
    {
      // Initialize default styles using fmt
      m_message_style.ts = fmt::text_style();

      m_warning_style.ts = fmt::fg( fmt::color::yellow );

      m_error_style.ts = fmt::emphasis::italic | fmt::fg( fmt::color::red );

      m_fatal_style.ts = fmt::emphasis::underline | fmt::fg( fmt::color::red );
    }

    //! Destructor.
    ~Console() = default;

    //! Change the message level.
    /*!
     * \param new_level New level for message output.
     *
     * The admissible levels runs from -1 up to 4.
     * - Level -1 means that all messages are suppressed.
     * - Level 4 means that all messages are printed
     */
    void change_level( int new_level )
    {
      if ( new_level < -1 || new_level > 4 )
      {
        fmt::print( "Console::change_level( new_level = {})\nnew_level must be in the range [-1,4]\n", new_level );
        return;
      }
      m_level = new_level;
    }

    //! Change the message \ref console_level "level".
    /*!
     * \param new_level New level for message output.
     * \deprecated
     */
    void changeLevel( int new_level ) { this->change_level( new_level ); }

    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream.
     */
    void change_stream( ostream_type * new_stream ) { m_stream = new_stream; }
    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream.
     * \deprecated
     */
    void changeStream( ostream_type * new_stream ) { m_stream = new_stream; }

    //! Get the current message \ref console_level "level".
    /*!
     * \return Current message level.
     */
    int get_level() const { return m_level; }
    //! Get the current message level.
    /*!
     * \return Current message level.
     * \deprecated
     */
    int getLevel() const { return m_level; }

    //! Get the current output stream.
    /*!
     * \return Pointer to the current output stream.
     */
    ostream_type * get_stream() const { return m_stream; }
    //! Get the current output stream.
    /*!
     * \return Pointer to the current output stream.
     * \deprecated
     */
    ostream_type * getStream() const { return m_stream; }

    //! Flush the output stream.
    void flush() const { m_stream->flush(); }

    //! Output a message at a specified \ref console_level "level".
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 4).
     */
    void message( string_view const msg, int const msg_level = 4 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( m_message_style.ts, "{}", msg ); }
    }

    //! Output a semaphore message.
    /*!
     * \param ryg Semaphore indicator (red, yellow, green).
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void semaphore( unsigned const ryg, string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard             lock_access( m_message_mutex );
      static constexpr fmt::color ryg_color[3]{ fmt::color::red, fmt::color::yellow, fmt::color::green };
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( ryg_color[ryg % 3] ), "{}", msg ); }
    }

    //! Output a message with specified colors.
    /*!
     * \param c Color code.
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void colors( unsigned const c, string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard             lock_access( m_message_mutex );
      static constexpr fmt::color rvg_color[5]{ fmt::color::red,
                                                fmt::color::magenta,
                                                fmt::color::yellow,
                                                fmt::color::cyan,
                                                fmt::color::green };
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( rvg_color[c % 5] ), "{}", msg ); }
    }

    //! Output a warning message.
    /*!
     * \param msg The warning message to output.
     */
    void warning( string_view const msg ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level >= 2 ) { *m_stream << fmt::format( m_warning_style.ts, "{}", msg ); }
    }

    //! Output an error message.
    /*!
     * \param msg The error message to output.
     */
    void error( string_view const msg ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( m_level >= 1 ) { *m_stream << fmt::format( m_error_style.ts, "{}", msg ); }
    }

    //! Output a fatal message.
    /*!
     * \param msg The fatal message to output.
     */
    void fatal( string_view const msg ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      *m_stream << fmt::format( m_fatal_style.ts, "{}", msg );
    }

    //! Output a message in black color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void black( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::black ), "{}", msg ); }
    }

    //! Output a message in red color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void red( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::red ), "{}", msg ); }
    }

    //! Output a message in green color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void green( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::green ), "{}", msg ); }
    }

    //! Output a message in yellow color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void yellow( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::yellow ), "{}", msg ); }
    }

    //! Output a message in blue color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void blue( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::blue ), "{}", msg ); }
    }

    //! Output a message in magenta color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void magenta( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::magenta ), "{}", msg ); }
    }

    //! Output a message in cyan color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void cyan( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::cyan ), "{}", msg ); }
    }

    //! Output a message in gray color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void gray( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level ) { *m_stream << fmt::format( fmt::fg( fmt::color::gray ), "{}", msg ); }
    }

    //! Output a message in black reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void black_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::black ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in red reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void red_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::red ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in green reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void green_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::green ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in yellow reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void yellow_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::yellow ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in blue reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void blue_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::blue ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in magenta reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void magenta_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::magenta ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in cyan reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void cyan_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::cyan ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Output a message in gray reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void gray_reversed( string_view const msg, int const msg_level = 0 ) const
    {
      std::lock_guard lock_access( m_message_mutex );
      if ( msg_level <= m_level )
      {
        *m_stream << fmt::format( fmt::fg( fmt::color::gray ) | fmt::emphasis::reverse, "{}", msg );
      }
    }

    //! Sets the message style.
    //! \param s The text style (e.g., bold, underline).
    //! \param f The foreground color of the text.
    //! \param b The background color of the text.
    void set_message_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      m_message_style.ts = s | fmt::fg( f ) | fmt::bg( b );
    }

    //! Sets the message style using a higher-level function.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_message_style`
    void setMessageStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_message_style( s, f, b ); }

    //! Sets the warning style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_warning_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      m_warning_style.ts = s | fmt::fg( f ) | fmt::bg( b );
    }

    //! Sets the warning style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_warning_style`
    void setWarningStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_warning_style( s, f, b ); }

    //! Sets the error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_error_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      m_error_style.ts = s | fmt::fg( f ) | fmt::bg( b );
    }

    //! Sets the error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_error_style`
    void setErrorStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_error_style( s, f, b ); }

    //! Sets the fatal error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void set_fatal_style( fmt::emphasis s, fmt::color f, fmt::color b )
    {
      m_fatal_style.ts = s | fmt::fg( f ) | fmt::bg( b );
    }

    //! Sets the fatal error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_fatal_style`
    void setFatalStyle( fmt::emphasis s, fmt::color f, fmt::color b ) { this->set_fatal_style( s, f, b ); }

    //! Disables coloring.
    //! Turns off color control in non-Windows terminals.
    void set_off() const
    {
      // fmt automatically handles color support detection
      // No equivalent to rang::setControlMode needed
    }

    //! Disables coloring.
    //! Turns off color control in non-Windows terminals.
    //! \deprecated use `set_off`
    void setOff() const { this->set_off(); }

    //! Sets coloring to automatic mode.
    //! Enables automatic color control in terminals depending on the operating
    //! system.
    void set_auto() const
    {
      // fmt automatically handles color support detection
      // No equivalent to rang::setControlMode or rang::setWinTermMode needed
    }

    //! Sets coloring to automatic mode.
    //! Enables automatic color control in terminals depending on the operating
    //! system.
    //! \deprecated use `set_auto`
    void setAuto() const { this->set_auto(); }
  };

}  // namespace Utils

#endif  // UTILS_CONSOLE_HXX
