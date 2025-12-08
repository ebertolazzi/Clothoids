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

namespace Utils
{

  using istream_type = std::basic_istream<char>;  //!< Type for input stream
  using ostream_type = std::basic_ostream<char>;  //!< Type for output stream
  using string       = std::string;               //!< Type for string

  //! Get the base name of a file.
  /*!
   * \param filename C-style string representing the file name.
   * \return Base name of the file.
   */
  string basename( string_view filename );

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  public:
    //! Structure representing console style.
    class Console_style
    {
    public:
      rang::style s;  //!< Text style
      rang::fg    f;  //!< Foreground color
      rang::bg    b;  //!< Background color
    };
#endif

  private:
    ostream_type * m_stream{ nullptr };  //!< Output stream pointer
    int            m_level{ 4 };         //!< Message level threshold

    Console_style m_message_style = { rang::style::reset, rang::fg::reset, rang::bg::reset };    //!< Message style
    Console_style m_warning_style = { rang::style::reset, rang::fg::yellow, rang::bg::reset };   //!< Warning style
    Console_style m_error_style   = { rang::style::italic, rang::fg::red, rang::bg::reset };     //!< Error style
    Console_style m_fatal_style   = { rang::style::underline, rang::fg::red, rang::bg::reset };  //!< Fatal style

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
    explicit Console( ostream_type * stream = &std::cout, int level = 4 );

    //! Destructor.
    ~Console() {}

    //! Change the message level.
    /*!
     * \param new_level New level for message output.
     *
     * The admissible levels runs from -1 up to 4.
     * - Level -1 means that all messages are suppressed.
     * - Level 4 means that all messages are printed
     */
    void change_level( int new_level );
    //! Change the message \ref console_level "level".
    /*!
     * \param new_level New level for message output.
     * \deprecated
     */
    void
    changeLevel( int new_level )
    {
      this->change_level( new_level );
    }

    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream.
     */
    void
    change_stream( ostream_type * new_stream )
    {
      m_stream = new_stream;
    }
    //! Change the output stream.
    /*!
     * \param new_stream Pointer to the new output stream.
     * \deprecated
     */
    void
    changeStream( ostream_type * new_stream )
    {
      m_stream = new_stream;
    }

    //! Get the current message \ref console_level "level".
    /*!
     * \return Current message level.
     */
    int
    get_level() const
    {
      return m_level;
    }
    //! Get the current message level.
    /*!
     * \return Current message level.
     * \deprecated
     */
    int
    getLevel() const
    {
      return m_level;
    }

    //! Get the current output stream.
    /*!
     * \return Pointer to the current output stream.
     */
    ostream_type *
    get_stream() const
    {
      return m_stream;
    }
    //! Get the current output stream.
    /*!
     * \return Pointer to the current output stream.
     * \deprecated
     */
    ostream_type *
    getStream() const
    {
      return m_stream;
    }

    //! Flush the output stream.
    void
    flush() const
    {
      m_stream->flush();
    }

    //! Output a message at a specified \ref console_level "level".
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 4).
     */
    void message( string_view msg, int msg_level = 4 ) const;

    //! Output a semaphore message.
    /*!
     * \param ryg Semaphore indicator (red, yellow, green).
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void semaphore( unsigned ryg, string_view msg, int msg_level = 0 ) const;

    //! Output a message with specified colors.
    /*!
     * \param c Color code.
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void colors( unsigned c, string_view msg, int msg_level = 0 ) const;

    //! Output a warning message.
    /*!
     * \param msg The warning message to output.
     */
    void warning( string_view msg ) const;  // level >= 2

    //! Output an error message.
    /*!
     * \param msg The error message to output.
     */
    void error( string_view msg ) const;  // level >= 1

    //! Output a fatal message.
    /*!
     * \param msg The fatal message to output.
     */
    void fatal( string_view msg ) const;  // level >= 0

    //! Output a message in black color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void black( string_view msg, int msg_level = 0 ) const;

    //! Output a message in red color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void red( string_view msg, int msg_level = 0 ) const;

    //! Output a message in green color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void green( string_view msg, int msg_level = 0 ) const;

    //! Output a message in yellow color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void yellow( string_view msg, int msg_level = 0 ) const;

    //! Output a message in blue color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void blue( string_view msg, int msg_level = 0 ) const;

    //! Output a message in magenta color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void magenta( string_view msg, int msg_level = 0 ) const;

    //! Output a message in cyan color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void cyan( string_view msg, int msg_level = 0 ) const;

    //! Output a message in gray color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void gray( string_view msg, int msg_level = 0 ) const;

    //! Output a message in black reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void black_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in red reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void red_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in green reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void green_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in yellow reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void yellow_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in blue reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void blue_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in magenta reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void magenta_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in cyan reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void cyan_reversed( string_view msg, int msg_level = 0 ) const;

    //! Output a message in gray reversed color.
    /*!
     * \param msg       The message to output.
     * \param msg_level The \ref console_level "level of the message" (default
     * is 0).
     */
    void gray_reversed( string_view msg, int msg_level = 0 ) const;

    //! Sets the message style.
    //! \param s The text style (e.g., bold, underline).
    //! \param f The foreground color of the text.
    //! \param b The background color of the text.
    void
    set_message_style( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      m_message_style.s = s;
      m_message_style.f = f;
      m_message_style.b = b;
    }

    //! Sets the message style using a higher-level function.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_message_style`
    void
    setMessageStyle( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      this->set_message_style( s, f, b );
    }

    //! Sets the warning style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void
    set_warning_style( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      m_warning_style.s = s;
      m_warning_style.f = f;
      m_warning_style.b = b;
    }

    //! Sets the warning style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_warning_style`
    void
    setWarningStyle( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      this->set_warning_style( s, f, b );
    }

    //! Sets the error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void
    set_error_style( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      m_error_style.s = s;
      m_error_style.f = f;
      m_error_style.b = b;
    }

    //! Sets the error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_error_style`
    void
    setErrorStyle( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      this->set_error_style( s, f, b );
    }

    //! Sets the fatal error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    void
    set_fatal_style( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      m_fatal_style.s = s;
      m_fatal_style.f = f;
      m_fatal_style.b = b;
    }

    //! Sets the fatal error style.
    //! \param s The text style.
    //! \param f The foreground color.
    //! \param b The background color.
    //! \deprecated use `set_fatal_style`
    void
    setFatalStyle( rang::style const & s, rang::fg const & f, rang::bg const & b )
    {
      this->set_fatal_style( s, f, b );
    }

    //! Disables coloring.
    //! Turns off color control in non-Windows terminals.
    void
    set_off() const
    {
#ifndef UTILS_OS_WINDOWS
      rang::setControlMode( rang::control::Off );
#endif
    }

    //! Disables coloring.
    //! Turns off color control in non-Windows terminals.
    //! \deprecated use `set_off`
    void
    setOff() const
    {
      this->set_off();
    }

    //! Sets coloring to automatic mode.
    //! Enables automatic color control in terminals depending on the operating
    //! system.
    void
    set_auto() const
    {
#ifdef UTILS_OS_WINDOWS
      rang::setWinTermMode( rang::winTerm::Auto );
#else
      rang::setControlMode( rang::control::Auto );
#endif
    }

    //! Sets coloring to automatic mode.
    //! Enables automatic color control in terminals depending on the operating
    //! system.
    //! \deprecated use `set_auto`
    void
    setAuto() const
    {
      this->set_auto();
    }
  };

}  // namespace Utils
