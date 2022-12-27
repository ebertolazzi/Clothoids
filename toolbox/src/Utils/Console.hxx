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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Console.hxx
///

namespace Utils {

  using istream_type = std::basic_istream<char>;
  using ostream_type = std::basic_ostream<char>;
  using string       = std::string;

  string basename( char const * filename );

  class Console {

    mutable std::mutex m_message_mutex; // mutex for critical section

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  public:
    class Console_style {
    public:
      rang::style s;
      rang::fg    f;
      rang::bg    b;
    };
  #endif

  private:

    ostream_type * m_stream;

    // 0 only fatal, error
    // 1 + warning
    // 2
    int m_level;

    Console_style m_message_style = { rang::style::reset,     rang::fg::reset,  rang::bg::reset };
    Console_style m_warning_style = { rang::style::reset,     rang::fg::yellow, rang::bg::reset };
    Console_style m_error_style   = { rang::style::italic,    rang::fg::red,    rang::bg::reset };
    Console_style m_fatal_style   = { rang::style::underline, rang::fg::red,    rang::bg::reset };

  public:

    Console() = delete;
    Console( Console const & ) = delete;

    explicit
    Console( ostream_type * stream = &std::cout, int level = 4 );

    void change_level( int new_level );
    void changeLevel( int new_level ) { this->change_level(new_level); }

    void change_stream( ostream_type * new_stream ) { m_stream = new_stream; }
    void changeStream ( ostream_type * new_stream ) { m_stream = new_stream; }

    int  get_level() const { return m_level; }
    int  getLevel()  const { return m_level; }

    ostream_type * get_stream() const { return m_stream; }
    ostream_type * getStream()  const { return m_stream; }

    Console const & flush() const { m_stream->flush(); return *this; }

    Console const & message  ( char const * msg, int msg_level = 4 ) const;
    Console const & semaphore( unsigned ryg, char const * msg, int msg_level = 0 ) const;
    Console const & colors   ( unsigned c,   char const * msg, int msg_level = 0 ) const;
    Console const & warning  ( char const * msg ) const; // level >= 2
    Console const & error    ( char const * msg ) const; // level >= 1
    Console const & fatal    ( char const * msg ) const; // level >= 0

    Console const & black   ( char const * msg, int msg_level = 0 ) const;
    Console const & red     ( char const * msg, int msg_level = 0 ) const;
    Console const & green   ( char const * msg, int msg_level = 0 ) const;
    Console const & yellow  ( char const * msg, int msg_level = 0 ) const;
    Console const & blue    ( char const * msg, int msg_level = 0 ) const;
    Console const & magenta ( char const * msg, int msg_level = 0 ) const;
    Console const & cyan    ( char const * msg, int msg_level = 0 ) const;
    Console const & gray    ( char const * msg, int msg_level = 0 ) const;

    Console const & black_reversed   ( char const * msg, int msg_level = 0 ) const;
    Console const & red_reversed     ( char const * msg, int msg_level = 0 ) const;
    Console const & green_reversed   ( char const * msg, int msg_level = 0 ) const;
    Console const & yellow_reversed  ( char const * msg, int msg_level = 0 ) const;
    Console const & blue_reversed    ( char const * msg, int msg_level = 0 ) const;
    Console const & magenta_reversed ( char const * msg, int msg_level = 0 ) const;
    Console const & cyan_reversed    ( char const * msg, int msg_level = 0 ) const;
    Console const & gray_reversed    ( char const * msg, int msg_level = 0 ) const;

    Console const &
    message( string const & msg, int msg_level = 4 ) const
    { return this->message( msg.c_str(), msg_level); }

    Console const &
    semaphore( unsigned ryg, string const & msg, int msg_level = 0 ) const
    { return this->semaphore( ryg, msg.c_str(), msg_level); }

    Console const &
    colors( unsigned c, string const & msg, int msg_level = 0 ) const
    { return this->colors( c, msg.c_str(), msg_level); }

    Console const &
    warning( string const & msg ) const // level >= 2
    { return this->warning( msg.c_str() ); }

    Console const &
    error( string const & msg ) const // level >= 1
    { return this->error( msg.c_str() ); }

    Console const &
    fatal( string const & msg ) const // level >= 0
    { return this->fatal( msg.c_str() ); }

    Console const &
    black( string const & msg, int msg_level = 0 ) const
    { return this->black( msg.c_str(), msg_level ); }

    Console const &
    red( string const & msg, int msg_level = 0 ) const
    { return this->red( msg.c_str(), msg_level ); }

    Console const &
    green( string const & msg, int msg_level = 0 ) const
    { return this->green( msg.c_str(), msg_level ); }

    Console const &
    yellow( string const & msg, int msg_level = 0 ) const
    { return this->yellow( msg.c_str(), msg_level ); }

    Console const &
    blue( string const & msg, int msg_level = 0 ) const
    { return this->blue( msg.c_str(), msg_level ); }

    Console const &
    magenta( string const & msg, int msg_level = 0 ) const
    { return this->magenta( msg.c_str(), msg_level ); }

    Console const &
    cyan( string const & msg, int msg_level = 0 ) const
    { return this->cyan( msg.c_str(), msg_level ); }

    Console const &
    gray( string const & msg, int msg_level = 0 ) const
    { return this->gray( msg.c_str(), msg_level ); }

    Console const &
    black_reversed( string const & msg, int msg_level = 0 ) const
    { return this->black_reversed( msg.c_str(), msg_level ); }

    Console const &
    red_reversed( string const & msg, int msg_level = 0 ) const
    { return this->red_reversed( msg.c_str(), msg_level ); }

    Console const &
    green_reversed( string const & msg, int msg_level = 0 ) const
    { return this->green_reversed( msg.c_str(), msg_level ); }

    Console const &
    yellow_reversed( string const & msg, int msg_level = 0 ) const
    { return this->yellow_reversed( msg.c_str(), msg_level ); }

    Console const &
    blue_reversed( string const & msg, int msg_level = 0 ) const
    { return this->blue_reversed( msg.c_str(), msg_level ); }

    Console const &
    magenta_reversed( string const & msg, int msg_level = 0 ) const
    { return this->magenta_reversed( msg.c_str(), msg_level ); }

    Console const &
    cyan_reversed( string const & msg, int msg_level = 0 ) const
    { return this->cyan_reversed( msg.c_str(), msg_level ); }

    Console const &
    gray_reversed( string const & msg, int msg_level = 0 ) const
    { return this->gray_reversed( msg.c_str(), msg_level ); }

    void
    set_message_style(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      m_message_style.s = s; m_message_style.f = f; m_message_style.b = b;
    }

    void
    setMessageStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      this->set_message_style( s, f, b );
    }

    void
    set_warning_style(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      m_warning_style.s = s;
      m_warning_style.f = f;
      m_warning_style.b = b;
    }

    void
    setWarningStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      this->set_warning_style( s, f, b );
    }

    void
    set_error_style(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      m_error_style.s = s;
      m_error_style.f = f;
      m_error_style.b = b;
    }

    void
    setErrorStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      this->set_error_style( s, f, b );
    }

    void
    set_fatal_style(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      m_fatal_style.s = s;
      m_fatal_style.f = f;
      m_fatal_style.b = b;
    }

    void
    setFatalStyle(
      rang::style const & s,
      rang::fg    const & f,
      rang::bg    const & b
    ) {
      this->set_fatal_style( s, f, b );
    }

    //! set off coloring
    void
    set_off() const {
      #ifndef UTILS_OS_WINDOWS
      rang::setControlMode( rang::control::Off );
      #endif
    }

    //! set off coloring
    void setOff() const { this->set_off(); }

    //! set coloring automatic
    void
    set_auto() const {
      #ifdef UTILS_OS_WINDOWS
      rang::setWinTermMode( rang::winTerm::Auto );
      #else
      rang::setControlMode( rang::control::Auto );
      #endif
    }

    //! set coloring automatic
    void setAuto() const { this->set_auto(); }
  };

}

///
/// eof: Console.hxx
///
