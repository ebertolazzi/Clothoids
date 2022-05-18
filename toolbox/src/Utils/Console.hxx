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

  typedef std::basic_istream<char> istream_type;
  typedef std::basic_ostream<char> ostream_type;

  std::string basename( char const * filename );

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
    void change_stream( ostream_type * new_stream );
    int  get_level() const { return m_level; }

    ostream_type * get_stream() const { return m_stream; }

    void changeLevel( int new_level ) { this->change_level(new_level); }
    void changeStream( ostream_type * new_stream ) { this->change_stream( new_stream ); }
    int  getLevel() const { return m_level; }

    ostream_type * getStream() const { return m_stream; }

    Console const & flush() const { m_stream->flush(); return *this; }

    Console const &
    message( std::string const & msg, int msg_level = 4 ) const;

    Console const &
    semaphore( unsigned ryg, std::string const & msg, int msg_level = 0 ) const;

    Console const &
    warning( std::string const & msg ) const; // level >= 2

    Console const &
    error( std::string const & msg ) const; // level >= 1

    Console const &
    fatal( std::string const & msg ) const; // level >= 0

    Console const & black   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & red     ( std::string const & msg, int msg_level = 0 ) const;
    Console const & green   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & yellow  ( std::string const & msg, int msg_level = 0 ) const;
    Console const & blue    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & magenta ( std::string const & msg, int msg_level = 0 ) const;
    Console const & cyan    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & gray    ( std::string const & msg, int msg_level = 0 ) const;

    Console const & black_reversed   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & red_reversed     ( std::string const & msg, int msg_level = 0 ) const;
    Console const & green_reversed   ( std::string const & msg, int msg_level = 0 ) const;
    Console const & yellow_reversed  ( std::string const & msg, int msg_level = 0 ) const;
    Console const & blue_reversed    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & magenta_reversed ( std::string const & msg, int msg_level = 0 ) const;
    Console const & cyan_reversed    ( std::string const & msg, int msg_level = 0 ) const;
    Console const & gray_reversed    ( std::string const & msg, int msg_level = 0 ) const;

    void
    set_message_style(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      m_message_style.s = s; m_message_style.f = f; m_message_style.b = b;
    }

    void
    setMessageStyle(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      this->set_message_style( s, f, b );
    }

    void
    set_warning_style(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      m_warning_style.s = s; m_warning_style.f = f; m_warning_style.b = b;
    }

    void
    setWarningStyle(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      this->set_warning_style( s, f, b );
    }

    void
    set_error_style(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      m_error_style.s = s; m_error_style.f = f; m_error_style.b = b;
    }

    void
    setErrorStyle(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      this->set_error_style( s, f, b );
    }

    void
    set_fatal_style(
      rang::style const & s, rang::fg const & f, rang::bg const & b
    ) {
      m_fatal_style.s = s; m_fatal_style.f = f; m_fatal_style.b = b;
    }

    void
    setFatalStyle(
      rang::style const & s, rang::fg const & f, rang::bg const & b
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
