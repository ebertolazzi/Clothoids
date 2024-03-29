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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: TicToc.hxx
///

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_OS_WINDOWS
  #include <chrono>
  #include <thread>
#endif
#endif

namespace Utils {

#ifdef UTILS_OS_WINDOWS
  //!
  //! Class for timing code execution.
  //!
  class TicToc {

    using real_type = double;
    int64_t   m_frequency;   // ticks per second
    int64_t   m_t1, m_t2;    // ticks
    real_type m_elapsed_time;

    TicToc( TicToc const & ) = delete;
    TicToc const & operator = ( TicToc const & ) const = delete;

  public:

    TicToc();
    ~TicToc() {}

    //!
    //! Start timing
    //!
    void tic();

    //!
    //! End timing
    //!
    void toc();

    //!
    //! Return elapsed time (between tic-toc) in seconds
    //!
    real_type elapsed_s()  const { return 1e-3*m_elapsed_time; }

    //!
    //! Return elapsed time (between tic-toc) in milliseconds
    //!
    real_type elapsed_ms() const { return m_elapsed_time; }

    //!
    //! Return elapsed time (between tic-toc) in microseconds
    //!
    real_type elapsed_mus() const { return 1000*m_elapsed_time; }

    //!
    //! Return elapsed time (between tic-toc) in nanoseconds
    //!
    real_type elapsed_ns() const { return 1e6*m_elapsed_time; }

  };

  void sleep_for_seconds( unsigned s );
  void sleep_for_milliseconds( unsigned ms );
  void sleep_for_microseconds( unsigned mus );
  void sleep_for_nanoseconds( unsigned ns );

#else

  //!
  //! Class for timing code execution.
  //!
  class TicToc {

    using real_type = double;
    #ifdef TIC_TOC_USE_HRC
    using clock = std::chrono::high_resolution_clock;
    #else
    using clock = std::chrono::steady_clock;
    #endif

    using elapsed_resolution = std::chrono::microseconds;

    clock::time_point m_start_time;
    clock::time_point m_stop_time;

    elapsed_resolution m_elapsed_time;

   public:

    TicToc( TicToc const & ) = delete;
    TicToc const & operator = ( TicToc const & ) const = delete;

    TicToc()
    : m_elapsed_time(0)
    { this->tic(); }

    ~TicToc() UTILS_DEFAULT;

    //!
    //! Start timing
    //!
    void tic();

    //!
    //! End timing
    //!
    void toc();

    //!
    //! Return elapsed time (between tic-toc) in seconds
    //!
    real_type elapsed_s() const;

    //!
    //! Return elapsed time (between tic-toc) in milliseconds
    //!
    real_type elapsed_ms() const;

    //!
    //! Return elapsed time (between tic-toc) in microseconds
    //!
    real_type elapsed_mus() const;

    //!
    //! Return elapsed time (between tic-toc) in nanoseconds
    //!
    real_type elapsed_ns() const;
  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { std::this_thread::sleep_for(std::chrono::seconds(s)); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  {  std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }

  inline
  void
  sleep_for_microseconds( unsigned mus )
  {  std::this_thread::sleep_for(std::chrono::microseconds(mus)); }

  inline
  void
  sleep_for_nanoseconds( unsigned ns )
  {  std::this_thread::sleep_for(std::chrono::nanoseconds(ns)); }

#endif

}

///
/// eof: TicToc.hxx
///
