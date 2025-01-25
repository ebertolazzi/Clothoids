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

//
// file: TicToc.hxx
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_OS_WINDOWS
  #include <chrono>
  #include <thread>
#endif

#endif

namespace Utils {

  //!
  //! \brief A class for timing code execution.
  //!
  //! The `TicToc` class provides a simple and efficient way to measure
  //! the execution time of code segments. It uses high-resolution clocks
  //! for precise timing, allowing users to track performance and optimize
  //! code execution times.
  //!
  //! **Usage**
  //!
  //! To measure execution time, instantiate a `TicToc` object, call
  //! `tic()` to start timing, and `toc()` to stop timing.
  //! You can then retrieve the elapsed time in various units such
  //! as seconds, milliseconds, microseconds, or nanoseconds.
  //!
  //! **Example**
  //!
  //! \code{cpp}
  //! TicToc timer;
  //! timer.tic();
  //! // ... code to be timed ...
  //! timer.toc();
  //! std::cout << "Elapsed time: " << timer.elapsed_ms() << " ms" << std::endl;
  //! \endcode
  //!

#ifdef UTILS_OS_WINDOWS
  class TicToc {

    using real_type = double;
    int64_t   m_frequency;   // ticks per second
    int64_t   m_t1, m_t2;    // ticks
    real_type m_elapsed_time;

    TicToc( TicToc const & ) = delete;
    TicToc const & operator = ( TicToc const & ) const = delete;

  public:

    //!
    //! \brief Constructs a TicToc object and starts the timer.
    //!
    //! The constructor initializes the elapsed time to zero and calls
    //! the `tic()` function to start the timing.
    //!
    TicToc();
    ~TicToc() {}

    //!
    //! \brief Start timing.
    //!
    //! This function captures the current time point, marking the start of the timing.
    //!
    void tic();

    //!
    //! \brief End timing.
    //!
    //! This function captures the current time point, marking the end of the timing
    //! and calculates the elapsed time.
    //!
    void toc();

    //!
    //! \brief Return elapsed time (between tic-toc) in seconds.
    //!
    //! \return The elapsed time in seconds as a double.
    //!
    real_type elapsed_s()  const { return 1e-3*m_elapsed_time; }

    //!
    //! \brief Return elapsed time (between tic-toc) in milliseconds.
    //!
    //! \return The elapsed time in milliseconds as a double.
    //!
    real_type elapsed_ms() const { return m_elapsed_time; }

    //!
    //! \brief Return elapsed time (between tic-toc) in microseconds.
    //!
    //! \return The elapsed time in microseconds as a double.
    //!
    real_type elapsed_mus() const { return 1000*m_elapsed_time; }

    //!
    //! \brief Return elapsed time (between tic-toc) in nanoseconds.
    //!
    //! \return The elapsed time in nanoseconds as a double.
    //!
    real_type elapsed_ns() const { return 1e6*m_elapsed_time; }

  };

  //!
  //! \brief Sleep for a specified number of seconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of seconds.
  //!
  //! \param s The number of seconds to sleep.
  //!
  void sleep_for_seconds( unsigned s );

  //!
  //! \brief Sleep for a specified number of milliseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of milliseconds.
  //!
  //! \param ms The number of milliseconds to sleep.
  //!
  void sleep_for_milliseconds( unsigned ms );

  //!
  //! \brief Sleep for a specified number of microseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of microseconds.
  //!
  //! \param mus The number of microseconds to sleep.
  //!
  void sleep_for_microseconds( unsigned mus );

  //!
  //! \brief Sleep for a specified number of nanoseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of nanoseconds.
  //!
  //! \param ns The number of nanoseconds to sleep.
  //!
  void sleep_for_nanoseconds( unsigned ns );

#else
  class TicToc {

    using real_type          = double;
    using clock              = std::chrono::high_resolution_clock;
    using elapsed_resolution = std::chrono::microseconds;

    clock::time_point m_start_time;
    clock::time_point m_stop_time;

    elapsed_resolution m_elapsed_time{0};

   public:

    TicToc( TicToc const & ) = delete;
    TicToc const & operator = ( TicToc const & ) const = delete;

    //!
    //! \brief Constructs a TicToc object and starts the timer.
    //!
    //! The constructor initializes the elapsed time to zero and calls
    //! the `tic()` function to start the timing.
    //!
    TicToc() { this->tic(); }

    ~TicToc() = default;

    //!
    //! \brief Start timing.
    //!
    //! This function captures the current time point, marking the start of the timing.
    //!
    void
    tic()
    { m_start_time = clock::now(); }

    //!
    //! \brief End timing.
    //!
    //! This function captures the current time point, marking the end of the timing
    //! and calculates the elapsed time.
    //!
    void
    toc() {
      m_stop_time    = clock::now();
      m_elapsed_time = std::chrono::duration_cast<elapsed_resolution>(m_stop_time - m_start_time);
    }

    //!
    //! \brief Return elapsed time (between tic-toc) in seconds.
    //!
    //! \return The elapsed time in seconds as a double.
    //!
    real_type elapsed_s() const { return real_type(1e-6*m_elapsed_time.count()); }

    //!
    //! \brief Return elapsed time (between tic-toc) in milliseconds.
    //!
    //! \return The elapsed time in milliseconds as a double.
    //!
    real_type elapsed_ms() const { return real_type(1e-3*m_elapsed_time.count()); }

    //!
    //! \brief Return elapsed time (between tic-toc) in microseconds.
    //!
    //! \return The elapsed time in microseconds as a double.
    //!
    real_type elapsed_mus() const { return real_type(m_elapsed_time.count()); }

    //!
    //! \brief Return elapsed time (between tic-toc) in nanoseconds.
    //!
    //! \return The elapsed time in nanoseconds as a double.
    //!
    real_type elapsed_ns() const { return real_type(1e3*m_elapsed_time.count()); }
  };

  //!
  //! \brief Sleep for a specified number of seconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of seconds.
  //!
  //! \param s The number of seconds to sleep.
  //!
  inline
  void
  sleep_for_seconds( unsigned s )
  { std::this_thread::sleep_for(std::chrono::seconds(s)); }

  //!
  //! \brief Sleep for a specified number of milliseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of milliseconds.
  //!
  //! \param ms The number of milliseconds to sleep.
  //!
  inline
  void
  sleep_for_milliseconds( unsigned ms )
  {  std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }

  //!
  //! \brief Sleep for a specified number of microseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of microseconds.
  //!
  //! \param mus The number of microseconds to sleep.
  //!
  inline
  void
  sleep_for_microseconds( unsigned mus )
  {  std::this_thread::sleep_for(std::chrono::microseconds(mus)); }

  //!
  //! \brief Sleep for a specified number of nanoseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of nanoseconds.
  //!
  //! \param ns The number of nanoseconds to sleep.
  //!
  inline
  void
  sleep_for_nanoseconds( unsigned ns )
  {  std::this_thread::sleep_for(std::chrono::nanoseconds(ns)); }

#endif

}

//
// eof: TicToc.hxx
//
