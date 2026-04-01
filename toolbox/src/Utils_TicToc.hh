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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: TicToc.hh
//

#pragma once
#ifndef UTILS_TICTOC_HXX
#define UTILS_TICTOC_HXX

#ifdef UTILS_OS_WINDOWS
#include <windows.h>
#else
#include <chrono>
#include <thread>
#endif

namespace Utils
{

  //!
  //! \brief A class for timing code execution with nanosecond precision.
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
  //! std::cout << "Elapsed time: " << timer.elapsed_ns() << " ns" << std::endl;
  //! std::cout << "Elapsed time: " << timer.elapsed_ms() << " ms" << std::endl;
  //! \endcode
  //!

#ifdef UTILS_OS_WINDOWS

  // Windows-specific helper function for high-resolution sleep
  inline BOOLEAN nanosleep( LONGLONG ns100 )
  {
    LARGE_INTEGER li;  // Time definition
    // Create timer
    HANDLE timer = CreateWaitableTimerW( NULL, TRUE, NULL );
    if ( !timer ) return FALSE;
    // Set timer properties
    li.QuadPart = -ns100;
    if ( !SetWaitableTimer( timer, &li, 0, NULL, NULL, FALSE ) )
    {
      CloseHandle( timer );
      return FALSE;
    }
    WaitForSingleObject( timer, INFINITE );  // Start & wait for timer
    CloseHandle( timer );                    // Clean resources
    return TRUE;                             // Slept without problems
  }

  class TicToc
  {
    using real_type = double;
    int64_t m_frequency{ 0 };      // ticks per second
    int64_t m_t1{ 0 }, m_t2{ 0 };  // ticks
    int64_t m_elapsed_ticks{ 0 };

    TicToc( TicToc const & )                         = delete;
    TicToc const & operator=( TicToc const & ) const = delete;

  public:
    //!
    //! \brief Constructs a TicToc object and starts the timer.
    //!
    //! The constructor initializes the elapsed time to zero and calls
    //! the `tic()` function to start the timing.
    //!
    TicToc() { tic(); }

    ~TicToc() {}

    //!
    //! \brief Initialize the timer frequency.
    //!
    //! This should be called once before using timing functions.
    //!
    static void init()
    {
      static bool initialized = false;
      if ( !initialized )
      {
        // QueryPerformanceFrequency is already called in constructor
        initialized = true;
      }
    }

    //!
    //! \brief Start timing.
    //!
    //! This function captures the current time point, marking the start of the
    //! timing.
    //!
    void tic()
    {
      LARGE_INTEGER t1;
      QueryPerformanceCounter( &t1 );
      m_t1 = t1.QuadPart;
    }

    //!
    //! \brief End timing.
    //!
    //! This function captures the current time point, marking the end of the
    //! timing and calculates the elapsed time in ticks.
    //!
    void toc()
    {
      LARGE_INTEGER t2;
      QueryPerformanceCounter( &t2 );
      m_t2            = t2.QuadPart;
      m_elapsed_ticks = m_t2 - m_t1;
    }

    //!
    //! \brief Return elapsed time in ticks.
    //!
    //! \return The elapsed time in ticks as int64_t.
    //!
    int64_t elapsed_ticks() const { return m_elapsed_ticks; }

    //!
    //! \brief Return elapsed time (between tic-toc) in seconds.
    //!
    //! \return The elapsed time in seconds as a double.
    //!
    real_type elapsed_s() const
    {
      LARGE_INTEGER frequency;
      QueryPerformanceFrequency( &frequency );
      return static_cast<real_type>( m_elapsed_ticks ) / static_cast<real_type>( frequency.QuadPart );
    }

    //!
    //! \brief Return elapsed time (between tic-toc) in milliseconds.
    //!
    //! \return The elapsed time in milliseconds as a double.
    //!
    real_type elapsed_ms() const { return elapsed_s() * 1000.0; }

    //!
    //! \brief Return elapsed time (between tic-toc) in microseconds.
    //!
    //! \return The elapsed time in microseconds as a double.
    //!
    real_type elapsed_mus() const { return elapsed_s() * 1000000.0; }

    //!
    //! \brief Return elapsed time (between tic-toc) in nanoseconds.
    //!
    //! \return The elapsed time in nanoseconds as a double.
    //!
    real_type elapsed_ns() const { return elapsed_s() * 1000000000.0; }
  };

  //!
  //! \brief Sleep for a specified number of seconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of seconds.
  //!
  //! \param s The number of seconds to sleep.
  //!
  inline void sleep_for_seconds( unsigned s )
  {
    Sleep( DWORD( s ) * 1000 );
  }

  //!
  //! \brief Sleep for a specified number of milliseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of milliseconds.
  //!
  //! \param ms The number of milliseconds to sleep.
  //!
  inline void sleep_for_milliseconds( unsigned ms )
  {
    Sleep( DWORD( ms ) );
  }

  //!
  //! \brief Sleep for a specified number of microseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of microseconds.
  //!
  //! \param mus The number of microseconds to sleep.
  //!
  inline void sleep_for_microseconds( unsigned mus )
  {
    nanosleep( LONGLONG( mus * 10 ) );
  }

  //!
  //! \brief Sleep for a specified number of nanoseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of nanoseconds.
  //!
  //! \param ns The number of nanoseconds to sleep.
  //!
  inline void sleep_for_nanoseconds( unsigned ns )
  {
    nanosleep( LONGLONG( ns / 100 ) );
  }

#else

  class TicToc
  {
    using real_type             = double;
    using clock                 = std::chrono::high_resolution_clock;
    using nanosecond_resolution = std::chrono::nanoseconds;

    clock::time_point m_start_time;
    clock::time_point m_stop_time;

    nanosecond_resolution m_elapsed_time{ 0 };

    TicToc( TicToc const & )                         = delete;
    TicToc const & operator=( TicToc const & ) const = delete;

  public:
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
    //! This function captures the current time point, marking the start of the
    //! timing.
    //!
    void tic() { m_start_time = clock::now(); }

    //!
    //! \brief End timing.
    //!
    //! This function captures the current time point, marking the end of the
    //! timing and calculates the elapsed time in nanoseconds.
    //!
    void toc()
    {
      m_stop_time    = clock::now();
      m_elapsed_time = std::chrono::duration_cast<nanosecond_resolution>( m_stop_time - m_start_time );
    }

    //!
    //! \brief Return elapsed time in nanoseconds as integer.
    //!
    //! \return The elapsed time in nanoseconds as int64_t.
    //!
    int64_t elapsed_ns_int() const { return m_elapsed_time.count(); }

    //!
    //! \brief Return elapsed time (between tic-toc) in seconds.
    //!
    //! \return The elapsed time in seconds as a double.
    //!
    real_type elapsed_s() const { return m_elapsed_time.count() * 1e-9; }

    //!
    //! \brief Return elapsed time (between tic-toc) in milliseconds.
    //!
    //! \return The elapsed time in milliseconds as a double.
    //!
    real_type elapsed_ms() const { return m_elapsed_time.count() * 1e-6; }

    //!
    //! \brief Return elapsed time (between tic-toc) in microseconds.
    //!
    //! \return The elapsed time in microseconds as a double.
    //!
    real_type elapsed_mus() const { return m_elapsed_time.count() * 1e-3; }

    //!
    //! \brief Return elapsed time (between tic-toc) in nanoseconds.
    //!
    //! \return The elapsed time in nanoseconds as a double.
    //!
    real_type elapsed_ns() const { return static_cast<real_type>( m_elapsed_time.count() ); }
  };

  //!
  //! \brief Sleep for a specified number of seconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of seconds.
  //!
  //! \param s The number of seconds to sleep.
  //!
  inline void sleep_for_seconds( unsigned s )
  {
    std::this_thread::sleep_for( std::chrono::seconds( s ) );
  }

  //!
  //! \brief Sleep for a specified number of milliseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of milliseconds.
  //!
  //! \param ms The number of milliseconds to sleep.
  //!
  inline void sleep_for_milliseconds( unsigned ms )
  {
    std::this_thread::sleep_for( std::chrono::milliseconds( ms ) );
  }

  //!
  //! \brief Sleep for a specified number of microseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of microseconds.
  //!
  //! \param mus The number of microseconds to sleep.
  //!
  inline void sleep_for_microseconds( unsigned mus )
  {
    std::this_thread::sleep_for( std::chrono::microseconds( mus ) );
  }

  //!
  //! \brief Sleep for a specified number of nanoseconds.
  //!
  //! This function pauses the execution of the current thread for the specified
  //! number of nanoseconds.
  //!
  //! \param ns The number of nanoseconds to sleep.
  //!
  inline void sleep_for_nanoseconds( unsigned ns )
  {
    std::this_thread::sleep_for( std::chrono::nanoseconds( ns ) );
  }

#endif  // UTILS_OS_WINDOWS

}  // namespace Utils

#endif  // UTILS_TICTOC_HXX

//
// eof: TicToc.hh
//
