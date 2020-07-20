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

#ifndef TIC_TOC_HH
#define TIC_TOC_HH

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
  #define TIC_TOC_USE_CXX11
#endif

#ifdef G2LIB_OS_WINDOWS
  #include <windows.h>
  class TicToc {

    typedef double real_type;
    LARGE_INTEGER frequency;        // ticks per second
    LARGE_INTEGER t1, t2;           // ticks
    real_type elapsed_time;

    TicToc( TicToc const & );
    TicToc const & operator = ( TicToc const & ) const;

  public:

    TicToc()
    : elapsed_time(0)
    { QueryPerformanceFrequency(&frequency); tic(); }

    ~TicToc() {}

    void
    tic()
    { QueryPerformanceCounter(&t1); }

    void
    toc() {
      QueryPerformanceCounter(&t2);
      elapsed_time = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;;
    }

    real_type
    elapsed_s() const
    { return 1e-3*elapsed_time; }

    real_type
    elapsed_ms() const
    { return elapsed_time; }

  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { Sleep(DWORD(s)*1000); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  { Sleep(DWORD(ms)); }

#else

  #ifdef TIC_TOC_USE_CXX11

    #include <chrono>
    #include <thread>

    class TicToc {

      typedef double real_type;
      typedef std::chrono::high_resolution_clock clock;

      using elapsed_resolution = std::chrono::microseconds;

      clock::time_point start_time;
      clock::time_point stop_time;

      elapsed_resolution elapsed_time;

      TicToc( TicToc const & );
      TicToc const & operator = ( TicToc const & ) const;

    public:

      TicToc()
      : elapsed_time(0)
      { tic(); }

      ~TicToc() {}

      void
      tic()
      { start_time = clock::now(); }

      void
      toc() {
        stop_time    = clock::now();
        elapsed_time = std::chrono::duration_cast<elapsed_resolution>(stop_time - start_time);
      }

      real_type
      elapsed_s() const
      { return 1e-6*elapsed_time.count(); }

      real_type
      elapsed_ms() const
      { return 1e-3*elapsed_time.count(); }
    };

    inline
    void
    sleep_for_seconds( unsigned s )
    { std::this_thread::sleep_for(std::chrono::seconds(s)); }

    inline
    void
    sleep_for_milliseconds( unsigned ms )
    { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }

  #else

    #include <sys/time.h>
    #include <unistd.h>

    inline
    bool
    getTime( long & sec, long & usec ) {
      struct timeval now ;
      bool ok = gettimeofday(&now, NULL) == 0 ;
      if ( ok ) {
        sec  = now . tv_sec;
        usec = now . tv_usec;
      } else {
        sec = usec = 0 ;
      }
      return ok ;
    }

    class TicToc {

      typedef double real_type;
      long sec, usec ;
      real_type elapsed_time;

      TicToc( TicToc const & );
      TicToc const & operator = ( TicToc const & ) const;

    public:

      TicToc()
      : elapsed_time(0)
      { tic(); }

      ~TicToc() {}

      void
      tic()
      { getTime( sec, usec ); }

      void
      toc() {
        long new_sec, new_usec ;
        getTime( new_sec, new_usec ) ;
        elapsed_time = 1e3*(new_sec-sec)+1e-3*(new_usec-usec);
      }

      real_type
      elapsed_s() const
      { return 1e-3*elapsed_time; }

      real_type
      elapsed_ms() const
      { return elapsed_time; }
    };

    inline
    void
    sleep_for_seconds( unsigned s ) {
      useconds_t us = s; us *= 1000000;
      usleep(us);
    }

    inline
    void
    sleep_for_milliseconds( unsigned ms ) {
      useconds_t us = ms; us *= 1000;
      usleep(us);
    }

  #endif

#endif

#endif
