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
/// file: TicToc.cc
///

#include "Utils.hh"

#ifdef UTILS_OS_WINDOWS

#include <windows.h>

namespace Utils {

  #define TOINT64(A) (static_cast<int64_t>(A.HighPart) << 32) | A.LowPart

  TicToc::TicToc() : m_elapsed_time(0) {
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    m_frequency = frequency.QuadPart;
    tic();
  }

  void
  TicToc::tic() {
    LARGE_INTEGER t1;
    QueryPerformanceCounter(&t1);
    m_t1 = t1.QuadPart;
  }

  void
  TicToc::toc() {
    LARGE_INTEGER t2;
    QueryPerformanceCounter(&t2);
    m_t2 = t2.QuadPart;
    m_elapsed_time = (m_t2 - m_t1) * 1000.0 / m_frequency;
    ;
  }

  void
  sleep_for_seconds( unsigned s ) {
    Sleep(DWORD(s) * 1000);
  }

  void
  sleep_for_milliseconds( unsigned ms ) {
    Sleep(DWORD(ms));
  }
}

#else

namespace Utils {

  void
  TicToc::tic()
  { m_start_time = clock::now(); }

  void
  TicToc::toc() {
    m_stop_time    = clock::now();
    m_elapsed_time = std::chrono::duration_cast<elapsed_resolution>(m_stop_time - m_start_time);
  }

  typename TicToc::real_type
  TicToc::elapsed_s() const
  { return 1e-6*m_elapsed_time.count(); }

  typename TicToc::real_type
  TicToc::elapsed_ms() const
  { return 1e-3*m_elapsed_time.count(); }
}

#endif

///
/// eof: TicToc.cc
///
