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
// file: TicToc.cc
//

#include "Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifdef UTILS_OS_WINDOWS

#include <windows.h>

namespace Utils {

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
    m_t2           = t2.QuadPart;
    m_elapsed_time = (m_t2 - m_t1) * 1000.0 / m_frequency;
  }

  BOOLEAN
  nanosleep( LONGLONG ns100 ) {
    LARGE_INTEGER li;	   // Time defintion
    // Create timer
    HANDLE timer = CreateWaitableTimerW(NULL, TRUE, NULL);
    if ( !timer ) return FALSE;
    // Set timer properties
    li.QuadPart = -ns100;
    if ( !SetWaitableTimer(timer, &li, 0, NULL, NULL, FALSE) ) {
      CloseHandle(timer);
      return FALSE;
    }
    WaitForSingleObject(timer, INFINITE); // Start & wait for timer
    CloseHandle(timer);                   // Clean resources
    return TRUE;                          // Slept without problems
  }

  void
  sleep_for_seconds( unsigned s )
  { Sleep(DWORD(s) * 1000); }

  void
  sleep_for_milliseconds( unsigned ms )
  { Sleep(DWORD(ms)); }

  void
  sleep_for_microseconds( unsigned mus ) {
    nanosleep( LONGLONG(mus*10) );
  }

  void
  sleep_for_nanoseconds( unsigned ns ) {
    nanosleep( LONGLONG(ns/100) );
  }

}

#endif

#endif

//
// eof: TicToc.cc
//
