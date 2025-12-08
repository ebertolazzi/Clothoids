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

//
// file: ThreadUtils.cc
//

#if defined( __llvm__ ) || defined( __clang__ )
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#include "Utils.hh"
#include "Utils_fmt.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace Utils
{

  WorkerLoop::WorkerLoop()
    : m_active( true )
    , m_running( false )  // block
    , m_do_job( false )
    , m_job( []() -> void {} )
  {
    m_running_thread = std::thread( &WorkerLoop::worker_loop, this );
  }

  WorkerLoop::~WorkerLoop()
  {
    m_active  = false;
    m_running = true;
    m_do_job  = true;  // force loop exit
    m_cv.notify_one();
    m_running_thread.join();
  }

  void
  WorkerLoop::worker_loop()
  {
    while ( m_active )
    {
      std::unique_lock lk( m_mutex );
      m_cv.wait( lk, [this] { return this->m_do_job; } );
      if ( !m_active ) break;
      m_running = true;
      m_job();
      m_do_job  = false;
      m_running = false;
      m_cv.notify_one();
    }
  }

  void
  WorkerLoop::exec( std::function<void()> & fun )
  {
    std::unique_lock lk( m_mutex );
    m_cv.wait( lk, [this] { return !this->m_do_job; } );   // another job
    m_cv.wait( lk, [this] { return !this->m_running; } );  // still running
    m_job    = fun;
    m_do_job = true;
    m_cv.notify_one();
  }

  void
  WorkerLoop::exec()
  {
    std::unique_lock lk( m_mutex );
    m_cv.wait( lk, [this] { return !this->m_do_job; } );   // another job
    m_cv.wait( lk, [this] { return !this->m_running; } );  // still running
    m_do_job = true;
    m_cv.notify_one();
  }

  void
  WorkerLoop::wait()
  {
    std::unique_lock lk( m_mutex );
    m_cv.wait( lk, [this] { return !this->m_do_job; } );
    m_cv.notify_one();
  }

#ifdef UTILS_OS_WINDOWS
  WinMutex::WinMutex() : m_mutex( NULL )
  {
    m_mutex = CreateMutex( NULL,   // no security descriptor
                           FALSE,  // mutex not owned
                           NULL    // object name
    );
    UTILS_ASSERT( m_mutex != NULL, "WinMutex(): error: {}.\n", GetLastError() );
  }

  void
  WinMutex::lock()
  {
    DWORD res = WaitForSingleObject( m_mutex, INFINITE );
    UTILS_ASSERT0( res == WAIT_OBJECT_0, "WinMutex::lock, WAIT_TIMEOUT" );
  }

  void
  WinMutex::unlock()
  {
    DWORD res = ReleaseMutex( m_mutex );
    UTILS_ASSERT0( res == WAIT_OBJECT_0, "WinMutex::lock, WAIT_TIMEOUT" );
  }
#endif

}  // namespace Utils

#endif

//
// eof: ThreadUtils.cc
//
