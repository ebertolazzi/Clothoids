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
/// file: ThreadUtils.cc
///

#include "Utils.hh"

namespace Utils {

  WorkerLoop::WorkerLoop()
  : m_active(true)
  , m_running(false) // block
  , m_do_job(false)
  {
    m_job            = []()->void{}; // empty job
    m_running_thread = std::thread( &WorkerLoop::worker_loop, this );
  }

  WorkerLoop::~WorkerLoop() {
    m_active  = false;
    m_running = true;
    m_do_job  = true; // force loop exit
    m_cv.notify_one();
    m_running_thread.join();
  }

  void
  WorkerLoop::worker_loop() {
    while ( m_active ) {
      std::unique_lock<std::mutex> lk(m_mutex);
      m_cv.wait(lk, [this]{ return this->m_do_job;} );
      if ( !m_active ) break;
      m_running = true;
      m_job();
      m_do_job  = false;
      m_running = false;
      m_cv.notify_one();
    }
  }

  void
  WorkerLoop::exec( std::function<void()> & fun ) {
    std::unique_lock<std::mutex> lk(m_mutex);
    m_cv.wait(lk, [this]{ return !this->m_do_job;}  ); // another job
    m_cv.wait(lk, [this]{ return !this->m_running;} ); // still running
    m_job    = fun;
    m_do_job = true;
    m_cv.notify_one();
  }

  void
  WorkerLoop::exec() {
    std::unique_lock<std::mutex> lk(m_mutex);
    m_cv.wait(lk, [this]{ return !this->m_do_job;}  ); // another job
    m_cv.wait(lk, [this]{ return !this->m_running;} ); // still running
    m_do_job = true;
    m_cv.notify_one();
  }

  void
  WorkerLoop::wait() {
    std::unique_lock<std::mutex> lk(m_mutex);
    m_cv.wait(lk, [this]{ return !this->m_do_job;} );
    m_cv.notify_one();
  }

}

///
/// eof: ThreadUtils.cc
///
