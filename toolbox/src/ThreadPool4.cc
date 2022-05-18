/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
/// file: ThreadPool4.cc
///

#include "Utils.hh"

namespace Utils {

  ThreadPool4::ThreadPool4( unsigned thread_count, unsigned queue_capacity )
  : m_done(false)
  , m_running_task(0)
  , m_running_thread(0)
  , m_work_queue( queue_capacity == 0 ? 4 * (thread_count+1) : queue_capacity )
  , m_pop_waiting(0)
  , m_push_waiting(0)
  {
    create_workers( thread_count );
  }

  void
  ThreadPool4::push_task( TaskData * task ) {
    m_tm.tic();
    // --------------------------
    ++m_push_waiting;
    m_work_on_queue_mutex.lock();
    m_queue_push_cv.wait( m_work_on_queue_mutex, [&]()->bool { return !m_work_queue.is_full(); } );
    m_work_queue.push( task );
    --m_push_waiting;
    m_work_on_queue_mutex.unlock();
    if ( m_pop_waiting > 0 ) m_queue_pop_cv.notify_one();
    if ( m_push_waiting > 0 && !m_work_queue.is_full() ) m_queue_push_cv.notify_one();
    // --------------------------
    m_tm.toc();
    m_push_ms += m_tm.elapsed_ms();
  }

  tp::Queue::TaskData *
  ThreadPool4::pop_task() {
    ++m_pop_waiting;
    m_work_on_queue_mutex.lock();
    m_queue_pop_cv.wait( m_work_on_queue_mutex, [&]()->bool { return !m_work_queue.empty(); } );
    TaskData * task = m_work_queue.pop();
    ++m_running_task; // must be incremented in the locked part
    --m_pop_waiting;
    m_work_on_queue_mutex.unlock();
    if ( m_push_waiting > 0 ) m_queue_push_cv.notify_one();
    if ( m_pop_waiting  > 0 && !m_work_queue.empty() ) m_queue_pop_cv.notify_one();
    return task;
  }

  void
  ThreadPool4::worker_thread(
    real_type & pop_ms,
    real_type & job_ms,
    unsigned  & n_job
  ) {
    TicToc tm;
    ++m_running_thread;
    while ( !m_done ) {
      // ---------------------------- POP
      tm.tic();
      TaskData * task = pop_task();
      tm.toc();
      pop_ms += tm.elapsed_ms();
      // ---------------------------- RUN
      tm.tic();
      (*task)(); // run and delete task;
      tm.toc();
      job_ms += tm.elapsed_ms();
      // ---------------------------- UPDATE
      --m_running_task; ++n_job;
    }
    --m_running_thread;
  }

  void
  ThreadPool4::create_workers( unsigned thread_count ) {
    m_worker_threads.clear();
    m_worker_threads.reserve(thread_count);
    m_job_ms.resize( std::size_t(thread_count) );
    m_pop_ms.resize( std::size_t(thread_count) );
    m_n_job.resize( std::size_t(thread_count) );
    std::fill( m_job_ms.begin(), m_job_ms.end(), 0 );
    std::fill( m_pop_ms.begin(), m_pop_ms.end(), 0 );
    std::fill( m_n_job.begin(), m_n_job.end(), 0 );
    m_push_ms = 0;
    m_done    = false;
    try {
      for ( unsigned i=0; i<thread_count; ++i )
        m_worker_threads.push_back(
          std::thread(
            &ThreadPool4::worker_thread, this,
            std::ref(m_pop_ms[i]),
            std::ref(m_job_ms[i]),
            std::ref(m_n_job[i])
          )
        );
    } catch(...) {
      m_done = true;
      throw;
    }
  }

  void
  ThreadPool4::resize( unsigned thread_count, unsigned queue_capacity ) {
    join();
    if ( queue_capacity == 0 ) queue_capacity = 4 * (thread_count+1);
    m_work_queue.resize( queue_capacity );
    create_workers( thread_count );
  }

  void
  ThreadPool4::wait()
  { while ( !m_work_queue.empty() || m_running_task > 0 ) nano_sleep(); }

  void
  ThreadPool4::join() {
    this->wait(); // finish all the running task
    m_done = true;
    unsigned i = m_running_thread;
    while ( i-- > 0 ) push_task( new TaskData([](){}) );
    while ( m_running_thread > 0 ) nano_sleep();
    m_work_queue.clear(); // remove spurious (null task) remained
    for ( std::thread & w : m_worker_threads ) { if (w.joinable()) w.join(); }
    m_worker_threads.clear(); // destroy the workers threads vector
  }

  void
  ThreadPool4::info( ostream_type & s ) const {
    unsigned nw = unsigned(m_pop_ms.size());
    for ( unsigned i = 0; i < nw; ++i ) {
      unsigned njob = m_n_job[i];
      fmt::print( s,
        "Worker {:2}, #job = {:4}, [job {:10}, POP {:10}]\n",
        i, njob,
        fmt::format( "{:.3} mus", 1000*m_job_ms[i]/njob ),
        fmt::format( "{:.3} mus", 1000*m_pop_ms[i]/njob )
      );
    }
    fmt::print( s, "PUSH {:10.6} ms\n\n", m_push_ms );
  }

}

///
/// eof: ThreadPool4.cc
///
