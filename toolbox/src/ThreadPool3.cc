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
/// file: ThreadPool3.cc
///

#include "Utils.hh"

namespace Utils {

  ThreadPool3::ThreadPool3( unsigned thread_count, unsigned queue_capacity )
  : m_done(false)
  , m_running_task(0)
  , m_running_thread(0)
  , m_work_queue( queue_capacity == 0 ? std::max( 10 * (thread_count+1), unsigned(4096) ) : queue_capacity )
  , m_push_waiting(0)
  , m_pop_waiting(0)
  {
    create_workers( thread_count );
  }

  //
  // https://stackoverflow.com/questions/48936591/is-the-performance-of-notify-one-really-this-bad
  //

  void
  ThreadPool3::push_task( TaskData * task ) {
    m_queue_push_mutex.lock();
    ++m_push_waiting;
    m_queue_push_cv.wait( m_queue_push_mutex, [&]()->bool { return !m_work_queue.is_full(); } );
    //-----------------
    m_queue_spin.lock();
    m_work_queue.push( task );
    --m_push_waiting;
    m_queue_spin.unlock();
    //-----------------
    m_queue_push_mutex.unlock();
    // push done
    if ( m_pop_waiting > 0 ) m_queue_pop_cv.notify_one();
    if ( m_push_waiting > 0 ) {
      m_queue_spin.lock();
      if ( !m_work_queue.is_full() ) m_queue_push_cv.notify_one();
      m_queue_spin.unlock();
    }
  }

  tp::Queue::TaskData *
  ThreadPool3::pop_task() {
    m_queue_pop_mutex.lock();
    ++m_pop_waiting;
    m_queue_pop_cv.wait( m_queue_pop_mutex, [&]()->bool { return !m_work_queue.empty(); } );
    //-----------------
    m_queue_spin.lock();
    TaskData * task = m_work_queue.pop();
    --m_pop_waiting;
    m_queue_spin.unlock();
    ++m_running_task; // must be incremented in the locked part
    //-----------------
    m_queue_pop_mutex.unlock();
    if ( m_push_waiting > 0 ) m_queue_push_cv.notify_one();
    if ( m_pop_waiting  > 0 ) {
      m_queue_spin.lock();
      if ( !m_work_queue.empty() ) m_queue_pop_cv.notify_one();
      m_queue_spin.unlock();
    }
    return task;
  }

  void
  ThreadPool3::worker_thread(
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
  ThreadPool3::create_workers( unsigned thread_count ) {
    m_worker_threads.clear();
    m_worker_threads.reserve(thread_count);
    m_job_ms.resize( std::size_t(thread_count) );
    m_pop_ms.resize( std::size_t(thread_count) );
    m_n_job.resize( std::size_t(thread_count) );
    std::fill( m_job_ms.begin(), m_job_ms.end(), 0 );
    std::fill( m_pop_ms.begin(), m_pop_ms.end(), 0 );
    std::fill( m_n_job.begin(), m_n_job.end(), 0 );
    m_push_ms      = 0;
    m_done         = false;
    m_push_waiting = 0;
    m_pop_waiting  = 0;
    try {
      for ( unsigned i=0; i<thread_count; ++i )
        m_worker_threads.emplace_back(
          std::thread(
            &ThreadPool3::worker_thread, this,
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
  ThreadPool3::join() {
    wait();
    m_done = true;
    { // send null task until all the workers stopped
      std::function<void()> null_job = [](){};
      for ( unsigned i = m_running_thread; i > 0; --i )
        push_task( new TaskData(null_job) );
      while ( m_running_thread > 0 ) std::this_thread::yield();
    }
    m_work_queue.clear();
    for ( std::thread & w : m_worker_threads ) { if (w.joinable()) w.join(); }
    m_worker_threads.clear();
  }

  void
  ThreadPool3::resize( unsigned thread_count, unsigned queue_capacity ) {
    join();
    if ( queue_capacity == 0 ) queue_capacity = 4 * (thread_count+1);
    m_work_queue.resize( queue_capacity );
    create_workers( thread_count );
  }

  void
  ThreadPool3::info( ostream_type & s ) const {
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
/// eof: ThreadPool3.cc
///
