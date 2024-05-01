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
/// file: ThreadPool2.cc
///

#include "Utils.hh"

namespace Utils {

  void
  HQueue::FixedCapacityQueue::reserve( unsigned capacity ) {
    assert(empty()); // Copying / moving of Fun not supported.
    if ( capacity != m_capacity ) {
      m_size     = capacity+1;
      m_capacity = capacity;
      m_fun_vec.resize( m_size );
    }
  }

  HQueue::FixedCapacityQueue::~FixedCapacityQueue() {
    while (!empty()) pop();
  }

  void
  HQueue::help( int return_if_idle ) {

    unsigned min_queue_size = return_if_idle < 0 ? 0 : return_if_idle;

    // Increment total worker count, decrement again on scope exit
    { std::lock_guard<std::mutex> lock(m_push_mutex); ++m_total_workers; }

    // execute at exit
    auto x1 = at_scope_exit([this](){
      std::lock_guard<std::mutex> lock(this->m_push_mutex);
      if (--this->m_total_workers == this->m_idle_workers)
        this->m_waiters_cond.notify_all();
    });

    FixedCapacityQueue function_queue(1);

    while ( true ) {
      std::unique_lock<std::mutex> lock(m_pop_mutex);
      unsigned queue_size;

      // Try to get the next task(s)
      while ((queue_size = m_queue.size()) <= min_queue_size) {
        if ( int(queue_size) <= return_if_idle) return;
        if ( queue_size > 0 ) break;
        // The queue is empty, wait for more tasks to be put()
        lock.unlock();
        {
          std::unique_lock<std::mutex> lock2(m_push_mutex);
          while (m_queue.empty() && !m_shutting_down) {
            if ( ++m_idle_workers == m_total_workers ) m_waiters_cond.notify_all();
            m_waiting_workers_cond.wait(lock2); // Wait for task to be queued
            m_wakeup_is_pending = false;
            --m_idle_workers;
          }
        }
        if (m_shutting_down) return;
        lock.lock();
      }

      // There is at least one task in the queue and the back is locked.

      unsigned stride = (m_maxpart == 0) ? 1 : queue_size / m_maxpart;
      if (stride <= 0) stride = 1;
      if (stride > function_queue.capacity()) function_queue.reserve(2 * stride);
      while (stride--) function_queue.push(m_queue.pop());
      lock.unlock();

      if ( m_idle_workers && !m_wakeup_is_pending && queue_size )
        m_waiting_workers_cond.notify_one();

      while (!function_queue.empty()) function_queue.pop()();
    }
  }

  void
  HQueue::shutdown() {
    std::unique_lock<std::mutex> push_lock(m_push_mutex);
    std::unique_lock<std::mutex> pop_lock(m_pop_mutex);
    m_shutting_down = true;
    while (!m_queue.empty()) m_queue.pop();
    m_waiting_workers_cond.notify_all();
    m_waiters_cond.notify_all();
  }

  void
  HQueue:: wait() {
    if ( std::uncaught_exception() ) shutdown();
    std::exception_ptr e;
    std::unique_lock<std::mutex> lock(m_push_mutex);
    while ( !m_queue.empty() || m_idle_workers != m_total_workers ) {
      while ( !m_queue.empty() ) {
        lock.unlock();
        try {
          try_help(0);
        } catch (...) {
          if (e == nullptr) e = std::current_exception();
        }
        lock.lock();
      }
      while ( m_idle_workers != m_total_workers ) m_waiters_cond.wait(lock);
    }
    if (e != nullptr && !std::uncaught_exception())
      std::rethrow_exception(std::move(e));
  }

  ThreadPool2::ThreadPool2( unsigned thread_count )
  : m_queue( 50 * (thread_count+1), 3 * (thread_count+1) )
  , m_worker_threads( thread_count )
  {
    for ( std::thread & w : m_worker_threads )
      w = std::thread( std::bind( &HQueue::work, &m_queue, false ) );
  }

  void
  ThreadPool2::wait() {
    m_queue.work( true ); // Help out instead of sitting around idly.
    m_queue.wait();
  }

  void
  ThreadPool2::join() {
    m_queue.work( true ); // Help out instead of sitting around idly.
    m_queue.wait();
    m_queue.shutdown();
    m_queue.work( false ); // Instead of hanging around, help the workers!
    for ( std::thread & w : m_worker_threads )
      { if (w.joinable()) w.join(); }
  }

  void
  ThreadPool2::resize(
    unsigned thread_count,
    unsigned queue_size,
    unsigned maxpart
  ) {
    join();
    if ( queue_size == 0 ) queue_size = 50 * (thread_count+1);
    if ( maxpart    == 0 ) maxpart    = 3 * (thread_count+1);

    new (&m_queue) HQueue( queue_size, maxpart );

    m_worker_threads.clear();
    m_worker_threads.resize( thread_count );
    for ( std::thread & w : m_worker_threads )
      w = std::thread( std::bind( &HQueue::work, &m_queue, false ) );
  }

  ThreadPool2::~ThreadPool2() {
    // Abort processing if destructor runs during exception handling.
    if ( std::uncaught_exception() ) m_queue.shutdown();
    join(); // Running threads would continue to access the destructed pool.
  }

}

///
/// eof: ThreadPool2.cc
///
