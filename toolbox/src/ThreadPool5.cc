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
/// file: ThreadPool5.cc
///

#include "Utils.hh"

namespace Utils {
  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  ThreadPool5::ThreadPool5( unsigned nthread )
  : ThreadPoolBase()
  {
    resize_workers( nthread );
    //info( std::cout );
  }

  ThreadPool5::~ThreadPool5() {
    join();
    m_workers.clear();
    m_stack.clear();
  }

  void
  ThreadPool5::resize_workers( unsigned numThreads ) {
    m_stack.clear(); // empty stack
    m_stack.reserve( size_t(numThreads) );
    m_workers.resize( size_t(numThreads) );
    unsigned id = 0;
    for ( Worker & w : m_workers ) { w.setup( this, id ); ++id; }
    while ( id-- > 0 ) push_worker( id );
    setup();
  }

  void
  ThreadPool5::push_worker( unsigned id ) {
    std::unique_lock<std::mutex> lock(m_stack_mutex);
    m_stack.push_back(id);
    m_stack_cond.notify_one();
  }

  unsigned
  ThreadPool5::pop_worker() {
    std::unique_lock<std::mutex> lock(m_stack_mutex);
    m_stack_cond.wait( lock, [&]()->bool { return !m_stack.empty(); } );
    unsigned id = m_stack.back(); m_stack.pop_back();
    return id;
  }

  void
  ThreadPool5::info_stack( ostream_type & s ) const {
    fmt::print( s, "STACK[{}]: ", m_stack.size() );
    for ( unsigned const & id : m_stack )
      fmt::print( s, "{}, ", id );
    s << '\n';
  }

  void
  ThreadPool5::info( ostream_type & s ) const {
    for ( Worker const & w : m_workers ) w.info(s);
    fmt::print( s,
      "LAUNCH {} ms\n"
      "POP    {} ms\n",
      m_exec_ms, m_pop_ms
    );
    info_stack( s );
    fmt::print( s, "\n" );
  }

  /*\
   |  __        __         _
   |  \ \      / /__  _ __| | _____ _ __
   |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
   |    \ V  V / (_) | |  |   <  __/ |
   |     \_/\_/ \___/|_|  |_|\_\___|_|
  \*/

  void
  ThreadPool5::Worker::start() {
    if ( !m_active ) {
      m_active = true;
      m_running_thread = std::thread( &Worker::worker_loop, this );
    }
  }

  void
  ThreadPool5::Worker::stop() {
    if ( m_active ) {
      wait();               // if running task wait it terminate
      m_active = false;     // deactivate computation
      m_job = [](){};       // dummy task
      m_is_running.green(); // start computation (exiting loop)
      if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
      m_is_running.red();   // end of computation (for double stop);
    }
    //fmt::print( "worker_loop {} stopped\n", m_worker_id );
  }

  void
  ThreadPool5::Worker::worker_loop() {
    m_is_running.red(); // block computation
    while ( m_active ) {
      m_tm.tic();
      m_is_running.wait(); // wait signal to start computation
      m_tm.toc();
      m_wait_ms += m_tm.elapsed_ms();
      // ----------------------------------------
      if ( !m_active ) break; // if finished exit
      m_tm.tic();
      m_job();
      m_tm.toc();
      m_job_ms += m_tm.elapsed_ms();
      // ----------------------------------------
      m_tm.tic();
      m_is_running.red();     // block computation
      ++m_job_done_counter;
      m_tp->push_worker( m_worker_id ); // worker ready for a new computation
      m_tm.toc();
      m_sync_ms += m_tm.elapsed_ms();
      std::this_thread::yield();
    }
  }

  void
  ThreadPool5::Worker::info( ostream_type & s ) const {
    fmt::print( s,
      "Worker {:2}, #job = {:4}, [job {:10}, sync {:10}, wait {:10}]\n",
      m_worker_id, m_job_done_counter,
      fmt::format( "{:.3} mus", 1000*elapsed_job_ms()/m_job_done_counter),
      fmt::format( "{:.3} mus", 1000*elapsed_sync_ms()/m_job_done_counter),
      fmt::format( "{:.3} mus", 1000*elapsed_wait_ms()/m_job_done_counter)
    );
  }

}

///
/// eof: ThreadPool5.cc
///
