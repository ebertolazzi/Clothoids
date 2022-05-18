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
/// file: ThreadPool1.cc
///

#include "Utils.hh"

namespace Utils {

  /*\
   |  __        __         _
   |  \ \      / /__  _ __| | _____ _ __
   |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
   |    \ V  V / (_) | |  |   <  __/ |
   |     \_/\_/ \___/|_|  |_|\_\___|_|
  \*/

  void
  ThreadPool1::Worker::worker_loop() {
    TicToc tm;
    m_running.wait(); // wait to start the first job
    while ( m_active ) {
      // ----------------------------
      tm.tic();
      m_job(); ++m_n_job; // run a job
      tm.toc();
      m_job_ms += tm.elapsed_ms();
      // ----------------------------
      tm.tic();
      m_running.red(); // job done
      m_running.wait(); // wait to start a new job
      tm.toc();
      m_wait_ms += tm.elapsed_ms();
    }
  }

  void
  ThreadPool1::Worker::start() {
    if ( !m_active ) {
      m_active = true;
      m_running.red();
      m_running_thread = std::thread( &Worker::worker_loop, this );
    }
  }

  void
  ThreadPool1::Worker::stop() {
    if ( m_active ) {
      m_active = false;     // deactivate computation
      m_running.wait_red(); // se gia occupato in task aspetta
      m_job = [](){};       // dummy task
      m_running.green();    // start computation
      if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
    }
  }

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  ThreadPool1::ThreadPool1( unsigned nthread )
  : ThreadPoolBase()
  {
    m_workers.resize( size_t( nthread ) );
    setup();
  }

  ThreadPool1::~ThreadPool1()
  { join(); m_workers.clear(); }

  void
  ThreadPool1::wait() {
    m_thread_to_send = 0;
    for ( auto && w : m_workers ) w.wait();
  }

  void
  ThreadPool1::resize( unsigned numThreads ) {
    wait();
    stop();
    m_workers.resize( size_t(numThreads) );
    setup();
  }

  void
  ThreadPool1::start()
  { m_thread_to_send = 0; for ( auto && w : m_workers ) w.start(); }

  void
  ThreadPool1::stop()
  { m_thread_to_send = 0; for ( auto && w : m_workers ) w.stop(); }

  void
  ThreadPool1::info( ostream_type & s ) const {
    unsigned i = 0;
    for ( auto const & w : m_workers )
      fmt::print( s,
        "Worker {:2}, #job = {:4}, [job {:10}, WAIT {:10}]"
        " PUSH = {:10}\n",
        i++, w.n_job(),
        fmt::format( "{:.3} mus", 1000*w.job_ms()/w.n_job() ),
        fmt::format( "{:.3} mus", 1000*w.wait_ms()/w.n_job() ),
        fmt::format( "{:.3} mus", 1000*w.push_ms()/w.n_job() )
      );
  }

}

///
/// eof: ThreadPool1.cc
///
