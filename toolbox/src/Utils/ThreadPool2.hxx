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

//
// file: ThreadPool2.hxx
//

namespace Utils {

  /*!
   * \addtogroup THREAD
   * @{
   */

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  //!
  //! \brief Manages a pool of worker threads for concurrent task execution.
  //!
  //! The `ThreadPool2` class is derived from `ThreadPoolBase` and provides
  //! functionality to execute tasks concurrently using a pool of worker threads.
  //! Based on tghe code from https://github.com/nbsdx/ThreadPool/tree/master
  //!
  class ThreadPool2 : public ThreadPoolBase {

    using TYPE = std::function<void(void)>;

    std::atomic_int          m_jobs_left{0};
    std::atomic_bool         m_bailout{false};
    std::atomic_bool         m_finished{false};

    std::condition_variable  m_job_available_var;
    std::condition_variable  m_wait_var;
    std::mutex               m_wait_mutex;
    std::mutex               m_queue_mutex;

    std::list<TYPE>          m_queue;
    std::vector<std::thread> m_threads;

    /*!
     *  Take the next job in the queue and run it.
     *  Notify the main thread that a job has completed.
     */
    void
    Task() {
      TYPE job;
      while( !m_bailout ) {
        {
          std::unique_lock<std::mutex> job_lock( m_queue_mutex );
          while ( m_queue.empty() && !m_bailout ) m_job_available_var.wait( job_lock );
          if( !m_bailout ) { job = m_queue.front(); m_queue.pop_front(); }
        }
        if ( !m_bailout ) { job(); --m_jobs_left; }
        {
          std::lock_guard<std::mutex> wait_lock( m_wait_mutex );
          m_wait_var.notify_one();
        }
      }
    }


  public:

    //!
    //! \brief Constructs a new ThreadPool1 with a specified number of threads.
    //!
    //! \param nthread Number of threads to create (default: hardware concurrency - 1).
    //!
    explicit
    ThreadPool2( unsigned nthread = std::max( unsigned(1), unsigned(std::thread::hardware_concurrency()-1) )) {
      m_threads.clear();
      m_threads.reserve(nthread);
      for( unsigned i{0}; i < nthread; ++i )
        m_threads.emplace_back( [this]{ this->Task(); } );
    }

    //!
    //! \brief Destroys the ThreadPool1 and stops all worker threads.
    //!
    virtual ~ThreadPool2() { join(); }

    //!
    //! \brief Executes a task in the thread pool.
    //!
    //! \param job The function to be executed.
    //!
    void
    exec( TYPE && job ) override {
      std::lock_guard<std::mutex> guard( m_queue_mutex );
      ++m_jobs_left;
      m_queue.emplace_back( job );
      m_job_available_var.notify_one();
    }

    void
    wait() override {
      std::unique_lock<std::mutex> lk( m_wait_mutex );
      while ( m_jobs_left > 0 ) m_wait_var.wait( lk );
    }

    void
    join() override {
      if( m_finished ) return;
      wait();
      // note that we're done, and wake up any thread that's
      // waiting for a new job
      m_bailout = true;
      m_job_available_var.notify_all();
      for( auto & t : m_threads )
        if( t.joinable() )
          t.join();
      m_finished = true;
    }

    //!
    //! \brief Returns the number of threads in the pool.
    //!
    //! \return The number of threads.
    //!
    unsigned thread_count() const override { return unsigned(m_threads.size()); }

    void
    resize( unsigned nthread ) override {
      join();
      m_jobs_left = 0;
      m_bailout   = false;
      m_finished  = false;
      m_threads.clear();
      m_threads.reserve(nthread);
      for( unsigned i{0}; i < nthread; ++i )
        m_threads.emplace_back( [this]{ this->Task(); } );
    };

    char const * name() const override { return "ThreadPool2"; } //!< Returns the name of the thread pool

    //!
    //! \brief Returns the ID of the specified worker thread.
    //!
    //! \param i Index of the worker.
    //! \return The thread ID.
    //!
    std::thread::id get_id( unsigned i ) const { return m_threads[size_t(i)].get_id(); }

    //!
    //! \brief Returns the thread object of the specified worker.
    //!
    //! \param i Index of the worker.
    //! \return The thread object.
    //!
    std::thread const & get_thread( unsigned i ) const { return m_threads[size_t(i)]; }
    std::thread       & get_thread( unsigned i )       { return m_threads[size_t(i)]; }
  };

  /*! @} */

}

//
// eof: ThreadPool2.hxx
//
