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
// file: ThreadPool1.hxx
//

namespace Utils
{

  /*!
   * \addtogroup THREAD
   * @{
   */

  class ThreadPool1 : public ThreadPoolBase
  {
    using FUN = ThreadPoolBase::FUN;

    class Worker
    {
      std::mutex              m_task_mutex;
      std::condition_variable m_task_cv;

      bool m_active{ true };     //!< Indicates if the worker is active
      bool m_task_done{ true };  //!< Indicates if the worker has completed its task
      FUN  m_task;               //!< Function to be executed by the worker

      // Startup synchronization
      std::promise<void>       m_startup_promise;
      std::shared_future<void> m_startup_future;

      std::thread m_running_thread;  //!< Thread executing the worker's tasks

      //! The main loop for the worker thread
      void worker_loop()
      {
        // Wait until constructor signals we're ready to proceed
        m_startup_future.wait();
        while ( true )
        {
          FUN task_to_run;
          {  // scope lock
            std::unique_lock<std::mutex> lock( m_task_mutex );
            // Wait for a new task or for deactivation
            m_task_cv.wait( lock, [&] { return ( !m_task_done || !m_active ); } );
            if ( !m_active ) break;
            // Move the task out to local, mark as not done
            task_to_run = std::move( m_task );
            m_task_done = false;
          }
          // Run the task outside the lock
          task_to_run();

          // After completion, mark as done and notify
          {
            std::unique_lock<std::mutex> lock( m_task_mutex );
            m_task_done = true;
            m_task_cv.notify_one();
          }
        }
      }

    public:
      // A Worker owns a std::thread started with &Worker::worker_loop and this.
      // Therefore its address must remain stable for the whole thread lifetime.
      // Moving a live Worker would leave the running thread referring to the old
      // object address.
      Worker( Worker const & )             = delete;
      Worker & operator=( Worker const & ) = delete;
      Worker( Worker && )                  = delete;
      Worker & operator=( Worker && )      = delete;

      Worker() : m_active( true ), m_task_done( true ), m_startup_future( m_startup_promise.get_future().share() )
      {
        // Now that everything is initialized, start the thread
        m_running_thread = std::thread( &Worker::worker_loop, this );
        // Signal that initialization is done and the thread can start the main
        // loop
        m_startup_promise.set_value();
      }

      ~Worker()
      {
        auto dummy_task = []() -> void {};
        {
          std::unique_lock<std::mutex> lock( m_task_mutex );
          m_active    = false;  // deactivate computation
          m_task_done = false;  // Wake worker if waiting
          m_task      = dummy_task;
          m_task_cv.notify_one();
        }
        if ( m_running_thread.joinable() ) m_running_thread.join();
      }


      void wait()
      {
        std::unique_lock<std::mutex> lock( m_task_mutex );
        while ( !m_task_done ) m_task_cv.wait( lock );
      }

      void exec( FUN && fun )
      {
        std::unique_lock<std::mutex> lock( m_task_mutex );
        while ( !m_task_done ) m_task_cv.wait( lock );
        m_task_done = false;
        m_task      = std::move( fun );
        m_task_cv.notify_one();
      }
    };

    std::size_t         m_thread_to_send{ 0 };
    std::vector<Worker> m_workers;

  public:
    explicit ThreadPool1(
      unsigned nthread = std::max( unsigned( 1 ), unsigned( std::thread::hardware_concurrency() - 1 ) ) )
      : ThreadPoolBase(), m_workers( size_t( nthread ) )
    {
    }

    virtual ~ThreadPool1() override { m_workers.clear(); }

    void exec( FUN && fun ) override
    {
      m_workers[m_thread_to_send].exec( std::move( fun ) );
      if ( ++m_thread_to_send >= m_workers.size() ) m_thread_to_send = 0;
    }

    void wait() override
    {
      for ( auto && w : m_workers ) w.wait();
    }

    unsigned thread_count() const override { return unsigned( m_workers.size() ); }

    static char const * Name() { return "ThreadPool1"; }

    char const * name() const override { return Name(); }
  };

  /*! @} */

}  // namespace Utils

//
// eof: ThreadPool1.hxx
//
