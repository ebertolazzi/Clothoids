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
// file: ThreadPool1.hxx
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
  //! The `ThreadPool1` class is derived from `ThreadPoolBase` and provides
  //! functionality to execute tasks concurrently using a pool of worker threads.
  //!
  class ThreadPool1 : public ThreadPoolBase {

    using FUN = ThreadPoolBase::FUN;
    //using PFUN = ThreadPoolBase::PFUN;

    /*\
     |  __        __         _
     |  \ \      / /__  _ __| | _____ _ __
     |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
     |    \ V  V / (_) | |  |   <  __/ |
     |     \_/\_/ \___/|_|  |_|\_\___|_|
    \*/
    //!
    //! \brief Represents a worker thread in the thread pool.
    //!
    //! The `Worker` class manages an individual thread that can execute tasks.
    //! It provides methods to start, stop, and execute functions on the worker thread.
    //!
    class Worker {
      
      std::mutex              m_task_mutex;
      std::condition_variable m_task_cv;

      bool        m_active{true};     //!< Indicates if the worker is active
      std::thread m_running_thread;   //!< Thread executing the worker's tasks
      bool        m_task_done{true};  //!< Indicates if the worker is active
      FUN         m_task;             //!< Function to be executed by the worker

      //! The main loop for the worker thread
      void
      worker_loop() {
        std::unique_lock<std::mutex> lock(m_task_mutex);
        while ( m_active ) {
          while ( m_task_done ) m_task_cv.wait( lock );
          m_task();
          m_task_done = true;
          m_task_cv.notify_one();
        }
      }

    public:

      //disable copy
      Worker( Worker const & )              = delete;
      //Worker( Worker && )                   = delete;
      Worker& operator = ( Worker const & ) = delete;
      Worker& operator = ( Worker && )      = delete;

      //!
      //! \brief Constructs a new Worker and starts the thread.
      //!
      Worker() : m_running_thread( std::thread( &Worker::worker_loop, this ) ) {}

      //!
      //! \brief Destroys the Worker and stops the thread.
      //!
      ~Worker() {
        auto dummy_task = []()->void{};
        m_active = false; // deactivate computation
        exec( dummy_task );
        if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
      }
      
      //!
      //! \brief Moves the Worker from one instance to another.
      //!
      Worker( Worker && rhs ) noexcept
      : m_active(rhs.m_active)
      , m_running_thread(std::move(rhs.m_running_thread))
      , m_task(std::move(rhs.m_task))
      {}

      //!
      //! \brief Waits for the current task to finish.
      //!
      void
      wait() {
        std::unique_lock<std::mutex> lock(m_task_mutex);
        while ( !m_task_done ) m_task_cv.wait( lock );
      }

      //!
      //! \brief Executes a function in the worker thread.
      //!
      //! \param fun The function to execute.
      //!
      void
      exec( FUN && fun ) {
        // aspetto che pointer si liberi se occupato
        std::unique_lock<std::mutex> lock(m_task_mutex);
        while ( !m_task_done ) m_task_cv.wait( lock );
        m_task_done = false;
        m_task = std::move(fun); // cambia funzione da eseguire
        m_task_cv.notify_one();
      }
    };

    //! Index of the next thread to send a task to
    std::size_t m_thread_to_send{0};

    //! Vector of worker threads
    std::vector<Worker> m_workers;

  public:

    //!
    //! \brief Constructs a new ThreadPool1 with a specified number of threads.
    //!
    //! \param nthread Number of threads to create (default: hardware concurrency - 1).
    //!
    explicit
    ThreadPool1(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    ) : ThreadPoolBase(), m_workers( size_t( nthread ) ) { }

    //!
    //! \brief Destroys the ThreadPool1 and stops all worker threads.
    //!
    virtual
    ~ThreadPool1() { m_workers.clear(); }

    //!
    //! \brief Executes a task in the thread pool.
    //!
    //! \param fun The function to be executed.
    //!
    void
    exec( FUN && fun ) override {
      m_workers[m_thread_to_send].exec( std::move(fun) );
      if ( ++m_thread_to_send >= m_workers.size() ) m_thread_to_send = 0;
    }

    //! Waits for all tasks to finish
    void wait() override { for ( auto && w : m_workers ) w.wait(); }

    //!
    //! \brief Returns the number of threads in the pool.
    //!
    //! \return The number of threads.
    //!
    unsigned thread_count() const override { return unsigned(m_workers.size()); }

    static char const * Name() { return "ThreadPool1"; } //!< Returns the name of the thread pool

    char const * name() const override { return Name(); }

  };

  /*! @} */

}

//
// eof: ThreadPool1.hxx
//
