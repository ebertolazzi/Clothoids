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

      using real_type = double;

      bool                  m_active;         //!< Indicates if the worker is active
      UTILS_SEMAPHORE       m_running;        //!< Semaphore for controlling task execution
      std::thread           m_running_thread; //!< Thread executing the worker's tasks
      std::function<void()> m_job;            //!< Function to be executed by the worker

      TicToc m_tm; //!< Timing object for measuring task execution time

      void worker_loop(); //!< The main loop for the worker thread

    public:

      //disable copy
      Worker( Worker const & )              = delete;
      //Worker( Worker && )                   = delete;
      Worker& operator = ( Worker const & ) = delete;
      Worker& operator = ( Worker && )      = delete;

      //!
      //! \brief Constructs a new Worker and starts the thread.
      //!
      Worker() : m_active(false) { start(); }

      //!
      //! \brief Destroys the Worker and stops the thread.
      //!
      ~Worker() { stop(); }

      //!
      //! \brief Moves the Worker from one instance to another.
      //!
      Worker( Worker && rhs ) noexcept
      : m_active(rhs.m_active)
      , m_running_thread(std::move(rhs.m_running_thread))
      , m_job(std::move(rhs.m_job))
      {}

      void start(); //!< Starts the worker thread
      void stop();  //!< Stops the worker thread

      //!
      //! \brief Waits for the current task to finish.
      //!
      void wait() { m_running.wait_red(); }

      //!
      //! \brief Executes a function in the worker thread.
      //!
      //! \param fun The function to execute.
      //!
      void
      exec( std::function<void()> & fun ) {
        m_running.wait_red(); // se gia occupato in task aspetta
        m_job = fun;          // cambia funzione da eseguire
        m_running.green();    // activate computation
      }

      std::thread::id     get_id()     const { return m_running_thread.get_id(); }
      std::thread const & get_thread() const { return m_running_thread; }
      std::thread &       get_thread()       { return m_running_thread; }
    };

    //! Index of the next thread to send a task to
    std::size_t m_thread_to_send{0};

    //! Vector of worker threads
    std::vector<Worker> m_workers;

    //! Initializes and starts all workers
    void setup() { for ( auto & w: m_workers ) w.start(); }

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
    );

    //!
    //! \brief Destroys the ThreadPool1 and stops all worker threads.
    //!
    virtual ~ThreadPool1();

    //!
    //! \brief Executes a task in the thread pool.
    //!
    //! \param fun The function to be executed.
    //!
    void
    exec( std::function<void()> && fun ) override {
      m_workers[m_thread_to_send].exec( fun );
      if ( ++m_thread_to_send >= m_workers.size() ) m_thread_to_send = 0;
    }

    void wait() override;             //!< Waits for all tasks to finish
    void join() override { stop(); }  //!< Stops and joins all threads

    //!
    //! \brief Returns the number of threads in the pool.
    //!
    //! \return The number of threads.
    //!
    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    void resize( unsigned numThreads ) override; //!< Resizes the thread pool
    char const * name() const override { return "ThreadPool1"; } //!< Returns the name of the thread pool

    // EXTRA
    void start(); //!< Starts all worker threads
    void stop();  //!< Stops all worker threads

    //!
    //! \brief Returns the ID of the specified worker thread.
    //!
    //! \param i Index of the worker.
    //! \return The thread ID.
    //!
    std::thread::id
    get_id( unsigned i ) const
    { return m_workers[size_t(i)].get_id(); }

    //!
    //! \brief Returns the thread object of the specified worker.
    //!
    //! \param i Index of the worker.
    //! \return The thread object.
    //!
    std::thread const &
    get_thread( unsigned i ) const
    { return m_workers[size_t(i)].get_thread(); }

    std::thread &
    get_thread( unsigned i )
    { return m_workers[size_t(i)].get_thread(); }

    // ALIAS
    void wait_all()  { this->wait();  } //!< Alias for wait()
    void start_all() { this->start(); } //!<  Alias for start()
    void stop_all()  { this->stop();  } //!<  Alias for stop()
    unsigned size() const { return this->thread_count(); } //!< Returns the number of threads
  };

  using ThreadPool = ThreadPool1;

  /*! @} */

}

//
// eof: ThreadPool1.hxx
//
