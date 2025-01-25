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
// file: ThreadPool5.hxx
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
  //! \brief A thread pool for concurrent task execution with worker management.
  //!
  //! The `ThreadPool5` class provides an implementation of a thread pool that
  //! allows for efficient concurrent execution of tasks. This class manages a
  //! collection of worker threads that can be dynamically resized, and utilizes
  //! a semaphore-based mechanism for worker synchronization.
  //!
  //! This class extends `ThreadPoolBase` and supports task execution, joining,
  //! and resizing of worker threads.
  //!
  class ThreadPool5 : public ThreadPoolBase {

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
    //! \brief Represents an individual worker thread in the thread pool.
    //!
    //! Each `Worker` is responsible for executing tasks and reporting their
    //! completion status. Workers are managed by the `ThreadPool5` class and
    //! can wait for tasks to be assigned to them.
    //!
    class Worker {
    
    
      FUN             m_push_worker;
      bool            m_active{true};    //!< Indicates if the worker is active.
      UTILS_SEMAPHORE m_is_running;      //!< Semaphore to manage task execution.
      std::thread     m_running_thread;  //!< The thread that runs the worker loop.
      FUN             m_job;             //!< Function to be executed by the worker.

      //!
      //! \brief The main loop for the worker thread, continuously fetching and executing tasks.
      //!
      void
      worker_loop() {
        m_is_running.red(); // block computation
        while ( m_active ) {
          m_is_running.wait(); // wait signal to start computation
          // ----------------------------------------
          m_job();
          // ----------------------------------------
          m_is_running.red(); // block computation
          m_push_worker();    // worker ready for a new computation
        }
      }

    public:

      Worker() { m_is_running.red(); }

      //!
      //! \brief Constructs a new worker and starts its thread.
      //!
      //! \param tp Pointer to the parent ThreadPool5.
      //! \param id The unique ID of the worker.
      //!
      explicit
      Worker( ThreadPool5 * tp, unsigned id ) {
        m_is_running.red();
        m_push_worker    = [tp,id]()->void { tp->push_worker( id ); };
        m_running_thread = std::thread( &Worker::worker_loop, this );
      }

      //!
      //! \brief Destructor for the Worker class.
      //!
      //! Ensures the worker stops execution before destruction.
      //!
      ~Worker() {
        wait();               // if running task wait it terminate
        m_active = false;     // deactivate computation
        m_job    = [](){};    // dummy task
        m_is_running.green(); // start computation (exiting loop)
        if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
      }

      //!
      //! \brief Waits for the currently running task to complete.
      //!
      void wait() { m_is_running.wait_red(); }

      //!
      //! \brief Executes a task.
      //!
      //! \param fun The function to be executed as a task.
      //!
      void
      exec( FUN && fun ) {
        m_is_running.wait_red();
        m_job = std::move(fun); // cambia funzione da eseguire
        m_is_running.green();   // activate computation
      }

    };

    // =========================================================================
    // =========================================================================
    // =========================================================================

    // Vector of workers managed by the thread pool.
    std::vector<Worker>     m_workers;
    std::list<unsigned>     m_queue;
    std::mutex              m_queue_mutex; //!< Mutex for accessing the worker stack.
    std::condition_variable m_queue_cond;  //!< Condition variable for worker availability.

    //!
    //! \brief Pushes a worker ID back onto the stack of available workers.
    //!
    //! \param id The ID of the worker to be pushed.
    //!
    void
    push_worker( unsigned id ) {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      m_queue.emplace_back(id);
      m_queue_cond.notify_one();
    }

    //!
    //! \brief Pops a worker ID from the stack of available workers.
    //!
    //! \return The ID of the popped worker.
    //!
    unsigned
    pop_worker() {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      while ( m_queue.empty() ) m_queue_cond.wait( lock );
      unsigned id{ m_queue.back() }; m_queue.pop_back();
      return id;
    }

  public:

    //!
    //! \brief Constructs a new ThreadPool5 instance with a specified number of threads.
    //!
    //! \param nthread The number of threads to create in the pool. Defaults to the maximum hardware threads available.
    //!
    ThreadPool5(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    )
    : ThreadPoolBase()
    , m_workers( size_t(nthread) )
    {
      m_queue.clear();
      unsigned id{0};
      for ( Worker & w : m_workers ) { new (&w) Worker( this, id ); ++id; }
      while ( id-- > 0 ) push_worker( id );
    }

    //!
    //! \brief Destructor for the ThreadPool5 class.
    //!
    //! Ensures all workers are stopped and joined before destruction.
    //!
    virtual
    ~ThreadPool5() {
      wait();
      m_workers.clear();
    }

    //!
    //! \brief Executes a task and assigns it to an available worker.
    //!
    //! \param fun The function to be executed as a task.
    //!
    void
    exec( FUN && fun ) override {
      // cerca prima thread libera
      m_workers[pop_worker()].exec( std::move(fun) );
    }

    //!
    //! \brief Waits for all tasks to be completed.
    //!
    void
    wait() override
    { for ( auto & w : m_workers ) w.wait(); }

    //!
    //! \brief Gets the current number of threads in the pool.
    //!
    //! \return The number of threads in the pool.
    //!
    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    static char const * Name() { return "ThreadPool5"; }

    char const * name() const override { return Name(); }

  };

  /*! @} */

}

//
// eof: ThreadPool5.hxx
//
