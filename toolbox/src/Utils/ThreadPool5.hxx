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

    using real_type = double;

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

      bool                  m_active           = false;   //!< Indicates if the worker is active.
      unsigned              m_job_done_counter = 0;       //!< Counter for completed jobs.
      unsigned              m_worker_id        = 0;       //!< Unique ID for the worker.
      ThreadPool5 *         m_tp               = nullptr; //!< Pointer to the parent thread pool.
      UTILS_SEMAPHORE       m_is_running;                 //!< Semaphore to manage task execution.
      std::thread           m_running_thread;             //!< The thread that runs the worker loop.
      std::function<void()> m_job;                        //!< Function to be executed by the worker.
      TicToc                m_tm;                         //!< Timing utility for measuring execution durations.
      real_type             m_job_ms  = 0;                //!< Time taken for job execution.
      real_type             m_sync_ms = 0;                //!< Time spent in synchronization.
      real_type             m_wait_ms = 0;                //!< Time spent waiting for tasks.

      //!
      //! \brief The main loop for the worker thread, continuously fetching and executing tasks.
      //!
      void worker_loop();

    public:

      //!
      //! \brief Constructs a new worker and starts its thread.
      //!
      explicit Worker() { start(); }

      //!
      //! \brief Destructor for the Worker class.
      //!
      //! Ensures the worker stops execution before destruction.
      //!
      ~Worker() { stop(); }

      // dummy copy constructor for resize
      Worker( Worker && ) {}

      //!
      //! \brief Sets up the worker with its thread pool and ID.
      //!
      //! \param tp Pointer to the parent ThreadPool5.
      //! \param id The unique ID of the worker.
      //!
      void
      setup( ThreadPool5 * tp, unsigned id ) {
        m_worker_id        = id;
        m_job_done_counter = 0;
        m_tp               = tp;
      }

      //!
      //! \brief Starts the worker thread.
      //!
      void start();

      //!
      //! \brief Waits for the currently running task to complete.
      //!
      void wait() { m_is_running.wait_red(); }

      //!
      //! \brief Stops the worker and waits for the current task to finish.
      //!
      void stop();

      //!
      //! \brief Executes a task.
      //!
      //! \param fun The function to be executed as a task.
      //!
      void
      exec( std::function<void()> & fun ) {
        m_is_running.wait_red();
        m_job = fun;          // cambia funzione da eseguire
        m_is_running.green(); // activate computation
      }

      //!
      //! \brief Gets the number of jobs completed by this worker.
      //!
      //! \return The count of completed jobs.
      //!
      unsigned job_done_counter() const { return m_job_done_counter; }

      //!
      //! \brief Gets the elapsed time spent executing jobs.
      //!
      //! \return The elapsed time in milliseconds.
      //!
      real_type elapsed_job_ms()  const { return m_job_ms; }

      //!
      //! \brief Gets the elapsed time spent in synchronization.
      //!
      //! \return The elapsed time in milliseconds.
      //!
      real_type elapsed_sync_ms() const { return m_sync_ms; }

      //!
      //! \brief Gets the elapsed time spent waiting for tasks.
      //!
      //! \return The elapsed time in milliseconds.
      //!
      real_type elapsed_wait_ms() const { return m_wait_ms; }

      //!
      //! \brief Outputs information about the worker's performance.
      //!
      //! \param s The output stream to which information will be written.
      //!
      void info( ostream_type & s ) const;
    };

    // =========================================================================
    // =========================================================================
    // =========================================================================

    // Vector of workers managed by the thread pool.
    std::vector<Worker>     m_workers;
    // Stack of available worker IDs.
    std::vector<unsigned>   m_stack;
    std::mutex              m_stack_mutex; //!< Mutex for accessing the worker stack.
    std::condition_variable m_stack_cond;  //!< Condition variable for worker availability.
    TicToc                  m_tm;          //!< Timing utility for measuring performance.
    real_type               m_exec_ms = 0; //!< Total execution time of tasks.
    real_type               m_pop_ms  = 0; //!< Total time spent popping tasks.

    //!
    //! \brief Initializes and starts all workers.
    //!
    void setup() { for ( auto & w: m_workers ) w.start(); }

    //!
    //! \brief Resizes the number of workers in the thread pool.
    //!
    //! \param numThreads The new number of worker threads.
    //!
    void resize_workers( unsigned numThreads );

    //!
    //! \brief Pushes a worker ID back onto the stack of available workers.
    //!
    //! \param id The ID of the worker to be pushed.
    //!
    void push_worker( unsigned id );

    //!
    //! \brief Pops a worker ID from the stack of available workers.
    //!
    //! \return The ID of the popped worker.
    //!
    unsigned pop_worker();

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
    );

    //!
    //! \brief Destructor for the ThreadPool5 class.
    //!
    //! Ensures all workers are stopped and joined before destruction.
    //!
    virtual ~ThreadPool5();

    //!
    //! \brief Executes a task and assigns it to an available worker.
    //!
    //! \param fun The function to be executed as a task.
    //!
    void
    exec( std::function<void()> && fun ) override {
      // cerca prima thread libera

      m_tm.tic();
      Worker & w = m_workers[pop_worker()];
      m_tm.toc();
      m_pop_ms += m_tm.elapsed_ms();

      m_tm.tic();
      w.exec( fun );
      m_tm.toc();
      m_exec_ms += m_tm.elapsed_ms();
    }

    //!
    //! \brief Waits for all tasks to be completed.
    //!
    void
    wait() override
    { for ( auto & w : m_workers ) w.wait(); }

    //!
    //! \brief Joins all worker threads and ensures they are stopped.
    //!
    void
    join() override
    { stop(); }

    //!
    //! \brief Gets the current number of threads in the pool.
    //!
    //! \return The number of threads in the pool.
    //!
    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    //!
    //! \brief Resizes the thread pool to a new number of workers.
    //!
    //! \param numThreads The new number of threads to maintain in the pool.
    //!
    void
    resize( unsigned numThreads ) override
    { wait(); stop(); resize_workers( numThreads ); }

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    char const * name() const override { return "ThreadPool5"; }

    //!
    //! \brief Starts all workers in the pool.
    //!
    void start() { for ( auto && w : m_workers ) w.start(); }

    //!
    //! \brief Stops all workers in the pool.
    //!
    void stop() { for ( auto && w : m_workers ) w.stop(); }

    //!
    //! \brief Outputs information about the stack of available workers.
    //!
    //! \param s The output stream to which information will be written.
    //!
    void info_stack( ostream_type & s ) const;

    //!
    //! \brief Outputs information about the thread pool's performance.
    //!
    //! \param s The output stream to which information will be written.
    //!
    void info( ostream_type & s ) const override;
  };

  /*! @} */

}

//
// eof: ThreadPool5.hxx
//
