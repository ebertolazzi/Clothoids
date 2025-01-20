//
// file: ThreadPool4.hxx
//

namespace Utils {

  /*!
   * \addtogroup THREAD
   * @{
   */

  //!
  //! \brief A thread pool for concurrent task execution with enhanced synchronization.
  //!
  //! The `ThreadPool4` class provides a flexible and efficient thread pool implementation
  //! that allows for concurrent execution of tasks. It includes mechanisms for dynamic
  //! resizing of the thread pool and effective management of task synchronization using
  //! atomic operations and condition variables.
  //!
  //! This class extends `ThreadPoolBase` and manages a collection of worker threads
  //! that execute tasks from a queue.
  //!
  class ThreadPool4 : public ThreadPoolBase {


    using real_type = double;              //!< Type used for timing measurements.
    using TaskData  = tp::Queue::TaskData; //!< Type representing a task in the queue.

    std::atomic<bool>        m_done;                //!< Flag indicating if the thread pool has completed execution.
    std::atomic<unsigned>    m_running_task;        //!< Count of tasks currently running.
    std::atomic<unsigned>    m_running_thread;      //!< Count of threads currently active.
    std::vector<std::thread> m_worker_threads;      //!< Collection of worker threads.
    tp::Queue                m_work_queue;          //!< Queue for holding tasks; not thread-safe.
    std::mutex               m_work_on_queue_mutex; //!< Mutex for safe access to the task queue.

    // -----------------------------------------
    std::condition_variable_any m_queue_pop_cv; //!< Condition variable for task popping notifications.
    std::atomic<unsigned>       m_pop_waiting;  //!< Count of threads waiting to pop tasks from the queue.
    // -----------------------------------------
    std::condition_variable_any m_queue_push_cv; //!< Condition variable for task pushing notifications.
    std::atomic<unsigned>       m_push_waiting;  //!< Count of threads waiting to push tasks into the queue.
    // -----------------------------------------

    TicToc                   m_tm;      //!< Timing utility for measuring execution durations.
    std::vector<real_type>   m_job_ms;  //!< Duration of executed jobs.
    std::vector<real_type>   m_pop_ms;  //!< Duration for popping tasks from the queue.
    std::vector<unsigned>    m_n_job;   //!< Count of jobs executed.
    real_type                m_push_ms; //!< Duration for pushing tasks into the queue.

    //!
    //! \brief Sleep for a brief period to yield control.
    //!
    //! On Windows, it calls Sleep(0); otherwise, it sleeps for one nanosecond.
    //!
    inline
    void
    nano_sleep() const
    #ifdef UTILS_OS_WINDOWS
    { Sleep(0); }
    #else
    { sleep_for_nanoseconds(1); }
    //{ std::this_thread::yield(); }
    #endif

    /*!
     * \brief Pops a task from the task queue.
     *
     * \return A pointer to the popped task data, or nullptr if the queue is empty.
     */
    TaskData * pop_task();

    //!
    //! \brief Pushes a task onto the task queue.
    //!
    //! \param task A pointer to the task data to be pushed onto the queue.
    //!
    void push_task( TaskData * task );

    //!
    //! \brief The main function executed by worker threads.
    //!
    //! This function continuously pops tasks from the queue and executes them.
    //!
    //! \param pop_ms Reference to a variable for measuring the duration of task popping.
    //! \param job_ms Reference to a variable for measuring the duration of job execution.
    //! \param n_job Reference to a variable tracking the number of executed jobs.
    //!
    void
    worker_thread(
      real_type & pop_ms,
      real_type & job_ms,
      unsigned  & n_job
    );

    //!
    //! \brief Creates worker threads for the thread pool.
    //!
    //! \param thread_count The number of worker threads to create.
    //!
    void create_workers( unsigned thread_count );

  public:

    //!
    //! \brief Constructor for the ThreadPool4 class.
    //!
    //! \param thread_count The number of threads in the pool. Defaults to the number of hardware threads available.
    //! \param queue_capacity The capacity of the task queue. Defaults to 0, allowing unlimited tasks.
    //!
    explicit
    ThreadPool4(
      unsigned thread_count   = std::thread::hardware_concurrency(),
      unsigned queue_capacity = 0
    );

    //!
    //! \brief Destructor for the ThreadPool4 class.
    //!
    //! Ensures all tasks are completed before destruction by calling join().
    //!
    virtual ~ThreadPool4() { join(); }

    //!
    //! \brief Resizes the thread pool to a new thread count.
    //!
    //! \param thread_count The new number of threads to maintain in the pool.
    //!
    void resize( unsigned thread_count ) override { resize( thread_count, 0 ); }

    //!
    //! \brief Resizes the thread pool and the task queue.
    //!
    //! \param thread_count The new number of threads to maintain in the pool.
    //! \param queue_capacity The new capacity for the task queue.
    //!
    void resize( unsigned thread_count, unsigned queue_capacity );

    //!
    //! \brief Executes a task and adds it to the task queue.
    //!
    //! \param fun The function to be executed as a task.
    //!
    void
    exec( std::function<void()> && fun ) override
    { push_task( new TaskData(std::move(fun)) ); }

    //!
    //! \brief Waits for all tasks in the queue to complete.
    //!
    //! This function ensures that all queued tasks have finished executing.
    //!
    void wait() override;

    //!
    //! \brief Joins all worker threads and waits for them to finish executing.
    //!
    void join() override;

    //!
    //! \brief Provides information about the thread pool's performance.
    //!
    //! \param s The output stream to which information will be written.
    //!
    void info( ostream_type & s ) const override;

    //!
    //! \brief Gets the current number of threads in the pool.
    //!
    //! \return The number of threads in the pool.
    //!
    unsigned thread_count() const override { return unsigned(m_worker_threads.size()); }

    //!
    //! \brief Gets the current capacity of the task queue.
    //!
    //! \return The capacity of the task queue.
    //!
    unsigned queue_capacity() const { return m_work_queue.capacity(); }

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    char const * name() const override { return "ThreadPool4"; }
  };

  /*! @} */

}

//
// eof: ThreadPool4.hxx
//
