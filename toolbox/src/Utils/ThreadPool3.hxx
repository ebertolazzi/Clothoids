//
// file: ThreadPool3.hxx
//

namespace Utils {

  /*!
   * \addtogroup THREAD
   * @{
   */

  //!
  //! \brief A thread pool implementation for concurrent task execution.
  //!
  //! The `ThreadPool3` class provides a high-level interface for managing a pool of threads
  //! that can execute tasks concurrently. It supports dynamic resizing of the thread pool
  //! and provides mechanisms for task synchronization and performance measurement.
  //!
  //! This class is built on top of the `ThreadPoolBase` and utilizes atomic operations
  //! to manage the state of the thread pool safely across multiple threads.
  //!
  class ThreadPool3 : public ThreadPoolBase {

    using real_type = double;              //!< Type used for timing measurements.
    using TaskData  = tp::Queue::TaskData; //!< Type representing a task in the queue.

    std::atomic<bool>           m_done{false};       //!< Flag indicating if the pool is finished.
    std::atomic<unsigned>       m_running_task{0};   //!< Number of tasks currently being executed.
    std::atomic<unsigned>       m_running_thread{0}; //!< Number of threads currently running tasks.
    std::vector<std::thread>    m_worker_threads;    //!< Vector of worker threads.
    tp::Queue                   m_work_queue;        //!< Queue for tasks; not thread-safe by itself.
    // -----------------------------------------
    std::mutex                  m_queue_push_mutex; //!< Mutex for managing concurrent task pushes.
    std::condition_variable_any m_queue_push_cv;    //!< Condition variable for notifying task pushes.
    std::atomic<unsigned>       m_push_waiting{0};  //!< Count of threads waiting to push tasks.
    // -----------------------------------------
    std::mutex                  m_queue_pop_mutex; //!< Mutex for managing concurrent task pops.
    std::condition_variable_any m_queue_pop_cv;    //!< Condition variable for notifying task pops.
    std::atomic<unsigned>       m_pop_waiting{0};  //!< Count of threads waiting to pop tasks.
    // -----------------------------------------
    UTILS_SPINLOCK              m_queue_spin; //!< Spinlock for quick access to the queue.

    TicToc                 m_tm;      //!< Timing utility for measuring task durations.
    std::vector<real_type> m_job_ms;  //!< Duration of job executions.
    std::vector<real_type> m_pop_ms;  //!< Duration of task popping.
    std::vector<unsigned>  m_n_job;   //!< Number of jobs executed.
    real_type              m_push_ms; //!< Duration of task pushing.

    //!
    //! \brief Pops a task from the task queue.
    //!
    //! \return Pointer to the popped task data, or nullptr if the queue is empty.
    //!
    TaskData * pop_task();

    //!
    //! \brief Pushes a task onto the task queue.
    //!
    //! \param task Pointer to the task data to be pushed onto the queue.
    //!
    void push_task( TaskData * task );

    //!
    //! \brief The main function for worker threads.
    //!
    //! This function continuously pops tasks from the queue and executes them.
    //!
    //! \param pop_ms Reference to the variable for measuring pop duration.
    //! \param job_ms Reference to the variable for measuring job execution duration.
    //! \param n_job Reference to the variable tracking the number of executed jobs.
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
    //! \brief Constructor for the ThreadPool3 class.
    //!
    //! \param thread_count The number of threads in the pool. Defaults to the number of hardware threads available.
    //! \param queue_capacity The capacity of the task queue. Defaults to 0, which allows unlimited tasks.
    //!
    explicit
    ThreadPool3(
      unsigned thread_count   = std::thread::hardware_concurrency(),
      unsigned queue_capacity = 0
    );

    //!
    //! \brief Destructor for the ThreadPool3 class.
    //!
    //! Cleans up resources and ensures that all tasks are completed before destruction by calling join().
    //!
    virtual ~ThreadPool3() { join(); }

    //!
    //! \brief Executes a task and adds it to the task queue.
    //!
    //! \param fun The function to be executed as a task.
    //!
    void
    exec( std::function<void()> && fun ) override {
      m_tm.tic();
      push_task( new TaskData(std::move(fun)) );
      m_tm.toc();
      m_push_ms += m_tm.elapsed_ms();
    }

    //!
    //! \brief Waits for all tasks in the queue to complete.
    //!
    //! This function yields control to allow other threads to proceed until all tasks are done.
    //!
    void
    wait() override
    { while ( !m_work_queue.empty() || m_running_task > 0 ) std::this_thread::yield(); }

    //!
    //! \brief Joins all worker threads and waits for them to finish.
    //!
    void join() override;

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
    unsigned thread_count()   const override { return unsigned(m_worker_threads.size()); }

    //!
    //! \brief Gets the current capacity of the task queue.
    //!
    //! \return The capacity of the task queue.
    //!
    unsigned queue_capacity() const          { return m_work_queue.capacity(); }

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    char const * name() const override { return "ThreadPool3"; }
  };

  /*! @} */

}

//
// eof: ThreadPool3.hxx
//
