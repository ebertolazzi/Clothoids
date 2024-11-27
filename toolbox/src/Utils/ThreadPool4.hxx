/**
 */

///
/// file: ThreadPool4.hxx
///

namespace Utils {

  class ThreadPool4 : public ThreadPoolBase {

    using real_type = double;
    using TaskData  = tp::Queue::TaskData;

    std::atomic<bool>        m_done;
    std::atomic<unsigned>    m_running_task;
    std::atomic<unsigned>    m_running_thread;
    std::vector<std::thread> m_worker_threads;
    tp::Queue                m_work_queue; // not thread safe
    SpinLock                 m_work_on_queue_mutex;

    // -----------------------------------------
    std::condition_variable_any m_queue_pop_cv;
    std::atomic<unsigned>       m_pop_waiting;
    // -----------------------------------------
    std::condition_variable_any m_queue_push_cv;
    std::atomic<unsigned>       m_push_waiting;
    // -----------------------------------------

    TicToc                   m_tm;
    std::vector<real_type>   m_job_ms;
    std::vector<real_type>   m_pop_ms;
    std::vector<unsigned>    m_n_job;
    real_type                m_push_ms;

    inline
    void
    nano_sleep() const
    #ifdef UTILS_OS_WINDOWS
    { Sleep(0); }
    #else
    { sleep_for_nanoseconds(1); }
    //{ std::this_thread::yield(); }
    #endif

    TaskData * pop_task();
    void push_task( TaskData * task );

    void
    worker_thread(
      real_type & pop_ms,
      real_type & job_ms,
      unsigned  & n_job
    );

    void create_workers( unsigned thread_count );

  public:

    explicit
    ThreadPool4(
      unsigned thread_count   = std::thread::hardware_concurrency(),
      unsigned queue_capacity = 0
    );

    virtual ~ThreadPool4() { join(); }

    void resize( unsigned thread_count ) override { resize( thread_count, 0 ); }
    void resize( unsigned thread_count, unsigned queue_capacity );

    void
    exec( std::function<void()> && fun ) override
    { push_task( new TaskData(std::move(fun)) ); }

    void wait() override;

    void join() override;
    void info( ostream_type & s ) const override;

    unsigned thread_count()   const override { return unsigned(m_worker_threads.size()); }
    unsigned queue_capacity() const          { return m_work_queue.capacity(); }

    char const * name() const override { return "ThreadPool4"; }
  };

}

///
/// eof: ThreadPool4.hxx
///
