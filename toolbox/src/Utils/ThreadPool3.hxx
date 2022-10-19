/*!
 */

///
/// file: ThreadPool3.hxx
///

namespace Utils {

  class ThreadPool3 : public ThreadPoolBase {

    using real_type = double;
    using TaskData  = tp::Queue::TaskData;

    std::atomic<bool>           m_done;
    std::atomic<unsigned>       m_running_task;
    std::atomic<unsigned>       m_running_thread;
    std::vector<std::thread>    m_worker_threads;
    tp::Queue                   m_work_queue; // not thread safe
    // -----------------------------------------
    std::mutex                  m_queue_push_mutex;
    std::condition_variable_any m_queue_push_cv;
    unsigned                    m_push_waiting;
    // -----------------------------------------
    std::mutex                  m_queue_pop_mutex;
    std::condition_variable_any m_queue_pop_cv;
    unsigned                    m_pop_waiting;
    // -----------------------------------------
    UTILS_SPINLOCK              m_queue_spin;

    TicToc                 m_tm;
    std::vector<real_type> m_job_ms;
    std::vector<real_type> m_pop_ms;
    std::vector<unsigned>  m_n_job;
    real_type              m_push_ms;

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
    ThreadPool3(
      unsigned thread_count   = std::thread::hardware_concurrency(),
      unsigned queue_capacity = 0
    );

    virtual ~ThreadPool3() { join(); }

    void
    exec( std::function<void()> && fun ) override {
      m_tm.tic();
      push_task( new TaskData(std::move(fun)) );
      m_tm.toc();
      m_push_ms += m_tm.elapsed_ms();
    }

    void
    wait() override
    { while ( !m_work_queue.empty() || m_running_task > 0 ) std::this_thread::yield(); }

    void join() override;
    void resize( unsigned thread_count ) override { resize( thread_count, 0 ); }
    void resize( unsigned thread_count, unsigned queue_capacity );

    void info( ostream_type & s ) const override;

    unsigned thread_count()   const override { return unsigned(m_worker_threads.size()); }
    unsigned queue_capacity() const          { return m_work_queue.capacity(); }

    char const * name() const override { return "ThreadPool3"; }
  };

}

///
/// eof: ThreadPool3.hxx
///
