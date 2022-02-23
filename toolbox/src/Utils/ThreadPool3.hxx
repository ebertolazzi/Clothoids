/*!
 */

///
/// file: ThreadPool3.hxx
///

namespace Utils {

  class ThreadPool3 : public ThreadPoolBase {

    typedef tp::Queue::TaskData TaskData;

    std::atomic<bool>        m_done;
    std::atomic<unsigned>    m_running_task;
    std::atomic<unsigned>    m_running_thread;
    std::vector<std::thread> m_worker_threads;
    tp::Queue                m_work_queue; // not thread safe
    SpinLock                 m_spin_push;
    SpinLock                 m_spin_pop;

    inline
    void
    nano_sleep() const {
      std::this_thread::sleep_for(std::chrono::nanoseconds(10));
    }

    TaskData *
    pop_task() {
      m_spin_pop.lock();
      while ( m_work_queue.empty() ) nano_sleep();
      TaskData * task = m_work_queue.pop();
      ++m_running_task; // must be incremented in the locked part
      m_spin_pop.unlock();
      return task;
    }

    void
    push_task( TaskData * task ) {
      m_spin_push.lock();
      while ( m_work_queue.is_full() ) nano_sleep();
      m_work_queue.push( task );
      m_spin_push.unlock();
    }

    void
    worker_thread() {
      ++m_running_thread;
      while ( !m_done ) {
        (*pop_task())();
        --m_running_task;
      }
      --m_running_thread;
    }

    void
    create_workers( unsigned thread_count ) {
      m_worker_threads.clear();
      m_worker_threads.reserve(thread_count);
      m_done = false;
      try {
        for ( unsigned i=0; i<thread_count; ++i )
          m_worker_threads.push_back(
            std::thread(&ThreadPool3::worker_thread,this)
          );
      } catch(...) {
        m_done = true;
        throw;
      }
    }

  public:

    explicit
    ThreadPool3(
      unsigned thread_count   = std::thread::hardware_concurrency(),
      unsigned queue_capacity = 0
    )
    : m_done(false)
    , m_running_task(0)
    , m_running_thread(0)
    , m_work_queue( queue_capacity == 0 ? 4 * (thread_count+1) : queue_capacity )
    {
      create_workers( thread_count );
    }

    virtual ~ThreadPool3() { join(); }

    void
    exec( std::function<void()> && fun ) override {
      push_task( new TaskData(std::move(fun)) );
    }

    void
    wait() override
    { while ( !m_work_queue.empty() || m_running_task > 0 ) nano_sleep(); }

    void
    join() {
      wait();
      m_done = true;
      { // send null task until all the workers stopped
        for ( unsigned i = m_running_thread; i > 0; --i ) {
          std::function<void()> null_job = [](){};
          push_task( new TaskData(std::move(null_job)) );
        }
        while ( m_running_thread > 0 ) nano_sleep();
      }
      m_work_queue.clear();
      for ( std::thread & w : m_worker_threads ) { if (w.joinable()) w.join(); }
      m_worker_threads.clear();
    }

    void resize( unsigned thread_count ) override { resize( thread_count, 0 ); }

    void
    resize( unsigned thread_count, unsigned queue_capacity = 0 ) {
      join();
      if ( queue_capacity == 0 ) queue_capacity = 4 * (thread_count+1);
      m_work_queue.resize( queue_capacity );
      create_workers( thread_count );
    }

    unsigned thread_count()   const override { return unsigned(m_worker_threads.size()); }
    unsigned queue_capacity() const          { return m_work_queue.capacity(); }

    char const * name() const override { return "ThreadPool3"; }
  };

}

///
/// eof: ThreadPool3.hxx
///
