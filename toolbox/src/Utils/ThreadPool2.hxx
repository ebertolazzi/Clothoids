/**
 * Threadpool for C++11
 *
 * \copyright 2021 Enrico Bertolazzi
 *
 * based on the work of Ruediger Helsch, Ruediger.Helsch@t-online.de (2014)
 * version 2.0 (https://github.com/RuedigerHelsch/ThreadPool)
 */

///
/// file: ThreadPool2.hxx
///

namespace Utils {

  /**
   * Queue for functions with signature void()
   *
   * This queue is dependent on template parameter
   * Function. Only one type of functions can be queued. For
   * example if this queue is instantiated for a lambda
   * function, only this exact lambda function can be
   * queued, since each lambda function has its own type
   * separate from all other lambda functions.
   *
   * To make this queue more flexible, instantiate it either
   * with std::function<void()>, or use the virtual thread
   * pool VirtualThreadPool.
   *
   * \tparam Function
   *         The function type to queue.
   */
  template <class Function>
  class HQueue {

    /*\
      This queue requires attention for protection against
      concurrent access. Protect against:

      - Concurrent access by two worker threads both wanting to get()
        a task from the queue at the same time.

      - Concurrent access by two threads both wanting to
        put() a task into the queue at the same time.

      - A worker thread having determined that the queue is empty,
        while at the same time a new task is put() into the queue.

      - A task wanting to put() a task into the queue having found the
        queue full, while at the same time the queues fill level decreases.
    \*/

    unsigned const                   m_maxpart;
    bool                             m_shutting_down;
    unsigned                         m_idle_workers;
    unsigned                         m_total_workers;
    bool                             m_wakeup_is_pending;
    tp::FixedCapacityQueue<Function> m_queue;
    std::mutex                       m_pop_mutex;
    std::mutex                       m_push_mutex;
    std::condition_variable          m_waiting_workers_cond;
    std::condition_variable          m_waiters_cond;

    /**
     * Get tasks and execute them. Return as soon as the queue
     * shrinks to `return_if_idle` tasks.
     */
    void
    help(int return_if_idle) {

      unsigned min_queue_size = return_if_idle < 0 ? 0 : return_if_idle;

      // Increment total worker count, decrement again on scope exit
      { std::lock_guard<std::mutex> lock(m_push_mutex); ++m_total_workers; }

      // execute at exit
      auto x1 = at_scope_exit([this](){
        std::lock_guard<std::mutex> lock(this->m_push_mutex);
        if (--this->m_total_workers == this->m_idle_workers)
          this->m_waiters_cond.notify_all();
      });

      tp::FixedCapacityQueue<Function> function_queue(1);

      for (;;) {
        std::unique_lock<std::mutex> lock(m_pop_mutex);
        unsigned queue_size;

        // Try to get the next task(s)
        while ((queue_size = m_queue.size()) <= min_queue_size) {
          if ( int(queue_size) <= return_if_idle) return;
          if ( queue_size > 0 ) break;
          // The queue is empty, wait for more tasks to be put()
          lock.unlock();
          {
            std::unique_lock<std::mutex> lock2(m_push_mutex);
            while (m_queue.empty() && !m_shutting_down) {
              if ( ++m_idle_workers == m_total_workers ) m_waiters_cond.notify_all();
              m_waiting_workers_cond.wait(lock2); // Wait for task to be queued
              m_wakeup_is_pending = false;
              --m_idle_workers;
            }
          }
          if (m_shutting_down) return;
          lock.lock();
        }

        // There is at least one task in the queue and the back is locked.

        unsigned stride = (m_maxpart == 0) ? 1 : queue_size / m_maxpart;
        if (stride <= 0) stride = 1;
        if (stride > function_queue.capacity()) function_queue.reserve(2 * stride);
        while (stride--) function_queue.push(m_queue.pop());
        lock.unlock();

        if ( m_idle_workers && !m_wakeup_is_pending && queue_size )
          m_waiting_workers_cond.notify_one();

        while (!function_queue.empty()) function_queue.pop()();
      }
    }

    /**
     * Help, and shut down if an exception escapes.
     */
    void
    try_help(int return_if_idle) {
      try {
        help(return_if_idle);
      } catch (...) {
        shutdown();
        throw;
      }
    }

  public:

    HQueue( unsigned queue_size, unsigned maxpart )
    : m_maxpart(maxpart)
    , m_shutting_down(false)
    , m_idle_workers(0)
    , m_total_workers(0)
    , m_wakeup_is_pending(false)
    , m_queue(queue_size)
    { }

    /**
     * Get tasks and execute them. If `return_if_idle`, return
     * instead of idly waiting.
     */
    void
    work( bool return_if_idle ) {
      help( return_if_idle ? 0 : -1 );
    }

    /**
       Enqueue a task.
    */
    template<class C>
    void
    put( C&& c ) {
      std::unique_lock<std::mutex> lock(m_push_mutex);
      while ( m_queue.is_full() ) {
        // No space in the queue. Must wait for workers to advance.
        lock.unlock();
        try_help( m_queue.capacity() / 2 );
        lock.lock();
      }
      // Now there is space in the queue and we have locked the back.

      // Enqueue function.
      if ( m_shutting_down ) {
        Function fun(std::forward<C>(c)); // Run Function destructor
      } else {
        // Here we have exclusive access to the head of the queue.
        m_queue.push(std::forward<C>(c));

        if ( m_idle_workers && !m_wakeup_is_pending ) {
          m_wakeup_is_pending = true;
          m_waiting_workers_cond.notify_one(); // wake up a worker
        }
      }
    }

    void
    shutdown() {
      std::unique_lock<std::mutex> push_lock(m_push_mutex);
      std::unique_lock<std::mutex> pop_lock(m_pop_mutex);
      m_shutting_down = true;
      while (!m_queue.empty()) m_queue.pop();
      m_waiting_workers_cond.notify_all();
      m_waiters_cond.notify_all();
    }

    void
    wait() {
      if ( std::uncaught_exception() ) shutdown();
      std::exception_ptr e;
      std::unique_lock<std::mutex> lock(m_push_mutex);
      while ( !m_queue.empty() || m_idle_workers != m_total_workers ) {
        while ( !m_queue.empty() ) {
          lock.unlock();
          try {
            try_help(0);
          } catch (...) {
            if (e == nullptr) e = std::current_exception();
          }
          lock.lock();
        }
        while ( m_idle_workers != m_total_workers ) m_waiters_cond.wait(lock);
      }
      if (e != nullptr && !std::uncaught_exception())
        std::rethrow_exception(std::move(e));
    }

    unsigned queue_capacity() const { return m_queue.capacity(); }
    unsigned maxpart()        const { return m_maxpart; }

  };

  /**
   * VirtualTask type for thread pool
   *
   * The thread pool wraps functions into object of type
   * VirtualTask and enqueues them. The actual tasks can be
   * heterogenous and must only support the VirtualTask interface.
   */
  class VirtualTask {
  public:
    /**
     * The payload, users function to be run.
     *
     * Operator() is run from the thread pool.
     * Is responsible for deleting the task object once it is done.
     */
    virtual void operator()() = 0;

    /**
     * Destroy the task object
     */
    virtual ~VirtualTask() = default;
  };

  /**
   * A thread pool reading its tasks from a generic queue.
   *
   * \tparam class Queue
   *         The queue delivering the tasks.
   *
   *	Class Queue must provide the following members:
   *
   *	- void work()
   *    Gets tasks and works until the end of
   *    the queue is reached.
   *
   *	- void shutdown()
   *    Causes the queue to tell the threads
   *    asking for work to return.
   *
   */
  template<class Queue>
  class GenericThreadPool {

    std::mutex               m_mutex;
    Queue *                  m_queue;
    std::vector<std::thread> m_worker_threads;

    //! The main function of the thread.
    void work() { do_work(false); }

  public:

    // Copying and moving are not supported.
    GenericThreadPool( GenericThreadPool const & )             = delete;
    GenericThreadPool( GenericThreadPool && )                  = delete;
    GenericThreadPool& operator = (GenericThreadPool const & ) = delete;
    GenericThreadPool& operator = (GenericThreadPool && )      = delete;

    /**
     * Generic thread pool.
     *
     * As soon as the pool is created the threads start running
     * tasks from the queue until the queue is empty. When the
     * queue is empty the threads return and are ready to be
     * collected by join() or by the destructor.
     *
     * \param queue
     *        The queue delivering the tasks.
     *
     * \param thread_count
     *        The number of threads to use.
     *        If the thread count is not specified it defaults
     *        to the number of available hardware threads
     *        std::thread::hardware_concurrency().
     *
     */
    GenericThreadPool( Queue * queue, int thread_count )
    : m_queue(queue)
    , m_worker_threads(thread_count)
    {
      for ( std::thread & w : m_worker_threads )
        w = std::thread( std::bind( &GenericThreadPool::work, this ) );
    }

    /**
     * Help with the work.
     *
     * \param return_if_idle
     *        Never wait for work, return instead.
     */
    void do_work( bool return_if_idle ) { m_queue->work( return_if_idle ); }

    /**
     * Wait for all threads to finish and collect them.
     *
     * Leaves the thread pool ready for destruction.
     */
    void
    join() {
      work(); // Instead of hanging around, help the workers!
      for ( std::thread & w : m_worker_threads )
        { if (w.joinable()) w.join(); }
    }

    /**
     * Destroy the thread pool.
     *
     * Generally a destructor should not wait for a long time,
     * and it should not throw any exceptions. Unfortunately
     * threads are not abortable in C++.  The only way to make
     * sure the threads have terminated is to wait for them to
     * return by themselves.  If we would allow the
     * destruction of the thread pool before the threads have
     * returned, the threads would continue to access the
     * memory of the destroyed thread pool, potentially
     * clobbering other objects residing in the recycled
     * memory.  We could allocate parts of the memory with
     * new, and leave it behind for the threads after the
     * thread pool is destructed.  But even then, the user
     * supplied functions run by the threads might access
     * memory that gets destroyed if the function that
     * constructed the thread pool terminates.  The danger of
     * undetected and undebuggable memory corruption is just
     * too big.
     *
     * With regard to the exceptions rethrown in the
     * destructor, it is better to signal the exception than
     * to ignore it silently.
     *
     * If it is not acceptable for the destructor to wait or
     * to throw an exception, just call join() before the pool
     * is destructed.  After join() the destructor is
     * guaranteed to run fast and without exceptions.
     *
     * If it should really be necessary to keep threads
     * running after the function that created the thread pool
     * returns, just create the thread pool on the heap with
     * new. And if you want to make sure nobody destroys the
     * thread pool, feel free to throw away the handle.
     */
    ~GenericThreadPool() {
      // Abort processing if destructor runs during exception handling.
      if ( std::uncaught_exception() ) m_queue->shutdown();
      join(); // Running threads would continue to access the destructed pool.
    }

    unsigned
    thread_count() const
    { return unsigned(m_worker_threads.size()); }
  };

  /**
   * Store pointers into the queue. Decorate the pointers
   * with an operator() to make them callable as needed by
   * ThreadPool.
   */
  class QueueElement {
    VirtualTask * m_task;

  public:

    QueueElement()                                     = delete;
    QueueElement( QueueElement const & )               = delete;
    QueueElement & operator = ( QueueElement const & ) = delete;
    QueueElement & operator = ( QueueElement && )      = delete;

    QueueElement( VirtualTask * t ) : m_task(t) { }
    QueueElement( QueueElement && x ) noexcept : m_task(x.m_task) { x.m_task = nullptr; }
    void operator()() { (*m_task)(); m_task = nullptr; }
    ~QueueElement() { if (m_task) delete m_task; }
  };

  /**
   * Implementation of virtual thread pool.
   *
   * Implements the functionality of the virtual thread
   * pool. Only provides an interface to run a generic
   * VirtualTask. The convenience functions to run
   * different types of callable objects should be implemented
   * in a subclass.
   *
   * The template parameter is not used, only serves to make
   * this a class template which can be instantiated in multiple
   * compilation units without giving multiply defined symbol
   * errors.
   */
  class ThreadPool2 : public ThreadPoolBase  {

    using QUEUE = HQueue<QueueElement>;
    using POOL  = GenericThreadPool<QUEUE>;

    QUEUE * m_queue;
    POOL  * m_pool;

  public:

    explicit
    ThreadPool2(
      unsigned thread_count = std::thread::hardware_concurrency(),
      unsigned queue_size   = 0,
      unsigned maxpart      = 0
    ) {
      if ( queue_size == 0 ) queue_size = 50 * (thread_count+1);
      if ( maxpart    == 0 ) maxpart    = 3 * (thread_count+1);
      m_queue = new QUEUE( queue_size, maxpart );
      m_pool  = new POOL( m_queue, thread_count );
    }

    void
    run_task( std::unique_ptr<VirtualTask>&& t )
    { m_queue->put(t.release()); }

    void
    run_task( VirtualTask * t )
    { m_queue->put(t); }

    void
    exec( std::function<void()> && fun ) override {
      class WrappedFunction : public VirtualTask {
        std::function<void()> m_f;
      public:
        explicit
        WrappedFunction( std::function<void()> && f ) : m_f(std::move(f)) { }
        virtual void operator()() override { m_f(); delete this; }
      };
      run_task(new WrappedFunction( std::move(fun) ) );
    }

    void
    wait() override {
      m_pool->do_work(true); // Help out instead of sitting around idly.
      m_queue->wait();
    }

    /**
     * Discard all tasks from the queue that have not yet
     * started and wait for all threads to return.
     *
     * Also throws an exception if one of the tasks has
     * encountered an uncatched exception.
     *
     * Leaves the pool in a shutdown state not ready to run
     * tasks, but ready for destruction.
     */

    void
    join() override {
      wait();
      m_queue->shutdown();
      m_pool->join();
    }

    /**
     * Destroy the thread pool.
     *
     * Does the equivalent of wait() and join() before the
     * thread pool is destructed. This means, the destructor
     * can hang a long time and can throw an exception (unless
     * wait() or join() have been called before the
     * destructor).
     */

    unsigned thread_count()   const override { return m_pool->thread_count(); }
    unsigned queue_capacity() const          { return m_queue->queue_capacity(); }
    unsigned maxpart()        const          { return m_queue->maxpart(); }

    void
    resize( unsigned thread_count ) override {
      this->resize( thread_count, 0, 0 );
    }

    void
    resize(
      unsigned thread_count,
      unsigned queue_size,
      unsigned maxpart
    ) {
      this->~ThreadPool2();
      if ( queue_size == 0 ) queue_size = 50 * (thread_count+1);
      if ( maxpart    == 0 ) maxpart    = 3 * (thread_count+1);
      m_queue = new QUEUE( queue_size, maxpart );
      m_pool  = new POOL( m_queue, thread_count );
    }

    char const * name() const override { return "ThreadPool2"; }

    virtual
    ~ThreadPool2() {
      wait(); join();
      delete m_pool; delete m_queue;
    }
  };

}

///
/// eof: ThreadPool2.hxx
///
