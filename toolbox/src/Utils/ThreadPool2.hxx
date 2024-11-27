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
   * VirtualTask type for thread pool
   *
   * The thread pool wraps functions into object of type VirtualTask and enqueues them.
   * The actual tasks can be heterogenous and must only support the VirtualTask interface.
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
   * Queue for functions with signature void()
   *
   * This queue is dependent on template parameter Function. Only one type of functions can be queued.
   * For example if this queue is instantiated for a lambda function, only this exact lambda function
   * can be queued, since each lambda function has its own type separate from all other lambda functions.
   *
   * To make this queue more flexible, instantiate it either with std::function<void()>,
   * or use the virtual thread pool VirtualThreadPool.
   *
   * \tparam Function
   *         The function type to queue.
   */
  class HQueue {

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

    /*\
     |   ___ _            _  ___                    _ _         ___
     |  | __(_)_ _____ __| |/ __|__ _ _ __  __ _ __(_) |_ _  _ / _ \ _  _ ___ _  _ ___
     |  | _|| \ \ / -_) _` | (__/ _` | '_ \/ _` / _| |  _| || | (_) | || / -_) || / -_)
     |  |_| |_/_\_\___\__,_|\___\__,_| .__/\__,_\__|_|\__|\_, |\__\_\\_,_\___|\_,_\___|
     |                               |_|                  |__/
    \*/
    /*\
     * If we would use a deque, we would have to protect against overlapping accesses to
     * the front and the back. The standard containers do not allow this.
     * Use a vector instead.  With a vector it is possible to access both ends of the
     * queue at the same time, as push()ing and pop()ing does not modify the container
     * itself but only its elements.
    \*/
    class FixedCapacityQueue {

      union Fun {
        QueueElement m_fun; // Only used between pop_ptr and push_ptr
        Fun() noexcept { }
        Fun( Fun const & ) noexcept { }
        Fun( Fun && ) noexcept { }
        ~Fun() noexcept { }
      };

      std::vector<Fun> m_fun_vec;
      unsigned m_size     = 0;
      unsigned m_capacity = 0;
      unsigned m_push_ptr = 0;
      unsigned m_pop_ptr  = 0;

    public:

      FixedCapacityQueue( FixedCapacityQueue const & )              = delete;
      FixedCapacityQueue( FixedCapacityQueue && )                   = delete;
      FixedCapacityQueue& operator = ( FixedCapacityQueue const & ) = delete;
      FixedCapacityQueue& operator = ( FixedCapacityQueue && )      = delete;

      explicit
      FixedCapacityQueue( unsigned capacity )
      : m_fun_vec( size_t( capacity+1 ) )
      , m_size( capacity+1 )
      , m_capacity( capacity )
      { }

      void
      push( QueueElement && f ) {
        new (&m_fun_vec[m_push_ptr].m_fun) QueueElement(std::forward<QueueElement>(f));
        if ( ++m_push_ptr == m_size ) m_push_ptr = 0;
      }

      QueueElement
      pop() {
        QueueElement r{ std::move(m_fun_vec[m_pop_ptr].m_fun) };
        m_fun_vec[m_pop_ptr].m_fun.~QueueElement();
        if ( ++m_pop_ptr == m_size ) m_pop_ptr = 0;
        return r;
      }

      unsigned size()     const { return ((m_push_ptr + m_size) - m_pop_ptr) % m_size; }
      bool     empty()    const { return m_push_ptr == m_pop_ptr; }
      bool     is_full()  const { return this->size() >= m_capacity; }
      unsigned capacity() const { return m_capacity; }

      void reserve( unsigned capacity );

      ~FixedCapacityQueue();
    };

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

    unsigned const          m_maxpart;
    bool                    m_shutting_down{false};
    unsigned                m_idle_workers{0};
    unsigned                m_total_workers{0};
    bool                    m_wakeup_is_pending{false};
    FixedCapacityQueue      m_queue;
    std::mutex              m_pop_mutex;
    std::mutex              m_push_mutex;
    std::condition_variable m_waiting_workers_cond;
    std::condition_variable m_waiters_cond;

    /**
     * Get tasks and execute them. Return as soon as the queue shrinks to `return_if_idle` tasks.
     */
    void help( int return_if_idle );

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
        QueueElement fun(std::forward<C>(c)); // Run Function destructor
      } else {
        // Here we have exclusive access to the head of the queue.
        m_queue.push(std::forward<C>(c));

        if ( m_idle_workers && !m_wakeup_is_pending ) {
          m_wakeup_is_pending = true;
          m_waiting_workers_cond.notify_one(); // wake up a worker
        }
      }
    }

    void shutdown();
    void wait();

    unsigned queue_capacity() const { return m_queue.capacity(); }
    unsigned maxpart()        const { return m_maxpart; }

  };

  /**
   * Implementation of thread pool.
   */
  class ThreadPool2 : public ThreadPoolBase  {

    HQueue                   m_queue;
    std::vector<std::thread> m_worker_threads;

  public:

    explicit
    ThreadPool2( unsigned thread_count = std::thread::hardware_concurrency() );

    void
    exec( std::function<void()> && fun ) override {
      class WrappedFunction : public VirtualTask {
        std::function<void()> m_f;
      public:
        explicit
        WrappedFunction( std::function<void()> && f ) : m_f(std::move(f)) { }
        virtual void operator()() override { m_f(); delete this; }
      };
      m_queue.put( new WrappedFunction( std::move(fun) ) );
    }

    void wait() override;

    /**
     * Discard all tasks from the queue that have not yet started and wait for all threads to return.
     * Also throws an exception if one of the tasks has encountered an uncatched exception.
     * Leaves the pool in a shutdown state not ready to run tasks, but ready for destruction.
     */

    void join() override;

    /**
     * Destroy the thread pool.
     *
     * Does the equivalent of wait() and join() before the thread pool is destructed.
     * This means, the destructor can hang a long time and can throw an exception
     * (unless wait() or join() have been called before the destructor).
     */

    unsigned thread_count()   const override { return unsigned(m_worker_threads.size()); }
    unsigned queue_capacity() const          { return m_queue.queue_capacity(); }
    unsigned maxpart()        const          { return m_queue.maxpart(); }

    void resize( unsigned thread_count ) override { this->resize( thread_count, 0, 0 ); }
    void resize( unsigned thread_count, unsigned queue_size, unsigned maxpart );

    char const * name() const override { return "ThreadPool2"; }

    virtual ~ThreadPool2();

  };

}

///
/// eof: ThreadPool2.hxx
///
