/**
 * Threadpool for C++11
 *
 * \copyright 2021 Enrico Bertolazzi
 *
 * based on the work of Ruediger Helsch, Ruediger.Helsch@t-online.de (2014)
 * version 2.0 (https://github.com/RuedigerHelsch/ThreadPool)
 */

//
// file: ThreadPool2.hxx
//

namespace Utils {

  /*!
   * \addtogroup THREAD
   * @{
   */

  //!
  //! VirtualTask type for thread pool
  //!
  //! The thread pool wraps functions into object of type VirtualTask and enqueues them.
  //! The actual tasks can be heterogeneous and must only support the VirtualTask interface.
  //!
  class VirtualTask {
  public:
    //!
    //! The payload, users function to be run.
    //!
    //! Operator() is run from the thread pool.
    //! Is responsible for deleting the task object once it is done.
    //!
    virtual void operator()() = 0;

    //! Destroy the task object
    virtual ~VirtualTask() = default;
  };

  //!
  //! Queue for functions with signature void()
  //!
  //! This queue is dependent on template parameter Function. Only one type of functions can be queued.
  //! For example, if this queue is instantiated for a lambda function, only this exact lambda function
  //! can be queued, since each lambda function has its own type separate from all other lambda functions.
  //!
  //! To make this queue more flexible, instantiate it either with std::function<void()>,
  //! or use the virtual thread pool VirtualThreadPool.
  //!
  //! \tparam Function
  //!         The function type to queue.
  //!
  class HQueue {

    //!
    //! Store pointers into the queue.
    //! Decorate the pointers with an operator() to make them
    //! callable as needed by ThreadPool.
    //!
    class QueueElement {
      VirtualTask * m_task;

    public:

      QueueElement()                                     = delete;
      QueueElement( QueueElement const & )               = delete;
      QueueElement & operator = ( QueueElement const & ) = delete;
      QueueElement & operator = ( QueueElement && )      = delete;

      //!
      //! Constructor that initializes the QueueElement with a VirtualTask pointer
      //!
      //! \param t Pointer to the VirtualTask to be stored.
      //!
      QueueElement( VirtualTask * t ) : m_task(t) { }

      //!
      //! Move constructor for QueueElement
      //!
      //! \param x Rvalue reference to another QueueElement
      //!
      QueueElement( QueueElement && x ) noexcept : m_task(x.m_task) { x.m_task = nullptr; }

      //! Call operator to execute the stored task
      void operator()() { (*m_task)(); m_task = nullptr; }

      //! Destructor that deletes the task if it exists
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

      //!
      //! Constructor to initialize the FixedCapacityQueue
      //!
      //! \param capacity The maximum number of tasks that can be stored in the queue.
      //!
      explicit
      FixedCapacityQueue( unsigned capacity )
      : m_fun_vec( size_t( capacity+1 ) )
      , m_size( capacity+1 )
      , m_capacity( capacity )
      { }

      //!
      //! Push a QueueElement into the queue
      //!
      //! \param f The QueueElement to be pushed into the queue.
      //!
      void
      push( QueueElement && f ) {
        new (&m_fun_vec[m_push_ptr].m_fun) QueueElement(std::forward<QueueElement>(f));
        if ( ++m_push_ptr == m_size ) m_push_ptr = 0;
      }

      //!
      //! Pop a QueueElement from the queue
      //!
      //! \return The popped QueueElement.
      //!
      QueueElement
      pop() {
        QueueElement r{ std::move(m_fun_vec[m_pop_ptr].m_fun) };
        m_fun_vec[m_pop_ptr].m_fun.~QueueElement();
        if ( ++m_pop_ptr == m_size ) m_pop_ptr = 0;
        return r;
      }

      //!
      //! Get the number of elements currently in the queue
      //!
      //! \return The size of the queue.
      //!
      unsigned
      size() const
      { return ((m_push_ptr + m_size) - m_pop_ptr) % m_size; }

      //!
      //! Check if the queue is empty
      //!
      //! \return True if the queue is empty, otherwise false.
      //!
      bool
      empty() const
      { return m_push_ptr == m_pop_ptr; }

      //!
      //! Check if the queue is full
      //!
      //! \return True if the queue is full, otherwise false.
      //!
      bool
      is_full() const
      { return this->size() >= m_capacity; }

      //!
      //! Get the maximum capacity of the queue
      //!
      //! \return The capacity of the queue.
      //!
      unsigned
      capacity() const
      { return m_capacity; }

      void reserve( unsigned capacity );

      //! Destructor for FixedCapacityQueue
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

    unsigned const          m_maxpart;                   //!< Maximum number of partitions
    bool                    m_shutting_down{false};      //!< Indicates if the queue is shutting down
    unsigned                m_idle_workers{0};           //!< Number of idle worker threads
    unsigned                m_total_workers{0};          //!< Total number of worker threads
    bool                    m_wakeup_is_pending{false};  //!< Indicates if a wake-up is pending
    FixedCapacityQueue      m_queue;                     //!< Queue to store tasks
    std::mutex              m_pop_mutex;                 //!< Mutex for popping tasks
    std::mutex              m_push_mutex;                //!< Mutex for pushing tasks
    std::condition_variable m_waiting_workers_cond;      //!< Condition variable for waiting workers
    std::condition_variable m_waiters_cond;              //!< Condition variable for waiters

    //! Get tasks and execute them. Return as soon as the queue shrinks to `return_if_idle` tasks.
    void help( int return_if_idle );

    //! Help, and shut down if an exception escapes.
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

    //!
    //! Constructor for HQueue
    //!
    //! \param queue_size The size of the queue.
    //! \param maxpart The maximum number of partitions.
    //!
    HQueue( unsigned queue_size, unsigned maxpart )
    : m_maxpart(maxpart)
    , m_queue(queue_size)
    { }

    //! Get tasks and execute them. If `return_if_idle`, return instead of idly waiting.
    void
    work( bool return_if_idle ) {
      help( return_if_idle ? 0 : -1 );
    }

    //!
    //! Enqueue a task.
    //!
    //! \tparam C The type of the callable object to enqueue.
    //!
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

    //! Shut down the queue and all associated worker threads
    void shutdown();

    //! Wait for all tasks to complete
    void wait();

    //!
    //! Get the capacity of the queue
    //!
    //! \return The queue's capacity.
    //!
    unsigned queue_capacity() const { return m_queue.capacity(); }

    //!
    //! Get the maximum number of partitions
    //!
    //! \return The maximum number of partitions.
    //!
    unsigned maxpart()        const { return m_maxpart; }

  };

  //!
  //! Implementation of thread pool.
  //!
  //! This class manages a pool of threads that can execute tasks concurrently.
  //!
  class ThreadPool2 : public ThreadPoolBase  {

    HQueue                   m_queue;          //!< The queue holding tasks
    std::vector<std::thread> m_worker_threads; //!< Vector of worker threads

  public:

    //!
    //! Constructor to initialize the ThreadPool2
    //!
    //! \param thread_count The number of threads to create (defaults to the number of hardware threads).
    //!
    explicit
    ThreadPool2( unsigned thread_count = std::thread::hardware_concurrency() );

    //! Execute a task
    //!
    //! \param fun The function to execute, wrapped in std::function<void()>.
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

    //! Wait for all tasks to finish executing
    void wait() override;

    //!
    //! Discard all tasks from the queue that have not yet started and wait for all threads to return.
    //! Also throws an exception if one of the tasks has encountered an uncatched exception.
    //! Leaves the pool in a shutdown state not ready to run tasks, but ready for destruction.
    //!
    void join() override;

    /**
     * Destroy the thread pool.
     *
     * Does the equivalent of wait() and join() before the thread pool is destructed.
     * This means, the destructor can hang a long time and can throw an exception
     * (unless wait() or join() have been called before the destructor).
     */

    //!
    //! Get the number of worker threads
    //!
    //! \return The number of threads in the pool.
    //!
    unsigned
    thread_count() const override
    { return unsigned(m_worker_threads.size()); }

    //!
    //! Get the capacity of the queue
    //!
    //! \return The capacity of the queue.
    //!
    unsigned
    queue_capacity() const
    { return m_queue.queue_capacity(); }

    //!
    //! Get the maximum number of partitions
    //!
    //! \return The maximum number of partitions.
    //!
    unsigned maxpart() const
    { return m_queue.maxpart(); }

    //!
    //! Resize the thread pool
    //!
    //! \param thread_count The new number of threads in the pool.
    //!
    void
    resize( unsigned thread_count ) override
    { this->resize( thread_count, 0, 0 ); }

    //!
    //! Resize the thread pool with specified queue size and maximum partitions
    //!
    //! \param thread_count The new number of threads in the pool.
    //! \param queue_size The new size of the queue.
    //! \param maxpart The new maximum number of partitions.
    //!
    void resize( unsigned thread_count, unsigned queue_size, unsigned maxpart );

    //!
    //! Get the name of the thread pool
    //!
    //! \return The name of the thread pool.
    //!
    char const * name() const override { return "ThreadPool2"; }

    //! Destructor for ThreadPool2
    virtual ~ThreadPool2();

  };

  /*! @} */

}

//
// eof: ThreadPool2.hxx
//
