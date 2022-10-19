namespace threadpool {

  /*
    The thread pool for arbitrary functions works fine, and can be
    used to process elements of a container. But this means queuing
    a task for each element, with each task executing the same
    function on another element of the container. Certainly there
    is the possibility for optimization.
    Define a queue calling the same function for each object in an
    iterator range. This means we do not need a true queue but just
    an object using incremental values of an iterator until the end
    of the range is reached.
  */

  /**
   * Queue calling the function on single objects.
   *
   * \relates ForEach_ThreadPool
   * Conceptually ForEach_Queue is a member
   * of class ForEach_ThreadPool, but the standard does
   * not allow template specialization inside classes. I
   * had to move it out of the class.
   */
  template<class Iterator, class Last, class Function, bool forward_iterator>
  class ForEach_Queue {
  protected:
    Iterator   & m_current;
    Last const & m_last;
    Function   & m_fun;
    std::mutex   m_mutex; // Make sure threads do not access concurrently

  public:

    ForEach_Queue(
      Iterator   & first,
      Last const & last,
      Function   & fun,
      unsigned = 0 /*ignored*/
    )
    : m_current(first)
    , m_last(last)
    , m_fun(fun)
    { }

    void
    work(bool /*ignored*/) {
      using IT = iterval_traits<Iterator>;
      Last const & l(m_last);
      for (;;) {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (m_current == l) break;
        typename IT::type v(IT::copy(m_current));
        ++m_current;
        lock.unlock();
        m_fun(IT::pass(std::move(v)));
      }
    }

    /**
     * Shut the queue down, stop returning values
     */
    void
    shutdown() {
      std::lock_guard<std::mutex> lock(m_mutex);
      m_current = m_last;
    }
  };

  /*
    The queue just implemented would work fine. But if there are
    a lot of tasks with each task taking a very short time, it may
    cause a lot of overhead because each object is dequeued
    separately. Wouldn't it be nice if we could deliver larger
    tasks?
  */

  /**
   * Run a function on objects from a container.
   *
   * Queue with `forward_iterator` == false takes groups of
   * objects from the queue.
   *
   * This works only for random access iterators. The
   * specialization is selected with template parameter
   * forward_iterator = true. For all other iterators, use the
   * general case of the template above.
   *
   * \relates ForEach_ThreadPool
   * Conceptually ForEach_Queue is a member
   * of class ForEach_ThreadPool, but the standard does
   * not allow template specialization inside classes. I
   * had to move it out of the class.
   */
  template<class Iterator, class Last, class Function>
  class ForEach_Queue<Iterator, Last, Function, true> : public ForEach_Queue<Iterator, Last, Function, false> {
    using Base            = ForEach_Queue<Iterator, Last, Function, false>;
    using difference_type = typename std::iterator_traits<Iterator>::difference_type;
    unsigned const  m_maxpart;
    difference_type m_remaining;
  public:

    ForEach_Queue(
      Iterator   & first,
      Last const & last,
      Function   & fun,
      unsigned     maxpart
    )
    : Base(first, last, fun)
    , m_maxpart(maxpart)
    , m_remaining(std::distance(first, last))
    { }

    void
    work(bool /*ignored*/) {
      Last const & last = this->m_last; // Does never change
      for (;;) {
        Iterator c, l;
        {
          std::lock_guard<std::mutex> lock(this->m_mutex);
          if ((c = this->m_current) == last) break;
          difference_type stride = (m_maxpart == 0) ? 1 : m_remaining / m_maxpart;
          if (stride <= 0) stride = 1;
          l = c;
          std::advance(l, stride);
          this->m_current = l;
          m_remaining -= stride;
        }
        while (c != l) { this->m_fun(*c); ++c; }
      }
    }
  };

  /*
    Now write a parallel version of std::transform().
    We could reuse a few parts from the parallel for_each
    implementation, but this would not save much and create a
    dependence. Better to keep parallel_transform self-standing.
  */

  /**
   * Run a function on objects from a container, store the return
   * values in a result container.
   *
   * General case for arbitrary iterators. This is the difficult
   * case to implement, because the worker threads must synchronize
   * after they have done their work, and must make sure they write
   * the results in the order of the input objects and not in the
   * order the threads are finished. But this is the thing that
   * makes parallel_transform interesting for the users. They don't
   * need to synchronize, the algorithm makes it for them.
   *
   * \tparam InputIterator
   *         Type of the input iterator. In this
   *         specialization with forward_iterator = false,
   *         an arbitrary input iterator type.
   *
   * \tparam OutputIterator
   *         Type of the result iterator. In this
   *         specialization with forward_iterator = false,
   *         an arbitrary output iterator type.
   * \tparam Function
   *         Type of the function to be called with
   *         successive elements from the input iterator.
   *         The function must return a result
   *         which is stored through the result iterator.
   *
   * \tparam forward_iterator
   *         A bool selecting the specialization.
   *         The general case for arbitrary input and output
   *         iterators which is implemented here is selected
   *         with forward_iterator = false. The specialization
   *         for forward iterators follows below.
   *
   * \relates Transform_ThreadPool
   *          Transform_Queue is conceptually a member of class
   *          Transform_ThreadPool, but the standard does not
   *          allow template specialization inside classes.
   *          I had to move it out of the class.
   */
  template<
    class InputIterator,
    class Last,
    class OutputIterator,
    class Function,
    bool forward_iterator
  >
  class Transform_Queue {
    struct Results {
      typename std::remove_reference<decltype(std::declval<Function&>()(*std::declval<InputIterator>()))>::type result;
      std::unique_ptr<Results> next;
    };

    InputIterator  & m_current;
    Last const     & m_last;
    OutputIterator & m_result;
    Function       & m_fun;
    std::mutex       m_mutex;
    bool             m_do_shutdown = false;

    using counter_type = unsigned long long int;

    counter_type            m_input_counter           = 1; // Counter of objects got from the queue
    counter_type            m_output_counter          = 1; // Counter of objects written
    Results *               m_previous_results        = nullptr;
    counter_type            m_max_output_queue_length = 1000; // This should be configurable
    std::mutex              m_output_mutex;
    std::condition_variable m_output_queue_cond;
    unsigned                m_output_queue_waiters = 0;

  public:

    Transform_Queue(
      InputIterator  & first,
      Last const     & last,
      OutputIterator & result,
      Function       & fun,
      unsigned
    )
    : m_current(first)
    , m_last(last)
    , m_result(result)
    , m_fun(fun)
    { }

    void
    work(bool return_if_idle) {
      using IT = iterval_traits<InputIterator>;

      std::unique_ptr<Results> results;
      Last const & last = m_last; // Does never change.
      for (;;) {
        if (!results) results = std::unique_ptr<Results>(new Results);
        counter_type ctr;
        Results* prvres;
        {
          std::unique_lock<std::mutex> lock(m_mutex);
          if (m_current == last) break;
          ctr                = m_input_counter;
          prvres             = m_previous_results;
          m_previous_results = &*results;
          typename IT::type v(IT::copy(m_current));
          ++m_current;
          m_input_counter = ctr + 1;
          lock.unlock();
          results->result = fun(IT::pass(std::move(v)));
        }
        {
          /*
            We must store the results in the order they had
            in the input sequence, not in the order the
            tasks finish. Just work together: whoever is
            ready before his predecessor just leaves his
            work for the predecessor to clean up.
          */
          std::unique_lock<std::mutex> lock(m_output_mutex);
          while (ctr - m_output_counter > m_max_output_queue_length) {
            if (m_do_shutdown) return;
            if (return_if_idle) {
              prvres->next = std::move(results);
              return;
            }
            ++m_output_queue_waiters;
            m_output_queue_cond.wait(lock);
            --m_output_queue_waiters;
          }
          if (m_output_counter == ctr) {
            // Predecessor is done, we can store our things.
            lock.unlock();
            *m_result = std::move(results->result);
            ++m_result;
            ++ctr;
            lock.lock();
            // Now look whether our successors have left us their work.
            while (results->next) {
              results = std::move(results->next);
              lock.unlock();
              *m_result = std::move(results->result);
              ++m_result;
              ++ctr;
              lock.lock();
            }
            m_output_counter = ctr;
            if (m_output_queue_waiters) m_output_queue_cond.notify_all(); // All because we do not know who is the right one.
          } else {
            // Predecessor still running, let him clean up.
            prvres->next = std::move(results);
          }
        }
      }
    }

    /**
     * Shut the queue down, stop returning values
     */
    void
    shutdown() {
      std::lock_guard<std::mutex> lock1(m_mutex);
      std::lock_guard<std::mutex> lock2(m_output_mutex);
      m_current     = m_last;
      m_do_shutdown = true;
      m_output_queue_cond.notify_all();
    }
  };

  /**
   * Run a function on objects from a container, store the return
   * values in a result container.
   *
   * Specialization for forward iterators. It is used when
   * template argument forward_iterator is true. For all other
   * iterators, use the generic version above.
   *
   * \tparam InputIterator
   * Type of the input iterator.
   * In this specialization with forward_iterator = false,
   * an arbitrary input iterator type.
   *
   * \tparam OutputIterator
   * Type of the result iterator.
   * In this specialization with forward_iterator = false,
   * an arbitrary output iterator type.
   *
   * \tparam Function
   * Type of the function to be called with
   * successive elements from the input iterator.
   * The function must return a result
   * which is stored through the result iterator.
   *
   * \relates Transform_ThreadPool
   * Transform_Queue is conceptually a member
   * of class Transform_ThreadPool, but the standard
   * does not allow template specialization inside classes.
   * I had to move it out of the class.
   */
  template<class InputIterator, class Last,class OutputIterator, class Function>
  class Transform_Queue<InputIterator, Last, OutputIterator, Function, true> {
    using difference_type = typename std::iterator_traits<InputIterator>::difference_type;
    InputIterator   & m_current;
    Last const      & m_last;
    OutputIterator  & m_result;
    Function        & m_fun;
    std::mutex        m_mutex;
    unsigned const    m_maxpart;
    difference_type   m_remaining;
  public:
    Transform_Queue(
      InputIterator  & first,
      Last const     & last,
      OutputIterator & result,
      Function       & fun,
      unsigned         maxpart
    )
    : m_current(first)
    , m_last(last)
    , m_result(result)
    , m_fun(fun)
    , m_maxpart(maxpart)
    , m_remaining(std::distance(first, last))
    { }

    void
    work(bool /* ignored */) {
      Last const & last = m_last; // Does never change
      for (;;) {
        InputIterator  c, l;
        OutputIterator r;
        {
          std::lock_guard<std::mutex> lock(m_mutex);
          if ((c = m_current) == last) break;
          difference_type stride = (m_maxpart == 0) ? 1 : m_remaining / m_maxpart;
          if (stride <= 0) stride = 1;
          l = c;
          std::advance(l, stride);
          r = m_result;
          m_current = l;
          std::advance(m_result, stride);
          m_remaining -= stride;
        }
        while (c != l) { *r = m_fun(*c); ++c; ++r; }
      }
    }

    /**
     * Shut the queue down, stop returning values
     */
    void
    shutdown() {
      std::lock_guard<std::mutex> lock(m_mutex);
      m_current = m_last;
    }
  };

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

    /*
      If we would use a deque, we would have to protect
      against overlapping accesses to the front and the
      back. The standard containers do not allow this. Use a
      vector instead.  With a vector it is possible to access
      both ends of the queue at the same time, as push()ing
      and pop()ing does not modify the container itself but
      only its elements.
    */
    class Queue {
      union Fun {
        Function m_fun; // Only used between pop_ptr and push_ptr
        Fun() noexcept { }
        Fun( Fun const & ) noexcept { }
        Fun( Fun && ) noexcept { }
        ~Fun() noexcept { }
      };
      std::vector<Fun> m_fun_vec;
      unsigned m_size     = 0;
      unsigned m_push_ptr = 0;
      unsigned m_pop_ptr  = 0;

      Queue( Queue const & )              = delete;
      Queue( Queue && )                   = delete;
      Queue& operator = ( Queue const & ) = delete;
      Queue& operator = ( Queue && )      = delete;

    public:

      Queue( unsigned s )
      : m_fun_vec( std::size_t(s+1) )
      , m_size(s+1) { }

      template<class F>
      void
      push( F && f ) {
        new (&m_fun_vec[m_push_ptr].m_fun) Function(std::forward<F>(f));
        if (++m_push_ptr == m_size) m_push_ptr = 0;
      }

      Function
      pop() {
        Function r = std::move(m_fun_vec[m_pop_ptr].m_fun);
        m_fun_vec[m_pop_ptr].m_fun.~Function();
        if (++m_pop_ptr == m_size) m_pop_ptr = 0;
        return r;
      }

      unsigned
      size() const
      { return ((m_push_ptr + m_size) - m_pop_ptr) % m_size; }

      bool     empty()    const { return m_push_ptr == m_pop_ptr; }
      unsigned capacity() const { return m_size - 1; }

      void
      reserve( unsigned s ) {
        assert(empty()); // Copying / moving of Fun not supported.
        if ( s >= m_size ) {
          m_fun_vec.resize(s + 1);
          m_size = s+1;
        }
      }

      ~Queue() { while (!empty()) pop(); }
    };

    /*
      This queue requires attention for protection against
      concurrent access. Protect against:
      - Concurrent access by two worker threads both
        wanting to get() a task from the queue at the same
        time.
      - Concurrent access by two threads both wanting to
        put() a task into the queue at the same time.
      - A worker thread having determined that the queue
        is empty, while at the same time a new task is put()
        into the queue.
      - A task wanting to put() a task into the queue
        having found the queue full, while at the same time
        the queues fill level decreases.
    */

    unsigned const          m_maxpart;
    bool                    m_shutting_down;
    unsigned                m_idle_workers;
    unsigned                m_total_workers;
    bool                    m_wakeup_is_pending;
    Queue                   m_queue;
    std::mutex              m_pop_mutex;
    std::mutex              m_push_mutex;
    std::condition_variable m_waiting_workers_cond;
    std::condition_variable m_waiters_cond;

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

      Queue functions(1);

      for (;;) {
        std::unique_lock<std::mutex> lock(m_pop_mutex);
        unsigned queue_size;

        // Try to get the next task(s)
        while ((queue_size = m_queue.size()) <= min_queue_size) {
          if (static_cast<int>(queue_size) <= return_if_idle) return;
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
        if (stride > functions.capacity()) functions.reserve(2 * stride);
        while (stride--) functions.push(m_queue.pop());
        lock.unlock();

        if ( m_idle_workers && !m_wakeup_is_pending && queue_size )
          m_waiting_workers_cond.notify_one();

        while (!functions.empty()) functions.pop()();
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
    work(bool return_if_idle) {
      help(return_if_idle ? 0 : -1);
    }

    /**
       Enqueue a task.
    */
    template<class C>
    void
    put(C&& c) {
      std::unique_lock<std::mutex> lock(m_push_mutex);
      while ( m_queue.size() >= m_queue.capacity() ) {
        // No space in the queue. Must wait for workers to advance.
        lock.unlock();
        try_help( m_queue.capacity() / 2 );
        lock.lock();
      }
      // Now there is space in the queue and we have locked the back.

      // Enqueue function.
      if (m_shutting_down) {
        Function fun(std::forward<C>(c)); // Run Function destructor
      } else {
        // Here we have exclusive access to the head of the queue.
        m_queue.push(std::forward<C>(c));

        if ( m_idle_workers && !m_wakeup_is_pending ) {
          m_wakeup_is_pending = true;
          m_waiting_workers_cond.notify_one();
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
      if (std::uncaught_exception())
      shutdown();
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

    unsigned queue_size() const { return m_queue.capacity(); }
    unsigned maxpart()    const { return m_maxpart; }

  };

} // End of namespace threadpool
