/**
 * \copyright 2021 Enrico Bertolazzi
 *
 * based on the work of Ruediger Helsch, Ruediger.Helsch@t-online.de (2014)
 * version 2.0 (https://github.com/RuedigerHelsch/ThreadPool)
 */

namespace threadpool {

  namespace parallel {

    /*
      Now that the queue implementation is done, the definition of a
      thread pool for container processing is easy.
    */

    /**
     * A thread pool to be used as base for container processing.
     *
     * \tparam Iterator
     *	       The iterator type to be used to traverse the container.
     *
     * \tparam Last
     *         The iterator type for the last element.
     *
     * \tparam Function
     *         The function type to execute. Must be callable
     *         with a reference to a container value_type as
     *         parameter, e.g.  as void(Element& e).
     *
     */
    template<class Iterator, class Last, class Function>
    class ForEach_ThreadPool {
      typedef ForEach_Queue<
        Iterator,
        Last,
        Function,
        is_forward_iterator<Iterator>::value
      > Queue;
      Queue                    m_queue;
      GenericThreadPool<Queue> m_pool;

    public:

      /**
       * Run a function on all objects in an iterator range.
       *
       * \param first
       *        The start of the range to be processed.
       *
       * \param last
       *        One past the end of the range to be processed.
       *
       * \param fun
       *        The function to apply to all objects
       *        in the container.
       *
       * \param thread_count
       *        The number of threads to use. If the
       *        thread count is specified as -1 it
       *        defaults to the number of available
       *        hardware threads
       *        std::thread::hardware_concurrency().
       *
       * \param maxpart
       *        The maximum part of the remaining
       *        input range that one thread is allowed
       *        to take.  If maxpart is for example 5
       *        and 100 elements remain to be
       *        processed, then a task will take 100 /
       *        5 = 20 elements and process them. If a
       *        large value is chosen for maxpart,
       *        each thread will take small chunks of
       *        work and will look for more work
       *        frequently, causing increased
       *        synchronization overhead. If a small
       *        value is chosen for maxpart, each
       *        thread will take huge chunks of work,
       *        possibly leaving the remaining threads
       *        out of work at the end. A good value
       *        might be three times the number of
       *        threads. A value of 0 enforces
       *        single-object processing.
       */
      ForEach_ThreadPool(
        Iterator   & first,
        Last const & last,
        Function   & fun,
        unsigned     thread_count,
        unsigned     maxpart
      )
      : m_queue( first, last, fun, maxpart )
      , m_pool( &m_queue, thread_count )
      { }

      /**
       * Collect threads, throw any pending exceptions.
       *
       * Use this function if waiting in the desctructor or throwing
       * exceptions from the destructor is undesirable.  After
       * join() returned, the thread pool can be destroyed without
       * delay and without throwing.
       */
       void join() { m_pool.join(); }

       unsigned thread_count() const { return m_pool.thread_count(); }

    };

    /**
     * Run a function on all objects in a range of iterators.
     *
     * \param first
     *        The start of the range to be processed.
     *
     * \param last
     *        One past the end of the range to be processed.
     *
     * \param fun
     *        The function to apply to all objects in the range.
     *
     * \returns
     *        The final value of the function
     *
     * \tparam thread_count
     *         The number of threads to spawn. If the default
     *         value of -1 is specified, the thread count is
     *         determined based on the number of available
     *         hardware threads. A value of 1 selects the
     *         single-threaded algorithm.
     *
     * \tparam maxpart The maximum part of the remaining input
     *         range that one thread is allowed to take.  If
     *         maxpart is for example 5 and 100 elements
     *         remain to be processed, then a task will take
     *         100 / 5 = 20 elements and process them. If a
     *         large value is chosen for maxpart, each thread
     *         will take small chunks of work and will look
     *         for more work frequently, causing increased
     *         synchronization overhead. If a small value is
     *         chosen for maxpart, each thread will take huge
     *         chunks of work, possibly leaving the remaining
     *         threads out of work at the end. A good value
     *         might be three times the number of
     *         threads. The default value of 1 lets the
     *         system determine a value, which is three times
     *         the number of threads. A value of 0 enforces
     *         single-object processing.
     */
    template<
      int   thread_count,
      class Iterator,
      class Last,
      class Function,
      class = typename std::enable_if<std::is_same<Iterator,Last>::value ||
              !std::is_integral<Iterator>::value ||
              !std::is_integral<Last>::value
              >::type
    >
    typename std::decay<Function>::type
    for_each(
      Iterator     first,
      Last const & last,
      Function  && fun
    ) {
      if (thread_count <= 1) {
        return std::for_each(first, last, fun);
      } else {
        ForEach_ThreadPool<Iterator, Last, Function>(
          first, last, fun, thread_count, 3 * (thread_count + 1)
        );
        return std::forward<Function>(fun);
      }
    }

    /**
     * Run a function with each of a range of integral values.
     *
     * \param first
     *        The start of the range to be processed.
     *
     * \param last
     *        One past the end of the range to be processed.
     *
     * \param fun
     *        The function to call with all numbers in the range.
     *
     * \returns
     *        The final value of the function
     *
     * \tparam thread_count
     *         The number of threads to spawn. If the default
     *         value of -1 is specified, the thread count is
     *         determined based on the number of available
     *         hardware threads. A value of 1 selects the
     *         single-threaded algorithm.
     *
     * \tparam maxpart The maximum part of the remaining input
     *         range that one thread is allowed to take.  If
     *         maxpart is for example 5 and 100 elements
     *         remain to be processed, then a task will take
     *         100 / 5 = 20 elements and process them. If a
     *         large value is chosen for maxpart, each thread
     *         will take small chunks of work and will look
     *         for more work frequently, causing increased
     *         synchronization overhead. If a small value is
     *         chosen for maxpart, each thread will take huge
     *         chunks of work, possibly leaving the remaining
     *         threads out of work at the end. A good value
     *         might be three times the number of
     *         threads. The default value of 1 lets the
     *         system determine a value, which is three times
     *         the number of threads. A value of 0 enforces
     *         single-object processing.
     */
    template<
      int   thread_count,
      class Iterator,
      class Last,
      class Function,
      class = typename std::enable_if<!std::is_same<Iterator,Last>::value &&
              std::is_integral<Iterator>::value &&
              std::is_integral<Last>::value
              >::type
    >
    typename std::decay<Function>::type
    for_each(
      Iterator  && first,
      Last const & last,
      Function  && fun
    )
    {
      /*
        We can not use the generic function because the user
        might specify `first` as 0 which makes type `Iterator'
        become `int`, and `last` as something of type
        `std::size_t` not representable in an `int`. This loop
        would run forever. Just extend type `Iterator`.
       */
      typedef typename std::common_type<Iterator, Last>::type common_type;
      typedef IntegralIterator<common_type> CommonIterator;

      return for_each<thread_count>(
        CommonIterator(std::forward<Iterator>(first)),
        CommonIterator(last),
        std::forward<Function>(fun)
      );
    }

    /**
     * Run a function on all objects in a container.
     *
     * Version for lvalue containers. The objects in the container
     * are passed to `fun` by reference.
     *
     * \param container
     *        The container.
     *
     * \param fun
     *        The function to apply to all objects in the
     *        container.
     *
     * \returns
     *        The final value of the function
     *
     * \tparam thread_count
     *         The number of threads to spawn. If the default
     *         value of -1 is specified, the thread count is
     *         determined based on the number of available
     *         hardware threads. A value of 1 selects the
     *         single-threaded algorithm.
     *
     * \tparam maxpart The maximum part of the remaining input
     *         range that one thread is allowed to take.  If
     *         maxpart is for example 5 and 100 elements
     *         remain to be processed, then a task will take
     *         100 / 5 = 20 elements and process them. If a
     *         large value is chosen for maxpart, each thread
     *         will take small chunks of work and will look
     *         for more work frequently, causing increased
     *         synchronization overhead. If a small value is
     *         chosen for maxpart, each thread will take huge
     *         chunks of work, possibly leaving the remaining
     *         threads out of work at the end. A good value
     *         might be three times the number of
     *         threads. The default value of 1 lets the
     *         system determine a value, which is three times
     *         the number of threads. A value of 0 enforces
     *         single-object processing.
     */
    template<
      int   thread_count,
      class Container,
      class Function
    >
    typename std::decay<Function>::type
    for_each( Container & container, Function&& fun ) {
      return for_each<thread_count>(
        std::begin(container),
        std::end(container),
        std::forward<Function>(fun)
      );
    }

    /**
     * Run a function on all objects in a container.
     *
     * Version for rvalue containers. The objects in the container
     * are passed to `fun` by rvalue reference, so they can be
     * move()d.
     *
     * \param container
     *        The container.
     *
     * \param fun
     *        The function to apply to all objects in the
     *        container.
     *
     * \returns
     *        The final value of the function
     *
     * \tparam thread_count
     *         The number of threads to spawn. If the default
     *         value of -1 is specified, the thread count is
     *         determined based on the number of available
     *         hardware threads. A value of 1 selects the
     *         single-threaded algorithm.
     *
     * \tparam maxpart The maximum part of the remaining input
     *         range that one thread is allowed to take.  If
     *         maxpart is for example 5 and 100 elements
     *         remain to be processed, then a task will take
     *         100 / 5 = 20 elements and process them. If a
     *         large value is chosen for maxpart, each thread
     *         will take small chunks of work and will look
     *         for more work frequently, causing increased
     *         synchronization overhead. If a small value is
     *         chosen for maxpart, each thread will take huge
     *         chunks of work, possibly leaving the remaining
     *         threads out of work at the end. A good value
     *         might be three times the number of
     *         threads. The default value of 1 lets the
     *         system determine a value, which is three times
     *         the number of threads. A value of 0 enforces
     *         single-object processing.
     */
    template<
      int thread_count,
      class Container,
      class Function,
      class = typename std::enable_if<!std::is_lvalue_reference<Container>::value>::type
    >
    typename std::decay<Function>::type
    for_each( Container&& container, Function&& fun ) {
      return for_each<thread_count>(
        std::make_move_iterator(std::begin(container)),
        std::make_move_iterator(std::end(container)),
        std::forward<Function>(fun));
    }
  } // End of namespace parallel
} // End of namespace threadpool
