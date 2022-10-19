/**
 * \copyright 2021 Enrico Bertolazzi
 *
 * based on the work of Ruediger Helsch, Ruediger.Helsch@t-online.de (2014)
 * version 2.0 (https://github.com/RuedigerHelsch/ThreadPool)
 */

namespace threadpool {

  namespace parallel {

    /**
     * A thread pool to be used as base for container processing.
     *
     * \tparam InputIterator
     *         The iterator type to be used to get input values.
     *
     * \tparam OutputIterator
     *         The iterator type to be used to store the results.
     *
     * \tparam Function
     *         Type of the function to be called with
     *         successive elements from the input
     *         iterator. The function must return a result
     *         which is stored through the result iterator.
     */
    template<class InputIterator, class Last,class OutputIterator, class Function>
    class Transform_ThreadPool {
      using Queue = Transform_Queue<
        InputIterator,
        Last,
        OutputIterator,
        Function,
        is_forward_iterator<InputIterator>::value &&
        is_forward_iterator<OutputIterator>::value
      >;
      Queue                    m_queue;
      GenericThreadPool<Queue> m_pool;

    public:

      /**
       * Run a function on all objects in an input iterator range,
       * store the return values through an output iterator.
       *
       * \param first
       *        The start of the range to be processed.
       *
       * \param last
       *        One past the end of the range to be processed.
       *
       * \param result
       *        The iterator receiving the results.
       *
       * \param fun
       *        The function to apply to all objects in the container.
       *
       * \param thread_count
       *        The number of threads to spawn. If the
       *        default value of -1 is specified, the
       *        thread count is determined based on
       *        the number of available hardware
       *        threads. A value of 1 selects the
       *        single-threaded algorithm.
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
       *        might be four times the number of
       *        threads. A value of 0 enforces
       *        single-object processing.
       */
      Transform_ThreadPool(
        InputIterator  & first,
        Last const     & last,
        OutputIterator & result,
        Function       & fun,
        unsigned         thread_count,
        unsigned         maxpart
      )
      : m_queue( first, last, result, fun, maxpart )
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
      void
      join() { m_pool.join(); }
    };

    /**
     * Run a function on all objects in a range of iterators, store
     * the return values through an output iterator.
     *
     * \param first
     *        The start of the range to be processed.
     *
     * \param last
     *        One past the end of the range to be processed.
     *
     * \param result
     *        The iterator receiving the return values
     *
     * \param fun
     *        The function to apply to all objects in the range.
     *
     * \returns
     *        The final value of the result iterator
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
      class InputIterator,
      class Last,
      class OutputIterator,
      class Function,
      class = typename std::enable_if<std::is_same<InputIterator,Last>::value ||
              !std::is_integral<InputIterator>::value ||
              !std::is_integral<Last>::value
              >::type
    >
    typename std::decay<OutputIterator>::type
    transform(
      InputIterator  first,
      Last const &   last,
      OutputIterator result,
      Function&&     fun
    ) {
      if (thread_count <= 1) {
        return std::transform(first, last, result, fun);
      } else {
        Transform_ThreadPool<InputIterator, Last, OutputIterator, Function>(
          first, last, result, fun, thread_count, 3 * (thread_count + 1)
        );
        return result;
      }
    }

    /**
     * Run a function with each of a range of integral values, store
     * the return values through an output iterator.
     *
     * \param first
     *        The start of the range to be processed.
     *
     * \param last
     *        One past the end of the range to be processed.
     *
     * \param result
     *        The iterator receiving the return values
     *
     * \param fun
     *        The function to call with all numbers in the range.
     *
     * \returns
     *        The final value of the result iterator
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
      class InputIterator,
      class Last,
      class OutputIterator,
      class Function,
      class = typename std::enable_if<!std::is_same<InputIterator,Last>::value &&
              std::is_integral<InputIterator>::value &&
              std::is_integral<Last>::value
              >::type
    >
    typename std::decay<OutputIterator>::type
    transform(
      InputIterator  && first,
      Last const      & last,
      OutputIterator && result,
      Function&& fun
    ) {
      /*
        We can not use the generic function because the user
        might specify `first` as 0 which makes type
        `InputIterator' become `int`, and `last` as something of
        type `std::size_t` not representable in an `int`. This
        loop would run forever. Just extend type `InputIterator`.
       */
      using common_type    = typename std::common_type<InputIterator, Last>::type;
      using CommonIterator = IntegralIterator<common_type>;

      return transform<thread_count>(
        CommonIterator(std::forward<InputIterator>(first)),
        CommonIterator(last),
        std::forward<OutputIterator>(result),
        std::forward<Function>(fun)
      );
    }

    /**
     * Run a function on all objects in a container, store the return
     * values through an output iterator.
     *
     * Version for lvalue containers
     *
     * \param container
     *        The container.
     *
     * \param result
     *        The iterator receiving the return values
     *
     * \param fun
     *        The function to apply to all objects in the container.
     *
     * \returns
     *        The final value of the result iterator
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
      class OutputIterator,
      class Function
    >
    typename std::decay<OutputIterator>::type
    transform(
      Container      &  container,
      OutputIterator && result,
      Function       && fun
    ) {
      return transform<thread_count>(
        std::begin(container),
        std::end(container),
        std::forward<OutputIterator>(result),
        std::forward<Function>(fun)
      );
    }

    /**
     * Run a function on all objects in a container, store the return
     * values through an output iterator.
     *
     * Version for rvalue containers
     *
     * \param container
     *        The container.
     *
     * \param result
     *        The iterator receiving the return values
     *
     * \param fun
     *        The function to apply to all objects in the container.
     *
     * \returns
     *        The final value of the result iterator
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
      class OutputIterator,
      class Function,
      class = typename std::enable_if<!std::is_lvalue_reference<Container>::value>::type
    >
    typename std::decay<OutputIterator>::type
    transform(
      Container      && container,
      OutputIterator && result,
      Function       && fun
    ) {
      return transform<thread_count>(
        std::make_move_iterator(std::begin(container)),
        std::make_move_iterator(std::end(container)),
        std::forward<OutputIterator>(result),
        std::forward<Function>(fun)
      );
    }

  } // End of namespace parallel

} // End of namespace threadpool
