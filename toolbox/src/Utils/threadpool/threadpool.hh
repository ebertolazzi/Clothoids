/**
 * \file threadpool/threadpool.h
 *
 * Threadpool for C++11, header for thread pool
 *
 * \copyright 2021 Enrico Bertolazzi
 *
 * based on the work of Ruediger Helsch, Ruediger.Helsch@t-online.de (2014)
 * version 2.0 (https://github.com/RuedigerHelsch/ThreadPool)
 */

#ifndef THREADPOOL_dot_HH
#define THREADPOOL_dot_HH

#include <cstddef>
#include <memory>
#include <type_traits>  // For std::remove_reference()
#include <functional>		// For std::bind()
#include <cassert>
#include <iterator>
#include <utility>	    // For std::move(), std::forward()

#include <vector>
#include <limits>
#include <limits>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "threadpool_utils.hxx"
#include "threadpool_queue.hxx"

namespace threadpool {

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
    virtual ~VirtualTask() {}
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
    std::exception_ptr       m_pending_exception;
    Queue *                  m_queue;
    std::vector<std::thread> m_worker_threads;
    bool                     m_ignore_thread_pool_exceptions = true;

    //! The main function of the thread.
    void work() { help(false); }

    //! Wait for all workers to finish.
    void
    join_workers() {
      work(); // Instead of hanging around, help the workers!
      for ( std::thread & w : m_worker_threads )
        { if (w.joinable()) w.join(); }
    }

    // Copying and moving are not supported.
    GenericThreadPool( GenericThreadPool const & )             = delete;
    GenericThreadPool( GenericThreadPool && )                  = delete;
    GenericThreadPool& operator = (GenericThreadPool const & ) = delete;
    GenericThreadPool& operator = (GenericThreadPool && )      = delete;

	public:

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
    : m_pending_exception(nullptr)
    , m_queue(queue)
    , m_worker_threads(thread_count)
    {
      for ( std::thread & w : m_worker_threads )
        //w = std::move(std::thread(std::bind(&GenericThreadPool::work, this)));
        w = std::thread(std::bind(&GenericThreadPool::work, this));
    }

    /**
     * Help with the work.
     *
     * \param return_if_idle
     *        Never wait for work, return instead.
     */
    void
    help( bool return_if_idle ) {
      if ( m_ignore_thread_pool_exceptions ) {
        m_queue->work( return_if_idle );
      } else {
        try {
          m_queue->work( return_if_idle );
        } catch (...) {
          {
            std::exception_ptr e = std::current_exception();
            std::lock_guard<std::mutex> lock(m_mutex);
            if (!m_pending_exception) m_pending_exception = std::move(e);
          }
          m_queue->shutdown();
        }
      }
    }

    /**
     * Rethrow a potentially pending exception from a worker thread.
     */
    void
    rethrow_exception() {
      if ( m_pending_exception && !std::uncaught_exception() ) {
        m_queue->shutdown();
        join_workers();
        if (!std::uncaught_exception()) {
          std::exception_ptr e = m_pending_exception;
          m_pending_exception = nullptr;
          std::rethrow_exception(std::move(e));
        }
      }
    }

    /**
     * Wait for all threads to finish and collect them.
     *
     * Leaves the thread pool ready for destruction.
     */
    void
    join() {
      join_workers();
      rethrow_exception();
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

    /**
     * Switch exception handling on/off
     */
    void
    ignore_thread_pool_exceptions( bool flg = true )
    { m_ignore_thread_pool_exceptions = flg; }

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

    QueueElement()                                     = delete;
    QueueElement( QueueElement const & )               = delete;
    QueueElement & operator = ( QueueElement const & ) = delete;
    QueueElement & operator = ( QueueElement && )      = delete;

  public:

    QueueElement( VirtualTask * t ) : m_task(t) { }
    QueueElement( QueueElement && x ) : m_task(x.m_task) { x.m_task = nullptr; }
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
  class ThreadPool {

    typedef HQueue<QueueElement>     QUEUE;
    typedef GenericThreadPool<QUEUE> POOL;

    QUEUE * m_queue;
    POOL  * m_pool;

  public:

    explicit
    ThreadPool(
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

    /**
     * Wrap void functions in a task and run them without
     * exception handling.
     */
    template<class Function>
    void
    run( Function && f ) {
      typedef typename std::remove_reference<Function>::type function_type;
      class WrappedFunction : public VirtualTask {
        function_type m_f;
      public:
        WrappedFunction( function_type && f ) : m_f(std::move(f)) { }
        virtual void operator()() override { m_f(); delete this; }
      };
      run_task(new WrappedFunction(std::forward<Function>(f)));
    }

    /**
     * Wrap functions with arguments in a task and run them without
     * exception handling.
     */
    template <class Function, class... Args>
    void
    run( Function && f, Args && ... args ) {
      run( std::bind(std::forward<Function>(f),std::forward<Args>(args)...) );
    }

    /**
     * Run a function on all objects in an iterator range
     *
     * \param first
     *        Start of the range
     *
     * \param last
     *        End of the range
     *
     * \param fun
     *        The function taking one parameter
     *        by reference and returning void.
     *
     * Does not wait for all tasks to finish! Caller is
     * responsible for wait()ing on the pool if necessary.
     */
    template <class Iterator, class Function>
    void
    for_each( Iterator first, const Iterator& last, Function&& fun ) {
      while (first != last) {
        typedef iterval_traits<Iterator> IT;
        Wrap<typename IT::type> e(IT::copy(first));
        ++first;
        run([&fun,e](){ fun(IT::pass(std::move(e.value))); });
      }
    }

    /**
     * Run a function on all members of a container
     *
     * \param container
     *        The container to process
     *
     * \param fun
     *        The function taking one parameter
     *        by reference and returning void.
     *
     * Does not wait for all tasks to finish! Caller is
     * responsible for wait()ing on the pool if necessary.
     */
    template<class Container, class Function>
    void
    for_each( Container&& container, Function&& fun ) {
      for ( auto & e: container ) run([&fun,&e](){ fun(e); });
    }

    /**
     * Wait for all active tasks to finish.
     *
     * Also throws an exception if one of the tasks has
     * encountered an uncatched exception.
     *
     * Leaves the pool in a valid state ready to run more
     * tasks, unless an exception has been thrown.
     */
    void
    wait() {
      m_pool->help(true); // Help out instead of sitting around idly.
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
    join() {
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

    unsigned thread_count() const { return m_pool->thread_count(); }
    unsigned queue_size()   const { return m_queue->queue_size(); }
    unsigned maxpart()      const { return m_queue->maxpart(); }

    void
    resize(
      unsigned thread_count = std::thread::hardware_concurrency(),
      unsigned queue_size   = 0,
      unsigned maxpart      = 0
    ) {
      this->~ThreadPool();
      if ( queue_size == 0 ) queue_size = 50 * (thread_count+1);
      if ( maxpart    == 0 ) maxpart    = 3 * (thread_count+1);
      m_queue = new QUEUE( queue_size, maxpart );
      m_pool  = new POOL( m_queue, thread_count );
    }

    ~ThreadPool() { wait(); join(); delete m_pool; delete m_queue; }
  };

} // End of namespace threadpool

#include "threadpool_parallel_for_each.hxx"
#include "threadpool_parallel_transform.hxx"
#include "threadpool_make_iterator.hxx"

#endif // !defined(THREADPOOL_THREADPOOL_H)
