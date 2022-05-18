/*--------------------------------------------------------------------------*\
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/
/*
  Original code by Thomas Nagler




 */




// Copyright 2021 Thomas Nagler (MIT License)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include <exception>
#include <memory>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#ifdef UTILS_QP_USE_FUTURE
#include <future>
#endif

#if (defined __linux__ || defined AFFINITY)
#include <pthread.h>
#endif

// Layout of quickpool.hpp
//
// 1. Memory related utilities.
//    - Memory order aliases
//    - Class for padding bytes
//    - Memory aligned allocation utilities
//    - Class for cache aligned atomics
//    - Class for load/assign atomics with relaxed order
// 2. Loop related utilities.
//    - Worker class for parallel for loops
// 3. Scheduling utilities.
//    - Ring buffer
//    - Task queue
//    - Task manager
// 4. Thread pool class
// 5. Free-standing functions (main API)

//! quickpool namespace
namespace quickpool {

  // 1. --------------------------------------------------------------------------

  //! Memory related utilities.
  namespace mem {

    //! convenience definitions
    static constexpr std::memory_order relaxed = std::memory_order_relaxed;
    static constexpr std::memory_order acquire = std::memory_order_acquire;
    static constexpr std::memory_order release = std::memory_order_release;
    static constexpr std::memory_order seq_cst = std::memory_order_seq_cst;

    //! Padding char[]s always must hold at least one char. If the size of the
    //! object ends at an alignment point, we don't want to pad one extra byte
    //! however. The construct below ensures that padding bytes are only added if
    //! necessary.
    namespace padding_impl {

      //! constexpr modulo operator.
      constexpr size_t mod( size_t a, size_t b ) { return a - b * (a / b); }

      // Padding bytes from end of aligned object until next alignment point. char[]
      // must hold at least one byte.
      template < typename T, size_t Align >
      struct padding_bytes {
        static constexpr size_t free_space = Align - mod(sizeof(std::atomic<T>), Align);
        static constexpr size_t required   = free_space > 1 ? free_space : 1;
        char padding_[required];
      };

      struct empty_struct{};

      //! Class holding padding bytes is necessary. Classes can inherit from this
      //! to automically add padding if necessary.
      template < typename T, size_t Align >
      struct padding : std::conditional<
                         mod(sizeof(std::atomic<T>), Align) != 0,
                         padding_bytes<T, Align>,
                         empty_struct
                       >::type
      {};

    } // end namespace padding_impl

    namespace aligned {

      // The alloc/dealloc mechanism is pretty much
      // https://www.boost.org/doc/libs/1_76_0/boost/align/detail/aligned_alloc.hpp

      inline
      void*
      alloc( size_t alignment, size_t size ) noexcept
      { return aligned_alloc( alignment, size ); }

      inline
      void
      free( void * & ptr ) noexcept
      { if ( ptr ) std::free( ptr ); ptr = nullptr; }

      //! short version of
      //! https://www.boost.org/doc/libs/1_65_0/boost/align/aligned_allocator.hpp
      template < typename T, std::size_t Alignment = 64 >
      class allocator : public std::allocator<T> {
      private:
        static constexpr size_t min_align =
          (Alignment >= alignof(void*)) ? Alignment : alignof(void*);
      public:
        template < typename U >
        struct rebind {
          typedef allocator<U, Alignment> other;
        };

        T*
        allocate( size_t size, const void * = nullptr ) {
          if ( size == 0 ) return 0;
          void * p = mem::aligned::alloc(min_align, sizeof(T) * size);
          if (!p) throw std::bad_alloc();
          return static_cast<T*>(p);
        }

        void deallocate( T* ptr, size_t ) { mem::aligned::free(ptr); }

        template < typename U, class... Args>
        void
        construct( U* ptr, Args&&... args )
        { ::new ((void*)ptr) U(std::forward<Args>(args)...); }

        template < typename U >
        void
        destroy( U * ptr ) {
          (void)ptr;
          ptr->~U();
        }
      };

      //! Memory-aligned atomic `std::atomic<T>`. Behaves like `std::atomic<T>`, but
      //! overloads operators `new` and `delete` to align its memory location. Padding
      //! bytes are added if necessary.
      template < typename T, size_t Align = 64 >
      struct alignas(Align) atomic : public std::atomic<T>,
                                     private padding_impl::padding<T, Align> {
      public:
        atomic() noexcept = default;
        atomic(T desired) noexcept : std::atomic<T>(desired) {}

        // Assignment operators have been deleted, must redefine.
        T operator=(T x) noexcept { return std::atomic<T>::operator=(x); }
        T operator=(T x) volatile noexcept { return std::atomic<T>::operator=(x); }

        static
        void*
        operator new(size_t count) noexcept
        { return mem::aligned::alloc(Align, count); }

        static
        void
        operator delete(void* ptr)
        { mem::aligned::free(ptr); }
      };

      //! Fast and simple load/assign atomic with no memory ordering guarantees.
      template <typename T>
      struct relaxed_atomic : public mem::aligned::atomic<T> {
        explicit relaxed_atomic(T value) : mem::aligned::atomic<T>(value) {}
        operator T() const noexcept { return this->load(mem::relaxed); }

        T operator=(T desired) noexcept {
          this->store(desired, mem::relaxed);
          return desired;
        }
      };

      //! vector class for aligned types.
      template <typename T, size_t Alignment = 64>
      using vector = std::vector<T, mem::aligned::allocator<T, Alignment>>;

    } // end namespace aligned

  } // end namespace mem

  // 2. --------------------------------------------------------------------------

  //! Loop related utilities.
  namespace loop {

    //! Worker state.
    struct State {
      int pos; //!< position in the loop range
      int end; //!< end of range assigned to worker
    };

    //! Worker class for parallel loops.
    //!
    //! When a worker completes its own range, it steals half of the remaining range
    //! of another worker. The number of steals (= only source of contention) is
    //! therefore only logarithmic in the number of tasks. The algorithm uses
    //! double-width compare-and-swap, which is lock-free on most modern processor
    //! architectures.
    //!
    //! @tparam type of function processing the loop (required to account for
    //! type-erasure).
    template <typename Function>
    struct Worker {
      Worker() {}

      Worker(int begin, int end, Function fun)
      : state{ State{ begin, end } }, f{ fun }
      {}

      Worker(Worker&& other)
      : state{ other.state.load() }
      , f{ std::forward<Function>(other.f) }
      {}

      size_t
      tasks_left() const {
        State s = state.load();
        return s.end - s.pos;
      }

      bool done() const { return (tasks_left() == 0); }

      //! @param others pointer to the vector of all workers.
      void
      run( std::shared_ptr<mem::aligned::vector<Worker>> others ) {
        State s, s_old; // temporary state variables
        do {
          s = state.load();
          if ( s.pos < s.end ) {
            // Protect slot by trying to advance position before doing work.
            s_old = s;
            s.pos++;

            // Another worker might have changed the end of the range in
            // the meanwhile. Check atomically if the state is unaltered
            // and, if so, replace by advanced state.
            if ( state.compare_exchange_weak(s_old, s) ) {
              f(s_old.pos); // succeeded, do work
            } else {
              continue; // failed, try again
            }
          }
          if ( s.pos == s.end ) {
            // Reached end of own range, steal range from others. Range
            // remains empty if all work is done, so we can leave the loop.
            this->steal_range(*others);
          }
        } while (!this->done());
      }

      //! @param workers vector of all workers.
      void
      steal_range( mem::aligned::vector<Worker> & workers ) {
        do {
          Worker& other = find_victim(workers);
          State s = other.state.load();
          if ( s.pos >= s.end ) continue; // other range is empty by now

          // Remove second half of the range. Check atomically if the
          // state is unaltered and, if so, replace with reduced range.
          auto s_old = s;
          s.end -= (s.end - s.pos + 1) / 2;
          if ( other.state.compare_exchange_weak(s_old, s) ) {
            // succeeded, update own range
            state = State{ s.end, s_old.end };
            break;
          }
        } while (!all_done(workers)); // failed steal, try again
      }

      //! @param workers vector of all workers.
      bool
      all_done( mem::aligned::vector<Worker> const & workers ) {
        for ( auto const & worker : workers ) {
          if ( !worker.done() ) return false;
        }
        return true;
      }

      //! targets the worker with the largest remaining range to minimize
      //! number of steal events.
      //! @param others vector of all workers.
      Worker &
      find_victim( mem::aligned::vector<Worker> & workers ) {
        std::vector<size_t> tasks_left;
        tasks_left.reserve(workers.size());
        for ( auto const & worker : workers ) {
          tasks_left.push_back(worker.tasks_left());
        }
        auto max_it = std::max_element(tasks_left.begin(), tasks_left.end());
        auto idx    = std::distance(tasks_left.begin(), max_it);
        return workers[idx];
      }

      mem::aligned::relaxed_atomic<State> state; //!< worker state `{pos, end}`
      Function                            f;     //< function applied to the loop index
    };

    //! creates loop workers. They must be passed to each worker using a shared
    //! pointer, so that they persist if an inner `parallel_for()` in a nested
    //! loop exits.
    template <typename Function>
    std::shared_ptr<mem::aligned::vector<Worker<Function>>>
    create_workers( Function const & f, int begin, int end, size_t num_workers ) {
      auto num_tasks = std::max(end - begin, static_cast<int>(0));
      num_workers = std::max(num_workers, static_cast<size_t>(1));
      auto workers = new mem::aligned::vector<Worker<Function>>;
      workers->reserve(num_workers);
      for ( size_t i = 0; i < num_workers; ++i ) {
        workers->emplace_back(
          begin + num_tasks * i / num_workers,
          begin + num_tasks * (i + 1) / num_workers,
          f
        );
      }
      return std::shared_ptr<mem::aligned::vector<Worker<Function>>>(std::move(workers));
    }
  } // end namespace loop

  // 3. -------------------------------------------------------------------------

  //! Task management utilities.
  namespace sched {

    //! A simple ring buffer class.
    template <typename T>
    class RingBuffer {
    public:
      explicit RingBuffer(size_t capacity)
      : buffer_{ std::unique_ptr<T[]>(new T[capacity]) }
      , capacity_{ capacity }
      , mask_{ capacity - 1 }
      {}

      size_t capacity() const { return capacity_; }

      void set_entry( size_t i, T val ) { buffer_[i & mask_] = val; }

      T get_entry(size_t i) const { return buffer_[i & mask_]; }

      RingBuffer<T>*
      enlarged_copy(size_t bottom, size_t top) const {
        RingBuffer<T>* new_buffer = new RingBuffer{ 2 * capacity_ };
        for (size_t i = top; i != bottom; ++i)
          new_buffer->set_entry(i, this->get_entry(i));
        return new_buffer;
      }

    private:
      std::unique_ptr<T[]> buffer_;
      size_t               capacity_;
      size_t               mask_;
    };

    //! A multi-producer, multi-consumer queue; pops are lock free.
    class TaskQueue {
      using Task = std::function<void()>;

    public:
      //! @param capacity must be a power of two.
      TaskQueue(size_t capacity = 256)
      : buffer_{ new RingBuffer<Task*>(capacity) }
      {}

      ~TaskQueue() noexcept {
        // Must free memory allocated by push(), but not freed by try_pop().
        auto buf_ptr = buffer_.load();
        for (int i = top_; i < bottom_.load(mem::relaxed); ++i)
          delete buf_ptr->get_entry(i);
        delete buf_ptr;
      }

      TaskQueue(TaskQueue const& other) = delete;
      TaskQueue& operator=(TaskQueue const& other) = delete;

      //! checks if queue is empty.
      bool
      empty() const {
        return (bottom_.load(mem::relaxed) <= top_.load(mem::relaxed));
      }

      //! pushes a task to the bottom of the queue; returns false if queue is
      //! currently locked; enlarges the queue if full.
      void
      push( Task&& task ) {
        // Must hold lock in case of multiple producers.
        std::unique_lock<std::mutex> lk(mutex_);
        auto b = bottom_.load(mem::relaxed);
        auto t = top_.load(mem::acquire);
        RingBuffer<Task*>* buf_ptr = buffer_.load(mem::relaxed);

        if ( static_cast<int>(buf_ptr->capacity()) < (b - t) + 1 ) {
          // Buffer is full, create enlarged copy before continuing.
          auto old_buf = buf_ptr;
          buf_ptr = std::move(buf_ptr->enlarged_copy(b, t));
          old_buffers_.emplace_back(old_buf);
          buffer_.store(buf_ptr, mem::relaxed);
        }

        //! Store pointer to new task in ring buffer.
        buf_ptr->set_entry(b, new Task{ std::forward<Task>(task) });
        bottom_.store(b + 1, mem::release);

        lk.unlock(); // can release before signal
        cv_.notify_one();
      }

      //! pops a task from the top of the queue; returns false if lost race.
      bool
      try_pop( Task & task ) {
        auto t = top_.load(mem::acquire);
        std::atomic_thread_fence(mem::seq_cst);
        auto b = bottom_.load(mem::acquire);

        if ( t < b ) {
          // Must load task pointer before acquiring the slot, because it
          // could be overwritten immediately after.
          auto task_ptr = buffer_.load(mem::acquire)->get_entry(t);

          // Atomically try to advance top.
          if ( top_.compare_exchange_strong( t, t + 1, mem::seq_cst, mem::relaxed) ) {
            task = std::move(*task_ptr); // won race, get task
            delete task_ptr;             // fre memory allocated in push()
            return true;
          }
        }
        return false; // queue is empty or lost race
      }

      //! waits for tasks or stop signal.
      void
      wait() {
        std::unique_lock<std::mutex> lk(mutex_);
        cv_.wait(lk, [this] { return !this->empty() || stopped_; });
      }

      //! stops the queue and wakes up all workers waiting for jobs.
      void
      stop() {
        {
          std::lock_guard<std::mutex> lk(mutex_);
          stopped_ = true;
        }
        cv_.notify_one();
      }

      void
      wake_up() {
        {
          std::lock_guard<std::mutex> lk(mutex_);
        }
       cv_.notify_one();
      }

    private:

      //! queue indices
      mem::aligned::atomic<int> top_{ 0 };
      mem::aligned::atomic<int> bottom_{ 0 };

      //! ring buffer holding task pointers
      std::atomic<RingBuffer<Task*>*> buffer_{ nullptr };

      //! pointers to buffers that were replaced by enlarged buffer
      std::vector<std::unique_ptr<RingBuffer<Task*>>> old_buffers_;

      //! synchronization variables
      std::mutex mutex_;
      std::condition_variable cv_;
      bool stopped_{ false };
    };

    //! Task manager based on work stealing.
    class TaskManager {
    public:

      explicit
      TaskManager(size_t num_queues)
      : queues_(num_queues)
      , num_queues_(num_queues)
      , owner_id_(std::this_thread::get_id())
      {}

      TaskManager &
      operator = ( TaskManager&& other ) {
        std::swap(queues_, other.queues_);
        num_queues_  = other.num_queues_;
        status_      = other.status_.load();
        num_waiting_ = other.num_waiting_.load();
        push_idx_    = other.push_idx_.load();
        todo_        = other.todo_.load();
        return *this;
      }

      void
      resize( size_t num_queues ) {
        num_queues_ = std::max(num_queues, static_cast<size_t>(1));
        if ( num_queues > queues_.size() ) {
          queues_ = mem::aligned::vector<TaskQueue>(num_queues);
          // thread pool must have stopped the manager, reset
          num_waiting_ = 0;
          todo_        = 0;
          status_      = Status::running;
        }
      }

      template <typename Task>
      void
      push(Task&& task) {
        rethrow_exception(); // push() throws if a task has errored.
        if (is_running()) {
          todo_.fetch_add(1, mem::release);
          queues_[push_idx_++ % num_queues_].push(task);
        }
      }

      template <typename Task>
      bool
      try_pop( Task& task, size_t worker_id = 0 ) {
        // Always start pop cycle at own queue to avoid contention.
        for ( size_t k = 0; k <= num_queues_; ++k ) {
          if ( queues_[(worker_id + k) % num_queues_].try_pop(task) ) {
            return is_running(); // Throw away task if pool has stopped or errored.
          }
        }
        return false;
      }

      void
      wake_up_all_workers() {
        for ( auto & q : queues_ ) q.wake_up();
      }

      void
      wait_for_jobs( size_t id ) {
        if ( has_errored() ) {
          // Main thread may be waiting to reset the pool.
          std::lock_guard<std::mutex> lk(mtx_);
          if ( ++num_waiting_ == queues_.size() ) cv_.notify_all();
        } else {
          ++num_waiting_;
        }
        queues_[id].wait();
        --num_waiting_;
      }

      //! @param millis if > 0: stops waiting after millis ms
      void
      wait_for_finish( size_t millis = 0 ) {
        if ( called_from_owner_thread() && is_running() ) {
          auto wake_up = [this] { return (todo_ <= 0) || !is_running(); };
          std::unique_lock<std::mutex> lk(mtx_);
          if ( millis == 0 ) {
            cv_.wait(lk, wake_up);
          } else {
            cv_.wait_for(lk, std::chrono::milliseconds(millis), wake_up);
          }
        }
        rethrow_exception();
      }

      bool
      called_from_owner_thread() const {
        return (std::this_thread::get_id() == owner_id_);
      }

      void
      report_success() {
        auto n = todo_.fetch_sub(1, mem::release) - 1;
        if (n == 0) {
          // all jobs are done; lock before signal to prevent spurious failure
          {
            std::lock_guard<std::mutex> lk{ mtx_ };
          }
          cv_.notify_all();
        }
      }

      void
      report_fail( std::exception_ptr err_ptr ) {
        std::lock_guard<std::mutex> lk(mtx_);
        if ( has_errored() ) return; // only catch first exception

        err_ptr_ = err_ptr;
        status_  = Status::errored;

        // Some threads may change todo_ after we stop. The large
        // negative number forces them to exit the processing loop.
        todo_.store(std::numeric_limits<int>::min() / 2);
        cv_.notify_all();
      }

      void
      stop() {
        {
          std::lock_guard<std::mutex> lk(mtx_);
          status_ = Status::stopped;
        }
        // Worker threads wait on queue-specific mutex -> notify all queues.
        for ( auto & q : queues_ ) q.stop();
      }

      void
      rethrow_exception() {
        // Exceptions are only thrown from the owner thread, not in workers.
        if ( called_from_owner_thread() && has_errored() ) {
          {
            // Wait for all threads to idle so we can clean up after them.
            std::unique_lock<std::mutex> lk(mtx_);
            cv_.wait(lk, [this] { return num_waiting_ == queues_.size(); });
          }
          // Before throwing: restore defaults for potential future use of
          // the task manager.
          todo_ = 0;
          auto current_exception = err_ptr_;
          err_ptr_ = nullptr;
          status_  = Status::running;

          std::rethrow_exception(current_exception);
        }
      }

      bool
      is_running() const {
        return status_.load(mem::relaxed) == Status::running;
      }

      bool
      has_errored() const {
        return status_.load(mem::relaxed) == Status::errored;
      }

      bool
      stopped() const {
        return status_.load(mem::relaxed) == Status::stopped;
      }

      bool done() const { return (todo_.load(mem::relaxed) <= 0); }

    private:

      //! worker queues
      mem::aligned::vector<TaskQueue> queues_;
      size_t num_queues_;

      //! task management
      mem::aligned::relaxed_atomic<size_t> num_waiting_{ 0 };
      mem::aligned::relaxed_atomic<size_t> push_idx_{ 0 };
      mem::aligned::atomic<int> todo_{ 0 };

      //! synchronization variables
      std::thread::id const owner_id_;
      enum class Status { running, errored, stopped };
      mem::aligned::atomic<Status> status_{ Status::running };
      std::mutex mtx_;
      std::condition_variable cv_;
      std::exception_ptr err_ptr_{ nullptr };
    };

  } // end namespace sched

  // 4. ------------------------------------------------------------------------

  //! A work stealing thread pool.
  class ThreadPool {
  public:
    //! @brief constructs a thread pool.
    //! @param threads number of worker threads to create; defaults to
    //! number of available (virtual) hardware cores.
    explicit
    ThreadPool( std::size_t threads = std::thread::hardware_concurrency() )
    : task_manager_{ threads }
    {
      set_active_threads(threads);
    }

    ~ThreadPool() {
      task_manager_.stop();
      join_threads();
    }

    ThreadPool( ThreadPool&& ) = delete;
    ThreadPool( ThreadPool const & ) = delete;
    ThreadPool& operator=( ThreadPool const &) = delete;
    ThreadPool& operator=( ThreadPool && other ) = delete;

    //! @brief returns a reference to the global thread pool instance.
    static
    ThreadPool & global_instance() {
#ifdef _WIN32
      // Must leak resource, because windows + R deadlock otherwise.
      // Memory is released on shutdown.
      static auto ptr = new ThreadPool;
      return *ptr;
#else
      static ThreadPool instance_;
      return instance_;
#endif
    }

    //! @brief sets the number of active worker threads in the thread pool.
    //! @param threads the number of worker threads.
    //! Has no effect when not called from owner thread.
    void
    set_active_threads( size_t threads ) {
      if ( !task_manager_.called_from_owner_thread() ) return;

      if ( threads <= workers_.size() ) {
        task_manager_.resize(threads);
      } else {
        if ( workers_.size() > 0 ) {
          task_manager_.stop();
          join_threads();
        }
        workers_      = std::vector<std::thread>{ threads };
        task_manager_ = quickpool::sched::TaskManager{ threads };
        for ( size_t id = 0; id < threads; ++id) add_worker(id);
      }
      active_threads_ = threads;
    }

    //! @brief retrieves the number of active worker threads in the thread pool.
    size_t get_active_threads() const { return active_threads_; }

    //! @brief pushes a job to the thread pool.
    //! @param f a function.
    //! @param args (optional) arguments passed to `f`.
    template <typename Function, class... Args>
    void
    push( Function&& f, Args&&... args ) {
      if (active_threads_ == 0) return f(args...);
      task_manager_.push(std::bind(std::forward<Function>(f),std::forward<Args>(args)...));
    }

    #ifdef UTILS_QP_USE_FUTURE
    //! @brief executes a job asynchronously on the global thread pool.
    //! @param f a function.
    //! @param args (optional) arguments passed to `f`.
    //! @return A `std::future` for the task. Call `future.get()` to
    //! retrieve the results at a later point in time (blocking).
    template <typename Function, class... Args>
    auto
    async(Function&& f, Args&&... args) -> std::future<decltype(f(args...))> {
      auto  pack     = std::bind(std::forward<Function>(f), std::forward<Args>(args)...);
      using pack_t   = std::packaged_task<decltype(f(args...))()>;
      auto  task_ptr = std::make_shared<pack_t>(std::move(pack));
      this->push([task_ptr] { (*task_ptr)(); });
      return task_ptr->get_future();
    }
    #endif

    //! @brief computes an index-based parallel for loop.
    //!
    //! Waits until all tasks have finished, unless called from a thread
    //! that didn't create the pool. If this is taken into account, parallel
    //! loops can be nested.
    //!
    //! @param begin first index of the loop.
    //! @param end the loop runs in the range `[begin, end)`.
    //! @param f a function taking `int` argument (the 'loop body').
    template <typename UnaryFunction>
    void
    parallel_for( int begin, int end, UnaryFunction f ) {
      // each worker has its dedicated range, but can steal part of
      // another worker's ranges when done with own
      auto n       = std::min(end - begin, static_cast<int>(1));
      auto workers = loop::create_workers<UnaryFunction>(f, begin, end, n);
      for ( int k = 0; k < n; ++k )
        this->push([=] { workers->at(k).run(workers); });
      this->wait();
    }

    //! @brief computes a iterator-based parallel for loop.
    //!
    //! Waits until all tasks have finished, unless called from a thread
    //! that didn't create the pool. If this is taken into account, parallel
    //! loops can be nested.
    //!
    //! @param items an object allowing for `std::begin()` and `std::end()`.
    //! @param f function to be applied as `f(*it)` for the iterator in the
    //! range `[begin, end)` (the 'loop body').
    template <typename Items, typename UnaryFunction>
    inline
    void
    parallel_for_each( Items & items, UnaryFunction f ) {
      auto begin = std::begin(items);
      auto size  = std::distance(begin, std::end(items));
      this->parallel_for(0, size, [=](int i) { f(begin[i]); });
    }

    //! @brief waits for all jobs currently running on the thread
    //! pool. Has no effect when called from threads other than the one that
    //! created the pool.
    //! @param millis if > 0: stops waiting after millis ms.
    void wait(size_t millis = 0) { task_manager_.wait_for_finish(millis); }

    //! @brief checks whether all jobs are done.
    bool done() const { return task_manager_.done(); }

    //! @brief allocator respecting memory alignment.
    static
    void*
    operator new(size_t count) {
      return mem::aligned::alloc(alignof(ThreadPool), count);
    }

    //! @brief deallocator respecting memory alignment.
    static
    void
    operator delete(void* ptr) { mem::aligned::free(ptr); }

  private:
    //! joins all worker threads.
    void
    join_threads() {
      for ( auto & worker : workers_)
        if (worker.joinable())
          worker.join();
    }

    //! adds one worker thread to the thread pool.
    //! @param id worker id (used for matching threads with queues and cores)
    void
    add_worker( size_t id ) {
      workers_[id] = std::thread([&, id] {
        std::function<void()> task;
        while ( !task_manager_.stopped() ) {
          task_manager_.wait_for_jobs(id);
          do {
            // inner while to save some time calling done()
            while ( task_manager_.try_pop(task, id) ) this->execute_safely(task);
          } while (!task_manager_.done());
        }
      });

      // set thread affinity on linux
      this->set_thread_affinity(id);
    }

    //! sets thread affinity (if there are as less workers than cores).
    //! This works on linux by default. In OSX compile with -pthread -DAFFINITY.
    void
    set_thread_affinity( size_t id ) {
#if (defined __linux__ || defined AFFINITY)
      auto hardware_cores = std::thread::hardware_concurrency();
      if ( id >= hardware_cores ) id = id % hardware_cores;
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(id, &cpuset);
      int rc = pthread_setaffinity_np(
        workers_.at(id).native_handle(), sizeof(cpu_set_t), &cpuset
      );
      if (rc != 0) {
        throw std::runtime_error("Error calling pthread_setaffinity_np");
      }
#endif
    }

    void
    execute_safely( std::function<void()>& task ) {
      try {
        task();
        task_manager_.report_success();
      } catch (...) {
        task_manager_.report_fail(std::current_exception());
      }
    }

    sched::TaskManager task_manager_;
    std::vector<std::thread> workers_;
    std::atomic_size_t active_threads_;
  };

  // 5. ---------------------------------------------------

  //! Free-standing functions (main API)

  //! @brief push a job to the global thread pool.
  //! @param f a function.
  //! @param args (optional) arguments passed to `f`.
  template <typename Function, class... Args>
  inline
  void
  push( Function&& f, Args&&... args ) {
    ThreadPool::global_instance().push( std::forward<Function>(f), std::forward<Args>(args)...);
  }

  #ifdef UTILS_QP_USE_FUTURE
  //! @brief executes a job asynchronously the global thread pool.
  //! @param f a function.
  //! @param args (optional) arguments passed to `f`.
  //! @return A `std::future` for the task. Call `future.get()` to retrieve
  //! the results at a later point in time (blocking).
  template<class Function, class... Args>
  inline
  auto
  async(Function&& f, Args&&... args) -> std::future<decltype(f(args...))> {
    return ThreadPool::global_instance().async(std::forward<Function>(f), std::forward<Args>(args)...);
  }
  #endif

  //! @brief waits for all jobs currently running on the global thread pool.
  //! Has no effect when not called from main thread.
  inline
  void
  wait() { ThreadPool::global_instance().wait(); }

  //! @brief checks whether all globel jobs are done.
  inline
  bool
  done() { return ThreadPool::global_instance().done(); }

  //! @brief sets the number of active worker threads in the global thread pool.
  //! @param threads the number of worker threads.
  //! Has no effect when not called from main thread.
  inline
  void
  set_active_threads( size_t threads )
  { ThreadPool::global_instance().set_active_threads(threads); }

  //! @brief retrieves the number of active worker threads in the global thread
  //! pool.
  inline
  size_t
  get_active_threads()
  { return ThreadPool::global_instance().get_active_threads(); }

  //! @brief computes an index-based parallel for loop.
  //!
  //! Waits until all tasks have finished, unless called from a thread that
  //! didn't create the pool. If this is taken into account, parallel loops
  //! can be nested.
  //!
  //! @param begin first index of the loop.
  //! @param end the loop runs in the range `[begin, end)`.
  //! @param f a function taking `int` argument (the 'loop body').
  template <typename UnaryFunction>
  inline
  void
  parallel_for( int begin, int end, UnaryFunction&& f ) {
    ThreadPool::global_instance().parallel_for( begin, end, std::forward<UnaryFunction>(f) );
  }

  //! @brief computes a iterator-based parallel for loop.
  //!
  //! Waits until all tasks have finished, unless called from a thread that
  //! didn't create the pool. If this is taken into account, parallel loops
  //! can be nested.
  //!
  //! @param items an object allowing for `std::begin()` and `std::end()`.
  //! @param f function to be applied as `f(*it)` for the iterator in the
  //! range `[begin, end)` (the 'loop body').
  template <typename Items, typename UnaryFunction>
  inline
  void
  parallel_for_each( Items& items, UnaryFunction&& f ) {
    ThreadPool::global_instance().parallel_for_each( items, std::forward<UnaryFunction>(f));
  }

} // end namespace quickpool
