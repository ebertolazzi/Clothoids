/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
 |                                                                          |
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

///
/// eof: ThreadUtils.hxx
///

namespace Utils {

  #ifdef UTILS_OS_WINDOWS
    #define UTILS_SEMAPHORE Utils::WinSemaphore
    #define UTILS_MUTEX     Utils::WinCriticalSection
    #define UTILS_SPINLOCK  Utils::WinCriticalSection
    #define UTILS_BARRIER   Utils::WinBarrier
  #else
    #define UTILS_SEMAPHORE Utils::SimpleSemaphore
    #define UTILS_MUTEX     std::mutex
    #define UTILS_SPINLOCK  Utils::SpinLock
    #define UTILS_BARRIER   Utils::Barrier
  #endif

  #ifdef UTILS_OS_WINDOWS

  class WinMutex {
    HANDLE m_mutex;

  public:

    WinMutex() : m_mutex(NULL) {
      m_mutex = CreateMutex(
        NULL,  // no security descriptor
        FALSE, // mutex not owned
        NULL   // object name
      );
      UTILS_ASSERT(
        m_mutex != NULL,
        "WinMutex(): error: {}.\n", GetLastError()
      );
    }

    ~WinMutex() { CloseHandle(m_mutex); }

    void
    lock() {
    	DWORD res = WaitForSingleObject(m_mutex, INFINITE);
      UTILS_ASSERT0( res == WAIT_OBJECT_0, "WinMutex::lock, WAIT_TIMEOUT" );
    }

    void
    unlock() {
    	DWORD res = ReleaseMutex(m_mutex);
      UTILS_ASSERT0( res == WAIT_OBJECT_0, "WinMutex::lock, WAIT_TIMEOUT" );
    }

  };

  class WinCriticalSection {
    CRITICAL_SECTION m_critical;
  public:
    WinCriticalSection()
    { InitializeCriticalSection(&m_critical); }
    ~WinCriticalSection() { DeleteCriticalSection(&m_critical); }
    void lock() {	EnterCriticalSection(&m_critical); }
    void unlock() {	LeaveCriticalSection(&m_critical); }
    bool try_lock() {	return TRUE == TryEnterCriticalSection(&m_critical); }
    void wait() {	lock(); unlock(); }
    CRITICAL_SECTION const & data() const { return m_critical; }
    CRITICAL_SECTION       & data()       { return m_critical; }
  };

  class WinSemaphore {
    bool               m_is_red;
    WinCriticalSection m_critical;
    CONDITION_VARIABLE m_condition;
    unsigned           m_waiting_green;
    unsigned           m_waiting_red;

    //void notify_one() noexcept { ReleaseSemaphore( m_semaphore, 1, NULL ); }
    void notify_one() noexcept { WakeConditionVariable( &m_condition ); }
    void notify_all() noexcept { WakeAllConditionVariable( &m_condition ); }
    void wait_cond() noexcept { SleepConditionVariableCS( &m_condition, &m_critical.data(), INFINITE ); }

  public:

    WinSemaphore()
    : m_is_red(false)
    , m_waiting_green(0)
    , m_waiting_red(0)
    { InitializeConditionVariable( &m_condition ); }

    ~WinSemaphore() {}

    //!
    //! unblock semaphore
    //!
    void
    green() noexcept {
      m_critical.lock();
      m_is_red = false;
      m_critical.unlock();
      if      ( m_waiting_green > 1 ) notify_all();
      else if ( m_waiting_green > 0 ) notify_one();
    }

    //!
    //! block semaphore
    //!
    void
    red() noexcept {
      m_critical.lock();
      m_is_red = true;
      m_critical.unlock();
      if      ( m_waiting_red > 1 ) notify_all();
      else if ( m_waiting_red > 0 ) notify_one();
    }

    //!
    //! wait until m_count <= 0
    //!
    void
    wait() {
      m_critical.lock();
      ++m_waiting_green;
      while ( m_is_red ) wait_cond();
      --m_waiting_green;
      m_critical.unlock();
    }

    //!
    //! wait until m_count > 0
    //!
    void
    wait_red() {
      m_critical.lock();
      ++m_waiting_red;
      while ( !m_is_red ) wait_cond();
      --m_waiting_red;
      m_critical.unlock();
    }

  };

  class WinBarrier {
    int                m_to_be_done;
    WinCriticalSection m_critical;
    CONDITION_VARIABLE m_condition;

    void notify_one() noexcept { WakeConditionVariable( &m_condition ); }
    void notify_all() noexcept { WakeAllConditionVariable( &m_condition ); }
    void wait_cond() noexcept { SleepConditionVariableCS( &m_condition, &m_critical.data(), INFINITE ); }

  public:
    WinBarrier() : m_to_be_done(0) {}

    void
    setup( int nthreads )
    { m_to_be_done = nthreads ; }

    void
    count_down() {
      m_critical.lock();
      if ( --m_to_be_done <= 0 ) notify_all() ; // wake up all tread
      m_critical.unlock();
    }

    void
    wait() {
      m_critical.lock();
      while( m_to_be_done > 0 ) wait_cond();
      m_critical.unlock();
    }

    void
    count_down_and_wait() {
      m_critical.lock();
      if ( --m_to_be_done <= 0 ) {
        notify_all() ; // wake up all tread
      } else {
        while( m_to_be_done > 0 ) wait_cond();
      }
      m_critical.unlock();
    }
  };

  #endif

  /*\
   |   ___      _      _            _
   |  / __|_ __(_)_ _ | |   ___  __| |__
   |  \__ \ '_ \ | ' \| |__/ _ \/ _| / /
   |  |___/ .__/_|_||_|____\___/\__|_\_\
   |      |_|
  \*/

  class SpinLock {
  private:
    std::atomic<bool> m_lock;
  public:
    SpinLock() : m_lock(false) { }
    SpinLock( SpinLock const & ) = delete;
    ~SpinLock() = default;

    void
    wait() {
      while( m_lock.load(std::memory_order_acquire) )
        std::this_thread::yield();
    }

    void
    lock() {
      while( m_lock.exchange(true, std::memory_order_acquire) )
        std::this_thread::yield();
    }

    bool
    try_lock() {
      return !m_lock.exchange(true, std::memory_order_acquire);
    }

    void
    unlock() {
      m_lock.store(false, std::memory_order_release);
    }
  };

  /*\
   |   ___      _      _            _     _                  _
   |  / __|_ __(_)_ _ | |   ___  __| |__ | |__  __ _ _ _ _ _(_)___ _ _
   |  \__ \ '_ \ | ' \| |__/ _ \/ _| / / | '_ \/ _` | '_| '_| / -_) '_|
   |  |___/ .__/_|_||_|____\___/\__|_\_\_|_.__/\__,_|_| |_| |_\___|_|
   |      |_|                         |___|
  \*/

  class SpinLock_barrier {
  private:
    std::atomic<unsigned> m_count;
    std::atomic<unsigned> m_generation;
    unsigned              m_count_reset_value;
  public:
    SpinLock_barrier( SpinLock_barrier const & ) = delete;
    SpinLock_barrier& operator=( SpinLock_barrier const & ) = delete;

    explicit
    SpinLock_barrier()
    : m_count(0)
    , m_generation(0)
    , m_count_reset_value(0)
    {}

    void
    setup( unsigned count ) {
      m_count_reset_value = m_count = count ;
    }

    void
    count_down() {
      unsigned gen = m_generation.load();
      if ( --m_count == 0 ) {
        if ( m_generation.compare_exchange_weak(gen, gen + 1) )
          m_count = m_count_reset_value;
        return;
      }
    }

    void
    wait() {
      unsigned gen = m_generation.load();
      while ((gen == m_generation) && (m_count != 0))
        std::this_thread::yield();
    }

    void
    count_down_and_wait() {
      unsigned gen = m_generation.load();
      if ( --m_count == 0 ) {
        if ( m_generation.compare_exchange_weak(gen, gen + 1) )
          m_count = m_count_reset_value;
        return;
      }
      while ((gen == m_generation) && (m_count != 0))
        std::this_thread::yield();
    }
  };

  /*\
   |   ___               _
   |  | _ ) __ _ _ _ _ _(_)___ _ _
   |  | _ \/ _` | '_| '_| / -_) '_|
   |  |___/\__,_|_| |_| |_\___|_|
  \*/

  class Barrier {
    int                     m_to_be_done;
    std::mutex              m_mtx;
    std::condition_variable m_cond;
  public:
    Barrier() : m_to_be_done(0) {}

    void
    setup( int nthreads )
    { m_to_be_done = nthreads ; }

    void
    count_down() {
      std::unique_lock<std::mutex> lck(m_mtx);
      if ( --m_to_be_done <= 0 ) m_cond.notify_all() ; // wake up all tread
    }

    void
    wait() {
      std::unique_lock<std::mutex> lck(m_mtx);
      m_cond.wait( lck, [&]()->bool { return m_to_be_done <= 0; } );
    }

    void
    count_down_and_wait() {
      std::unique_lock<std::mutex> lck(m_mtx);
      if ( --m_to_be_done <= 0 ) {
        m_cond.notify_all() ; // wake up all tread
      } else {
        m_cond.wait( lck, [&]()->bool { return m_to_be_done <= 0; } );
      }
    }
  };

  /*\
   |   ___ _            _     ___                      _
   |  / __(_)_ __  _ __| |___/ __| ___ _ __  __ _ _ __| |_  ___ _ _ ___
   |  \__ \ | '  \| '_ \ / -_)__ \/ -_) '  \/ _` | '_ \ ' \/ _ \ '_/ -_)
   |  |___/_|_|_|_| .__/_\___|___/\___|_|_|_\__,_| .__/_||_\___/_| \___|
   |              |_|                            |_|
  \*/

  class SimpleSemaphore {
  private:
    bool                        m_is_red;
    unsigned                    m_waiting_green;
    unsigned                    m_waiting_red;
    std::mutex                  m_mutex;
    std::condition_variable_any m_cv_red;
    std::condition_variable_any m_cv_green;
  public:
    SimpleSemaphore() noexcept
    : m_is_red(false)
    , m_waiting_green(0)
    , m_waiting_red(0)
    { }

    //!
    //! unblock semaphore
    //!
    void
    green() noexcept {
      m_mutex.lock();
      m_is_red = false;
      m_mutex.unlock();
      if      ( m_waiting_green > 1 ) m_cv_green.notify_all();
      else if ( m_waiting_green > 0 ) m_cv_green.notify_one();
    }

    //!
    //! block semaphore
    //!
    void
    red() noexcept {
      m_mutex.lock();
      m_is_red = true;
      m_mutex.unlock();
      if      ( m_waiting_red > 1 ) m_cv_red.notify_all();
      else if ( m_waiting_red > 0 ) m_cv_red.notify_one();
    }

    //!
    //! wait until m_count <= 0
    //!
    void
    wait() noexcept {
      m_mutex.lock();
      ++m_waiting_green;
      m_cv_green.wait( m_mutex, [&]()->bool { return !m_is_red; } );
      --m_waiting_green;
      m_mutex.unlock();
    }

    //!
    //! wait until m_count > 0
    //!
    void
    wait_red() noexcept {
      m_mutex.lock();
      ++m_waiting_red;
      m_cv_red.wait( m_mutex, [&]()->bool { return m_is_red; } );
      --m_waiting_red;
      m_mutex.unlock();
    }

  };

  /*\
   |  __        __         _             _
   |  \ \      / /__  _ __| | _____ _ __| |    ___   ___  _ __
   |   \ \ /\ / / _ \| '__| |/ / _ \ '__| |   / _ \ / _ \| '_ \
   |    \ V  V / (_) | |  |   <  __/ |  | |__| (_) | (_) | |_) |
   |     \_/\_/ \___/|_|  |_|\_\___|_|  |_____\___/ \___/| .__/
   |                                                     |_|
  \*/
  class WorkerLoop {

    bool                    m_active;
    bool                    m_running;
    bool                    m_do_job;
    std::thread             m_running_thread;
    std::function<void()>   m_job;

    std::mutex              m_mutex;
    std::condition_variable m_cv;

    void worker_loop();

  public:

    WorkerLoop( WorkerLoop && )                   = delete;
    WorkerLoop( WorkerLoop const & )              = delete;
    WorkerLoop& operator = ( WorkerLoop const & ) = delete;
    WorkerLoop& operator = ( WorkerLoop && )      = delete;

    WorkerLoop();
    ~WorkerLoop();

    void exec( std::function<void()> & fun );
    void exec();
    void wait();

    template <typename Func, typename... Args>
    void
    run( Func && func, Args && ... args ) {
      std::function<void()> f = std::bind( func, std::forward<Args>(args)... );
      this->exec( f );
    }

    std::thread::id     get_id()     const { return m_running_thread.get_id(); }
    std::thread const & get_thread() const { return m_running_thread; }
    std::thread &       get_thread()       { return m_running_thread; }
  };


  /*\
   |  __      __    _ _ __      __       _
   |  \ \    / /_ _(_) |\ \    / /__ _ _| |_____ _ _
   |   \ \/\/ / _` | |  _\ \/\/ / _ \ '_| / / -_) '_|
   |    \_/\_/\__,_|_|\__|\_/\_/\___/_| |_\_\___|_|
  \*/

  class WaitWorker {
  private:
    #ifdef UTILS_OS_WINDOWS
    std::atomic<int> n_worker;
    #else
    std::atomic<int> n_worker = {0};
    #endif
  public:
    #ifdef UTILS_OS_WINDOWS
    WaitWorker() { n_worker = 0; }
    #else
    WaitWorker() = default;
    #endif

    void
    wait()
    { while (n_worker.load(std::memory_order_relaxed) != 0 ){} }

    void enter() { ++n_worker; }
    void leave() { --n_worker; }
  };

  /*\
   |   ___ _                    ___                  _
   |  | _ |_)_ _  __ _ _ _ _  _/ __| ___ __ _ _ _ __| |_
   |  | _ \ | ' \/ _` | '_| || \__ \/ -_) _` | '_/ _| ' \
   |  |___/_|_||_\__,_|_|  \_, |___/\___\__,_|_| \__|_||_|
   |                       |__/
  \*/

  template <typename DATA>
  class BinarySearch {
  private:
    typedef std::pair<std::thread::id,DATA*> DATA_TYPE;
    mutable std::vector<DATA_TYPE>           m_data;
    mutable UTILS_SPINLOCK                   m_spin_write;
    mutable WaitWorker                       m_worker_read;

  public:

    BinarySearch() {
      m_data.clear();
      m_data.reserve(64);
    }

    ~BinarySearch() {
      m_spin_write.wait();
      for ( auto & a : m_data ) delete a.second;
      m_data.clear();
      m_spin_write.wait();
    }

    void
    clear() {
      m_spin_write.wait();
      for ( auto & a : m_data ) delete a.second;
      m_data.clear(); m_data.reserve(64);
      m_spin_write.wait();
    }

    DATA *
    search( std::thread::id const & id, bool & ok ) const {
      m_spin_write.wait(); // wait writing finished
      m_worker_read.enter();
      ok = true;

      size_t U = m_data.size();

      if ( U == 0 ) {
        m_worker_read.leave();
        m_spin_write.lock();
        m_worker_read.wait(); // wait all read finished
        ok = false;
        U  = m_data.size(); // MAI USATO
        m_data.resize(1);
        DATA_TYPE & dL = m_data[0];
        dL.first = id;
        DATA * res = dL.second = new DATA();
        m_spin_write.unlock();
        return res;
      }

      size_t L = 0;
      while ( U-L > 1 ) {
        size_t pos = (L+U)>>1;
        std::thread::id const & id_pos = m_data[pos].first;
        if ( id_pos < id ) L = pos; else U = pos;
      }
      DATA_TYPE & dL = m_data[L];
      if ( dL.first == id ) { m_worker_read.leave(); return dL.second; }
      DATA_TYPE & dU = m_data[U];
      if ( dU.first == id ) { m_worker_read.leave(); return dU.second; }
      m_worker_read.leave();

      // not found must insert
      m_spin_write.lock();
      m_worker_read.wait(); // wait all read finished
      ok = false;
      if ( dL.first < id ) ++L;
      U = m_data.size();
      m_data.resize(U+1);
      while ( U > L ) {
        --U;
        m_data[U+1].first  = m_data[U].first;
        m_data[U+1].second = m_data[U].second;
      }
      DATA_TYPE & dL1 = m_data[L];
      dL1.first = id;
      DATA * res = dL1.second = new DATA();
      m_spin_write.unlock();
      return res;
    }

  };

  /*\
   |         _                                             _ _
   |    __ _| |_     ___  ___ ___  _ __   ___     _____  _(_) |_
   |   / _` | __|   / __|/ __/ _ \| '_ \ / _ \   / _ \ \/ / | __|
   |  | (_| | |_    \__ \ (_| (_) | |_) |  __/  |  __/>  <| | |_
   |   \__,_|\__|___|___/\___\___/| .__/ \___|___\___/_/\_\_|\__|
   |           |_____|            |_|       |_____|
  \*/

  /**
   * Call some function at the end of the current block
   */
  template <class Destructor>
  class at_scope_exit_impl {
    Destructor m_destructor;
    bool       m_active;

  public:

    at_scope_exit_impl( at_scope_exit_impl const & ) = delete;
    at_scope_exit_impl & operator=( at_scope_exit_impl const & ) = delete;

    at_scope_exit_impl() : m_active(false) { }

    explicit
    at_scope_exit_impl( Destructor&& destructor )
    : m_destructor(std::forward<Destructor>(destructor))
    , m_active(true)
    { }

    explicit
    at_scope_exit_impl( Destructor const & destructor )
    : m_destructor(destructor)
    , m_active(true)
    { }

    at_scope_exit_impl(at_scope_exit_impl&& x) noexcept
    : m_destructor(std::move(x.m_destructor))
    , m_active(x.m_active)
    { x.m_active = false; }

    at_scope_exit_impl&
    operator=(at_scope_exit_impl&& x) {
      m_destructor = std::move(x.m_destructor);
      m_active     = x.m_active;
      x.m_active   = false;
    }

    ~at_scope_exit_impl() { if (m_active) m_destructor(); }
  };

  /**
   * Create a variable that when destructed at the end of the scope
   * executes a destructor function.
   *
   * \tparam Destructor&& destructor
   *         The destructor function, maybe a lambda function.
   *
   * Use like this:
   *
   * static int a = 0;
   *
   * { // Enter scope
   *     ++a;
   *     auto x1 = at_scope_exit([&](){ --a; }
   *     // Do something, possibly throwing an exception
   * } // x1 goes out of scope, 'delete a' is called.
   */
  template<class Function>
  auto at_scope_exit(Function&& fun) -> at_scope_exit_impl<Function>
  { return at_scope_exit_impl<Function>(std::forward<Function>(fun)); }

  template<class Function>
  auto at_scope_exit(Function const & fun) -> at_scope_exit_impl<Function const &>
  { return at_scope_exit_impl<Function const &>(fun); }

}

///
/// eof: ThreadUtils.hxx
///
