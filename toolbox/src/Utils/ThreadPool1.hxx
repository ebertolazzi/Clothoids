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
/// file: ThreadPool1.hxx
///

#ifdef UTILS_OS_LINUX
  #include <pthread.h>
#endif

namespace Utils {

  /*\
   |  __        __         _
   |  \ \      / /__  _ __| | _____ _ __
   |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
   |    \ V  V / (_) | |  |   <  __/ |
   |     \_/\_/ \___/|_|  |_|\_\___|_|
  \*/

  class Worker {
    //friend class ThreadPool1;

    bool                  m_active;
    SimpleSemaphore       m_is_running;
    std::thread           m_running_thread;
    std::function<void()> m_job;

    //disable copy
    Worker( Worker const & )              = delete;
    //Worker( Worker && )                   = delete;
    Worker& operator = ( Worker const & ) = delete;
    Worker& operator = ( Worker && )      = delete;

    void
    worker_loop() {
      while ( m_active ) {
        m_is_running.wait();
        m_job();
        m_is_running.red();
      }
    }

  public:

    Worker() : m_active(false) { start(); }
    ~Worker() { stop(); }

    Worker( Worker && rhs ) {
      m_active         = rhs.m_active;
      m_job            = std::move(rhs.m_job);
      m_running_thread = std::move(rhs.m_running_thread);
    }

    void
    start() {
      if ( !m_active ) {
        m_active = true;
        m_is_running.red();
        m_running_thread = std::thread(&Worker::worker_loop,this);
      }
    }

    void
    stop() {
      if ( m_active ) {
        m_active = false;        // deactivate computation
        m_is_running.wait_red(); // se gia occupato in task aspetta
        m_job = [](){};          // dummy task
        m_is_running.green();    // start computation
        if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
      }
    }

    //!
    //! wait task is done
    //!
    void wait() { m_is_running.wait_red(); }

    void
    exec( std::function<void()> & fun ) {
      m_is_running.wait_red(); // se gia occupato in task aspetta
      m_job = fun;             // cambia funzione da eseguire
      m_is_running.green();    // activate computation
    }

    template < class Func, class... Args >
    void
    run( Func && func, Args && ... args ) {
      m_is_running.wait_red(); // se gia occupato in task aspetta
      m_job = std::bind(std::forward<Func>(func),std::forward<Args>(args)...);
      m_is_running.green();    // activate computation
    }

    std::thread::id     get_id()     const { return m_running_thread.get_id(); }
    std::thread const & get_thread() const { return m_running_thread; }
    std::thread &       get_thread()       { return m_running_thread; }

  };

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  class ThreadPool1 : public ThreadPoolBase {

    std::size_t m_thread_to_send = 0;

    // need to keep track of threads so we can join them
    std::vector<Worker> m_workers;

    #ifdef UTILS_OS_LINUX
    void
    setup() {
      sched_param sch;
      int         policy;
      for ( auto & w: m_workers ) {
        w.start();
        std::thread & t = w.get_thread();
        pthread_getschedparam( t.native_handle(), &policy, &sch );
        sch.sched_priority = sched_get_priority_max( SCHED_RR );
        pthread_setschedparam( t.native_handle(), SCHED_RR, &sch );
      }
    }
    #else
    void
    setup() {
      for ( auto & w: m_workers ) w.start();
    }
    #endif

  public:

    ThreadPool1(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    )
    : ThreadPoolBase()
    {
      m_workers.resize( size_t( nthread ) );
      setup();
    }

    virtual ~ThreadPool1() { join(); m_workers.clear(); }

    void
    exec( std::function<void()> && fun ) override {
      m_workers[m_thread_to_send].exec( fun );
      if ( ++m_thread_to_send >= m_workers.size() ) m_thread_to_send = 0;
    }

    void
    wait() override {
      m_thread_to_send = 0;
      for ( auto && w : m_workers ) w.wait();
    }

    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    void
    resize( unsigned numThreads ) override {
      wait();
      stop();
      m_workers.resize( size_t(numThreads) );
      setup();
    }

    char const * name() const override { return "ThreadPool1"; }

    // EXTRA

    void start() { m_thread_to_send = 0; for ( auto && w : m_workers ) w.start(); }
    void stop()  { m_thread_to_send = 0; for ( auto && w : m_workers ) w.stop(); }
    void join()  { stop(); }

    std::thread::id
    get_id( unsigned i ) const
    { return m_workers[size_t(i)].get_id(); }

    std::thread const &
    get_thread( unsigned i ) const
    { return m_workers[size_t(i)].get_thread(); }

    std::thread &
    get_thread( unsigned i )
    { return m_workers[size_t(i)].get_thread(); }

    // ALIAS
    void wait_all()  { this->wait();  }
    void start_all() { this->start(); }
    void stop_all()  { this->stop();  }
    unsigned size() const { return this->thread_count(); }

  };

  using ThreadPool = ThreadPool1;

}

///
/// eof: ThreadPool1.hxx
///
