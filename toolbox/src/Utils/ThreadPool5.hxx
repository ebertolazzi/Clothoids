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
/// file: ThreadPool5.hxx
///

#ifdef UTILS_OS_LINUX
  #include <pthread.h>
#endif

namespace Utils {

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  class ThreadPool5 : public ThreadPoolBase {

    /*\
     |  __        __         _
     |  \ \      / /__  _ __| | _____ _ __
     |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
     |    \ V  V / (_) | |  |   <  __/ |
     |     \_/\_/ \___/|_|  |_|\_\___|_|
    \*/

    class Worker {

      bool                  m_active           = false;
      unsigned              m_job_done_counter = 0;
      unsigned              m_worker_id        = 0;
      ThreadPool5 *         m_tp               = nullptr;
      SimpleSemaphore       m_is_running;
      std::thread           m_running_thread;
      std::function<void()> m_job;

      void
      worker_loop() {
        while ( m_active ) {
          m_is_running.red();     // block computation
          m_is_running.wait();    // wait signal to start computation
          if ( !m_active ) break; // if finished exit
          m_job();
          ++m_job_done_counter;
          m_tp->push_worker( m_worker_id ); // worker ready for a new computation
        }
        //fmt::print( "worker_loop {} exiting\n", m_worker_id );
      }

    public:

      explicit Worker() { start(); }
      ~Worker() { stop(); }

      // dummy copy constructor for resize
      Worker( Worker && ) {}

      void
      setup( ThreadPool5 * tp, unsigned id ) {
        m_worker_id        = id;
        m_job_done_counter = 0;
        m_tp               = tp;
      }

      void
      start() {
        if ( !m_active ) {
          m_active = true;
          m_running_thread = std::thread( &Worker::worker_loop, this );
        }
      }

      //!
      //! wait task is done
      //!
      void wait() { m_is_running.wait_red(); }

      void
      stop() {
        if ( m_active ) {
          wait();               // if running task wait it terminate
          m_active = false;     // deactivate computation
          m_job = [](){};       // dummy task
          m_is_running.green(); // start computation (exiting loop)
          if ( m_running_thread.joinable() ) m_running_thread.join(); // wait thread for exiting
          m_is_running.red();   // end of computation (for double stop);
        }
        //fmt::print( "worker_loop {} stopped\n", m_worker_id );
      }

      void
      exec( std::function<void()> & fun ) {
        m_is_running.wait_red();
        m_job = fun;             // cambia funzione da eseguire
        m_is_running.green();    // activate computation
      }

      unsigned job_done_counter() const { return m_job_done_counter; }

      void
      info( ostream_type & s ) const {
        fmt::print(
          s,"Worker {}, #job = {}\n", m_worker_id, m_job_done_counter
        );
      }
    };

    // =========================================================================
    // =========================================================================
    // =========================================================================

    // need to keep track of threads so we can join them
    std::vector<Worker>     m_workers;
    // stack of available workers
    std::vector<unsigned>   m_stack;
    std::mutex              m_stack_mutex;
    std::condition_variable m_stack_cond;

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
    void setup() { for ( auto & w: m_workers ) w.start(); }
    #endif

    void
    resize_workers( unsigned numThreads ) {
      m_stack.clear(); // empty stack
      m_stack.reserve( size_t(numThreads) );
      m_workers.resize( size_t(numThreads) );
      unsigned id = 0;
      for ( Worker & w : m_workers ) { w.setup( this, id ); ++id; }
      while ( id-- > 0 ) push_worker( id );
      setup();
    }

    void
    push_worker( unsigned id ) {
      std::unique_lock<std::mutex> lock(m_stack_mutex);
      m_stack.push_back(id);
      m_stack_cond.notify_one();
    }

    unsigned
    pop_worker() {
      std::unique_lock<std::mutex> lock(m_stack_mutex);
      m_stack_cond.wait( lock, [&]()->bool { return !m_stack.empty(); } );
      unsigned id = m_stack.back(); m_stack.pop_back();
      return id;
    }

  public:

    ThreadPool5(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    )
    : ThreadPoolBase()
    {
      resize_workers( nthread );
      //info( std::cout );
    }

    virtual
    ~ThreadPool5() {
      join();
      m_workers.clear();
      m_stack.clear();
    }

    void
    exec( std::function<void()> && fun ) override {
      // cerca prima thread libera
      unsigned id = pop_worker();
      assert( id >= 0 && id < m_workers.size() );
      m_workers[id].exec( fun );
    }

    void
    wait() override
    { for ( auto & w : m_workers ) w.wait(); }

    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    void
    resize( unsigned numThreads ) override
    { wait(); stop(); resize_workers( numThreads ); }

    char const * name() const override { return "ThreadPool5"; }

    void start() { for ( auto && w : m_workers ) w.start(); }
    void stop()  { for ( auto && w : m_workers ) w.stop(); }
    void join()  { stop(); }

    void
    info_stack( ostream_type & s ) const {
      fmt::print( s, "STACK[{}]: ", m_stack.size() );
      for ( unsigned const & id : m_stack )
        fmt::print( s, "{}, ", id );
      s << '\n';
    }

    void
    info( ostream_type & s ) const override {
      for ( Worker const & w : m_workers ) w.info(s);
      info_stack( s );
    }
  };
}

///
/// eof: ThreadPool5.hxx
///
