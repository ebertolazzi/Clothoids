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

namespace Utils {

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  class ThreadPool5 : public ThreadPoolBase {

    typedef double real_type;

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
      UTILS_SEMAPHORE       m_is_running;
      std::thread           m_running_thread;
      std::function<void()> m_job;
      TicToc                m_tm;
      real_type             m_job_ms  = 0;
      real_type             m_sync_ms = 0;
      real_type             m_wait_ms = 0;

      void worker_loop();

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

      void start();

      //!
      //! wait task is done
      //!
      void wait() { m_is_running.wait_red(); }
      void stop();

      void
      exec( std::function<void()> & fun ) {
        m_is_running.wait_red();
        m_job = fun;          // cambia funzione da eseguire
        m_is_running.green(); // activate computation
      }

      unsigned job_done_counter() const { return m_job_done_counter; }

      real_type elapsed_job_ms()  const { return m_job_ms; }
      real_type elapsed_sync_ms() const { return m_sync_ms; }
      real_type elapsed_wait_ms() const { return m_wait_ms; }

      void info( ostream_type & s ) const;
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
    TicToc                  m_tm;
    real_type               m_exec_ms = 0;
    real_type               m_pop_ms  = 0;

    void setup() { for ( auto & w: m_workers ) w.start(); }

    void     resize_workers( unsigned numThreads );
    void     push_worker( unsigned id );
    unsigned pop_worker();

  public:

    ThreadPool5(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    );

    virtual ~ThreadPool5();

    void
    exec( std::function<void()> && fun ) override {
      // cerca prima thread libera

      m_tm.tic();
      Worker & w = m_workers[pop_worker()];
      m_tm.toc();
      m_pop_ms += m_tm.elapsed_ms();

      m_tm.tic();
      w.exec( fun );
      m_tm.toc();
      m_exec_ms += m_tm.elapsed_ms();
    }

    void
    wait() override
    { for ( auto & w : m_workers ) w.wait(); }

    void
    join() override
    { stop(); }

    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    void
    resize( unsigned numThreads ) override
    { wait(); stop(); resize_workers( numThreads ); }

    char const * name() const override { return "ThreadPool5"; }

    void start() { for ( auto && w : m_workers ) w.start(); }
    void stop()  { for ( auto && w : m_workers ) w.stop(); }

    void info_stack( ostream_type & s ) const;
    void info( ostream_type & s ) const override;
  };
}

///
/// eof: ThreadPool5.hxx
///
