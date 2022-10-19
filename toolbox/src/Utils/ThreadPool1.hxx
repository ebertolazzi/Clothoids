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

namespace Utils {

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  class ThreadPool1 : public ThreadPoolBase {

    /*\
     |  __        __         _
     |  \ \      / /__  _ __| | _____ _ __
     |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
     |    \ V  V / (_) | |  |   <  __/ |
     |     \_/\_/ \___/|_|  |_|\_\___|_|
    \*/

    class Worker {

      using real_type = double;

      bool                  m_active;
      UTILS_SEMAPHORE       m_running;
      std::thread           m_running_thread;
      std::function<void()> m_job;

      TicToc m_tm;

      unsigned  m_n_job   = 0;
      real_type m_job_ms  = 0;
      real_type m_wait_ms = 0;
      real_type m_push_ms = 0;

      void worker_loop();

    public:

      //disable copy
      Worker( Worker const & )              = delete;
      //Worker( Worker && )                   = delete;
      Worker& operator = ( Worker const & ) = delete;
      Worker& operator = ( Worker && )      = delete;

      Worker() : m_active(false) { start(); }
      ~Worker() { stop(); }

      Worker( Worker && rhs ) noexcept
      : m_active(rhs.m_active)
      , m_running_thread(std::move(rhs.m_running_thread))
      , m_job(std::move(rhs.m_job))
      {}

      void start();
      void stop();

      //!
      //! wait task is done
      //!
      void wait() { m_running.wait_red(); }

      void
      exec( std::function<void()> & fun ) {
        m_tm.tic();
        m_running.wait_red(); // se gia occupato in task aspetta
        m_job = fun;          // cambia funzione da eseguire
        m_running.green();    // activate computation
        m_tm.toc();
        m_push_ms += m_tm.elapsed_ms();
      }

      std::thread::id     get_id()     const { return m_running_thread.get_id(); }
      std::thread const & get_thread() const { return m_running_thread; }
      std::thread &       get_thread()       { return m_running_thread; }

      unsigned  n_job()   const { return m_n_job; }
      real_type job_ms()  const { return m_job_ms; }
      real_type wait_ms() const { return m_wait_ms; }
      real_type push_ms() const { return m_push_ms; }

    };
    std::size_t m_thread_to_send = 0;

    // need to keep track of threads so we can join them
    std::vector<Worker> m_workers;

    void setup() { for ( auto & w: m_workers ) w.start(); }

  public:

    explicit
    ThreadPool1(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    );

    virtual ~ThreadPool1();

    void
    exec( std::function<void()> && fun ) override {
      m_workers[m_thread_to_send].exec( fun );
      if ( ++m_thread_to_send >= m_workers.size() ) m_thread_to_send = 0;
    }

    void wait() override;
    void join() override { stop(); }

    unsigned
    thread_count() const override
    { return unsigned(m_workers.size()); }

    void resize( unsigned numThreads ) override;
    char const * name() const override { return "ThreadPool1"; }

    // EXTRA
    void start();
    void stop();

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

    void info( ostream_type & s ) const override;
  };

  using ThreadPool = ThreadPool1;

}

///
/// eof: ThreadPool1.hxx
///
