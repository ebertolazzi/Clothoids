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
/// eof: ThreadPoolBase.hxx
///

namespace Utils {

  /*\
   |   _____ _                    _ ___          _ ___
   |  |_   _| |_  _ _ ___ __ _ __| | _ \___  ___| | _ ) __ _ ___ ___
   |    | | | ' \| '_/ -_) _` / _` |  _/ _ \/ _ \ | _ \/ _` (_-</ -_)
   |    |_| |_||_|_| \___\__,_\__,_|_| \___/\___/_|___/\__,_/__/\___|
  \*/

  class ThreadPoolBase {

  public:

    //disable copy
    ThreadPoolBase( ThreadPoolBase const & )               = delete;
    ThreadPoolBase( ThreadPoolBase && )                    = delete;
    ThreadPoolBase & operator = ( ThreadPoolBase const & ) = delete;
    ThreadPoolBase & operator = ( ThreadPoolBase && )      = delete;

    ThreadPoolBase() = default;

    virtual
    void
    exec( std::function<void()> && ) = 0;

    template <typename Func, typename... Args>
    void
    run( Func && func, Args && ... args ) {
      this->exec(
        std::bind(
          std::forward<Func>(func),
          std::forward<Args>(args)...
        )
      );
    }

    virtual void         wait() = 0;
    virtual void         join() = 0;
    virtual unsigned     thread_count() const = 0;
    virtual void         resize( unsigned numThreads ) = 0;
    virtual char const * name() const = 0;
    virtual void         info( ostream_type & ) const { }
  };

  namespace tp {

    /*\
     |   ___ _            _  ___                    _ _         ___
     |  | __(_)_ _____ __| |/ __|__ _ _ __  __ _ __(_) |_ _  _ / _ \ _  _ ___ _  _ ___
     |  | _|| \ \ / -_) _` | (__/ _` | '_ \/ _` / _| |  _| || | (_) | || / -_) || / -_)
     |  |_| |_/_\_\___\__,_|\___\__,_| .__/\__,_\__|_|\__|\_, |\__\_\\_,_\___|\_,_\___|
     |                               |_|                  |__/
    \*/
    /*\
        If we would use a deque, we would have to protect
        against overlapping accesses to the front and the
        back. The standard containers do not allow this. Use a
        vector instead.  With a vector it is possible to access
        both ends of the queue at the same time, as push()ing
        and pop()ing does not modify the container itself but
        only its elements.
    \*/
    template <class Function>
    class FixedCapacityQueue {

      union Fun {
        Function m_fun; // Only used between pop_ptr and push_ptr
        Fun() noexcept { }
        Fun( Fun const & ) noexcept { }
        Fun( Fun && ) noexcept { }
        ~Fun() noexcept { }
      };

      std::vector<Fun> m_fun_vec;
      unsigned m_size     = 0;
      unsigned m_capacity = 0;
      unsigned m_push_ptr = 0;
      unsigned m_pop_ptr  = 0;

    public:

      FixedCapacityQueue( FixedCapacityQueue const & )              = delete;
      FixedCapacityQueue( FixedCapacityQueue && )                   = delete;
      FixedCapacityQueue& operator = ( FixedCapacityQueue const & ) = delete;
      FixedCapacityQueue& operator = ( FixedCapacityQueue && )      = delete;

      explicit
      FixedCapacityQueue( unsigned capacity )
      : m_fun_vec( size_t( capacity+1 ) )
      , m_size( capacity+1 )
      , m_capacity( capacity )
      { }

      void
      push( Function && f ) {
        new (&m_fun_vec[m_push_ptr].m_fun) Function(std::forward<Function>(f));
        if ( ++m_push_ptr == m_size ) m_push_ptr = 0;
      }

      Function
      pop() {
        Function r = std::move(m_fun_vec[m_pop_ptr].m_fun);
        m_fun_vec[m_pop_ptr].m_fun.~Function();
        if ( ++m_pop_ptr == m_size ) m_pop_ptr = 0;
        return r;
      }

      unsigned size()     const { return ((m_push_ptr + m_size) - m_pop_ptr) % m_size; }
      bool     empty()    const { return m_push_ptr == m_pop_ptr; }
      bool     is_full()  const { return this->size() >= m_capacity; }
      unsigned capacity() const { return m_capacity; }

      void
      reserve( unsigned capacity ) {
        assert(empty()); // Copying / moving of Fun not supported.
        if ( capacity != m_capacity ) {
          m_size     = capacity+1;
          m_capacity = capacity;
          m_fun_vec.resize( m_size );
        }
      }

      ~FixedCapacityQueue() { while (!empty()) pop(); }
    };

    /*\
     |    ___
     |   / _ \ _   _  ___ _   _  ___
     |  | | | | | | |/ _ \ | | |/ _ \
     |  | |_| | |_| |  __/ |_| |  __/
     |   \__\_\\__,_|\___|\__,_|\___|
    \*/

    class Queue {
    public:

      class TaskData {
        std::function<void()> m_fun;
      public:
        TaskData( std::function<void()> && f ) : m_fun(std::move(f)) { }
        TaskData( std::function<void()> & f ) : m_fun(f) { }
        void operator()() { m_fun(); delete this; }
        ~TaskData() = default;
      };

    private:

      std::vector<TaskData*> m_queue_data;

      unsigned m_size, m_capacity, m_push_ptr, m_pop_ptr;

    public:

      Queue( Queue const & )              = delete;
      Queue( Queue && )                   = delete;
      Queue& operator = ( Queue const & ) = delete;
      Queue& operator = ( Queue && )      = delete;

      explicit
      Queue( unsigned capacity )
      : m_queue_data( std::size_t( capacity+1 ) )
      , m_size( capacity+1 )
      , m_capacity( capacity )
      , m_push_ptr( 0 )
      , m_pop_ptr( 0 )
      { }

      void
      push( TaskData * task ) {
        m_queue_data[m_push_ptr] = task;
        if ( ++m_push_ptr == m_size ) m_push_ptr = 0;
      }

      TaskData *
      pop() {
        unsigned ipos = m_pop_ptr;
        if ( ++m_pop_ptr == m_size ) m_pop_ptr = 0;
        return m_queue_data[ipos];
      }

      unsigned size()     const { return ((m_push_ptr + m_size) - m_pop_ptr) % m_size; }
      bool     empty()    const { return m_push_ptr == m_pop_ptr; }
      bool     is_full()  const { return this->size() >= m_capacity; }
      unsigned capacity() const { return m_capacity; }

      //! clear queue and delete tasks
      void clear() { while( !empty() ) delete pop(); }

      void
      resize( unsigned capacity ) {
        this->clear();
        m_size     = capacity+1;
        m_capacity = capacity;
        m_queue_data.resize( m_size );
      }

      ~Queue() = default;
    };
  }
}

///
/// eof: ThreadPoolBase.hxx
///
