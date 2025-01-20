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

//
// file: ThreadPool0.hxx
//

#ifdef UTILS_OS_LINUX
  #include <pthread.h>
#endif

namespace Utils {

  /*!
   * \addtogroup THREAD
   * @{
   */

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  //!
  //!  \brief Fake thread pool class.
  //!
  //!  This class simulates a thread pool by executing tasks
  //!  in the calling thread, providing a simple interface for
  //!  task execution without actual multi-threading capabilities.
  //!
  //!  It is primarily used for testing or environments where
  //!  threading is not needed.
  //!
  class ThreadPool0 : public ThreadPoolBase {
    using Func = std::function<void()>;

  public:

    //!
    //!  \brief Constructs a fake thread pool with a specified number of threads.
    //!
    //!  This constructor initializes the thread pool with a specified number
    //!  of threads, although it does not actually create any threads.
    //!  The specified number is ignored.
    //!
    explicit
    ThreadPool0( unsigned )
    : ThreadPoolBase()
    {}

    //! Destructor.
    virtual ~ThreadPool0() = default;

    //!
    //! \brief Executes a given function in the current thread.
    //!
    //! This method takes a function and executes it immediately
    //! in the calling thread. No actual thread pooling is performed.
    //!
    //! \param fun The function to be executed.
    //!
    void exec( Func && fun ) override { fun(); }

    //!
    //! \brief No-op method to wait for all tasks to complete.
    //!
    //! This method does not perform any actions, as there are no
    //! threads or tasks to wait for.
    //!
    void wait() override { }

    //!
    //! \brief No-op method to join all threads.
    //!
    //! This method does not perform any actions, as there are no
    //! threads to join.
    //!
    void join() override { }

    //!
    //! \brief Returns the number of threads in the pool.
    //!
    //! This method returns the number of threads in the pool,
    //! which is always 1 for this fake thread pool.
    //!
    //! \return Always returns 1.
    //!
    unsigned thread_count() const override { return 1; }

    //!
    //! \brief Resizes the pool to a specified number of threads.
    //!
    //! This method does not perform any actions, as the fake
    //! thread pool does not support resizing.
    //!
    void resize( unsigned ) override { }

    //!
    //! \brief Returns the name of the thread pool.
    //!
    //! This method returns a string describing the type of thread pool.
    //!
    //! \return A string indicating that this is a "ThreadPool0 (fake thread)".
    //!
    char const * name() const override { return "ThreadPool0 (fake thread)"; }

    //!
    //! \brief Returns the size of the thread pool.
    //!
    //! This method returns the size of the thread pool, which is always 1.
    //!
    //! \return Always returns 1.
    //!
    unsigned size() const { return 1; }
  };

  /*! @} */

}

//
// eof: ThreadPool0.hxx
//
