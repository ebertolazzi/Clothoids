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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: ThreadPool0.hxx
//

#ifdef UTILS_OS_LINUX
#include <pthread.h>
#endif

namespace Utils
{

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
  class ThreadPool0 : public ThreadPoolBase
  {
    using ThreadPoolBase::FUN;
    unsigned n_thread{ 1 };

  public:
    //!
    //!  \brief Constructs a fake thread pool with a specified number of
    //!  threads.
    //!
    //!  This constructor initializes the thread pool with a specified number
    //!  of threads, although it does not actually create any threads.
    //!  The specified number is ignored.
    //!
    explicit ThreadPool0( unsigned n ) : ThreadPoolBase(), n_thread{ n } {}

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
    void
    exec( FUN && fun ) override
    {
      fun();
    }

    //!
    //! \brief No-op method to wait for all tasks to complete.
    //!
    //! This method does not perform any actions, as there are no
    //! threads or tasks to wait for.
    //!
    void
    wait() override
    {
    }

    //!
    //! \brief Returns the number of threads in the pool.
    //!
    //! This method returns the number of threads in the pool,
    //! which is always 1 for this fake thread pool.
    //!
    //! \return Always returns 1.
    //!
    unsigned
    thread_count() const override
    {
      return n_thread;
    }

    //!
    //! \brief Returns the name of the thread pool.
    //!
    //! This method returns a string describing the type of thread pool.
    //!
    //! \return A string indicating that this is a "ThreadPool0 (fake thread)".
    //!
    static char const *
    Name()
    {
      return "ThreadPool0 (fake thread)";
    }

    char const *
    name() const override
    {
      return Name();
    }
  };

  /*! @} */

}  // namespace Utils

//
// eof: ThreadPool0.hxx
//
