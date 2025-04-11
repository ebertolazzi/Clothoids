/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
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
// file: ThreadPool3.hxx
//

#include "3rd/BS_thread_pool.hpp"

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
  //! \brief A thread pool for concurrent task execution with worker management.
  //!
  //! The `ThreadPool3` class provides an implementation of a thread pool that
  //! allows for efficient concurrent execution of tasks. This class manages a
  //! collection of worker threads that can be dynamically resized, and utilizes
  //! a semaphore-based mechanism for worker synchronization.
  //!
  //! This class extends `ThreadPoolBase` and supports task execution, joining,
  //! and resizing of worker threads.
  //!
  class ThreadPool3 : public ThreadPoolBase {

    BS::light_thread_pool m_pool;

  public:

    //!
    //! \brief Constructs a new ThreadPool5 instance with a specified number of threads.
    //!
    //! \param nthread The number of threads to create in the pool. Defaults to the maximum hardware threads available.
    //!
    ThreadPool3(
      unsigned nthread = std::max(
        unsigned(1),
        unsigned(std::thread::hardware_concurrency()-1)
      )
    ) : m_pool( nthread ) {}

    //!
    //! \brief Destructor for the ThreadPool5 class.
    //!
    //! Ensures all workers are stopped and joined before destruction.
    //!
    virtual ~ThreadPool3() {}

    //!
    //! \brief Executes a task and assigns it to an available worker.
    //!
    //! \param fun The function to be executed as a task.
    //!
    void exec( FUN && fun ) override { m_pool.detach_task( std::move(fun) ); }

    //!
    //! \brief Waits for all tasks to be completed.
    //!
    void wait() override { m_pool.wait(); }

    //!
    //! \brief Gets the current number of threads in the pool.
    //!
    //! \return The number of threads in the pool.
    //!
    unsigned thread_count() const override { return unsigned( m_pool.get_thread_count() ); }

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    static char const * Name() { return "ThreadPool3 [BS]"; }

    char const * name() const override { return Name(); }

  };

  /*! @} */

}

//
// eof: ThreadPool3.hxx
//
