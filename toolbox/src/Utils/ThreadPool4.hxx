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
// file: ThreadPool4.hxx
//

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
  //! \brief A thread pool for concurrent task execution with worker management.
  //!
  //! The `ThreadPool4` class provides an implementation of a thread pool that
  //! allows for efficient concurrent execution of tasks. This class manages a
  //! collection of worker threads that can be dynamically resized, and utilizes
  //! a semaphore-based mechanism for worker synchronization.
  //!
  //! This class extends `ThreadPoolBase` and supports task execution, joining,
  //! and resizing of worker threads.
  //!
  class ThreadPool4 : public ThreadPool3
  {
  public:
    //!
    //! \brief Constructs a new ThreadPool5 instance with a specified number of
    //! threads.
    //!
    //! \param nthread The number of threads to create in the pool. Defaults to
    //! the maximum hardware threads available.
    //!
    using ThreadPool3::ThreadPool3;

    //!
    //! \brief Destructor for the ThreadPool5 class.
    //!
    //! Ensures all workers are stopped and joined before destruction.
    //!
    virtual ~ThreadPool4() override {}

    //!
    //! \brief Gets the name of the thread pool implementation.
    //!
    //! \return A constant character pointer to the name of the thread pool.
    //!
    static char const * Name() { return "ThreadPool4 [BS compat]"; }

    char const * name() const override { return Name(); }
  };

  /*! @} */

}  // namespace Utils

//
// eof: ThreadPool4.hxx
//
