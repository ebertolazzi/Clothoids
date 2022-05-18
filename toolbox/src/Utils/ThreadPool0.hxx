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
/// file: ThreadPool0.hxx
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

  //! Fake thread pool!
  class ThreadPool0 : public ThreadPoolBase {
    typedef std::function<void()> Func;

  public:

    explicit
    ThreadPool0( unsigned )
    : ThreadPoolBase()
    {}

    virtual ~ThreadPool0() = default;

    void         exec( Func && fun )  override { fun(); }
    void         wait()               override { }
    void         join()               override { }
    unsigned     thread_count() const override { return 1; }
    void         resize( unsigned )   override { }
    char const * name() const         override { return "ThreadPool0 (fake thread)"; }
    unsigned     size() const                  { return 1; }
  };
}

///
/// eof: ThreadPool0.hxx
///
