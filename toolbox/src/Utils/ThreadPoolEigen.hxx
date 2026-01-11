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

#pragma once

#include "3rd/Eigen/ThreadPool"

namespace Utils
{

  class ThreadPoolEigen : public ThreadPoolBase
  {
  public:
    using FUN = ThreadPoolBase::FUN;

  private:
    Eigen::ThreadPool   m_pool;
    std::atomic<size_t> m_pending;
    unsigned            m_nthread;

  public:
    explicit ThreadPoolEigen( unsigned nthread ) : m_pool( int( nthread ) ), m_pending( 0 ), m_nthread( nthread ) {}

    virtual ~ThreadPoolEigen() override { wait(); }

    static char const * Name() { return "ThreadPoolEigen"; }

    char const * name() const override { return Name(); }

    unsigned thread_count() const override { return m_nthread; }

    void exec( FUN && fun ) override
    {
      m_pending.fetch_add( 1, std::memory_order_relaxed );

      m_pool.Schedule(
        [this, f = std::move( fun )]() mutable
        {
          try
          {
            f();
          }
          catch ( ... )
          {
            // opzionale: logging
          }
          m_pending.fetch_sub( 1, std::memory_order_release );
        } );
    }

    void wait() override
    {
      while ( m_pending.load( std::memory_order_acquire ) != 0 ) { std::this_thread::yield(); }
    }
  };

}  // namespace Utils
