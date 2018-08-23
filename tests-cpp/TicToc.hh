/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 \*--------------------------------------------------------------------------*/

#ifndef TIC_TOC_HH
#define TIC_TOC_HH

#include <chrono>
#include <thread>

#define TIC_TOC_USE_HRC

class TicToc {

  #ifdef TIC_TOC_USE_HRC
  std::chrono::high_resolution_clock::time_point start_time;
  std::chrono::high_resolution_clock::time_point stop_time;
  #else
  std::chrono::system_clock::time_point start_time;
  std::chrono::system_clock::time_point stop_time;
  #endif
  
  std::chrono::microseconds elapsedPartial;
  std::chrono::microseconds elapsedTotal;

  TicToc( TicToc const & );
  TicToc const & operator = ( TicToc const & ) const;

public:

  TicToc()
  : elapsedPartial(0)
  , elapsedTotal(0) { tic(); }

  ~TicToc() {}

  void reset() {
    elapsedTotal   = std::chrono::microseconds(0);
    elapsedPartial = std::chrono::microseconds(0);
  }

  void
  tic()
  #ifdef TIC_TOC_USE_HRC
  { start_time = std::chrono::high_resolution_clock::now(); }
  #else
  { start_time = std::chrono::system_clock::now(); }
  #endif

  void
  toc() {
    #ifdef TIC_TOC_USE_HRC
    stop_time = std::chrono::high_resolution_clock::now();
    #else
    stop_time = std::chrono::system_clock::now();
    #endif
    elapsedPartial = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
    elapsedTotal  += elapsedPartial;
  }

  double totalElapsedSeconds()      const { return 1e-6*elapsedTotal.count(); }
  double totalElapsedMilliseconds() const { return 1e-3*elapsedTotal.count(); }

  double elapsedSeconds()      const { return 1e-6*elapsedPartial.count(); }
  double elapsedMilliseconds() const { return 1e-3*elapsedPartial.count(); }
  
  void sleep_for_seconds( unsigned s ) {
    std::this_thread::sleep_for(std::chrono::seconds(s));
  }

  void sleep_for_milliseconds( unsigned ms ) {
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
  }

};

#endif
