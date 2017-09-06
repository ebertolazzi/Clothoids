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
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Triangle2D.hh
///

#ifndef TRIANGLE2D_HH
#define TRIANGLE2D_HH

#include <cmath>

//! Clothoid computations routine
namespace Triangle2D {

  /*\
   |   _____     _                   _      ____  ____
   |  |_   _| __(_) __ _ _ __   __ _| | ___|___ \|  _ \
   |    | || '__| |/ _` | '_ \ / _` | |/ _ \ __) | | | |
   |    | || |  | | (_| | | | | (_| | |  __// __/| |_| |
   |    |_||_|  |_|\__,_|_| |_|\__, |_|\___|_____|____/
   |                           |___/
  \*/
  //! \brief Class to manage Triangle for BB of clothoid curve

  template <typename T>
  class Triangle2D {
    T p1[2], p2[2], p3[2] ;

  public:

    Triangle2D( ) {
      p1[0] = p1[1] =
      p2[0] = p2[1] =
      p3[0] = p3[1] = 0 ;
    }

    Triangle2D( T x1, T y1,
                T x2, T y2,
                T x3, T y3 ) {
      p1[0] = x1; p1[1] = y1;
      p2[0] = x2; p2[1] = y2;
      p3[0] = x3; p3[1] = y3;
    }

    Triangle2D( T const _p1[2],
                T const _p2[2],
                T const _p3[2] ) {
      p1[0] = _p1[0] ; p1[1] = _p1[1] ;
      p2[0] = _p2[0] ; p2[1] = _p2[1] ;
      p3[0] = _p3[0] ; p3[1] = _p3[1] ;
    }

    ~Triangle2D() {}
    
    void
    setup( T const _p1[2],
           T const _p2[2],
           T const _p3[2] ) {
      p1[0] = _p1[0] ; p1[1] = _p1[1] ;
      p2[0] = _p2[0] ; p2[1] = _p2[1] ;
      p3[0] = _p3[0] ; p3[1] = _p3[1] ;
    }
    
    T x1() const { return p1[0] ; }
    T y1() const { return p1[1] ; }
    T x2() const { return p2[0] ; }
    T y2() const { return p2[1] ; }
    T x3() const { return p3[0] ; }
    T y3() const { return p3[1] ; }

    T const * P1() const { return p1 ; }
    T const * P2() const { return p2 ; }
    T const * P3() const { return p3 ; }

    bool intersect( Triangle2D<T> const & ) const ;
    bool overlap( Triangle2D<T> const & ) const ;

  };

  extern template class Triangle2D<float> ;
  extern template class Triangle2D<double> ;

}

#endif

///
/// eof: Triangle2D.hh
///
