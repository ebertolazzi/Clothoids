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
#include <vector>
#include "G2lib.hh"

//! Clothoid computations routine
namespace G2lib {

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

    void
    maxmin3( T const & a, T const & b, T const & c, T & vmin, T & vmax) const ;

    class AABBTree {
      int numChildren ;
      T   xmin, ymin, xmax, ymax ;
      union {
        Triangle2D<T> const * pTriangle;
        AABBTree      const * pChildren[2];
      } data ;

      AABBTree( typename std::vector<Triangle2D<T> const *>::iterator & begin,
                typename std::vector<Triangle2D<T> const *>::iterator & end ) ;

    public:

      AABBTree( std::vector<Triangle2D<T> > & triangles );
      AABBTree( Triangle2D<T> const & triangle );
      AABBTree( AABBTree const * pTree, Triangle2D<T> const & triangle );
      AABBTree( AABBTree const * pTreeL, AABBTree const * pTreeR );

      ~AABBTree();

      void
      bbox( T & _xmin, T & _ymin, T & _xmax, T & _ymax ) const
      { _xmin = xmin ; _ymin = ymin ; _xmax = xmax ; _ymax = ymax ; }

      bool overlap( Triangle2D<T> const & triangle ) const ;
      bool overlap( AABBTree const * pTree ) const ;

      Triangle2D<T> const & getTriangle() const {
        G2LIB_ASSERT( numChildren == 0,
                      "Triangle2D::AABBTree::getTriangle() not a leaf" ) ;
        return *data.pTriangle ;
      }
    };

  public:

    Triangle2D( Triangle2D<T> const & t ) {
      p1[0] = t.p1[0] ; p1[1] = t.p1[1] ;
      p2[0] = t.p2[0] ; p2[1] = t.p2[1] ;
      p3[0] = t.p3[0] ; p3[1] = t.p3[1] ;
    }

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

    void
    bbox( T & xmin, T & ymin, T & xmax, T & ymax ) const {
      maxmin3( p1[0], p2[0], p3[0], xmin, xmax ) ;
      maxmin3( p1[1], p2[1], p3[1], ymin, ymax ) ;
    }

    T baricenterX() const { return (p1[0]+p2[0]+p3[0])/3 ; }
    T baricenterY() const { return (p1[1]+p2[1]+p3[1])/3 ; }

    T const * P1() const { return p1 ; }
    T const * P2() const { return p2 ; }
    T const * P3() const { return p3 ; }

    bool intersect( Triangle2D<T> const & ) const ;
    bool overlap( Triangle2D<T> const & ) const ;

  };
  
  // explicit instantiation declaration to suppress warnings
  #ifdef G2LIB_USE_CXX11

  #ifdef __GCC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class Triangle2D<float> ;
  extern template class Triangle2D<double> ;
  
  #endif

}

#endif

///
/// eof: Triangle2D.hh
///
