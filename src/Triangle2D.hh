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

  class Triangle2D {
    real_type p1[2], p2[2], p3[2];

    void
    maxmin3( real_type   a,
             real_type   b,
             real_type   c,
             real_type & vmin,
             real_type & vmax) const;

    class AABBTree {
      int       numChildren;
      real_type xmin, ymin, xmax, ymax;
      union {
        Triangle2D const * pTriangle;
        AABBTree   const * pChildren[2];
      } data;

      AABBTree( std::vector<Triangle2D const *>::iterator & begin,
                std::vector<Triangle2D const *>::iterator & end );

    public:

      AABBTree( std::vector<Triangle2D > & triangles );
      AABBTree( Triangle2D const & triangle );
      AABBTree( AABBTree const * pTree, Triangle2D const & triangle );
      AABBTree( AABBTree const * pTreeL, AABBTree const * pTreeR );

      ~AABBTree();

      void
      bbox( real_type & _xmin,
            real_type & _ymin,
            real_type & _xmax,
            real_type & _ymax ) const
      { _xmin = xmin; _ymin = ymin; _xmax = xmax; _ymax = ymax; }

      bool overlap( Triangle2D const & triangle ) const;
      bool overlap( AABBTree const * pTree ) const;

      Triangle2D const & getTriangle() const {
        G2LIB_ASSERT( numChildren == 0,
                      "Triangle2D::AABBTree::getTriangle() not a leaf" );
        return *data.pTriangle;
      }
    };

  public:

    Triangle2D( Triangle2D const & t ) {
      p1[0] = t.p1[0]; p1[1] = t.p1[1];
      p2[0] = t.p2[0]; p2[1] = t.p2[1];
      p3[0] = t.p3[0]; p3[1] = t.p3[1];
    }

    Triangle2D( ) {
      p1[0] = p1[1] =
      p2[0] = p2[1] =
      p3[0] = p3[1] = 0;
    }

    Triangle2D( real_type x1, real_type y1,
                real_type x2, real_type y2,
                real_type x3, real_type y3 ) {
      p1[0] = x1; p1[1] = y1;
      p2[0] = x2; p2[1] = y2;
      p3[0] = x3; p3[1] = y3;
    }

    Triangle2D( real_type const _p1[2],
                real_type const _p2[2],
                real_type const _p3[2] ) {
      p1[0] = _p1[0]; p1[1] = _p1[1];
      p2[0] = _p2[0]; p2[1] = _p2[1];
      p3[0] = _p3[0]; p3[1] = _p3[1];
    }

    ~Triangle2D() {}

    void
    build( real_type const _p1[2],
           real_type const _p2[2],
           real_type const _p3[2] ) {
      p1[0] = _p1[0]; p1[1] = _p1[1];
      p2[0] = _p2[0]; p2[1] = _p2[1];
      p3[0] = _p3[0]; p3[1] = _p3[1];
    }

    void
    build( real_type x1, real_type y1,
           real_type x2, real_type y2,
           real_type x3, real_type y3 ) {
      p1[0] = x1; p1[1] = y1;
      p2[0] = x2; p2[1] = y2;
      p3[0] = x3; p3[1] = y3;
    }

    real_type x1() const { return p1[0]; }
    real_type y1() const { return p1[1]; }

    real_type x2() const { return p2[0]; }
    real_type y2() const { return p2[1]; }

    real_type x3() const { return p3[0]; }
    real_type y3() const { return p3[1]; }

    void
    translate( real_type tx, real_type ty ) {
      p1[0] += tx; p2[0] += tx; p3[0] += tx;
      p1[1] += ty; p2[1] += ty; p3[1] += ty;
    }

    void
    rotate( real_type angle, real_type cx, real_type cy );

    void
    scale( real_type sc ) {
      p1[0] *= sc; p1[1] *= sc;
      p2[0] *= sc; p2[1] *= sc;
      p3[0] *= sc; p3[1] *= sc;
    }

    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const {
      maxmin3( p1[0], p2[0], p3[0], xmin, xmax );
      maxmin3( p1[1], p2[1], p3[1], ymin, ymax );
    }

    real_type baricenterX() const { return (p1[0]+p2[0]+p3[0])/3; }
    real_type baricenterY() const { return (p1[1]+p2[1]+p3[1])/3; }

    real_type const * P1() const { return p1; }
    real_type const * P2() const { return p2; }
    real_type const * P3() const { return p3; }

    bool intersect( Triangle2D const & ) const;
    bool overlap( Triangle2D const & ) const;

    /*!
    //  return +1 = CounterClockwise
    //  return -1 = Clockwise
    //  return  0 = flat
    */
    int_type
    isCounterClockwise() const {
      return G2lib::isCounterClockwise( p1, p2, p3 );
    }

    /*!
    //  return +1 = inside
    //  return -1 = outside
    //  return  0 = on the border
    */
    int_type
    isInside( real_type x, real_type y ) const {
      real_type const pt[2] = {x,y};
      return isPointInTriangle( pt, p1, p2, p3 );
    }

    int_type
    isInside( real_type const pt[2] ) const {
      return isPointInTriangle( pt, p1, p2, p3 );
    }

    real_type
    distMin( real_type x, real_type y ) const;

    real_type
    distMax( real_type x, real_type y ) const;

    void
    info( ostream_type & stream ) const
    { stream << "Triangle2D\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, Triangle2D const & c );

  };

}

#endif

///
/// eof: Triangle2D.hh
///
