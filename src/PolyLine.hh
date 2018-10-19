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
/// file: PolyLine.hh
///

#ifndef POLY_LINE_HH
#define POLY_LINE_HH

#include "Line.hh"
#include <iterator>     // std::back_inserter
#include <vector>

namespace G2lib {

  /*\
   |  ____       _       _     _
   | |  _ \ ___ | |_   _| |   (_)_ __   ___
   | | |_) / _ \| | | | | |   | | '_ \ / _ \
   | |  __/ (_) | | |_| | |___| | | | |  __/
   | |_|   \___/|_|\__, |_____|_|_| |_|\___|
   |               |___/
  \*/

  class CircleArc;
  class Biarc;
  class ClothoidCurve;
  class ClothoidList;

  class PolyLine {

    std::vector<LineSegment> lvec;
    std::vector<real_type>   s0;
    real_type                xe, ye;

    mutable int_type isegment;
    void search( real_type s ) const;

  public:

    PolyLine()
    : isegment(0)
    {}

    void
    copy( PolyLine const & l ) {
      lvec.clear();
      lvec.reserve(l.lvec.size());
      std::copy( l.lvec.begin(),
                 l.lvec.end(),
                 back_inserter(lvec) );
    }

    PolyLine( PolyLine const & s ) { copy(s); }

    PolyLine const & operator = ( PolyLine const & s )
    { copy(s); return *this; }

    LineSegment const &
    getSegment( int_type n ) const;

    int_type
    numSegment() const
    { return int_type(lvec.size()); }

    int_type
    numPoints() const
    { return int_type(s0.size()); }

    void
    polygon( real_type x[], real_type y[]) const;

    void
    init( real_type x0, real_type y0 );

    void
    push_back( real_type x, real_type y );

    void
    push_back( LineSegment const & C );

    void
    push_back( CircleArc const & C, real_type tol );

    void
    push_back( Biarc const & C, real_type tol );

    void
    push_back( ClothoidCurve const & C, real_type tol );

    void
    push_back( ClothoidList const & L, real_type tol );

    void
    build( real_type const x[],
           real_type const y[],
           int_type npts );

    void
    build( LineSegment const & L );

    void
    build( CircleArc const & C, real_type tol );

    void
    build( Biarc const & C, real_type tol );

    void
    build( ClothoidCurve const & C, real_type tol );

    void
    build( ClothoidList const & L, real_type tol );

    real_type length() const { return s0.back(); }

    real_type xBegin() const { return lvec.front().xBegin(); }
    real_type yBegin() const { return lvec.front().yBegin(); }
    real_type xEnd()   const { return lvec.back().xEnd(); }
    real_type yEnd()   const { return lvec.back().yEnd(); }

    real_type
    X( real_type s ) const {
      this->search( s );
      return lvec[size_t(isegment)].X( s-s0[size_t(isegment)] );
    }

    real_type
    Y( real_type s ) const {
      this->search( s );
      return lvec[size_t(isegment)].Y(s-s0[size_t(isegment)]);
    }

    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const {
      this->search( s );
      lvec[size_t(isegment)].eval(s-s0[size_t(isegment)],x,y);
    }

    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const {
      this->search( s );
      lvec[size_t(isegment)].eval_D( s-s0[size_t(isegment)], x_D, y_D );
    }

    void
    eval_DD( real_type,
             real_type & x_DD,
             real_type & y_DD ) const
    { x_DD = y_DD = 0; }

    void
    eval_DDD( real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const
    { x_DDD = y_DDD = 0; }

    // ---

    void
    eval( real_type   s,
          real_type   t,
          real_type & x,
          real_type & y ) const {
      this->search( s );
      lvec[size_t(isegment)].eval( s-s0[size_t(isegment)], t, x, y );
    }

    void
    eval_D( real_type   s,
            real_type   t,
            real_type & x_D,
            real_type & y_D ) const {
      this->search( s );
      lvec[size_t(isegment)].eval_D( s-s0[size_t(isegment)], t, x_D, y_D );
    }

    void
    eval_DD( real_type,
             real_type,
             real_type & x_DD,
             real_type & y_DD ) const
    { x_DD = y_DD = 0; }

    void
    eval_DDD( real_type,
              real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const
    { x_DDD = y_DDD = 0; }

    void
    translate( real_type tx, real_type ty ) {
      std::vector<LineSegment>::iterator il;
      for ( il = lvec.begin(); il != lvec.end(); ++il )
        il->translate( tx, ty );
    }

    void
    rotate( real_type angle, real_type cx, real_type cy ) {
      std::vector<LineSegment>::iterator il;
      for ( il = lvec.begin(); il != lvec.end(); ++il )
        il->rotate( angle, cx, cy );
    }

    void reverse();

    /*!
     * \brief compute the point at minimum distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param X x-coordinate of the closest point
     * \param Y y-coordinate of the closest point
     * \param S param of the closest point
     * \return the distance point-segment
    \*/
    real_type
    closestPoint( real_type   x,
                  real_type   y,
                  real_type & X,
                  real_type & Y,
                  real_type & S ) const;

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \param S param at minimum distance
     * \return the distance point-segment
    \*/
    real_type
    distance( real_type   x,
              real_type   y,
              real_type & S ) const {
      real_type X, Y;
      return closestPoint( x, y, X, Y, S );
    }

    /*!
     * \brief compute the distance from a point `[x,y]` and the line segment
     *
     * \param x x-coordinate
     * \param y y-coordinate
     * \return the distance point-segment
    \*/
    real_type
    distance( real_type x, real_type y ) const {
      real_type ss;
      return distance( x, y, ss );
    }

    void
    intersect( PolyLine const         & pl,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2 ) const;

    bool
    intersect( PolyLine const & pl ) const;

    void
    info( ostream_type & stream ) const
    { stream << "PolyLine\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, PolyLine const & P );

  };

}

#endif

///
/// eof: PolyLine.hh
///
