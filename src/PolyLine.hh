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

    mutable int_type lastSegment;

    int_type
    search( real_type s ) const {
      int_type npts = int_type(s0.size());
      G2LIB_ASSERT( npts > 1,
                    "PolyLine::search(" << s << ") empty PolyLine" );
      real_type sl = s0.front();
      real_type sr = s0.back();
      G2LIB_ASSERT( s >= sl && s <= sr,
                    "PolyLine::search( " << s <<
                    " ) out of range: [" << sl << ", " << sr << "]" );
      updateInterval( lastSegment, s, &s0.front(), npts );
      return lastSegment;
    }

  public:

    PolyLine()
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
    getSegment( int_type n ) const
    { return lvec[n]; }

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
           int_type npts ) {
      init( x[0], y[0] );
      for ( int_type k = 1; k < npts; ++k )
        push_back( x[k], y[k] );
    }

    void
    build( LineSegment const & C );

    void
    build( CircleArc const & C, real_type tol );

    void
    build( Biarc const & C, real_type tol );

    void
    build( ClothoidCurve const & C, real_type tol );

    void
    build( ClothoidList const & L, real_type tol );

    int_type
    numPoints() const
    { return int_type(s0.size()); }

    void
    polygon( real_type x[], real_type y[]) const {
      int_type n = int_type(lvec.size());
      for ( int_type k = 0; k < n; ++k ) {
        x[k] = lvec[k].xBegin();
        y[k] = lvec[k].yBegin();
      }
      x[n] = lvec[n-1].xEnd();
      y[n] = lvec[n-1].yEnd();
    }

    real_type length() const { return s0.back(); }

    real_type xBegin() const { return lvec.front().xBegin(); }
    real_type yBegin() const { return lvec.front().yBegin(); }
    real_type xEnd()   const { return lvec.back().xEnd(); }
    real_type yEnd()   const { return lvec.back().yEnd(); }

    real_type
    X( real_type s ) const {
      int_type ns = search( s );
      return lvec[ns].X(s-s0[ns]);
    }

    real_type
    Y( real_type s ) const {
      int_type ns = search( s );
      return lvec[ns].Y(s-s0[ns]);
    }

    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const {
      int_type ns = search( s );
      lvec[ns].eval(s-s0[ns],x,y);
    }

    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const {
      int_type ns = search( s );
      lvec[ns].eval_D(s-s0[ns],x_D,y_D);
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
      int_type ns = search( s );
      lvec[ns].eval(s-s0[ns],t,x,y);
    }

    void
    eval_D( real_type   s,
            real_type   t,
            real_type & x_D,
            real_type & y_D ) const {
      int_type ns = search( s );
      lvec[ns].eval_D(s-s0[ns],t,x_D,y_D);
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

  };

}

#endif

///
/// eof: PolyLine.hh
///
