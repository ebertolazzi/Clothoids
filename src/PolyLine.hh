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
#include "AABBtree.hh"

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

  class PolyLine : public BaseCurve {
  private:
    vector<LineSegment> polylineList;
    vector<real_type>   s0;
    real_type           xe, ye;

    mutable int_type isegment;
    void search( real_type s ) const;

    mutable bool     aabb_done;
    mutable AABBtree aabb;

  public:

    explicit
    PolyLine()
    : BaseCurve(G2LIB_POLYLINE)
    , isegment(0)
    , aabb_done(false)
    {}

    void
    copy( PolyLine const & l ) {
      polylineList.clear();
      polylineList.reserve(l.polylineList.size());
      std::copy( l.polylineList.begin(),
                 l.polylineList.end(),
                 back_inserter(polylineList) );
      aabb_done = false;
    }

    explicit
    PolyLine( PolyLine const & PL )
    : BaseCurve(G2LIB_POLYLINE)
    { copy(PL); }

    explicit
    PolyLine( BaseCurve const & C );

    PolyLine const & operator = ( PolyLine const & s )
    { copy(s); return *this; }

    LineSegment const &
    getSegment( int_type n ) const;

    int_type
    numSegment() const
    { return int_type(polylineList.size()); }

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

    virtual
    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const G2LIB_OVERRIDE;

    virtual
    void
    bbox( real_type   /* offs */,
          real_type & /* xmin */,
          real_type & /* ymin */,
          real_type & /* xmax */,
          real_type & /* ymax */ ) const G2LIB_OVERRIDE {
      G2LIB_ASSERT( false, "PolyLine::bbox( offs ... ) not available!");
    }

    virtual
    real_type
    length() const G2LIB_OVERRIDE
    { return s0.back(); }

    virtual
    real_type
    length( real_type ) const G2LIB_OVERRIDE {
      G2LIB_ASSERT( false, "PolyLine::length( offs ) not available!");
      return 0;
    }

    virtual
    real_type
    xBegin() const G2LIB_OVERRIDE
    { return polylineList.front().xBegin(); }

    virtual
    real_type
    yBegin() const G2LIB_OVERRIDE
    { return polylineList.front().yBegin(); }

    virtual
    real_type
    xEnd() const G2LIB_OVERRIDE
    { return polylineList.back().xEnd(); }

    virtual
    real_type
    yEnd() const G2LIB_OVERRIDE
    { return polylineList.back().yEnd(); }

    virtual
    real_type
    X( real_type s ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      return polylineList[size_t(isegment)].X(s-ss);
    }

    virtual
    real_type
    X_D( real_type s ) const G2LIB_OVERRIDE {
      this->search( s );
      return polylineList[size_t(isegment)].c0;
    }

    virtual
    real_type
    X_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    X_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y( real_type s ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      return polylineList[size_t(isegment)].Y(s-ss);
    }

    virtual
    real_type
    Y_D( real_type s ) const G2LIB_OVERRIDE {
      this->search( s );
      return polylineList[size_t(isegment)].s0;
    }

    virtual
    real_type
    Y_DD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    Y_DDD( real_type ) const G2LIB_OVERRIDE
    { return 0; }

    virtual
    real_type
    theta( real_type s ) const G2LIB_OVERRIDE;

    virtual
    real_type
    theta_D( real_type s ) const G2LIB_OVERRIDE;

    virtual
    real_type
    theta_DD( real_type s ) const G2LIB_OVERRIDE;

    virtual
    real_type
    theta_DDD( real_type s ) const G2LIB_OVERRIDE;

    virtual
    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      polylineList[size_t(isegment)].eval( s-ss, x, y );
    }

    virtual
    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      polylineList[size_t(isegment)].eval_D( s-ss, x_D, y_D );
    }

    virtual
    void
    eval_DD( real_type,
             real_type & x_DD,
             real_type & y_DD ) const G2LIB_OVERRIDE
    { x_DD = y_DD = 0; }

    virtual
    void
    eval_DDD( real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const G2LIB_OVERRIDE
    { x_DDD = y_DDD = 0; }

    // ---

    virtual
    void
    eval( real_type   s,
          real_type   t,
          real_type & x,
          real_type & y ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      polylineList[size_t(isegment)].eval( s-ss, t, x, y );
    }

    virtual
    void
    eval_D( real_type   s,
            real_type   t,
            real_type & x_D,
            real_type & y_D ) const G2LIB_OVERRIDE {
      this->search( s ); real_type ss = s0[size_t(isegment)];
      polylineList[size_t(isegment)].eval_D( s-ss, t, x_D, y_D );
    }

    virtual
    void
    eval_DD( real_type,
             real_type,
             real_type & x_DD,
             real_type & y_DD ) const G2LIB_OVERRIDE
    { x_DD = y_DD = 0; }

    virtual
    void
    eval_DDD( real_type,
              real_type,
              real_type & x_DDD,
              real_type & y_DDD ) const G2LIB_OVERRIDE
    { x_DDD = y_DDD = 0; }

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    virtual
    void
    translate( real_type tx, real_type ty ) G2LIB_OVERRIDE {
      std::vector<LineSegment>::iterator il;
      for ( il = polylineList.begin(); il != polylineList.end(); ++il )
        il->translate( tx, ty );
    }

    virtual
    void
    rotate( real_type angle,
            real_type cx,
            real_type cy ) G2LIB_OVERRIDE {
      std::vector<LineSegment>::iterator il;
      for ( il = polylineList.begin(); il != polylineList.end(); ++il )
        il->rotate( angle, cx, cy );
    }

    virtual
    void
    reverse() G2LIB_OVERRIDE;

    virtual
    void
    scale( real_type sc ) G2LIB_OVERRIDE;

    virtual
    void
    changeOrigin( real_type newx0, real_type newy0 ) G2LIB_OVERRIDE;

    virtual
    void
    trim( real_type s_begin, real_type s_end ) G2LIB_OVERRIDE;

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
    virtual
    int_type
    closestPoint( real_type   x,
                  real_type   y,
                  real_type & X,
                  real_type & Y,
                  real_type & S,
                  real_type & T,
                  real_type & DST ) const G2LIB_OVERRIDE;

    virtual
    int_type
    closestPoint( real_type   /* x    */,
                  real_type   /* y    */,
                  real_type   /* offs */,
                  real_type & /* X    */,
                  real_type & /* Y    */,
                  real_type & /* S    */,
                  real_type & /* T    */,
                  real_type & /* DST */ ) const G2LIB_OVERRIDE {
      G2LIB_ASSERT( false, "PolyLine::closestPoint( ... offs ... ) not available!");
    }

    void
    intersect( PolyLine const         & pl,
               std::vector<real_type> & s1,
               std::vector<real_type> & s2 ) const;

    bool
    intersect( PolyLine const & pl ) const;

    void
    intersect( PolyLine const & pl,
               IntersectList  & ilist,
               bool             swap_s_vals ) {
      std::vector<real_type> s1, s2;
      this->intersect( pl, s1, s2 );
      ilist.reserve( ilist.size() + s1.size() );
      for ( size_t i=0; i < s1.size(); ++i ) {
        real_type ss1 = s1[i];
        real_type ss2 = s2[i];
        if ( swap_s_vals ) std::swap( ss1, ss2 );
        ilist.push_back( Ipair(ss1,ss2) );
      }
    }

    void
    intersect( real_type        /* offs        */,
               PolyLine const & /* pl          */,
               real_type        /* offs_pl     */,
               IntersectList  & /* ilist       */,
               bool             /* swap_s_vals */ ) {
      G2LIB_ASSERT( false, "PolyLine::intersect( offs ... ) not available!");
    }

    virtual
    void
    info( ostream_type & stream ) const G2LIB_OVERRIDE
    { stream << "PolyLine\n" << *this << '\n'; }

    friend
    ostream_type &
    operator << ( ostream_type & stream, PolyLine const & P );

    void
    build_AABBtree( AABBtree & aabb ) const;

    void
    build_AABBtree() const {
      if ( !aabb_done ) {
        this->build_AABBtree( aabb );
        aabb_done = true;
      }
    }

  };

}

#endif

///
/// eof: PolyLine.hh
///
