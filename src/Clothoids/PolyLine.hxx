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
/// file: PolyLine.hxx
///

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
  class BiarcList;
  class ClothoidCurve;
  class ClothoidList;

  //! Class to manage a collection of straight segment
  class PolyLine : public BaseCurve {
    friend class ClothoidList;
    friend class BiarcList;
  private:
    vector<LineSegment> m_polyline_list;
    vector<real_type>   m_s0;
    real_type           m_xe;
    real_type           m_ye;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable std::mutex                                         m_last_interval_mutex;
    mutable std::map<std::thread::id,std::shared_ptr<integer>> m_last_interval;
    #else
    mutable integer m_last_interval{0};
    #endif

    mutable bool       m_aabb_done{false};
    mutable AABB_TREE  m_aabb_tree;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable std::mutex m_aabb_mutex;
    #endif

    void
    reset_last_interval() {
      #ifdef CLOTHOIDS_USE_THREADS
      std::unique_lock<std::mutex> lock(m_last_interval_mutex);
      auto id = std::this_thread::get_id();
      auto it = m_last_interval.find(id);
      if ( it == m_last_interval.end() ) it = m_last_interval.insert( {id,std::make_shared<integer>()} ).first;
      integer & last_interval{ *it->second.get() };
      #else
      integer & last_interval = m_last_interval;
      #endif
      last_interval = 0;
    }

  public:

    PolyLine() = delete;

    //explicit
    PolyLine( string const & name ) : BaseCurve( name )
    { this->reset_last_interval(); }

    void setup( GenericContainer const & gc ) override;

    void init();

    void copy( PolyLine const & l );

    //explicit
    PolyLine( PolyLine const & PL ) : BaseCurve( PL.name() )
    { this->reset_last_interval(); this->copy(PL); }

    integer find_at_s( real_type & s ) const;

    explicit PolyLine( LineSegment const & LS );
    explicit PolyLine( CircleArc const & C, real_type tol );
    explicit PolyLine( Biarc const & B, real_type tol );
    explicit PolyLine( ClothoidCurve const & B, real_type tol );
    explicit PolyLine( ClothoidList const & B, real_type tol );
    explicit PolyLine( BaseCurve const * pC );

    CurveType type() const override { return CurveType::POLYLINE; }

    PolyLine const & operator = ( PolyLine const & s )
    { this->copy(s); return *this; }

    LineSegment const &
    getSegment( integer n ) const;

    integer
    num_segments() const
    { return integer(m_polyline_list.size()); }

    integer
    numPoints() const
    { return integer(m_s0.size()); }

    void polygon( real_type x[], real_type y[] ) const;
    void init( real_type x0, real_type y0 );
    void push_back( real_type x, real_type y );
    void push_back( LineSegment const & C );
    void push_back( CircleArc const & C, real_type tol );
    void push_back( Biarc const & C, real_type tol );
    void push_back( ClothoidCurve const & C, real_type tol );
    void push_back( ClothoidList const & L, real_type tol );

    void
    build(
      integer         npts,
      real_type const x[],
      real_type const y[]
    );

    void build( LineSegment const & L );
    void build( CircleArc const & C, real_type tol );
    void build( Biarc const & B, real_type tol );
    void build( ClothoidCurve const & C, real_type tol );
    void build( ClothoidList const & CL, real_type tol );

    void build( CircleArc const & C );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( BiarcList const & );
    void build( PolyLine const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    void
    bbox(
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
    ) const override;

    void
    bbox_ISO(
      real_type   /* offs */,
      real_type & /* xmin */,
      real_type & /* ymin */,
      real_type & /* xmax */,
      real_type & /* ymax */
    ) const override;

    /*\
     |  _    _   _____    _                _
     | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
     | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
     | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
     |                               |___/
    \*/

    void
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override;

    void
    bb_triangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override;

    void
    bb_triangles_SAE(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      this->bb_triangles_ISO( -offs, tvec, max_angle, max_size, icurve );
    }

    real_type
    length() const override
    { return m_s0.back(); }

    real_type
    length_ISO( real_type ) const override;

    real_type
    x_begin() const override
    { return m_polyline_list.front().x_begin(); }

    real_type
    y_begin() const override
    { return m_polyline_list.front().y_begin(); }

    real_type
    x_end() const override
    { return m_polyline_list.back().x_end(); }

    real_type
    y_end() const override
    { return m_polyline_list.back().y_end(); }

    real_type
    X( real_type s ) const override {
      integer idx = this->find_at_s( s );
      real_type ss = m_s0[idx];
      return m_polyline_list[size_t(idx)].X(s-ss);
    }

    real_type
    X_D( real_type s ) const override {
      integer idx = this->find_at_s( s );
      return m_polyline_list.at(size_t(idx)).m_c0;
    }

    real_type
    X_DD( real_type ) const override
    { return 0; }

    real_type
    X_DDD( real_type ) const override
    { return 0; }

    real_type
    Y( real_type s ) const override {
      integer idx = this->find_at_s( s );
      real_type ss = m_s0[idx];
      return m_polyline_list[size_t(idx)].Y(s-ss);
    }

    real_type
    Y_D( real_type s ) const override {
      integer idx = this->find_at_s( s );
      return m_polyline_list[size_t(idx)].m_s0;
    }

    real_type
    Y_DD( real_type ) const override
    { return 0; }

    real_type
    Y_DDD( real_type ) const override
    { return 0; }

    real_type theta    ( real_type s ) const override;
    real_type theta_D  ( real_type s ) const override;
    real_type theta_DD ( real_type s ) const override;
    real_type theta_DDD( real_type s ) const override;

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const override {
      integer idx = this->find_at_s( s );
      real_type ss = m_s0[idx];
      m_polyline_list[size_t(idx)].eval( s-ss, x, y );
    }

    void
    eval_D(
      real_type   s,
      real_type & x_D,
      real_type & y_D
    ) const override {
      integer idx = this->find_at_s( s );
      real_type ss = m_s0[idx];
      m_polyline_list[size_t(idx)].eval_D( s-ss, x_D, y_D );
    }

    void
    eval_DD(
      real_type,
      real_type & x_DD,
      real_type & y_DD
    ) const override
    { x_DD = y_DD = 0; }

    void
    eval_DDD(
      real_type,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override
    { x_DDD = y_DDD = 0; }

    // ---

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const override {
      integer idx{ this->find_at_s( s ) };
      real_type ss{ m_s0[idx] };
      m_polyline_list[size_t(idx)].eval_ISO( s-ss, offs, x, y );
    }

    void
    eval_ISO_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const override {
      integer idx{ this->find_at_s( s ) };
      real_type ss{ m_s0[idx] };
      m_polyline_list[size_t(idx)].eval_ISO_D( s-ss, offs, x_D, y_D );
    }

    void
    eval_ISO_DD(
      real_type,
      real_type,
      real_type & x_DD,
      real_type & y_DD
    ) const override
    { x_DD = y_DD = 0; }

    void
    eval_ISO_DDD(
      real_type,
      real_type,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override
    { x_DDD = y_DDD = 0; }

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void
    translate( real_type tx, real_type ty ) override {
      for ( auto & il : m_polyline_list ) il.translate( tx, ty );
    }

    void
    rotate(
      real_type angle,
      real_type cx,
      real_type cy
    ) override {
      for ( auto & il : m_polyline_list ) il.rotate( angle, cx, cy );
    }

    void reverse() override;

    void scale( real_type sc ) override;

    void change_origin( real_type newx0, real_type newy0 ) override;

    void trim( real_type s_begin, real_type s_end ) override;

    void trim( real_type s_begin, real_type s_end, PolyLine & newPL ) const;

    //!
    //! Compute the point at minimum distance from a point `[x,y]` and the line segment
    //!
    //! \param[in]  x   \f$x\f$-coordinate
    //! \param[in]  y   \f$y\f$-coordinate
    //! \param[out] X   \f$x\f$-coordinate of the closest point
    //! \param[out] Y   \f$y\f$-coordinate of the closest point
    //! \param[out] S   s-param of the closest point
    //! \param[out] T   t-param of the closest point
    //! \param[out] DST the distance point-segment
    //! \return the distance point-segment
    //!
    integer
    closest_point_ISO(
      real_type   x,
      real_type   y,
      real_type & X,
      real_type & Y,
      real_type & S,
      real_type & T,
      real_type & DST
    ) const override;

    integer
    closest_point_ISO(
      real_type   /* x    */,
      real_type   /* y    */,
      real_type   /* offs */,
      real_type & /* X    */,
      real_type & /* Y    */,
      real_type & /* S    */,
      real_type & /* T    */,
      real_type & /* DST  */
    ) const override;

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    bool
    collision( PolyLine const & C ) const;

    bool
    collision_ISO(
      real_type        offs,
      PolyLine const & CL,
      real_type        offs_CL
    ) const;

    bool
    collision( BaseCurve const * pC ) const override;

    bool
    collision_ISO(
      real_type         offs,
      BaseCurve const * pC,
      real_type         offs_C
    ) const override;

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    //!
    //! Intersect PolyLine with another PolyLine
    //!
    //! \param[in]  pl  other PolyLine
    //! \param[out] ss0 list of the paramter of intersection
    //! \param[out] ss1 list of the paramter of intersection of the other PolyLine
    //!
    void
    intersect(
      PolyLine const    & pl,
      vector<real_type> & ss0,
      vector<real_type> & ss1
    ) const;

    //!
    //! Intersect PolyLine with another PolyLine
    //!
    //! \param[in]  pl    other PolyLine
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      PolyLine const & pl,
      IntersectList  & ilist
    ) const;

    //!
    //! Intersect PolyLine with another PolyLine (not yet available)
    //!
    //! \param[in]  offs    PolyLine offset
    //! \param[in]  pl      other PolyLine
    //! \param[in]  offs_pl Other Poliline offset
    //! \param[out] ilist   list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type        offs,
      PolyLine const & pl,
      real_type        offs_pl,
      IntersectList  & ilist
    ) const;

    void
    intersect(
      BaseCurve const * pC,
      IntersectList   & ilist
    ) const override;

    void
    intersect_ISO(
      real_type         offs,
      BaseCurve const * pC,
      real_type         offs_LS,
      IntersectList   & ilist
    ) const override;

    string info() const;

    void
    info( ostream_type & stream ) const override
    { stream << info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, PolyLine const & P );

    void
    build_AABBtree() const;

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "PolyLine_compatibility.hxx"
#endif

  };

}

///
/// eof: PolyLine.hxx
///
