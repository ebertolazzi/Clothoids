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
/// file: BiarcList.hh
///

namespace G2lib {

  using std::vector;

  class ClothoidList;

  /*\
   |  ____  _                _     _     _
   | | __ )(_) __ _ _ __ ___| |   (_)___| |_
   | |  _ \| |/ _` | '__/ __| |   | / __| __|
   | | |_) | | (_| | | | (__| |___| \__ \ |_
   | |____/|_|\__,_|_|  \___|_____|_|___/\__|
  \*/
  //!
  //! Class to manage a list of biarc Curve (not necessarily \f$ G^2 \f$ or \f$ G^1 \f$ connected)
  //!
  //! @html_image{biarc_list.png,width=60%}
  //!
  class BiarcList : public BaseCurve {

    friend class ClothoidList;

    vector<real_type> m_s0;
    vector<Biarc>     m_biarc_list;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable std::mutex                                         m_last_interval_mutex;
    mutable std::map<std::thread::id,std::shared_ptr<integer>> m_last_interval;
    #else
    mutable integer m_last_interval{0};
    #endif

    mutable bool               m_aabb_done{false};
    mutable AABB_TREE          m_aabb_tree;
    mutable real_type          m_aabb_offs{real_type(0)};
    mutable real_type          m_aabb_max_angle{real_type(0)};
    mutable real_type          m_aabb_max_size{real_type(0)};
    mutable vector<Triangle2D> m_aabb_triangles;

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

    integer
    closest_point_internal(
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & dst
    ) const;

  public:

    #include "BaseCurve_using.hxx"

    //!
    //! Build an empty biarc spline.
    //!
    BiarcList( string const & name ) : BaseCurve( name )
    { this->reset_last_interval(); }

    ~BiarcList() override {
      m_s0.clear();
      m_biarc_list.clear();
      m_aabb_triangles.clear();
    }

    void setup( GenericContainer const & gc ) override;

    //!
    //! Build a copy of another biarc spline.
    //!
    BiarcList( BiarcList const & s ) : BiarcList( s.name() )
    { this->copy(s); }

    //!
    //! Empty the the biarc list.
    //!
    void init();

    //!
    //! Reserve memory for `n` biarcs.
    //!
    void reserve( integer n );

    //!
    //! Copy another biarc spline.
    //!
    void copy( BiarcList const & L );

    CurveType type() const override { return CurveType::BIARC_LIST; }

    //!
    //! Copy another biarc spline.
    //!
    BiarcList const & operator = ( BiarcList const & s )
    { this->copy(s); return *this; }

    //!
    //! Build a biarc list from a line segment.
    //!
    explicit BiarcList( LineSegment const & LS );

    //!
    //! Build a biarc list from a single circle arc.
    //!
    explicit BiarcList( CircleArc const & C );

    //!
    //! Build a biarc list from a single biarc.
    //!
    explicit BiarcList( Biarc const & C );

    //!
    //! Build a biarc list from a single polyline.
    //!
    explicit BiarcList( PolyLine const & pl );

    //!
    //! Build a biarc list from another curve.
    //!
    explicit BiarcList( BaseCurve const * pC );

    void build( LineSegment const & );
    void build( CircleArc const & );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( BiarcList const & );
    void build( PolyLine const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    //!
    //! Append a line segment to the biarc list
    //! (transformed to a degenerate biarc).
    //!
    void push_back( LineSegment const & c );

    //!
    //! Append a line circle to the biarc
    //! list (transformed to a degenerate biarc).
    //!
    void push_back( CircleArc const & c );

    //!
    //! Append a biarc to the biarc list.
    //!
    void push_back( Biarc const & c );

    //!
    //! Append a polyline to the biarc list
    //! (transformed to a list of degenerate biarc).
    //!
    void push_back( PolyLine const & c );

    //!
    //! Construct a biarc passing from the points
    //! \f$ (x_0,y_0) \f$ to the point \f$ (x_1,y_1) \f$
    //! with initial angle \f$ \theta_0 \f$ and final angle \f$ \theta_1 \f$
    //! and append the biarc to the tail of biarc list.
    //! The initial point and angle is taken from the tail of the biarc list.
    //!
    //! \param[in] x1      \f$ x_1      \f$
    //! \param[in] y1      \f$ y_1      \f$
    //! \param[in] theta1  \f$ \theta_1 \f$
    //!
    void push_back_G1( real_type x1, real_type y1, real_type theta1 );

    //!
    //! Construct a biarc passing from the points
    //! \f$ (x_0,y_0) \f$ to the point \f$ (x_1,y_1) \f$
    //! with initial angle \f$ \theta_0 \f$ and final angle \f$ \theta_1 \f$
    //! and append the biarc to the tail of biarc list.
    //!
    //! \param[in] x0      \f$ x_0      \f$
    //! \param[in] y0      \f$ y_0      \f$
    //! \param[in] theta0  \f$ \theta_0 \f$
    //! \param[in] x1      \f$ x_1      \f$
    //! \param[in] y1      \f$ y_1      \f$
    //! \param[in] theta1  \f$ \theta_1 \f$
    //!
    void
    push_back_G1(
      real_type x0, real_type y0, real_type theta0,
      real_type x1, real_type y1, real_type theta1
    );

    //!
    //! Construct a biarc list passing to the points \f$ (x_i,y_i) \f$.
    //!
    //! \param[in] n number of points
    //! \param[in] x \f$x\f$-coordinates
    //! \param[in] y \f$y\f$-coordinates
    //!
    bool
    build_G1(
      integer         n,
      real_type const x[],
      real_type const y[]
    );

    //!
    //! Construct a biarc list passing to the points \f$ (x_i,y_i) \f$
    //! with angles  \f$ \theta_i \f$.
    //!
    //! \param[in] n     number of points
    //! \param[in] x     \f$x\f$-coordinates
    //! \param[in] y     \f$y\f$-coordinates
    //! \param[in] theta angles at nodes
    //!
    bool
    build_G1(
      integer         n,
      real_type const x[],
      real_type const y[],
      real_type const theta[]
    );

    //!
    //! Get the `idx`-th biarc.
    //!
    Biarc const & get( integer idx ) const;

    //!
    //! Get the biarc that contain the curvilinear coordinate \f$s\f$.
    //!
    Biarc const & get_at_s( real_type s ) const;

    //!
    //! Return the number of biarc in the biarc list.
    //!
    integer num_segments() const { return integer(m_biarc_list.size()); }

    //!
    //! Get the of the biarc that contain
    //! the curvilinear coordinate \f$s\f$.
    //!
    integer find_at_s( real_type & s ) const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type length() const override;
    real_type length_ISO( real_type offs ) const override;

    //!
    //! The length of the `nseg`-th biarc.
    //!
    real_type
    segment_length( integer nseg ) const;

    //!
    //! The length of the `nseg`-th biarc with offset `offs`.
    //!
    real_type
    segment_length_ISO( integer nseg, real_type offs ) const;

    /*\
     |  _    _   _____    _                _
     | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
     | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
     | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
     |                               |___/
    \*/

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
    ) const override;

    void
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override;

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    //!
    //! Build the internal AABB tree of the biarc list with offset (ISO)
    //!
    //! \param[out] offs      curve offset
    //! \param[out] max_angle maximum angle variation of the arc covered by a triangle
    //! \param[out] max_size  maximum admissible size of the covering tirnagles
    //!
    void
    build_AABBtree_ISO(
      real_type offs,
      real_type max_angle = Utils::m_pi/6, // 30 degree
      real_type max_size  = 1e100
    ) const;

    //!
    //! Build the internal AABB tree of the biarc list with offset (SAE)
    //!
    //! \param[out] offs      curve offset
    //! \param[out] max_angle maximum angle variation of the arc covered by a triangle
    //! \param[out] max_size  maximum admissible size of the covering tirnagles
    //!
    void
    build_AABBtree_SAE(
      real_type offs,
      real_type max_angle = Utils::m_pi/6, // 30 degree
      real_type max_size  = 1e100
    ) const {
      build_AABBtree_ISO( -offs, max_angle, max_size );
    }
    #endif

    /*\
     |   _     _
     |  | |__ | |__   _____  __
     |  | '_ \| '_ \ / _ \ \/ /
     |  | |_) | |_) | (_) >  <
     |  |_.__/|_.__/ \___/_/\_\
    \*/

    void
    bbox(
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
    ) const override {
      bbox_ISO( 0, xmin, ymin, xmax, ymax );
    }

    void
    bbox_ISO(
      real_type   offs,
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
    ) const override;

    /*\
     |   ____             _          _______           _
     |  | __ )  ___  __ _(_)_ __    / / ____|_ __   __| |
     |  |  _ \ / _ \/ _` | | '_ \  / /|  _| | '_ \ / _` |
     |  | |_) |  __/ (_| | | | | |/ / | |___| | | | (_| |
     |  |____/ \___|\__, |_|_| |_/_/  |_____|_| |_|\__,_|
     |              |___/
    \*/

    real_type
    theta_begin() const override
    { return m_biarc_list.front().theta_begin(); }

    real_type
    theta_end() const override
    { return m_biarc_list.back().theta_end(); }

    real_type
    x_begin() const override
    { return m_biarc_list.front().x_begin(); }

    real_type
    y_begin() const override
    { return m_biarc_list.front().y_begin(); }

    real_type
    x_end() const override
    { return m_biarc_list.back().x_end(); }

    real_type
    y_end() const override
    { return m_biarc_list.back().y_end(); }

    real_type
    x_begin_ISO( real_type offs ) const override
    { return m_biarc_list.front().x_begin_ISO( offs ); }

    real_type
    y_begin_ISO( real_type offs ) const override
    { return m_biarc_list.front().y_begin_ISO( offs ); }

    real_type
    x_end_ISO( real_type offs ) const override
    { return m_biarc_list.back().x_end_ISO( offs ); }

    real_type
    y_end_ISO( real_type offs ) const override
    { return m_biarc_list.back().y_end_ISO( offs ); }

    real_type
    tx_begin() const override
    { return m_biarc_list.front().tx_begin(); }

    real_type
    ty_begin() const override
    { return m_biarc_list.front().ty_begin(); }

    real_type
    tx_end() const override
    { return m_biarc_list.back().tx_end(); }

    real_type
    ty_end() const override
    { return m_biarc_list.back().ty_end(); }

    real_type
    nx_begin_ISO() const override
    { return m_biarc_list.front().nx_begin_ISO(); }

    real_type
    ny_begin_ISO() const override
    { return m_biarc_list.front().ny_begin_ISO(); }

    real_type
    nx_end_ISO() const override
    { return m_biarc_list.back().nx_end_ISO(); }

    real_type
    ny_end_ISO() const override
    { return m_biarc_list.back().ny_end_ISO(); }

    /*\
     |  _   _          _
     | | |_| |__   ___| |_ __ _
     | | __| '_ \ / _ \ __/ _` |
     | | |_| | | |  __/ || (_| |
     |  \__|_| |_|\___|\__\__,_|
    \*/

    real_type theta    ( real_type ) const override;
    real_type theta_D  ( real_type ) const override;
    real_type theta_DD ( real_type ) const override;
    real_type theta_DDD( real_type ) const override;

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type ) const override;
    real_type ty    ( real_type ) const override;
    real_type tx_D  ( real_type ) const override;
    real_type ty_D  ( real_type ) const override;
    real_type tx_DD ( real_type ) const override;
    real_type ty_DD ( real_type ) const override;
    real_type tx_DDD( real_type ) const override;
    real_type ty_DDD( real_type ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg(
      real_type   s,
      real_type & tg_x,
      real_type & tg_y
    ) const override;

    void
    tg_D(
      real_type   s,
      real_type & tg_x_D,
      real_type & tg_y_D
    ) const override;

    void
    tg_DD(
      real_type   s,
      real_type & tg_x_DD,
      real_type & tg_y_DD
    ) const override;

    void
    tg_DDD(
      real_type   s,
      real_type & tg_x_DDD,
      real_type & tg_y_DDD
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate(
      real_type   s,
      real_type & th,
      real_type & k,
      real_type & x,
      real_type & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate_ISO(
      real_type   s,
      real_type   offs,
      real_type & th,
      real_type & k,
      real_type & x,
      real_type & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X    ( real_type ) const override;
    real_type Y    ( real_type ) const override;
    real_type X_D  ( real_type ) const override;
    real_type Y_D  ( real_type ) const override;
    real_type X_DD ( real_type ) const override;
    real_type Y_DD ( real_type ) const override;
    real_type X_DDD( real_type ) const override;
    real_type Y_DDD( real_type ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const override;

    void
    eval_D(
      real_type   s,
      real_type & x_D,
      real_type & y_D
    ) const override;

    void
    eval_DD(
      real_type   s,
      real_type & x_DD,
      real_type & y_DD
    ) const override;

    void
    eval_DDD(
      real_type   s,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override;

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    real_type X_ISO    ( real_type s, real_type offs ) const override;
    real_type Y_ISO    ( real_type s, real_type offs ) const override;
    real_type X_ISO_D  ( real_type s, real_type offs ) const override;
    real_type Y_ISO_D  ( real_type s, real_type offs ) const override;
    real_type X_ISO_DD ( real_type s, real_type offs ) const override;
    real_type Y_ISO_DD ( real_type s, real_type offs ) const override;
    real_type X_ISO_DDD( real_type s, real_type offs ) const override;
    real_type Y_ISO_DDD( real_type s, real_type offs ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const override;

    void
    eval_ISO_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const override;

    void
    eval_ISO_DD(
      real_type   s,
      real_type   offs,
      real_type & x_DD,
      real_type & y_DD
    ) const override;

    void
    eval_ISO_DDD(
      real_type   s,
      real_type   offs,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override;

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void translate( real_type tx, real_type ty ) override;
    void rotate( real_type angle, real_type cx, real_type cy ) override;
    void scale( real_type sc ) override;
    void reverse() override;
    void change_origin( real_type newx0, real_type newy0 ) override;
    void trim( real_type s_begin, real_type s_end ) override;

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    integer
    closest_point_ISO(
      real_type   qx,
      real_type   qy,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst
    ) const override;

    integer // true if projection is unique and orthogonal
    closest_point_ISO(
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & t,
      real_type & dst
    ) const override;

    string info() const;

    void
    info( ostream_type & stream ) const override
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, BiarcList const & CL );

    //!
    //! Return the biarc as a list of nodes angles and curvatures.
    //!
    //! \param[out] s     nodes
    //! \param[out] theta angles
    //! \param[out] kappa curvature
    //!
    void
    get_STK(
      real_type s[],
      real_type theta[],
      real_type kappa[]
    ) const;

    //!
    //! Return the biarc XY nodes
    //!
    //! \param[out] x \f$ x \f$-nodes
    //! \param[out] y \f$ y \f$-nodes
    //!
    void
    get_XY( real_type x[], real_type y[] ) const;

    //!
    //! Find parametric coordinate.
    //!
    //! \param  x    \f$x\f$-coordinate point
    //! \param  y    \f$y\f$-coordinate point
    //! \param  s    value \f$ s \f$
    //! \param  t    value \f$ t \f$
    //! \return idx  the segment with point at minimal distance, otherwise
    //!              `-(idx+1)` if \f$(x,y)\f$ cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      real_type   x,
      real_type   y,
      real_type & s,
      real_type & t
    ) const;

    //!
    //! Find parametric coordinate.
    //!
    //! \param  ibegin initial segment to compute the distance
    //! \param  iend   final segment to compute the distance
    //! \param  x      \f$x\f$-coordinate point
    //! \param  y      \f$y\f$-coordinate point
    //! \param  s      value \f$ s \f$
    //! \param  t      value \f$ t \f$
    //! \return idx    the segment with point at minimal distance, otherwise
    //!                `-(idx+1)` if \f$(x,y)\f$ cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      integer     ibegin,
      integer     iend,
      real_type   x,
      real_type   y,
      real_type & s,
      real_type & t
    ) const;

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    //!
    //! Detect a collision with another biarc list.
    //!
    bool
    collision( BiarcList const & BL ) const {
      return collision_ISO( 0, BL, 0 );
    }

    //!
    //! Detect a collision with another biarc list with offset.
    //!
    //! \param[in] offs   offset of first biarc
    //! \param[in] BL     second biarc
    //! \param[in] offs_C offset of second biarc
    //!
    bool
    collision_ISO(
      real_type         offs,
      BiarcList const & BL,
      real_type         offs_C
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
    //! intersect a biarc list with another biarc list
    //!
    //! \param[in]  BL    second biarc
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      BiarcList const & BL,
      IntersectList   & ilist
    ) const {
      this->intersect_ISO( 0, BL, 0, ilist );
    }

    //!
    //! Intersect a biarc list with another biarc list with offset (ISO).
    //!
    //! \param[in]  offs     offset of first biarc
    //! \param[in]  BL       second biarc
    //! \param[in]  offs_obj offset of second biarc
    //! \param[out] ilist    list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type         offs,
      BiarcList const & BL,
      real_type         offs_obj,
      IntersectList   & ilist
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

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "BiarcList_compatibility.hxx"
#endif

  };

}

///
/// eof: BiarcList.hxx
///
