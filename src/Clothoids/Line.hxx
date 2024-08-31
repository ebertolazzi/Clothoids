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
/// file: Line.hxx
///

namespace G2lib {

  /*\
   |   _     _
   |  | |   (_)_ __   ___
   |  | |   | | '_ \ / _ \
   |  | |___| | | | |  __/
   |  |_____|_|_| |_|\___|
  \*/

  //!
  //! Class to manage a straight segment
  //!
  class LineSegment : public BaseCurve {

    friend class CircleArc;
    friend class PolyLine;

    real_type m_x0{0};     //!< initial \f$x\f$-coordinate of the line
    real_type m_y0{0};     //!< initial \f$y\f$-coordinate of the line
    real_type m_theta0{0}; //!< angle of the line

    real_type m_c0{1};     //!< `cos(theta0)`
    real_type m_s0{0};     //!< `sin(theta0)`
    real_type m_L{0};      //!< length of the segment

  public:

    #include "BaseCurve_using.hxx"

    LineSegment() = delete;
    LineSegment( string const & name ) : BaseCurve( name ) {};

    void setup( GenericContainer const & gc ) override;

    LineSegment( LineSegment const & s ) : BaseCurve( s.name() )
    { this->copy(s); }

    explicit
    LineSegment( BaseCurve const * pC );

    //!
    //! Construct a circle curve with the standard parameters
    //!
    explicit
    LineSegment(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type L,
      string const & name
    )
    : BaseCurve( name )
    , m_x0(x0)
    , m_y0(y0)
    , m_theta0(theta0)
    , m_c0(cos(theta0))
    , m_s0(sin(theta0))
    , m_L(L)
    {}

    void
    copy( LineSegment const & c ) {
      m_x0     = c.m_x0;
      m_y0     = c.m_y0;
      m_theta0 = c.m_theta0;
      m_c0     = c.m_c0;
      m_s0     = c.m_s0;
      m_L      = c.m_L;
    }

    LineSegment const & operator = ( LineSegment const & s )
    { this->copy(s); return *this; }

    CurveType type() const override { return CurveType::LINE; }

    real_type
    length() const override
    { return m_L; }

    real_type
    length_ISO( real_type ) const override
    { return m_L; }

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
    ) const override;

    void
    bbox_ISO(
      real_type   offs,
      real_type & xmin,
      real_type & ymin,
      real_type & xmax,
      real_type & ymax
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
      real_type            max_size  = 1e100, // unused
      integer              icurve    = 0
    ) const override;

    void
    bb_triangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100, // unused
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

    /*\
     |   ____             _          _______           _
     |  | __ )  ___  __ _(_)_ __    / / ____|_ __   __| |
     |  |  _ \ / _ \/ _` | | '_ \  / /|  _| | '_ \ / _` |
     |  | |_) |  __/ (_| | | | | |/ / | |___| | | | (_| |
     |  |____/ \___|\__, |_|_| |_/_/  |_____|_| |_|\__,_|
     |              |___/
    \*/

    real_type tx_begin() const override { return m_c0; }
    real_type ty_begin() const override { return m_s0; }
    real_type tx_end()   const override { return m_c0; }
    real_type ty_end()   const override { return m_s0; }

    real_type nx_begin_ISO() const override { return -m_s0; }
    real_type ny_begin_ISO() const override { return m_c0; }
    real_type nx_end_ISO()   const override { return -m_s0; }
    real_type ny_end_ISO()   const override { return m_c0; }

    real_type x_begin() const override { return m_x0; }
    real_type y_begin() const override { return m_y0; }
    real_type x_end()   const override { return m_x0+m_L*m_c0; }
    real_type y_end()   const override { return m_y0+m_L*m_s0; }

    real_type
    x_begin_ISO( real_type offs ) const override
    { return m_x0+offs*nx_begin_ISO(); }

    real_type
    y_begin_ISO( real_type offs ) const override
    { return m_y0+offs*ny_begin_ISO(); }

    real_type
    x_end_ISO( real_type offs ) const override
    { return x_end()+offs*nx_begin_ISO(); }

    real_type
    y_end_ISO( real_type offs ) const override
    { return y_end()+offs*ny_begin_ISO(); }

    /*\
     |  _   _          _
     | | |_| |__   ___| |_ __ _
     | | __| '_ \ / _ \ __/ _` |
     | | |_| | | |  __/ || (_| |
     |  \__|_| |_|\___|\__\__,_|
    \*/

    real_type theta    ( real_type ) const override { return m_theta0; }
    real_type theta_D  ( real_type ) const override { return 0; }
    real_type theta_DD ( real_type ) const override { return 0; }
    real_type theta_DDD( real_type ) const override { return 0; }


    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type ) const override { return m_c0; }
    real_type ty    ( real_type ) const override { return m_s0; }
    real_type tx_D  ( real_type ) const override { return 0; }
    real_type ty_D  ( real_type ) const override { return 0; }
    real_type tx_DD ( real_type ) const override { return 0; }
    real_type ty_DD ( real_type ) const override { return 0; }
    real_type tx_DDD( real_type ) const override { return 0; }
    real_type ty_DDD( real_type ) const override { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg( real_type, real_type & tx, real_type & ty ) const override
    { tx = m_c0; ty = m_s0; }

    void
    tg_D( real_type, real_type & tx_D, real_type & ty_D ) const override
    { tx_D = ty_D = 0; }

    void
    tg_DD( real_type, real_type & tx_DD, real_type & ty_DD ) const override
    { tx_DD = ty_DD = 0; }

    void
    tg_DDD( real_type, real_type & tx_DDD, real_type & ty_DDD ) const override
    { tx_DDD = ty_DDD = 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X( real_type s ) const override { return m_x0+s*m_c0; }
    real_type Y( real_type s ) const override { return m_y0+s*m_s0; }

    real_type X_D( real_type ) const override { return m_c0; }
    real_type Y_D( real_type ) const override { return m_s0; }

    real_type X_DD( real_type ) const override { return 0; }
    real_type Y_DD( real_type ) const override { return 0; }

    real_type X_DDD( real_type ) const override { return 0; }
    real_type Y_DDD( real_type ) const override { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const override {
      x = m_x0+s*m_c0;
      y = m_y0+s*m_s0;
    }

    void
    eval_D(
      real_type,
      real_type & x_D,
      real_type & y_D
    ) const override {
      x_D = m_c0;
      y_D = m_s0;
    }

    void
    eval_DD(
      real_type,
      real_type & x_DD,
      real_type & y_DD
    ) const override {
      x_DD = 0;
      y_DD = 0;
    }

    void
    eval_DDD(
      real_type,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override {
      x_DDD = 0;
      y_DDD = 0;
    }

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    real_type
    X_ISO( real_type s, real_type offs ) const override
    { return m_x0 + s*m_c0 + offs*nx_begin_ISO(); }

    real_type
    Y_ISO( real_type s, real_type offs ) const override
    { return m_y0 + s*m_s0 + offs*ny_begin_ISO(); }

    real_type
    X_ISO_D( real_type, real_type ) const override
    { return m_c0; }

    real_type
    Y_ISO_D( real_type, real_type ) const override
    { return m_s0; }

    real_type
    X_ISO_DD( real_type, real_type ) const override
    { return 0; }

    real_type
    Y_ISO_DD( real_type, real_type ) const override
    { return 0; }

    real_type
    X_ISO_DDD( real_type, real_type ) const override
    { return 0; }

    real_type
    Y_ISO_DDD( real_type, real_type ) const override
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const override {
      x = m_x0 + s*m_c0 + offs*nx_begin_ISO();
      y = m_y0 + s*m_s0 + offs*ny_begin_ISO();
    }

    void
    eval_ISO_D(
      real_type,
      real_type,
      real_type & x_D,
      real_type & y_D
    ) const override {
      x_D = m_c0;
      y_D = m_s0;
    }

    void
    eval_ISO_DD(
      real_type,
      real_type,
      real_type & x_DD,
      real_type & y_DD
    ) const override {
      x_DD = y_DD = 0;
    }

    void
    eval_ISO_DDD(
      real_type,
      real_type,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override {
      x_DDD = y_DDD = 0;
    }

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void
    translate( real_type tx, real_type ty ) override
    { m_x0 += tx; m_y0 += ty; }

    void
    rotate( real_type angle, real_type cx, real_type cy ) override;

    void
    reverse() override;

    void
    change_origin( real_type newx0, real_type newy0 ) override
    { m_x0 = newx0; m_y0 = newy0; }

    void
    scale( real_type sc ) override
    { m_L *= sc; }

    void
    trim( real_type s_begin, real_type s_end ) override {
      m_x0 += m_c0 * s_begin;
      m_y0 += m_s0 * s_begin;
      m_L   = s_end - s_begin;
    }

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    //!
    //! Compute the point at minimum distance from a point `[x,y]` and the line segment
    //!
    //! \param[in]  qx  \f$x\f$-coordinate
    //! \param[in]  qy  \f$y\f$-coordinate
    //! \param[out] x   \f$x\f$-coordinate of the closest point
    //! \param[out] y   \f$y\f$-coordinate of the closest point
    //! \param[out] s   param of the closest point
    //! \param[out] t   signed distance if projection is orthogonal to segment
    //! \param[out] dst signed distance from the segment
    //! \return 1 = point is projected orthogonal
    //!         0 = more than one projection (first returned)
    //!        -1 = minimum point is not othogonal projection to curve
    //!

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

    integer
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
    operator << ( ostream_type & stream, LineSegment const & c );

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type L
    ) {
      m_x0     = x0;
      m_y0     = y0;
      m_theta0 = theta0;
      m_c0     = cos(theta0);
      m_s0     = sin(theta0);
      m_L      = L;
    }

    //!
    //! Construct a clothoid with the standard parameters
    //!
    void
    build_2P(
      real_type _x0,
      real_type _y0,
      real_type _x1,
      real_type _y1
    );

    //!
    //! Construct a clothoid with the standard parameters
    //!
    void
    build_2P( real_type const p0[2], real_type const p1[2] )
    { build_2P( p0[0], p0[1], p1[0], p1[1] ); }

    void
    p1p2( real_type p1[2], real_type p2[2] ) const {
      p1[0] = m_x0;
      p1[1] = m_y0;
      p2[0] = m_x0+m_L*m_c0;
      p2[1] = m_y0+m_L*m_s0;
    }

    void build( LineSegment const & LS );
    void build( CircleArc const & );
    void build( Biarc const & );
    void build( ClothoidCurve const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    /*
    //             _ _ _     _
    //    ___ ___ | | (_)___(_) ___  _ __
    //   / __/ _ \| | | / __| |/ _ \| '_ \
    //  | (_| (_) | | | \__ \ | (_) | | | |
    //   \___\___/|_|_|_|___/_|\___/|_| |_|
    */

    bool
    collision( LineSegment const & S ) const;

    bool
    collision_ISO(
      real_type           offs,
      LineSegment const & S,
      real_type           S_offs
    ) const;

    bool
    collision( BaseCurve const * pC ) const override;

    bool
    collision_ISO(
      real_type         offs,
      BaseCurve const * pC,
      real_type         offs_C
    ) const override;

    /*
    //   _       _                          _
    //  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
    //  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
    //  | | | | | ||  __/ |  \__ \  __/ (__| |_
    //  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    */

    bool
    intersect(
      LineSegment const & S,
      real_type         & s1,
      real_type         & s2
    ) const;

    bool
    intersect_ISO(
      real_type           offs,
      LineSegment const & S,
      real_type           S_offs,
      real_type         & s1,
      real_type         & s2
    ) const;

    void
    intersect(
      LineSegment const & LS,
      IntersectList     & ilist
    ) const;

    void
    intersect_ISO(
      real_type           offs,
      LineSegment const & LS,
      real_type           offs_LS,
      IntersectList     & ilist
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

    /*\
     |   _   _ _   _ ____  ____ ____
     |  | \ | | | | |  _ \| __ ) ___|
     |  |  \| | | | | |_) |  _ \___ \
     |  | |\  | |_| |  _ <| |_) |__) |
     |  |_| \_|\___/|_| \_\____/____/
    \*/

    void
    paramNURBS( integer & n_knots, integer & n_pnts ) const;

    void
    toNURBS( real_type knots[], real_type Poly[][3] ) const;

    virtual
    void
    toBS( real_type knots[], real_type Poly[][2] ) const;

    friend class ClothoidCurve;

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "Line_compatibility.hxx"
#endif

  };

}

///
/// eof: Line.hxx
///
