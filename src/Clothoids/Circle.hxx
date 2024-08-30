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
/// file: Circle.hxx
///

namespace G2lib {

  /*\
   |    ____ _          _         _
   |   / ___(_)_ __ ___| | ___   / \   _ __ ___
   |  | |   | | '__/ __| |/ _ \ / _ \ | '__/ __|
   |  | |___| | | | (__| |  __// ___ \| | | (__
   |   \____|_|_|  \___|_|\___/_/   \_\_|  \___|
  \*/

  //!
  //! Class to manage a circle arc
  //!
  class CircleArc : public BaseCurve {

    friend class Biarc;

    real_type m_x0{0};     //!< initial \f$x\f$-coordinate of the clothoid
    real_type m_y0{0};     //!< initial \f$y\f$-coordinate of the clothoid
    real_type m_theta0{0}; //!< initial angle of the clothoid
    real_type m_c0{1};     //!< initial \f$\cos(\f$ angle \f$)\f$ of the clothoid
    real_type m_s0{0};     //!< initial \f$\sin(\f$ angle \f$)\f$ of the clothoid
    real_type m_k{0};      //!< curvature

    real_type m_L{0};      //!< length of the circle segment

  public:

    #include "BaseCurve_using.hxx"

    //!
    //! Build an empty circle
    //!
    CircleArc() = delete;
    CircleArc( string const & name ) : BaseCurve( name ) {};

    void setup( GenericContainer const & gc ) override;

    //!
    //! Build a copy of an existing circle arc.
    //!
    CircleArc( CircleArc const & s ) : BaseCurve( s.name() )
    { this->copy(s); }

    //!
    //! Construct a circle arc with the standard parameters.
    //!
    //! \param[in] x0     starting position \f$x\f$-coordinate
    //! \param[in] y0     starting position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] k      curvature
    //! \param[in] L      length
    //! \param[in] name   name of the circle arc
    //!
    explicit
    CircleArc(
      real_type      x0,
      real_type      y0,
      real_type      theta0,
      real_type      k,
      real_type      L,
      string const & name
    )
    : BaseCurve( name )
    , m_x0(x0)
    , m_y0(y0)
    , m_theta0(theta0)
    , m_c0(cos(theta0))
    , m_s0(sin(theta0))
    , m_k(k)
    , m_L(L)
    {}

    //!
    //! Construct a circle arc from a line
    //! segment (degenerate circle).
    //!
    explicit
    CircleArc( LineSegment const & LS )
    : BaseCurve( LS.name() )
    , m_x0(LS.x_begin())
    , m_y0(LS.y_begin())
    , m_theta0(LS.m_theta0)
    , m_c0(LS.m_c0)
    , m_s0(LS.m_s0)
    , m_k(0)
    , m_L(LS.length())
    {}

    //!
    //! Make a copy of an existing circle arc.
    //!
    void
    copy( CircleArc const & c ) {
      m_x0     = c.m_x0;
      m_y0     = c.m_y0;
      m_theta0 = c.m_theta0;
      m_c0     = c.m_c0;
      m_s0     = c.m_s0;
      m_k      = c.m_k;
      m_L      = c.m_L;
    }

    //!
    //! Build a circle arc from a generic curve (if possibile).
    //!
    explicit
    CircleArc( BaseCurve const * pC );

    CurveType type() const override { return CurveType::CIRCLE; }

    //!
    //! Make a copy of an existing circle arc.
    //!
    CircleArc const &
    operator = ( CircleArc const & s )
    { this->copy(s); return *this; }

    //!
    //! Construct a circle arc with the standard parameters.
    //!
    //! \param[in] x0     starting position \f$x\f$-coordinate
    //! \param[in] y0     starting position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] k      curvature
    //! \param[in] L      length
    //!
    void
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type k,
      real_type L
    ) {
      m_x0     = x0;
      m_y0     = y0;
      m_theta0 = theta0;
      m_k      = k;
      m_L      = L;
    }

    //!
    //! Build a circle by solving the hermite \f$ G^1 \f$ problem.
    //!
    //! \param[in] x0     starting position \f$x\f$-coordinate
    //! \param[in] y0     starting position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position \f$x\f$-coordinate
    //! \param[in] y1     final position \f$y\f$-coordinate
    //! \return true if success
    //!
    bool
    build_G1(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1
    );

    //!
    //! Build a circle passing by 3 points.
    //!
    //! \param[in] x0 starting point \f$x\f$-coordinate
    //! \param[in] y0 starting point \f$y\f$-coordinate
    //! \param[in] x1 intermediate point \f$x\f$-coordinate
    //! \param[in] y1 intermediate point \f$y\f$-coordinate
    //! \param[in] x2 final point \f$x\f$-coordinate
    //! \param[in] y2 final point \f$y\f$-coordinate
    //! \return true if success
    //!
    bool
    build_3P(
      real_type x0,
      real_type y0,
      real_type x1,
      real_type y1,
      real_type x2,
      real_type y2
    );

    //!
    //! Construct a circle arc from a line
    //! segment (degenerate circle).
    //!
    void build( LineSegment const & );
    void build( CircleArc const & );
    void build( Biarc const & );
    void build( ClothoidCurve const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    //!
    //! Get the bounding box triangle (if angle variation less that pi/3).
    //!
    //! \param[out] x0 first triangle point \f$x\f$-coordinate
    //! \param[out] y0 first triangle point \f$y\f$-coordinate
    //! \param[out] x1 second triangle point \f$x\f$-coordinate
    //! \param[out] y1 second triangle point \f$y\f$-coordinate
    //! \param[out] x2 third triangle point \f$x\f$-coordinate
    //! \param[out] y2 third triangle point \f$y\f$-coordinate
    //! \return true if success
    //!
    bool
    bbTriangle(
      real_type & x0, real_type & y0,
      real_type & x1, real_type & y1,
      real_type & x2, real_type & y2
    ) const;

    //!
    //! Get the bounding box triangle of the circle arc with offset
    //! (if angle variation less that \f$ \pi/3 \f$).
    //!
    //! \param[in]  offs offset
    //! \param[out] x0   first triangle point \f$x\f$-coordinate
    //! \param[out] y0   first triangle point \f$y\f$-coordinate
    //! \param[out] x1   second triangle point \f$x\f$-coordinate
    //! \param[out] y1   second triangle point \f$y\f$-coordinate
    //! \param[out] x2   third triangle point \f$x\f$-coordinate
    //! \param[out] y2   third triangle point \f$y\f$-coordinate
    //! \return true if success
    //!
    bool
    bbTriangle_ISO(
      real_type   offs,
      real_type & x0, real_type & y0,
      real_type & x1, real_type & y1,
      real_type & x2, real_type & y2
    ) const;

    //!
    //! Get the bounding box triangle of the circle arc with offset
    //! if angle variation less that \f$ \pi/3 \f$).
    //!
    //! \param[in]  offs offset
    //! \param[out] x0   first triangle point \f$x\f$-coordinate
    //! \param[out] y0   first triangle point \f$y\f$-coordinate
    //! \param[out] x1   second triangle point \f$x\f$-coordinate
    //! \param[out] y1   second triangle point \f$y\f$-coordinate
    //! \param[out] x2   third triangle point \f$x\f$-coordinate
    //! \param[out] y2   third triangle point \f$y\f$-coordinate
    //! \return true if success
    //!
    bool
    bbTriangle_SAE(
      real_type   offs,
      real_type & x0, real_type & y0,
      real_type & x1, real_type & y1,
      real_type & x2, real_type & y2
    ) const {
      return this->bbTriangle_ISO( -offs, x0, y0, x1, y1, x2, y2 );
    }

    //!
    //! Get the bounding box triangle of the circle
    //! (if angle variation less that \f$ \pi/3 \f$).
    //!
    //! \param[out] p0 first triangle point
    //! \param[out] p1 second triangle point
    //! \param[out] p2 third triangle point
    //! \return true if success
    //!
    bool
    bbTriangle(
      real_type p0[],
      real_type p1[],
      real_type p2[]
    ) const {
      return bbTriangle( p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );
    }

    //!
    //! Get the bounding box triangle of the circle arc with offset
    //! (if angle variation less that \f$ \pi/3 \f$).
    //!
    //! \param[in]  offs offset
    //! \param[out] p0   first triangle point
    //! \param[out] p1   second triangle point
    //! \param[out] p2   third triangle point
    //! \return true if success
    //!
    bool
    bbTriangle_ISO(
      real_type   offs,
      real_type p0[],
      real_type p1[],
      real_type p2[]
    ) const {
      return bbTriangle_ISO( offs, p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );
    }

    //!
    //! Get the bounding box triangle of the circle arc with offset
    //! (if angle variation less that \f$ \pi/3 \f$).
    //!
    //! \param[in]  offs offset
    //! \param[out] p0   first triangle point
    //! \param[out] p1   second triangle point
    //! \param[out] p2   third triangle point
    //! \return true if success
    //!
    bool
    bbTriangle_SAE(
      real_type offs,
      real_type p0[],
      real_type p1[],
      real_type p2[]
    ) const {
      return bbTriangle_SAE( offs, p0[0], p0[1], p1[0], p1[1], p2[0], p2[1] );
    }

    //!
    //! Get the bounding box triangle list of the circle arc.
    //!
    //! \param[out] t      the bounding triangle
    //! \param[in]  ss0    `s0` stored in the triangle class
    //! \param[in]  ss1    `s1` stored in the triangle class
    //! \param[in]  icurve `id` stored in the triangle class
    //! \return true if success
    //!
    bool
    bbTriangle(
      Triangle2D & t,
      real_type    ss0    = 0,
      real_type    ss1    = 0,
      integer      icurve = 0
    ) const {
      real_type p0[2], p1[2], p2[2];
      bool ok = bbTriangle( p0, p1, p2 );
      if ( ok ) t.build( p0, p1, p2, ss0, ss1, icurve );
      return ok;
    }

    //!
    //! Get the bounding box triangle list
    //! of the circle arc with offset.
    //!
    //! \param[in]  offs   offset
    //! \param[out] t      the bounding triangle
    //! \param[in]  ss0    `s0` stored in the triangle class
    //! \param[in]  ss1    `s1` stored in the triangle class
    //! \param[in]  icurve `id` stored in the triangle class
    //! \return true if success
    //!
    bool
    bbTriangle_ISO(
      real_type    offs,
      Triangle2D & t,
      real_type    ss0    = 0,
      real_type    ss1    = 0,
      integer      icurve = 0
    ) const {
      real_type p0[2], p1[2], p2[2];
      bool ok = bbTriangle_ISO( offs, p0, p1, p2 );
      if ( ok ) t.build( p0, p1, p2, ss0, ss1, icurve );
      return ok;
    }

    //!
    //! Get the bounding box triangle list of
    //! the circle arc with offset.
    //!
    //! \param[in]  offs   offset
    //! \param[out] t      the bounding triangle
    //! \param[in]  ss0    `s0` stored in the triangle class
    //! \param[in]  ss1    `s1` stored in the triangle class
    //! \param[in]  icurve `id` stored in the triangle class
    //! \return true if success
    //!
    bool
    bbTriangle_SAE(
      real_type    offs,
      Triangle2D & t,
      real_type    ss0    = 0,
      real_type    ss1    = 0,
      integer      icurve = 0
    ) const {
      return this->bbTriangle_ISO( -offs, t, ss0, ss1, icurve );
    }

    //!
    //! Get the bounding box triangle list of the circle arc with offset.
    //!
    //! \param[out] tvec      the bounding triangle list
    //! \param[in]  max_angle maximum angle variation admitted for all splitted segment
    //! \param[in]  max_size  maximum size admitted for all splitted segment
    //! \param[in]  icurve    `id` stored in the triangles
    //!
    void
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type max_angle = Utils::m_pi/18,
      real_type max_size  = 1e100,
      integer   icurve    = 0
    ) const override; // 10 degree

    //!
    //! Get the bounding box triangle list of the circle arc with offset.
    //!
    //! \param[in]  offs      offset
    //! \param[out] tvec      the bounding triangle list
    //! \param[in]  max_angle maximum angle variation admitted for all splitted segment
    //! \param[in]  max_size  maximum size admitted for all splitted segment
    //! \param[in]  icurve    `id` stored in the triangles
    //!
    void
    bb_triangles_ISO(
      real_type offs,
      vector<Triangle2D> & tvec,
      real_type max_angle = Utils::m_pi/18,
      real_type max_size  = 1e100,
      integer   icurve    = 0
    ) const override; // 10 degree

    //!
    //! Get the bounding box triangle list of the circle arc with offset.
    //!
    //! \param[in]  offs      offset
    //! \param[out] tvec      the bounding triangle list
    //! \param[in]  max_angle maximum angle variation admitted for all splitted segment
    //! \param[in]  max_size  maximum size admitted for all splitted segment
    //! \param[in]  icurve    `id` stored in the triangles
    //!
    void
    bb_triangles_SAE(
      real_type offs,
      vector<Triangle2D> & tvec,
      real_type max_angle = Utils::m_pi/18,
      real_type max_size  = 1e100,
      integer   icurve    = 0
    ) const override {
      this->bb_triangles_ISO( -offs, tvec, max_angle, max_size, icurve );
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type
    length() const override
    { return m_L; }

    real_type
    length_ISO( real_type offs ) const override
    { return m_L*(1+m_k*offs); }

    real_type theta_begin()  const override { return m_theta0; }
    real_type kappa_begin()  const override { return m_k; }
    real_type kappa_end()    const override { return m_k; }
    real_type x_begin()      const override { return m_x0; }
    real_type y_begin()      const override { return m_y0; }
    real_type tx_begin()     const override { return m_c0; }
    real_type ty_begin()     const override { return m_s0; }
    real_type nx_begin_ISO() const override { return m_s0; }
    real_type ny_begin_ISO() const override { return -m_c0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type theta    ( real_type s ) const override { return m_theta0 + s*m_k; }
    real_type theta_D  ( real_type   ) const override { return m_k; }
    real_type theta_DD ( real_type   ) const override { return 0; }
    real_type theta_DDD( real_type   ) const override { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate(
      real_type   s,
      real_type & th,
      real_type & kappa,
      real_type & x,
      real_type & y
    ) const override {
      eval( s, x, y );
      th     = m_theta0 + s*m_k;
      kappa  = m_k;
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X( real_type s ) const override;
    real_type Y( real_type s ) const override;

    real_type X_D( real_type ) const override;
    real_type Y_D( real_type ) const override;

    real_type X_DD( real_type ) const override;
    real_type Y_DD( real_type ) const override;

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
      real_type,
      real_type & x_D,
      real_type & y_D
    ) const override;

    void
    eval_DD(
      real_type,
      real_type & x_DD,
      real_type & y_DD
    ) const override;

    void
    eval_DDD(
      real_type,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override;

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx( real_type s ) const override { return cos(theta(s)); }
    real_type ty( real_type s ) const override { return sin(theta(s)); }

    real_type tx_D( real_type s ) const override { return -sin(theta(s))*m_k; }
    real_type ty_D( real_type s ) const override { return cos(theta(s))*m_k; }

    real_type tx_DD( real_type s ) const override { return -cos(theta(s))*m_k*m_k; }
    real_type ty_DD( real_type s ) const override { return -sin(theta(s))*m_k*m_k; }

    real_type tx_DDD( real_type s ) const override { return sin(theta(s))*m_k*m_k*m_k; }
    real_type ty_DDD( real_type s ) const override { return -cos(theta(s))*m_k*m_k*m_k; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void tg( real_type s, real_type & tx, real_type & ty ) const override;
    void tg_D( real_type s, real_type & tx_D, real_type & ty_D ) const override;
    void tg_DD( real_type s, real_type & tx_DD, real_type & ty_DD ) const override;
    void tg_DDD( real_type s, real_type & tx_DDD, real_type & ty_DDD ) const override;

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

    void rotate( real_type angle, real_type cx, real_type cy ) override;
    void reverse() override;

    void
    change_origin( real_type newx0, real_type newy0 ) override
    { m_x0 = newx0; m_y0 = newy0; }

    void scale( real_type s ) override;
    void trim( real_type s_begin, real_type s_end ) override;

    /*\
     |        _                     _   ____       _       _
     |    ___| | ___  ___  ___  ___| |_|  _ \ ___ (_)_ __ | |_
     |   / __| |/ _ \/ __|/ _ \/ __| __| |_) / _ \| | '_ \| __|
     |  | (__| | (_) \__ \  __/\__ \ |_|  __/ (_) | | | | | |_
     |   \___|_|\___/|___/\___||___/\__|_|   \___/|_|_| |_|\__|
     |
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
    operator << ( ostream_type & stream, CircleArc const & bi );

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    //!
    //! Detect a collision with another circle arc.
    //!
    bool
    collision( CircleArc const & ) const;

    //!
    //! Detect a collision with another circle arc with offset.
    //!
    //! \param[in] offs     offset of first circle arc
    //! \param[in] C        second circle arc
    //! \param[in] offs_obj offset of second circle arc
    //!
    bool
    collision_ISO(
      real_type         offs,
      CircleArc const & C,
      real_type         offs_obj
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
    //! Intersect a circle arc with another circle arc.
    //!
    //! \param[in]  obj   second biarc
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      CircleArc const & obj,
      IntersectList   & ilist
    ) const;

    //!
    //! Intersect a circle arc with another circle arc with offset (ISO).
    //!
    //! \param[in]  offs     offset of first circle arc
    //! \param[in]  C        second circle arc
    //! \param[in]  offs_obj offset of second circle arc
    //! \param[out] ilist    list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type         offs,
      CircleArc const & C,
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

    //!
    //! Return \f$ \sin \theta_0 \f$ where
    //! \f$ \theta_0 \f$ is the initial tangent angle.
    //!
    real_type sin_theta0() const { return sin(m_theta0); }

    //!
    //! Return \f$ \cos \theta_0 \f$ where
    //! \f$ \theta_0 \f$ is the initial tangent angle.
    //!
    real_type cos_theta0() const { return cos(m_theta0); }

    //!
    //! Return curvature of the circle arc.
    //!
    real_type curvature() const { return m_k; }

    //!
    //! Return the length of the arc that
    //! can approximated by a line segment.
    //!
    real_type len_tolerance( real_type tol ) const;

    //!
    //! Return the tangent angle variation in the circle arc.
    //!
    real_type delta_theta() const { return std::abs(m_L*m_k); }

    //!
    //! Return the absolute value of the tangent
    //! angle variation in the circle arc.
    //!
    real_type theta_total_variation() const { return delta_theta(); }

    //!
    //! Minimum and maximum tangent angle.
    //!
    //! \param[out] thMin mimimum tangent angle
    //! \param[out] thMax maximum tangent angle
    //! \return `thMax`-`thMin`
    //!
    real_type
    theta_min_max( real_type & thMin, real_type & thMax ) const;

    //!
    //! Change the origin of the circle arc at \f$ s_0 \f$
    //! and the length of the arc to  \f$ L \f$.
    //!
    //! \param[in] s0   \f$ s_0 \f$
    //! \param[in] newL \f$ L   \f$
    //!
    void
    change_curvilinear_origin( real_type s0, real_type newL );

    //!
    //! Get the center of the circle arc \f$ (c_x,c_y) \f$.
    //!
    //! \param[in] cx \f$ c_x \f$
    //! \param[in] cy \f$ c_y \f$
    //!
    void
    center( real_type & cx, real_type & cy ) const;

    //!
    //! Get the ray of the circle arc.
    //!
    real_type ray() const { return 1/std::abs(m_k); }

    /*\
     |   _   _ _   _ ____  ____ ____
     |  | \ | | | | |  _ \| __ ) ___|
     |  |  \| | | | | |_) |  _ \___ \
     |  | |\  | |_| |  _ <| |_) |__) |
     |  |_| \_|\___/|_| \_\____/____/
    \*/

    //!
    //! Get the parameters to build a NURBS for the circle ars.
    //!
    //! \param[out] n_knots number of knots for the NURBS
    //! \param[out] n_pnts  number of point of the polygon of the NURBS
    //!
    void
    paramNURBS( integer & n_knots, integer & n_pnts ) const;

    //!
    //! Get the parameters to build a NURBS for the circle ars.
    //!
    //! \param[out] knots vector of the knots
    //! \param[out] Poly  points of the polygon of the NURBS
    //!
    void
    toNURBS( real_type knots[], real_type Poly[][3] ) const;

    friend class ClothoidCurve;

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "Circle_compatibility.hxx"
#endif

  };

}

///
/// eof: Circle.hxx
///
