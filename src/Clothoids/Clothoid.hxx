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
/// file: Clothoid.hxx
///

namespace G2lib {

  using std::vector;

  /*\
   |    ____ _       _   _           _     _  ____
   |   / ___| | ___ | |_| |__   ___ (_) __| |/ ___|   _ _ ____   _____
   |  | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |  | | | | '__\ \ / / _ \
   |  | |___| | (_) | |_| | | | (_) | | (_| | |__| |_| | |   \ V /  __/
   |   \____|_|\___/ \__|_| |_|\___/|_|\__,_|\____\__,_|_|    \_/ \___|
  \*/
  //!
  //! Class to manage Clothoid Curve.
  //! A clothoid curve is described by the following generalized Fresnel integrals
  //!
  //! \f[
  //!   \begin{cases}
  //!      x(s) = x_0 + \displaystyle\int_0^s \cos(as^2+bs+c)\,\mathrm{d}t \\[1em]
  //!      y(s) = y_0 + \displaystyle\int_0^s \sin(as^2+bs+c)\,\mathrm{d}t
  //!   \end{cases}
  //! \f]
  //!
  //! @html_image{G1problem.png,width=60%}
  //!
  class ClothoidCurve : public BaseCurve {
    friend class ClothoidList;
  private:

    ClothoidData m_CD;  //!< clothoid data
    real_type    m_L;   //!< length of clothoid segment

    void
    optimized_sample_internal_ISO(
      real_type           s_begin,
      real_type           s_end,
      real_type           offs,
      real_type           ds,
      real_type           max_angle,
      vector<real_type> & s
    ) const;

    void
    bb_triangles_internal_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            s0,
      real_type            s1,
      real_type            max_angle,
      real_type            max_size,
      integer              icurve
    ) const;

    void
    closest_point_internal(
      real_type   s_begin,
      real_type   s_end,
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & dst
    ) const;

    void
    closest_point_internal(
      real_type   qx,
      real_type   qy,
      real_type   offs,
      real_type & x,
      real_type & y,
      real_type & s,
      real_type & dst
    ) const;

    static integer   m_max_iter;
    static real_type m_tolerance;

    mutable bool               m_aabb_done{false};
    mutable AABB_TREE          m_aabb_tree;
    mutable real_type          m_aabb_offs{real_type(0)};
    mutable real_type          m_aabb_max_angle{real_type(0)};
    mutable real_type          m_aabb_max_size{real_type(0)};
    mutable vector<Triangle2D> m_aabb_triangles;

    #ifdef CLOTHOIDS_USE_THREADS
    mutable std::mutex m_aabb_mutex;
    #endif

    bool
    aabb_intersect_ISO(
      Triangle2D    const & T1,
      real_type             offs,
      ClothoidCurve const * pC,
      Triangle2D    const & T2,
      real_type             C_offs,
      real_type           & ss1,
      real_type           & ss2
    ) const;

  public:

    #include "BaseCurve_using.hxx"

    //!
    //! Build an empty clothoid curve
    //!
    ClothoidCurve( string const & name );

    //!
    //! Build a copy of an existing clothoid curve
    //!
    ClothoidCurve( ClothoidCurve const & s );

    void setup( GenericContainer const & gc ) override;

    //!
    //! Construct a clothoid with the standard parameters.
    //!
    //! \param[in] x0     starting position \f$x\f$-coordinate
    //! \param[in] y0     starting position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] k      curvature
    //! \param[in] dk     curvature derivative
    //! \param[in] L      length
    //! \param[in] name   name of the clothoid curve
    //!
    explicit
    ClothoidCurve(
      real_type      x0,
      real_type      y0,
      real_type      theta0,
      real_type      k,
      real_type      dk,
      real_type      L,
      string const & name
    );

    //!
    //! Construct a clothoid \f$ G(s) \f$ solving the \f$ G^1 \f$ problem.
    //!
    //! \f[
    //!   \begin{array}{rcl}
    //!     G(0)  &=& \mathbf{p}_0 \\[1em]
    //!     G(L)  &=& \mathbf{p}_1 \\[1em]
    //!     G'(0) &=& (\cos\theta_0,\sin\theta_0)^T \\[1em]
    //!     G'(L) &=& (\cos\theta_1,\sin\theta_1)^T \\[1em]
    //!   \end{array}
    //! \f]
    //!
    //! \param[in] P0     initial point \f$ \mathbf{p}_0 \f$
    //! \param[in] theta0 initial angle \f$ \theta_0 \f$
    //! \param[in] P1     final point \f$ \mathbf{p}_1 \f$
    //! \param[in] theta1 final angle \f$ \theta_1 \f$
    //! \param[in] name   name of the clothoid curve
    //!
    explicit
    ClothoidCurve(
      real_type const   P0[],
      real_type         theta0,
      real_type const   P1[],
      real_type         theta1,
      string    const & name
    );

    //!
    //! Build a clothoid copying an existing one.
    //!
    void copy( ClothoidCurve const & c );

    //!
    //! Build a clothoid copying an existing line segment.
    //!
    explicit
    ClothoidCurve( LineSegment const & LS );

    //!
    //! Build a clothoid copying an existing circle arc.
    //!
    explicit
    ClothoidCurve( CircleArc const & C );

    //!
    //! Build a clothoid copying an existing curve.
    //!
    explicit
    ClothoidCurve( BaseCurve const * pC );

    //!
    //! Copy an existing clothoid.
    //!
    ClothoidCurve const & operator = ( ClothoidCurve const & s )
    { this->copy(s); return *this; }

    CurveType type() const override { return CurveType::CLOTHOID; }

    /*\
     |  _         _ _    _
     | | |__ _  _(_) |__| |
     | | '_ \ || | | / _` |
     | |_.__/\_,_|_|_\__,_|
    \*/
    //!
    //! Build a clothoid with the standard parameters
    //!
    //! \param[in] x0     starting position \f$x\f$-coordinate
    //! \param[in] y0     starting position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] k      curvature
    //! \param[in] dk     curvature derivative
    //! \param[in] L      length
    //!
    void
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type k,
      real_type dk,
      real_type L
    );

    //!
    //! Build a clothoid by solving the hermite \f$ G^1 \f$ problem.
    //!
    //! \param[in] x0     initial x position \f$ x_0      \f$
    //! \param[in] y0     initial y position \f$ y_0      \f$
    //! \param[in] theta0 initial angle      \f$ \theta_0 \f$
    //! \param[in] x1     final x position   \f$ x_1      \f$
    //! \param[in] y1     final y position   \f$ y_1      \f$
    //! \param[in] theta1 final angle        \f$ \theta_1 \f$
    //! \param[in] tol    tolerance
    //! \return number of iteration performed
    //!
    int
    build_G1(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type tol = 1e-12
    );

    //!
    //! Build a clothoid by solving the hermite \f$ G^1 \f$ problem.
    //!
    //! \param[in]  x0     initial x position \f$ x_0      \f$
    //! \param[in]  y0     initial y position \f$ y_0      \f$
    //! \param[in]  theta0 initial angle      \f$ \theta_0 \f$
    //! \param[in]  x1     final x position   \f$ x_1      \f$
    //! \param[in]  y1     final y position   \f$ y_1      \f$
    //! \param[in]  theta1 final angle        \f$ \theta_1 \f$
    //! \param[out] L_D    derivative of the length \f$ L(\theta_0,\theta_1) \f$
    //! \param[out] k_D    derivative of the curvature \f$ \kappa(\theta_0,\theta_1) \f$
    //! \param[out] dk_D   derivative of the curvature variation \f$ \kappa'(\theta_0,\theta_1) \f$
    //! \param[out] tol = \f$10^{-12}\f$
    //! \return number of iteration performed
    //!
    int
    build_G1_D(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type L_D[2],
      real_type k_D[2],
      real_type dk_D[2],
      real_type tol = 1e-12
    );

    //!
    //! Build a clothoid by solving the forward problem.
    //!
    //! \param[in] x0     initial \f$x\f$-position \f$ x_0 \f$
    //! \param[in] y0     initial \f$y\f$-position \f$ y_0 \f$
    //! \param[in] theta0 initial angle \f$ \theta_0 \f$
    //! \param[in] kappa0 initial curvature \f$ \kappa_0 \f$
    //! \param[in] x1     final \f$x\f$-position \f$ x_1 \f$
    //! \param[in] y1     final \f$y\f$-position \f$ y_1 \f$
    //! \param[in] tol    tolerance of the forward problem
    //!
    bool
    build_forward(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type kappa0,
      real_type x1,
      real_type y1,
      real_type tol = 1e-12
    );

    //!
    //! Build a clothoid from a line segment.
    //!
    void build( LineSegment const & LS );

    //!
    //! Build a clothoid from a circle arc.
    //!
    void build( CircleArc const & );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    //!
    //! Return the point at infinity of the clothoids \f$ P(s) \f$.
    //!
    //! \param[out] x    \f$x\f$-coordinate of the \f$\infty\f$ point
    //! \param[out] y    \f$y\f$-coordinate of the \f$\infty\f$ point
    //! \param[out] plus it true return \f$ \lim_{s\to+\infty} P(s) \f$
    //!                  otherwise return \f$ \lim_{s\to-\infty} P(s) \f$
    //!
    //! @html_image{Pinfinity.png,width=60%}
    //!
    void
    Pinfinity( real_type & x, real_type & y, bool plus = true ) const
    { m_CD.Pinfinity( x, y, plus ); }

    //!
    //! Derivative of the curvature of the clothoid.
    //!
    real_type dkappa() const { return m_CD.m_dk; }

    //!
    //! Clothoid curve total variation of the angle.
    //!
    real_type theta_total_variation() const;

    //!
    //! Max and min angle of the curve.
    //!
    real_type theta_min_max( real_type & thMin, real_type & thMax ) const;

    //!
    //! Clothoid angle range.
    //!
    real_type
    delta_theta() const
    { real_type thMin, thMax; return theta_min_max( thMin, thMax ); }

    //!
    //! Max and min of the curvatire of the clothoid curve.
    //!
    real_type curvature_min_max( real_type & kMin, real_type & kMax ) const;

    //!
    //! Clothoid total curvature variation.
    //!
    real_type curvature_total_variation() const;

    //!
    //! Given the clothoid curve \f$ P(s) \f$ compute.
    //!
    //! \f[
    //!    \int_0^L |P''(s)|^2 \mathrm{d}s
    //! \f]
    //!
    real_type integral_curvature2() const;

    //!
    //! Given the clothoid curve \f$ P(s) \f$ compute.
    //!
    //! \f[
    //!    \int_0^L |P'''(s)|^2 \mathrm{d}s
    //! \f]
    //!
    real_type integral_jerk2() const;

    //!
    //! Given the clothoid curve \f$ P(s) \f$ compute.
    //!
    //! \f[
    //!    \int_0^L |P''''(s)|^2 \mathrm{d}s
    //! \f]
    //!
    real_type integral_snap2() const;

    //!
    //! Return a vector of optimized sample parameters for plotting.
    //!
    //! \param offs      offset of the sampled curve
    //! \param npts      suggested minimum number of sampled points
    //! \param max_angle maximum angle variation between two sampled points
    //! \param s         vector of computed parameters
    //!
    void
    optimized_sample_ISO(
      real_type           offs,
      integer             npts,
      real_type           max_angle,
      vector<real_type> & s
    ) const;

    //!
    //! Return a vector of optimized sample parameters for plotting.
    //!
    //! \param offs      offset of the sampled curve
    //! \param npts      suggested minimum number of sampled points
    //! \param max_angle maximum angle variation between two sampled points
    //! \param s         vector of computed parameters
    //!
    void
    optimized_sample_SAE(
      real_type           offs,
      integer             npts,
      real_type           max_angle,
      vector<real_type> & s
    ) const {
      optimized_sample_ISO( -offs, npts, max_angle, s );
    }

    /*\
     |     _ _    _
     |  __| (_)__| |_ __ _ _ _  __ ___
     | / _` | (_-<  _/ _` | ' \/ _/ -_)
     | \__,_|_/__/\__\__,_|_||_\__\___|
    \*/
    //!
    //! Compute the point on clothoid at minimal distance from a given point.
    //!
    //! \param  ds sampling step
    //! \param  qx \f$ x \f$-coordinate of the given point
    //! \param  qy \f$ y \f$-coordinate of the given point
    //! \param  X  \f$ x \f$-coordinate of the point on clothoid at minimal distance
    //! \param  Y  \f$ y \f$-coordinate of the point on clothoid at minimal distance
    //! \param  S  curvilinear coordinate of the point (X,Y) on the clothoid
    //! \return the distance of the point from the clothoid
    //!
    real_type
    closest_point_by_sample(
      real_type   ds,
      real_type   qx,
      real_type   qy,
      real_type & X,
      real_type & Y,
      real_type & S
    ) const;

    //!
    //! Approximate the point on clothoid at minimal distance from a given point
    //! using simple sampling.
    //!
    //! \param  ds sampling step
    //! \param  qx \f$ x \f$-coordinate of the given point
    //! \param  qy \f$ y \f$-coordinate of the given point
    //! \param  S  curvilinear coordinate of the point (X,Y) on the clothoid
    //! \return the distance of the point from the clothoid
    //!
    real_type
    distance_by_sample(
      real_type   ds,
      real_type   qx,
      real_type   qy,
      real_type & S
    ) const {
      real_type X, Y;
      return closest_point_by_sample( ds, qx, qy, X, Y, S );
    }

    //!
    //! Approximate the point on clothoid at minimal distance from a given point
    //! using simple sampling.
    //!
    //! \param  ds sampling step
    //! \param  qx \f$x\f$-coordinate of the given point
    //! \param  qy \f$y\f$-coordinate of the given point
    //! \return the distance of the point from the clothoid
    //!
    real_type
    distance_by_sample(
      real_type ds,
      real_type qx,
      real_type qy
    ) const {
      real_type X, Y, S;
      return closest_point_by_sample( ds, qx, qy, X, Y, S );
    }

    /*\
     |  _    _   _____    _                _
     | | |__| |_|_   _| _(_)__ _ _ _  __ _| |___
     | | '_ \ '_ \| || '_| / _` | ' \/ _` | / -_)
     | |_.__/_.__/|_||_| |_\__,_|_||_\__, |_\___|
     |                               |___/
    \*/

    //!
    //! Get the triangle bounding box
    //! (if angle variation less that \f$ \pi/2 \f$ )
    //!
    bool
    bbTriangle(
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const {
      return m_CD.bbTriangle( m_L, xx0, yy0, xx1, yy1, xx2, yy2 );
    }

    //!
    //! Get the triangle bounding box
    //! (if angle variation less that \f$ \pi/2 \f$)
    //!
    bool
    bbTriangle_ISO(
      real_type   offs,
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const {
      return m_CD.bbTriangle_ISO( m_L, offs, xx0, yy0, xx1, yy1, xx2, yy2 );
    }

    //!
    //! Get the triangle bounding box
    //! (if angle variation less that \f$ \pi/2 \f$)
    //!
    bool
    bbTriangle_SAE(
      real_type   offs,
      real_type & xx0, real_type & yy0,
      real_type & xx1, real_type & yy1,
      real_type & xx2, real_type & yy2
    ) const {
      return m_CD.bbTriangle_SAE( m_L, offs, xx0, yy0, xx1, yy1, xx2, yy2 );
    }

    bool
    bbTriangle( Triangle2D & t, integer icurve = 0 ) const {
      real_type x0, y0, x1, y1, x2, y2;
      bool ok = m_CD.bbTriangle( m_L, x0, y0, x1, y1, x2, y2 );
      if ( ok ) t.build( x0, y0, x1, y1, x2, y2, 0, 0, icurve );
      return ok;
    }

    bool
    bbTriangle_ISO( real_type offs, Triangle2D & t, integer icurve = 0 ) const {
      real_type x0, y0, x1, y1, x2, y2;
      bool ok = m_CD.bbTriangle_ISO( m_L, offs, x0, y0, x1, y1, x2, y2 );
      if ( ok ) t.build( x0, y0, x1, y1, x2, y2, 0, 0, icurve );
      return ok;
    }

    bool
    bbTriangle_SAE( real_type offs, Triangle2D & t, integer icurve = 0 ) const {
      real_type x0, y0, x1, y1, x2, y2;
      bool ok = m_CD.bbTriangle_SAE( m_L, offs, x0, y0, x1, y1, x2, y2 );
      if ( ok ) t.build( x0, y0, x1, y1, x2, y2, 0, 0, icurve );
      return ok;
    }

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

    void
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/6, // 30 degree
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      this->bb_triangles_ISO( 0, tvec, max_angle, max_size, icurve );
    }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

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

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type length() const override { return m_L; }

    real_type length_ISO( real_type ) const override;

    real_type theta_begin()  const override { return m_CD.m_theta0; }
    real_type kappa_begin()  const override { return m_CD.m_kappa0; }
    real_type x_begin()      const override { return m_CD.m_x0; }
    real_type x_end()        const override { return m_CD.X(m_L); }
    real_type y_begin()      const override { return m_CD.m_y0; }
    real_type y_end()        const override { return m_CD.Y(m_L); }
    real_type tx_begin()     const override { return m_CD.tg0_x(); }
    real_type ty_begin()     const override { return m_CD.tg0_y(); }
    real_type nx_begin_ISO() const override { return m_CD.nor0_x_ISO(); }
    real_type ny_begin_ISO() const override { return m_CD.nor0_y_ISO(); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type s ) const override { return m_CD.tg_x( s ); }
    real_type ty    ( real_type s ) const override { return m_CD.tg_y( s ); }
    real_type tx_D  ( real_type s ) const override { return m_CD.tg_x_D( s ); }
    real_type ty_D  ( real_type s ) const override { return m_CD.tg_y_D( s ); }
    real_type tx_DD ( real_type s ) const override { return m_CD.tg_x_DD( s ); }
    real_type ty_DD ( real_type s ) const override { return m_CD.tg_y_DD( s ); }
    real_type tx_DDD( real_type s ) const override { return m_CD.tg_x_DDD( s ); }
    real_type ty_DDD( real_type s ) const override { return m_CD.tg_y_DDD( s ); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg(
      real_type   s,
      real_type & tx,
      real_type & ty
    ) const override
    { m_CD.tg( s, tx, ty ); }

    void
    tg_D(
      real_type   s,
      real_type & tx_D,
      real_type & ty_D
    ) const override
    { m_CD.tg_D( s, tx_D, ty_D ); }

    void
    tg_DD(
      real_type   s,
      real_type & tx_DD,
      real_type & ty_DD
    ) const override
    { m_CD.tg_DD( s, tx_DD, ty_DD ); }

    void
    tg_DDD(
      real_type   s,
      real_type & tx_DDD,
      real_type & ty_DDD
    ) const override
    { m_CD.tg_DDD( s, tx_DDD, ty_DDD ); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    //!
    //! Get clothoid angle at curvilinear cooordinate `s`.
    //!
    //! \param  s curvilinear cooordinate
    //! \return angle (radiant) at curvilinear cooordinate `s`
    //!
    real_type
    theta( real_type s ) const override
    { return m_CD.theta(s); }

    //!
    //! Get clothoid angle derivative (=curvature)
    //! at curvilinear cooordinate `s`.
    //!
    //! \param  s curvilinear cooordinate
    //! \return angle derivative (radiant/s) at curvilinear cooordinate `s`
    //!
    real_type
    theta_D( real_type s ) const override
    { return m_CD.kappa(s); }

    //!
    //! Get clothoid angle second derivative
    //! at curvilinear cooordinate `s`.
    //!
    //! \return angle second derivative (radiant/s^2) at curvilinear cooordinate `s`
    //!
    real_type
    theta_DD( real_type ) const override
    { return m_CD.m_dk; }

    //!
    //! Get clothoid angle third derivative
    //! at curvilinear cooordinate `s`.
    //!
    //! \return angle third derivative (radiant/s^3) at curvilinear cooordinate `s`
    //!
    real_type
    theta_DDD( real_type ) const override
    { return 0; }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    //!
    //! Clothoid \f$x\f$ coordinate at curvilinear coordinate \f$s\f$.
    //!
    //! \param s curvilinear coordinate
    //! \return clothoid \f$x\f$ coordinate
    //!
    real_type X    ( real_type s ) const override { return m_CD.X(s); }
    real_type X_D  ( real_type s ) const override { return m_CD.X_D(s); }
    real_type X_DD ( real_type s ) const override { return m_CD.X_DD(s); }
    real_type X_DDD( real_type s ) const override { return m_CD.X_DDD(s); }

    //!
    //! Clothoid \f$y\f$ coordinate at curvilinear coordinate \f$s\f$.
    //!
    //! \param s curvilinear coordinate
    //! \return clothoid \f$y\f$ coordinate
    //!
    real_type Y    ( real_type s ) const override { return m_CD.Y(s); }
    real_type Y_D  ( real_type s ) const override { return m_CD.Y_D(s); }
    real_type Y_DD ( real_type s ) const override { return m_CD.Y_DD(s); }
    real_type Y_DDD( real_type s ) const override { return m_CD.Y_DDD(s); }

    //!
    //! Clothoid \f$x\f$ coordinate at curvilinear coordinate \f$s\f$.
    //!
    //! \param s    curvilinear coordinate
    //! \param offs lateral offset
    //! \return     clothoid \f$x\f$ coordinate
    //!
    real_type
    X_ISO( real_type s, real_type offs ) const override
    { return m_CD.X_ISO(s,offs); }

    real_type
    X_ISO_D( real_type s, real_type offs ) const override
    { return m_CD.X_ISO_D(s,offs); }

    real_type
    X_ISO_DD( real_type s, real_type offs ) const override
    { return m_CD.X_ISO_DD(s,offs); }

    real_type
    X_ISO_DDD( real_type s, real_type offs ) const override
    { return m_CD.X_ISO_DDD(s,offs); }

    //!
    //! Clothoid \f$y\f$ coordinate at curvilinear coordinate \f$s\f$.
    //!
    //! \param s curvilinear coordinate
    //! \param offs lateral offset
    //! \return clothoid \f$y\f$ coordinate
    //!
    real_type
    Y_ISO( real_type s, real_type offs ) const override
    { return m_CD.Y_ISO(s,offs); }

    real_type
    Y_ISO_D( real_type s, real_type offs ) const override
    { return m_CD.Y_ISO_D(s,offs); }

    real_type
    Y_ISO_DD( real_type s, real_type offs ) const override
    { return m_CD.Y_ISO_DD(s,offs); }

    real_type
    Y_ISO_DDD( real_type s, real_type offs ) const override
    { return m_CD.Y_ISO_DDD(s,offs); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval(
      real_type   s,
      real_type & x,
      real_type & y
    ) const override
    { m_CD.eval( s, x, y ); }

    void
    eval_D(
      real_type   s,
      real_type & x_D,
      real_type & y_D
    ) const override
    { m_CD.eval_D( s, x_D, y_D ); }

    void
    eval_DD(
      real_type   s,
      real_type & x_DD,
      real_type & y_DD
    ) const override
    { m_CD.eval_DD( s, x_DD, y_DD ); }

    void
    eval_DDD(
      real_type   s,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override
    { m_CD.eval_DDD( s, x_DDD, y_DDD ); }

    void
    eval_ISO(
      real_type   s,
      real_type   offs,
      real_type & x,
      real_type & y
    ) const override
    { m_CD.eval_ISO( s, offs, x, y ); }

    void
    eval_ISO_D(
      real_type   s,
      real_type   offs,
      real_type & x_D,
      real_type & y_D
    ) const override
    { m_CD.eval_ISO_D( s, offs, x_D, y_D ); }

    void
    eval_ISO_DD(
      real_type   s,
      real_type   offs,
      real_type & x_DD,
      real_type & y_DD
    ) const override
    { m_CD.eval_ISO_DD( s, offs, x_DD, y_DD ); }

    void
    eval_ISO_DDD(
      real_type   s,
      real_type   offs,
      real_type & x_DDD,
      real_type & y_DDD
    ) const override
    { m_CD.eval_ISO_DDD( s, offs, x_DDD, y_DDD ); }

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void
    translate( real_type tx, real_type ty ) override
    { m_CD.m_x0 += tx; m_CD.m_y0 += ty; }

    void
    rotate( real_type angle, real_type cx, real_type cy ) override
    { m_CD.rotate( angle, cx, cy ); }

    void
    scale( real_type s ) override {
      m_CD.m_kappa0 /= s;
      m_CD.m_dk     /= s*s;
      m_L           *= s;
    }

    void
    reverse() override
    { m_CD.reverse(m_L); }

    void
    change_origin( real_type newx0, real_type newy0 ) override
    { m_CD.m_x0 = newx0; m_CD.m_y0 = newy0; }

    void
    trim( real_type s_begin, real_type s_end ) override {
      m_CD.origin_at( s_begin );
      m_L = s_end - s_begin;
    }

    //!
    //! change the origin of the clothoid at \f$ s_0 \f$
    //! and the length to  \f$ L \f$.
    //!
    //! \param[in] s0   \f$ s_0 \f$
    //! \param[in] newL \f$ L   \f$
    //!
    void
    change_curvilinear_origin( real_type s0, real_type newL ) {
      m_CD.origin_at( s0 );
      m_L = newL;
    }

    /*\
     |        _                     _   ____       _       _
     |    ___| | ___  ___  ___  ___| |_|  _ \ ___ (_)_ __ | |_
     |   / __| |/ _ \/ __|/ _ \/ __| __| |_) / _ \| | '_ \| __|
     |  | (__| | (_) \__ \  __/\__ \ |_|  __/ (_) | | | | | |_
     |   \___|_|\___/|___/\___||___/\__|_|   \___/|_|_| |_|\__|
    \*/
    //!
    //! Compute the point on clothoid at minimal distance from a given point
    //! using the optimized algorithm described in the publication:
    //!
    //! - **E.Bertolazzi, M.Frego**, Point-Clothoid distance and projection computation
    //!   SIAM J. Scientific Computing, Vol. 41, No. 5, pp. A3326-A3353
    //!
    //! \param  qx  \f$x\f$-coordinate of the given point
    //! \param  qy  \f$y\f$-coordinate of the given point
    //! \param  x   \f$x\f$-coordinate of the point on clothoid at minimal distance
    //! \param  y   \f$y\f$-coordinate of the point on clothoid at minimal distance
    //! \param  s   curvilinear coordinate of the point (X,Y) on the clothoid
    //! \param  t   normal coordinate of the point (X,Y) on the clothoid
    //! \param  dst the distance of the point from the clothoid
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

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    void
    build_AABBtree_ISO(
      real_type offs,
      real_type max_angle = Utils::m_pi/18, // 10 degree
      real_type max_size  = 1e100
    ) const;
    #endif

    // collision detection
    bool
    approximate_collision_ISO(
      real_type             offs,
      ClothoidCurve const & c,
      real_type             c_offs,
      real_type             max_angle,
      real_type             max_size
    ) const;

    bool
    collision( ClothoidCurve const & C ) const {
      return collision_ISO( 0, C, 0 );
    }

    bool
    collision_ISO(
      real_type             offs,
      ClothoidCurve const & C,
      real_type             offs_C
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

    void
    intersect_ISO(
      real_type             offs,
      ClothoidCurve const & C,
      real_type             offs_C,
      IntersectList       & ilist
    ) const;

    void
    intersect(
      ClothoidCurve const & C,
      IntersectList       & ilist
    ) const {
      this->intersect_ISO( 0, C, 0, ilist );
    }

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
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidCurve const & c );

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "Clothoid_compatibility.hxx"
#endif

  };

}

///
/// eof: Clothoid.hxx
///
