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
/// file: Biarc.hxx
///

namespace G2lib {

  /*\
   |   ____  _
   |  | __ )(_) __ _ _ __ ___
   |  |  _ \| |/ _` | '__/ __|
   |  | |_) | | (_| | | | (__
   |  |____/|_|\__,_|_|  \___|
  \*/

  //!
  //! Compute biarc fitting by Hemite data
  //!
  //! @html_image{biarc.png,width=60%}
  //!
  class Biarc : public BaseCurve {
    CircleArc m_C0{"Biarc_C0"};
    CircleArc m_C1{"Biarc_C1"};

    void
    gfun( real_type alpha, real_type g[3] ) const {
      real_type so  = sin(alpha);
      real_type co  = cos(alpha);
      real_type oco = alpha*co;
      g[0] = so + oco;
      g[1] = 2*co - alpha*so;
      g[2] = -3*so - oco;
    }

  public:

    #include "BaseCurve_using.hxx"

    //!
    //! Construct and empty biarc
    //!
    Biarc() = delete;
    Biarc( string const & name ) : BaseCurve( name ) {};

    ~Biarc() override = default;

    void setup( GenericContainer const & gc ) override;

    //!
    //! Make a copy of an existing biarc
    //!
    Biarc( Biarc const & ba ) : BaseCurve( ba.name() )
    { this->copy(ba); }

    //!
    //! Construct a biarc passing from the points
    //! \f$(x_0,y_0)\f$ to the point \f$(x_1,y_1)\f$
    //! with initial angle \f$\theta_0\f$ and final angle \f$\theta_1\f$
    //!
    //! \param[in] x0      \f$ x_0      \f$
    //! \param[in] y0      \f$ y_0      \f$
    //! \param[in] theta0  \f$ \theta_0 \f$
    //! \param[in] x1      \f$ x_1      \f$
    //! \param[in] y1      \f$ y_1      \f$
    //! \param[in] theta1  \f$ \theta_1 \f$
    //! \param[in] name    name of the biarc
    //!
    explicit
    Biarc(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      string const & name
    );

    explicit
    Biarc( LineSegment const & LS ) : BaseCurve( LS.name() )
    { this->build( LS ); }

    explicit
    Biarc( CircleArc const & C ) : BaseCurve( C.name() )
    { this->build( C ); }

    explicit
    Biarc( BaseCurve const * pC );

    //!
    //! Make a copy of an existing biarc.
    //!
    void
    copy( Biarc const & c ) {
      m_C0.copy(c.m_C0);
      m_C1.copy(c.m_C1);
    }

    CurveType type() const override { return CurveType::BIARC; }

    //!
    //! Make a copy of an existing biarc.
    //!
    Biarc const & operator = ( Biarc const & ba )
    { this->copy(ba); return *this; }

    //!
    //! Return the first circle arc of the biarc.
    //!
    CircleArc const & C0() const { return m_C0; }

    //!
    //! Return the first circle arc of the biarc.
    //!
    CircleArc const & C1() const { return m_C1; }

    //!
    //! Construct a biarc passing from the points
    //! \f$ (x_0,y_0) \f$ to the point \f$ (x_1,y_1) \f$
    //! with initial angle \f$ \theta_0 \f$ and final angle \f$ \theta_1 \f$
    //!
    //! \param[in] x0      \f$ x_0      \f$
    //! \param[in] y0      \f$ y_0      \f$
    //! \param[in] theta0  \f$ \theta_0 \f$
    //! \param[in] x1      \f$ x_1      \f$
    //! \param[in] y1      \f$ y_1      \f$
    //! \param[in] theta1  \f$ \theta_1 \f$
    //! \return false if biarc cannot be computed
    //!
    bool
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1
    );

    //!
    //! Construct a biarc by 3 point at "minimum energy"
    //!
    //! - Planar point set fairing and fitting by arc splines
    //! - Xunnian Yang and Guozhao Wang
    //! - Computer-Aided Design, vol 33, 2001
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
    //! Build a biarc from a line segment
    //!
    void build( LineSegment const & LS );

    //!
    //! Build a biarc from a circle arc.
    //!
    void build( CircleArc const & C );
    void build( Biarc const & );
    void build( ClothoidCurve const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

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
    { return m_C0.length()+m_C1.length(); }

    real_type
    length_ISO( real_type offs ) const override
    { return m_C0.length_ISO(offs)+m_C1.length_ISO(offs); }

    real_type theta_begin()  const override { return m_C0.theta_begin(); }
    real_type theta_end()    const override { return m_C1.theta_end(); }
    real_type kappa_begin()  const override { return m_C0.kappa_begin(); }
    real_type kappa_end()    const override { return m_C1.kappa_end(); }
    real_type x_begin()      const override { return m_C0.x_begin(); }
    real_type x_end()        const override { return m_C1.x_end(); }
    real_type y_begin()      const override { return m_C0.y_begin(); }
    real_type y_end()        const override { return m_C1.y_end(); }
    real_type tx_begin()     const override { return m_C0.tx_begin(); }
    real_type tx_end()       const override { return m_C1.tx_end(); }
    real_type ty_begin()     const override { return m_C0.ty_begin(); }
    real_type ty_end()       const override { return m_C1.ty_end(); }
    real_type nx_begin_ISO() const override { return m_C0.nx_begin_ISO(); }
    real_type nx_end_ISO()   const override { return m_C1.nx_end_ISO(); }
    real_type ny_begin_ISO() const override { return m_C0.ny_begin_ISO(); }
    real_type ny_end_ISO()   const override { return m_C1.ny_end_ISO(); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type theta    ( real_type ) const override;
    real_type theta_D  ( real_type ) const override;
    real_type theta_DD ( real_type ) const override { return 0; }
    real_type theta_DDD( real_type ) const override { return 0; }

    void
    evaluate(
      real_type   s,
      real_type & th,
      real_type & k,
      real_type & x,
      real_type & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X( real_type s ) const override;
    real_type Y( real_type s ) const override;

    real_type X_D( real_type ) const override;
    real_type Y_D( real_type ) const override;

    real_type X_DD( real_type ) const override;
    real_type Y_DD( real_type ) const override;

    real_type X_DDD( real_type ) const override;
    real_type Y_DDD( real_type ) const override;

    real_type X_ISO( real_type s, real_type offs ) const override;
    real_type Y_ISO( real_type s, real_type offs ) const override;

    real_type X_ISO_D( real_type, real_type offs ) const override;
    real_type Y_ISO_D( real_type, real_type offs ) const override;

    real_type X_ISO_DD( real_type, real_type offs ) const override;
    real_type Y_ISO_DD( real_type, real_type offs ) const override;

    real_type X_ISO_DDD( real_type, real_type offs ) const override;
    real_type Y_ISO_DDD( real_type, real_type offs ) const override;

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
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type s ) const override;
    real_type tx_D  ( real_type s ) const override;
    real_type tx_DD ( real_type s ) const override;
    real_type tx_DDD( real_type s ) const override;
    real_type ty    ( real_type s ) const override;
    real_type ty_D  ( real_type s ) const override;
    real_type ty_DD ( real_type s ) const override;
    real_type ty_DDD( real_type s ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg(
      real_type   s,
      real_type & tx,
      real_type & ty
    ) const override;

    void
    tg_D(
      real_type   s,
      real_type & tx_D,
      real_type & ty_D
    ) const override;

    void
    tg_DD(
      real_type   s,
      real_type & tx_DD,
      real_type & ty_DD
    ) const override;

    void
    tg_DDD(
      real_type   s,
      real_type & tx_DDD,
      real_type & ty_DDD
    ) const override;

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void
    translate( real_type tx, real_type ty ) override
    { m_C0.translate(tx,ty); m_C1.translate(tx,ty); }

    void
    rotate( real_type angle, real_type cx, real_type cy ) override
    { m_C0.rotate(angle,cx,cy); m_C1.rotate(angle,cx,cy); }

    void
    reverse() override;

    void
    change_origin( real_type newx0, real_type newy0 ) override;

    void
    trim( real_type s_begin, real_type s_end ) override;

    void
    scale( real_type s ) override;

    /*\
     |        _                     _   ____       _       _
     |    ___| | ___  ___  ___  ___| |_|  _ \ ___ (_)_ __ | |_
     |   / __| |/ _ \/ __|/ _ \/ __| __| |_) / _ \| | '_ \| __|
     |  | (__| | (_) \__ \  __/\__ \ |_|  __/ (_) | | | | | |_
     |   \___|_|\___/|___/\___||___/\__|_|   \___/|_|_| |_|\__|
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

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    //!
    //! Return the \f$x\f$-coordinate of the juction point of the biarc.
    //!
    real_type x_middle() const { return m_C1.x_begin(); }

    //!
    //! Return the \f$y\f$-coordinate of the juction point of the biarc.
    //!
    real_type y_middle() const { return m_C1.y_begin(); }

    //!
    //! Return the angle of the juction point of the biarc.
    //!
    real_type theta_middle() const { return m_C1.theta_begin(); }

    //!
    //! Return the curvature of the first arc of the biarc.
    //!
    real_type kappa0() const { return m_C0.curvature(); }

    //!
    //! Intial angle of the biarc.
    //!
    real_type length0() const { return m_C0.length(); }

    //!
    //! Return the curvature of the first arc of the biarc.
    //!
    real_type kappa1() const { return m_C1.curvature(); }

    //!
    //! Final angle of the biarc.
    //!
    real_type length1() const { return m_C1.length(); }

    //!
    //! Change of the angle of the biarc (\f$ \theta_1-\theta_0 \f$).
    //!
    real_type delta_theta() const { return m_C0.delta_theta() + m_C1.delta_theta(); }

    void
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bb_triangles( tvec, max_angle, max_size, icurve );
      m_C1.bb_triangles( tvec, max_angle, max_size, icurve );
    }

    void
    bb_triangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bb_triangles_ISO( offs, tvec, max_angle, max_size, icurve );
      m_C1.bb_triangles_ISO( offs, tvec, max_angle, max_size, icurve );
    }

    void
    bb_triangles_SAE(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bb_triangles_SAE( offs, tvec, max_angle, max_size, icurve );
      m_C1.bb_triangles_SAE( offs, tvec, max_angle, max_size, icurve );
    }

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    //!
    //! Detect a collision with another biarc.
    //!
    bool
    collision( Biarc const & B ) const {
      return m_C0.collision( B.m_C0 ) || m_C0.collision( B.m_C1 ) ||
             m_C1.collision( B.m_C0 ) || m_C1.collision( B.m_C1 );
    }

    //!
    //! Detect a collision with another biarc with offset.
    //!
    //! \param[in] offs   offset of first biarc
    //! \param[in] B      second biarc
    //! \param[in] offs_B offset of second biarc
    //!
    bool
    collision_ISO(
      real_type     offs,
      Biarc const & B,
      real_type     offs_B
    ) const {
      return m_C0.collision_ISO( offs, B.m_C0, offs_B ) ||
             m_C0.collision_ISO( offs, B.m_C1, offs_B ) ||
             m_C1.collision_ISO( offs, B.m_C0, offs_B ) ||
             m_C1.collision_ISO( offs, B.m_C1, offs_B );
    }

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
    //! Intersect a biarc with another biarc.
    //!
    //! \param[in]  B     second biarc
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      Biarc const   & B,
      IntersectList & ilist
    ) const;

    //!
    //! Intersect a biarc with another biarc with offset (ISO).
    //!
    //! \param[in]  offs   offset of first biarc
    //! \param[in]  B      second biarc
    //! \param[in]  offs_B offset of second biarc
    //! \param[out] ilist  list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type       offs,
      Biarc const   & B,
      real_type       offs_B,
      IntersectList & ilist
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
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, Biarc const & bi );

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "Biarc_compatibility.hxx"
#endif

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Given a list of points \f$ (x_i,y_i) \f$
  //! build a guess of angles for a spline of biarc.
  //!
  //! \param[in]  n     number of points
  //! \param[in]  x     \f$x\f$-coordinates
  //! \param[in]  y     \f$y\f$-coordinates
  //! \param[out] theta guessed angles
  //!
  bool
  build_guess_theta(
    integer         n,
    real_type const x[],
    real_type const y[],
    real_type       theta[]
  );

}

///
/// eof: Biarc.hxx
///
