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
/// file: Dubins.hxx
///

namespace G2lib {

  /*\
   |   ____        _     _
   |  |  _ \ _   _| |__ (_)_ __  ___
   |  | | | | | | | '_ \| | '_ \/ __|
   |  | |_| | |_| | |_) | | | | \__ \
   |  |____/ \__,_|_.__/|_|_| |_|___/
  \*/

  //!
  //! Class to manage a circle arc
  //!
  class Dubins : public BaseCurve {
  public:
    using DubinsType = enum class DubinsType : integer
    { LSL, RSR, LSR, RSL, LRL, RLR };

  private:
    DubinsType m_solution_type;

    CircleArc m_C0, m_C1, m_C2; //! Three arc solution of DUBINS problem

  public:

    //!
    //! Build an empty circle
    //!
    Dubins() = default;

    //!
    //! Build a copy of an existing Dubins problem.
    //!
    Dubins( Dubins const & s )
    { this->copy(s); }

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in] x0     initial position x-coordinate
    //! \param[in] y0     initial position y-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position x-coordinate
    //! \param[in] y1     final position y-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] kmax   max curvature
    //!
    explicit
    Dubins(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type k_max
    ) {
      this->build( x0, y0, theta0, x1, y1, theta1, k_max );
    }

    //!
    //! Make a copy of an existing Dubins solution.
    //!
    void
    copy( Dubins const & d ) {
      m_C0 = d.m_C0;
      m_C1 = d.m_C1;
      m_C2 = d.m_C2;
      m_solution_type = d.m_solution_type;
    }

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in] x0     initial position x-coordinate
    //! \param[in] y0     initial position y-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position x-coordinate
    //! \param[in] y1     final position y-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] kmax   max curvature
    //!
    bool
    build(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type k_max
    );

    void build( LineSegment const & L );
    void build( CircleArc const & C );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( BiarcList const & );
    void build( PolyLine const & );
    void build( ClothoidList const & );
    void build( Dubins const & );

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

    //!
    //! Return the first cicle of the Dubins solution
    //!
    CircleArc const & C0() const { return m_C0; }

    //!
    //! Return the second cicle of the Dubins solution
    //!
    CircleArc const & C1() const { return m_C1; }

    //!
    //! Return the third cicle of the Dubins solution
    //!
    CircleArc const & C2() const { return m_C2; }

    void
    get_solution( ClothoidList & CL ) const {
      CL.init();
      CL.reserve(3);
      CL.push_back( m_C0 );
      CL.push_back( m_C1 );
      CL.push_back( m_C2 );
    }

    real_type length() const override;
    real_type length_ISO( real_type offs ) const override;

    DubinsType solution_type() const { return m_solution_type; }

    real_type length0() const { return m_C0.length(); }
    real_type length1() const { return m_C1.length(); }
    real_type length2() const { return m_C2.length(); }

    real_type kappa0() const { return m_C0.kappa_begin(); }
    real_type kappa1() const { return m_C1.kappa_begin(); }
    real_type kappa2() const { return m_C2.kappa_begin(); }

    real_type X0( real_type s ) const { return m_C0.X(s); }
    real_type Y0( real_type s ) const { return m_C0.Y(s); }

    real_type X1( real_type s ) const { return m_C1.X(s); }
    real_type Y1( real_type s ) const { return m_C1.Y(s); }

    real_type X2( real_type s ) const { return m_C2.X(s); }
    real_type Y2( real_type s ) const { return m_C2.Y(s); }

    real_type theta0( real_type s ) const { return m_C0.theta(s); }
    real_type theta1( real_type s ) const { return m_C1.theta(s); }
    real_type theta2( real_type s ) const { return m_C2.theta(s); }

    real_type theta_begin()  const override { return m_C0.theta_begin(); }
    real_type theta_end()    const override { return m_C2.theta_end(); }
    real_type kappa_begin()  const override { return m_C0.kappa_begin(); }
    real_type kappa_end()    const override { return m_C2.kappa_end(); }
    real_type x_begin()      const override { return m_C0.x_begin(); }
    real_type y_begin()      const override { return m_C0.y_begin(); }
    real_type x_end()        const override { return m_C2.x_end(); }
    real_type y_end()        const override { return m_C2.y_end(); }
    real_type tx_begin()     const override { return m_C0.tx_begin(); }
    real_type ty_begin()     const override { return m_C0.ty_begin(); }
    real_type nx_begin_ISO() const override { return m_C0.nx_begin_ISO(); }
    real_type ny_begin_ISO() const override { return m_C0.ny_begin_ISO(); }

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type theta    ( real_type s ) const override;
    real_type theta_D  ( real_type   ) const override;
    real_type theta_DD ( real_type   ) const override { return 0; }
    real_type theta_DDD( real_type   ) const override { return 0; }

    real_type X( real_type s ) const override;
    real_type Y( real_type s ) const override;

    real_type X_D( real_type s ) const override;
    real_type Y_D( real_type s ) const override;

    real_type X_DD( real_type s ) const override;
    real_type Y_DD( real_type s ) const override;

    real_type X_DDD( real_type s ) const override;
    real_type Y_DDD( real_type s ) const override;

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void
    translate( real_type tx, real_type ty ) override
    { m_C0.translate(tx,ty); m_C1.translate(tx,ty); m_C2.translate(tx,ty); }

    void
    rotate( real_type angle, real_type cx, real_type cy ) override
    { m_C0.rotate(angle,cx,cy); m_C1.rotate(angle,cx,cy); m_C2.rotate(angle,cx,cy); }

    void reverse() override;

    void change_origin( real_type newx0, real_type newy0 ) override;

    void
    trim( real_type, real_type ) override {
      UTILS_ERROR0( "Dubins::trim not defined, convert to ClothoidList to trim the curve!");
    }

    void scale( real_type s ) override;

    void
    eval(
      real_type   s,
      real_type & theta,
      real_type & kappa,
      real_type & x,
      real_type & y
    ) const;

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

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    bbTriangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bbTriangles( tvec, max_angle, max_size, icurve );
      m_C1.bbTriangles( tvec, max_angle, max_size, icurve );
      m_C2.bbTriangles( tvec, max_angle, max_size, icurve );
    }

    void
    bbTriangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bbTriangles_ISO( offs, tvec, max_angle, max_size, icurve );
      m_C1.bbTriangles_ISO( offs, tvec, max_angle, max_size, icurve );
      m_C2.bbTriangles_ISO( offs, tvec, max_angle, max_size, icurve );
    }

    void
    bbTriangles_SAE(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bbTriangles_SAE( offs, tvec, max_angle, max_size, icurve );
      m_C1.bbTriangles_SAE( offs, tvec, max_angle, max_size, icurve );
      m_C2.bbTriangles_SAE( offs, tvec, max_angle, max_size, icurve );
    }

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
    bool collision( Dubins const & B ) const;

    //!
    //! Detect a collision with another biarc with offset.
    //!
    //! \param[in] offs   offset of first biarc
    //! \param[in] B      second biarc
    //! \param[in] offs_B offset of second biarc
    //!
    bool
    collision_ISO(
      real_type      offs,
      Dubins const & B,
      real_type      offs_B
    ) const;

    bool collision( BaseCurve const * pC ) const override;

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
      Dubins const  & B,
      IntersectList & ilist
    ) const;

    //!
    //! Intersect a biarc with another biarc with offset (ISO).
    //!
    //! \param[in]  offs   offset of first biarc
    //! \param[in]  B      second Dubins
    //! \param[in]  offs_B offset of second biarc
    //! \param[out] ilist  list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type       offs,
      Dubins const  & B,
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

    string
    info() const
    { return fmt::format( "Dubins\n{}\n", *this ); }

    void
    info( ostream_type & stream ) const override
    { stream << this->info(); }

    //!
    //! Pretty print of Dubins class.
    //!
    friend
    ostream_type &
    operator << ( ostream_type & stream, Dubins const & bi );

    CurveType type() const override { return CurveType::DUBINS; }

  };

  inline
  string
  to_string( Dubins::DubinsType n ) {
    string res = "";
    switch ( n ) {
    case Dubins::DubinsType::LSL: res = "LSL"; break;
    case Dubins::DubinsType::RSR: res = "RSR"; break;
    case Dubins::DubinsType::LSR: res = "LSR"; break;
    case Dubins::DubinsType::RSL: res = "RSL"; break;
    case Dubins::DubinsType::LRL: res = "LRL"; break;
    case Dubins::DubinsType::RLR: res = "RLR"; break;
    }
    return res;
  };
}

///
/// eof: Dubins.hxx
///
