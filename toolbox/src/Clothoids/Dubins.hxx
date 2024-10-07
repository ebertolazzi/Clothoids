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
  //! Type of Dubins solution
  //!
  //! - LSL left-stright-left solution
  //! - RSR right-stright-right solution
  //! - LSR left-stright-right solution
  //! - RSL right-stright-left solution
  //! - LRL left-right-left solution
  //! - RLR right-left-right solution
  //!
  using DubinsType = enum class DubinsType : integer
  { LSL, RSR, LSR, RSL, LRL, RLR, DUBINS_ERROR };

  //!
  //! Convert Dubins type solution to an integer.
  //!
  //! \param d type of curve
  //! \return a number from `0` to `5`
  //!
  integer to_integer( DubinsType d );

  bool
  Dubins_build(
    real_type    x0,
    real_type    y0,
    real_type    theta0,
    real_type    x1,
    real_type    y1,
    real_type    theta1,
    real_type    k_max,
    DubinsType & type,
    real_type  & L1,
    real_type  & L2,
    real_type  & L3,
    real_type    grad[2]
  );

  //!
  //! Class to manage a Dubins curve
  //!
  class Dubins : public BaseCurve {
  private:
    DubinsType m_solution_type{DubinsType::DUBINS_ERROR};
    real_type  m_length{0};
    real_type  m_length_Dalpha{0};
    real_type  m_length_Dbeta{0};

    CircleArc m_C0{"Dubins_C0"}; //! Three arc solution of DUBINS problem
    CircleArc m_C1{"Dubins_C1"}; //! Three arc solution of DUBINS problem
    CircleArc m_C2{"Dubins_C2"}; //! Three arc solution of DUBINS problem

  public:

    //!
    //! Build an empty circle
    //!
    Dubins() = delete;
    Dubins( string const & name ) : BaseCurve( name ) {};

    void setup( GenericContainer const & gc ) override;

    //!
    //! Build a copy of an existing Dubins problem.
    //!
    Dubins( Dubins const & s ) : BaseCurve( s.name() )
    { this->copy(s); }

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in] x0     initial position \f$x\f$-coordinate
    //! \param[in] y0     initial position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position \f$x\f$-coordinate
    //! \param[in] y1     final position \f$y\f$-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] k_max  max curvature
    //! \param[in] name   name of the Dubins object
    //!
    explicit
    Dubins(
      real_type      x0,
      real_type      y0,
      real_type      theta0,
      real_type      x1,
      real_type      y1,
      real_type      theta1,
      real_type      k_max,
      string const & name
    ) : BaseCurve( name ) {
      this->build( x0, y0, theta0, x1, y1, theta1, k_max );
    }

    //!
    //! Make a copy of an existing Dubins solution.
    //!
    void
    copy( Dubins const & d ) {
      m_C0            = d.m_C0;
      m_C1            = d.m_C1;
      m_C2            = d.m_C2;
      m_length        = d.m_length;
      m_length_Dalpha = d.m_length_Dalpha;
      m_length_Dbeta  = d.m_length_Dbeta;
      m_solution_type = d.m_solution_type;
    }

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in] x0     initial position \f$x\f$-coordinate
    //! \param[in] y0     initial position \f$y\f$-coordinate
    //! \param[in] theta0 initial angle
    //! \param[in] x1     final position \f$x\f$-coordinate
    //! \param[in] y1     final position \f$y\f$-coordinate
    //! \param[in] theta1 final angle
    //! \param[in] k_max  max curvature
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

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in]  x0     initial position \f$x\f$-coordinate
    //! \param[in]  y0     initial position \f$y\f$-coordinate
    //! \param[in]  x1     final position \f$x\f$-coordinate
    //! \param[in]  y1     final position \f$y\f$-coordinate
    //! \param[in]  theta1 final angle
    //! \param[in]  k_max  max curvature
    //! \param[out] angles range points
    //! \return     number of range points
    //!
    integer
    get_range_angles_begin(
      real_type x0,
      real_type y0,
      real_type x1,
      real_type y1,
      real_type theta1,
      real_type k_max,
      real_type angles[]
    ) const;

    //!
    //! Construct a Dubins solution
    //!
    //! \param[in]  x0     initial position \f$x\f$-coordinate
    //! \param[in]  y0     initial position \f$y\f$-coordinate
    //! \param[in]  theta0 initial angle
    //! \param[in]  x1     final position \f$x\f$-coordinate
    //! \param[in]  y1     final position \f$y\f$-coordinate
    //! \param[in]  k_max  max curvature
    //! \param[out] angles range points
    //! \return     number of range points
    //!
    integer
    get_range_angles_end(
      real_type x0,
      real_type y0,
      real_type theta0,
      real_type x1,
      real_type y1,
      real_type k_max,
      real_type angles[]
    ) const;

    void build( LineSegment const & L );
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

    real_type length()        const override { return m_length; }
    real_type length_Dalpha() const          { return m_length_Dalpha; }
    real_type length_Dbeta()  const          { return m_length_Dbeta; }
    real_type length_ISO( real_type offs ) const override;

    void
    length_grad( real_type grad[2] ) const {
      grad[0] = m_length_Dalpha;
      grad[1] = m_length_Dbeta;
    }

    DubinsType solution_type() const { return m_solution_type; }

    string solution_type_string() const;
    string solution_type_string_short() const;

    integer icode() const { return to_integer(m_solution_type); }

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
    real_type x_end()        const override { return m_C2.x_end(); }

    real_type y_begin()      const override { return m_C0.y_begin(); }
    real_type y_end()        const override { return m_C2.y_end(); }

    real_type tx_begin()     const override { return m_C0.tx_begin(); }
    real_type tx_end()       const override { return m_C2.tx_end(); }

    real_type ty_begin()     const override { return m_C0.ty_begin(); }
    real_type ty_end()       const override { return m_C2.ty_end(); }

    real_type nx_begin_ISO() const override { return m_C0.nx_begin_ISO(); }
    real_type nx_end_ISO()   const override { return m_C2.nx_end_ISO(); }

    real_type ny_begin_ISO() const override { return m_C0.ny_begin_ISO(); }
    real_type ny_end_ISO()   const override { return m_C2.ny_end_ISO(); }

    real_type theta0_begin() const { return m_C0.theta_begin(); }
    real_type theta0_end()   const { return m_C0.theta_end(); }

    real_type theta1_begin() const { return m_C1.theta_begin(); }
    real_type theta1_end()   const { return m_C1.theta_end(); }

    real_type theta2_begin() const { return m_C2.theta_begin(); }
    real_type theta2_end()   const { return m_C2.theta_end(); }

    real_type x0_begin() const { return m_C0.x_begin(); }
    real_type x0_end()   const { return m_C0.x_end(); }

    real_type y0_begin() const { return m_C0.y_begin(); }
    real_type y0_end()   const { return m_C0.y_end(); }

    real_type x1_begin() const { return m_C1.x_begin(); }
    real_type x1_end()   const { return m_C1.x_end(); }

    real_type y1_begin() const { return m_C1.y_begin(); }
    real_type y1_end()   const { return m_C1.y_end(); }

    real_type x2_begin() const { return m_C2.x_begin(); }
    real_type x2_end()   const { return m_C2.x_end(); }

    real_type y2_begin() const { return m_C2.y_begin(); }
    real_type y2_end()   const { return m_C2.y_end(); }

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

    void trim( real_type, real_type ) override;

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
    bb_triangles(
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_C0.bb_triangles( tvec, max_angle, max_size, icurve );
      m_C1.bb_triangles( tvec, max_angle, max_size, icurve );
      m_C2.bb_triangles( tvec, max_angle, max_size, icurve );
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
      m_C2.bb_triangles_ISO( offs, tvec, max_angle, max_size, icurve );
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
      m_C2.bb_triangles_SAE( offs, tvec, max_angle, max_size, icurve );
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

    string info() const;

    void
    info( ostream_type & stream ) const override
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, Dubins const & bi );

    CurveType type() const override { return CurveType::DUBINS; }

    friend class Dubins3p;
  };

  //!
  //! \param[in] n Dubins type solution
  //! \return  the string with the name of the solution
  //!
  inline
  string
  to_string( DubinsType n ) {
    string res{""};
    switch ( n ) {
    case DubinsType::LSL:          res = "LSL";   break;
    case DubinsType::RSR:          res = "RSR";   break;
    case DubinsType::LSR:          res = "LSR";   break;
    case DubinsType::RSL:          res = "RSL";   break;
    case DubinsType::LRL:          res = "LRL";   break;
    case DubinsType::RLR:          res = "RLR";   break;
    case DubinsType::DUBINS_ERROR: res = "ERROR"; break;
    }
    return res;
  };

}

///
/// eof: Dubins.hxx
///
