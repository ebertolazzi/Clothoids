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
/// file: Dubins3p.hxx
///

namespace G2lib {

  /*
  //   ____        _     _           _____
  //  |  _ \ _   _| |__ (_)_ __  ___|___ / _ __
  //  | | | | | | | '_ \| | '_ \/ __| |_ \| '_ \
  //  | |_| | |_| | |_) | | | | \__ \___) | |_) |
  //  |____/ \__,_|_.__/|_|_| |_|___/____/| .__/
  //                                      |_|
  */

  //!
  //! Type of Dubins for three points construction algorithm
  //!
  //! - SAMPLE_ONE_DEGREE
  //!   search solution by sampling a fixed angles
  //!
  //! - PATTERN_SEARCH
  //!   search solution by using pattern search
  //!
  //! - PATTERN_TRICHOTOMY
  //!   search solution by using pattern search based on tricotomy
  //!
  //! - PATTERN_SEARCH_WITH_ALGO748
  //!   search solution by using pattern search and refinement using ALGO748
  //!
  //! - PATTERN_TRICHOTOMY_WITH_ALGO748
  //!   search solution by using pattern search with tricotomy and refinement using ALGO748
  //!
  //! - ELLIPSE
  //!   search solution by using ellipse geometric construction
  //!
  //! - POLYNOMIAL_SYSTEM
  //!   search solution by solving polinomial system
  //!
  using Dubins3pBuildType = enum class Dubins3pBuildType : integer {
    SAMPLE_ONE_DEGREE,
    PATTERN_SEARCH,
    PATTERN_TRICHOTOMY,
    PATTERN_SEARCH_WITH_ALGO748,
    PATTERN_TRICHOTOMY_WITH_ALGO748,
    ELLIPSE,
    POLYNOMIAL_SYSTEM
  };

  //!
  //!  \param[in] str name of the `Dubins3pBuildType` type
  //!  \return the `Dubins3pBuildType` enumerator
  //!
  Dubins3pBuildType string_to_Dubins3pBuildType( string const & str );

  //!
  //! Class to manage a circle arc
  //!
  class Dubins3p : public BaseCurve {
  private:

    Dubins m_Dubins0{"Dubins0"};
    Dubins m_Dubins1{"Dubins1"};

    real_type m_tolerance{Utils::m_pi/180}; // one degree
    real_type m_sample_angle{10*Utils::m_pi/180}; // ten degree
    real_type m_sample_points{360};
    integer   m_max_evaluation{1000};
    integer   m_evaluation;

    bool
    build_sample(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max
    );

    bool
    build_pattern_search(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max,
      real_type tolerance      = 1e-8,
      bool      use_trichotomy = true,
      bool      use_748        = true
    );

    bool
    build_poly_system(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max
    );

    bool
    build_ellipse(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max
    );

  public:

    //!
    //! Build an empty circle
    //!
    Dubins3p() = delete;
    Dubins3p( string const & name ) : BaseCurve( name ) {};

    void setup( GenericContainer const & gc ) override;

    //!
    //! Build a copy of an existing Dubins problem.
    //!
    Dubins3p( Dubins3p const & s ) : BaseCurve( s.name() )
    { this->copy(s); }

    //!
    //! Construct a Dubins3p solution
    //!
    //! \param[in] xi     initial position \f$x\f$-coordinate
    //! \param[in] yi     initial position \f$y\f$-coordinate
    //! \param[in] thetai initial angle
    //! \param[in] xm     intermediate position \f$x\f$-coordinate
    //! \param[in] ym     intermediate position \f$y\f$-coordinate
    //! \param[in] xf     final position \f$x\f$-coordinate
    //! \param[in] yf     final position \f$y\f$-coordinate
    //! \param[in] thetaf final angle
    //! \param[in] k_max  max curvature
    //! \param[in] method construction method
    //! \param[in] name   name of the 3 points Dubins object
    //!
    explicit
    Dubins3p(
      real_type         xi,
      real_type         yi,
      real_type         thetai,
      real_type         xm,
      real_type         ym,
      real_type         xf,
      real_type         yf,
      real_type         thetaf,
      real_type         k_max,
      Dubins3pBuildType method,
      string const &    name
    ) : BaseCurve( name ) {
      this->build( xi, yi, thetai, xm, ym, xf, yf, thetaf, k_max, method );
    }

    //!
    //! Make a copy of an existing Dubins solution.
    //!
    void
    copy( Dubins3p const & d3p ) {
      m_Dubins0.copy( d3p.m_Dubins0 );
      m_Dubins1.copy( d3p.m_Dubins1 );
    }

    //!
    //! Construct a Dubins3p solution
    //!
    //! \param[in] xi     initial position \f$x\f$-coordinate
    //! \param[in] yi     initial position \f$y\f$-coordinate
    //! \param[in] thetai initial angle
    //! \param[in] xm     intermediate position \f$x\f$-coordinate
    //! \param[in] ym     intermediate position \f$y\f$-coordinate
    //! \param[in] xf     final position \f$x\f$-coordinate
    //! \param[in] yf     final position \f$y\f$-coordinate
    //! \param[in] thetaf final angle
    //! \param[in] k_max  max curvature
    //! \param[in] method construction method
    //!
    bool
    build(
      real_type         xi,
      real_type         yi,
      real_type         thetai,
      real_type         xm,
      real_type         ym,
      real_type         xf,
      real_type         yf,
      real_type         thetaf,
      real_type         k_max,
      Dubins3pBuildType method
    );

    void set_tolerance( real_type tol );
    void set_sample_angle( real_type ang );
    void set_max_evaluation( integer max_eval );
    void set_sample_points( integer npts );

    real_type  tolerance()          const { return m_tolerance; }
    real_type  sample_angle()       const { return m_sample_angle; }
    DubinsType solution_type0()     const { return m_Dubins0.solution_type(); }
    DubinsType solution_type1()     const { return m_Dubins1.solution_type(); }
    integer    icode()              const { return m_Dubins0.icode()+16*m_Dubins1.icode(); }
    integer    icode0()             const { return m_Dubins0.icode(); }
    integer    icode1()             const { return m_Dubins1.icode(); }
    integer    num_evaluation()     const { return m_evaluation; }
    integer    max_num_evaluation() const { return m_max_evaluation; }

    string solution_type_string() const;
    string solution_type_string_short() const;

    void build( LineSegment const & L );
    void build( CircleArc const & C );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( BiarcList const & );
    void build( PolyLine const & );
    void build( ClothoidList const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    //!
    //! Get possible point of discontinuity for the length
    //!
    //! \param[in]  xi     initial position \f$x\f$-coordinate
    //! \param[in]  yi     initial position \f$y\f$-coordinate
    //! \param[in]  thetai initial angle
    //! \param[in]  xm     intermediate position \f$x\f$-coordinate
    //! \param[in]  ym     intermediate position \f$y\f$-coordinate
    //! \param[in]  xf     final position \f$x\f$-coordinate
    //! \param[in]  yf     final position \f$y\f$-coordinate
    //! \param[in]  thetaf final angle
    //! \param[in]  k_max  max curvature
    //! \param[out] angles angles of possibile dicontinuity
    //! \return     number of point computed
    //!
    integer
    get_range_angles(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max,
      real_type angles[]
    ) const;

    //!
    //! Get sample point for length minimization
    //!
    //! \param[in]  xi        initial position \f$x\f$-coordinate
    //! \param[in]  yi        initial position \f$y\f$-coordinate
    //! \param[in]  thetai    initial angle
    //! \param[in]  xm        intermediate position \f$x\f$-coordinate
    //! \param[in]  ym        intermediate position \f$y\f$-coordinate
    //! \param[in]  xf        final position \f$x\f$-coordinate
    //! \param[in]  yf        final position \f$y\f$-coordinate
    //! \param[in]  thetaf    final angle
    //! \param[in]  k_max     max curvature
    //! \param[in]  tolerance tolerance angle used in point approximaton
    //! \param[out] angles    sample points
    //!
    void
    get_sample_angles(
      real_type xi,
      real_type yi,
      real_type thetai,
      real_type xm,
      real_type ym,
      real_type xf,
      real_type yf,
      real_type thetaf,
      real_type k_max,
      real_type tolerance,
      vector<real_type> & angles
    ) const;

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

    //! Return the first cicle of the Dubins solution
    CircleArc const & C0() const { return m_Dubins0.m_C0; }
    //! Return the second cicle of the Dubins solution
    CircleArc const & C1() const { return m_Dubins0.m_C1; }
    //! Return the third cicle of the Dubins solution
    CircleArc const & C2() const { return m_Dubins0.m_C2; }
    //! Return the first cicle of the Dubins solution
    CircleArc const & C3() const { return m_Dubins1.m_C0; }
    //! Return the seco4d cicle of the Dubins solution
    CircleArc const & C4() const { return m_Dubins1.m_C1; }
    //! Return the third cicle of the Dubins solution
    CircleArc const & C5() const { return m_Dubins1.m_C2; }

    void
    get_solution( ClothoidList & CL ) const {
      CL.init();
      CL.reserve(6);
      CL.push_back( m_Dubins0.m_C0 );
      CL.push_back( m_Dubins0.m_C1 );
      CL.push_back( m_Dubins0.m_C2 );
      CL.push_back( m_Dubins1.m_C0 );
      CL.push_back( m_Dubins1.m_C1 );
      CL.push_back( m_Dubins1.m_C2 );
    }

    real_type length() const override;
    real_type length_ISO( real_type offs ) const override;

    real_type length0() const { return m_Dubins0.m_C0.length(); }
    real_type length1() const { return m_Dubins0.m_C1.length(); }
    real_type length2() const { return m_Dubins0.m_C2.length(); }
    real_type length3() const { return m_Dubins1.m_C0.length(); }
    real_type length4() const { return m_Dubins1.m_C1.length(); }
    real_type length5() const { return m_Dubins1.m_C2.length(); }

    real_type kappa0() const { return m_Dubins0.m_C0.kappa_begin(); }
    real_type kappa1() const { return m_Dubins0.m_C1.kappa_begin(); }
    real_type kappa2() const { return m_Dubins0.m_C2.kappa_begin(); }
    real_type kappa3() const { return m_Dubins1.m_C0.kappa_begin(); }
    real_type kappa4() const { return m_Dubins1.m_C1.kappa_begin(); }
    real_type kappa5() const { return m_Dubins1.m_C2.kappa_begin(); }

    real_type X0( real_type s ) const { return m_Dubins0.m_C0.X(s); }
    real_type Y0( real_type s ) const { return m_Dubins0.m_C0.Y(s); }

    real_type X1( real_type s ) const { return m_Dubins0.m_C1.X(s); }
    real_type Y1( real_type s ) const { return m_Dubins0.m_C1.Y(s); }

    real_type X2( real_type s ) const { return m_Dubins0.m_C2.X(s); }
    real_type Y2( real_type s ) const { return m_Dubins0.m_C2.Y(s); }

    real_type X3( real_type s ) const { return m_Dubins1.m_C0.X(s); }
    real_type Y3( real_type s ) const { return m_Dubins1.m_C0.Y(s); }

    real_type X4( real_type s ) const { return m_Dubins1.m_C1.X(s); }
    real_type Y4( real_type s ) const { return m_Dubins1.m_C1.Y(s); }

    real_type X5( real_type s ) const { return m_Dubins1.m_C2.X(s); }
    real_type Y5( real_type s ) const { return m_Dubins1.m_C2.Y(s); }

    real_type theta0( real_type s ) const { return m_Dubins0.m_C0.theta(s); }
    real_type theta1( real_type s ) const { return m_Dubins0.m_C1.theta(s); }
    real_type theta2( real_type s ) const { return m_Dubins0.m_C2.theta(s); }
    real_type theta3( real_type s ) const { return m_Dubins1.m_C0.theta(s); }
    real_type theta4( real_type s ) const { return m_Dubins1.m_C1.theta(s); }
    real_type theta5( real_type s ) const { return m_Dubins1.m_C2.theta(s); }

    real_type theta_begin()  const override { return m_Dubins0.m_C0.theta_begin(); }
    real_type theta_end()    const override { return m_Dubins1.m_C2.theta_end(); }

    real_type kappa_begin()  const override { return m_Dubins0.m_C0.kappa_begin(); }
    real_type kappa_end()    const override { return m_Dubins1.m_C2.kappa_end(); }

    real_type x_begin()      const override { return m_Dubins0.m_C0.x_begin(); }
    real_type x_end()        const override { return m_Dubins1.m_C2.x_end(); }

    real_type y_begin()      const override { return m_Dubins0.m_C0.y_begin(); }
    real_type y_end()        const override { return m_Dubins1.m_C2.y_end(); }

    real_type tx_begin()     const override { return m_Dubins0.m_C0.tx_begin(); }
    real_type tx_end()       const override { return m_Dubins1.m_C2.tx_end(); }

    real_type ty_begin()     const override { return m_Dubins0.m_C0.ty_begin(); }
    real_type ty_end()       const override { return m_Dubins1.m_C2.ty_end(); }

    real_type nx_begin_ISO() const override { return m_Dubins0.m_C0.nx_begin_ISO(); }
    real_type nx_end_ISO()   const override { return m_Dubins1.m_C2.nx_end_ISO(); }

    real_type ny_begin_ISO() const override { return m_Dubins0.m_C0.ny_begin_ISO(); }
    real_type ny_end_ISO()   const override { return m_Dubins1.m_C2.ny_end_ISO(); }

    real_type theta0_begin() const { return m_Dubins0.m_C0.theta_begin(); }
    real_type theta0_end()   const { return m_Dubins0.m_C0.theta_end(); }

    real_type theta1_begin() const { return m_Dubins0.m_C1.theta_begin(); }
    real_type theta1_end()   const { return m_Dubins0.m_C1.theta_end(); }

    real_type theta2_begin() const { return m_Dubins0.m_C2.theta_begin(); }
    real_type theta2_end()   const { return m_Dubins0.m_C2.theta_end(); }

    real_type theta3_begin() const { return m_Dubins1.m_C0.theta_begin(); }
    real_type theta3_end()   const { return m_Dubins1.m_C0.theta_end(); }

    real_type theta4_begin() const { return m_Dubins1.m_C1.theta_begin(); }
    real_type theta4_end()   const { return m_Dubins1.m_C1.theta_end(); }

    real_type theta5_begin() const { return m_Dubins1.m_C2.theta_begin(); }
    real_type theta5_end()   const { return m_Dubins1.m_C2.theta_end(); }

    real_type x0_begin() const { return m_Dubins0.m_C0.x_begin(); }
    real_type x0_end()   const { return m_Dubins0.m_C0.x_end(); }

    real_type y0_begin() const { return m_Dubins0.m_C0.y_begin(); }
    real_type y0_end()   const { return m_Dubins0.m_C0.y_end(); }

    real_type x1_begin() const { return m_Dubins0.m_C1.x_begin(); }
    real_type x1_end()   const { return m_Dubins0.m_C1.x_end(); }

    real_type y1_begin() const { return m_Dubins0.m_C1.y_begin(); }
    real_type y1_end()   const { return m_Dubins0.m_C1.y_end(); }

    real_type x2_begin() const { return m_Dubins0.m_C2.x_begin(); }
    real_type x2_end()   const { return m_Dubins0.m_C2.x_end(); }

    real_type y2_begin() const { return m_Dubins0.m_C2.y_begin(); }
    real_type y2_end()   const { return m_Dubins0.m_C2.y_end(); }

    real_type x3_begin() const { return m_Dubins1.m_C0.x_begin(); }
    real_type x3_end()   const { return m_Dubins1.m_C0.x_end(); }

    real_type y3_begin() const { return m_Dubins1.m_C0.y_begin(); }
    real_type y3_end()   const { return m_Dubins1.m_C0.y_end(); }

    real_type x4_begin() const { return m_Dubins1.m_C1.x_begin(); }
    real_type x4_end()   const { return m_Dubins1.m_C1.x_end(); }

    real_type y4_begin() const { return m_Dubins1.m_C1.y_begin(); }
    real_type y4_end()   const { return m_Dubins1.m_C1.y_end(); }

    real_type x5_begin() const { return m_Dubins1.m_C2.x_begin(); }
    real_type x5_end()   const { return m_Dubins1.m_C2.x_end(); }

    real_type y5_begin() const { return m_Dubins1.m_C2.y_begin(); }
    real_type y5_end()   const { return m_Dubins1.m_C2.y_end(); }

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
    translate( real_type tx, real_type ty ) override {
      m_Dubins0.translate(tx,ty);
      m_Dubins1.translate(tx,ty);
    }

    void
    rotate( real_type angle, real_type cx, real_type cy ) override {
      m_Dubins0.rotate(angle,cx,cy);
      m_Dubins1.rotate(angle,cx,cy);
    }

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
      m_Dubins0.bb_triangles( tvec, max_angle, max_size, icurve );
      m_Dubins1.bb_triangles( tvec, max_angle, max_size, icurve );
    }

    void
    bb_triangles_ISO(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_Dubins0.bb_triangles_ISO( offs, tvec, max_angle, max_size, icurve );
      m_Dubins1.bb_triangles_ISO( offs, tvec, max_angle, max_size, icurve );
    }

    void
    bb_triangles_SAE(
      real_type            offs,
      vector<Triangle2D> & tvec,
      real_type            max_angle = Utils::m_pi/18,
      real_type            max_size  = 1e100,
      integer              icurve    = 0
    ) const override {
      m_Dubins0.bb_triangles_SAE( offs, tvec, max_angle, max_size, icurve );
      m_Dubins1.bb_triangles_SAE( offs, tvec, max_angle, max_size, icurve );
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
    bool collision( Dubins3p const & B ) const;

    //!
    //! Detect a collision with another biarc with offset.
    //!
    //! \param[in] offs   offset of first biarc
    //! \param[in] B      second biarc
    //! \param[in] offs_B offset of second biarc
    //!
    bool
    collision_ISO(
      real_type        offs,
      Dubins3p const & B,
      real_type        offs_B
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
      Dubins3p const & B,
      IntersectList  & ilist
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
      real_type        offs,
      Dubins3p const & B,
      real_type        offs_B,
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
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, Dubins3p const & bi );

    CurveType type() const override { return CurveType::DUBINS; }

  };

}

///
/// eof: Dubins3p.hxx
///
