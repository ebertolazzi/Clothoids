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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: ClothoidList.hxx
///

namespace G2lib {

  using std::vector;

  //!
  //! Type of Spline list to be build
  //!
  using TargetType = enum class TargetType : integer
  { P1, P2, P3, P4, P5, P6, P7, P8, P9 };

  //!
  //! Convert `TargetType` to string
  //!
  //! \param[in] n the target type
  //! \return a string with the name of the target
  //!
  static
  inline
  string_view
  to_string( TargetType n ) {
    switch ( n ) {
    case TargetType::P1: return "P1";
    case TargetType::P2: return "P2";
    case TargetType::P3: return "P3";
    case TargetType::P4: return "P4";
    case TargetType::P5: return "P5";
    case TargetType::P6: return "P6";
    case TargetType::P7: return "P7";
    case TargetType::P8: return "P8";
    case TargetType::P9: return "P9";
    }
    return "";
  };
}

#include "ClothoidList_G2solve2arc.hxx"
#include "ClothoidList_G2solveCLC.hxx"
#include "ClothoidList_G2solve3arc.hxx"
#include "ClothoidList_ClothoidSplineG2.hxx"

namespace G2lib {

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/
  //!
  //! Manage a piecewise clothoids \f$ G(s) \f$ composed by
  //! \f$ n \f$ clothoids (not necessarily \f$ G^2 \f$ or \f$ G^1 \f$ connected)
  //!
  //! @html_image{G2problem3arc.png,width=60%}
  //!
  class ClothoidList : public BaseCurve {

    bool                  m_curve_is_closed{false};
    vector<real_type>     m_s0;
    vector<ClothoidCurve> m_clothoid_list;

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
      real_type const qx,
      real_type const qy,
      real_type const offs,
      real_type      & x,
      real_type      & y,
      real_type      & s,
      real_type      & DST
    ) const;

    // mette a posto gli angoli delle clotoidi in modo che differiscano meno di 2*pi!
    void adjust_angles();

  public:

    #include "BaseCurve_using.hxx"

    ClothoidList() = delete;

    //!
    //! Build an empty clothoid list
    //!
    explicit
    ClothoidList( string_view name ) : BaseCurve( name )
    { this->reset_last_interval(); }

    ~ClothoidList() override {
      m_s0.clear();
      m_clothoid_list.clear();
      m_aabb_triangles.clear();
    }

    void setup( GenericContainer const & gc ) override;

    //!
    //! Build a copy of an existing clothoid list
    //!
    ClothoidList( ClothoidList const & s ) : BaseCurve( s.name() )
    { this->reset_last_interval(); this->copy(s); }

    CurveType type() const override { return CurveType::CLOTHOID_LIST; }

    //!
    //! Initialize the clothoid list
    //!
    void init();

    //!
    //! Reserve memory for `n` clothoid
    //!
    void reserve( integer const n );

    //!
    //! Build a clothoid list copying an existing one
    //!
    void copy( ClothoidList const & L );
    
    //!
    //! Return the internal list of clothoids
    //!
    vector<ClothoidCurve> const & get_list() const { return m_clothoid_list; }

    //!
    //! Copy an existing clothoid list
    //!
    ClothoidList const & operator = ( ClothoidList const & s )
    { this->copy(s); return *this; }

    //!
    //! Build a clothoid from a line segment
    //!
    explicit ClothoidList( LineSegment const & LS );

    //!
    //! Build a clothoid from a circle arc
    //!
    explicit ClothoidList( CircleArc const & C );

    //!
    //! Build a clothoid from a biarc
    //!
    explicit ClothoidList( Biarc const & B );

    //!
    //! Build a clothoid from a list of biarc
    //!
    explicit ClothoidList( BiarcList const & BL );

    //!
    //! Build a clothoid from a clothoid curve
    //!
    explicit ClothoidList( ClothoidCurve const & CL );

    //!
    //! Build a clothoid from a list line segment
    //!
    explicit ClothoidList( PolyLine const & PL );

    //!
    //! Build a clothoid from G2solve2arc
    //!
    explicit ClothoidList( G2solve2arc const & C, string_view name );

    //!
    //! Build a clothoid from G2solve3arc
    //!
    explicit ClothoidList( G2solve3arc const & C, string_view name  );

    //!
    //! Build a clothoid from G2solveCLC
    //!
    explicit ClothoidList( G2solveCLC const & C, string_view name  );

    //!
    //! Build a clothoid from a curve
    //!
    explicit ClothoidList( BaseCurve const * pC );

    void build( LineSegment const & );
    void build( CircleArc const & );
    void build( ClothoidCurve const & );
    void build( Biarc const & );
    void build( PolyLine const & );
    void build( BiarcList const & );
    void build( ClothoidList const & );
    void build( G2solve2arc const & );
    void build( G2solve3arc const & );
    void build( G2solveCLC const & );
    void build( Dubins const & );
    void build( Dubins3p const & );

    //!
    //! Add a line segment to the tail of clothoid list
    //!
    void push_back( LineSegment const & c );

    //!
    //! Add a circle arc to the tail of clothoid list
    //!
    void push_back( CircleArc const & c );

    //!
    //! Add a biarc to the tail of clothoid list
    //!
    void push_back( Biarc const & c );

    //!
    //! Add a biarc list to the tail of clothoid list
    //!
    void push_back( BiarcList const & c );

    //!
    //! Add a clothoid curve to the tail of clothoid list
    //!
    void push_back( ClothoidCurve const & c );

    //!
    //! Add a clothoid list to the tail of clothoid list
    //!
    void push_back( ClothoidList const & c );

    //!
    //! Add a dubins 3 arc curve to the tail of clothoid list
    //!
    void push_back( Dubins const & c );

    //!
    //! Add a dubins 6 arc curve to the tail of clothoid list
    //!
    void push_back( Dubins3p const & c );

    //!
    //! Add a G2solve2arc to the tail of clothoid list
    //!
    void push_back( G2solve2arc const & c );

    //!
    //! Add a G2solve3arc list to the tail of clothoid list
    //!
    void push_back( G2solve3arc const & c );

    //!
    //! Add a clothoid list to the tail of clothoid list
    //!
    void push_back( G2solveCLC const & c );

    //!
    //! Add a list of line segment to the tail of clothoid list
    //!
    void push_back( PolyLine const & c );

    //!
    //! Add a clothoid to the tail of the clothoid list.
    //!
    //! \param kappa0 initial curvature
    //! \param dkappa derivative of the curvature
    //! \param L      length of the segment
    //!
    void
    push_back(
      real_type kappa0,
      real_type dkappa,
      real_type L
    );

    //!
    //! Add a clothoid to the tail of the clothoid list.
    //! The builded clothoid is translated to the tail of the clothioid list.
    //!
    //! \param x0     initial x
    //! \param y0     initial y
    //! \param theta0 initial angle
    //! \param kappa0 initial curvature
    //! \param dkappa derivative of the curvature
    //! \param L      length of the segment
    //!
    void
    push_back(
      real_type const x0,
      real_type const y0,
      real_type const theta0,
      real_type const kappa0,
      real_type const dkappa,
      real_type const L
    );

    //!
    //! Add a clothoid to the tail of the clothoid list solving the \f$ G^1 \f$ problem.
    //! The initial point and angle are taken from the tail of the clothoid list.
    //!
    //! \param x1     final x
    //! \param y1     final y
    //! \param theta1 final angle
    //!
    void
    push_back_G1(
      real_type const x1,
      real_type const y1,
      real_type const theta1
    );

    //!
    //! Add a clothoid to the tail of the clothoid list solving the \f$ G^1 \f$ problem.
    //! The initial point and angle are taken from the tail of the clothoid list.
    //! The builded clothoid is translated to the tail of the clothioid list.
    //!
    //! \param x0     initial x
    //! \param y0     initial y
    //! \param theta0 initial angle
    //! \param x1     final x
    //! \param y1     final y
    //! \param theta1 final angle
    //!
    void
    push_back_G1(
      real_type const x0,
      real_type const y0,
      real_type const theta0,
      real_type const x1,
      real_type const y1,
      real_type const theta1
    );

    //!
    //! True if curve is closed
    //!
    bool is_closed() const { return m_curve_is_closed; }

    //!
    //! Set clousure flag to true
    //!
    void make_closed() { m_curve_is_closed = true; }

    //!
    //! Set clousure flag to false
    //!
    void make_open() { m_curve_is_closed = false; }

    //!
    //! Difference initial final point x component
    //!
    real_type closure_gap_x() const { return this->x_end() - this->x_begin(); }

    //!
    //! Difference initial final point y component
    //!
    real_type closure_gap_y() const { return this->y_end() - this->y_begin(); }

    //!
    //! Difference initial final tangent x component
    //!
    real_type closure_gap_tx() const { return this->tx_end() - this->tx_begin(); }

    //!
    //! Difference initial final tangent y component
    //!
    real_type closure_gap_ty() const { return this->ty_end() - this->ty_begin(); }

    //!
    //! check if clothoid list is closed
    //!
    //! \param[in] tol_xy position tolerance
    //! \param[in] tol_tg angle (tangent) tolerance
    //!
    //! \return true if curve is closed
    //!
    bool
    closure_check( real_type const tol_xy = 1e-6, real_type const tol_tg = 1e-6 ) const {
      return std::abs(closure_gap_x())  < tol_xy &&
             std::abs(closure_gap_y())  < tol_xy &&
             std::abs(closure_gap_tx()) < tol_tg &&
             std::abs(closure_gap_ty()) < tol_tg;
    }

    //!
    //! Build clothoid list passing to a list of points
    //! solving a series of \f$ G^1 \f$ fitting problems.
    //! The angle at points are estimated using the routine `xy_to_guess_angle`
    //!
    //! \param[in] n number of points
    //! \param[in] x \f$x\f$-coordinates
    //! \param[in] y \f$y\f$-coordinates
    //!
    //! \return false if routine fails
    //!
    bool
    build_G1(
      integer   const n,
      real_type const x[],
      real_type const y[]
    );

    //!
    //! Build clothoid list passing to a list of points
    //! solving a series of \f$ G^1 \f$ fitting problems.
    //!
    //! \param[in] n     number of points
    //! \param[in] x     \f$x\f$-coordinates
    //! \param[in] y     \f$y\f$-coordinates
    //! \param[in] theta angles at the points
    //!
    //! \return false if routine fails
    //!
    bool
    build_G1(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const theta[]
    );

    //!
    //! \brief
    //! Builds a sequence of clothoid (Cornu) arcs interpolating a given set of
    //! points with \f$ G^2 \f$ (curvature) continuity between segments.
    //!
    //! \details
    //! Constructs a smooth open curve passing through the specified points,
    //! starting with a given orientation and curvature, and ending with the
    //! prescribed final orientation and curvature. The resulting clothoid list
    //! ensures continuous position, tangent, and curvature along the entire path.
    //!
    //! \param[in] n           Number of points to interpolate.
    //! \param[in] x           Array of point \f$x_i\f$ coordinates.
    //! \param[in] y           Array of point \f$y_i\f$ coordinates.
    //! \param[in] theta_init  Initial angle \f$\theta_i\f$ (radians).
    //! \param[in] kappa_init  Initial curvature \f$\kappa_i\f$ (1/length units).
    //! \param[in] theta_end   Final angle \f$\theta_e\f$ (radians).
    //! \param[in] kappa_end   Final curvature \f$\kappa_e\f$ (1/length units).
    //!
    //! \return
    //! `true` if the \f$ G^2 \f$ clothoid construction succeeds,
    //! `false` otherwise (e.g., degenerate or inconsistent input data).
    //!
    //! \note
    //! - Arrays `x[]` and `y[]` must have the same length `n`.
    //! - The function does not modify the input arrays.
    //!
    bool
    build_G2(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const theta_init,
      real_type const kappa_init,
      real_type const theta_end,
      real_type const kappa_end
    );

    //!
    //! \brief
    //! Builds a sequence of clothoid (Cornu) arcs interpolating a given set of
    //! points with \f$ G^2 \f$ (curvature) continuity between segments.
    //!
    //! \details
    //! Constructs a smooth open curve passing through the specified points,
    //! starting with a given orientation and curvature, and ending with the
    //! prescribed final orientation and curvature. The resulting clothoid list
    //! ensures continuous position, tangent, and curvature along the entire path.
    //!
    //! \param[in] n           Number of points to interpolate.
    //! \param[in] x           Array of point \f$x_i\f$ coordinates.
    //! \param[in] y           Array of point \f$y_i\f$ coordinates.
    //! \param[in] theta_init  Initial angle \f$\theta_i\f$ (radians).
    //! \param[in] theta_end   Final angle \f$\theta_e\f$ (radians).
    //!
    //! \return
    //! `true` if the \f$ G^2 \f$ clothoid construction succeeds,
    //! `false` otherwise (e.g., degenerate or inconsistent input data).
    //!
    //! \note
    //! - Arrays `x[]` and `y[]` must have the same length `n`.
    //! - The function does not modify the input arrays.
    //!
    bool
    build_G2(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const theta_init,
      real_type const theta_end
    );

    bool
    build_G2_with_target(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const theta[],
      real_type const wL[],
      real_type const wR[],
      real_type const theta_init,
      real_type const theta_end,
      std::function<real_type(ClothoidList const & lst)> const & target
    );

    //!
    //! \brief
    //! Builds a cyclic list of clothoid (Cornu) arcs interpolating a set of points
    //! with \f$ G^2 \f$ geometric continuity.
    //!
    //! \details
    //! Creates a closed sequence of clothoid segments ensuring continuous
    //! position, tangent, and curvature between consecutive points.
    //!
    //! \param[in] n  Number of points (must be ≥ 3 for a closed curve).
    //! \param[in] x  Array of \f$x_i\f$ coordinates.
    //! \param[in] y  Array of \f$y_i\f$ coordinates.
    //!
    //! \return
    //! `true` if the cyclic \f$ G^2 \f$ clothoid construction succeeds,
    //! `false` otherwise.
    //!
    //! \note
    //! The input arrays `x[]` and `y[]` must have the same length and are not modified.
    //!
    bool
    build_G2_cyclic(
      integer   const n,
      real_type const x[],
      real_type const y[]
    );

    bool
    build_G2_cyclic_with_target(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const theta[],
      real_type const wL[],
      real_type const wR[],
      std::function<real_type(ClothoidList const & lst)> const & target
    );

    //!
    //! \brief Smooths a G1 spline of clothoids to minimize jumps in curvature.
    //!
    //! This function applies a smoothing algorithm to reduce discontinuities
    //! in curvature, enhancing the overall continuity and visual quality of the spline.
    //! The smoothing process iterates up to a specified maximum number of iterations.
    //!
    //! \param max_iter Maximum number of iterations for the smoothing process.
    //!                 This parameter controls how many times the smoothing
    //!                 algorithm will attempt to refine the spline.
    //! \param epsi     Tolerance for the smoothing process. It determines the
    //!                 threshold below which curvature adjustments are considered
    //!                 negligible and thus terminates the algorithm.
    //! \param max_dK   Reference to a variable that will hold the maximum change
    //!                 in curvature observed during the smoothing process. This
    //!                 output can be used to assess the effectiveness of the
    //!                 smoothing operation.
    //!
    //! \return Returns true if the smoothing was successful; otherwise, false.
    //!
    //! \note This function assumes that the input spline is already in a G1
    //!       state and is defined using clothoids. Ensure that the input parameters
    //!       are properly initialized before calling this function.
    //!
    //! \note If curve is closed and cyclic algorithm try to reduce the jump in curvature
    //!       between first and last point.
    //!
    bool
    smooth_quasi_G2(
      integer   const max_iter,
      real_type const epsi,
      real_type     & max_dK
    );

    //!
    //! Build clothoid list with \f$ G^2 \f$ continuity.
    //! The vector `s` contains the breakpoints of the curve.
    //! Between two breakpoint the curvature change linearly (is a clothoid)
    //!
    //! \param[in] x0     initial x
    //! \param[in] y0     initial y
    //! \param[in] theta0 initial angle
    //! \param[in] n      number of segments
    //! \param[in] s      break point of the piecewise curve
    //! \param[in] kappa  curvature at the break point
    //!
    //! \return true if curve is closed
    //!
    bool
    build(
      real_type const x0,
      real_type const y0,
      real_type const theta0,
      integer   const n,
      real_type const s[],
      real_type const kappa[]
    );

    //!
    //! Build clothoid list with \f$ G^2 \f$ continuity.
    //! The vector `s` contains the breakpoints of the curve.
    //! Between two breakpoint the curvature change linearly (is a clothoid)
    //!
    //! \param[in] x0     initial x
    //! \param[in] y0     initial y
    //! \param[in] theta0 initial angle
    //! \param[in] s      break point of the piecewise curve
    //! \param[in] kappa  curvature at the break point
    //!
    //! \return true if curve is closed
    //!
    bool
    build(
      real_type         const   x0,
      real_type         const   y0,
      real_type         const   theta0,
      vector<real_type> const & s,
      vector<real_type> const & kappa
    ) {
      if ( s.size() != kappa.size() ) return false;
      return build(
        x0, y0, theta0,
        integer(s.size()),
        &s.front(), &kappa.front()
      );
    }

    //!
    //! Build clothoid listy using raw data.
    //!
    //! \param[in] n        number of points
    //! \param[in] x        \f$x\f$-coordinates
    //! \param[in] y        \f$y\f$-coordinates
    //! \param[in] abscissa break point of the piecewise curve
    //! \param[in] theta    angles at breakpoints
    //! \param[in] kappa    curvature at the break point
    //!
    //! \return false if fails
    //!
    bool
    build_raw(
      integer   const n,
      real_type const x[],
      real_type const y[],
      real_type const abscissa[],
      real_type const theta[],
      real_type const kappa[]
    );

    //!
    //! Build clothoid listy using raw data.
    //!
    //! \param[in] x        \f$x\f$-coordinates
    //! \param[in] y        \f$y\f$-coordinates
    //! \param[in] abscissa break point of the piecewise curve
    //! \param[in] theta    angles at breakpoints
    //! \param[in] kappa    curvature at the break point
    //!
    //! \return false if fails
    //!
    bool
    build_raw(
      vector<real_type> const & x,
      vector<real_type> const & y,
      vector<real_type> const & abscissa,
      vector<real_type> const & theta,
      vector<real_type> const & kappa
    ) {
      integer n = integer(x.size());
      if ( n != integer(y.size())        ||
           n != integer(abscissa.size()) ||
           n != integer(theta.size())    ||
           n != integer(kappa.size()) ) return false;
      return build_raw(
        n, &x.front(), &y.front(),
        &abscissa.front(), &theta.front(), &kappa.front()
      );
    }

    //!
    //! Get the `idx`-th clothoid of the list
    //!
    ClothoidCurve const & get( integer const idx ) const;

    //!
    //! Get the `idx`-th clothoid of the list where `idx` is the clothoid at parameter `s`
    //!
    ClothoidCurve const & get_at_s( real_type const s ) const;

    //!
    //! Return the numbber of clothoid of the list
    //!
    integer num_segments() const { return integer(m_clothoid_list.size()); }

    //!
    //! The list of clothoid has total length \f$ L \f$
    //! the parameter \f$ s \f$ us recomputed as \f$ s+kL \f$ in such a way
    //! \f$ s+kL\in[0,L) \f$ with \f$ k\in\mathbb{Z} \f$.
    //!
    void wrap_in_range( real_type & s ) const;

    //!
    //! Find the clothoid segment whose definiton range contains `s`
    //!
    integer find_at_s( real_type & s ) const;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type length() const override;
    real_type length_ISO( real_type const offs ) const override;

    //!
    //! Return the length of the `nseg`-th clothoid of the list
    //!
    real_type
    segment_length( integer const nseg ) const;

    //!
    //! Return the length of the `nseg`-th clothoid of the list with offset
    //!
    real_type
    segment_length_ISO( integer const nseg, real_type const offs ) const;

    //!
    //! Return the length of the `nseg`-th clothoid of the list with offset
    //!
    real_type
    segment_length_SAE( integer const nseg, real_type const offs ) const
    { return segment_length_ISO( nseg, -offs ); }

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
      real_type const      max_angle = Utils::m_pi/6, // 30 degree
      real_type const      max_size  = 1e100,
      integer   const      icurve    = 0
    ) const override;

    void
    bb_triangles_ISO(
      real_type const      offs,
      vector<Triangle2D> & tvec,
      real_type const      max_angle = Utils::m_pi/6, // 30 degree
      real_type const      max_size  = 1e100,
      integer   const      icurve    = 0
    ) const override;

    void
    bb_triangles_SAE(
      real_type const      offs,
      vector<Triangle2D> & tvec,
      real_type const      max_angle = Utils::m_pi/6, // 30 degree
      real_type const      max_size  = 1e100,
      integer   const      icurve    = 0
    ) const override {
      this->bb_triangles_ISO( -offs, tvec, max_angle, max_size, icurve );
    }

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    void
    build_AABBtree_ISO(
      real_type const offs,
      real_type const max_angle = Utils::m_pi/6, // 30 degree
      real_type const max_size  = 1e100
    ) const;
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
      real_type const offs,
      real_type     & xmin,
      real_type     & ymin,
      real_type     & xmax,
      real_type     & ymax
    ) const override;

    /*\
     |   ____             _          _______           _
     |  | __ )  ___  __ _(_)_ __    / / ____|_ __   __| |
     |  |  _ \ / _ \/ _` | | '_ \  / /|  _| | '_ \ / _` |
     |  | |_) |  __/ (_| | | | | |/ / | |___| | | | (_| |
     |  |____/ \___|\__, |_|_| |_/_/  |_____|_| |_|\__,_|
     |              |___/
    \*/

    real_type theta_begin() const override { return m_clothoid_list.front().theta_begin(); }
    real_type theta_end()   const override { return m_clothoid_list.back().theta_end(); }
    real_type x_begin()     const override { return m_clothoid_list.front().x_begin(); }
    real_type y_begin()     const override { return m_clothoid_list.front().y_begin(); }
    real_type x_end()       const override { return m_clothoid_list.back().x_end(); }
    real_type y_end()       const override { return m_clothoid_list.back().y_end(); }
    real_type x_begin_ISO( real_type const offs ) const override { return m_clothoid_list.front().x_begin_ISO( offs ); }
    real_type y_begin_ISO( real_type const offs ) const override { return m_clothoid_list.front().y_begin_ISO( offs ); }
    real_type x_end_ISO  ( real_type const offs ) const override { return m_clothoid_list.back().x_end_ISO( offs ); }
    real_type y_end_ISO  ( real_type const offs ) const override { return m_clothoid_list.back().y_end_ISO( offs ); }
    real_type tx_begin()     const override { return m_clothoid_list.front().tx_begin(); }
    real_type ty_begin()     const override { return m_clothoid_list.front().ty_begin(); }
    real_type tx_end()       const override { return m_clothoid_list.back().tx_end(); }
    real_type ty_end()       const override { return m_clothoid_list.back().ty_end(); }
    real_type nx_begin_ISO() const override { return m_clothoid_list.front().nx_begin_ISO(); }
    real_type ny_begin_ISO() const override { return m_clothoid_list.front().ny_begin_ISO(); }
    real_type nx_end_ISO()   const override { return m_clothoid_list.back().nx_end_ISO(); }
    real_type ny_end_ISO()   const override { return m_clothoid_list.back().ny_end_ISO(); }

    /*\
     |  _   _          _
     | | |_| |__   ___| |_ __ _
     | | __| '_ \ / _ \ __/ _` |
     | | |_| | | |  __/ || (_| |
     |  \__|_| |_|\___|\__\__,_|
    \*/

    real_type theta    ( real_type const s ) const override;
    real_type theta_D  ( real_type const s ) const override;
    real_type theta_DD ( real_type const s ) const override;
    real_type theta_DDD( real_type const s ) const override;

    G2LIB_DEFINE_1ARG_AUTODIFF( theta )

    /*\
     |  _____                   _   _   _
     | |_   _|   __ _ _ __   __| | | \ | |
     |   | |    / _` | '_ \ / _` | |  \| |
     |   | |   | (_| | | | | (_| | | |\  |
     |   |_|    \__,_|_| |_|\__,_| |_| \_|
    \*/

    real_type tx    ( real_type const s ) const override;
    real_type ty    ( real_type const s ) const override;
    real_type tx_D  ( real_type const s ) const override;
    real_type ty_D  ( real_type const s ) const override;
    real_type tx_DD ( real_type const s ) const override;
    real_type ty_DD ( real_type const s ) const override;
    real_type tx_DDD( real_type const s ) const override;
    real_type ty_DDD( real_type const s ) const override;

    G2LIB_DEFINE_1ARG_AUTODIFF( tx )
    G2LIB_DEFINE_1ARG_AUTODIFF( ty )

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    tg(
      real_type const s,
      real_type     & tg_x,
      real_type     & tg_y
    ) const override;

    void
    tg_D(
      real_type const s,
      real_type     & tg_x_D,
      real_type     & tg_y_D
    ) const override;

    void
    tg_DD(
      real_type const s,
      real_type     & tg_x_DD,
      real_type     & tg_y_DD
    ) const override;

    void
    tg_DDD(
      real_type const s,
      real_type     & tg_x_DDD,
      real_type     & tg_y_DDD
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate(
      real_type const s,
      real_type     & th,
      real_type     & k,
      real_type     & x,
      real_type     & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    evaluate_ISO(
      real_type const s,
      real_type const offs,
      real_type     & th,
      real_type     & k,
      real_type     & x,
      real_type     & y
    ) const override;

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    real_type X    ( real_type const s ) const override;
    real_type Y    ( real_type const s ) const override;
    real_type X_D  ( real_type const s ) const override;
    real_type Y_D  ( real_type const s ) const override;
    real_type X_DD ( real_type const s ) const override;
    real_type Y_DD ( real_type const s ) const override;
    real_type X_DDD( real_type const s ) const override;
    real_type Y_DDD( real_type const s ) const override;

    G2LIB_DEFINE_1ARG_AUTODIFF( X )
    G2LIB_DEFINE_1ARG_AUTODIFF( Y )

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval(
      real_type const s,
      real_type     & x,
      real_type     & y
    ) const override;

    void
    eval_D(
      real_type const s,
      real_type     & x_D,
      real_type     & y_D
    ) const override;

    void
    eval_DD(
      real_type const s,
      real_type     & x_DD,
      real_type     & y_DD
    ) const override;

    void
    eval_DDD(
      real_type const s,
      real_type     & x_DDD,
      real_type     & y_DDD
    ) const override;

    /*\
     |         __  __          _
     |   ___  / _|/ _|___  ___| |_
     |  / _ \| |_| |_/ __|/ _ \ __|
     | | (_) |  _|  _\__ \  __/ |_
     |  \___/|_| |_| |___/\___|\__|
    \*/

    real_type X_ISO    ( real_type const s, real_type const offs ) const override;
    real_type Y_ISO    ( real_type const s, real_type const offs ) const override;
    real_type X_ISO_D  ( real_type const s, real_type const offs ) const override;
    real_type Y_ISO_D  ( real_type const s, real_type const offs ) const override;
    real_type X_ISO_DD ( real_type const s, real_type const offs ) const override;
    real_type Y_ISO_DD ( real_type const s, real_type const offs ) const override;
    real_type X_ISO_DDD( real_type const s, real_type const offs ) const override;
    real_type Y_ISO_DDD( real_type const s, real_type const offs ) const override;

    G2LIB_DEFINE_1ARG_1PAR_AUTODIFF( X_ISO )
    G2LIB_DEFINE_1ARG_1PAR_AUTODIFF( Y_ISO )

    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    void
    eval_ISO(
      real_type const s,
      real_type const offs,
      real_type     & x,
      real_type     & y
    ) const override;

    void
    eval_ISO_D(
      real_type const s,
      real_type const offs,
      real_type     & x_D,
      real_type     & y_D
    ) const override;

    void
    eval_ISO_DD(
      real_type const s,
      real_type const offs,
      real_type     & x_DD,
      real_type     & y_DD
    ) const override;

    void
    eval_ISO_DDD(
      real_type const s,
      real_type const offs,
      real_type     & x_DDD,
      real_type     & y_DDD
    ) const override;

    /*\
     |  _                        __
     | | |_ _ __ __ _ _ __  ___ / _| ___  _ __ _ __ ___
     | | __| '__/ _` | '_ \/ __| |_ / _ \| '__| '_ ` _ \
     | | |_| | | (_| | | | \__ \  _| (_) | |  | | | | | |
     |  \__|_|  \__,_|_| |_|___/_|  \___/|_|  |_| |_| |_|
    \*/

    void translate( real_type const tx, real_type const ty ) override;
    void rotate( real_type const angle, real_type const cx, real_type const cy ) override;
    void scale( real_type const sc ) override;
    void reverse() override;
    void change_origin( real_type const newx0, real_type const newy0 ) override;
    void trim( real_type const s_begin, real_type const s_end ) override;
    void trim( real_type const s_begin, real_type const s_end, ClothoidList & newCL ) const;

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    //!
    //! \param  qx     \f$x\f$-coordinate of the point
    //! \param  qy     \f$y\f$-coordinate of the point
    //! \param  x      \f$x\f$-coordinate of the projected point on the curve
    //! \param  y      \f$y\f$-coordinate of the projected point on the curve
    //! \param  s      parameter on the curve of the projection
    //! \param  t      curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst    distance point projected point
    //! \return n >= 0 point is projected orthogonal, n is the number of the segment at minimum distance<br>
    //!        -(n+1)  minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_ISO(
      real_type const qx,
      real_type const qy,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst
    ) const override;

    //!
    //! \param  qx    \f$x\f$-coordinate of the point
    //! \param  qy    \f$y\f$-coordinate of the point
    //! \param  offs  offset of the curve
    //! \param  x     \f$x\f$-coordinate of the projected point on the curve
    //! \param  y     \f$y\f$-coordinate of the projected point on the curve
    //! \param  s     parameter on the curve of the projection
    //! \param  t     curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst   distance point projected point
    //! \return n > 0 point is projected orthogonal, n-1 is the number of the segment at minimum distance<br>
    //!        -(n+1) minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_ISO(
      real_type const qx,
      real_type const qy,
      real_type const offs,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst
    ) const override;

    //!
    //! Compute the point on clothoid at minimal distance from a given point
    //! using the optimized algorithm described in the publication:
    //!
    //! \param  ds sampling step
    //! \param  qx \f$x\f$-coordinate of the given point
    //! \param  qy \f$y\f$-coordinate of the given point
    //! \param  X  \f$x\f$-coordinate of the point on clothoid at minimal distance
    //! \param  Y  \f$y\f$-coordinate of the point on clothoid at minimal distance
    //! \param  S  curvilinear coordinate of the point (X,Y) on the clothoid
    //! \return the distance of the point from the clothoid
    //!
    real_type
    closest_point_by_sample(
      real_type const ds,
      real_type const qx,
      real_type const qy,
      real_type     & X,
      real_type     & Y,
      real_type     & S
    ) const;

    /*\
     |      _ _     _
     |   __| (_)___| |_ __ _ _ __   ___ ___
     |  / _` | / __| __/ _` | '_ \ / __/ _ \
     | | (_| | \__ \ || (_| | | | | (_|  __/
     |  \__,_|_|___/\__\__,_|_| |_|\___\___|
    \*/

    //!
    //! \param  qx  \f$x\f$-coordinate of the point
    //! \param  qy  \f$y\f$-coordinate of the point
    //! \return the segment at minimal distance from point (qx,qy)
    //!
    integer
    closest_segment( real_type qx, real_type qy ) const;

    //!
    //! \param  qx           \f$x\f$-coordinate of the point
    //! \param  qy           \f$y\f$-coordinate of the point
    //! \param  icurve_begin index of the initial segment
    //! \param  icurve_end   index of the past to the last segment
    //! \param  x            \f$x\f$-coordinate of the projected point on the curve
    //! \param  y            \f$y\f$-coordinate of the projected point on the curve
    //! \param  s            parameter on the curve of the projection
    //! \param  t            curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst          distance point projected point
    //! \param  icurve       number of the segment with the projected point
    //! \return 1            point is projected orthogonal<br>
    //!         0 =          more than one projection (first returned)<br>
    //!        -1 =          minimum point is not othogonal projection to curve
    //!
    integer
    closest_point_in_range_ISO(
      real_type const qx,
      real_type const qy,
      integer   const icurve_begin,
      integer   const icurve_end,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst,
      integer       & icurve
    ) const;

    //!
    //! \param  qx           \f$x\f$-coordinate of the point
    //! \param  qy           \f$y\f$-coordinate of the point
    //! \param  icurve_begin index of the initial segment
    //! \param  icurve_end   index of the past to the last segment
    //! \param  x            \f$x\f$-coordinate of the projected point on the curve
    //! \param  y            \f$y\f$-coordinate of the projected point on the curve
    //! \param  s            parameter on the curve of the projection
    //! \param  t            curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst          distance point projected point
    //! \param  icurve       number of the segment with the projected point
    //! \return 1            point is projected orthogonal<br>
    //!         0            = more than one projection (first returned)<br>
    //!        -1            = minimum point is not othogonal projection to curve<br>
    //!
    integer
    closest_point_in_range_SAE(
      real_type const qx,
      real_type const qy,
      integer   const icurve_begin,
      integer   const icurve_end,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst,
      integer       & icurve
    ) const {
      integer res = this->closest_point_in_range_ISO(
        qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve
      );
      t = -t;
      return res;
    }

    //!
    //! \param  qx      \f$x\f$-coordinate of the point
    //! \param  qy      \f$y\f$-coordinate of the point
    //! \param  s_begin initial curvilinear coordinate of the search range
    //! \param  s_end   final curvilinear coordinate of the search range
    //! \param  x       \f$x\f$-coordinate of the projected point on the curve
    //! \param  y       \f$y\f$-coordinate of the projected point on the curve
    //! \param  s       parameter on the curve of the projection
    //! \param  t       curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst     distance point projected point
    //! \param  icurve  number of the segment with the projected point
    //!
    //! \return 1 ok -1 projection failed
    //!
    integer
    closest_point_in_s_range_ISO(
      real_type const qx,
      real_type const qy,
      real_type const s_begin,
      real_type const s_end,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst,
      integer       & icurve
    ) const;

    //!
    //! \param  qx      \f$x\f$-coordinate of the point
    //! \param  qy      \f$y\f$-coordinate of the point
    //! \param  s_begin initial curvilinear coordinate of the search range
    //! \param  s_end   final curvilinear coordinate of the search range
    //! \param  x       \f$x\f$-coordinate of the projected point on the curve
    //! \param  y       \f$y\f$-coordinate of the projected point on the curve
    //! \param  s       parameter on the curve of the projection
    //! \param  t       curvilinear coordinate of the point x,y (if orthogonal projection)
    //! \param  dst     distance point projected point
    //! \param  icurve  number of the segment with the projected point
    //!
    //! \return 1 ok -1 projection failed
    //!
    integer
    closest_point_in_s_range_SAE(
      real_type const qx,
      real_type const qy,
      real_type const s_begin,
      real_type const s_end,
      real_type     & x,
      real_type     & y,
      real_type     & s,
      real_type     & t,
      real_type     & dst,
      integer       & icurve
    ) const {
      integer res = this->closest_point_in_s_range_ISO(
        qx, qy, s_begin, s_end, x, y, s, t, dst, icurve
      );
      t = -t;
      return res;
    }

    string info() const;

    void
    info( ostream_type & stream ) const override
    { stream << this->info(); }

    friend
    ostream_type &
    operator << ( ostream_type & stream, ClothoidList const & CL );

    //!
    //! Return the clothoid list as a list of nodes and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] kappa curvature
    //!
    void
    get_SK( real_type s[], real_type kappa[] ) const;

    //!
    //! Return the clothoid list as a list of nodes and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] kappa curvature
    //!
    void
    get_SK(
      vector<real_type> & s,
      vector<real_type> & kappa
    ) const {
      s.resize( m_clothoid_list.size()+1 );
      kappa.resize( m_clothoid_list.size()+1 );
      get_SK( &s.front(), &kappa.front() );
    }

    //!
    //! Return the clothoid list as a list of nodes angles and curvatures
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
    //! Return the clothoid list as a list of nodes angles and curvatures
    //!
    //! \param[out] s     nodes
    //! \param[out] theta angles
    //! \param[out] kappa curvature
    //!
    void
    get_STK(
      vector<real_type> & s,
      vector<real_type> & theta,
      vector<real_type> & kappa
    ) const {
      s.resize( m_clothoid_list.size()+1 );
      theta.resize( m_clothoid_list.size()+1 );
      kappa.resize( m_clothoid_list.size()+1 );
      get_STK( &s.front(), &theta.front(), &kappa.front() );
    }

    //!
    //! Return the points of the clothoid list at breakpoints
    //!
    //! \param[out] x \f$x\f$-coordinates
    //! \param[out] y \f$y\f$-coordinates
    //!
    void
    get_XY( real_type x[], real_type y[] ) const;

    void
    get_delta_theta( real_type delta_theta[] ) const;

    void
    get_delta_kappa( real_type deltaKappa[] ) const;

    //!
    //! Find parametric coordinate.
    //!
    //! \param  x    \f$x\f$-coordinate point
    //! \param  y    \f$y\f$-coordinate point
    //! \param  s    value \f$ s \f$
    //! \param  t    value \f$ t \f$
    //! \return idx  the segment with point at minimal distance, otherwise
    //!              -(idx+1) if \f$ (x,y) \f$ cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      real_type const x,
      real_type const y,
      real_type     & s,
      real_type     & t
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
    //!                `-(idx+1)` if \f$ (x,y) \f$ cannot be projected orthogonally on the segment
    //!
    integer
    findST1(
      integer   const ibegin,
      integer   const iend,
      real_type const x,
      real_type const y,
      real_type     & s,
      real_type     & t
    ) const;

    /*\
     |             _ _ _     _
     |    ___ ___ | | (_)___(_) ___  _ __
     |   / __/ _ \| | | / __| |/ _ \| '_ \
     |  | (_| (_) | | | \__ \ | (_) | | | |
     |   \___\___/|_|_|_|___/_|\___/|_| |_|
    \*/

    //!
    //! Detect a collision with another clothoid list
    //!
    bool
    collision( ClothoidList const & CL ) const {
      return collision_ISO( 0, CL, 0 );
    }

    //!
    //! Detect a collision with another clothoid list with offset
    //!
    //! \param[in] offs   offset of first clothoid list
    //! \param[in] CL     second clothoid list
    //! \param[in] offs_C offset of second clothoid list
    //!
    bool
    collision_ISO(
      real_type    const   offs,
      ClothoidList const & CL,
      real_type    const   offs_C
    ) const;

    bool
    collision( BaseCurve const * pC ) const override;

    bool
    collision_ISO(
      real_type const   offs,
      BaseCurve const * pC,
      real_type const   offs_C
    ) const override;

    /*\
     |   _       _                          _
     |  (_)_ __ | |_ ___ _ __ ___  ___  ___| |_
     |  | | '_ \| __/ _ \ '__/ __|/ _ \/ __| __|
     |  | | | | | ||  __/ |  \__ \  __/ (__| |_
     |  |_|_| |_|\__\___|_|  |___/\___|\___|\__|
    \*/

    //!
    //! Intersect a clothoid list with another clothoid list
    //!
    //! \param[in]  CL    second clothoid list
    //! \param[out] ilist list of the intersection (as parameter on the curves)
    //!
    void
    intersect(
      ClothoidList const & CL,
      IntersectList      & ilist
    ) const {
      this->intersect_ISO( 0, CL, 0, ilist );
    }

    //!
    //! Intersect a clothoid list with another clothoid list with offset (ISO)
    //!
    //! \param[in]  offs     offset of first clothoid list
    //! \param[in]  CL       second clothoid list
    //! \param[in]  offs_obj offset of second clothoid list
    //! \param[out] ilist    list of the intersection (as parameter on the curves)
    //!
    void
    intersect_ISO(
      real_type    const   offs,
      ClothoidList const & CL,
      real_type    const   offs_obj,
      IntersectList      & ilist
    ) const;

    void
    intersect(
      BaseCurve const * pC,
      IntersectList   & ilist
    ) const override;

    void
    intersect_ISO(
      real_type const   offs,
      BaseCurve const * pC,
      real_type const   offs_LS,
      IntersectList   & ilist
    ) const override;

    //!
    //! Save Clothoid list to a stream
    //!
    //!\param stream stream to save
    //!
    void
    export_table( ostream_type & stream ) const;

    //!
    //! Save Clothoid list to a stream
    //!
    //!\param stream streamstream to save
    //!
    void
    export_ruby( ostream_type & stream ) const;

    //!
    //! Save the clothoid list on a stream.
    //! The data is saved as follows
    //!
    //!       # x y theta kappa
    //!       x0 y0 theta0 kappa0
    //!       x1 y1 theta1 kappa1
    //!       ...
    //!       xn yn thetan kappan
    //!
    void save( ostream_type & stream ) const;

    //!
    //! Read the clothoid list from a stream.
    //! The data is assumed to be saved as follows
    //!
    //!       # x y theta kappa
    //!       x0 y0 theta0 kappa0
    //!       x1 y1 theta1 kappa1
    //!       ...
    //!       xn yn thetan kappan
    //!
    void load( istream_type & stream, real_type const epsi = 1e-8 );

#ifdef CLOTHOIDS_BACK_COMPATIBILITY
#include "ClothoidList_compatibility.hxx"
#endif

  };

}

///
/// eof: ClothoidList.hxx
///
