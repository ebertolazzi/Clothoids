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
/// file: Clothoid.hh
///

#ifndef CLOTHOID_HH
#define CLOTHOID_HH

#include "Fresnel.hh"
#include "Triangle2D.hh"
#include <vector>

//! Clothoid computations routine
namespace G2lib {

  using std::vector ;

  typedef Triangle2D<valueType> T2D ;

  /*\
   |    ____ _       _   _           _     _
   |   / ___| | ___ | |_| |__   ___ (_) __| |
   |  | |   | |/ _ \| __| '_ \ / _ \| |/ _` |
   |  | |___| | (_) | |_| | | | (_) | | (_| |
   |   \____|_|\___/ \__|_| |_|\___/|_|\__,_|
  \*/
  /*! \brief Compute the clothoid by Hemite data
   *
   * \param x0     initial x position            \f$ x_0      \f$
   * \param y0     initial y position            \f$ y_0      \f$
   * \param theta0 initial angle                 \f$ \theta_0 \f$
   * \param x1     final x position              \f$ x_1      \f$
   * \param y1     final y position              \f$ y_1      \f$
   * \param theta1 final angle                   \f$ \theta_1 \f$
   * \param k      computed curvature            \f$ K        \f$
   * \param dk     computed curvature derivative \f$ K'       \f$
   * \param L      computed length of the curve
   */
  indexType
  buildClothoid( valueType   x0,
                 valueType   y0,
                 valueType   theta0,
                 valueType   x1,
                 valueType   y1,
                 valueType   theta1,
                 valueType & k,
                 valueType & dk,
                 valueType & L ) ;

  indexType
  buildClothoid( valueType   x0,
                 valueType   y0,
                 valueType   theta0,
                 valueType   x1,
                 valueType   y1,
                 valueType   theta1,
                 valueType & k,
                 valueType & dk,
                 valueType & L,
                 valueType & k_1,
                 valueType & dk_1,
                 valueType & L_1,
                 valueType & k_2,
                 valueType & dk_2,
                 valueType & L_2 ) ;

  //! Compute Lommel function
  valueType
  LommelReduced( valueType mu, valueType nu, valueType z ) ;

  /*\
   |    ____ _       _   _           _     _  ____
   |   / ___| | ___ | |_| |__   ___ (_) __| |/ ___|   _ _ ____   _____
   |  | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |  | | | | '__\ \ / / _ \
   |  | |___| | (_) | |_| | | | (_) | | (_| | |__| |_| | |   \ V /  __/
   |   \____|_|\___/ \__|_| |_|\___/|_|\__,_|\____\__,_|_|    \_/ \___|
  \*/
  //! \brief Class to manage Clothoid Curve
  class ClothoidCurve {
  public:

    typedef struct {
      valueType    s0 ;
      valueType    L  ;
      ClothoidData cd ;
      T2D          t  ;
    } bbData ;

    typedef struct {
      valueType    split_angle ;
      valueType    split_size ;
      valueType    split_offs ;
      valueType    s0 ;
      valueType    L ;
      ClothoidData cd ;
    } bbData2 ;

  private:
    ClothoidData CD ; //!< clothoid data
    valueType    L ;  //!< lenght of clothoid segment

    int
    build( valueType x1, valueType y1, valueType theta1 );

    void
    bbSplit_internal( bbData2 const & data, vector<bbData> & bbV ) const ;

    //! Use newton and bisection to intersect two small clothoid segment
    bool
    intersect_internal( bbData const & c1, valueType c1_offs, valueType & s1,
                        bbData const & c2, valueType c2_offs, valueType & s2,
                        indexType max_iter,
                        valueType tolerance ) const ;

  public:

    ClothoidCurve()
    {
      CD.x0     = 0 ;
      CD.y0     = 0 ;
      CD.theta0 = 0 ;
      CD.kappa0 = 0 ;
      CD.dk     = 0 ;
      L         = 0 ;
    }

    //! construct a clothoid with the standard parameters
    ClothoidCurve( valueType _x0,
                   valueType _y0,
                   valueType _theta0,
                   valueType _k,
                   valueType _dk,
                   valueType _L )
    {
      CD.x0     = _x0 ;
      CD.y0     = _y0 ;
      CD.theta0 = _theta0 ;
      CD.kappa0 = _k ;
      CD.dk     = _dk ;
      L         = _L ;
    }

    //! construct a clothoid by solving the hermite G1 problem
    ClothoidCurve( valueType const _P0[],
                   valueType       _theta0,
                   valueType const _P1[],
                   valueType       _theta1 )
    {
      CD.x0     = _P0[0] ;
      CD.y0     = _P0[1] ;
      CD.theta0 = _theta0 ;
      build( _P1[0], _P1[1], _theta1 ) ;
    }

    void
    copy( ClothoidCurve const & c ) {
      CD = c.CD ;
      L  = c.L ;
    }

    ClothoidCurve( ClothoidCurve const & s ) { copy(s) ; }

    ClothoidCurve const & operator = ( ClothoidCurve const & s )
    { copy(s) ; return *this ; }

    valueType getX0()      const { return CD.x0 ; }
    valueType getY0()      const { return CD.y0 ; }
    valueType getTheta0()  const { return CD.theta0 ; }

    valueType getKappa()   const { return CD.kappa0 ; }
    valueType getKappa_D() const { return CD.dk ; }
    valueType getL()       const { return L ; }

    valueType getThetaBegin() const { return CD.theta0 ; }
    valueType getThetaEnd()   const { return CD.theta(L) ; }
    valueType getKappaBegin() const { return CD.kappa0 ; }
    valueType getKappaEnd()   const { return CD.kappa(L) ; }

    //! construct a clothoid with the standard parameters
    void
    build( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _k,
           valueType _dk,
           valueType _L ) {
      CD.x0     = _x0 ;
      CD.y0     = _y0 ;
      CD.theta0 = _theta0 ;
      CD.kappa0 = _k ;
      CD.dk     = _dk ;
      L         = _L ;
    }

    /*! \brief build a clothoid by solving the hermite G1 problem
     *
     * \param x0     initial x position            \f$ x_0      \f$
     * \param y0     initial y position            \f$ y_0      \f$
     * \param theta0 initial angle                 \f$ \theta_0 \f$
     * \param x1     final x position              \f$ x_1      \f$
     * \param y1     final y position              \f$ y_1      \f$
     * \param theta1 final angle                   \f$ \theta_1 \f$
     * \return number of iteration performed
     */
    int
    build_G1( valueType x0,
              valueType y0,
              valueType theta0,
              valueType x1,
              valueType y1,
              valueType theta1,
              valueType tol = 1e-12 ) ;

    /*! \brief build a clothoid by solving the hermite G1 problem
     *
     * \param x0     initial x position            \f$ x_0      \f$
     * \param y0     initial y position            \f$ y_0      \f$
     * \param theta0 initial angle                 \f$ \theta_0 \f$
     * \param x1     final x position              \f$ x_1      \f$
     * \param y1     final y position              \f$ y_1      \f$
     * \param theta1 final angle                   \f$ \theta_1 \f$
     * \return number of iteration performed
     */
    int
    build_G1_D( valueType x0,
                valueType y0,
                valueType theta0,
                valueType x1,
                valueType y1,
                valueType theta1,
                valueType L_D[2],
                valueType k_D[2],
                valueType dk_D[2],
                valueType tol = 1e-12 ) ;

    /*! \brief build a clothoid by solving the forward problem
     *
     * \param x0     initial x position            \f$ x_0      \f$
     * \param y0     initial y position            \f$ y_0      \f$
     * \param theta0 initial angle                 \f$ \theta_0 \f$
     * \param k      initial curvature             \f$ \kappa_0 \f$
     * \param x1     final x position              \f$ x_1      \f$
     * \param y1     final y position              \f$ y_1      \f$
     */
    bool
    build_forward( valueType x0,
                   valueType y0,
                   valueType theta0,
                   valueType k,
                   valueType x1,
                   valueType y1,
                   valueType tol = 1e-8 ) ;

    /*! \brief get clothoid angle at curvilinear cooordinate `s`
     *
     * \param  s curvilinear cooordinate
     * \return angle (radiant) at curvilinear cooordinate `s`
     */
    valueType
    theta( valueType s ) const
    { return CD.theta(s) ; }

    /*! \brief get clothoid angle derivative (=curvature) at curvilinear cooordinate `s`
     *
     * \param  s curvilinear cooordinate
     * \return angle derivative (radiant/s) at curvilinear cooordinate `s`
     */
    valueType
    theta_D( valueType s ) const
    { return CD.kappa(s) ; }

    /*! \brief get clothoid angle second derivative at curvilinear cooordinate `s`
     *
     * \return angle second derivative (radiant/s^2) at curvilinear cooordinate `s`
     */
    valueType
    theta_DD( valueType ) const
    { return CD.dk ; }

    /*! \brief get clothoid angle third derivative at curvilinear cooordinate `s`
     *
     * \return angle third derivative (radiant/s^3) at curvilinear cooordinate `s`
     */
    valueType
    theta_DDD( valueType ) const { return 0 ; }

    /*! \return clothoid total variation
     */
    valueType
    thetaTotalVariation() const ;

    valueType
    thetaMinMax( valueType & thMin, valueType & thMax ) const ;

    /*! \return clothoid angle range
     */
    valueType
    deltaTheta() const
    { valueType thMin, thMax ; return thetaMinMax( thMin, thMax ) ; }

    valueType
    curvatureMinMax( valueType & kMin, valueType & kMax ) const ;

    /*! \return clothoid total curvature variation
     */
    valueType
    curvatureTotalVariation() const ;

    valueType
    integralCurvature2() const ;

    valueType
    integralJerk2() const ;

    valueType
    integralSnap2() const ;

    /*! \brief clothoid X coordinate at curvilinear coordinate `s`
     * \param s curvilinear coordinate
     * \return clothoid X coordinate
     */
    valueType X( valueType s ) const { return CD.X(s) ; }

    /*! \brief clothoid Y coordinate at curvilinear coordinate `s`
     * \param s curvilinear coordinate
     * \return clothoid Y coordinate
     */
    valueType Y( valueType s ) const { return CD.Y(s) ; }

    valueType Xbegin()     const { return CD.x0 ; }
    valueType Ybegin()     const { return CD.y0 ; }
    valueType KappaBegin() const { return CD.kappa0 ; }
    valueType ThetaBegin() const { return CD.theta0 ; }

    valueType Xend()     const { return X(L) ; }
    valueType Yend()     const { return Y(L) ; }
    valueType KappaEnd() const { return theta_D(L) ; }
    valueType ThetaEnd() const { return theta(L) ; }

    void
    Pinfinity( valueType & x, valueType & y, bool plus = true ) const
    { CD.Pinfinity( x, y, plus ); }

    void
    eval( valueType   s,
          valueType & theta,
          valueType & kappa,
          valueType & x,
          valueType & y ) const
    { CD.eval( s, theta, kappa, x, y ); }

    void
    eval( valueType   s,
          valueType & x,
          valueType & y ) const {
      CD.eval( s, x, y ) ;
    }

    void
    eval_D( valueType   s,
            valueType & x_D,
            valueType & y_D ) const {
      CD.eval_D( s, x_D, y_D ) ;
    }

    void
    eval_DD( valueType   s,
             valueType & x_DD,
             valueType & y_DD ) const {
      CD.eval_DD( s, x_DD, y_DD ) ;
    }

    void
    eval_DDD( valueType   s,
              valueType & x_DDD,
              valueType & y_DDD ) const {
      CD.eval_DDD( s, x_DDD, y_DDD ) ;
    }

    // offset curve
    void
    eval( valueType   s,
          valueType   offs,
          valueType & x,
          valueType & y ) const {
      CD.eval( s, offs, x, y ) ;
    }

    void
    eval_D( valueType   s,
            valueType   offs,
            valueType & x_D,
            valueType & y_D ) const {
      CD.eval_D( s, offs, x_D, y_D ) ;
    }

    void
    eval_DD( valueType   s,
             valueType   offs,
             valueType & x_DD,
             valueType & y_DD ) const {
      CD.eval_DD( s, offs, x_DD, y_DD ) ;
    }

    void
    eval_DDD( valueType   s,
              valueType   offs,
              valueType & x_DDD,
              valueType & y_DDD ) const {
      CD.eval_DDD( s, offs, x_DDD, y_DDD ) ;
    }

    valueType
    closestPoint( valueType   qx,
                  valueType   qy,
                  valueType & X,
                  valueType & Y,
                  valueType & S ) const ;

    valueType
    distance( valueType qx, valueType qy, valueType & S ) const {
      valueType X, Y;
      return closestPoint( qx, qy, X, Y, S );
    }

    valueType
    distance( valueType qx, valueType qy ) const {
      valueType X, Y, S ;
      return closestPoint( qx, qy, X, Y, S );
    }

    valueType
    closestPointBySample( valueType   ds,
                          valueType   qx,
                          valueType   qy,
                          valueType & X,
                          valueType & Y,
                          valueType & S ) const ;

    valueType
    distanceBySample( valueType   ds,
                      valueType   qx,
                      valueType   qy,
                      valueType & S ) const {
      valueType X, Y;
      return closestPointBySample( ds, qx, qy, X, Y, S );
    }

    valueType
    distanceBySample( valueType ds,
                      valueType qx,
                      valueType qy ) const {
      valueType X, Y, S ;
      return closestPointBySample( ds, qx, qy, X, Y, S );
    }

    void
    trim( valueType s_begin, valueType s_end ) {
      valueType xx, yy ;
      CD.eval( s_begin, xx, yy ) ;
      CD.kappa0 += s_begin * CD.dk ;
      CD.theta0 += s_begin * ( CD.kappa0 + 0.5*s_begin * CD.dk ) ;
      L          = s_end - s_begin ;
      CD.x0 = xx ;
      CD.y0 = yy ;
    }

    //! get the bounding box triangle (if angle variation less that pi/2)
    bool
    bbTriangle( valueType offs,
                valueType p0[2],
                valueType p1[2],
                valueType p2[2] ) const {
      return CD.bbTriangle( L, offs, p0, p1, p2 ) ;
    }

    bool
    bbTriangle( valueType offs, T2D & t ) const {
      valueType p0[2], p1[2], p2[2] ;
      bool ok = CD.bbTriangle( L, offs, p0, p1, p2 ) ;
      if ( ok ) t.setup( p0, p1, p2 ) ;
      return ok ;
    }

    /*! \brief split clothois in smaller segments
     *
     * \param split_angle maximum angle variation
     * \param split_size  maximum height of the triangle
     * \param split_offs  curve offset
     * \param bb          splitting data structures vector
     */
    void
    bbSplit( valueType        split_angle,
             valueType        split_size,
             valueType        split_offs,
             vector<bbData> & bb ) const ;

    // intersect computation
    void
    intersect( ClothoidCurve const & c,
               vector<valueType>   & s1,
               vector<valueType>   & s2,
               indexType             max_iter,
               valueType             tolerance ) const {
      intersect( 0, c, 0, s1, s2, max_iter, tolerance ) ;
    }

    void
    intersect( valueType             offs,
               ClothoidCurve const & c,
               valueType             c_offs,
               vector<valueType>   & s1,
               vector<valueType>   & s2,
               indexType             max_iter,
               valueType             tolerance ) const ;

    // collision detection
    bool
    approximate_collision( valueType             offs,
                           ClothoidCurve const & c,
                           valueType             c_offs,
                           valueType             max_angle,         //!< maximum angle variation
                           valueType             max_size ) const ; //!< curve offset

    void
    rotate( valueType angle, valueType cx, valueType cy ) ;

    void
    translate( valueType tx, valueType ty )
    { CD.x0 += tx ; CD.y0 += ty ; }

    void
    changeCurvilinearOrigin( valueType s0, valueType newL ) ;

    void
    changeOrigin( valueType newx0, valueType newy0 )
    { CD.x0 = newx0 ; CD.y0 = newy0 ; }

    void
    scale( valueType sfactor ) ;

    void
    reverse() ;

    friend
    std::ostream &
    operator << ( std::ostream & stream, ClothoidCurve const & c ) ;

  } ;

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/
  // Clothoid-clothoid
  class G2solve2arc {

    valueType tolerance ;
    int       maxIter ;

    valueType x0 ;
    valueType y0 ;
    valueType theta0 ;
    valueType kappa0 ;

    valueType x1 ;
    valueType y1 ;
    valueType theta1 ;
    valueType kappa1 ;

    // standard problem
    valueType lambda, phi, xbar, ybar ;
    valueType th0, th1 ;
    valueType k0, k1 ;
    valueType DeltaK ;
    valueType DeltaTheta ;

    ClothoidCurve S0, S1 ;

    void
    evalA( valueType   alpha,
           valueType   L,
           valueType & A ) const ;

    void
    evalA( valueType   alpha,
           valueType   L,
           valueType & A,
           valueType & A_1,
           valueType & A_2 ) const ;

    void
    evalG( valueType alpha,
           valueType L,
           valueType th,
           valueType k,
           valueType G[2] ) const ;

    void
    evalG( valueType alpha,
           valueType L,
           valueType th,
           valueType k,
           valueType G[2],
           valueType G_1[2],
           valueType G_2[2] ) const ;

    void
    evalF( valueType const vars[2], valueType F[2] ) const ;

    void
    evalFJ( valueType const vars[2],
            valueType       F[2],
            valueType       J[2][2] ) const ;

    void
    buildSolution( valueType alpha, valueType L ) ;

  public:

    G2solve2arc()
    : tolerance(1e-10)
    , maxIter(20)
    , x0(0)
    , y0(0)
    , theta0(0)
    , kappa0(0)
    , x1(0)
    , y1(0)
    , theta1(0)
    , kappa1(0)
    , lambda(0)
    , phi(0)
    , xbar(0)
    , ybar(0)
    , th0(0)
    , th1(0)
    , k0(0)
    , k1(0)
    {}

    ~G2solve2arc() {}

    int
    build( valueType x0, valueType y0, valueType theta0, valueType kappa0,
           valueType x1, valueType y1, valueType theta1, valueType kappa1 ) ;

    void
    setTolerance( valueType tol ) ;

    void
    setMaxIter( int tol ) ;

    int
    solve() ;

    ClothoidCurve const & getS0() const { return S0 ; }
    ClothoidCurve const & getS1() const { return S1 ; }

  } ;

  /*\
   |    ____ ____            _            ____ _     ____
   |   / ___|___ \ ___  ___ | |_   _____ / ___| |   / ___|
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |   | |  | |
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/ |___| |__| |___
   |   \____|_____|___/\___/|_| \_/ \___|\____|_____\____|
  \*/

  // Clothoid-line-clothoid
  class G2solveCLC {

    valueType tolerance ;
    int       maxIter ;

    valueType x0 ;
    valueType y0 ;
    valueType theta0 ;
    valueType kappa0 ;
    valueType x1 ;
    valueType y1 ;
    valueType theta1 ;
    valueType kappa1 ;

    // standard problem
    valueType lambda, phi, xbar, ybar ;
    valueType th0, th1 ;
    valueType k0, k1 ;

    ClothoidCurve S0, SM, S1 ;

    bool
    buildSolution( valueType sM, valueType thM ) ;

  public:

    G2solveCLC()
    : tolerance(1e-10)
    , maxIter(20)
    , x0(0)
    , y0(0)
    , theta0(0)
    , kappa0(0)
    , x1(0)
    , y1(0)
    , theta1(0)
    , kappa1(0)
    , lambda(0)
    , phi(0)
    , xbar(0)
    , ybar(0)
    , th0(0)
    , th1(0)
    , k0(0)
    , k1(0)
    {}

    ~G2solveCLC() {}

    int
    build( valueType x0, valueType y0, valueType theta0, valueType kappa0,
           valueType x1, valueType y1, valueType theta1, valueType kappa1 ) ;

    void
    setTolerance( valueType tol ) ;

    void
    setMaxIter( int tol ) ;

    int
    solve() ;

    ClothoidCurve const & getS0() const { return S0 ; }
    ClothoidCurve const & getSM() const { return SM ; }
    ClothoidCurve const & getS1() const { return S1 ; }

  } ;

  /*\
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/
  // Clothoid-clothoid-clothoid with G2 continuity
  class G2solve3arc {

    ClothoidCurve S0, SM, S1 ;

    valueType tolerance ;
    int       maxIter ;

    // G2 interpolation data
    valueType x0 ;
    valueType y0 ;
    valueType theta0 ;
    valueType kappa0 ;
    valueType x1 ;
    valueType y1 ;
    valueType theta1 ;
    valueType kappa1 ;

    // standard scaled problem
    valueType phi, Lscale ;
    valueType th0, th1 ;
    valueType s0, s1 ;

    // precomputed values
    valueType K0, K1, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14 ;

    void
    evalFJ( valueType const vars[2],
            valueType       F[2],
            valueType       J[2][2] ) const ;

    void
    evalF( valueType const vars[2], valueType F[2] ) const ;

    void
    buildSolution( valueType sM, valueType thM ) ;

    int
    solve( valueType sM_guess, valueType thM_guess ) ;

  public:

    G2solve3arc()
    : tolerance(1e-10)
    , maxIter(100)
    {}

    ~G2solve3arc() {}

    void setTolerance( valueType tol ) ;
    void setMaxIter( int miter ) ;

    /*!
     | Compute the 3 arc clothoid spline that fit the data
     |
     | \param[in] x0      initial `x` position
     | \param[in] y0      initial `y` position
     | \param[in] theta0  initial angle
     | \param[in] kappa0  initial curvature
     | \param[in] x1      final `x` position
     | \param[in] y1      final `y` position
     | \param[in] theta1  final angle
     | \param[in] kappa1  final curvature
     | \param[in] Dmax    rough desidered maximum angle variation, if 0 computed automatically
     | \param[in] dmax    rough desidered maximum angle divergence from guess, if 0 computed automatically
     | \return
     |
    \*/
    int
    build( valueType x0,
           valueType y0,
           valueType theta0,
           valueType kappa0,
           valueType x1,
           valueType y1,
           valueType theta1,
           valueType kappa1,
           valueType Dmax = 0,
           valueType dmax = 0 ) ;


    /*!
     | \return get the first clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getS0() const { return S0 ; }
    /*!
     | \return get the last clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getS1() const { return S1 ; }
    /*!
     | \return get the middle clothoid for the 3 arc G2 fitting
    \*/
    ClothoidCurve const & getSM() const { return SM ; }

    /*!
     | \return get the length of the 3 arc G2 fitting
    \*/
    valueType
    totalLength() const {
      return S0.getL() + S1.getL() + SM.getL() ;
    }

    /*!
     | \return get the total angle variation of the 3 arc G2 fitting
    \*/
    valueType
    thetaTotalVariation() const {
      return S0.thetaTotalVariation() +
             S1.thetaTotalVariation() +
             SM.thetaTotalVariation() ;
    }

    /*!
     | \return get the total curvature variation of the 3 arc G2 fitting
    \*/
    valueType
    curvatureTotalVariation() const {
      return S0.curvatureTotalVariation() +
             S1.curvatureTotalVariation() +
             SM.curvatureTotalVariation() ;
    }

    /*!
     | \return get the integral of the curvature squared of the 3 arc G2 fitting
    \*/
    valueType
    integralCurvature2() const {
      return S0.integralCurvature2() +
             S1.integralCurvature2() +
             SM.integralCurvature2() ;
    }

    /*!
     | \return get the integral of the jerk squared of the 3 arc G2 fitting
    \*/
    valueType
    integralJerk2() const {
      return S0.integralJerk2() +
             S1.integralJerk2() +
             SM.integralJerk2() ;
    }

    /*!
     | \return get the integral of the snap squared of the 3 arc G2 fitting
    \*/
    valueType
    integralSnap2() const {
      return S0.integralSnap2() +
             S1.integralSnap2() +
             SM.integralSnap2() ;
    }

    /*!
     | \param[out] thMin minimum angle in the 3 arc G2 fitting curve
     | \param[out] thMax maximum angle in the 3 arc G2 fitting curve
     | \return the difference of `thMax` and `thMin`
    \*/
    valueType
    thetaMinMax( valueType & thMin, valueType & thMax ) const ;

    /*!
     | \return the difference of maximum-minimum angle in the 3 arc G2 fitting curve
    \*/
    valueType
    deltaTheta() const
    { valueType thMin, thMax ; return thetaMinMax( thMin, thMax ) ; }

    /*!
     | \param[out] kMin minimum curvature in the 3 arc G2 fitting curve
     | \param[out] kMax maximum curvature in the 3 arc G2 fitting curve
     | \return the difference of `kMax` and `kMin`
    \*/
    valueType
    curvatureMinMax( valueType & kMin, valueType & kMax ) const ;

    /*!
     | \return angle as a function of curvilinear coordinate
    \*/
    valueType theta( valueType s ) const ;

    /*!
     | \return angle derivative (curvature) as a function of curvilinear coordinate
    \*/
    valueType theta_D( valueType s ) const ;

    /*!
     | \return angle second derivative (curvature derivative) as a function of curvilinear coordinate
    \*/
    valueType theta_DD( valueType s ) const ;

    /*!
     | \return angle third derivative as a function of curvilinear coordinate
    \*/
    valueType theta_DDD( valueType s ) const ;

    /*!
     | \return x coordinate of the3 arc clothoid as a function of curvilinear coordinate
    \*/
    valueType X( valueType s ) const ;

    /*!
     | \return y coordinate of the3 arc clothoid as a function of curvilinear coordinate
    \*/
    valueType Y( valueType s ) const ;

    /*!
     | \return initial x coordinate of the 3 arc clothoid
    \*/
    valueType Xbegin()     const { return S0.Xbegin() ; }

    /*!
     | \return initial y coordinate of the 3 arc clothoid
    \*/
    valueType Ybegin()     const { return S0.Ybegin() ; }

    /*!
     | \return initial curvature of the 3 arc clothoid
    \*/
    valueType KappaBegin() const { return S0.KappaBegin() ; }

    /*!
     | \return initial angle of the 3 arc clothoid
    \*/
    valueType ThetaBegin() const { return S0.ThetaBegin() ; }

    /*!
     | \return final x coordinate of the 3 arc clothoid
    \*/
    valueType Xend()     const { return S1.Xend() ; }

    /*!
     | \return final y coordinate of the 3 arc clothoid
    \*/
    valueType Yend()     const { return S1.Yend() ; }

    /*!
     | \return final curvature of the 3 arc clothoid
    \*/
    valueType KappaEnd() const { return S1.KappaEnd() ; }

    /*!
     | \return final angle of the 3 arc clothoid
    \*/
    valueType ThetaEnd() const { return S1.ThetaEnd() ; }

    /*!
     | Compute parameters of 3 arc clothoid at curvilinear coordinate `s`
     |
     | \param[in]  s     curvilinear coordinate of where curve is computed
     | \param[out] theta the curve angle
     | \param[out] kappa the curve curvature
     | \param[out] x     the curve x-coordinate
     | \param[out] y     the curve y-coordinate
    \*/
    void
    eval( valueType   s,
          valueType & theta,
          valueType & kappa,
          valueType & x,
          valueType & y ) const ;

    void eval( valueType s, valueType & x, valueType & y ) const ;
    void eval_D( valueType s, valueType & x_D, valueType & y_D ) const ;
    void eval_DD( valueType s, valueType & x_DD, valueType & y_DD ) const ;
    void eval_DDD( valueType s, valueType & x_DDD, valueType & y_DDD ) const ;

    // offset curve
    void eval( valueType s, valueType offs, valueType & x, valueType & y ) const ;
    void eval_D( valueType s, valueType offs, valueType & x_D, valueType & y_D ) const ;
    void eval_DD( valueType s, valueType offs, valueType & x_DD, valueType & y_DD ) const ;
    void eval_DDD( valueType s, valueType offs, valueType & x_DDD, valueType & y_DDD ) const ;

    void
    rotate( valueType angle, valueType cx, valueType cy ) {
      S0.rotate( angle, cx, cy ) ;
      S1.rotate( angle, cx, cy ) ;
      SM.rotate( angle, cx, cy ) ;
    }

    void
    translate( valueType tx, valueType ty ){
      S0.translate( tx, ty ) ;
      S1.translate( tx, ty ) ;
      SM.translate( tx, ty ) ;
    }

    void
    reverse() {
      std::swap( S0, S1 ) ;
      S0.reverse() ;
      S1.reverse() ;
      SM.reverse() ;
    }

    friend
    std::ostream &
    operator << ( std::ostream & stream, ClothoidCurve const & c ) ;

  } ;

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/
  //! \brief Class to manage a list Clothoid Curve (not necessarily G2 or G1 connected)
  class ClothoidList {

    std::vector<valueType>     s0 ;
    std::vector<ClothoidCurve> clotoidList ;

  public:

    ClothoidList() {}
    ~ClothoidList() ;

    ClothoidList( ClothoidList const & s ) { copy(s) ; }
    ClothoidList const & operator = ( ClothoidList const & s )
    { copy(s) ; return *this ; }

    void reserve( indexType n );

    void copy( ClothoidList const & L );

    void add( ClothoidCurve const & c );

    ClothoidCurve const & get( indexType idx ) const;
    ClothoidCurve const & getAtS( valueType s, indexType & last_idx ) const;

    bool findAtS( valueType s, indexType & last_idx ) const;

    valueType theta( valueType s, indexType & last_idx ) const;
    valueType theta_D( valueType s, indexType & last_idx ) const;
    valueType theta_DD( valueType s, indexType & last_idx ) const;
    valueType theta_DDD( valueType, indexType & ) const { return 0 ; }

    valueType totalLength() const {
      if ( s0.empty() ) return 0;
      return s0.back() - s0.front() ;
    }

    valueType X( valueType s, indexType & last_idx ) const ;
    valueType Y( valueType s, indexType & last_idx ) const ;

    void
    eval( valueType   s,
          indexType & last_idx,
          valueType & theta,
          valueType & kappa,
          valueType & x,
          valueType & y ) const ;

    void
    eval( valueType   s,
          indexType & last_idx,
          valueType & x,
          valueType & y ) const ;
    void
    eval_D( valueType   s,
            indexType & last_idx,
            valueType & x_D,
            valueType & y_D ) const ;
    void
    eval_DD( valueType   s,
             indexType & last_idx,
             valueType & x_DD,
             valueType & y_DD ) const ;
    void
    eval_DDD( valueType   s,
              indexType & last_idx,
              valueType & x_DDD,
              valueType & y_DDD ) const ;

    // offset curve
    void
    eval( valueType   s,
          indexType & last_idx,
          valueType   offs,
          valueType & x,
          valueType & y ) const ;
    void
    eval_D( valueType   s,
            indexType & last_idx,
            valueType   offs,
            valueType & x_D,
            valueType & y_D ) const ;
    void
    eval_DD( valueType   s,
             indexType & last_idx,
             valueType   offs,
             valueType & x_DD,
             valueType & y_DD ) const ;
    void
    eval_DDD( valueType   s,
              indexType & last_idx,
              valueType   offs,
              valueType & x_DDD,
              valueType & y_DDD ) const ;

    valueType
    distance( valueType   qx,
              valueType   qy,
              valueType & S ) const ;

    valueType
    closestPoint( valueType   qx,
                  valueType   qy,
                  valueType & X,
                  valueType & Y,
                  valueType & S ) const ;

    void rotate( valueType angle, valueType cx, valueType cy ) ;
    void translate( valueType tx, valueType ty ) ;
    void changeOrigin( valueType newx0, valueType newy0 ) ;
    void scale( valueType sfactor ) ;
    void reverse() ;

  };
}

#endif

///
/// eof: Clothoid.hh
///
