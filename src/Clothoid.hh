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
namespace Clothoid {

  using std::vector ;
  using namespace Fresnel ;
  using namespace Triangle2D ;

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
  
  class ClothoidCurve ; // forward declaration

  /*\
   |    ____ _       _   _           _     _  ____
   |   / ___| | ___ | |_| |__   ___ (_) __| |/ ___|   _ _ ____   _____
   |  | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |  | | | | '__\ \ / / _ \
   |  | |___| | (_) | |_| | | | (_) | | (_| | |__| |_| | |   \ V /  __/
   |   \____|_|\___/ \__|_| |_|\___/|_|\__,_|\____\__,_|_|    \_/ \___|
  \*/
  //! \brief Class to manage Clothoid Curve
  class ClothoidCurve {

    valueType x0,       //!< initial x coordinate of the clothoid
              y0,       //!< initial y coordinate of the clothoid
              theta0,   //!< initial angle of the clothoid
              k,        //!< initial curvature
              dk ;      //!< curvature derivative

    valueType s_min,    //!< initial curvilinear coordinate of the clothoid segment
              s_max ;   //!< final curvilinear coordinate of the clothoid segment

    int
    build( valueType x1, valueType y1, valueType theta1 );

    void
    bbSplit_internal( valueType                        split_angle,
                      valueType                        split_size,
                      valueType                        split_offs,
                      vector<ClothoidCurve>          & c,
                      vector<Triangle2D<valueType> > & t ) const ;

    //! Use newton and bisection to intersect two small clothoid segment
    bool
    intersect_internal( ClothoidCurve & c1, valueType c1_offs, valueType & s1,
                        ClothoidCurve & c2, valueType c2_offs, valueType & s2,
                        indexType max_iter,
                        valueType tolerance ) const ;

  public:

    ClothoidCurve()
    : x0(0)
    , y0(0)
    , theta0(0)
    , k(0)
    , dk(0)
    , s_min(0)
    , s_max(0)
    {}

    //! construct a clothoid with the standard parameters
    ClothoidCurve( valueType _x0,
                   valueType _y0,
                   valueType _theta0,
                   valueType _k,
                   valueType _dk,
                   valueType _L )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , k(_k)
    , dk(_dk)
    , s_min(0)
    , s_max(_L)
    {}

    ClothoidCurve( valueType _x0,
                   valueType _y0,
                   valueType _theta0,
                   valueType _k,
                   valueType _dk,
                   valueType _smin ,
                   valueType _smax )
    : x0(_x0)
    , y0(_y0)
    , theta0(_theta0)
    , k(_k)
    , dk(_dk)
    , s_min(_smin)
    , s_max(_smax)
    {}

    //! construct a clothoid by solving the hermite G1 problem
    ClothoidCurve( valueType const _P0[],
                   valueType       _theta0,
                   valueType const _P1[],
                   valueType       _theta1 )
    : x0(_P0[0])
    , y0(_P0[1])
    , theta0(_theta0)
    , s_min(0) {
      build( _P1[0], _P1[1], _theta1 ) ;
      //buildClothoid( x0, y0, theta0, _P1[0], _P1[1], _theta1, k, dk, s_max ) ;
    }

    void
    copy( ClothoidCurve const & c ) {
      x0     = c.x0 ;
      y0     = c.y0 ;
      theta0 = c.theta0 ;
      k      = c.k ;
      dk     = c.dk ;
      s_min  = c.s_min ;
      s_max  = c.s_max ;
    }

    ClothoidCurve( ClothoidCurve const & s ) { copy(s) ; }

    ClothoidCurve const & operator = ( ClothoidCurve const & s )
    { copy(s) ; return *this ; }

    valueType getX0()      const { return x0 ; }
    valueType getY0()      const { return y0 ; }
    valueType getTheta0()  const { return theta0 ; }

    valueType getKappa()   const { return k ; }
    valueType getKappa_D() const { return dk ; }
    valueType getSmin()    const { return s_min ; }
    valueType getSmax()    const { return s_max ; }
    valueType getL()       const { return s_max-s_min ; }

    valueType getThetaBegin() const { return theta0 + s_min * ( k + s_min * dk / 2 ) ; }
    valueType getThetaEnd()   const { return theta0 + s_max * ( k + s_max * dk / 2 ) ; }
    valueType getKappaBegin() const { return k + s_min * dk ; }
    valueType getKappaEnd()   const { return k + s_max * dk ; }

    //! construct a clothoid with the standard parameters
    void
    setup( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _k,
           valueType _dk,
           valueType _L ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      k      = _k ;
      dk     = _dk ;
      s_min  = 0 ;
      s_max  = _L ;
    }

    void
    setup( valueType _x0,
           valueType _y0,
           valueType _theta0,
           valueType _k,
           valueType _dk,
           valueType _smin,
           valueType _smax ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      k      = _k ;
      dk     = _dk ;
      s_min  = _smin ;
      s_max  = _smax ;
    }

    /*! \brief build a clothoid by solving the hermite G1 problem
     *
     * \param _x0     initial x position            \f$ x_0      \f$
     * \param _y0     initial y position            \f$ y_0      \f$
     * \param _theta0 initial angle                 \f$ \theta_0 \f$
     * \param _x1     final x position              \f$ x_1      \f$
     * \param _y1     final y position              \f$ y_1      \f$
     * \param _theta1 final angle                   \f$ \theta_1 \f$
     */
    void
    setup_G1( valueType _x0,
              valueType _y0,
              valueType _theta0,
              valueType _x1,
              valueType _y1,
              valueType _theta1 ) {
      x0     = _x0 ;
      y0     = _y0 ;
      theta0 = _theta0 ;
      build( _x1, _y1, _theta1 ) ;
      //buildClothoid( _x0, _y0, _theta0, _x1, _y1, _theta1, k, dk, s_max ) ;
      //s_min  = 0 ;
    }


    /*! \brief build a clothoid by solving the forward problem
     *
     * \param _x0     initial x position            \f$ x_0      \f$
     * \param _y0     initial y position            \f$ y_0      \f$
     * \param _theta0 initial angle                 \f$ \theta_0 \f$
     * \param _k      initial curvature             \f$ \kappa_0 \f$
     * \param _x1     final x position              \f$ x_1      \f$
     * \param _y1     final y position              \f$ y_1      \f$
     */
    bool
    setup_forward( valueType _x0,
                   valueType _y0,
                   valueType _theta0,
                   valueType _k,
                   valueType _x1,
                   valueType _y1,
                   valueType tol = 1e-8 ) ;

    valueType
    theta( valueType s ) const { return theta0 + s*(k + 0.5*s*dk) ; }

    valueType
    theta_D( valueType s ) const { return k + s*dk ; }

    valueType
    theta_DD( valueType ) const { return dk ; }

    valueType
    theta_DDD( valueType ) const { return 0 ; }

    valueType
    totalLength() const { return s_max-s_min ; }

    valueType
    thetaTotalVariation() const ;

    valueType
    thetaMinMax( valueType & thMin, valueType & thMax ) const ;

    valueType
    deltaTheta() const
    { valueType thMin, thMax ; return thetaMinMax( thMin, thMax ) ; }

    valueType
    curvatureMinMax( valueType & kMin, valueType & kMax ) const ;

    valueType
    curvatureTotalVariation() const ;

    valueType
    integralCurvature2() const ;

    valueType
    integralJerk2() const ;

    valueType
    integralSnap2() const ;

    valueType X( valueType s ) const ;
    valueType Y( valueType s ) const ;

    valueType Xbegin()     const { return X(s_min) ; }
    valueType Ybegin()     const { return Y(s_min) ; }
    valueType KappaBegin() const { return theta_D(s_min) ; }
    valueType ThetaBegin() const { return theta(s_min) ; }

    valueType Xend()     const { return X(s_max) ; }
    valueType Yend()     const { return Y(s_max) ; }
    valueType KappaEnd() const { return theta_D(s_max) ; }
    valueType ThetaEnd() const { return theta(s_max) ; }

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
    trim( valueType s_begin, valueType s_end ) {
      s_min = s_begin ;
      s_max = s_end ;
    }

    //! set the origin of the clothoid to the curvilinear abscissa s0
    void change_origin( valueType s0 ) ;

    //! get the bounding box triangle (if angle variation less that pi/2)
    bool
    bbTriangle( valueType offs,
                valueType p0[2],
                valueType p1[2],
                valueType p2[2] ) const ;

    bool
    bbTriangle( valueType offs, Triangle2D<valueType> & t ) const {
      valueType p0[2], p1[2], p2[2] ;
      bool ok = bbTriangle( offs, p0, p1, p2 ) ;
      if ( ok ) t.setup( p0, p1, p2 ) ;
      return ok ;
    }

    void
    bbSplit( valueType                        split_angle, //!< maximum angle variation
             valueType                        split_size,  //!< maximum height of the triangle
             valueType                        split_offs,  //!< curve offset
             vector<ClothoidCurve>          & c,           //!< clothoid segments
             vector<Triangle2D<valueType> > & t ) const ;  //!< clothoid bounding box

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
    approsimate_collision( valueType             offs,
                           ClothoidCurve const & c,
                           valueType             c_offs,
                           valueType             max_angle,         //!< maximum angle variation
                           valueType             max_size ) const ; //!< curve offset

    friend
    std::ostream &
    operator << ( std::ostream & stream, ClothoidCurve const & c ) ;

    void
    rotate( valueType angle, valueType cx, valueType cy ) ;

    void
    translate( valueType tx, valueType ty )
    { x0 += tx ; y0 += ty ; }

    void
    moveOrigin( valueType newx0, valueType newy0 )
    { x0 = newx0 ; y0 = newy0 ; }

    void
    scale( valueType s ) ;

    void
    reverse() ;

  } ;

  /*\
   |    ____ ____            _           ____
   |   / ___|___ \ ___  ___ | |_   _____|___ \ __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ __) / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __// __/ (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|_____\__,_|_|  \___|
  \*/
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
           valueType & A,
           valueType & A_1,
           valueType & A_2 ) const ;

    void
    evalG( valueType alpha,
           valueType L,
           valueType th,
           valueType k,
           valueType G[2],
           valueType G_1[2],
           valueType G_2[2] ) const ;

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

    void
    setup( valueType x0, valueType y0, valueType theta0, valueType kappa0,
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
   |    ____ ____            _           _____
   |   / ___|___ \ ___  ___ | |_   _____|___ /  __ _ _ __ ___
   |  | |  _  __) / __|/ _ \| \ \ / / _ \ |_ \ / _` | '__/ __|
   |  | |_| |/ __/\__ \ (_) | |\ V /  __/___) | (_| | | | (__
   |   \____|_____|___/\___/|_| \_/ \___|____/ \__,_|_|  \___|
  \*/
  class G2solve3arc {

    ClothoidCurve S0, SM, S1, SG ;

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
    valueType phi, Lscale ;
    valueType th0, th1 ;
    valueType s0, s1 ;

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

    ClothoidCurve const & getS0()    const { return S0 ; }
    ClothoidCurve const & getS1()    const { return S1 ; }
    ClothoidCurve const & getSM()    const { return SM ; }
    ClothoidCurve const & getGuess() const { return SG ; }

    valueType
    totalLength() const {
      return S0.totalLength() +
             S1.totalLength() +
             SM.totalLength() ;
    }

    valueType
    thetaTotalVariation() const {
      return S0.thetaTotalVariation() +
             S1.thetaTotalVariation() +
             SM.thetaTotalVariation() ;
    }

    valueType
    curvatureTotalVariation() const {
      return S0.curvatureTotalVariation() +
             S1.curvatureTotalVariation() +
             SM.curvatureTotalVariation() ;
    }

    valueType
    integralCurvature2() const {
      return S0.integralCurvature2() +
             S1.integralCurvature2() +
             SM.integralCurvature2() ;
    }

    valueType
    integralJerk2() const {
      return S0.integralJerk2() +
             S1.integralJerk2() +
             SM.integralJerk2() ;
    }

    valueType
    integralSnap2() const {
      return S0.integralSnap2() +
             S1.integralSnap2() +
             SM.integralSnap2() ;
    }

    valueType
    thetaMinMax( valueType & thMin, valueType & thMax ) const ;

    valueType
    deltaTheta() const
    { valueType thMin, thMax ; return thetaMinMax( thMin, thMax ) ; }

    valueType
    curvatureMinMax( valueType & kMin, valueType & kMax ) const ;

    valueType getL0() const { return s0/Lscale; }
    valueType getL1() const { return s1/Lscale; }

  } ;

}

#endif

///
/// eof: Clothoid.hh
///
