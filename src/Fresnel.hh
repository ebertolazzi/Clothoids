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
/// file: Fresnel.hh
///

#ifndef FRESNEL_HH
#define FRESNEL_HH

#include "G2lib.hh"

#include <iostream>
#include <cmath>

//! Clothoid computations routine
namespace G2lib {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   |   _____                         _
   |  |  ___| __ ___  ___ _ __   ___| |
   |  | |_ | '__/ _ \/ __| '_ \ / _ \ |
   |  |  _|| | |  __/\__ \ | | |  __/ |
   |  |_|  |_|  \___||___/_| |_|\___|_|
  \*/
  //! Compute Fresnel integrals
  /*!
   * \f[ C(x) = \int_0^x \cos\left(\frac{\pi}{2}t^2\right) dt, \qquad
   *     S(x) = \int_0^x \sin\left(\frac{\pi}{2}t^2\right) dt \f]
   * \param x the input abscissa
   * \param S the value of \f$ S(x) \f$
   * \param C the value of \f$ C(x) \f$
   */
  void
  FresnelCS( real_type   x,
             real_type & C,
             real_type & S );

  //! Compute Fresnel integrals and its derivatives
  /*!
   * \f[ C(x) = \int_0^x \cos\left(\frac{\pi}{2}t^2\right) dt, \qquad
   *     S(x) = \int_0^x \sin\left(\frac{\pi}{2}t^2\right) dt \f]
   * \param x the input abscissa
   * \param S S[0]=\f$ S(x) \f$, S[1]=\f$ S'(x) \f$, S[2]=\f$ S''(x) \f$
   * \param C C[0]=\f$ C(x) \f$, C[1]=\f$ C'(x) \f$, C[2]=\f$ C''(x) \f$
   */
  void
  FresnelCS( int_type  nk,
             real_type x,
             real_type C[],
             real_type S[] );

  /*! \brief Compute the Fresnel integrals
   * \f[
   *   \int_0^1 t^k \cos\left(a\frac{t^2}{2} + b t + c\right) dt,\qquad
   *   \int_0^1 t^k \sin\left(a\frac{t^2}{2} + b t + c\right) dt
   * \f]
   * \param nk   number of momentae to compute
   * \param a    parameter \f$ a \f$
   * \param b    parameter \f$ b \f$
   * \param c    parameter \f$ c \f$
   * \param intC cosine integrals,
   * \param intS sine integrals
   */
  void
  GeneralizedFresnelCS( int_type  nk,
                        real_type a,
                        real_type b,
                        real_type c,
                        real_type intC[],
                        real_type intS[] );

  /*! \brief Compute the Fresnel integrals
   * \f[
   *   \int_0^1 t^k \cos\left(a\frac{t^2}{2} + b t + c\right) dt,\qquad
   *   \int_0^1 t^k \sin\left(a\frac{t^2}{2} + b t + c\right) dt
   * \f]
   * \param a      parameter \f$ a \f$
   * \param b      parameter \f$ b \f$
   * \param c      parameter \f$ c \f$
   * \param intC   cosine integrals,
   * \param intS   sine integrals
   */
  void
  GeneralizedFresnelCS( real_type   a,
                        real_type   b,
                        real_type   c,
                        real_type & intC,
                        real_type & intS );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class ClothoidData {
  public:

    real_type x0;     //!< initial x coordinate of the clothoid
    real_type y0;     //!< initial y coordinate of the clothoid
    real_type theta0; //!< initial angle of the clothoid
    real_type kappa0; //!< initial curvature
    real_type dk;     //!< curvature derivative

    ClothoidData()
    : x0(0)
    , y0(0)
    , theta0(0)
    , kappa0(0)
    , dk(0)
    {}

    real_type
    deltaTheta( real_type s ) const
    { return s*(kappa0 + 0.5*s*dk); }

    //! return angle at curvilinear coordinate `s`
    real_type theta
    ( real_type s ) const
    { return theta0 + s*(kappa0 + 0.5*s*dk); }

    real_type theta_D  ( real_type s ) const { return kappa0 + s*dk; }
    real_type theta_DD ( real_type   ) const { return dk; }
    real_type theta_DDD( real_type   ) const { return 0; }

    //! return curvature at curvilinear coordinate `s`
    real_type kappa    ( real_type s ) const { return kappa0 + s*dk; }
    real_type kappa_D  ( real_type   ) const { return dk; }
    real_type kappa_DD ( real_type   ) const { return 0; }
    real_type kappa_DDD( real_type   ) const { return 0; }

    real_type X( real_type s ) const;
    real_type Y( real_type s ) const;
    real_type X_D( real_type s ) const;
    real_type Y_D( real_type s ) const;
    real_type X_DD( real_type s ) const;
    real_type Y_DD( real_type s ) const;
    real_type X_DDD( real_type s ) const;
    real_type Y_DDD( real_type s ) const;

    real_type
    X( real_type s, real_type t ) const
    { return X(s) + t * nor_x(s); }

    real_type
    Y( real_type s, real_type t ) const
    { return Y(s) + t * nor_y(s); }

    real_type
    X_D( real_type s, real_type t ) const
    { return X_D(s) + t * nor_x_D(s); }

    real_type
    Y_D( real_type s, real_type t ) const
    { return Y_D(s) + t * nor_y_D(s); }

    real_type
    X_DD( real_type s, real_type t ) const
    { return X_DD(s) + t * nor_x_DD(s); }

    real_type
    Y_DD( real_type s, real_type t ) const
    { return Y_DD(s) + t * nor_y_DD(s); }

    real_type
    X_DDD( real_type s, real_type t ) const
    { return X_DDD(s) + t * nor_x_DDD(s); }

    real_type
    Y_DDD( real_type s, real_type t ) const
    { return Y_DDD(s) + t * nor_y_DDD(s); }

    real_type tg0_x() const { return cos(theta0); }
    real_type tg0_y() const { return sin(theta0); }

    real_type tg_x( real_type s ) const { return cos(theta(s)); }
    real_type tg_y( real_type s ) const { return sin(theta(s)); }

    real_type tg_x_D( real_type s ) const;
    real_type tg_y_D( real_type s ) const;

    real_type tg_x_DD( real_type s ) const;
    real_type tg_y_DD( real_type s ) const;

    real_type tg_x_DDD( real_type s ) const;
    real_type tg_y_DDD( real_type s ) const;

    real_type nor0_x() const { return -tg0_y(); }
    real_type nor0_y() const { return tg0_x(); }

    real_type nor_x( real_type s ) const { return -tg_y(s); }
    real_type nor_y( real_type s ) const { return tg_x(s); }

    real_type nor_x_D( real_type s ) const { return -tg_y_D(s); }
    real_type nor_y_D( real_type s ) const { return tg_x_D(s); }

    real_type nor_x_DD( real_type s ) const { return -tg_y_DD(s); }
    real_type nor_y_DD( real_type s ) const { return tg_x_DD(s); }

    real_type nor_x_DDD( real_type s ) const { return -tg_y_DDD(s); }
    real_type nor_y_DDD( real_type s ) const { return tg_x_DDD(s); }

    void tg( real_type s, real_type & tx, real_type & ty ) const;
    void tg_D( real_type s, real_type & tx, real_type & ty ) const;
    void tg_DD( real_type s, real_type & tx, real_type & ty ) const;
    void tg_DDD( real_type s, real_type & tx, real_type & ty ) const;

    void
    nor( real_type s, real_type & nx, real_type & ny ) const
    { tg( s, ny, nx ); ny = -ny; }

    void
    nor_D( real_type s, real_type & nx_D, real_type & ny_D ) const
    { tg_D( s, ny_D, nx_D ); nx_D = -nx_D; }

    void
    nor_DD( real_type s, real_type & nx_DD, real_type & ny_DD ) const
    { tg_DD( s, ny_DD, nx_DD ); nx_DD = -nx_DD; }

    void
    nor_DDD( real_type s, real_type & nx_DDD, real_type & ny_DDD ) const
    { tg_DDD( s, ny_DDD, nx_DDD ); nx_DDD = -nx_DDD; }

    void
    evaluate( real_type   s,
              real_type & theta,
              real_type & kappa,
              real_type & x,
              real_type & y ) const;

    void
    eval( real_type   s,
          real_type & x,
          real_type & y ) const;

    void
    eval_D( real_type   s,
            real_type & x_D,
            real_type & y_D ) const;

    void
    eval_DD( real_type   s,
             real_type & x_DD,
             real_type & y_DD ) const;

    void
    eval_DDD( real_type   s,
              real_type & x_DDD,
              real_type & y_DDD ) const;

    void
    eval( real_type   s,
          real_type   offs,
          real_type & x,
          real_type & y ) const;

    void
    eval_D( real_type   s,
            real_type   offs,
            real_type & x_D,
            real_type & y_D ) const;

    void
    eval_DD( real_type   s,
             real_type   offs,
             real_type & x_DD,
             real_type & y_DD ) const;

    void
    eval_DDD( real_type   s,
              real_type   offs,
              real_type & x_DDD,
              real_type & y_DDD ) const;

    void
    eval( real_type s, ClothoidData & C) const;

    real_type c0x() const { return x0 - (sin(theta0)/kappa0); }
    real_type c0y() const { return y0 + (cos(theta0)/kappa0); }

    void
    Pinfinity( real_type & x, real_type & y, bool plus ) const;

    void
    reverse( real_type L );

    void
    reverse( real_type L, ClothoidData & out) const;

    void
    rotate( real_type angle, real_type cx, real_type cy );

    void
    origin_at( real_type s_origin );

    real_type
    split_at_flex( ClothoidData & C0, ClothoidData & C1 ) const;

    real_type
    aplus( real_type dtheta ) const;

    bool
    bbTriangle( real_type   L,
                real_type & xx0, real_type & yy0,
                real_type & xx1, real_type & yy1,
                real_type & xx2, real_type & yy2 ) const;

    bool
    bbTriangle( real_type   L,
                real_type   offs,
                real_type & xx0, real_type & yy0,
                real_type & xx1, real_type & yy1,
                real_type & xx2, real_type & yy2 ) const;

    int
    build_G1( real_type   x0,
              real_type   y0,
              real_type   theta0,
              real_type   x1,
              real_type   y1,
              real_type   theta1,
              real_type   tol,
              real_type & L,
              bool        compute_deriv = false,
              real_type   L_D[2]        = nullptr,
              real_type   k_D[2]        = nullptr,
              real_type   dk_D[2]       = nullptr );

    bool
    build_forward( real_type   x0,
                   real_type   y0,
                   real_type   theta0,
                   real_type   kappa0,
                   real_type   x1,
                   real_type   y1,
                   real_type   tol,
                   real_type & L );

    void
    info( ostream_type & s ) const;

  };

}

#endif

///
/// eof: Fresnel.hh
///
