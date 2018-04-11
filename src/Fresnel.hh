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

  using namespace G2lib ;

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
  FresnelCS( valueType   x,
             valueType & C,
             valueType & S ) ;

  //! Compute Fresnel integrals and its derivatives
  /*!
   * \f[ C(x) = \int_0^x \cos\left(\frac{\pi}{2}t^2\right) dt, \qquad
   *     S(x) = \int_0^x \sin\left(\frac{\pi}{2}t^2\right) dt \f]
   * \param x the input abscissa
   * \param S S[0]=\f$ S(x) \f$, S[1]=\f$ S'(x) \f$, S[2]=\f$ S''(x) \f$
   * \param C C[0]=\f$ C(x) \f$, C[1]=\f$ C'(x) \f$, C[2]=\f$ C''(x) \f$
   */
  void
  FresnelCS( indexType nk,
             valueType x,
             valueType C[],
             valueType S[] ) ;

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
  GeneralizedFresnelCS( indexType nk,
                        valueType a,
                        valueType b,
                        valueType c,
                        valueType intC[],
                        valueType intS[] ) ;

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
  GeneralizedFresnelCS( valueType   a,
                        valueType   b,
                        valueType   c,
                        valueType & intC,
                        valueType & intS ) ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class ClothoidData {
  public:
    valueType x0 ;     //!< initial x coordinate of the clothoid
    valueType y0 ;     //!< initial y coordinate of the clothoid
    valueType theta0 ; //!< initial angle of the clothoid
    valueType kappa0 ;     //!< initial curvature
    valueType dk ;     //!< curvature derivative

    valueType deltaTheta( valueType s ) const
    { return s*(kappa0 + 0.5*s*dk) ; }

    valueType theta( valueType s ) const
    { return theta0 + s*(kappa0 + 0.5*s*dk) ; }

    valueType kappa( valueType s ) const
    { return kappa0 + s*dk ; }

    valueType X( valueType s ) const ;
    valueType Y( valueType s ) const ;

    void
    eval( valueType   s,
          valueType & theta,
          valueType & kappa,
          valueType & x,
          valueType & y ) const ;

    void
    eval( valueType   s,
          valueType & x,
          valueType & y ) const ;

    void
    eval_D( valueType   s,
            valueType & x_D,
            valueType & y_D ) const ;

    void
    eval_DD( valueType   s,
             valueType & x_DD,
             valueType & y_DD ) const ;

    void
    eval_DDD( valueType   s,
              valueType & x_DDD,
              valueType & y_DDD ) const ;

    void
    eval( valueType   s,
          valueType   offs,
          valueType & x,
          valueType & y ) const ;

    void
    eval_D( valueType   s,
            valueType   offs,
            valueType & x_D,
            valueType & y_D ) const ;

    void
    eval_DD( valueType   s,
             valueType   offs,
             valueType & x_DD,
             valueType & y_DD ) const ;

    void
    eval_DDD( valueType   s,
              valueType   offs,
              valueType & x_DDD,
              valueType & y_DDD ) const ;

    void
    eval( valueType s, ClothoidData & C) const ;

    valueType c0x() const { return x0 - (sin(theta0)/kappa0); }
    valueType c0y() const { return y0 + (cos(theta0)/kappa0); }

    void
    Pinfinity( valueType & x, valueType & y, bool plus ) const ;

    void
    reverse( valueType L ) ;

    void
    reverse( valueType L, ClothoidData & out) const ;

    valueType
    split_at_flex( ClothoidData & C0, ClothoidData & C1 ) const ;

    valueType
    aplus( valueType dtheta ) const ;

    valueType
    aminus( valueType dtheta ) const ;

    void
    info( std::ostream & s ) const ;

  } ;

}

#endif

///
/// eof: Fresnel.hh
///
