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

}

#endif

///
/// eof: Fresnel.hh
///
