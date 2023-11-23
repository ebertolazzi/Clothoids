/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
/// file: Utils_zeros.hh
///

#pragma once

#ifndef UTILS_ZEROS_dot_HH
#define UTILS_ZEROS_dot_HH

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "Utils.hh"

namespace Utils {

  using std::pow;
  using std::abs;

  /*
  //   __________ ____   ___  ____
  //  |__  / ____|  _ \ / _ \/ ___|
  //    / /|  _| | |_) | | | \___ \
  //   / /_| |___|  _ <| |_| |___) |
  //  /____|_____|_| \_\\___/|____/
  */

  //!
  //! Class function for `Zeros` class
  //!
  template <typename Real>
  class Zeros_base_fun {
  public:
    virtual Real eval    ( Real ) const = 0;
    virtual Real eval_D  ( Real ) const { UTILS_ERROR("Zeros_base_fun derivative not defined"); };
    virtual Real eval_DD ( Real ) const { UTILS_ERROR("Zeros_base_fun second derivative not defined"); };
    virtual Real eval_DDD( Real ) const { UTILS_ERROR("Zeros_base_fun third derivative not defined"); };
    Real operator () ( Real x ) const { return this->eval(x); }
    Real D           ( Real x ) const { return this->eval_D(x); }
    Real DD          ( Real x ) const { return this->eval_DD(x); }
    Real DDD         ( Real x ) const { return this->eval_DDD(x); }
  };

  //!
  //! Implementation of:
  //!
  //! - Newton
  //!
  template <typename Real>
  class Zeros {

    using Integer = int;

    Integer m_max_fun_evaluation{200}; // max number of function evaluations
    Integer m_max_iteration{100};       // max number of iterations
    Real    m_tolerance{pow(machine_eps<Real>(),Real(2./3.))};
    bool    m_converged{false};

    mutable Integer m_iteration_count{0};    // explore iteration counter
    mutable Integer m_fun_evaluation_count{0};

    //static Real Q( Real t ) { return 1/(1-2*t); }
    static Real Q( Real t ) { return 1+2*t; }

    static Real W( Real t, Real s ) {
      Real t2{t*t};
      return t2*(1-4*t)+(4*s+2)*t+s+1;
    }

    static Real H( Real t, Real s, Real u ) {
      Real t1  = t*t;
      Real t2  = t1*t1;
      Real t8  = s*s;
      Real t17 = s*t8;
      Real t23 = 2.0*u;
      return ( (8*u+6*t2+4)*s-(6*t8+4*(s+u+1))*t1 + 2*t8 - 4*t17 + t23 + 2 )*t +
               t1*(t8+s+u+1) + (1-3*t2+t23)*s + u - t17 + 1;
    }

    static Real J( Real t, Real s, Real u, Real v ) {
      Real t1  = s*s;
      Real t2  = t1*t1;
      Real t17 = t*t;
      Real t22 = u*u;
      Real t32 = t17*t17;
      Real t34 = t*t32;
      Real t37 = t*t17;
      Real t46 = 1.0+v;
      Real t65 = u+1.0+v;
      Real t76 = (-2.0*t22+u+4.0*v+2.0)*u;
      return (-1.0+2.0*t)*(2.0+5.0*t)*u*t*t2 +
             (4.0*t+1.0)*u*s*t2+
             (u*t22-2.0*u*v-u-v-1.0)*(4.0*t17+3.0*t+1.0)*(-1.0+t)-
             8.0*(t22*(t17/2.0-1.0/4.0)+u*(t17*t32-5.0/8.0*t34-3.0/4.0*t32+
             3.0/8.0*t37+3.0/4.0*t17-t/8.0-1.0/4.0)+3.0/4.0*t46*(t+1.0/2.0)*(t-2.0/3.0))*t*t1+
             4.0*(t22*(-3.0/2.0*t-1.0/4.0)+u*(t34-t32-3.0/2.0*t37+t17/4.0-t-1.0/4.0)-
             t46*(t+1.0/4.0))*s*t1+(1.0+v+t65*t17-4.0*t65*t37-3.0*t65*t32+6.0*t65*t34+t76+4.0*(1.0+v+t76)*t)*s;
    }

  public:

    Zeros() UTILS_DEFAULT;
    ~Zeros() UTILS_DEFAULT;

    void set_max_iterations( Integer mit );
    void set_max_fun_evaluation( Integer mfev );
    void set_tolerance( Real tol );

    Real solve_Newton( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Chebyshev( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Halley( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Order4( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Order8( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Order16( Real x_guess, Zeros_base_fun<Real> * fun );
    Real solve_Order32( Real x_guess, Zeros_base_fun<Real> * fun );

    Integer used_iter()    const { return m_iteration_count; }
    Integer num_fun_eval() const { return m_fun_evaluation_count; }
    Real    tolerance()    const { return m_tolerance; }
    bool    converged()    const { return m_converged; }

  };

  #ifndef UTILS_OS_WINDOWS
  extern template class Zeros<float>;
  extern template class Zeros<double>;
  #endif

}

#endif

///
/// EOF: Utils_zeros.hh
///
