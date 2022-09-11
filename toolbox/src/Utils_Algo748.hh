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
/// file: Utils_Algo748.hh
///

#pragma once

#ifndef UTILS_ALGO748_dot_HH
#define UTILS_ALGO748_dot_HH

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
  //      _    _           _____ _  _    ___
  //     / \  | | __ _  __|___  | || |  ( _ )
  //    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
  //   / ___ \| | (_| | (_) / / |__   _| (_) |
  //  /_/   \_\_|\__, |\___/_/     |_|  \___/
  //             |___/
  */

  //!
  //! Class function for `Algo748` class
  //!
  template <typename Real>
  class Algo748_base_fun {
  public:
    virtual Real eval( Real x ) const = 0;
    Real operator () ( Real x ) const { return this->eval(x); }
  };

  template <typename Real, typename PFUN>
  class Algo748_fun : public Algo748_base_fun<Real> {
    PFUN m_fun;
  public:
    explicit Algo748_fun( PFUN pfun ) : m_fun(pfun) {}
    Real eval( Real x ) const override { return m_fun(x); };
  };

  //!
  //! Implementation of:
  //!
  //! - **G. E. Alefeld, Florian A Potra, Yixun Shi**,
  //!   *Algorithm 748: enclosing zeros of continuous functions*,
  //!   ACM Transactions on Mathematical Software, vol 21, N.3, 1995
  //!
  template <typename Real>
  class Algo748 {

    typedef int Integer;

    Integer m_num_iter_done = 0;
    Integer m_num_fun_eval  = 0;
    Real    m_mu            = Real(0.5);
    Real    m_tolerance     = pow(machine_eps<Real>(),Real(2./3.));

    bool m_converged     = false;

    Real m_a  = 0;
    Real m_b  = 0;
    Real m_c  = 0;
    Real m_d  = 0;
    Real m_e  = 0;
    Real m_fa = 0;
    Real m_fb = 0;
    Real m_fc = 0;
    Real m_fd = 0;
    Real m_fe = 0;

    Algo748_base_fun<Real> * m_function = nullptr;

    bool bracketing();
    void set_tolerance( Real tol );
    Real pzero();
    Real newton_quadratic( Integer niter );

    Real evaluate( Real x ) { ++m_num_fun_eval; return m_function->eval(x); };

    bool all_different( Real a, Real b, Real c, Real d ) const;

    Real eval( Real a, Real b );

  public:

    Algo748() = default;
    ~Algo748() = default;

    Real
    eval( Real a, Real b, Algo748_base_fun<Real> * fun ) {
      m_function = fun;
      return this->eval(a,b);
    }

    template <typename PFUN>
    Real
    eval( Real a, Real b, PFUN pfun ) {
      Algo748_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->eval(a,b);
    }

    Integer used_iter()    const { return m_num_iter_done; }
    Integer num_fun_eval() const { return m_num_fun_eval; }
    bool converged() const { return m_converged; }

  };

  #ifndef UTILS_OS_WINDOWS
  extern template class Algo748<float>;
  extern template class Algo748<double>;
  #endif

}

#endif

///
/// EOF: Utils_Algo748.hh
///
