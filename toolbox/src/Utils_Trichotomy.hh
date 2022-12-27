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
/// file: Utils_Trichotomy.hh
///

#pragma once

#ifndef UTILS_TRICHOTOMY_dot_HH
#define UTILS_TRICHOTOMY_dot_HH

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "Utils.hh"

namespace Utils {

  using std::pow;
  using std::abs;

  //!
  //! Class function for `Trichotomy` class
  //!
  template <typename Real>
  class Trichotomy_base_fun {
  public:
    virtual Real eval( Real x ) const = 0;
    Real operator () ( Real x ) const { return this->eval(x); }
  };

  template <typename Real, typename PFUN>
  class Trichotomy_fun : public Trichotomy_base_fun<Real> {
    PFUN m_fun;
  public:
    explicit Trichotomy_fun( PFUN pfun ) : m_fun(pfun) {}
    Real eval( Real x ) const override { return m_fun(x); };
  };

  //!
  //!  A new zero-order 1-D optimization algorithm: trichotomy method
  //!  Alena Antonova, Olga Ibryaeva
  //!  arXiv:1903.07117
  //! https://doi.org/10.48550/arXiv.1903.07117
  //!
  template <typename Real>
  class Trichotomy {

    using Integer = int;

    Integer m_num_iter_done{0};
    Integer m_num_fun_eval{0};
    Integer m_max_iteration{100};
    Integer m_max_fun_evaluation{1000};
    Real    m_tolerance{pow(machine_eps<Real>(),Real(2./3.))};

    bool m_converged{false};

    Real m_a{0},  m_fa{0};
    Real m_x1{0}, m_f1{0};
    Real m_x2{0}, m_f2{0};
    Real m_x3{0}, m_f3{0};
    Real m_x4{0}, m_f4{0};
    Real m_x5{0}, m_f5{0};
    Real m_b{0},  m_fb{0};

    Trichotomy_base_fun<Real> * m_function{nullptr};

    bool bracketing();
    Real minimize();
    void set_tolerance( Real tol );
    Real evaluate( Real x ) { ++m_num_fun_eval; return m_function->eval(x); };
    Real eval( Real a, Real b );
    Real search( Real x, Real delta );

  public:

    Trichotomy() UTILS_DEFAULT;
    ~Trichotomy() UTILS_DEFAULT;

    Real
    eval( Real a, Real b, Trichotomy_base_fun<Real> * fun ) {
      m_function = fun;
      return this->eval( a, b );
    }

    template <typename PFUN>
    Real
    eval2( Real a, Real b, PFUN pfun ) {
      Trichotomy_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->eval(a,b);
    }

    Real
    search( Real x, Real delta, Trichotomy_base_fun<Real> * fun ) {
      m_function = fun;
      return this->search( x, delta );
    }

    template <typename PFUN>
    Real
    search2( Real x, Real delta, PFUN pfun ) {
      Trichotomy_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->search( x, delta );
    }

    Integer used_iter()    const { return m_num_iter_done; }
    Integer num_fun_eval() const { return m_num_fun_eval; }
    Real    tolerance()    const { return m_tolerance; }
    bool    converged()    const { return m_converged; }

    void set_max_iterations( Integer mit );
    void set_max_fun_evaluation( Integer mfev );

  };

  #ifndef UTILS_OS_WINDOWS
  extern template class Trichotomy<float>;
  extern template class Trichotomy<double>;
  #endif

}

#endif

///
/// EOF: Utils_Trichotomy.hh
///
