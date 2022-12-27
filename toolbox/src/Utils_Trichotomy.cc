 /*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#include "Utils_Trichotomy.hh"

#include <vector>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils {

  using std::isfinite;

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  Trichotomy<Real>::set_max_iterations( Integer mit ) {
    UTILS_ASSERT(
      mit > 0,
      "Trichotomy::set_max_iterations({}) argument must be >0\n", mit
    );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  Trichotomy<Real>::set_max_fun_evaluation( Integer mfev ) {
    UTILS_ASSERT(
      mfev > 0,
      "Trichotomy::set_max_fun_evaluation({}) argument must be >0\n", mfev
    );
    m_max_fun_evaluation = mfev;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Trichotomy<Real>::eval( Real a, Real b ) {
    m_a  = a;       m_fa = m_function->eval(m_a);
    m_b  = b;       m_fb = m_function->eval(m_b);
    m_x3 = (a+b)/2; m_f3 = m_function->eval(m_x3);
    m_num_iter_done = 0;
    m_num_fun_eval  = 3;
    return this->minimize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  bool
  Trichotomy<Real>::bracketing() {
    m_x2 = (m_a+2*m_x3)/3;
    m_f2 = this->evaluate(m_x2);
    if ( m_f2 <= m_f3 ) {
      m_x1 = (m_a+m_x2)/2;
      m_f1 = this->evaluate(m_x1);
      if ( m_f1 <= m_f2 ) {
        m_x3 = m_x1; m_f3 = m_f1;
        m_b  = m_x2; m_fb = m_f2;
      } else {
        m_a  = m_x1; m_fa = m_f1;
        m_b  = m_x3; m_fb = m_f3;
        m_x3 = m_x2; m_f3 = m_f2;
      }
    } else {
      m_x4 = (m_b+2*m_x3)/3;
      m_f4 = this->evaluate(m_x4);
      if ( m_f4 <= m_f3 ) {
        m_x5 = (2*m_b+m_x3)/3;
        m_f5 = this->evaluate(m_x5);
        if ( m_f5 <= m_f4 ) {
          m_a  = m_x4; m_fa = m_f4;
          m_x3 = m_x5; m_f3 = m_f5;
        } else {
          m_a  = m_x3; m_fa = m_f3;
          m_x3 = m_x4; m_f3 = m_f4;
          m_b  = m_x5; m_fb = m_f5;
        }
      } else {
        m_a = m_x2; m_fa = m_f2;
        m_b = m_x4; m_fb = m_f4;
      }
    }
    return m_b-m_a < m_tolerance;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Trichotomy<Real>::minimize() {
    m_num_iter_done = 0;
    m_converged     = false;
    while ( m_num_iter_done++ < m_max_iteration ) {
      m_converged = bracketing();
      if ( m_converged ) break;
      if ( m_num_fun_eval >= m_max_fun_evaluation ) break;
    }
    return m_x3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Trichotomy<Real>::search( Real x, Real delta ) {
    m_num_iter_done = 0;
    m_num_fun_eval  = 3;
    m_a  = x-delta; m_fa = m_function->eval(m_a);
    m_x3 = x;       m_f3 = m_function->eval(m_x3);
    m_b  = x+delta; m_fb = m_function->eval(m_b);
    if ( m_fa < m_fb ) {
      // enter fa < f3 <? fb
      while ( m_fa < m_f3 ) {
        m_b  = m_x3; m_fb = m_f3;
        m_x3 = m_a;  m_f3 = m_fa;
        m_a  = m_x3 - 2*(m_b-m_x3);
        m_fa = this->evaluate(m_a);
        ++m_num_iter_done;
      }
      // exit fa >= f3 < fb
      //      a ======== x3 == b
      m_x2 = (m_a+m_x3)/2;
      m_f2 = this->evaluate(m_x2);
      if ( m_f2 <= m_f3 ) {
        m_b  = m_x3; m_fb = m_f3;
        m_x3 = m_x2; m_f3 = m_f2; // fa >= f3 < fb
      } else {
        m_a = m_x2; m_fa = m_f2; // fa > f3 < fb
      }
    } else {
      // enter fa ?> f3 > fb
      while ( m_f3 > m_fb ) {
        m_a  = m_x3; m_fa = m_f3;
        m_x3 = m_b;  m_f3 = m_fb;
        m_b  = m_x3 + 2*(m_x3-m_a);
        m_fb = this->evaluate(m_b);
        ++m_num_iter_done;
      }
      // exit fa > f3 <= fb
      //      a == x3 ====== b
      m_x4 = (m_x3+m_b)/2;
      m_f4 = this->evaluate(m_x4);
      if ( m_f4 <= m_f3 ) {
        m_a  = m_x3; m_fa = m_f3;
        m_x3 = m_x4; m_f3 = m_f4; // fa >= f3 < fb
      } else {
        m_b = m_x4; m_fb = m_f4; // fa > f3 < fb
      }
    }
    //fmt::print(
    //  "a={} x3={} b={} d1={} d2={} iter={} #fun={}\n",
    //  m_a, m_x3, m_b, m_x3-m_a, m_b-m_x3, m_num_iter_done, m_num_fun_eval
    //);
    return this->minimize();
  }

  template class Trichotomy<float>;
  template class Trichotomy<double>;
}
