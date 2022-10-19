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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Trichotomy<Real>::eval( Real a, Real b ) {
    m_a  = a;       m_fa = m_function->eval(m_a);
    m_b  = b;       m_fb = m_function->eval(m_b);
    m_x3 = (a+b)/2; m_f3 = m_function->eval(m_x3);
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
    m_num_fun_eval  = 0;
    for ( m_num_iter_done = 0; m_num_iter_done <= m_max_iter; ++m_num_iter_done ) {
      m_converged = bracketing();
      if ( m_converged ) break;
      if ( m_num_fun_eval >= m_max_fun_eval ) break;
    }
    return m_x3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Trichotomy<Real>::search( Real x, Real delta ) {
    m_a  = x-delta; m_fa = m_function->eval(m_a);
    m_b  = x+delta; m_fb = m_function->eval(m_b);
    m_x3 = x;       m_f3 = m_function->eval(m_x3);
    while ( m_f3 > m_fa || m_f3 > m_fb ) {
      if ( m_fa < m_fb ) {
        m_b  = m_x3; m_fb = m_f3;
        m_x3 = m_a;  m_f3 = m_fa;
        m_a  = m_a - (m_b-m_x3);
        m_fa = m_function->eval(m_a);
      } else {
        m_a  = m_x3; m_fa = m_f3;
        m_x3 = m_b;  m_f3 = m_fb;
        m_b  = m_b + (m_a-m_x3);
        m_fb = m_function->eval(m_b);
      }
    }
    return this->minimize();
  }

  template class Trichotomy<float>;
  template class Trichotomy<double>;
}
