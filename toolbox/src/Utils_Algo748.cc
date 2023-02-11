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

#include "Utils_Algo748.hh"

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

  /*\
  :|:      _    _           _____ _  _    ___
  :|:     / \  | | __ _  __|___  | || |  ( _ )
  :|:    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
  :|:   / ___ \| | (_| | (_) / / |__   _| (_) |
  :|:  /_/   \_\_|\__, |\___/_/     |_|  \___/
  :|:             |___/
  \*/

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  Algo748<Real>::set_max_iterations( Integer mit ) {
    UTILS_ASSERT(
      mit > 0,
      "Algo748::set_max_iterations({}) argument must be >0\n", mit
    );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  Algo748<Real>::set_max_fun_evaluation( Integer mfev ) {
    UTILS_ASSERT(
      mfev > 0,
      "Algo748::set_max_fun_evaluation({}) argument must be >0\n", mfev
    );
    m_max_fun_evaluation = mfev;
  }

  template <typename Real>
  void
  Algo748<Real>::set_tolerance( Real B ) {
    Real eps = 2*Utils::machine_eps<Real>();
    m_tolerance = eps + 2*abs(B)*eps;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  bool
  Algo748<Real>::all_different( Real a, Real b, Real c, Real d ) const {
    return a != b &&
           a != c &&
           a != d &&
           b != c &&
           b != d &&
           c != d;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  bool
  Algo748<Real>::bracketing() {
    // Given current enclosing interval [a,b] and a number c in (a,b):
    //
    //  a) if f(c)=0 then sets the output a=c.
    //  b) Otherwise determines the new enclosing interval:
    //     [a,b]=[a,c] or [a,b]=[c,b].
    //     also updates the termination criterion corresponding
    //     to the new enclosing interval.
    //
    // Adjust c if (b-a) is very small or if c is very close to a or b.
    {
      Real tol = Real(0.7)*m_tolerance;
      Real hba = (m_b - m_a)/2;
      if      ( hba <= tol     ) m_c = m_a + hba;
      else if ( m_c <= m_a+tol ) m_c = m_a + tol;
      else if ( m_c >= m_b-tol ) m_c = m_b - tol;
    }

    m_fc = this->evaluate( m_c );

    // If f(c)=0, then set a=c and return.
    // This will terminate the procedure.

    if ( m_fc == 0 ) {
      m_a = m_c; m_fa = 0;
      m_d = 0;   m_fd = 0;
      return true;
    } else {
      // If f(c) is not zero, then determine the new enclosing interval.
      if ( m_fa*m_fc < 0 ) {
        // D <-- B <-- C
        m_d = m_b; m_fd = m_fb;
        m_b = m_c; m_fb = m_fc;
      } else {
        // D <-- A <-- C
        m_d = m_a; m_fd = m_fa;
        m_a = m_c; m_fa = m_fc;
      }
      // update the termination criterion according to the new enclosing interval.
      if ( abs(m_fb) <= abs(m_fa) ) this->set_tolerance(m_b);
      else                          this->set_tolerance(m_a);
      return false;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Algo748<Real>::pzero() {
    //
    // Uses cubic inverse interpolation of f(x) at a, b, d, and e to
    // get an approximate root of f(x).
    // Rewritten using divided difference.

    Real D1 = m_b-m_a;
    Real D2 = m_d-m_a;
    Real D3 = m_e-m_a;

    Real DD0 = D1/(m_fb-m_fa);
    Real DD1 = (D1-D2)/(m_fb-m_fd);
    Real DD2 = (D2-D3)/(m_fd-m_fe);

    Real DDD0  = (DD0-DD1)/(m_fa-m_fd);
    Real DDD1  = (DD1-DD2)/(m_fb-m_fe);

    Real DDDD0 = (DDD0-DDD1)/(m_fa-m_fe);

    Real c = m_a - m_fa*(DD0-m_fb*(DDD0-m_fd*DDDD0));

    Real tol = Real(0.7)*m_tolerance;
    if ( c <= m_a+tol || c >= m_b-tol ) c = (m_a+m_b)/2;

    // CALCULATE THE OUTPUT C.
    return c;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Algo748<Real>::newton_quadratic( Integer niter ) {
    // Uses `niter` newton steps to approximate the zero in (a,b) of the
    // quadratic polynomial interpolating f(x) at a, b, and d.
    // Safeguard is used to avoid overflow.

    Real A0 = m_fa;
    Real A1 = (m_fb-m_fa)/(m_b-m_a);
    Real A2 = ((m_fd-m_fb)/(m_d-m_b)-A1)/(m_d-m_a);

    // Safeguard to avoid overflow.
    if ( A2 == 0 ) return m_a-A0/A1;

    // Determine the starting point of newton steps.
    Real c = A2*m_fa > 0 ? m_a : m_b;

    // Start the safeguarded newton steps.
    bool ok = true;
    for ( Integer i=0; i < niter && ok ; ++i ) {
      Real PC  = A0+(A1+A2*(c-m_b))*(c-m_a);
      Real PDC = A1+A2*((2*c)-(m_a+m_b));
      ok = PDC != 0;
      if ( ok ) c -= PC/PDC;
    }

    if ( ok ) return c;
    else      return m_a-A0/A1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Algo748<Real>::eval( Real a, Real b ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    m_a = a; m_fa = this->evaluate(m_a);
    m_b = b; m_fb = this->evaluate(m_b);

    // check if solution can exists
    if ( m_fa*m_fb > 0 ) return 0;
    else                 return eval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Algo748<Real>::eval( Real a, Real b, Real amin, Real bmax ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    m_a = a; m_fa = this->evaluate(m_a);
    m_b = b; m_fb = this->evaluate(m_b);

    // try to enlarge interval
    while ( m_fa*m_fb > 0 ) {
      if ( Utils::is_finite(m_fa) && m_a > amin ) {
        m_a -= m_b - m_a;
        m_fa = this->evaluate(m_a);
      } else if ( Utils::is_finite(m_fb) && m_b < bmax ) {
        m_b += m_b - m_a;
        m_fb = this->evaluate(m_b);
      } else {
        break;
      }
    }
    return eval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Algo748<Real>::eval() {

    // check for trivial solution
    m_converged = m_fa == 0; if ( m_converged ) return m_a;
    m_converged = m_fb == 0; if ( m_converged ) return m_b;

    // Finds either an exact solution or an approximate solution
    // of the equation f(x)=0 in the interval [a,b].
    //
    // At the beginning of each iteration, the current enclosing interval
    // is recorded as [a0,b0].
    // The first iteration is simply a secant step.
    // Starting with the second iteration, three steps are taken in each iteration.
    // First two steps are either quadratic interpolation
    // or cubic inverse interpolation.
    // The third step is a double-size secant step.
    // If the diameter of the enclosing interval obtained after
    // those three steps is larger than 0.5*(b0-a0),
    // then an additional bisection step will be taken.

    m_e  = NaN<Real>();
    m_fe = NaN<Real>(); // Dumb values

    //
    // While f(left) or f(right) are infinite perform bisection
    //
    while ( !( isfinite(m_fa) && isfinite(m_fb) ) ) {
      ++m_iteration_count;
      m_c  = (m_a+m_b)/2;
      m_fc = this->evaluate(m_c);
      m_converged = m_fc == 0;
      if ( m_converged ) return m_c;
      if ( m_fa*m_fc < 0 ) {
        // --> [a,c]
        m_b = m_c; m_fb = m_fc;
      } else {
        // --> [c,b]
        m_a = m_c; m_fa = m_fc;
      }
      m_converged = (m_b-m_a) <= m_tolerance;
      if ( abs(m_fb) <= abs(m_fa) ) {
        this->set_tolerance(m_b);
        if ( m_converged ) return m_b;
      } else {
        this->set_tolerance(m_a);
        if ( m_converged ) return m_a;
      }
    }
    {
      Real ba  = m_b  - m_a;
      Real R   = ba/(m_fb - m_fa);
      m_c = abs(m_fb) < abs(m_fa) ? m_b+m_fb*R : m_a-m_fa*R;
      // impedisce m_c troppo vicino ad m_a o m_b
      Real delta = m_interval_shink*ba, amin, bmax;
      if      ( m_c < (amin=m_a+delta) ) m_c = amin;
      else if ( m_c > (bmax=m_b-delta) ) m_c = bmax;
    }
    //
    // Call "bracketing" to get a shrinked enclosing interval as
    // well as to update the termination criterion.
    // Stop the procedure if the criterion is satisfied or the
    // exact solution is obtained.
    //
    m_converged = bracketing();
    if ( m_converged ) return m_a;
    // Iteration starts.
    // The enclosing interval before executing the iteration is recorded as [a0, b0].
    m_converged = false;

    // ITERATION STARTS. THE ENCLOSING INTERVAL BEFORE EXECUTING THE
    // ITERATION IS RECORDED AS [A0, B0].
    while ( !m_converged ) {
      ++m_iteration_count;
      Real BA0 = m_b-m_a;

      // Calculates the termination criterion.
      // Stops the procedure if the criterion is satisfied.

      {
        Real abs_fa = abs(m_fa);
        Real abs_fb = abs(m_fb);
        Real c      = abs_fb <= abs_fa ? m_b : m_a;
        this->set_tolerance( c );
        m_converged = BA0 <= m_tolerance;
        if ( m_converged ) return c;
      }

      //
      // Starting with the second iteration, in the first two steps, either
      // quadratic interpolation is used by calling the subroutine "newtonquadratic"
      // or the cubic inverse interpolation is used by calling the subroutine
      // "pzero". in the following, if "prof" is not equal to 0, then the
      // four function values "fa", "fb", "fd", and "fe" are distinct, and
      // hence "pzero" will be called.
      //

      bool do_newton_quadratic = false;
      if ( !is_finite(m_fe) ) {
        do_newton_quadratic = true;
      } else if ( !this->all_different( m_fa, m_fb, m_fd, m_fe ) ) {
        do_newton_quadratic = true;
      } else {
        m_c = this->pzero();
        do_newton_quadratic = (m_c-m_a)*(m_c-m_b) >= 0;
      }
      if ( do_newton_quadratic ) m_c = this->newton_quadratic(2);

      m_e  = m_d;
      m_fe = m_fd;

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return abs(m_fa) < abs(m_fb) ? m_a : m_b;
      if ( !this->all_different( m_fa, m_fb, m_fd, m_fe ) ) {
        do_newton_quadratic = true;
      } else {
        m_c = this->pzero();
        do_newton_quadratic = (m_c-m_a)*(m_c-m_b) >= 0;
      }
      if ( do_newton_quadratic ) m_c = this->newton_quadratic(3);

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //

      {
        m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
        Real abs_fa = abs(m_fa);
        Real abs_fb = abs(m_fb);
        if ( m_converged ) return abs_fa < abs_fb ? m_a : m_b;

        m_e  = m_d;
        m_fe = m_fd;
        // Takes the double-size secant step.
        Real u, fu;
        if ( abs_fa < abs_fb ) { u = m_a; fu = m_fa; }
        else                   { u = m_b; fu = m_fb; }
        Real hba = (m_b-m_a)/2;
        m_c = u-4*(fu/(m_fb-m_fa))*hba;
        if ( abs(m_c-u) > hba ) m_c = m_a + hba;
      }

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return abs(m_fa) < abs(m_fb) ? m_a : m_b;
      //
      // Determines whether an additional bisection step is needed.
      // Takes it if necessary.
      //
      if ( (m_b-m_a) < m_mu*BA0 ) continue;

      m_e  = m_d;
      m_fe = m_fd;
      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      {
        Real ba = m_b-m_a;
        m_c = m_a + ba/2;
        m_converged = this->bracketing() || ba <= m_tolerance;
      }
    }
    // TERMINATES THE PROCEDURE AND RETURN THE "ROOT".
    return abs(m_fa) < abs(m_fb) ? m_a : m_b;
  }

  template class Algo748<float>;
  template class Algo748<double>;
}
