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

  /*\
  :|:      _    _           _____ _  _    ___
  :|:     / \  | | __ _  __|___  | || |  ( _ )
  :|:    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
  :|:   / ___ \| | (_| | (_) / / |__   _| (_) |
  :|:  /_/   \_\_|\__, |\___/_/     |_|  \___/
  :|:             |___/
  \*/

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
    if ( a == b ) return false;
    if ( a == c ) return false;
    if ( a == d ) return false;
    if ( b == c ) return false;
    if ( b == d ) return false;
    if ( c == d ) return false;
    return true;
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
      Real tol = 0.7*m_tolerance;
      Real ba  = m_b - m_a;
      if      ( ba  <= 2*tol   ) m_c = m_a+0.5*ba;
      else if ( m_c <= m_a+tol ) m_c = m_a+tol;
      else if ( m_c >= m_b-tol ) m_c = m_b-tol;
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
      if ( abs(m_fb) <= abs(m_fa) ) {
        this->set_tolerance(m_b);
      } else {
        this->set_tolerance(m_a);
      }
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
    // This procedure is a slight modification of Aitken-Neville
    // algorithm for interpolation described by Stoer and Bulirsch
    // in "Introduction to numerical analysis" springer-verlag. new york (1980).
    //
    //Real Q11 = (m_d-m_e)*m_fd/(m_fe-m_fd);
    //Real Q21 = (m_e-m_d)*m_fb/(m_fd-m_fb);
    //Real Q31 = (m_a-m_b)*m_fa/(m_fb-m_fa);
    //Real D21 = (m_b-m_d)*m_fd/(m_fd-m_fb);
    //Real D31 = (m_a-m_b)*m_fb/(m_fb-m_fa);
    //Real Q22 = (D21-Q11)*m_fb/(m_fe-m_fb);
    //Real Q32 = (D31-Q21)*m_fa/(m_fb-m_fa);
    //Real D32 = (D31-Q21)*m_fd/(m_fd-m_fa);
    //Real Q33 = (D32-Q22)*m_fa/(m_fe-m_fa);

    Real Q11 = (m_d-m_e)/(m_fe/m_fd-1);
    Real Q21 = (m_e-m_d)/(m_fd/m_fb-1);
    Real Q31 = (m_a-m_b)/(m_fb/m_fa-1);
    Real D21 = (m_b-m_d)/(1-m_fb/m_fd);
    Real D31 = (m_a-m_b)/(1-m_fa/m_fb);
    Real Q22 = (D21-Q11)/(m_fe/m_fb-1);
    Real Q32 = (D31-Q21)/(m_fb/m_fa-1);
    Real D32 = (D31-Q21)/(1-m_fa/m_fd);
    Real Q33 = (D32-Q22)/(m_fe/m_fa-1);

    // CALCULATE THE OUTPUT C.
    return m_a+(Q31+Q32+Q33);
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
    m_num_iter_done = 0;
    m_num_fun_eval  = 0;

    m_a = a; m_fa = this->evaluate(m_a);
    m_b = b; m_fb = this->evaluate(m_b);

    // check for trivial solution
    m_converged = m_fa == 0;
    if ( m_converged ) return m_a;
    m_converged = m_fb == 0;
    if ( m_converged ) return m_b;

    // check if solution can exists
    if ( m_fa*m_fb > 0 ) return 0;

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

    // Until f(left) or f(right) are infinite perform bisection
    while ( !( isfinite(m_fa) && isfinite(m_fb) ) ) {
      ++m_num_iter_done;
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
    Real ba  = m_b-m_a;
    Real fba = m_fb-m_fa;
    Real R   = ba/fba;
    m_c = abs(m_fb) < abs(m_fa) ? m_b+m_fb*R : m_a-m_fa*R;
    // impedisce m_c troppo vicino ad m_a o m_b
    Real bmax = m_b-0.1*ba;
    Real amin = m_a+0.1*ba;
    if      ( m_c < amin ) m_c = amin;
    else if ( m_c > bmax ) m_c = bmax;
    //
    // Call "bracketing" to get a shrinked enclosing interval as
    // well as to update the termination criterion.
    // Stop the procedure if the criterion is satisfied or the
    // exact solution is obtained.
    //
    if ( bracketing() ) return m_a;
    // Iteration starts.
    // The enclosing interval before executing the iteration is recorded as [a0, b0].
    m_converged = false;

    // ITERATION STARTS. THE ENCLOSING INTERVAL BEFORE EXECUTING THE
    // ITERATION IS RECORDED AS [A0, B0].
    while ( !m_converged ) {
      ++m_num_iter_done;
      Real A0 = m_a;
      Real B0 = m_b;

      // Calculates the termination criterion.
      // Stops the procedure if the criterion is satisfied.

      if ( abs(m_fb) <= abs(m_fa) ) this->set_tolerance(m_b);
      else                          this->set_tolerance(m_a);

      m_converged = (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return (m_a+m_b)/2;

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
      if ( m_converged ) return m_a;
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
      m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return m_a;

      m_e  = m_d;
      m_fe = m_fd;
      // Takes the double-size secant step.
      Real u, fu;
      if ( abs(m_fa) < abs(m_fb) ) { u = m_a; fu = m_fa; }
      else                         { u = m_b; fu = m_fb; }
      {
        Real ba = m_b-m_a;
        m_c = u-2*(fu/(m_fb-m_fa))*ba;
        if ( abs(m_c-u) > 0.5*ba ) m_c = m_a+0.5*ba;
      }

      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
      if ( m_converged ) return m_a;
      //
      // Determines whether an additional bisection step is needed.
      // Takes it if necessary.
      //
      if ( (m_b-m_a) < m_mu*(B0-A0) ) continue;

      m_e  = m_d;
      m_fe = m_fd;
      //
      // Call subroutine "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion. stop the procedure
      // if the criterion is satisfied or the exact solution is obtained.
      //
      m_c = m_a+0.5*(m_b-m_a);
      m_converged = this->bracketing() || (m_b-m_a) <= m_tolerance;
    }
    // TERMINATES THE PROCEDURE AND RETURN THE "ROOT".
    return m_a;
  }

  template class Algo748<float>;
  template class Algo748<double>;

}
