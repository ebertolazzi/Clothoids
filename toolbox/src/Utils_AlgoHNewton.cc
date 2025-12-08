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
|      Università degli Studi di Trento                                    |
|      email: enrico.bertolazzi@unitn.it                                   |
|                                                                          |
\*--------------------------------------------------------------------------*/

#include "Utils_AlgoHNewton.hh"

#include <vector>

#include "Utils_fmt.hh"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils
{

  using std::abs;
  using std::ceil;
  using std::log2;
  using std::max;
  using std::min;
  using std::signbit;
  using std::swap;

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  AlgoHNewton<Real>::set_max_iterations( Integer mit )
  {
    UTILS_ASSERT( mit > 0, "AlgoHNewton::set_max_iterations({}) argument must be >0\n", mit );
    m_max_iteration = mit;
  }

  template <typename Real>
  void
  AlgoHNewton<Real>::set_tolerance( Real B )
  {
    Real eps{ 2 * machine_eps<Real>() };
    m_tolerance = eps + 2 * abs( B ) * eps;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoHNewton<Real>::eval( Real a, Real b )
  {
    m_iteration_count        = 0;
    m_fun_evaluation_count   = 0;
    m_fun_D_evaluation_count = 0;

    m_a  = a;
    m_fa = this->evaluate( m_a );
    m_b  = b;
    m_fb = this->evaluate( m_b );

    // check if solution can exists
    if ( m_fa * m_fb > 0 ) return m_a;
    return eval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoHNewton<Real>::eval( Real a, Real b, Real amin, Real bmax )
  {
    m_iteration_count        = 0;
    m_fun_evaluation_count   = 0;
    m_fun_D_evaluation_count = 0;

    m_a  = a;
    m_fa = this->evaluate( m_a );
    m_b  = b;
    m_fb = this->evaluate( m_b );

    // try to enlarge interval
    while ( m_fa * m_fb > 0 )
    {
      if ( is_finite( m_fa ) && m_a > amin )
      {
        m_a -= m_b - m_a;
        m_fa = this->evaluate( m_a );
      }
      else if ( is_finite( m_fb ) && m_b < bmax )
      {
        m_b += m_b - m_a;
        m_fb = this->evaluate( m_b );
      }
      else
      {
        break;
      }
    }
    return eval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoHNewton<Real>::p_zero2() const
  {
    Real D1{ m_fb - m_fa };
    Real theta{ ( m_d - m_a ) / m_ba };
    Real D2{ m_fd - m_fa };
    Real A{ D2 / theta - D1 };
    Real B{ D1 * theta - D2 / theta };
    Real C{ m_fa * ( theta - 1 ) };
    Real M{ max( max( abs( A ), abs( B ) ), abs( C ) ) };
    if ( C < 0 ) M = -M;
    A /= M;
    B /= M;
    C /= M;
    Real D{ B * B - 4 * A * C };
    if ( D < 0 ) return -m_fa / D1;
    D = sqrt( D );
    return B < 0 ? 2 * C / ( D - B ) : ( D + B ) / ( 2 * A );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoHNewton<Real>::invp_zero2() const
  {
    //
    // Uses quadratic inverse interpolation of f(x) at a, b, and d to
    // get an approximate root of f(x).
    // Rewritten using divided difference.

    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ m_fd };

    Real D0{ 0 };
    Real D1{ 1 };
    Real D2{ ( m_d - m_a ) / m_ba };

    Real D01{ ( D0 - D1 ) / ( x0 - x1 ) };
    Real D12{ ( D1 - D2 ) / ( x1 - x2 ) };

    Real D012{ ( D01 - D12 ) / ( x0 - x2 ) };

    Real O1{ 0 - x0 };
    Real O2{ ( 0 - x1 ) * O1 };

    Real P0{ D0 };
    Real P1{ P0 + D01 * O1 };
    Real P2{ P1 + D012 * O2 };

    UTILS_ASSERT( is_finite( P2 ),
                  "AlgoHNewton<Real>::pzero(), compute NaN or Inf at\n"
                  "a={} f(a)={}\n"
                  "b={} f(b)={}\n"
                  "d={} f(d)={}\n",
                  m_a, m_fa, m_b, m_fb, m_d, m_fd );

    // CALCULATE THE OUTPUT C.
    return P2;
  }

#if 0

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoHNewton<Real>::invp_zero3() const {
    //
    // Uses cubic inverse interpolation of f(x) at a, b, d, and e to
    // get an approximate root of f(x).
    // Rewritten using divided difference.

    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ m_fd };
    Real x3{ m_fe };

    Real D0{ 0 };
    Real D1{ 1 };
    Real D2{ (m_d-m_a) / m_ba };
    Real D3{ (m_e-m_a) / m_ba };

    Real D01{ (D0-D1)/(x0-x1) };
    Real D12{ (D1-D2)/(x1-x2) };
    Real D23{ (D2-D3)/(x2-x3) };

    Real D012{ (D01-D12)/(x0-x2) };
    Real D123{ (D12-D23)/(x1-x3) };

    Real D0123{ (D012-D123)/(x0-x3) };

    Real O1{ 0-x0 };
    Real O2{ (0-x1)*O1 };
    Real O3{ (0-x2)*O2 };

    Real P0{ D0 };
    Real P1{ P0 + D01   * O1 };
    Real P2{ P1 + D012  * O2 };
    Real P3{ P2 + D0123 * O3 };

    UTILS_ASSERT(
      is_finite(P3),
      "AlgoHNewton<Real>::pzero(), compute NaN or Inf at\n"
      "a={} f(a)={}\n"
      "b={} f(b)={}\n"
      "d={} f(d)={}\n"
      "e={} f(e)={}\n",
      m_a, m_fa,
      m_b, m_fb,
      m_d, m_fd,
      m_e, m_fe
    );

    // CALCULATE THE OUTPUT C.
    return P3;
  }
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*!

    Finds either an exact solution or an approximate solution
    of the equation \f$ f(x)=0 \f$ in the interval \f$ [a,b] \f$.

    1. The first iteration is simply a secant step.

    2. Starting with the second iteration, three steps are taken in each
    iteration.

       a. First two steps are either quadratic interpolation or cubic inverse
    interpolation.

       b. The third step is a double-size secant step.

    If the diameter of the enclosing interval obtained after
    those three steps is larger than \f$ (b-a)/2 \f$,
    then an additional bisection step will be taken.


  */

#if 0

  template <typename Real>
  Real
  AlgoHNewton<Real>::eval() {

    // check for trivial solution
    m_converged = m_fa == 0; if ( m_converged ) return m_a;
    m_converged = m_fb == 0; if ( m_converged ) return m_b;

    //
    // While f(left) or f(right) are infinite perform bisection
    //
    while ( !( is_finite(m_fa) && is_finite(m_fb) ) ) {
      ++m_iteration_count;
      m_c  = (m_a+m_b)/2;
      m_fc = this->evaluate(m_c);
      m_converged = m_fc == 0;
      if ( m_converged ) return m_c;
      if ( m_fa*m_fc < 0 ) {
        // --> [a,c]
        m_b = m_c; m_fb = m_fc;
      } else if ( m_fb*m_fc < 0 ) {
        // --> [c,b]
        m_a = m_c; m_fa = m_fc;
      } else {
        UTILS_ERROR(
          "AlgoHNewton::eval() cannot determine if to choose left or right segment\n"
          "a={} fa={}\n"
          "b={} fb={}\n"
          "c={} fc={}\n",
          m_a, m_fa,
          m_c, m_fb,
          m_a, m_fc
        );
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

    Real abs_fa{ abs(m_fa) };
    Real abs_fb{ abs(m_fb) };

    while ( !m_converged ) {
      ++m_iteration_count;
      Real BA0{ m_b - m_a };

      // Calculates the termination criterion.
      // Stops the procedure if the criterion is satisfied.

      {
        Real c{ abs_fb <= abs_fa ? m_b : m_a };
        this->set_tolerance( c );
        m_converged = BA0 <= m_tolerance;
        if ( m_converged ) return c;
      }

      if ( abs_fa < abs_fb ) {
        Real fa_D{ evaluate_D( m_a )};
        m_c = m_a - m_fa/fa_D;
      } else {
        Real fb_D{ evaluate_D( m_b )};
        m_c = m_b - m_fb/fb_D;
      }

      if ( m_c <= m_a + m_tolerance || m_c >= m_b - m_tolerance ) m_c = m_a + BA0/2;

      m_fc = evaluate( m_c );

      UTILS_ASSERT(
        Utils::is_finite(m_fc),
        "AlgoHNewton::eval() fiund NaN in c\n"
        "a={} fa={}\n"
        "b={} fb={}\n"
        "c={} fc={}\n",
        m_a, m_fa,
        m_c, m_fb,
        m_a, m_fc
      );

      Real abs_fc{ abs( m_fc ) };
      m_converged = abs_fc < m_tolerance;
      if ( m_converged ) return m_c;

      if ( m_fa * m_fc < 0 ) { m_b = m_c; m_fb = m_fc; abs_fb = abs_fc; }
      else                   { m_a = m_c; m_fa = m_fc; abs_fa = abs_fc; }

    }
    // TERMINATES THE PROCEDURE AND RETURN THE "ROOT".
    return abs_fa < abs_fb ? m_a : m_b;
  }

#endif

  template <typename Real>
  Real
  AlgoHNewton<Real>::eval()
  {
    // check for trivial solution
    m_converged = m_fa == 0;
    if ( m_converged ) return m_a;
    m_converged = m_fb == 0;
    if ( m_converged ) return m_b;
    m_tolerance = Utils::machine_eps<Real>() + 2 * max( abs( m_a ), abs( m_b ) ) * Utils::machine_eps<Real>();

    //
    // While f(left) or f(right) are infinite perform bisection
    //
    while ( !( is_finite( m_fa ) && is_finite( m_fb ) ) )
    {
      ++m_iteration_count;
      m_c         = ( m_a + m_b ) / 2;
      m_fc        = this->evaluate( m_c );
      m_converged = m_fc == 0;
      if ( m_converged ) return m_c;
      if ( m_fa * m_fc < 0 )
      {
        // --> [a,c]
        m_b  = m_c;
        m_fb = m_fc;
      }
      else if ( m_fb * m_fc < 0 )
      {
        // --> [c,b]
        m_a  = m_c;
        m_fa = m_fc;
      }
      else
      {
        UTILS_ERROR(
            "AlgoHNewton::eval() cannot determine if to choose left or right "
            "segment\n"
            "a={} fa={}\n"
            "b={} fb={}\n"
            "c={} fc={}\n",
            m_a, m_fa, m_c, m_fb, m_a, m_fc );
      }
      m_converged = ( m_b - m_a ) <= m_tolerance;
      if ( m_converged ) return abs( m_fb ) < abs( m_fa ) ? m_b : m_a;
    }

    m_tolerance = 10 * Utils::machine_eps<Real>() * ( 1 + max( abs( m_a ), abs( m_b ) ) );

    Real abs_fa{ abs( m_fa ) };
    Real abs_fb{ abs( m_fb ) };

    while ( ++m_iteration_count < m_max_iteration )
    {
      // fmt::print( "a={}, b={}, b-a={} iter={}\n", m_a, m_b, m_b-m_a,
      // m_iteration_count );

      // Calculates the termination criterion.
      // Stops the procedure if the criterion is satisfied.

      Real epsi{ max( Real( 1 ), max( abs_fa, abs_fb ) ) * Utils::sqrt_machine_eps<Real>() };

      auto diff = [epsi]( Real a, Real b ) -> bool { return abs( a - b ) > epsi; };

      m_ba        = m_b - m_a;
      m_converged = m_ba <= m_tolerance;
      if ( m_converged ) break;

      if ( abs_fa < abs_fb )
      {
        Real fa_D{ evaluate_D( m_a ) };
        m_c = m_a - m_fa / fa_D;
      }
      else
      {
        Real fb_D{ evaluate_D( m_b ) };
        m_c = m_b - m_fb / fb_D;
      }

      Real delta{ m_ba * m_kappa };
      if ( m_c > m_b - delta || m_c < m_a + delta ) m_c = m_a + m_ba / 2;

      m_fc = evaluate( m_c );

      m_converged = m_fc == 0;
      if ( m_converged ) return m_c;

      if ( m_fa * m_fc < 0 )
      {
        if ( m_c - m_a <= m_b - m_c )
        {
          m_b    = m_c;
          m_fb   = m_fc;
          abs_fb = abs( m_fb );
          continue;
        }
      }
      else
      {
        if ( m_c - m_a >= m_b - m_c )
        {
          m_a    = m_c;
          m_fa   = m_fc;
          abs_fa = abs( m_fa );
          continue;
        }
      }

      bool all_diff{ diff( m_fa, m_fd ) && diff( m_fb, m_fd ) };
      Real x{ all_diff ? invp_zero2() : p_zero2() };
      if ( x < m_kappa )
        x = m_kappa;
      else if ( x > 1 - m_kappa )
        x = 1 - m_kappa;
      m_d  = m_a + x * m_ba;
      m_fd = evaluate( m_d );

      m_converged = m_fd == 0;
      if ( m_converged ) return m_d;

      if ( m_c > m_d )
      {
        swap( m_c, m_d );
        swap( m_fc, m_fd );
      }

      if ( m_fc * m_fd < 0 )
      {
        m_a    = m_c;
        m_fa   = m_fc;
        abs_fa = abs( m_fa );
        m_b    = m_d;
        m_fb   = m_fd;
        abs_fb = abs( m_fb );
      }
      else
      {
        if ( m_fa * m_fc > 0 )
        {
          m_a    = m_d;
          m_fa   = m_fd;
          abs_fa = abs( m_fa );
        }
        else
        {
          m_b    = m_c;
          m_fb   = m_fc;
          abs_fb = abs( m_fb );
        }
      }
    }
    return abs( m_fa ) < abs( m_fb ) ? m_a : m_b;
  }

  template class AlgoHNewton<float>;
  template class AlgoHNewton<double>;
}  // namespace Utils
