/*--------------------------------------------------------------------------*\
|                                                                          |
|  Copyright (C) 2025                                                      |
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

#include "Utils_AlgoBracket.hh"

#include <algorithm>
#include <cmath>
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

  /*
  //   ____                 _        _
  //  | __ ) _ __ __ _  ___| | _____| |_
  //  |  _ \| '__/ _` |/ __| |/ / _ \ __|
  //  | |_) | | | (_| | (__|   <  __/ |_
  //  |____/|_|  \__,_|\___|_|\_\___|\__|
  //
  */

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  AlgoBracket<Real>::set_max_iterations( Integer mit )
  {
    UTILS_ASSERT( mit > 0, "AlgoBracket::set_max_iterations({}) argument must be >0\n", mit );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  AlgoBracket<Real>::set_max_fun_evaluation( Integer mfev )
  {
    UTILS_ASSERT( mfev > 0, "AlgoBracket::set_max_fun_evaluation({}) argument must be >0\n", mfev );
    m_max_fun_evaluation = mfev;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename Real>
  Real
  AlgoBracket<Real>::p_zero2( Real d, Real fd ) const
  {
    Real D1{ m_fb - m_fa };
    Real ba{ m_b - m_a };
    Real theta{ ( d - m_a ) / ba };
    Real D2{ fd - m_fa };
    Real A{ D2 / theta - D1 };
    Real B{ D1 * theta - D2 / theta };
    Real C{ m_fa * ( theta - 1 ) };
    Real M{ max( max( abs( A ), abs( B ) ), abs( C ) ) };
    if ( C < 0 ) M = -M;
    A /= M;
    B /= M;
    C /= M;
    Real D{ B * B - 4 * A * C };
    if ( D >= 0 )
    {
      D = sqrt( D );
      return B < 0 ? 2 * C / ( D - B ) : ( D + B ) / ( 2 * A );
    }
    return -m_fa / D1;
  }

  template <typename Real>
  Real
  AlgoBracket<Real>::invp_zero2( Real d, Real fd ) const
  {
    //
    // Uses quadratic inverse interpolation of f(x) at a, b, and d to
    // get an approximate root of f(x).
    // Rewritten using divided difference.

    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ fd };

    Real D0{ 0 };
    Real D1{ 1 };
    Real D2{ ( d - m_a ) / ( m_b - m_a ) };

    Real D01{ ( D0 - D1 ) / ( x0 - x1 ) };
    Real D12{ ( D1 - D2 ) / ( x1 - x2 ) };

    Real D012{ ( D01 - D12 ) / ( x0 - x2 ) };

    Real O1{ 0 - x0 };
    Real O2{ ( 0 - x1 ) * O1 };

    Real P0{ D0 };
    Real P1{ P0 + D01 * O1 };
    Real P2{ P1 + D012 * O2 };

    UTILS_ASSERT( is_finite( P2 ),
                  "AlgoBracket<Real>::pzero(), compute NaN or Inf at\n"
                  "a={} f(a)={}\n"
                  "b={} f(b)={}\n"
                  "d={} f(d)={}\n",
                  m_a, m_fa, m_b, m_fb, d, fd );

    // CALCULATE THE OUTPUT C.
    return P2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoBracket<Real>::invp_zero3( Real d, Real fd, Real e, Real fe ) const
  {
    //
    // Uses cubic inverse interpolation of f(x) at a, b, d, and e to
    // get an approximate root of f(x).
    // Rewritten using divided difference.

    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ fd };
    Real x3{ fe };

    Real D0{ 0 };
    Real D1{ 1 };
    Real D2{ ( d - m_a ) / ( m_b - m_a ) };
    Real D3{ ( e - m_a ) / ( m_b - m_a ) };

    Real D01{ ( D0 - D1 ) / ( x0 - x1 ) };
    Real D12{ ( D1 - D2 ) / ( x1 - x2 ) };
    Real D23{ ( D2 - D3 ) / ( x2 - x3 ) };

    Real D012{ ( D01 - D12 ) / ( x0 - x2 ) };
    Real D123{ ( D12 - D23 ) / ( x1 - x3 ) };

    Real D0123{ ( D012 - D123 ) / ( x0 - x3 ) };

    Real O1{ 0 - x0 };
    Real O2{ ( 0 - x1 ) * O1 };
    Real O3{ ( 0 - x2 ) * O2 };

    Real P0{ D0 };
    Real P1{ P0 + D01 * O1 };
    Real P2{ P1 + D012 * O2 };
    Real P3{ P2 + D0123 * O3 };

    UTILS_ASSERT( is_finite( P3 ),
                  "AlgoBracket<Real>::pzero(), compute NaN or Inf at\n"
                  "a={} f(a)={}\n"
                  "b={} f(b)={}\n"
                  "d={} f(d)={}\n"
                  "e={} f(e)={}\n",
                  m_a, m_fa, m_b, m_fb, d, fd, e, fe );

    // CALCULATE THE OUTPUT C.
    return P3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoBracket<Real>::eval( Real a, Real b )
  {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

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
  AlgoBracket<Real>::eval( Real a, Real b, Real amin, Real bmax )
  {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

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

  template <typename Real>
  Real
  AlgoBracket<Real>::eval()
  {
    auto check = [this]( Real c, Real fc ) -> void
    {
      UTILS_ASSERT( is_finite( fc ),
                    "AlgoBracket::eval() found Inf or NaN\n"
                    "a={} fa={}\n"
                    "b={} fb={}\n"
                    "c={} fc={}\n",
                    m_a, m_fa, m_b, m_fb, c, fc );
    };

    {
      Real & a{ m_a };
      Real & b{ m_b };
      Real & fa{ m_fa };
      Real & fb{ m_fb };

      UTILS_ASSERT( !is_NaN( fa ) && !is_NaN( fb ),
                    "AlgoBracket::eval() bad initial interval\n"
                    "a = {}, fa = {}\n"
                    "b = {}, fb = {}\n",
                    a, fa, b, fb );

      // check for trivial solution
      m_converged = fa == 0;
      if ( m_converged ) return a;
      m_converged = fb == 0;
      if ( m_converged ) return b;

      //
      // While f(left) or f(right) are infinite perform bisection
      //
      bool ffa, ffb;
      while ( ( ffa = !is_finite( fa ) ) || ( ffb = !is_finite( fb ) ) )
      {
        UTILS_ASSERT( ++m_iteration_count <= m_max_iteration,
                      "AlgoBracket::eval() too many iteration\n"
                      "a={} fa={}\n"
                      "b={} fb={}\n",
                      a, fa, b, fb );

        // tale "mid-point" close to infinite point
        Real c;
        if ( ffa )
          c = ( 15 * a + b ) / 16;
        else if ( ffb )
          c = ( 15 * b + a ) / 16;
        else
          c = ( a + b ) / 2;

        Real fc{ this->evaluate( c ) };

        UTILS_ASSERT( !is_NaN( fc ), "AlgoBracket::eval()\nc = {}, fc = {}\n", c, fc );

        m_converged = fc == 0;
        if ( m_converged )
        {
          a = b = c;
          fa = fb = 0;
          return c;
        }

        check( c, fc );

        if ( fa * fc < 0 )
        {
          b  = c;
          fb = fc;
        }  // --> [a,c]
        else
        {
          a  = c;
          fa = fc;
        }  // --> [c,b]

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };

        m_converged = ( b - a ) < m_tolerance_x || abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) return abs_fb < abs_fa ? b : a;
      }
    }

    /*
    //  _     _               _   _
    // | |__ (_)___  ___  ___| |_(_) ___  _ __
    // | '_ \| / __|/ _ \/ __| __| |/ _ \| '_ \
    // | |_) | \__ \  __/ (__| |_| | (_) | | | |
    // |_.__/|_|___/\___|\___|\__|_|\___/|_| |_|
    */

    auto bisection = [this, check]() -> void
    {
      Real & a{ m_a };
      Real & b{ m_b };
      Real & fa{ m_fa };
      Real & fb{ m_fb };
      while ( ++m_iteration_count < m_max_iteration )
      {
        Real ba{ b - a };
        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = ba < m_tolerance_x || abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) break;

        Real c{ a + ba / 2 };
        Real fc{ evaluate( c ) };
        check( c, fc );

        m_converged = fc == 0;
        if ( m_converged )
        {
          a = b = c;
          fa = fb = 0;
          break;
        }

        if ( ( fb > 0 ) == ( fc > 0 ) )
        {
          b  = c;
          fb = fc;
        }
        else
        {
          a  = c;
          fa = fc;
        }
      }
    };

    /*
    //   _ _ _ _             _
    //  (_) | (_)_ __   ___ (_)___
    //  | | | | | '_ \ / _ \| / __|
    //  | | | | | | | | (_) | \__ \
    //  |_|_|_|_|_| |_|\___/|_|___/
    */

    auto illinois = [this, check]() -> void
    {
      // illinois
      Real   dd{ 2 };
      Real & a{ m_a };
      Real & b{ m_b };
      Real & fa{ m_fa };
      Real & fb{ m_fb };
      while ( ++m_iteration_count < m_max_iteration )
      {
        Real ba{ b - a };
        Real tol{ m_tolerance_x / ( 2 * abs( ba ) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) break;

        Real c{ abs_fa < abs_fb ? a + max( tol, fa / ( fa - fb ) ) * ba : b - max( tol, fb / ( fb - fa ) ) * ba };
        Real fc{ evaluate( c ) };

        check( c, fc );

        m_converged = fc == 0;
        if ( m_converged )
        {
          a = b = c;
          fa = fb = 0;
          break;
        }

        if ( ( fa > 0 ) == ( fc > 0 ) )
        {
          a  = b;
          fa = fb;
          dd = 2;
        }
        else
        {
          fa /= dd;
          dd *= 2;
        }

        b  = c;
        fb = fc;
      }
    };

    /*
    //    ____ _                     _                        _   _
    //   / ___| |__   __ _ _ __   __| |_ __ _   _ _ __   __ _| |_| | __ _
    //  | |   | '_ \ / _` | '_ \ / _` | '__| | | | '_ \ / _` | __| |/ _` |
    //  | |___| | | | (_| | | | | (_| | |  | |_| | |_) | (_| | |_| | (_| |
    //   \____|_| |_|\__,_|_| |_|\__,_|_|   \__,_| .__/ \__,_|\__|_|\__,_|
    //                                           |_|
    */
    auto Chandrupatla = [this, check]() -> void
    {
      /*
        Tirupathi Chandrupatla,
        A new hybrid quadratic/bisection algorithm for finding the zero of
        a nonlinear function without using derivatives,
        Advances in Engineering Software,
        Volume 28, Number 3, pages 145-149, 1997.
      */

      Real   t{ 0.5 };
      Real & a{ m_a };
      Real & b{ m_b };
      Real & fa{ m_fa };
      Real & fb{ m_fb };

      while ( ++m_iteration_count < m_max_iteration )
      {
        Real dir{ b - a };
        Real tol{ m_tolerance_x / ( 2 * abs( dir ) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) break;

        if ( t < tol )
          t = tol;
        else if ( t > 1 - tol )
          t = 1 - tol;

        Real c{ a + t * dir };
        Real fc{ evaluate( c ) };

        check( c, fc );

        m_converged = fc == 0;
        if ( m_converged )
        {
          a = b = c;
          fa = fb = 0;
          break;
        }

        // Arrange 2-1-3: 2-1 Interval, 1 Middle, 3 Discarded point.
        Real d, fd;
        if ( ( 0 < fc ) == ( 0 < fa ) )
        {
          //  a => d     ---> [d,a,b] = [a,c,b]
          d  = a;
          fd = fa;
        }
        else
        {
          //  b => d
          //  a => b    ----> [b,a,d] = [a,c,b]
          d  = b;
          fd = fb;
          b  = a;
          fb = fa;
        }
        //  c => a
        a  = c;
        fa = fc;

        // If inverse quadratic interpolation holds, use it.
        Real ba{ b - a };
        Real fba{ fb - fa };
        Real bd{ b - d };
        Real fbd{ fb - fd };

        Real xi{ ba / bd };
        Real ph{ fba / fbd };
        Real fl{ 1 - sqrt( 1 - xi ) };
        Real fh{ sqrt( xi ) };

        if ( fl < ph && ph < fh )
        {
          Real da{ d - a };
          Real fda{ fd - fa };

          t = ( fa / fba ) * ( fd / fbd ) - ( fa / fda ) * ( fb / fbd ) * ( da / ba );
        }
        else
        {
          t = 0.5;
        }
      }
    };

    /*
    //   ____                 _
    //  | __ ) _ __ ___ _ __ | |_
    //  |  _ \| '__/ _ \ '_ \| __|
    //  | |_) | | |  __/ | | | |_
    //  |____/|_|  \___|_| |_|\__|
    */

    auto Brent = [this, check]() -> void
    {
      /*
      // Richard Brent,
      // Algorithms for Minimization without Derivatives,
      // Dover, 2002, ISBN: 0-486-41998-3
      */

      Real & a{ m_a };
      Real & fa{ m_fa };
      Real & b{ m_b };
      Real & fb{ m_fb };

      Real c{ a };
      Real fc{ fa };
      Real e{ b - a };
      Real d{ e };

      while ( ++m_iteration_count < m_max_iteration )
      {
        if ( abs( fc ) < abs( fb ) )
        {
          a  = b;
          b  = c;
          c  = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }

        Real tol{ 2 * Utils::machine_eps<Real>() * abs( b ) + m_tolerance_x };
        Real m{ ( c - b ) / 2 };

        m_converged = abs( m ) <= tol || fb == 0;
        if ( m_converged )
        {
          m_a  = m_b;
          m_fa = 0;
          break;
        }

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) break;

        if ( abs( e ) < tol || abs( fa ) <= abs( fb ) ) { d = e = m; }
        else
        {
          Real s{ fb / fa };
          Real p, q;
          if ( a == c )
          {
            p = 2 * m * s;
            q = 1 - s;
          }
          else
          {
            q = fa / fc;
            Real r{ fb / fc };
            p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
            q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
          }
          if ( p > 0 )
            q = -q;
          else
            p = -p;

          s = e;
          e = d;

          if ( 2 * p < 3 * m * q - abs( tol * q ) && 2 * p < abs( s * q ) ) { d = p / q; }
          else
          {
            e = m;
            d = e;
          }
        }
        a  = b;
        fa = fb;

        if ( tol < abs( d ) )
          b += d;
        else if ( m > 0 )
          b += tol;
        else
          b -= tol;

        fb = evaluate( b );
        check( b, fb );

        if ( ( 0 < fb && 0 < fc ) || ( fb <= 0 && fc <= 0 ) )
        {
          c  = a;
          fc = fa;
          e  = b - a;
          d  = e;
        }
      }
    };

    /*
    //   ____  _     _     _
    //  |  _ \(_) __| | __| | ___ _ __
    //  | |_) | |/ _` |/ _` |/ _ \ '__|
    //  |  _ <| | (_| | (_| |  __/ |
    //  |_| \_\_|\__,_|\__,_|\___|_|
    */
    auto Ridder = [this, check]() -> void
    {
      /*
         C.F.J.Ridders, “A New Algorithm for Computing a Single Root of a Real
         Continuous Function.” IEEE Trans. Circuits Systems 26, 979-980, 1979.
      */

      Real & a{ m_a };
      Real & fa{ m_fa };
      Real & b{ m_b };
      Real & fb{ m_fb };

      while ( ++m_iteration_count < m_max_iteration )
      {
        Real d{ ( b - a ) / 2 };
        m_converged = 2 * abs( d ) <= m_tolerance_x;
        if ( m_converged ) break;

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) break;

        // Compute the improved root x from Ridder's formula
        Real c{ a + d };
        Real fc{ evaluate( c ) };

        check( c, fc );
        m_converged = fc == 0;
        if ( m_converged )
        {
          a = b = c;
          fa = fb = 0;
          break;
        }

        Real A{ fc / fa };
        Real B{ fb / fa };
        Real dx{ d * A / sqrt( A * A - B ) };
        Real x{ c + dx };
        Real fx{ evaluate( x ) };

        m_converged = fx == 0;
        if ( m_converged )
        {
          a = b = x;
          fa = fb = 0;
          break;
        }

        if ( c > x )
        {
          swap( x, c );
          swap( fx, fc );
        }

        // Re-bracket the root as tightly as possible
        if ( ( fc > 0 ) != ( fx > 0 ) )
        {
          a  = c;
          fa = fc;
          b  = x;
          fb = fx;
        }
        else
        {
          if ( ( fa > 0 ) != ( fc > 0 ) )
          {
            b  = c;
            fb = fc;
          }
          else
          {
            a  = x;
            fa = fx;
          }
        }
      }
    };

    /*
    //   __  __           _ _  __ _          _
    //  |  \/  | ___   __| (_)/ _(_) ___  __| |
    //  | |\/| |/ _ \ / _` | | |_| |/ _ \/ _` |
    //  | |  | | (_) | (_| | |  _| |  __/ (_| |
    //  |_|  |_|\___/ \__,_|_|_| |_|\___|\__,_|
    //
    //      _              _                                 ____  _ _
    //     / \   _ __   __| | ___ _ __ ___  ___  _ __       | __ )(_) ___  _ __|
    | __
    //    / _ \ | '_ \ / _` |/ _ \ '__/ __|/ _ \| '_ \ _____|  _ \| |/ _ \| '__|
    |/ /
    //   / ___ \| | | | (_| |  __/ |  \__ \ (_) | | | |_____| |_) | | (_) | |  |
    <
    //  /_/   \_\_| |_|\__,_|\___|_|  |___/\___/|_| |_|     |____// |\___/|_|
    |_|\_\
    //                                                          |__/
    */

    auto modified_AB = [this, check]() -> void
    {
      /*
      N. Ganchovski, A. Traykov
      Modified Anderson-Bjork’s method for solving nonlinear equations in
      structural mechanics IOP Conf. Series: Materials Science and Engineering
      1276 (2023) 012010 doi:10.1088/1757-899X/1276/1/012010
      */

      Integer side{ 0 };
      // Integer N{ Integer(floor(1-log2( m_tolerance )/2)) };

      Real & x1{ m_a };
      Real & f1{ m_fa };
      Real & x2{ m_b };
      Real & f2{ m_fb };

      bool bisection{ true };

      while ( ++m_iteration_count < m_max_iteration )
      {
        Real dir{ x2 - x1 };
        Real tol{ m_tolerance_x / ( 2 * abs( dir ) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;

        Real abs_f1{ abs( f1 ) };
        Real abs_f2{ abs( f2 ) };
        m_converged = abs_f1 < m_tolerance_f || abs_f2 < m_tolerance_f;
        if ( m_converged ) break;

        Real x3, f3;

        if ( bisection )
        {
          x3 = x1 + dir / 2;    // Midpoint abscissa
          f3 = evaluate( x3 );  // and function value

          check( x3, f3 );
          m_converged = f3 == 0;
          if ( m_converged )
          {
            x1 = x2 = x3;
            f1 = f2 = 0;
            break;
          }

          Real fm{ ( f1 + f2 ) / 2 };  // Ordinate of chord at midpoint
          bisection = abs( fm - f3 ) >= ( abs( fm ) + abs( f3 ) ) / 4;
        }
        else
        {
          Real t{ f1 / ( f1 - f2 ) };
          if ( t < tol )
            t = tol;
          else if ( t > 1 - tol )
            t = 1 - tol;

          x3 = x1 + t * dir;  // False-position step
          f3 = evaluate( x3 );

          check( x3, f3 );
          m_converged = f3 == 0;
          if ( m_converged )
          {
            x1 = x2 = x3;
            f1 = f2 = 0;
            break;
          }
        }

        switch ( side )
        {
          case 1:  // Apply Anderson-Bjork modification for side 1
          {
            Real m{ 1 - f3 / f1 };
            f2 *= m > Utils::machine_eps<Real>() ? m : static_cast<Real>( 0.5 );
          }
          break;
          case 2:  // Apply Anderson-Bjork modification for side 2
          {
            Real m{ 1 - f3 / f2 };
            f1 *= m > Utils::machine_eps<Real>() ? m : static_cast<Real>( 0.5 );
          }
          break;
        }

        if ( ( f1 > 0 ) == ( f3 > 0 ) )
        {                              // If the left interval does not change sign
          if ( !bisection ) side = 1;  // Store the side that move
          x1 = x3;
          f1 = f3;  // Move the left end
        }
        else
        {                              // If the right interval does not change sign
          if ( !bisection ) side = 2;  // Store the side that move
          x2 = x3;
          f2 = f3;  // Move the right end
        }

        // if ( (m_iteration_count % N) == 0 ) {
        //   bisection = true;
        //   side = 0;
        // }
      }
    };

    /*
    //      _    _           _____ _  _    ___
    //     / \  | | __ _  __|___  | || |  ( _ )
    //    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
    //   / ___ \| | (_| | (_) / / |__   _| (_) |
    //  /_/   \_\_|\__, |\___/_/     |_|  \___/
    //             |___/
    */
    auto algo748 = [this, check]() -> void
    {
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

      auto all_different = []( Real a, Real b, Real c, Real d ) -> bool
      { return a != b && a != c && a != d && b != c && b != d && c != d; };

      Real & a{ m_a };
      Real & fa{ m_fa };
      Real & b{ m_b };
      Real & fb{ m_fb };

      Real m_interval_shink{ Real( 0.025 ) };

      Real c{ NaN<Real>() };
      Real fc{ NaN<Real>() };  // Dumb values
      Real d{ NaN<Real>() };
      Real fd{ NaN<Real>() };  // Dumb values
      Real e{ NaN<Real>() };
      Real fe{ NaN<Real>() };  // Dumb values

      auto pzero = [this, &d, &fd, &e, &fe]() -> Real
      {
        // Uses cubic inverse interpolation of f(x) at a, b, d, and e to
        // get an approximate root of f(x). Rewritten using divided difference.

        Real & a{ m_a };
        Real & fa{ m_fa };
        Real & b{ m_b };
        Real & fb{ m_fb };

        Real D1{ b - a };
        Real D2{ d - a };
        Real D3{ e - a };

        Real DD0{ D1 / ( fb - fa ) };
        Real DD1{ ( D1 - D2 ) / ( fb - fd ) };
        Real DD2{ ( D2 - D3 ) / ( fd - fe ) };

        Real DDD0{ ( DD0 - DD1 ) / ( fa - fd ) };
        Real DDD1{ ( DD1 - DD2 ) / ( fb - fe ) };

        Real DDDD0{ ( DDD0 - DDD1 ) / ( fa - fe ) };

        Real c{ a - fa * ( DD0 - fb * ( DDD0 - fd * DDDD0 ) ) };

        Real tol{ Real( 0.7 ) * m_tolerance_x };
        if ( c <= a + tol || c >= b - tol ) c = ( a + b ) / 2;

        UTILS_ASSERT( is_finite( c ),
                      "AlgoBracket[748]::pzero(), compute NaN or Inf at\n"
                      "a={} f(a)={}\n"
                      "b={} f(b)={}\n"
                      "c={}\n"
                      "d={} f(d)={}\n"
                      "e={} f(e)={}\n",
                      a, fa, b, fb, c, d, fd, e, fe );

        return c;
      };

      auto newton_quadratic = [this, &d, &fd]( Integer const niter, Real & c ) -> bool
      {
        // Uses `niter` newton steps to approximate the zero in (a,b) of the
        // quadratic polynomial interpolating f(x) at a, b, and d.
        // Safeguard is used to avoid overflow.

        Real & a{ m_a };
        Real & fa{ m_fa };
        Real & b{ m_b };
        Real & fb{ m_fb };

        UTILS_ASSERT( a < b && a != d && b != d,
                      "AlgoBracket[748]::newton_quadratic() bad data\n"
                      "a={} f(a)={}\n"
                      "b={} f(b)={}\n"
                      "d={} f(d)={}\n",
                      a, fa, b, fb, d, fd );

        Real A0{ fa };
        Real A1{ ( fb - fa ) / ( b - a ) };
        Real A2{ ( ( fd - fb ) / ( d - b ) - A1 ) / ( d - a ) };

        UTILS_ASSERT( is_finite( A0 ) && is_finite( A1 ) && is_finite( A2 ),
                      "AlgoBracket[748]::newton_quadratic(), compute NaN or Inf at\n"
                      "a={} f(a)={}\n"
                      "b={} f(b)={}\n"
                      "d={} f(d)={}\n"
                      "A0={}\n"
                      "A1={}\n"
                      "A2={}\n",
                      a, fa, b, fb, d, fd, A0, A1, A2 );

        // Safeguard to avoid overflow.
        if ( A2 == 0 ) { c = a - A0 / A1; }
        else
        {
          // Determine the starting point of newton steps.
          c = A2 * fa > 0 ? a : b;

          // Start the safeguarded newton steps.
          bool ok{ true };
          for ( Integer i{ 0 }; i < niter && ok; ++i )
          {
            Real PC  = A0 + ( A1 + A2 * ( c - b ) ) * ( c - a );
            Real PDC = A1 + A2 * ( ( 2 * c ) - ( a + b ) );
            ok       = PDC != 0;
            if ( ok ) c -= PC / PDC;
          }
        }
        return is_finite( c ) && c > a && c < b;
      };

      auto bracketing = [this, check]( Real & c, Real & fc, Real & d, Real & fd ) -> bool
      {
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
          Real tol{ Real( 0.7 ) * m_tolerance_x };
          Real hba{ ( m_b - m_a ) / 2 };
          if ( hba <= tol )
            c = m_a + hba;
          else if ( c <= m_a + tol )
            c = m_a + tol;
          else if ( c >= m_b - tol )
            c = m_b - tol;
        }

        UTILS_ASSERT( is_finite( c ),
                      "AlgoBracket[748]::bracketing(), unexpected\n"
                      "c={} at [a,b] = [{},{}]\n",
                      c, m_a, m_b );

        fc = this->evaluate( c );
        check( c, fc );

        // If f(c)=0, then set a=b=c and return.
        // This will terminate the procedure.

        if ( fc == 0 )
        {
          m_a = m_b = c;
          m_fa = m_fb = d = fd = 0;
          return true;
        }
        // If f(c) is not zero, then determine the new enclosing interval.
        if ( m_fa * fc < 0 )
        {
          // D <-- B <-- C
          d    = m_b;
          fd   = m_fb;
          m_b  = c;
          m_fb = fc;
        }
        else
        {
          // D <-- A <-- C
          d    = m_a;
          fd   = m_fa;
          m_a  = c;
          m_fa = fc;
        }
        // update the termination criterion according to the new enclosing
        // interval.
        // if ( abs(m_fb) <= abs(m_fa) ) set_tolerance(m_b);
        // else                          set_tolerance(m_a);
        return false;
      };

      {
        using std::abs;
        using std::min;
        using std::max;
        Real ba{ b - a };
        Real R{ ba / ( fb - fa ) };
        c = abs( fb ) < abs( fa ) ? b + fb * R : a - fa * R;
        // impedisce m_c troppo vicino ad m_a o m_b
        Real delta{ m_interval_shink * ba };
        c = max( min( c, b - delta ), a + delta );
      }
      //
      // Call "bracketing" to get a shrinked enclosing interval as
      // well as to update the termination criterion.
      // Stop the procedure if the criterion is satisfied or the
      // exact solution is obtained.
      //
      m_converged = bracketing( c, fc, d, fd );
      if ( m_converged )
      {
        b  = a;
        fb = fa;
        return;
      }
      // Iteration starts.
      // The enclosing interval before executing the iteration is recorded as
      // [a0, b0].
      m_converged = false;

      // ITERATION STARTS. THE ENCLOSING INTERVAL BEFORE EXECUTING THE
      // ITERATION IS RECORDED AS [A0, B0].
      while ( !m_converged )
      {
        ++m_iteration_count;
        Real BA0{ b - a };

        // Calculates the termination criterion.
        // Stops the procedure if the criterion is satisfied.

        m_converged = BA0 <= m_tolerance_x;
        if ( m_converged ) return;

        Real abs_fa{ abs( fa ) };
        Real abs_fb{ abs( fb ) };
        m_converged = abs_fa < m_tolerance_f || abs_fb < m_tolerance_f;
        if ( m_converged ) return;

        //
        // Starting with the second iteration, in the first two steps, either
        // quadratic interpolation is used by calling the subroutine
        // "newtonquadratic" or the cubic inverse interpolation is used by
        // calling the subroutine "pzero". in the following, if "prof" is not
        // equal to 0, then the four function values "fa", "fb", "fd", and "fe"
        // are distinct, and hence "pzero" will be called.
        //

        bool do_newton_quadratic = false;
        if ( !is_finite( fe ) ) { do_newton_quadratic = true; }
        else if ( !all_different( fa, fb, fd, fe ) ) { do_newton_quadratic = true; }
        else
        {
          c                   = pzero();
          do_newton_quadratic = ( c - a ) * ( c - b ) >= 0;
        }
        if ( do_newton_quadratic )
        {
          if ( !newton_quadratic( 2, c ) ) c = a + ( b - a ) / 2;
        }

        e  = d;
        fe = fd;

        //
        // Call subroutine "bracketing" to get a shrinked enclosing interval as
        // well as to update the termination criterion. stop the procedure
        // if the criterion is satisfied or the exact solution is obtained.
        //
        m_converged = bracketing( c, fc, d, fd ) || ( b - a ) <= m_tolerance_x;
        if ( m_converged ) return;

        if ( !all_different( fa, fb, fd, fe ) ) { do_newton_quadratic = true; }
        else
        {
          c                   = pzero();
          do_newton_quadratic = ( c - a ) * ( c - b ) >= 0;
        }
        if ( do_newton_quadratic )
        {
          if ( !newton_quadratic( 3, c ) ) c = a + ( b - a ) / 2;
        }

        //
        // Call subroutine "bracketing" to get a shrinked enclosing interval as
        // well as to update the termination criterion. stop the procedure
        // if the criterion is satisfied or the exact solution is obtained.
        //

        {
          m_converged = bracketing( c, fc, d, fd ) || ( b - a ) <= m_tolerance_x;
          if ( m_converged ) return;

          e  = d;
          fe = fd;
          // Takes the double-size secant step.
          Real u, fu;
          if ( abs( fa ) < abs( fb ) )
          {
            u  = a;
            fu = fa;
          }
          else
          {
            u  = b;
            fu = fb;
          }
          Real hba{ ( b - a ) / 2 };
          c = u - 4 * ( fu / ( fb - fa ) ) * hba;
          if ( abs( c - u ) > hba ) c = a + hba;
        }

        //
        // Call subroutine "bracketing" to get a shrinked enclosing interval as
        // well as to update the termination criterion. stop the procedure
        // if the criterion is satisfied or the exact solution is obtained.
        //
        m_converged = bracketing( c, fc, d, fd ) || ( b - a ) <= m_tolerance_x;
        if ( m_converged ) return;
        //
        // Determines whether an additional bisection step is needed.
        // Takes it if necessary.
        //
        Real mu{ Real( 0.5 ) };
        if ( ( b - a ) < mu * BA0 ) continue;

        e  = d;
        fe = fd;
        //
        // Call subroutine "bracketing" to get a shrinked enclosing interval as
        // well as to update the termination criterion. stop the procedure
        // if the criterion is satisfied or the exact solution is obtained.
        //
        {
          Real ba{ b - a };
          c           = a + ba / 2;
          m_converged = bracketing( c, fc, d, fd ) || ba <= m_tolerance_x;
        }
      }
    };

    switch ( m_select )
    {
      case Method::BISECTION:
        bisection();
        break;
      case Method::ILLINOIS:
        illinois();
        break;
      case Method::CHANDRUPATLA:
        Chandrupatla();
        break;
      case Method::BRENT:
        Brent();
        break;
      case Method::RIDDER:
        Ridder();
        break;
      case Method::MODIFIED_AB:
        modified_AB();
        break;
      case Method::ALGO748:
        algo748();
        break;
    }

    if ( m_a > m_b )
    {
      swap( m_a, m_b );
      swap( m_fa, m_fb );
    }

    return abs( m_fa ) < abs( m_fb ) ? m_a : m_b;
  }

  template class AlgoBracket<float>;
  template class AlgoBracket<double>;
}  // namespace Utils
