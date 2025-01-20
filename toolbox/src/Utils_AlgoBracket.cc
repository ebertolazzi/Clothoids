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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Utils_AlgoBracket.hh"
#include "Utils_fmt.hh"

#include <vector>
#include <cmath>
#include <algorithm>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils {

  using std::max;
  using std::min;
  using std::abs;
  using std::signbit;
  using std::ceil;
  using std::log2;
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
  AlgoBracket<Real>::set_max_iterations( Integer mit ) {
    UTILS_ASSERT(
      mit > 0,
      "AlgoBracket::set_max_iterations({}) argument must be >0\n", mit
    );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  AlgoBracket<Real>::set_max_fun_evaluation( Integer mfev ) {
    UTILS_ASSERT(
      mfev > 0,
      "AlgoBracket::set_max_fun_evaluation({}) argument must be >0\n", mfev
    );
    m_max_fun_evaluation = mfev;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename Real>
  Real
  AlgoBracket<Real>::p_zero2( Real d, Real fd ) const {
    Real D1    { m_fb - m_fa };
    Real ba    { m_b - m_a };
    Real theta { (d - m_a)/ba };
    Real D2    { fd - m_fa };
    Real A     { D2/theta-D1 };
    Real B     { D1*theta-D2/theta };
    Real C     { m_fa*(theta-1) };
    Real M     { max(max(abs(A),abs(B)),abs(C)) };
    if ( C < 0 ) M = -M;
    A /= M; B /= M; C /= M;
    Real D{ B*B-4*A*C };
    if ( D >= 0 ) {
      D = sqrt(D);
      return B < 0 ? 2*C/(D-B) : (D+B)/(2*A);
    } else {
      return -m_fa / D1;
    }
  }

  template <typename Real>
  Real
  AlgoBracket<Real>::invp_zero2( Real d, Real fd ) const {
    //
    // Uses quadratic inverse interpolation of f(x) at a, b, and d to
    // get an approximate root of f(x).
    // Rewritten using divided difference.
    
    Real x0{ m_fa };
    Real x1{ m_fb };
    Real x2{ fd   };

    Real D0{ 0 };
    Real D1{ 1 };
    Real D2{ (d-m_a)/(m_b-m_a) };

    Real D01{ (D0-D1)/(x0-x1) };
    Real D12{ (D1-D2)/(x1-x2) };

    Real D012{ (D01-D12)/(x0-x2) };
    
    Real O1{ 0-x0 };
    Real O2{ (0-x1)*O1 };

    Real P0{ D0             };
    Real P1{ P0 + D01  * O1 };
    Real P2{ P1 + D012 * O2 };

    UTILS_ASSERT(
      is_finite(P2),
      "AlgoBracket<Real>::pzero(), compute NaN or Inf at\n"
      "a={} f(a)={}\n"
      "b={} f(b)={}\n"
      "d={} f(d)={}\n",
      m_a, m_fa,
      m_b, m_fb,
      d,   fd
    );

    // CALCULATE THE OUTPUT C.
    return P2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoBracket<Real>::invp_zero3( Real d, Real fd, Real e, Real fe ) const {
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
    Real D2{ (d-m_a) / (m_b-m_a) };
    Real D3{ (e-m_a) / (m_b-m_a) };

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
      "AlgoBracket<Real>::pzero(), compute NaN or Inf at\n"
      "a={} f(a)={}\n"
      "b={} f(b)={}\n"
      "d={} f(d)={}\n"
      "e={} f(e)={}\n",
      m_a, m_fa,
      m_b, m_fb,
      d,   fd,
      e,   fe
    );

    // CALCULATE THE OUTPUT C.
    return P3;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoBracket<Real>::eval( Real a, Real b ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    m_a = a; m_fa = this->evaluate(m_a);
    m_b = b; m_fb = this->evaluate(m_b);

    // check if solution can exists
    if ( m_fa*m_fb > 0 ) return m_a;
    else                 return eval();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  AlgoBracket<Real>::eval( Real a, Real b, Real amin, Real bmax ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    m_a = a; m_fa = this->evaluate(m_a);
    m_b = b; m_fb = this->evaluate(m_b);

    // try to enlarge interval
    while ( m_fa*m_fb > 0 ) {
      if ( is_finite(m_fa) && m_a > amin ) {
        m_a -= m_b - m_a;
        m_fa = this->evaluate(m_a);
      } else if ( is_finite(m_fb) && m_b < bmax ) {
        m_b += m_b - m_a;
        m_fb = this->evaluate(m_b);
      } else {
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

    2. Starting with the second iteration, three steps are taken in each iteration.

       a. First two steps are either quadratic interpolation or cubic inverse interpolation.

       b. The third step is a double-size secant step.

    If the diameter of the enclosing interval obtained after
    those three steps is larger than \f$ (b-a)/2 \f$,
    then an additional bisection step will be taken.


  */

  template <typename Real>
  Real
  AlgoBracket<Real>::eval() {

    Real & a  { m_a  };
    Real & b  { m_b  };
    Real & fa { m_fa };
    Real & fb { m_fb };

    UTILS_ASSERT(
      !is_NaN( fa ) && !is_NaN( fb ),
      "AlgoBracket::eval() bad initial interval\n"
      "a = {}, fa = {}\n"
      "b = {}, fb = {}\n",
      a, fa,
      b, fb
    );

    // check for trivial solution
    m_converged = fa == 0; if ( m_converged ) return a;
    m_converged = fb == 0; if ( m_converged ) return b;
    m_tolerance = 10*Utils::machine_eps<Real>()*( 1 + 2*max(abs(a),abs(b)) );

    auto check = [this] ( Real c, Real fc ) -> void {
      UTILS_ASSERT(
        is_finite(fc),
        "AlgoBracket::eval() found Inf or NaN\n"
        "a={} fa={}\n"
        "b={} fb={}\n"
        "c={} fc={}\n",
        m_a, m_fa,
        m_b, m_fb,
        c,   fc
      );
    };

    //
    // While f(left) or f(right) are infinite perform bisection
    //
    while ( !( is_finite(fa) && is_finite(fb) ) ) {
    
      UTILS_ASSERT(
        ++m_iteration_count <= m_max_iteration,
        "AlgoBracket::eval() too many iteration\n"
        "a={} fa={}\n"
        "b={} fb={}\n",
        a, fa,
        b, fb
      );
 
      Real c  { (a+b)/2           };
      Real fc { this->evaluate(c) };

      UTILS_ASSERT( !is_NaN( fc ), "AlgoBracket::eval()\nc = {}, fc = {}\n", c, fc );

      m_converged = fc == 0;
      if ( m_converged ) { a = b = c; fa = fb = 0; return c; }
      
      check( c, fc );

      if ( m_fa*fc < 0 ) { b = c; fb = fc; } // --> [a,c]
      else               { a = c; fa = fc; } // --> [c,b]

      m_converged = (b-a) <= m_tolerance;
      if ( m_converged ) return abs(fb) < abs(fa) ? b : a;
    }
    
    /*
    //  _     _               _   _
    // | |__ (_)___  ___  ___| |_(_) ___  _ __
    // | '_ \| / __|/ _ \/ __| __| |/ _ \| '_ \
    // | |_) | \__ \  __/ (__| |_| | (_) | | | |
    // |_.__/|_|___/\___|\___|\__|_|\___/|_| |_|
    */

    auto bisection = [this,check] () -> void {
      Real & a  { m_a  };
      Real & b  { m_b  };
      Real & fa { m_fa };
      Real & fb { m_fb };
      while ( ++m_iteration_count < m_max_iteration ) {

        Real ba{ b - a };
        m_converged = ba <= m_tolerance;
        if ( m_converged ) break;
        
        Real c  { a + ba/2      };
        Real fc { evaluate( c ) };
        check( c, fc );

        m_converged = fc == 0;
        if ( m_converged ) { a = b = c; fa = fb = 0; break; }
          
        if ( (fb > 0) == (fc > 0) ) { b = c; fb = fc; }
        else                        { a = c; fa = fc; }
      }
    };

    /*
    //   _ _ _ _             _
    //  (_) | (_)_ __   ___ (_)___
    //  | | | | | '_ \ / _ \| / __|
    //  | | | | | | | | (_) | \__ \
    //  |_|_|_|_|_| |_|\___/|_|___/
    */

    auto illinois = [this,check] () -> void {
      // illinois
      Real dd{2};
      Real & a  { m_a  };
      Real & b  { m_b  };
      Real & fa { m_fa };
      Real & fb { m_fb };
      while ( ++m_iteration_count < m_max_iteration ) {
      
        Real ba{ b - a };
        Real tol{ m_tolerance / ( 2 * abs(ba) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;

        Real c{ abs(fa) < abs(fb) ?
                a + max( tol, fa/(fa-fb) ) * ba :
                b - max( tol, fb/(fb-fa) ) * ba };
        Real fc{ evaluate( c ) };
      
        check( c, fc );
      
        m_converged = fc == 0;
        if ( m_converged ) { a = b = c; fa = fb = 0; break; }

        if ( (fa > 0) == (fc > 0) ) { a = b; fa = fb; dd = 2; }
        else                        { fa /= dd; dd *= 2; }

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
    auto Chandrupatla = [this,check] () -> void {

      /*
        Tirupathi Chandrupatla,
        A new hybrid quadratic/bisection algorithm for finding the zero of
        a nonlinear function without using derivatives,
        Advances in Engineering Software,
        Volume 28, Number 3, pages 145-149, 1997.
      */

      Real   t  { 0.5  };
      Real & a  { m_a  };
      Real & b  { m_b  };
      Real & fa { m_fa };
      Real & fb { m_fb };
      
      while ( ++m_iteration_count < m_max_iteration ) {

        Real dir{ b-a };
        Real tol{ m_tolerance / ( 2 * abs( dir ) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;
        
        if      ( t < tol   ) t = tol;
        else if ( t > 1-tol ) t = 1-tol;

        Real c  { a + t * dir };
        Real fc { evaluate( c ) };

        check( c, fc );

        m_converged = fc == 0;
        if ( m_converged ) { a = b = c; fa = fb = 0; break; }

        // Arrange 2-1-3: 2-1 Interval, 1 Middle, 3 Discarded point.
        Real d, fd;
        if ( ( 0 < fc ) == ( 0 < fa ) ) {
          //  a => d     ---> [d,a,b] = [a,c,b]
          d = a; fd = fa;
        } else {
          //  b => d
          //  a => b    ----> [b,a,d] = [a,c,b]
          d = b; fd = fb;
          b = a; fb = fa;
        }
        //  c => a
        a = c; fa = fc;

        // If inverse quadratic interpolation holds, use it.
        Real ba  { b  - a  };
        Real fba { fb - fa };
        Real bd  { b  - d  };
        Real fbd { fb - fd };

        Real xi { ba / bd };
        Real ph { fba / fbd };
        Real fl { 1 - sqrt( 1 - xi ) };
        Real fh { sqrt( xi ) };

        if ( fl < ph && ph < fh ) {
          
          Real da  { d  - a  };
          Real fda { fd - fa };
       
          t = (fa/fba) * (fd/fbd) - (fa/fda) * (fb/fbd) * (da/ba);

        } else {
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
    
    auto Brent = [this,check] () -> void {
      /*
      // Richard Brent,
      // Algorithms for Minimization without Derivatives,
      // Dover, 2002, ISBN: 0-486-41998-3
      */

      Real & a  { m_a };
      Real & fa { m_fa };
      Real & b  { m_b };
      Real & fb { m_fb };

      Real c  { a  };
      Real fc { fa };
      Real e  { b - a };
      Real d  { e };

      while ( ++m_iteration_count < m_max_iteration ) {

        if ( abs( fc ) < abs( fb ) ) {
          a  = b;  b  = c;  c  = a;
          fa = fb; fb = fc; fc = fa;
        }

        Real tol { 2 * Utils::machine_eps<Real>() * abs( b ) + m_tolerance };
        Real m   { ( c - b )/2 };
        
        m_converged = abs( m ) <= tol || fb == 0;
        if ( m_converged ) { m_a = m_b; m_fa = 0; break; }

        if ( abs( e ) < tol || abs( fa ) <= abs( fb ) ) {
          d = e = m;
        } else {
          Real s{ fb / fa };
          Real p, q;
          if ( a == c ) {
            p = 2 * m * s;
            q = 1 - s;
          } else {
            q = fa / fc;
            Real r{ fb / fc };
            p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
            q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
          }
          if ( p > 0 ) q = -q;
          else         p = -p;
   
          s = e;
          e = d;

          if ( 2 * p < 3 * m * q - abs( tol * q ) && 2 * p < abs( s * q ) ) {
            d = p / q;
          } else {
            e = m;
            d = e;
          }
        }
        a  = b;
        fa = fb;

        if      ( tol < abs(d) ) b += d;
        else if ( m   > 0      ) b += tol;
        else                     b -= tol;

        fb = evaluate( b );
        check( b, fb );

        if ( ( 0 < fb && 0 < fc ) || ( fb <= 0 && fc <= 0 ) ) {
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
    auto Ridder = [this,check] () -> void {
      /*
         C.F.J.Ridders, “A New Algorithm for Computing a Single Root of a Real Continuous Function.”
         IEEE Trans. Circuits Systems 26, 979-980, 1979.
      */

      Real & a  { m_a  };
      Real & fa { m_fa };
      Real & b  { m_b  };
      Real & fb { m_fb };

      while ( ++m_iteration_count < m_max_iteration ) {
      
        Real d{ (b - a)/2 };
        m_converged = 2*abs(d) <= m_tolerance;
        if ( m_converged ) break;

        // Compute the improved root x from Ridder's formula
        Real c  { a + d };
        Real fc { evaluate(c) };

        check( c, fc );
        m_converged = fc == 0;
        if ( m_converged ) { a = b = c; fa = fb = 0; break; }
        
        Real A  { fc / fa };
        Real B  { fb / fa };
        Real dx { d*A/sqrt( A*A - B ) };
        Real x  { c + dx };
        Real fx { evaluate(x) };

        m_converged = fx == 0;
        if ( m_converged ) { a = b = x; fa = fb = 0; break; }
        
        if ( c > x ) { swap( x, c ); swap( fx, fc ); }

        // Re-bracket the root as tightly as possible
        if ( ( fc > 0 ) != ( fx > 0 ) ) {
          a = c; fa = fc;
          b = x; fb = fx;
        } else {
          if ( ( fa > 0 ) != ( fc > 0 ) ) { b = c; fb = fc; }
          else                            { a = x; fa = fx; }
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
    //      _              _                                 ____  _            _
    //     / \   _ __   __| | ___ _ __ ___  ___  _ __       | __ )(_) ___  _ __| | __
    //    / _ \ | '_ \ / _` |/ _ \ '__/ __|/ _ \| '_ \ _____|  _ \| |/ _ \| '__| |/ /
    //   / ___ \| | | | (_| |  __/ |  \__ \ (_) | | | |_____| |_) | | (_) | |  |   <
    //  /_/   \_\_| |_|\__,_|\___|_|  |___/\___/|_| |_|     |____// |\___/|_|  |_|\_\
    //                                                          |__/
    */

    auto modified_AB = [this,check] () -> void {
    
      /*
      N. Ganchovski, A. Traykov
      Modified Anderson-Bjork’s method for solving nonlinear equations in structural mechanics
      IOP Conf. Series: Materials Science and Engineering
      1276 (2023) 012010 doi:10.1088/1757-899X/1276/1/012010
      */
    
      Integer side{0};
      // Integer N{ Integer(floor(1-log2( m_tolerance )/2)) };

      Real & x1{ m_a  };
      Real & f1{ m_fa };
      Real & x2{ m_b  };
      Real & f2{ m_fb };
      
      bool bisection{true};

      while ( ++m_iteration_count < m_max_iteration ) {

        Real dir{ x2 - x1 };
        Real tol{ m_tolerance / ( 2 * abs( dir ) ) };
        m_converged = tol >= 0.5;
        if ( m_converged ) break;
        
        Real x3, f3;
        
        if ( bisection ) {
          x3 = x1 + dir/2;   // Midpoint abscissa
          f3 = evaluate(x3); // and function value

          check( x3, f3 );
          m_converged = f3 == 0;
          if ( m_converged ) { x1 = x2 = x3; f1 = f2 = 0; break; }
          
          Real fm{ (f1+f2)/2 }; // Ordinate of chord at midpoint
          bisection = abs(fm-f3) >= (abs(fm)+abs(f3))/4;
        } else {

          Real t{ f1/(f1-f2) };
          if      ( t < tol   ) t = tol;
          else if ( t > 1-tol ) t = 1-tol;

          x3 = x1 + t*dir; // False-position step
          f3 = evaluate(x3);

          check( x3, f3 );
          m_converged = f3 == 0;
          if ( m_converged ) { x1 = x2 = x3; f1 = f2 = 0; break; }
        }

        switch ( side ) {
        case 1: // Apply Anderson-Bjork modification for side 1
          { Real m{ 1-f3/f1 }; f2 *= m > Utils::machine_eps<Real>() ? m : 0.5; }
          break;
        case 2: // Apply Anderson-Bjork modification for side 2
          { Real m{ 1-f3/f2 }; f1 *= m > Utils::machine_eps<Real>() ? m : 0.5; }
          break;
        }
        
        if ( (f1 > 0) == (f3 > 0 ) ) { // If the left interval does not change sign
          if ( !bisection ) side = 1;  // Store the side that move
          x1 = x3; f1 = f3;            // Move the left end
        } else {                       // If the right interval does not change sign
          if ( !bisection ) side = 2;  // Store the side that move
          x2 = x3; f2 = f3;            // Move the right end
        }
        
        //if ( (m_iteration_count % N) == 0 ) {
        //  bisection = true;
        //  side = 0;
        //}
      }
    };
    
    switch ( m_select ) {
      case Method::BISECTION:    bisection();    break;
      case Method::ILLINOIS:     illinois();     break;
      case Method::CHANDRUPATLA: Chandrupatla(); break;
      case Method::BRENT:        Brent();        break;
      case Method::RIDDER:       Ridder();       break;
      case Method::MODIFIED_AB:  modified_AB();  break;
    }

    if ( m_a > m_b ) { swap(m_a,m_b); swap(m_fa, m_fb); }

    return abs(m_fa) < abs(m_fb) ? m_a : m_b;
  }

  template class AlgoBracket<float>;
  template class AlgoBracket<double>;
}
