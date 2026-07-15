/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

#include "PolynomialRoots.hh"

#define MAX_ITER_SAFETY 50

namespace PolynomialRoots
{

  using std::abs;
  using std::cbrt;
  using std::sqrt;

  template <typename T_real> inline T_real thirdT()
  { return T_real( 1 ) / T_real( 3 ); }
  template <typename T_real> inline T_real one27thT()
  { return T_real( 1 ) / T_real( 27 ); }
  template <typename T_real> inline T_real two27thT()
  { return T_real( 2 ) / T_real( 27 ); }

  template <typename T_real, typename T_complex> integer CubicT<T_real, T_complex>::get_real_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( m_cplx )
    {
      if ( m_nrts > 2 ) r[nr++] = m_r2;
    }
    else
    {
      if ( m_nrts > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 ) r[nr++] = m_r1;
      if ( m_nrts > 2 ) r[nr++] = m_r2;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer CubicT<T_real, T_complex>::get_positive_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( m_cplx )
    {
      if ( m_nrts > 2 && m_r2 > 0 ) r[nr++] = m_r2;
    }
    else
    {
      if ( m_nrts > 0 && m_r0 > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > 0 ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 > 0 ) r[nr++] = m_r2;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer CubicT<T_real, T_complex>::get_negative_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( m_cplx )
    {
      if ( m_nrts > 2 && m_r2 < 0 ) r[nr++] = m_r2;
    }
    else
    {
      if ( m_nrts > 0 && m_r0 < 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 < 0 ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 < 0 ) r[nr++] = m_r2;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer CubicT<T_real, T_complex>::get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    integer nr = 0;
    if ( m_cplx )
    {
      if ( m_nrts > 2 && m_r2 >= a && m_r2 <= b ) r[nr++] = m_r2;
    }
    else
    {
      if ( m_nrts > 0 && m_r0 >= a && m_r0 <= b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 >= a && m_r1 <= b ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 >= a && m_r2 <= b ) r[nr++] = m_r2;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer CubicT<T_real, T_complex>::get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    integer nr = 0;
    if ( m_cplx )
    {
      if ( m_nrts > 2 && m_r2 > a && m_r2 < b ) r[nr++] = m_r2;
    }
    else
    {
      if ( m_nrts > 0 && m_r0 > a && m_r0 < b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > a && m_r1 < b ) r[nr++] = m_r1;
      if ( m_nrts > 2 && m_r2 > a && m_r2 < b ) r[nr++] = m_r2;
    }
    return nr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T_real> static T_real guess1( T_real const a[3] )
  {
    static T_real const p = 1.09574;
    static T_real const q = -3.239E-1;
    static T_real const r = -3.239E-1;
    static T_real const s = 9.57439E-2;
    return p + q * a[1] + r * a[2] + s * a[1] * a[2];
  }

  template <typename T_real> static T_real guess2( T_real const a[3] )
  {
    static T_real const p = -1.09574;
    static T_real const q = 3.239E-1;
    static T_real const r = -3.239E-1;
    static T_real const s = 9.57439E-2;
    return p + q * a[1] + r * a[2] + s * a[1] * a[2];
  }

  template <typename T_real> static T_real guess3( T_real const a[3] )
  {
    static T_real const p = 1.14413;
    static T_real const q = -2.75509E-1;
    static T_real const r = -4.45578E-1;
    static T_real const s = -2.59342E-2;
    if ( T_real const t = a[2] / 3; a[0] < t * ( 2 * t * t - 1 ) ) return p + q * a[0] + r * a[2] + s * a[0] * a[2];
    return -p + q * a[0] + r * a[2] - s * a[0] * a[2];
  }

  template <typename T_real> static T_real guess4( T_real const a[3] )
  {
    static T_real const q = -7.71845E-1;
    static T_real const s = -2.28155E-1;
    if ( a[0] > 0 ) return ( q + s * a[2] ) * a[0];
    return ( q - s * a[2] ) * a[0];
  }

  template <typename T_real> static T_real guess5( T_real const a[3] )
  {
    T_real       p, q, r, s;
    T_real const tmp = two27thT<T_real>() - a[1] / 3;
    if ( a[1] <= thirdT<T_real>() )
    {
      if ( a[0] < tmp )
      {
        p = 8.78558E-1;
        q = -5.71888E-1;
        r = -7.11154E-1;
        s = -3.22313E-1;
      }
      else
      {
        p = -1.92823E-1;
        q = -5.66324E-1;
        r = +5.05734E-1;
        s = -2.64881E-1;
      }
    }
    else
    {
      if ( a[0] < tmp )
      {
        p = 1.19748;
        q = -2.83772E-1;
        r = -8.37476E-1;
        s = -3.56228E-1;
      }
      else
      {
        p = -3.45219E-1;
        q = -4.01231E-1;
        r = 2.07216E-1;
        s = -4.45532E-3;
      }
    }
    return p + q * a[0] + r * a[1] + s * a[0] * a[1];
  }

  template <typename T_real> static T_real guess6( T_real const a[3] )
  {
    T_real       p, q, r, s;
    T_real const tmp = a[1] / 3 - two27thT<T_real>();
    if ( a[1] <= thirdT<T_real>() )
    {
      if ( a[0] > tmp )
      {
        p = -8.78558E-1;
        q = -5.71888E-1;
        r = 7.11154E-1;
        s = -3.22313E-1;
      }
      else
      {
        p = 1.92823E-1;
        q = -5.66324E-1;
        r = -5.05734E-1;
        s = -2.64881E-1;
      }
    }
    else
    {
      if ( a[0] > tmp )
      {
        p = -1.19748;
        q = -2.83772E-1;
        r = 8.37476E-1;
        s = -3.56228E-1;
      }
      else
      {
        p = 3.45219E-1;
        q = -4.01231E-1;
        r = -2.07216E-1;
        s = -4.45532E-3;
      }
    }
    return p + q * a[0] + r * a[1] + s * a[0] * a[1];
  }

  /*
  ||   _   _               _              ____  _               _   _
  ||  | \ | | _____      _| |_ ___  _ __ | __ )(_)___  ___  ___| |_(_) ___  _ __
  ||  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \|  _ \| / __|/ _ \/ __| __| |/ _ \| '_ \
  ||  | |\  |  __/\ V  V /| || (_) | | | | |_) | \__ \  __/ (__| |_| | (_) | | | |
  ||  |_| \_|\___| \_/\_/  \__\___/|_| |_|____/|_|___/\___|\___|\__|_|\___/|_| |_|
  */

  // x^3 + a * x^2 + b * x + c
  template <typename T_real>
  static integer NewtonBisection( T_real const a, T_real const b, T_real const c, T_real & x )
  {
    T_real p, dp;
    evalMonicCubic( x, a, b, c, p, dp );
    T_real t = p;  // save p(x) for sign comparison
    x -= p / dp;   // 1st improved root

    integer iter       = 1;
    integer oscillate  = 0;
    integer nconverged = 0;
    bool    bisection  = false;
    bool    converged  = false;
    T_real  s          = 0;
    T_real  u          = 0;  // to mute warning
    while ( !( nconverged > 1 || bisection ) && iter < MAX_ITER_SAFETY )
    {
      ++iter;
      evalMonicCubic( x, a, b, c, p, dp );
      if ( p * t < 0 )
      {  // does Newton start oscillating ?
        if ( p < 0 )
        {
          ++oscillate;  // increment oscillation counter
          s = x;        // save lower bisection bound
        }
        else
        {
          u = x;  // save upper bisection bound
        }
        t = p;  // save current p(x)
      }
      dp = p / dp;                                                      // Newton correction
      x -= dp;                                                          // new Newton root
      bisection = oscillate > 2;                                        // activate bisection
      converged = abs( dp ) <= ( 1 + abs( x ) ) * machepsiT<T_real>();  // Newton convergence indicator
      if ( converged )
        ++nconverged;
      else
        nconverged = 0;
    }
    if ( bisection )
    {
      t = u - s;  // initial bisection interval
      while ( abs( t ) > ( 1 + abs( x ) ) * machepsiT<T_real>() && iter < MAX_ITER_SAFETY )
      {  // bisection iterates
        ++iter;
        p = evalMonicCubic( x, a, b, c );
        if ( p < 0 )
          s = x;
        else
          u = x;            // keep bracket on root
        t = ( u - s ) / 2;  // new bisection interval
        x = s + t;          // new bisection root
      }
    }
    return iter;
  }

  /*\
   *  Calculate the zeros of the cubic a*z^3 + b*z^2 + c*z + d.
  \*/

  template <typename T_real, typename T_complex> void CubicT<T_real, T_complex>::find_roots()
  {
    T_real const & A = m_ABCD[0];
    T_real const & B = m_ABCD[1];
    T_real const & C = m_ABCD[2];
    T_real const & D = m_ABCD[3];
    m_nrts = m_iter = 0;
    m_cplx = m_dblx = m_trpx = false;
    // special cases
    if ( A == 0 )
    {
      QuadraticT<T_real, T_complex> const qsolve( B, C, D );
      m_nrts = qsolve.num_roots();
      m_cplx = qsolve.complex_root();
      m_dblx = qsolve.double_root();
      m_r0   = qsolve.real_root0();
      m_r1   = qsolve.real_root1();
      return;
    }
    if ( D == 0 )
    {
      QuadraticT<T_real, T_complex> const qsolve( A, B, C );
      m_nrts = qsolve.num_roots() + 1;
      m_cplx = qsolve.complex_root();
      m_r0   = qsolve.real_root0();
      m_r1   = qsolve.real_root1();
      m_r2   = 0;
      if ( !m_cplx )
      {  // reorder
        if ( m_r2 < m_r1 ) std::swap( m_r1, m_r2 );
        if ( m_r1 < m_r0 ) std::swap( m_r0, m_r1 );
        if ( m_r2 < m_r1 ) std::swap( m_r1, m_r2 );
      }
      return;
    }
    /*              _
    ||  ___ __ __ _| |___
    || (_-</ _/ _` | / -_)
    || /__/\__\__,_|_\___|
    */
    // x^3 + aa * x^2 + bb * x + cc
    T_real const aa = B / A;
    T_real const bb = C / A;
    T_real const cc = D / A;
    // scale Cubic Monic Polynomial
    T_real const absa = abs( aa );
    T_real const absb = sqrt( abs( bb ) );
    T_real const absc = cbrt( abs( cc ) );

    integer i_case = 0;  // c MAX
    if ( absa < absb )
    {
      if ( absc < absb ) i_case = 1;  // |a| < |b| and |b| < |c| --> b MAX
      // |a| < |b| <= |c| --> c MAX
    }
    else
    {
      if ( absc < absa ) i_case = 2;  // |b| <= |a| and |c| < |a| --> a MAX
      // |b| <= |a| < |c| --> c MAX
    }

    T_real scale = 0;
    T_real a[3];
    switch ( i_case )
    {
      case 0:
        scale = absc;
        a[2]  = aa / absc;
        a[1]  = ( bb / absc ) / absc;
        a[0]  = cc > 0 ? 1 : -1;
        break;
      case 1:
        scale = absb;
        a[2]  = aa / absb;
        a[1]  = bb > 0 ? 1 : -1;
        a[0]  = ( ( cc / absb ) / absb ) / absb;
        break;
      case 2:
        scale = absa;
        a[2]  = aa > 0 ? 1 : -1;
        a[1]  = ( bb / absa ) / absa;
        a[0]  = ( ( cc / absa ) / absa ) / absa;
        break;
    }

    /*
    ||   __ _ _  _ ___ ______
    ||  / _` | || / -_|_-<_-<
    ||  \__, |\_,_\___/__/__/
    ||  |___/
    */
    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    integer iclass = -1;
    switch ( i_case )
    {
      case 0: iclass = a[0] > 0 ? 2 : 1; break;
      case 1: iclass = a[1] > 0 ? 4 : 3; break;
      case 2: iclass = a[2] > 0 ? 6 : 5; break;
    }
    bool use_shifted = false;
    m_trpx           = false;
    switch ( iclass )
    {
      case 1: m_r2 = guess1<T_real>( a ); break;
      case 2: m_r2 = guess2<T_real>( a ); break;
      case 3: m_r2 = guess3<T_real>( a ); break;
      case 4: m_r2 = guess4<T_real>( a ); break;
      case 5:
        m_r0   = a[1] - thirdT<T_real>();
        m_r1   = a[0] + one27thT<T_real>();
        m_trpx = abs( m_r0 ) <= machepsiT<T_real>() && abs( m_r1 ) <= machepsiT<T_real>();  // check for triple root
        if ( m_trpx )
        {
          m_r0 = m_r1 = m_r2 = thirdT<T_real>() * scale;
          m_nrts             = 3;
          return;
        }
        use_shifted = abs( m_r0 ) <= 0.01 && abs( m_r1 ) <= 0.01;
        m_r2        = guess5( a );
        break;
      case 6:
        m_r0   = a[1] - thirdT<T_real>();
        m_r1   = a[0] - one27thT<T_real>();
        m_trpx = abs( m_r0 ) <= machepsiT<T_real>() && abs( m_r1 ) <= machepsiT<T_real>();  // check for triple root
        if ( m_trpx )
        {
          m_r0 = m_r1 = m_r2 = -thirdT<T_real>() * scale;
          m_nrts             = 3;
          return;
        }
        use_shifted = abs( m_r0 ) <= 0.01 && abs( m_r1 ) <= 0.01;
        m_r2        = guess6( a );
        break;
    }

    /*
    ||          _
    ||  ___ ___| |_ _____
    || (_-</ _ \ \ V / -_)
    || /__/\___/_|\_/\___|
    */
    m_iter = 0;
    if ( use_shifted )
    {
      if ( iclass == 5 )
      {
        // a[2] == -1
        // y^3 + (a[1]-1/3)* y + (a[0]+a[1]/3-2/27), x = y+1/3
        m_r2 -= thirdT<T_real>();  // shift guess
        m_iter = NewtonBisection<T_real>( 0, m_r0, a[0] + a[1] / 3 - two27thT<T_real>(), m_r2 );
        m_r2 += thirdT<T_real>();  // unshift solution
      }
      else
      {
        // a[2] == 1
        // y^3 + (a[1]-1/3)* y + (a[0]-a[1]/3+2/27), x = y+1/3
        m_r2 += thirdT<T_real>();  // shift guess
        m_r1 -= a[1] / 3 - one27thT<T_real>();
        // if ( std::abs(r3) < machepsi ) r3 = 0;
        m_iter = NewtonBisection<T_real>( 0, m_r0, a[0] - a[1] / 3 + two27thT<T_real>(), m_r2 );
        m_r2 -= thirdT<T_real>();  // unshift solution
      }
    }
    else
    {
      m_iter = NewtonBisection<T_real>( a[2], a[1], a[0], m_r2 );
    }

    // scale root
    m_r2 *= scale;

    /*
    // deflate
    // x^3 + aa*x^2 + bb*x + cc
    //    = (x-r2)*(x^2+b1*x+b0)
    //    = x^3 + x^2 * ( b1 - r2 ) + x * ( b0 - b1*r2 ) - r2 * b0
    //
    //    aa == b1 - r2
    //    bb == b0 - b1*r2
    //    cc == -r2 * b0
    //
    //  Solve the overdetermined linear system:
    //
    //  / 0    1  \            / aa + r2 \
    //  |         |  / b0 \    |         |
    //  | 1   -r2 |  |    |  = |   bb    |
    //  |         |  \ b1 /    |         |
    //  \ -r2  0  /            \   cc    /
    //
    //  if |r2| < 1 then solve
    //
    //  / 0    1  \  / b0 \    / aa + r2 \
    //  |         |  |    | =  |         |
    //  \ -r2  0  /  \ b1 /    \   cc    /
    //
    //  otherwise solve
    //
    //  / 1   -r2 \  / b0 \   / bb \
    //  |         |  |    | = |    |
    //  \ -r2  0  /  \ b1 /   \ cc /
    */
    T_real const b0 = -cc / m_r2;
    T_real const b1 = abs( m_r2 ) < 1 ? aa + m_r2 : -( cc / m_r2 + bb ) / m_r2;

    // solve quadratic polynomial
    QuadraticT<T_real, T_complex> const qsolve( 1.0, b1, b0 );
    m_nrts = qsolve.num_roots() + 1;
    m_cplx = qsolve.complex_root();
    m_dblx = qsolve.double_root();
    m_r0   = qsolve.real_root0();
    m_r1   = qsolve.real_root1();

    if ( !m_cplx )
    {  // if real roots sort it!
      if ( m_r1 > m_r2 ) std::swap( m_r1, m_r2 );
      if ( m_r0 > m_r1 ) std::swap( m_r0, m_r1 );
      if ( m_r1 > m_r2 ) std::swap( m_r1, m_r2 );
    }
  }

  template <typename T_real, typename T_complex> void CubicT<T_real, T_complex>::info( ostream_type & s ) const
  {
    T_real const & A = m_ABCD[0];
    T_real const & B = m_ABCD[1];
    T_real const & C = m_ABCD[2];
    T_real const & D = m_ABCD[3];
    s << std::format(
      "poly a={} b={} c={} d={}\n"
      "n. roots = {}\n"
      "complex  = {}\n"
      "triple   = {}\n"
      "double   = {}\n",
      A,
      B,
      C,
      D,
      m_nrts,
      ( m_cplx ? "YES" : "NO" ),
      ( m_trpx ? "YES" : "NO" ),
      ( m_dblx ? "YES" : "NO" ) );
    if ( m_cplx )
    {
      s << std::format(
        "x0 = ({},{})\n"
        "x1 = ({},{})\n",
        m_r0,
        m_r1,
        m_r0,
        -m_r1 );
      if ( m_nrts > 2 ) s << std::format( "x3 = {}\n", m_r2 );
    }
    else
    {
      if ( m_nrts > 0 ) s << std::format( "x0 = {}\n", m_r0 );
      if ( m_nrts > 1 ) s << std::format( "x1 = {}\n", m_r1 );
      if ( m_nrts > 2 ) s << std::format( "x2 = {}\n", m_r2 );
    }
    s << '\n';
  }

  template <typename T_real, typename T_complex> bool CubicT<T_real, T_complex>::check( ostream_type & s ) const
  {
    T_real const & A    = m_ABCD[0];
    T_real const & B    = m_ABCD[1];
    T_real const & C    = m_ABCD[2];
    T_real const & D    = m_ABCD[3];
    bool           ok   = true;
    T_real const   epsi = ( abs( A ) + abs( B ) + abs( C ) + abs( D ) ) * toleranceT<T_real>();
    if ( m_cplx )
    {
      T_real const z0 = abs( eval( root0() ) );
      T_real const z1 = abs( eval( root1() ) );
      T_real const z2 = abs( eval( root2() ) );
      T_real const zr = eval( real_root2() );
      s << std::format(
        "|p(r0)| = {}\n"
        "|p(r1)| = {}\n"
        "|p(r2)| = {}\n"
        "p(real_part(r2)) = {}\n",
        z0,
        z1,
        z2,
        zr );
      ok = z0 < epsi && z1 < epsi && z2 < epsi;
    }
    else if ( m_nrts == 1 )
    {
      T_real const z0 = eval( real_root0() );
      s << std::format( "p(r0) = {}\n", z0 );
      ok = abs( z0 ) < epsi;
    }
    else if ( m_nrts == 2 )
    {
      T_real const z0 = abs( eval( root0() ) );
      T_real const z1 = abs( eval( root1() ) );
      s << std::format(
        "p(r0) = {}\n"
        "p(r1) = {}\n",
        z0,
        z1 );
      ok = abs( z0 ) < epsi && abs( z1 ) < epsi;
    }
    else if ( m_nrts == 3 )
    {
      T_real const z0 = eval( real_root0() );
      T_real const z1 = eval( real_root1() );
      T_real const z2 = eval( real_root2() );
      s << std::format(
        "p(r0) = {}\n"
        "p(r1) = {}\n"
        "p(r2) = {}\n",
        z0,
        z1,
        z2 );
      ok = abs( z0 ) < epsi && abs( z1 ) < epsi && abs( z2 ) < epsi;
    }
    return ok;
  }

  template class CubicT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template class CubicT<quad_real, quad_complex>;
#endif

}  // namespace PolynomialRoots

// EOF: PolynomialRoots-2-Cubic.cc
