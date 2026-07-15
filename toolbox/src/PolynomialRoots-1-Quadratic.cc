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

namespace PolynomialRoots
{

  using std::abs;
  using std::sqrt;

  template <typename T_real, typename T_complex>
  integer QuadraticT<T_real, T_complex>::get_real_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( !m_cplx )
    {
      r[nr++] = m_r0;
      if ( m_nrts > 1 ) r[nr++] = m_r1;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuadraticT<T_real, T_complex>::get_positive_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( !m_cplx )
    {
      if ( m_nrts > 0 && m_r0 > 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > 0 ) r[nr++] = m_r1;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuadraticT<T_real, T_complex>::get_negative_roots( T_real r[] ) const
  {
    integer nr = 0;
    if ( !m_cplx )
    {
      if ( m_nrts > 0 && m_r0 < 0 ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 < 0 ) r[nr++] = m_r1;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuadraticT<T_real, T_complex>::get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    integer nr = 0;
    if ( !m_cplx )
    {
      if ( m_nrts > 0 && m_r0 >= a && m_r0 <= b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 >= a && m_r1 <= b ) r[nr++] = m_r1;
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuadraticT<T_real, T_complex>::get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    integer nr = 0;
    if ( !m_cplx )
    {
      if ( m_nrts > 0 && m_r0 > a && m_r0 < b ) r[nr++] = m_r0;
      if ( m_nrts > 1 && m_r1 > a && m_r1 < b ) r[nr++] = m_r1;
    }
    return nr;
  }

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/

  template <typename T_real, typename T_complex> void QuadraticT<T_real, T_complex>::find_roots()
  {
    T_real const & A = m_ABC[0];
    T_real const & B = m_ABC[1];
    T_real const & C = m_ABC[2];

    m_r0 = m_r1 = 0;
    m_nrts      = 0;
    m_cplx = m_dblx = false;

    if ( A == 0 )
    {  // less than two roots b*z + c = 0
      if ( B != 0 )
      {
        m_nrts = 1;
        m_r0   = -C / B;
      }
    }
    else if ( C == 0 )
    {  // a*z^2 + b*z  = 0
      m_nrts = 2;
      m_dblx = B == 0;
      if ( !m_dblx )
      {
        m_r0 = -B / A;
        if ( m_r0 < 0 ) std::swap( m_r0, m_r1 );
      }
    }
    else
    {                              // Compute discriminant avoiding overflow.
      T_real const hb    = B / 2;  // b now b/2
      T_real const abs_b = abs( hb );
      T_real const abs_c = abs( C );
      T_real       e;
      T_real       d;
      if ( abs_b < abs_c )
      {
        e = C < 0 ? -A : A;
        e = ( hb * hb ) - e * abs_c;
        d = abs( e );
        d = sqrt( d );
      }
      else
      {
        e = 1 - ( A / hb ) * ( C / hb );
        d = sqrt( abs( e ) ) * abs_b;
      }
      m_nrts = 2;
      m_cplx = e < 0;
      if ( m_cplx )
      {
        // complex conjugate zeros
        m_r0 = -hb / A;       // real part
        m_r1 = abs( d / A );  // immaginary part
      }
      else
      {
        // real zeros
        m_dblx = d == 0;
        if ( m_dblx ) { m_r0 = m_r1 = -hb / A; }
        else
        {
          if ( hb >= 0 ) d = -d;
          m_r0 = ( d - hb ) / A;
          // r1 = (-d-hb)/a;
          if ( m_r0 != 0 ) m_r1 = ( C / m_r0 ) / A;
          if ( m_r0 > m_r1 ) std::swap( m_r0, m_r1 );  // order roots
        }
      }
    }
  }

  template <typename T_real, typename T_complex> void QuadraticT<T_real, T_complex>::info( ostream_type & s ) const
  {
    T_real const & A = m_ABC[0];
    T_real const & B = m_ABC[1];
    T_real const & C = m_ABC[2];
    s << std::format(
      "poly A={} B={} C={}\n"
      "n. roots = {}\n"
      "complex  = {}\n"
      "double   = {}\n",
      A,
      B,
      C,
      m_nrts,
      ( m_cplx ? "YES" : "NO" ),
      ( m_dblx ? "YES" : "NO" ) );
    if ( m_cplx )
      s << std::format(
        "x0 = ({},{})\n"
        "x1 = ({},{})\n",
        m_r0,
        m_r1,
        m_r0,
        -m_r1 );
    else if ( m_dblx )
      s << std::format( "x0 = x1 = {}\n", m_r0 );
    else if ( m_nrts == 1 )
      s << std::format( "x0 = {}\n", m_r0 );
    else if ( m_nrts == 2 )
      s << std::format(
        "x0 = {}\n"
        "x1 = {}\n",
        m_r0,
        m_r1 );
  }

  template <typename T_real, typename T_complex> bool QuadraticT<T_real, T_complex>::check( ostream_type & s ) const
  {
    T_real const & A    = m_ABC[0];
    T_real const & B    = m_ABC[1];
    T_real const & C    = m_ABC[2];
    bool           ok   = true;
    T_real const   epsi = ( abs( A ) + abs( B ) + abs( C ) ) * toleranceT<T_real>();
    if ( m_cplx )
    {
      T_real const z0 = abs( eval( root0() ) );
      T_real const z1 = abs( eval( root1() ) );
      s << std::format(
        "|p(r0)| = {}\n"
        "|p(r1)| = {}\n",
        z0,
        z1 );
      ok = z0 < epsi && z1 < epsi;
    }
    else if ( m_nrts == 1 )
    {
      T_real const z0 = eval( real_root0() );
      s << std::format( "p(r0) = {}\n", z0 );
      ok = abs( z0 ) < epsi;
    }
    else if ( m_nrts == 2 )
    {
      T_real const z0 = eval( real_root0() );
      T_real const z1 = eval( real_root1() );
      s << std::format(
        "p(r0) = {}\n"
        "p(r1) = {}\n",
        z0,
        z1 );
      ok = abs( z0 ) < epsi && abs( z1 ) < epsi;
    }
    return ok;
  }

  template class QuadraticT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template class QuadraticT<quad_real, quad_complex>;
#endif

}  // namespace PolynomialRoots

// EOF: PolynomialRoots-1-Quadratic.cc
