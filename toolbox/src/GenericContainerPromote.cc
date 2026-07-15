/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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

//
// file: GenericContainer.cc
//

#ifdef _MSC_VER
#pragma warning( disable : 4661 )
#pragma warning( disable : 4244 )
#endif

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <cmath>
#include <algorithm>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#if __cplusplus >= 201103L || ( defined( _MSC_VER ) && _MSC_VER >= 1900 )
#else
#error This library needs at least a C++11 compliant compiler
#endif

namespace GC_namespace
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  using std::fpclassify;

  static bool isZero0( real_type const x )
  {
    int const c{ fpclassify( x ) };
    return FP_ZERO == c || FP_SUBNORMAL == c;
  }

#endif

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_int( vec_int_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t ne = get_num_elements();
    v.reserve( ne );
    long_type    lval;
    real_type    rval;
    complex_type cval;
    int_type     val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = _b() ? 1 : 0; break;
        case GC_type::INTEGER: val = _i(); break;
        case GC_type::LONG:
          lval = _l();
          GC_assert(
            std::in_range<int_type>( lval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            lval );
          val = static_cast<int_type>( lval );
          break;
        case GC_type::REAL:
          rval = _r();
          GC_assert(
            GC_details::real_fits_integral<int_type>( rval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            rval );
          val = static_cast<int_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = _v_b()[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = _v_i()[i]; break;
        case GC_type::VEC_LONG:
          lval = _v_l()[i];
          GC_assert(
            std::in_range<int_type>( lval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            lval );
          val = static_cast<int_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = _v_r()[i];
          GC_assert(
            GC_details::real_fits_integral<int_type>( rval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            rval );
          val = static_cast<int_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = _c();
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<int_type>( cval.real() ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = _v_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<int_type>( cval.real() ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER: val = _m_i()[i]; break;
        case GC_type::MAT_LONG:
          lval = _m_l()[i];
          GC_assert(
            std::in_range<int_type>( lval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            lval );
          val = static_cast<int_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = _m_r()[i];
          GC_assert(
            GC_details::real_fits_integral<int_type>( rval ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            rval );
          val = static_cast<int_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = _m_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<int_type>( cval.real() ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_int( "GenericContainer::copyto_vec_int " ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_int_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_uint( vec_uint_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t ne{ get_num_elements() };
    v.reserve( ne );
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    uint_type    val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = _b() ? 1 : 0; break;
        case GC_type::INTEGER:
          ival = _i();
          GC_assert(
            std::in_range<uint_type>( ival ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            ival );
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::LONG:
          lval = _l();
          GC_assert(
            std::in_range<uint_type>( lval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            lval );
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::REAL:
          rval = _r();
          GC_assert(
            GC_details::real_fits_integral<uint_type>( rval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            rval );
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = _v_b()[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER:
          ival = _v_i()[i];
          GC_assert(
            std::in_range<uint_type>( ival ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            ival );
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::VEC_LONG:
          lval = _v_l()[i];
          GC_assert(
            std::in_range<uint_type>( lval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            lval );
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = _v_r()[i];
          GC_assert(
            GC_details::real_fits_integral<uint_type>( rval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            rval );
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = _c();
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<uint_type>( cval.real() ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = _v_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<uint_type>( cval.real() ),
            "{} copyto_vec_int: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER:
          ival = _m_i()[i];
          GC_assert(
            std::in_range<uint_type>( ival ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            ival );
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::MAT_LONG:
          lval = _m_l()[i];
          GC_assert(
            std::in_range<uint_type>( lval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            lval );
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = _m_r()[i];
          GC_assert(
            GC_details::real_fits_integral<uint_type>( rval ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            rval );
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = _m_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<uint_type>( cval.real() ),
            "{} copyto_vec_uint: v[{}] = {} cannot be converted to `unsigned integer'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_uint( "GenericContainer::copyto_vec_uint " ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_uint_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_long( vec_long_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t ne{ get_num_elements() };
    v.reserve( ne );
    real_type    rval;
    complex_type cval;
    long_type    val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = _b() ? 1 : 0; break;
        case GC_type::INTEGER: val = static_cast<long_type>( _i() ); break;
        case GC_type::LONG: val = _l(); break;
        case GC_type::REAL:
          rval = _r();
          GC_assert(
            GC_details::real_fits_integral<long_type>( rval ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            rval );
          val = static_cast<long_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = _v_b()[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = static_cast<long_type>( _v_i()[i] ); break;
        case GC_type::VEC_LONG: val = _v_l()[i]; break;
        case GC_type::VEC_REAL:
          rval = _v_r()[i];
          GC_assert(
            GC_details::real_fits_integral<long_type>( rval ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            rval );
          val = static_cast<long_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = _c();
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<long_type>( cval.real() ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = _v_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<long_type>( cval.real() ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER: val = static_cast<long_type>( _m_i()[i] ); break;
        case GC_type::MAT_LONG: val = _m_l()[i]; break;
        case GC_type::MAT_REAL:
          rval = _m_r()[i];
          GC_assert(
            GC_details::real_fits_integral<long_type>( rval ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            rval );
          val = static_cast<long_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = _m_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<long_type>( cval.real() ),
            "{} copyto_vec_long: v[{}] = {} cannot be converted to `long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_long( "copyto_vec_long" ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_long_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_ulong( vec_ulong_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t ne{ get_num_elements() };
    v.reserve( ne );
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    ulong_type   val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = _b() ? 1 : 0; break;
        case GC_type::INTEGER:
          ival = _i();
          GC_assert(
            std::in_range<ulong_type>( ival ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            ival );
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::LONG:
          lval = _l();
          GC_assert(
            std::in_range<ulong_type>( lval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            lval );
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::REAL:
          rval = _r();
          GC_assert(
            GC_details::real_fits_integral<ulong_type>( rval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            rval );
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = _v_b()[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER:
          ival = _v_i()[i];
          GC_assert(
            std::in_range<ulong_type>( ival ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            ival );
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::VEC_LONG:
          lval = _v_l()[i];
          GC_assert(
            std::in_range<ulong_type>( lval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            lval );
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = _v_r()[i];
          GC_assert(
            GC_details::real_fits_integral<ulong_type>( rval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            rval );
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = _c();
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<ulong_type>( cval.real() ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = _v_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<ulong_type>( cval.real() ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER:
          ival = _m_i()[i];
          GC_assert(
            std::in_range<ulong_type>( ival ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            ival );
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::MAT_LONG:
          lval = _m_l()[i];
          GC_assert(
            std::in_range<ulong_type>( lval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            lval );
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = _m_r()[i];
          GC_assert(
            GC_details::real_fits_integral<ulong_type>( rval ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            rval );
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = _m_c()[i];
          GC_assert(
            isZero0( cval.imag() ) && GC_details::real_fits_integral<ulong_type>( cval.real() ),
            "{} copyto_vec_ulong: v[{}] = {} cannot be converted to `std::size_t long'",
            where,
            i,
            to_string( cval ) );
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_ulong( "copyto_vec_ulong" ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_ulong_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_real( vec_real_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t ne{ get_num_elements() };
    v.reserve( ne );
    complex_type cval;
    real_type    val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = _b() ? 1 : 0; break;
        case GC_type::INTEGER: val = static_cast<real_type>( _i() ); break;
        case GC_type::LONG: val = static_cast<real_type>( _l() ); break;
        case GC_type::REAL: val = _r(); break;
        case GC_type::VEC_BOOL: val = _v_b()[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = static_cast<real_type>( _v_i()[i] ); break;
        case GC_type::VEC_LONG: val = static_cast<real_type>( _v_l()[i] ); break;
        case GC_type::VEC_REAL: val = _v_r()[i]; break;
        case GC_type::COMPLEX:
          cval = _c();
          GC_assert(
            isZero0( cval.imag() ),
            "{} copyto_vec_real: v[{}] = {} cannot be converted to `real_type'",
            where,
            i,
            to_string( cval ) );
          val = cval.real();
          break;
        case GC_type::VEC_COMPLEX:
          cval = _v_c()[i];
          GC_assert(
            isZero0( cval.imag() ),
            "{} copyto_vec_real: v[{}] = {} cannot be converted to `real_type'",
            where,
            i,
            to_string( cval ) );
          val = cval.real();
          break;
        case GC_type::MAT_INTEGER: val = static_cast<real_type>( _m_i()[i] ); break;
        case GC_type::MAT_LONG: val = static_cast<real_type>( _m_l()[i] ); break;
        case GC_type::MAT_REAL: val = _m_r()[i]; break;
        case GC_type::MAT_COMPLEX:
          cval = _m_c()[i];
          GC_assert(
            isZero0( cval.imag() ),
            "{} copyto_vec_real: v[{}] = {} cannot be converted to `real_type'",
            where,
            i,
            to_string( cval ) );
          val = cval.real();
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_number(); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_real_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_complex( vec_complex_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t const ne{ get_num_elements() };
    v.reserve( ne );
    complex_type val{ 0 };
    v.reserve( ne );
    for ( std::size_t i{ 0 }; i < ne; ++i )
    {
      switch ( get_type() )
      {
        case GC_type::BOOL: val = { static_cast<real_type>( _b() ? 1 : 0 ), 0 }; break;
        case GC_type::INTEGER: val = { static_cast<real_type>( _i() ), 0 }; break;
        case GC_type::LONG: val = { static_cast<real_type>( _l() ), 0 }; break;
        case GC_type::REAL: val = { _r(), 0 }; break;
        case GC_type::VEC_BOOL: val = { static_cast<real_type>( _v_b()[i] ? 1 : 0 ), 0 }; break;
        case GC_type::VEC_INTEGER: val = { static_cast<real_type>( _v_i()[i] ), 0 }; break;
        case GC_type::VEC_LONG: val = { static_cast<real_type>( _v_l()[i] ), 0 }; break;
        case GC_type::VEC_REAL: val = { _v_r()[i], 0 }; break;
        case GC_type::COMPLEX: val = _c(); break;
        case GC_type::VEC_COMPLEX: val = _v_c()[i]; break;
        case GC_type::MAT_INTEGER: val = { static_cast<real_type>( _m_i()[i] ), 0 }; break;
        case GC_type::MAT_LONG: val = { static_cast<real_type>( _m_l()[i] ), 0 }; break;
        case GC_type::MAT_REAL: val = { _m_r()[i], 0 }; break;
        case GC_type::MAT_COMPLEX: val = _m_c()[i]; break;
        case GC_type::VECTOR: val = ( *this )( i ).get_complex(); break;

        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_assert(
            false,
            "{} bad data type: `{}' cannot be converted into `vec_complex_type'",
            where,
            to_string( get_type() ) );
      }
      v.emplace_back( val );
    }
  }

  void GenericContainer::copyto_vec_string( vec_string_type & v, string_view const where ) const
  {
    v.clear();
    std::size_t const ne{ get_num_elements() };
    switch ( get_type() )
    {
      case GC_type::STRING:
        v.reserve( ne );
        v.emplace_back( _s() );
        break;
      case GC_type::VEC_STRING:
        v.resize( ne );
        std::copy( _v_s().begin(), _v_s().end(), v.begin() );
        break;
      case GC_type::VECTOR:
        v.reserve( ne );
        for ( std::size_t i{ 0 }; i < ne; ++i )
        {
          GenericContainer const & gc = get_gc_at( i, where );
          v.emplace_back( gc.get_string( where ) );
        }
        break;
      case GC_type::NOTYPE:
      case GC_type::BOOL:
      case GC_type::INTEGER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::MAT_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::MAT_LONG:
      case GC_type::VEC_REAL:
      case GC_type::MAT_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::MAT_COMPLEX:
      case GC_type::POINTER:
      case GC_type::VEC_POINTER:
      case GC_type::MAP:
        GC_assert(
          false,
          "{} bad data type: `{}' cannot be converted into `vec_string_type'",
          where,
          to_string( get_type() ) );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_int( mat_int_type & m, string_view const where ) const
  {
    m.clear();
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _b() ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _i();
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ &_v_b() };
        m.resize( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ &_v_i() };
        m.resize( static_cast<std::size_t>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ &_m_i() };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        if ( auto nc{ static_cast<std::size_t>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( std::size_t j{ 1 }; j < nc; ++j )
          {
            GC_assert(
              v[j].get_num_elements() == nr,
              "{} copyto_mat_int() cannot promote vector of size {} to a column of mat_int_type of size {} x {}",
              where,
              v[j].get_num_elements(),
              nr,
              nc );
          }
          m.resize( nr, nc );
          for ( std::size_t j{ 0 }; j < nc; ++j )
          {
            vec_int_type vj;
            v[j].copyto_vec_int( vj, where );
            for ( std::size_t i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_assert( false, "{} copyto_mat_int() cannot promote vector of size {} to mat_int_type", where, v.size() );
        [[fallthrough]];
      }
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::LONG:
      case GC_type::COMPLEX:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP:
        GC_assert( false, "{} copyto_mat_int() cannot promote {} to mat_int_type", where, get_type_name() );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_long( mat_long_type & m, string_view const where ) const
  {
    m.clear();
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _b() ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<long_type>( _i() );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _l();
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ &_v_b() };
        m.resize( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ &_v_i() };
        m.resize( static_cast<std::size_t>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ &_v_l() };
        m.resize( static_cast<std::size_t>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ &_m_i() };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ &_m_l() };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        if ( auto nc{ static_cast<std::size_t>( v.size() ) }; nc > 0 )
        {
          std::size_t nr{ v[0].get_num_elements() };
          for ( std::size_t j{ 1 }; j < nc; ++j )
          {
            GC_assert(
              v[j].get_num_elements() == nr,
              "{} copyto_mat_long() cannot promote vector of size {} to a column of mat_long_type of size {} x {}",
              where,
              v[j].get_num_elements(),
              nr,
              nc );
          }
          m.resize( nr, nc );
          for ( std::size_t j{ 0 }; j < nc; ++j )
          {
            vec_long_type vj;
            v[j].copyto_vec_long( vj, where );
            for ( std::size_t i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_assert( false, "{} copyto_mat_long() cannot promote vector of size {} to mat_long_type", where, v.size() );
        [[fallthrough]];
      }
      case GC_type::MAT_REAL:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP:
        GC_assert( false, "{} copyto_mat_long() cannot promote {} to mat_long_type", where, get_type_name() );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_real( mat_real_type & m, string_view const where ) const
  {
    m.clear();
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _b() ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<real_type>( _i() );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<real_type>( _l() );
      }
      break;
      case GC_type::REAL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _r();
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ &_v_b() };
        m.resize( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ &_v_i() };
        m.resize( static_cast<std::size_t>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ &_v_l() };
        m.resize( static_cast<std::size_t>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ &_v_r() };
        m.resize( static_cast<std::size_t>( v_r->size() ), 1 );
        std::copy_n( v_r->data(), v_r->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ &_m_i() };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ &_m_l() };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const * m_r{ &_m_r() };
        m.resize( m_r->num_rows(), m_r->num_cols() );
        std::copy_n( m_r->data(), m_r->num_rows() * m_r->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        if ( auto nc{ static_cast<std::size_t>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( std::size_t j{ 1 }; j < nc; ++j )
          {
            GC_assert(
              v[j].get_num_elements() == nr,
              "{} copyto_mat_real() cannot promote vector of size {} to a column of mat_real_type of size {} x {}",
              where,
              v[j].get_num_elements(),
              nr,
              nc );
          }
          m.resize( nr, nc );
          for ( std::size_t j{ 0 }; j < nc; ++j )
          {
            vec_real_type vj;
            v[j].copyto_vec_real( vj, where );
            for ( std::size_t i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_assert( false, "{} copyto_mat_real() cannot promote vector of size {} to mat_real_type", where, v.size() );
        [[fallthrough]];
      }
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP:
        GC_assert( false, "{} copyto_mat_real() cannot promote {} to mat_real_type", where, get_type_name() );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_complex( mat_complex_type & m, string_view const where ) const
  {
    m.clear();
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( _b() ? 1.0 : 0.0, 0.0 );
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( static_cast<real_type>( _i() ), 0.0 );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( static_cast<real_type>( _l() ), 0.0 );
      }
      break;
      case GC_type::REAL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( _r(), 0.0 );
      }
      break;
      case GC_type::COMPLEX:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = _c();
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ &_v_b() };
        m.resize( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ &_v_i() };
        m.resize( static_cast<std::size_t>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ &_v_l() };
        m.resize( static_cast<std::size_t>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ &_v_r() };
        m.resize( static_cast<std::size_t>( v_r->size() ), 1 );
        std::copy_n( v_r->data(), v_r->size(), m.data() );
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        vec_complex_type const * v_c{ &_v_c() };
        m.resize( static_cast<std::size_t>( v_c->size() ), 1 );
        std::copy_n( v_c->data(), v_c->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ &_m_i() };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ &_m_l() };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const * m_r{ &_m_r() };
        m.resize( m_r->num_rows(), m_r->num_cols() );
        std::copy_n( m_r->data(), m_r->num_rows() * m_r->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_COMPLEX:
      {
        mat_complex_type const * m_c{ &_m_c() };
        m.resize( m_c->num_rows(), m_c->num_cols() );
        std::copy_n( m_c->data(), m_c->num_rows() * m_c->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ _v() };
        if ( auto nc{ static_cast<std::size_t>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( std::size_t j{ 1 }; j < nc; ++j )
          {
            GC_assert(
              v[j].get_num_elements() == nr,
              "{} copyto_mat_complex() cannot promote vector of size {} to a column of mat_complex_type of size {} x {}",
              where,
              v[j].get_num_elements(),
              nr,
              nc );
          }
          m.resize( nr, nc );
          for ( std::size_t j{ 0 }; j < nc; ++j )
          {
            vec_complex_type vj;
            v[j].copyto_vec_complex( vj, where );
            for ( std::size_t i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_assert(
          false,
          "{} copyto_mat_complex() cannot promote vector of size {} to mat_complex_type",
          where,
          v.size() );
        [[fallthrough]];
      }
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAP:
        GC_assert( false, "{} copyto_mat_complex() cannot promote {} to mat_complex_type", where, get_type_name() );
    }
  }

  /*
  //   ____                            _
  //  |  _ \ _ __ ___  _ __ ___   ___ | |_ ___
  //  | |_) | '__/ _ \| '_ ` _ \ / _ \| __/ _ \
  //  |  __/| | | (_) | | | | | | (_) | ||  __/
  //  |_|   |_|  \___/|_| |_| |_|\___/ \__\___|
  */

  GenericContainer const & GenericContainer::promote_to_int()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: set_int( 0 ); break;
      case GC_type::BOOL: set_int( _b() ? 1 : 0 ); break;
      case GC_type::INTEGER: break;
      case GC_type::POINTER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_int() cannot promote {} to int", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_long()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: set_long( 0 ); break;
      case GC_type::BOOL: set_long( _b() ? 1 : 0 ); break;
      case GC_type::INTEGER: set_long( _i() ); break;
      case GC_type::LONG: break;
      case GC_type::POINTER:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_LONG:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_long() cannot promote {} to long", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_real()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: set_real( 0 ); break;
      case GC_type::BOOL: set_real( _b() ? 1 : 0 ); break;
      case GC_type::INTEGER: set_real( static_cast<real_type>( _i() ) ); break;
      case GC_type::LONG: set_real( static_cast<real_type>( _l() ) ); break;
      case GC_type::REAL: break;
      case GC_type::POINTER:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_real() cannot promote {} to real", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_complex()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE: set_complex( 0, 0 ); break;
      case GC_type::BOOL: set_complex( _b() ? 1 : 0, 0 ); break;
      case GC_type::INTEGER: set_complex( static_cast<real_type>( _i() ), 0 ); break;
      case GC_type::LONG: set_complex( static_cast<real_type>( _l() ), 0 ); break;
      case GC_type::REAL: set_complex( static_cast<real_type>( _r() ), 0 ); break;
      case GC_type::COMPLEX: break;
      case GC_type::VEC_POINTER:
      case GC_type::VEC_BOOL:
      case GC_type::VEC_INTEGER:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL: promote_to_vec_complex(); break;
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_complex_type() cannot promote {} to real type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_int()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_vec_int( 1 );
        get_int_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        int_type const tmp{ _b() ? 1 : 0 };
        set_vec_int( 1 );
        get_int_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        int_type const tmp{ _i() };
        set_vec_int( 1 );
        get_int_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_vec_int( v_b->size() );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _v_i()[i] = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER: break;
      case GC_type::POINTER:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_LONG:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_vec_int() cannot promote {} to vec_int_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_long()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_vec_long( 1 );
        get_long_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        long_type const tmp{ _b() ? 1 : 0 };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        long_type const tmp{ static_cast<long_type>( _i() ) };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        long_type const tmp{ _l() };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_vec_long( v_b->size() );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _v_l()[i] = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_vec_long( v_i->size() );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _v_l()[i] = static_cast<int_type>( ( *v_i )[i] );
      }
      break;
      case GC_type::VEC_LONG:  // nothing to do
        break;
      case GC_type::POINTER:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_vec_long() cannot promote {} to vec_long_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_real()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_vec_real( 1 );
        get_real_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( _b() ? 1 : 0 ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( _i() ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( _l() ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ _r() };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_vec_real( v_b->size() );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _v_r()[i] = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_vec_real( v_i->size() );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _v_r()[i] = static_cast<real_type>( ( *v_i )[i] );
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_vec_real( v_l->size() );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i ) _v_r()[i] = static_cast<real_type>( ( *v_l )[i] );
      }
      break;
      case GC_type::VEC_REAL: break;
      case GC_type::POINTER:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_vec_real() cannot promote {} vec_real_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_complex()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_vec_complex( 1 );
        get_complex_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( _b() ? 1 : 0 ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( _i() ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( _l() ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ _r() };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::COMPLEX:
      {
        complex_type const tmp{ _c() };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_vec_complex( v_b->size() );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _v_c()[i] = complex_type( ( *v_b )[i] ? 1 : 0, 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_vec_complex( v_i->size() );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i )
          _v_c()[i] = complex_type( static_cast<real_type>( ( *v_i )[i] ), 0 );
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_vec_complex( v_l->size() );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i )
          _v_c()[i] = complex_type( static_cast<real_type>( ( *v_l )[i] ), 0 );
      }
      break;
      case GC_type::VEC_REAL:
      {
        auto const v_r{ take_box<vec_real_type>() };
        set_vec_complex( v_r->size() );
        for ( std::size_t i{ 0 }; i < v_r->size(); ++i ) _v_c()[i] = complex_type( ( *v_r )[i], 0 );
      }
      break;
      case GC_type::VEC_COMPLEX: break;
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_vec_real() cannot promote {} to vec_complex_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_int()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        int_type const tmp{ _b() ? 1 : 0 };
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        int_type const tmp{ _i() };
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_mat_int( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _m_i()( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_mat_int( static_cast<std::size_t>( v_i->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _m_i()( i, 0 ) = ( *v_i )[i];
      }
      break;
      case GC_type::MAT_INTEGER: break;
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::REAL:
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::LONG:
      case GC_type::COMPLEX:
      case GC_type::VEC_LONG:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP: GC_assert( false, "promote_to_mat_int() cannot promote {} to mat_int_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_long()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        long_type const tmp{ _b() ? 1 : 0 };
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        long_type const tmp{ _i() };
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_mat_long( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _m_l()( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_mat_long( static_cast<std::size_t>( v_i->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _m_l()( i, 0 ) = static_cast<long_type>( ( *v_i )[i] );
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_mat_long( static_cast<std::size_t>( v_l->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i ) _m_l()( i, 0 ) = ( *v_l )[i];
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        auto const m_i{ take_box<mat_int_type>() };
        set_mat_long( m_i->num_rows(), m_i->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_i->size() ); ++i )
          _m_l()[i] = static_cast<long_type>( ( *m_i )[i] );
      }
      break;
      case GC_type::MAT_LONG: break;
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::LONG:
      case GC_type::REAL:
      case GC_type::COMPLEX:
      case GC_type::VEC_REAL:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_mat_long() cannot promote {} to mat_long_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_real()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( _b() ? 1 : 0 ) };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( _i() ) };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ _r() };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_mat_real( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _m_r()( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_mat_real( static_cast<std::size_t>( v_i->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _m_r()( i, 0 ) = static_cast<real_type>( ( *v_i )[i] );
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_mat_real( static_cast<std::size_t>( v_l->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i ) _m_r()( i, 0 ) = static_cast<real_type>( ( *v_l )[i] );
      }
      break;
      case GC_type::VEC_REAL:
      {
        auto const v_r{ take_box<vec_real_type>() };
        set_mat_real( static_cast<std::size_t>( v_r->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_r->size(); ++i ) _m_r()( i, 0 ) = ( *v_r )[i];
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        auto const m_i{ take_box<mat_int_type>() };
        set_mat_real( m_i->num_rows(), m_i->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_i->size() ); ++i )
          _m_r()[i] = static_cast<real_type>( ( *m_i )[i] );
      }
      break;
      case GC_type::MAT_LONG:
      {
        auto const m_l{ take_box<mat_long_type>() };
        set_mat_real( m_l->num_rows(), m_l->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_l->size() ); ++i )
          _m_r()[i] = static_cast<real_type>( ( *m_l )[i] );
      }
      break;
      case GC_type::MAT_REAL: break;
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::LONG:
      case GC_type::COMPLEX:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_mat_real() cannot promote {} to mat_real_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_complex()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( _b() ? 1 : 0 ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( _i() ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( _l() ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ _r() };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_mat_complex( static_cast<std::size_t>( v_b->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _m_c()( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_mat_complex( static_cast<std::size_t>( v_i->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _m_c()( i, 0 ) = static_cast<real_type>( ( *v_i )[i] );
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_mat_complex( static_cast<std::size_t>( v_l->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i ) _m_c()( i, 0 ) = static_cast<real_type>( ( *v_l )[i] );
      }
      break;
      case GC_type::VEC_REAL:
      {
        auto const v_r{ take_box<vec_real_type>() };
        set_mat_complex( static_cast<std::size_t>( v_r->size() ), 1 );
        for ( std::size_t i{ 0 }; i < v_r->size(); ++i ) _m_c()( i, 0 ) = ( *v_r )[i];
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        auto const m_i{ take_box<mat_int_type>() };
        set_mat_complex( m_i->num_rows(), m_i->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_i->size() ); ++i )
          _m_c()[i] = complex_type( static_cast<real_type>( ( *m_i )[i] ), 0 );
      }
      break;
      case GC_type::MAT_LONG:
      {
        auto const m_l{ take_box<mat_long_type>() };
        set_mat_complex( m_l->num_rows(), m_l->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_l->size() ); ++i )
          _m_c()[i] = complex_type( static_cast<real_type>( ( *m_l )[i] ), 0 );
      }
      break;
      case GC_type::MAT_REAL:
      {
        auto const m_r{ take_box<mat_real_type>() };
        set_mat_complex( m_r->num_rows(), m_r->num_cols() );
        for ( std::size_t i{ 0 }; i < static_cast<std::size_t>( m_r->size() ); ++i )
          _m_c()[i] = complex_type( ( *m_r )[i], 0 );
      }
      break;
      case GC_type::MAT_COMPLEX: break;
      case GC_type::POINTER:
      case GC_type::COMPLEX:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_STRING:
      case GC_type::VECTOR:
      case GC_type::MAP:
        GC_assert( false, "promote_to_mat_real() cannot promote {} to mat_complex_type", get_type_name() );
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vector()
  {
    switch ( get_type() )
    {
      case GC_type::NOTYPE:
      {
        set_vector( 1 );
        ( *this )[0].clear();
      }  // set data to no type
      break;
      case GC_type::POINTER:
      {
        pointer_type const tmp{ _p() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::BOOL:
      {
        bool_type const tmp{ _b() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        int_type const tmp{ _i() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::LONG:
      {
        long_type const tmp{ _l() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ _r() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::COMPLEX:
      {
        complex_type const tmp{ _c() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::STRING:
      {
        string_type const tmp{ _s() };
        set_vector( 1 );
        ( *this )[0] = tmp;
      }
      break;
      case GC_type::VEC_POINTER:
      {
        auto const v_p{ take_box<vec_pointer_type>() };
        set_vector( static_cast<std::size_t>( v_p->size() ) );
        for ( std::size_t i{ 0 }; i < v_p->size(); ++i ) _v()[i] = ( *v_p )[i];
      }
      break;
      case GC_type::VEC_BOOL:
      {
        auto const v_b{ take_box<vec_bool_type>() };
        set_vector( static_cast<std::size_t>( v_b->size() ) );
        for ( std::size_t i{ 0 }; i < v_b->size(); ++i ) _v()[i] = ( *v_b )[i];
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        auto const v_i{ take_box<vec_int_type>() };
        set_vector( static_cast<std::size_t>( v_i->size() ) );
        for ( std::size_t i{ 0 }; i < v_i->size(); ++i ) _v()[i] = ( *v_i )[i];
      }
      break;
      case GC_type::VEC_LONG:
      {
        auto const v_l{ take_box<vec_long_type>() };
        set_vector( static_cast<std::size_t>( v_l->size() ) );
        for ( std::size_t i{ 0 }; i < v_l->size(); ++i ) _v()[i] = ( *v_l )[i];
      }
      break;
      case GC_type::VEC_REAL:
      {
        auto const v_r{ take_box<vec_real_type>() };
        set_vector( static_cast<std::size_t>( v_r->size() ) );
        for ( std::size_t i{ 0 }; i < v_r->size(); ++i ) _v()[i] = ( *v_r )[i];
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        auto const v_c{ take_box<vec_complex_type>() };
        set_vector( static_cast<std::size_t>( v_c->size() ) );
        for ( std::size_t i{ 0 }; i < v_c->size(); ++i ) _v()[i] = ( *v_c )[i];
      }
      break;
      case GC_type::VEC_STRING:
      {
        auto const v_s{ take_box<vec_string_type>() };
        set_vector( static_cast<std::size_t>( v_s->size() ) );
        for ( std::size_t i{ 0 }; i < v_s->size(); ++i ) _v()[i] = ( *v_s )[i];
      }
      break;
      case GC_type::VECTOR: break;
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP: GC_assert( false, "promote_to_vector() cannot promote {} to vector_type", get_type_name() );
    }
    return *this;
  }

  // instantate classes
  template class mat_type<int_type>;
  template class mat_type<long_type>;
  template class mat_type<real_type>;
  template class mat_type<complex_type>;

}  // namespace GC_namespace

//
// eof: GenericContainer.cc
//
