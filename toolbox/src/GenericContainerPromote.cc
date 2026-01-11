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

  static bool isInteger( real_type const x )
  {
    using std::round;
    return isZero0( x - round( x ) );
  }

  static bool isUnsigned( real_type const x )
  {
    return isInteger( x ) && x >= 0;
  }

#endif

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_int( vec_int_type & v, string_view const where ) const
  {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve( ne );
    long_type    lval;
    real_type    rval;
    complex_type cval;
    int_type     val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = m_data.b ? 1 : 0; break;
        case GC_type::INTEGER: val = m_data.i; break;
        case GC_type::LONG:
          lval = m_data.l;
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval,
            where << " copyto_vec_int: v[" << i << "] = " << lval << " cannot be converted to `integer'" );
          val = static_cast<int_type>( lval );
          break;
        case GC_type::REAL:
          rval = m_data.r;
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_int: v[" << i << "] = " << rval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = ( *m_data.v_b )[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = ( *m_data.v_i )[i]; break;
        case GC_type::VEC_LONG:
          lval = ( *m_data.v_l )[i];
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval,
            where << " copyto_vec_int: v[" << i << "] = " << lval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = ( *m_data.v_r )[i];
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_int: v[" << i << "] = " << rval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = *m_data.c;
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_int: v[" << i << "] = " << cval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = ( *m_data.v_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_int: v[" << i << "] = " << cval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER: val = ( *m_data.m_i )[i]; break;
        case GC_type::MAT_LONG:
          lval = ( *m_data.m_l )[i];
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval,
            where << " copyto_vec_int: v[" << i << "] = " << lval << " cannot be converted to `integer'" );
          val = static_cast<int_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = ( *m_data.m_r )[i];
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_int: v[" << i << "] = " << rval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = ( *m_data.m_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_int: v[" << i << "] = " << cval << " cannot be converted to `integer'" )
          val = static_cast<int_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_int( "GenericContainer::copyto_vec_int " ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_int_type'" )
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_uint( vec_uint_type & v, string_view const where ) const
  {
    v.clear();
    unsigned ne{ get_num_elements() };
    v.reserve( ne );
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    uint_type    val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = m_data.b ? 1 : 0; break;
        case GC_type::INTEGER:
          ival = m_data.i;
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_uint: value = " << ival << " cannot be converted to `unsigned integer'" );
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::LONG:
          lval = m_data.l;
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval && lval >= 0,
            where << " copyto_vec_uint: v[" << i << "] = " << lval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::REAL:
          rval = m_data.r;
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_uint: v[" << i << "] = " << rval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = ( *m_data.v_b )[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER:
          ival = ( *m_data.v_i )[i];
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_uint: value = " << ival << " cannot be converted to `unsigned integer'" );
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::VEC_LONG:
          lval = ( *m_data.v_l )[i];
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval && lval >= 0,
            where << " copyto_vec_uint: v[" << i << "] = " << lval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = ( *m_data.v_r )[i];
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_uint: v[" << i << "] = " << rval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = *m_data.c;
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_uint: v[" << i << "] = " << cval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = ( *m_data.v_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_int: v[" << i << "] = " << cval << " cannot be converted to `unsigned integer'" );
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER:
          ival = ( *m_data.m_i )[i];
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_uint: value = " << ival << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( ival );
          break;
        case GC_type::MAT_LONG:
          lval = ( *m_data.m_l )[i];
          GC_ASSERT(
            static_cast<int_type>( lval ) == lval && lval >= 0,
            where << " copyto_vec_uint: v[" << i << "] = " << lval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = ( *m_data.m_r )[i];
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_uint: v[" << i << "] = " << rval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = ( *m_data.m_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_uint: v[" << i << "] = " << cval << " cannot be converted to `unsigned integer'" )
          val = static_cast<uint_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_uint( "GenericContainer::copyto_vec_uint " ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_uint_type'" )
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_long( vec_long_type & v, string_view const where ) const
  {
    v.clear();
    unsigned ne{ get_num_elements() };
    v.reserve( ne );
    real_type    rval;
    complex_type cval;
    long_type    val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = m_data.b ? 1 : 0; break;
        case GC_type::INTEGER: val = static_cast<long_type>( m_data.i ); break;
        case GC_type::LONG: val = m_data.l; break;
        case GC_type::REAL:
          rval = m_data.r;
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_long: v[" << i << "] = " << rval << " cannot be converted to `long'" )
          val = static_cast<long_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = ( *m_data.v_b )[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = static_cast<long_type>( ( *m_data.v_i )[i] ); break;
        case GC_type::VEC_LONG: val = ( *m_data.v_l )[i]; break;
        case GC_type::VEC_REAL:
          rval = ( *m_data.v_r )[i];
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_long: v[" << i << "] = " << rval << " cannot be converted to `long'" )
          val = static_cast<long_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = *m_data.c;
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_long: v[" << i << "] = " << cval << " cannot be converted to `long'" )
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = ( *m_data.v_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_long: v[" << i << "] = " << cval << " cannot be converted to `long'" )
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER: val = static_cast<long_type>( ( *m_data.m_i )[i] ); break;
        case GC_type::MAT_LONG: val = ( *m_data.m_l )[i]; break;
        case GC_type::MAT_REAL:
          rval = ( *m_data.m_r )[i];
          GC_ASSERT(
            isInteger( rval ),
            where << " copyto_vec_long: v[" << i << "] = " << rval << " cannot be converted to `long'" )
          val = static_cast<long_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = ( *m_data.m_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isInteger( cval.real() ),
            where << " copyto_vec_long: v[" << i << "] = " << cval << " cannot be converted to `long'" )
          val = static_cast<long_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_long( "copyto_vec_long" ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_long_type'" )
      }
      v.emplace_back( val );
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_vec_ulong( vec_ulong_type & v, string_view const where ) const
  {
    v.clear();
    unsigned ne{ get_num_elements() };
    v.reserve( ne );
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    ulong_type   val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = m_data.b ? 1 : 0; break;
        case GC_type::INTEGER:
          ival = m_data.i;
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_ulong: value = " << ival << " cannot be converted to `unsigned long'" );
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::LONG:
          lval = m_data.l;
          GC_ASSERT(
            lval >= 0,
            where << " copyto_vec_ulong: v[" << i << "] = " << lval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::REAL:
          rval = m_data.r;
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_ulong: v[" << i << "] = " << rval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::VEC_BOOL: val = ( *m_data.v_b )[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER:
          ival = ( *m_data.v_i )[i];
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_ulong: value = " << ival << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::VEC_LONG:
          lval = ( *m_data.v_l )[i];
          GC_ASSERT(
            lval >= 0,
            where << " copyto_vec_ulong: v[" << i << "] = " << lval << " cannot be converted to `unsigned long'" );
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::VEC_REAL:
          rval = ( *m_data.v_r )[i];
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_ulong: v[" << i << "] = " << rval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::COMPLEX:
          cval = *m_data.c;
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_ulong: v[" << i << "] = " << cval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::VEC_COMPLEX:
          cval = ( *m_data.v_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_ulong: v[" << i << "] = " << cval << " cannot be converted to `unsigned long'" );
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::MAT_INTEGER:
          ival = ( *m_data.m_i )[i];
          GC_ASSERT(
            ival >= 0,
            where << " copyto_vec_ulong: value = " << ival << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( ival );
          break;
        case GC_type::MAT_LONG:
          lval = ( *m_data.m_l )[i];
          GC_ASSERT(
            lval >= 0,
            where << " copyto_vec_ulong: v[" << i << "] = " << lval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( lval );
          break;
        case GC_type::MAT_REAL:
          rval = ( *m_data.m_r )[i];
          GC_ASSERT(
            isUnsigned( rval ),
            where << " copyto_vec_ulong: v[" << i << "] = " << rval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( rval );
          break;
        case GC_type::MAT_COMPLEX:
          cval = ( *m_data.m_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ) && isUnsigned( cval.real() ),
            where << " copyto_vec_ulong: v[" << i << "] = " << cval << " cannot be converted to `unsigned long'" )
          val = static_cast<ulong_type>( cval.real() );
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_as_ulong( "copyto_vec_ulong" ); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_ulong_type'" )
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
    unsigned ne{ get_num_elements() };
    v.reserve( ne );
    complex_type cval;
    real_type    val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = m_data.b ? 1 : 0; break;
        case GC_type::INTEGER: val = static_cast<real_type>( m_data.i ); break;
        case GC_type::LONG: val = static_cast<real_type>( m_data.l ); break;
        case GC_type::REAL: val = m_data.r; break;
        case GC_type::VEC_BOOL: val = ( *m_data.v_b )[i] ? 1 : 0; break;
        case GC_type::VEC_INTEGER: val = static_cast<real_type>( ( *m_data.v_i )[i] ); break;
        case GC_type::VEC_LONG: val = static_cast<real_type>( ( *m_data.v_l )[i] ); break;
        case GC_type::VEC_REAL: val = ( *m_data.v_r )[i]; break;
        case GC_type::COMPLEX:
          cval = *m_data.c;
          GC_ASSERT(
            isZero0( cval.imag() ),
            where << " copyto_vec_real: v[" << i << "] = " << cval << " cannot be converted to `real_type'" )
          val = cval.real();
          break;
        case GC_type::VEC_COMPLEX:
          cval = ( *m_data.v_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ),
            where << " copyto_vec_real: v[" << i << "] = " << cval << " cannot be converted to `real_type'" )
          val = cval.real();
          break;
        case GC_type::MAT_INTEGER: val = static_cast<real_type>( ( *m_data.m_i )[i] ); break;
        case GC_type::MAT_LONG: val = static_cast<real_type>( ( *m_data.m_l )[i] ); break;
        case GC_type::MAT_REAL: val = ( *m_data.m_r )[i]; break;
        case GC_type::MAT_COMPLEX:
          cval = ( *m_data.m_c )[i];
          GC_ASSERT(
            isZero0( cval.imag() ),
            where << " copyto_vec_real: v[" << i << "] = " << cval << " cannot be converted to `real_type'" )
          val = cval.real();
          break;
        case GC_type::VECTOR: val = ( *this )( i ).get_number(); break;
        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_real_type'" )
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
    unsigned const ne{ get_num_elements() };
    v.reserve( ne );
    complex_type val{ 0 };
    v.reserve( ne );
    for ( unsigned i{ 0 }; i < ne; ++i )
    {
      switch ( m_data_type )
      {
        case GC_type::BOOL: val = { static_cast<real_type>( m_data.b ? 1 : 0 ), 0 }; break;
        case GC_type::INTEGER: val = { static_cast<real_type>( m_data.i ), 0 }; break;
        case GC_type::LONG: val = { static_cast<real_type>( m_data.l ), 0 }; break;
        case GC_type::REAL: val = { m_data.r, 0 }; break;
        case GC_type::VEC_BOOL: val = { static_cast<real_type>( ( *m_data.v_b )[i] ? 1 : 0 ), 0 }; break;
        case GC_type::VEC_INTEGER: val = { static_cast<real_type>( ( *m_data.v_i )[i] ), 0 }; break;
        case GC_type::VEC_LONG: val = { static_cast<real_type>( ( *m_data.v_l )[i] ), 0 }; break;
        case GC_type::VEC_REAL: val = { ( *m_data.v_r )[i], 0 }; break;
        case GC_type::COMPLEX: val = *m_data.c; break;
        case GC_type::VEC_COMPLEX: val = ( *m_data.v_c )[i]; break;
        case GC_type::MAT_INTEGER: val = { static_cast<real_type>( ( *m_data.m_i )[i] ), 0 }; break;
        case GC_type::MAT_LONG: val = { static_cast<real_type>( ( *m_data.m_l )[i] ), 0 }; break;
        case GC_type::MAT_REAL: val = { ( *m_data.m_r )[i], 0 }; break;
        case GC_type::MAT_COMPLEX: val = ( *m_data.m_c )[i]; break;
        case GC_type::VECTOR: val = ( *this )( i ).get_complex(); break;

        case GC_type::NOTYPE:
        case GC_type::POINTER:
        case GC_type::STRING:
        case GC_type::VEC_POINTER:
        case GC_type::VEC_STRING:
        case GC_type::MAP:
          GC_DO_ERROR(
            where << " bad data type: `" << to_string( m_data_type )
                  << "' cannot be converted into `vec_complex_type'" )
      }
      v.emplace_back( val );
    }
  }

  void GenericContainer::copyto_vec_string( vec_string_type & v, string_view const where ) const
  {
    v.clear();
    unsigned const ne{ get_num_elements() };
    switch ( m_data_type )
    {
      case GC_type::STRING:
        v.reserve( ne );
        v.emplace_back( *m_data.s );
        break;
      case GC_type::VEC_STRING:
        v.resize( ne );
        std::copy( m_data.v_s->begin(), m_data.v_s->end(), v.begin() );
        break;
      case GC_type::VECTOR:
        v.reserve( ne );
        for ( unsigned i{ 0 }; i < ne; ++i )
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
        GC_DO_ERROR(
          where << " bad data type: `" << to_string( m_data_type ) << "' cannot be converted into `vec_string_type'" )
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_int( mat_int_type & m, string_view const where ) const
  {
    m.clear();
    switch ( m_data_type )
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
        m( 0, 0 ) = m_data.b ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = m_data.i;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m.resize( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m.resize( static_cast<unsigned>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ *m_data.v };
        if ( auto nc{ static_cast<unsigned>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( unsigned j{ 1 }; j < nc; ++j )
          {
            GC_ASSERT(
              v[j].get_num_elements() == nr,
              where << " copyto_mat_int() cannot promote vector of size " << v[j].get_num_elements()
                    << " to a column of mat_int_type of size " << nr << " x " << nc )
          }
          m.resize( nr, nc );
          for ( unsigned j{ 0 }; j < nc; ++j )
          {
            vec_int_type vj;
            v[j].copyto_vec_int( vj, where );
            for ( unsigned i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_DO_ERROR( where << " copyto_mat_int() cannot promote vector of size " << v.size() << " to mat_int_type" )
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
        GC_DO_ERROR( where << " copyto_mat_int() cannot promote " << get_type_name() << " to mat_int_type" )
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_long( mat_long_type & m, string_view const where ) const
  {
    m.clear();
    switch ( m_data_type )
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
        m( 0, 0 ) = m_data.b ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<long_type>( m_data.i );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = m_data.l;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m.resize( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m.resize( static_cast<unsigned>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m.resize( static_cast<unsigned>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ m_data.m_l };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ *m_data.v };
        if ( auto nc{ static_cast<unsigned>( v.size() ) }; nc > 0 )
        {
          unsigned nr{ v[0].get_num_elements() };
          for ( unsigned j{ 1 }; j < nc; ++j )
          {
            GC_ASSERT(
              v[j].get_num_elements() == nr,
              where << " copyto_mat_long() cannot promote vector of size " << v[j].get_num_elements()
                    << " to a column of mat_long_type of size " << nr << " x " << nc )
          }
          m.resize( nr, nc );
          for ( unsigned j{ 0 }; j < nc; ++j )
          {
            vec_long_type vj;
            v[j].copyto_vec_long( vj, where );
            for ( unsigned i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_DO_ERROR( where << " copyto_mat_long() cannot promote vector of size " << v.size() << " to mat_long_type" )
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
        GC_DO_ERROR( where << "copyto_mat_long() cannot promote " << get_type_name() << " to mat_long_type" )
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_real( mat_real_type & m, string_view const where ) const
  {
    m.clear();
    switch ( m_data_type )
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
        m( 0, 0 ) = m_data.b ? 1 : 0;
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<real_type>( m_data.i );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = static_cast<real_type>( m_data.l );
      }
      break;
      case GC_type::REAL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = m_data.r;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m.resize( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m.resize( static_cast<unsigned>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m.resize( static_cast<unsigned>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m.resize( static_cast<unsigned>( v_r->size() ), 1 );
        std::copy_n( v_r->data(), v_r->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ m_data.m_l };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const * m_r{ m_data.m_r };
        m.resize( m_r->num_rows(), m_r->num_cols() );
        std::copy_n( m_r->data(), m_r->num_rows() * m_r->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ *m_data.v };
        if ( auto nc{ static_cast<unsigned>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( unsigned j{ 1 }; j < nc; ++j )
          {
            GC_ASSERT(
              v[j].get_num_elements() == nr,
              where << " copyto_mat_real() cannot promote vector of size " << v[j].get_num_elements()
                    << " to a column of mat_real_type of size " << nr << " x " << nc )
          }
          m.resize( nr, nc );
          for ( unsigned j{ 0 }; j < nc; ++j )
          {
            vec_real_type vj;
            v[j].copyto_vec_real( vj, where );
            for ( unsigned i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_DO_ERROR( where << " copyto_mat_real() cannot promote vector of size " << v.size() << " to mat_real_type" )
      }
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::COMPLEX:
      case GC_type::VEC_COMPLEX:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP:
        GC_DO_ERROR( where << " copyto_mat_real() cannot promote " << get_type_name() << " to mat_real_type" )
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void GenericContainer::copyto_mat_complex( mat_complex_type & m, string_view const where ) const
  {
    m.clear();
    switch ( m_data_type )
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
        m( 0, 0 ) = complex_type( m_data.b ? 1.0 : 0.0, 0.0 );
      }
      break;
      case GC_type::INTEGER:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( static_cast<real_type>( m_data.i ), 0.0 );
      }
      break;
      case GC_type::LONG:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( static_cast<real_type>( m_data.l ), 0.0 );
      }
      break;
      case GC_type::REAL:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = complex_type( m_data.r, 0.0 );
      }
      break;
      case GC_type::COMPLEX:
      {
        m.resize( 1, 1 );
        m( 0, 0 ) = *m_data.c;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m.resize( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) m( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m.resize( static_cast<unsigned>( v_i->size() ), 1 );
        std::copy_n( v_i->data(), v_i->size(), m.data() );
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m.resize( static_cast<unsigned>( v_l->size() ), 1 );
        std::copy_n( v_l->data(), v_l->size(), m.data() );
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m.resize( static_cast<unsigned>( v_r->size() ), 1 );
        std::copy_n( v_r->data(), v_r->size(), m.data() );
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        vec_complex_type const * v_c{ m_data.v_c };
        m.resize( static_cast<unsigned>( v_c->size() ), 1 );
        std::copy_n( v_c->data(), v_c->size(), m.data() );
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m.resize( m_i->num_rows(), m_i->num_cols() );
        std::copy_n( m_i->data(), m_i->num_rows() * m_i->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ m_data.m_l };
        m.resize( m_l->num_rows(), m_l->num_cols() );
        std::copy_n( m_l->data(), m_l->num_rows() * m_l->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const * m_r{ m_data.m_r };
        m.resize( m_r->num_rows(), m_r->num_cols() );
        std::copy_n( m_r->data(), m_r->num_rows() * m_r->num_cols(), m.data() );
      }
      break;
      case GC_type::MAT_COMPLEX:
      {
        mat_complex_type const * m_c{ m_data.m_c };
        m.resize( m_c->num_rows(), m_c->num_cols() );
        std::copy_n( m_c->data(), m_c->num_rows() * m_c->num_cols(), m.data() );
      }
      break;
      case GC_type::VECTOR:
      {
        vector_type const & v{ *m_data.v };
        if ( auto nc{ static_cast<unsigned>( v.size() ) }; nc > 0 )
        {
          auto nr{ v[0].get_num_elements() };
          for ( unsigned j{ 1 }; j < nc; ++j )
          {
            GC_ASSERT(
              v[j].get_num_elements() == nr,
              where << " copyto_mat_complex() cannot promote vector of size " << v[j].get_num_elements()
                    << " to a column of mat_complex_type of size " << nr << " x " << nc )
          }
          m.resize( nr, nc );
          for ( unsigned j{ 0 }; j < nc; ++j )
          {
            vec_complex_type vj;
            v[j].copyto_vec_complex( vj, where );
            for ( unsigned i{ 0 }; i < nr; ++i ) m( i, j ) = vj[i];
          }
          break;  // finito esco.
        }
        GC_DO_ERROR(
          where << " copyto_mat_complex() cannot promote vector of size " << v.size() << " to mat_complex_type" )
      }
      case GC_type::POINTER:
      case GC_type::STRING:
      case GC_type::VEC_POINTER:
      case GC_type::VEC_STRING:
      case GC_type::MAP:
        GC_DO_ERROR( where << "copyto_mat_complex() cannot promote " << get_type_name() << " to mat_complex_type" )
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
    switch ( m_data_type )
    {
      case GC_type::NOTYPE: set_int( 0 ); break;
      case GC_type::BOOL: set_int( m_data.b ? 1 : 0 ); break;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_int() cannot promote " << get_type_name() << " to int" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_long()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE: set_long( 0 ); break;
      case GC_type::BOOL: set_long( m_data.b ? 1 : 0 ); break;
      case GC_type::INTEGER: set_long( m_data.i ); break;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_long() cannot promote " << get_type_name() << " to long" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_real()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE: set_real( 0 ); break;
      case GC_type::BOOL: set_real( m_data.b ? 1 : 0 ); break;
      case GC_type::INTEGER: set_real( static_cast<real_type>( m_data.i ) ); break;
      case GC_type::LONG: set_real( static_cast<real_type>( m_data.l ) ); break;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_real() cannot promote " << get_type_name() << " to real" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_complex()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE: set_complex( 0, 0 ); break;
      case GC_type::BOOL: set_complex( m_data.b ? 1 : 0, 0 ); break;
      case GC_type::INTEGER: set_complex( static_cast<real_type>( m_data.i ), 0 ); break;
      case GC_type::LONG: set_complex( static_cast<real_type>( m_data.l ), 0 ); break;
      case GC_type::REAL: set_complex( static_cast<real_type>( m_data.r ), 0 ); break;
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
        GC_DO_ERROR( "promote_to_complex_type() cannot promote " << get_type_name() << " to real type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_int()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_vec_int( 1 );
        get_int_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        int_type const tmp{ m_data.b ? 1 : 0 };
        set_vec_int( 1 );
        get_int_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        int_type const tmp{ m_data.i };
        set_vec_int( 1 );
        get_int_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_vec_int( v_b->size() );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.v_i )[i] = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_vec_int() cannot promote " << get_type_name() << " to vec_int_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_long()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_vec_long( 1 );
        get_long_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        long_type const tmp{ m_data.b ? 1 : 0 };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        long_type const tmp{ static_cast<long_type>( m_data.i ) };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        long_type const tmp{ m_data.l };
        set_vec_long( 1 );
        get_long_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_vec_long( v_b->size() );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.v_l )[i] = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_vec_long( v_i->size() );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.v_l )[i] = static_cast<int_type>( ( *v_i )[i] );
        delete v_i;
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
        GC_DO_ERROR( "promote_to_vec_long() cannot promote " << get_type_name() << " to vec_long_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_real()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_vec_real( 1 );
        get_real_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( m_data.b ? 1 : 0 ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( m_data.i ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( m_data.l ) };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ m_data.r };
        set_vec_real( 1 );
        get_real_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };  // salva puntatore
        m_data_type = GC_type::NOTYPE;
        set_vec_real( v_b->size() );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.v_r )[i] = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };  // salva puntatore
        m_data_type = GC_type::NOTYPE;
        set_vec_real( v_i->size() );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.v_r )[i] = static_cast<real_type>( ( *v_i )[i] );
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };  // salva puntatore
        m_data_type = GC_type::NOTYPE;
        set_vec_real( v_l->size() );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i ) ( *m_data.v_r )[i] = static_cast<real_type>( ( *v_l )[i] );
        delete v_l;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_vec_real() cannot promote " << get_type_name() << " vec_real_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vec_complex()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_vec_complex( 1 );
        get_complex_at( 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( m_data.b ? 1 : 0 ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( m_data.i ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( m_data.l ) };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ m_data.r };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::COMPLEX:
      {
        complex_type const tmp{ *m_data.c };
        set_vec_complex( 1 );
        get_complex_at( 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_vec_complex( v_b->size() );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.v_c )[i] = complex_type( ( *v_b )[i] ? 1 : 0, 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_vec_complex( v_i->size() );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i )
          ( *m_data.v_c )[i] = complex_type( static_cast<real_type>( ( *v_i )[i] ), 0 );
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m_data_type = GC_type::NOTYPE;
        set_vec_complex( v_l->size() );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i )
          ( *m_data.v_c )[i] = complex_type( static_cast<real_type>( ( *v_l )[i] ), 0 );
        delete v_l;
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m_data_type = GC_type::NOTYPE;
        set_vec_complex( v_r->size() );
        for ( unsigned i{ 0 }; i < v_r->size(); ++i ) ( *m_data.v_c )[i] = complex_type( ( *v_r )[i], 0 );
        delete v_r;
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
        GC_DO_ERROR( "promote_to_vec_real() cannot promote " << get_type_name() << " to vec_complex_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_int()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        int_type const tmp{ m_data.b ? 1 : 0 };
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        int_type const tmp{ m_data.i };
        set_mat_int( 1, 1 );
        get_int_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_mat_int( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_int( static_cast<unsigned>( v_i->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.m_r )( i, 0 ) = static_cast<real_type>( ( *v_i )[i] );
        delete v_i;
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
      case GC_type::MAP: GC_DO_ERROR( "promote_to_mat_int() cannot promote " << get_type_name() << " to mat_int_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_long()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        long_type const tmp{ m_data.b ? 1 : 0 };
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        long_type const tmp{ m_data.i };
        set_mat_long( 1, 1 );
        get_long_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_long( static_cast<unsigned>( v_i->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.m_l )( i, 0 ) = static_cast<long_type>( ( *v_i )[i] );
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m_data_type = GC_type::NOTYPE;
        set_mat_long( static_cast<unsigned>( v_l->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i ) ( *m_data.m_l )( i, 0 ) = ( *v_l )[i];
        delete v_l;
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_long( m_i->num_rows(), m_i->num_cols() );
        for ( unsigned i{ 0 }; i < m_i->size(); ++i ) ( *m_data.m_l )[i] = static_cast<long_type>( ( *m_i )[i] );
        delete m_i;
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
        GC_DO_ERROR( "promote_to_mat_long() cannot promote " << get_type_name() << " to mat_long_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_real()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( m_data.b ? 1 : 0 ) };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( m_data.i ) };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ m_data.r };
        set_mat_real( 1, 1 );
        get_real_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( static_cast<unsigned>( v_i->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.m_r )( i, 0 ) = static_cast<real_type>( ( *v_i )[i] );
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( static_cast<unsigned>( v_l->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i ) ( *m_data.m_r )( i, 0 ) = static_cast<real_type>( ( *v_l )[i] );
        delete v_l;
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( static_cast<unsigned>( v_r->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_r->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( *v_r )[i];
        delete v_r;
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( m_i->num_rows(), m_i->num_cols() );
        for ( unsigned i{ 0 }; i < m_i->size(); ++i ) ( *m_data.m_r )[i] = static_cast<real_type>( ( *m_i )[i] );
        delete m_i;
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ m_data.m_l };
        m_data_type = GC_type::NOTYPE;
        set_mat_real( m_l->num_rows(), m_l->num_cols() );
        for ( unsigned i{ 0 }; i < m_l->size(); ++i ) ( *m_data.m_r )[i] = static_cast<real_type>( ( *m_l )[i] );
        delete m_l;
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
        GC_DO_ERROR( "promote_to_mat_real() cannot promote " << get_type_name() << " to mat_real_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_mat_complex()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = 0;
      }
      break;
      case GC_type::BOOL:
      {
        real_type const tmp{ static_cast<real_type>( m_data.b ? 1 : 0 ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::INTEGER:
      {
        real_type const tmp{ static_cast<real_type>( m_data.i ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::LONG:
      {
        real_type const tmp{ static_cast<real_type>( m_data.l ) };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::REAL:
      {
        real_type const tmp{ m_data.r };
        set_mat_complex( 1, 1 );
        get_complex_at( 0, 0 ) = tmp;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( static_cast<unsigned>( v_b->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( ( *v_b )[i] ? 1 : 0 );
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( static_cast<unsigned>( v_i->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.m_r )( i, 0 ) = static_cast<real_type>( ( *v_i )[i] );
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( static_cast<unsigned>( v_l->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i ) ( *m_data.m_r )( i, 0 ) = static_cast<real_type>( ( *v_l )[i] );
        delete v_l;
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( static_cast<unsigned>( v_r->size() ), 1 );
        for ( unsigned i{ 0 }; i < v_r->size(); ++i ) ( *m_data.m_r )( i, 0 ) = ( *v_r )[i];
        delete v_r;
      }
      break;
      case GC_type::MAT_INTEGER:
      {
        mat_int_type const * m_i{ m_data.m_i };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( m_i->num_rows(), m_i->num_cols() );
        for ( unsigned i{ 0 }; i < m_i->size(); ++i )
          ( *m_data.m_c )[i] = complex_type( static_cast<real_type>( ( *m_i )[i] ), 0 );
        delete m_i;
      }
      break;
      case GC_type::MAT_LONG:
      {
        mat_long_type const * m_l{ m_data.m_l };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( m_l->num_rows(), m_l->num_cols() );
        for ( unsigned i{ 0 }; i < m_l->size(); ++i )
          ( *m_data.m_c )[i] = complex_type( static_cast<real_type>( ( *m_l )[i] ), 0 );
        delete m_l;
      }
      break;
      case GC_type::MAT_REAL:
      {
        mat_real_type const * m_r{ m_data.m_r };
        m_data_type = GC_type::NOTYPE;
        set_mat_complex( m_r->num_rows(), m_r->num_cols() );
        for ( unsigned i{ 0 }; i < m_r->size(); ++i ) ( *m_data.m_c )[i] = complex_type( ( *m_r )[i], 0 );
        delete m_r;
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
        GC_DO_ERROR( "promote_to_mat_real() cannot promote " << get_type_name() << " to mat_complex_type" )
    }
    return *this;
  }

  GenericContainer const & GenericContainer::promote_to_vector()
  {
    switch ( m_data_type )
    {
      case GC_type::NOTYPE:
      {
        set_vector( 1 );
        ( *this )[0].clear();
      }  // set data to no type
      break;
      case GC_type::POINTER:
      {
        set_vector( 1 );
        ( *this )[0] = m_data.p;
      }
      break;
      case GC_type::BOOL:
      {
        set_vector( 1 );
        ( *this )[0] = m_data.b;
      }
      break;
      case GC_type::INTEGER:
      {
        set_vector( 1 );
        ( *this )[0] = m_data.i;
      }
      break;
      case GC_type::LONG:
      {
        set_vector( 1 );
        ( *this )[0] = m_data.l;
      }
      break;
      case GC_type::REAL:
      {
        set_vector( 1 );
        ( *this )[0] = m_data.r;
      }
      break;
      case GC_type::COMPLEX:
      {
        set_vector( 1 );
        ( *this )[0] = *m_data.c;
      }
      break;
      case GC_type::STRING:
      {
        set_vector( 1 );
        ( *this )[0] = *m_data.s;
      }
      break;
      case GC_type::VEC_POINTER:
      {
        vec_pointer_type const * v_p{ m_data.v_p };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_p->size() ) );
        for ( unsigned i{ 0 }; i < v_p->size(); ++i ) ( *m_data.v )[i] = ( *v_p )[i];
        delete v_p;
      }
      break;
      case GC_type::VEC_BOOL:
      {
        vec_bool_type const * v_b{ m_data.v_b };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_b->size() ) );
        for ( unsigned i{ 0 }; i < v_b->size(); ++i ) ( *m_data.v )[i] = ( *v_b )[i];
        delete v_b;
      }
      break;
      case GC_type::VEC_INTEGER:
      {
        vec_int_type const * v_i{ m_data.v_i };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_i->size() ) );
        for ( unsigned i{ 0 }; i < v_i->size(); ++i ) ( *m_data.v )[i] = ( *v_i )[i];
        delete v_i;
      }
      break;
      case GC_type::VEC_LONG:
      {
        vec_long_type const * v_l{ m_data.v_l };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_l->size() ) );
        for ( unsigned i{ 0 }; i < v_l->size(); ++i ) ( *m_data.v )[i] = ( *v_l )[i];
        delete v_l;
      }
      break;
      case GC_type::VEC_REAL:
      {
        vec_real_type const * v_r{ m_data.v_r };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_r->size() ) );
        for ( unsigned i{ 0 }; i < v_r->size(); ++i ) ( *m_data.v )[i] = ( *v_r )[i];
        delete v_r;
      }
      break;
      case GC_type::VEC_COMPLEX:
      {
        vec_complex_type const * v_c{ m_data.v_c };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_c->size() ) );
        for ( unsigned i{ 0 }; i < v_c->size(); ++i ) ( *m_data.v )[i] = ( *v_c )[i];
        delete v_c;
      }
      break;
      case GC_type::VEC_STRING:
      {
        vec_string_type const * v_s{ m_data.v_s };
        m_data_type = GC_type::NOTYPE;
        set_vector( static_cast<unsigned>( v_s->size() ) );
        for ( unsigned i{ 0 }; i < v_s->size(); ++i ) ( *m_data.v )[i] = ( *v_s )[i];
        delete v_s;
      }
      break;
      case GC_type::VECTOR: break;
      case GC_type::MAT_INTEGER:
      case GC_type::MAT_LONG:
      case GC_type::MAT_REAL:
      case GC_type::MAT_COMPLEX:
      case GC_type::MAP: GC_DO_ERROR( "promote_to_vector() cannot promote " << get_type_name() << " to vector_type" )
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
