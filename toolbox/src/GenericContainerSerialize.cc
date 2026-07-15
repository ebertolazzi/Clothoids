/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2021                                                      |
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
// file: GenericContainerSerialize.cc
//

#ifdef _MSC_VER
#pragma warning( disable : 4661 )
#endif

#include "GenericContainer/GenericContainer.hh"

#include <limits>
#include "GenericContainer/GenericContainerLibs.hh"
#include <cstring>
#include <cstdint>

namespace GC_namespace
{

  using std::memcpy;

  static constexpr char serialize_msg[]{
    "GenericContainer::serialize, memory exausted use method mem_size() to estimate memory requirement"
  };

  static constexpr char deserialize_msg[]{
    "GenericContainer::de_serialize, buffer exhausted or corrupted input"
  };

  static int32_t int8_to_buffer( int8_t in, uint8_t * buffer, int32_t const available )
  {
    GC_assert( available >= 1, serialize_msg );
    buffer[0] = static_cast<uint8_t>( in );
    return sizeof( int8_t );
  }

  static int32_t uint32_to_buffer( uint32_t in, uint8_t * buffer, int32_t const available )
  {
    GC_assert( available >= 4, serialize_msg );
    buffer[0] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[1] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[2] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[3] = static_cast<uint8_t>( in & 0xFF );
    return sizeof( uint32_t );
  }

  static int32_t int32_to_buffer( int32_t in, uint8_t * buffer, int32_t const available )
  {
    uint32_t tmp;
    memcpy( &tmp, &in, sizeof( tmp ) );
    return uint32_to_buffer( tmp, buffer, available );
  }

  static int32_t uint64_to_buffer( uint64_t in, uint8_t * buffer, int32_t const available )
  {
    GC_assert( available >= 8, serialize_msg );
    buffer[0] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[1] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[2] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[3] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[4] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[5] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[6] = static_cast<uint8_t>( in & 0xFF );
    in >>= 8;
    buffer[7] = static_cast<uint8_t>( in & 0xFF );
    return sizeof( int64_t );
  }

  static int32_t int64_to_buffer( int64_t in, uint8_t * buffer, int32_t const available )
  {
    uint64_t tmp;
    memcpy( &tmp, &in, sizeof( tmp ) );
    return uint64_to_buffer( tmp, buffer, available );
  }

  static int32_t double_to_buffer( double const in, uint8_t * buffer, int32_t const available )
  {
    uint64_t tmp;
    memcpy( &tmp, &in, sizeof( tmp ) );
    uint64_to_buffer( tmp, buffer, available );
    return sizeof( double );
  }

  /* ---------------------------------------------------------------------------- */

  static int32_t buffer_to_uint8( uint8_t const * buffer, int32_t const available, uint8_t * out )
  {
    GC_assert( available >= 1, deserialize_msg );
    *out = buffer[0];
    return sizeof( uint8_t );
  }

  static int32_t buffer_to_uint32( uint8_t const * buffer, int32_t const available, uint32_t * out )
  {
    GC_assert( available >= 4, deserialize_msg );
    uint32_t const tmp0{ buffer[0] };
    uint32_t const tmp1{ buffer[1] };
    uint32_t const tmp2{ buffer[2] };
    uint32_t const tmp3{ buffer[3] };
    *out = tmp0 | ( tmp1 << 8 ) | ( tmp2 << 16 ) | ( tmp3 << 24 );
    return sizeof( int32_t );
  }

  static int32_t buffer_to_int32( uint8_t const * buffer, int32_t const available, int32_t * out )
  {
    uint32_t tmp;
    int32_t  nb{ buffer_to_uint32( buffer, available, &tmp ) };
    memcpy( out, &tmp, sizeof( tmp ) );
    return nb;
  }

  static int32_t buffer_to_uint64( uint8_t const * buffer, int32_t const available, uint64_t * out )
  {
    GC_assert( available >= 8, deserialize_msg );
    uint64_t const tmp0{ buffer[0] };
    uint64_t const tmp1{ buffer[1] };
    uint64_t const tmp2{ buffer[2] };
    uint64_t const tmp3{ buffer[3] };
    uint64_t const tmp4{ buffer[4] };
    uint64_t const tmp5{ buffer[5] };
    uint64_t const tmp6{ buffer[6] };
    uint64_t const tmp7{ buffer[7] };
    *out = tmp0 | ( tmp1 << 8 ) | ( tmp2 << 16 ) | ( tmp3 << 24 ) | ( tmp4 << 32 ) | ( tmp5 << 40 ) | ( tmp6 << 48 ) |
           ( tmp7 << 56 );
    return sizeof( uint64_t );
  }

  static int32_t buffer_to_int64( uint8_t const * buffer, int32_t const available, int64_t * out )
  {
    uint64_t tmp;
    int32_t  nb{ buffer_to_uint64( buffer, available, &tmp ) };
    memcpy( out, &tmp, sizeof( tmp ) );
    return nb;
  }

  static int32_t buffer_to_double( uint8_t const * buffer, int32_t const available, double * out )
  {
    uint64_t tmp;
    int32_t  nb{ buffer_to_uint64( buffer, available, &tmp ) };
    memcpy( out, &tmp, sizeof( tmp ) );
    return nb;
  }

  //!
  //! Serialized size computed in 64 bits through the public accessors, so a
  //! container larger than the int32 wire-format limit is detected instead
  //! of silently truncated.
  //!
  static uint64_t gc_mem_size( GenericContainer const & gc )
  {
    constexpr uint64_t header_size{ sizeof( int32_t ) };
    constexpr uint64_t len_size{ sizeof( int32_t ) };
    constexpr uint64_t ptr_size{ 8 };
    uint64_t           res{ 0 };
    switch ( gc.get_type() )
    {
      case GC_type::NOTYPE: res = header_size; break;
      case GC_type::BOOL: res = header_size + 1; break;
      case GC_type::INTEGER: res = header_size + sizeof( int_type ); break;
      case GC_type::LONG: res = header_size + sizeof( long_type ); break;
      case GC_type::REAL: res = header_size + sizeof( real_type ); break;
      case GC_type::POINTER: res = header_size + ptr_size; break;
      case GC_type::STRING: res = header_size + gc.get_string().length() + 5; break;
      case GC_type::COMPLEX: res = header_size + sizeof( complex_type ); break;
      case GC_type::VEC_POINTER: res = header_size + len_size + ptr_size * gc.get_vec_pointer().size(); break;
      case GC_type::VEC_BOOL: res = header_size + len_size + gc.get_vec_bool().size(); break;
      case GC_type::VEC_INTEGER: res = header_size + len_size + sizeof( int_type ) * gc.get_vec_int().size(); break;
      case GC_type::VEC_LONG: res = header_size + len_size + sizeof( long_type ) * gc.get_vec_long().size(); break;
      case GC_type::VEC_REAL: res = header_size + len_size + sizeof( real_type ) * gc.get_vec_real().size(); break;
      case GC_type::VEC_COMPLEX:
        res = header_size + len_size + sizeof( complex_type ) * gc.get_vec_complex().size();
        break;
      case GC_type::MAT_INTEGER:
        res = header_size + 2 * len_size + sizeof( int_type ) * uint64_t( gc.get_mat_int().size() );
        break;
      case GC_type::MAT_LONG:
        res = header_size + 2 * len_size + sizeof( long_type ) * uint64_t( gc.get_mat_long().size() );
        break;
      case GC_type::MAT_REAL:
        res = header_size + 2 * len_size + sizeof( real_type ) * uint64_t( gc.get_mat_real().size() );
        break;
      case GC_type::MAT_COMPLEX:
        res = header_size + 2 * len_size + sizeof( complex_type ) * uint64_t( gc.get_mat_complex().size() );
        break;
      case GC_type::VEC_STRING:
        res = header_size + len_size;
        for ( auto & s : gc.get_vec_string() ) res += len_size + s.length() + 1;
        break;
      case GC_type::VECTOR:
        res = header_size + len_size;
        for ( auto & S : gc.get_vector() ) res += gc_mem_size( S );
        break;
      case GC_type::MAP:
        res = header_size + len_size;
        for ( auto & [fst, snd] : gc.get_map() ) res += len_size + fst.length() + 1 + gc_mem_size( snd );
        break;
    }
    return res;
  }

  int32_t GenericContainer::mem_size() const
  {
    uint64_t const res{ gc_mem_size( *this ) };
    GC_assert(
      res <= uint64_t( std::numeric_limits<int32_t>::max() ),
      "GenericContainer::mem_size() serialized size {} exceeds the int32 wire-format limit",
      res
    );
    return static_cast<int32_t>( res );
  }

  int32_t GenericContainer::serialize( std::vector<uint8_t> & buffer ) const
  {
    int32_t const sz{ static_cast<int32_t>( mem_size() ) };
    buffer.resize( static_cast<std::size_t>( sz ) );
    return this->serialize( sz, buffer.data() );
  }

  int32_t GenericContainer::serialize( int32_t buffer_dim, uint8_t * buffer ) const
  {
    int32_t sz, available{ buffer_dim };

    int32_t nb = int32_to_buffer( static_cast<int32_t>( get_type() ), buffer, available );
    available -= nb;
    buffer += nb;

    switch ( get_type() )
    {
      case GC_type::NOTYPE: break;
      case GC_type::BOOL:
        nb = int8_to_buffer( _b() ? 1 : 0, buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::INTEGER:
        nb = int32_to_buffer( _i(), buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::LONG:
        nb = int64_to_buffer( _l(), buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::REAL:
        nb = double_to_buffer( _r(), buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::POINTER:
        nb = int64_to_buffer( reinterpret_cast<int64_t>( _p() ), buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::STRING:
        sz = static_cast<int32_t>( _s().length() + 1 );
        nb = int32_to_buffer( sz, buffer, available );
        buffer += nb;
        available -= nb;
        GC_assert( sz <= available, serialize_msg );
        memcpy( buffer, _s().c_str(), static_cast<size_t>( sz ) );
        buffer += sz;
        available -= sz;
        break;
      case GC_type::COMPLEX:
        nb = double_to_buffer( _c().real(), buffer, available );
        buffer += nb;
        available -= nb;
        nb = double_to_buffer( _c().imag(), buffer, available );
        buffer += nb;
        available -= nb;
        break;
      case GC_type::VEC_POINTER:
        nb = int32_to_buffer( static_cast<int32_t>( _v_p().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto p : _v_p() )
        {
          nb = int64_to_buffer( reinterpret_cast<int64_t>( p ), buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_BOOL:
        nb = int32_to_buffer( static_cast<int32_t>( _v_b().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto b : _v_b() )
        {
          nb = int8_to_buffer( b ? 1 : 0, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_INTEGER:
        nb = int32_to_buffer( static_cast<int32_t>( _v_i().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & i : _v_i() )
        {
          nb = int32_to_buffer( i, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_LONG:
        nb = int32_to_buffer( static_cast<int32_t>( _v_l().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & i : _v_l() )
        {
          nb = int64_to_buffer( i, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_REAL:
        nb = int32_to_buffer( static_cast<int32_t>( _v_r().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & r : _v_r() )
        {
          nb = double_to_buffer( r, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_COMPLEX:
        nb = int32_to_buffer( static_cast<int32_t>( _v_c().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & c : _v_c() )
        {
          nb = double_to_buffer( c.real(), buffer, available );
          buffer += nb;
          available -= nb;
          nb = double_to_buffer( c.imag(), buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::MAT_INTEGER:
        nb = int32_to_buffer( static_cast<int32_t>( _m_i().num_rows() ), buffer, available );
        buffer += nb;
        available -= nb;
        nb = int32_to_buffer( static_cast<int32_t>( _m_i().num_cols() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & i : _m_i() )
        {
          nb = int32_to_buffer( i, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::MAT_LONG:
        nb = int32_to_buffer( static_cast<int32_t>( _m_l().num_rows() ), buffer, available );
        buffer += nb;
        available -= nb;
        nb = int32_to_buffer( static_cast<int32_t>( _m_l().num_cols() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & i : _m_l() )
        {
          nb = int64_to_buffer( i, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::MAT_REAL:
        nb = int32_to_buffer( static_cast<int32_t>( _m_r().num_rows() ), buffer, available );
        buffer += nb;
        available -= nb;
        nb = int32_to_buffer( static_cast<int32_t>( _m_r().num_cols() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & r : _m_r() )
        {
          nb = double_to_buffer( r, buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::MAT_COMPLEX:
        nb = int32_to_buffer( static_cast<int32_t>( _m_c().num_rows() ), buffer, available );
        buffer += nb;
        available -= nb;
        nb = int32_to_buffer( static_cast<int32_t>( _m_c().num_cols() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & c : _m_c() )
        {
          nb = double_to_buffer( c.real(), buffer, available );
          buffer += nb;
          available -= nb;
          nb = double_to_buffer( c.imag(), buffer, available );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::VEC_STRING:
        nb = int32_to_buffer( static_cast<int32_t>( _v_s().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & s : _v_s() )
        {
          sz = static_cast<int32_t>( s.length() + 1 );
          nb = int32_to_buffer( sz, buffer, available );
          buffer += nb;
          available -= nb;
          GC_assert( sz <= available, serialize_msg );
          memcpy( buffer, s.c_str(), static_cast<size_t>( sz ) );
          buffer += sz;
          available -= sz;
        }
        break;
      case GC_type::VECTOR:
        nb = int32_to_buffer( static_cast<int32_t>( _v().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & S : _v() )
        {
          nb = S.serialize( available, buffer );
          buffer += nb;
          available -= nb;
        }
        break;
      case GC_type::MAP:
        nb = int32_to_buffer( static_cast<int32_t>( _m().size() ), buffer, available );
        buffer += nb;
        available -= nb;
        for ( auto & [fst, snd] : _m() )
        {
          sz = static_cast<int32_t>( fst.length() + 1 );
          nb = int32_to_buffer( sz, buffer, available );
          buffer += nb;
          available -= nb;
          GC_assert( sz <= available, serialize_msg );
          memcpy( buffer, fst.c_str(), static_cast<size_t>( sz ) );
          buffer += sz;
          available -= sz;
          nb = snd.serialize( available, buffer );
          buffer += nb;
          available -= nb;
        }
        break;
    }
    return buffer_dim - available;
  }

  int32_t GenericContainer::de_serialize( std::vector<uint8_t> const & buffer )
  {
    int32_t const sz{ static_cast<int32_t>( buffer.size() ) };
    return this->de_serialize( sz, buffer.data() );
  }

  int32_t GenericContainer::de_serialize( int32_t buffer_dim, uint8_t const * buffer )
  {
    // int_type ptr_size = 8;
    int32_t sz, nb, nbyte;
    double  bf;

    this->clear();

    int32_t nr, nc, i32;
    nbyte = nb = buffer_to_int32( buffer, buffer_dim, &i32 );
    buffer += nb;

    GC_assert(
      i32 >= 0 && i32 <= static_cast<int32_t>( GC_type::MAP ),
      "GenericContainer::de_serialize, invalid type tag {}",
      i32
    );
    switch ( static_cast<TypeAllowed>( i32 ) )
    {
      case GC_type::NOTYPE: m_data.emplace<std::monostate>(); break;
      case GC_type::BOOL:
        m_data.emplace<bool_type>();
        uint8_t b;
        nb = buffer_to_uint8( buffer, buffer_dim - nbyte, &b );
        buffer += nb;
        nbyte += nb;
        _b() = b > 0;
        break;
      case GC_type::INTEGER:
        m_data.emplace<int_type>();
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &_i() );
        buffer += nb;
        nbyte += nb;
        break;
      case GC_type::LONG:
        m_data.emplace<long_type>();
        nb = buffer_to_int64( buffer, buffer_dim - nbyte, &_l() );
        buffer += nb;
        nbyte += nb;
        break;
      case GC_type::REAL:
        m_data.emplace<real_type>();
        nb = buffer_to_double( buffer, buffer_dim - nbyte, &_r() );
        buffer += nb;
        nbyte += nb;
        break;
      case GC_type::POINTER:
      {
        m_data.emplace<pointer_type>();
        int64_t ptr_value{ 0 };
        nb = buffer_to_int64( buffer, buffer_dim - nbyte, &ptr_value );
        buffer += nb;
        nbyte += nb;
        _p() = reinterpret_cast<pointer_type>( static_cast<intptr_t>( ptr_value ) );
        break;
      }
      case GC_type::STRING:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert( i32 > 0, "GenericContainer::de_serialize, invalid string length" );
        GC_assert( i32 <= buffer_dim - nbyte, deserialize_msg );
        allocate_string();
        _s() = string_type( reinterpret_cast<char const *>( buffer ), static_cast<size_t>( i32 - 1 ) );
        buffer += i32;
        nbyte += i32;
        break;
      case GC_type::COMPLEX:
        allocate_complex();
        nb = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
        _c().real( bf );
        buffer += nb;
        nbyte += nb;
        nb = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
        _c().imag( bf );
        buffer += nb;
        nbyte += nb;
        break;
      case GC_type::VEC_POINTER:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * 8u <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_pointer( static_cast<std::size_t>( i32 ) );
        for ( auto & p : _v_p() )
        {
          int64_t ptr_value{ 0 };
          sz = buffer_to_int64( buffer, buffer_dim - nbyte, &ptr_value );
          p  = reinterpret_cast<pointer_type>( static_cast<intptr_t>( ptr_value ) );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::VEC_BOOL:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * 1u <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_bool( 0 );
        _v_b().reserve( static_cast<size_t>( i32 ) );
        for ( int32_t i = 0; i < i32; ++i )
        {
          uint8_t i8;
          sz = buffer_to_uint8( buffer, buffer_dim - nbyte, &i8 );
          buffer += sz;
          nbyte += sz;
          _v_b().push_back( i8 > 0 );
        }
        break;
      case GC_type::VEC_INTEGER:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * sizeof( int_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_int( static_cast<std::size_t>( i32 ) );
        for ( auto & i : _v_i() )
        {
          sz = buffer_to_int32( buffer, buffer_dim - nbyte, &i );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::VEC_LONG:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * sizeof( long_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_long( static_cast<std::size_t>( i32 ) );
        for ( auto & i : _v_l() )
        {
          sz = buffer_to_int64( buffer, buffer_dim - nbyte, &i );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::VEC_REAL:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * sizeof( real_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_real( static_cast<std::size_t>( i32 ) );
        for ( auto & r : _v_r() )
        {
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &r );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::VEC_COMPLEX:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * sizeof( complex_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_complex( static_cast<std::size_t>( i32 ) );
        for ( auto & c : _v_c() )
        {
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
          c.real( bf );
          buffer += sz;
          nbyte += sz;
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
          c.imag( bf );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::MAT_INTEGER:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nr );
        buffer += nb;
        nbyte += nb;
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nc );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          nr >= 0 && nc >= 0 && uint64_t( nr ) * uint64_t( nc ) * sizeof( int_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized matrix dimensions" );
        allocate_mat_int( static_cast<std::size_t>( nr ), static_cast<std::size_t>( nc ) );
        for ( auto & i : _m_i() )
        {
          sz = buffer_to_int32( buffer, buffer_dim - nbyte, &i );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::MAT_LONG:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nr );
        buffer += nb;
        nbyte += nb;
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nc );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          nr >= 0 && nc >= 0 && uint64_t( nr ) * uint64_t( nc ) * sizeof( long_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized matrix dimensions" );
        allocate_mat_long( static_cast<std::size_t>( nr ), static_cast<std::size_t>( nc ) );
        for ( auto & i : _m_l() )
        {
          sz = buffer_to_int64( buffer, buffer_dim - nbyte, &i );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::MAT_REAL:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nr );
        buffer += nb;
        nbyte += nb;
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nc );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          nr >= 0 && nc >= 0 && uint64_t( nr ) * uint64_t( nc ) * sizeof( real_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized matrix dimensions" );
        allocate_mat_real( static_cast<std::size_t>( nr ), static_cast<std::size_t>( nc ) );
        for ( auto & r : _m_r() )
        {
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &r );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::MAT_COMPLEX:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nr );
        buffer += nb;
        nbyte += nb;
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &nc );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          nr >= 0 && nc >= 0 && uint64_t( nr ) * uint64_t( nc ) * sizeof( complex_type ) <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized matrix dimensions" );
        allocate_mat_complex( static_cast<std::size_t>( nr ), static_cast<std::size_t>( nc ) );
        for ( auto & c : _m_c() )
        {
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
          c.real( bf );
          buffer += sz;
          nbyte += sz;
          sz = buffer_to_double( buffer, buffer_dim - nbyte, &bf );
          c.imag( bf );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::VEC_STRING:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * 4u <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vec_string( static_cast<std::size_t>( i32 ) );
        for ( auto & s : _v_s() )
        {
          nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
          buffer += nb;
          nbyte += nb;
          GC_assert( i32 > 0, "GenericContainer::de_serialize, invalid string length" );
          GC_assert( i32 <= buffer_dim - nbyte, deserialize_msg );
          s = string_type( reinterpret_cast<char const *>( buffer ), static_cast<size_t>( i32 - 1 ) );
          buffer += i32;
          nbyte += i32;
        }
        break;
      case GC_type::VECTOR:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        buffer += nb;
        nbyte += nb;
        GC_assert(
          i32 >= 0 && uint64_t( i32 ) * 4u <= uint64_t( buffer_dim - nbyte ),
          "GenericContainer::de_serialize, invalid or oversized vector size" );
        allocate_vector( static_cast<std::size_t>( i32 ) );
        for ( auto & S : _v() )
        {
          sz = S.de_serialize( buffer_dim - nbyte, buffer );
          buffer += sz;
          nbyte += sz;
        }
        break;
      case GC_type::MAP:
        nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
        nr = i32;
        buffer += nb;
        nbyte += nb;
        GC_assert( nr >= 0, "GenericContainer::de_serialize, invalid map size" );
        allocate_map();
        for ( int32_t i = 0; i < nr; ++i )
        {
          nb = buffer_to_int32( buffer, buffer_dim - nbyte, &i32 );
          buffer += nb;
          nbyte += nb;
          GC_assert( i32 > 0, "GenericContainer::de_serialize, invalid key length" );
          GC_assert( i32 <= buffer_dim - nbyte, deserialize_msg );
          string_type key( reinterpret_cast<char const *>( buffer ), static_cast<size_t>( i32 - 1 ) );
          buffer += i32;
          nbyte += i32;
          // std::cout << "key[" << i << "] = " << key << '\n';
          sz = _m()[key].de_serialize( buffer_dim - nbyte, buffer );
          buffer += sz;
          nbyte += sz;
        }
        break;
    }
    GC_assert( nbyte <= buffer_dim, "GenericContainer::serialize, buffer overflow" );
    return nbyte;
  }

}  // namespace GC_namespace

//
// eof: GenericContainerSerialize.cc
//
