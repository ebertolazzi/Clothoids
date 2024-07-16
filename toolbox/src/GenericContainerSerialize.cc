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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: GenericContainerSerialize.cc
//

#include "GenericContainer/GenericContainer.hh"
#include "GenericContainer/GenericContainerLibs.hh"
#include <cstring>

namespace GC_namespace {

  using std::memcpy;

  static
  uint32_t
  int8_to_buffer( int8_t in, uint8_t * buffer ) {
    buffer[0] = *reinterpret_cast<uint8_t*>(&in);
    return sizeof(int8_t);
  }

  static
  uint32_t
  uint32_to_buffer( uint32_t in, uint8_t * buffer ) {
    buffer[0] = uint8_t(in&0xFF); in >>= 8;
    buffer[1] = uint8_t(in&0xFF); in >>= 8;
    buffer[2] = uint8_t(in&0xFF); in >>= 8;
    buffer[3] = uint8_t(in&0xFF);
    return sizeof(uint32_t);
  }

  static
  uint32_t
  int32_to_buffer( int32_t in, uint8_t * buffer ) {
    return uint32_to_buffer( *reinterpret_cast<uint32_t*>(&in), buffer );
  }

  static
  uint32_t
  uint64_to_buffer( uint64_t in, uint8_t * buffer ) {
    buffer[0] = uint8_t(in&0xFF); in >>= 8;
    buffer[1] = uint8_t(in&0xFF); in >>= 8;
    buffer[2] = uint8_t(in&0xFF); in >>= 8;
    buffer[3] = uint8_t(in&0xFF); in >>= 8;
    buffer[4] = uint8_t(in&0xFF); in >>= 8;
    buffer[5] = uint8_t(in&0xFF); in >>= 8;
    buffer[6] = uint8_t(in&0xFF); in >>= 8;
    buffer[7] = uint8_t(in&0xFF);
    return sizeof(int64_t);
  }

  static
  uint32_t
  int64_to_buffer( int64_t in, uint8_t * buffer ) {
    return uint64_to_buffer( *reinterpret_cast<uint64_t*>(&in), buffer );
  }

  static
  uint32_t
  double_to_buffer( double in, uint8_t * buffer ) {
    union { double f; uint64_t i; } tmp;
    tmp.f = in;
    uint64_to_buffer( tmp.i, buffer );
    return sizeof(double);
  }

  /* ---------------------------------------------------------------------------- */

  static
  uint32_t
  buffer_to_uint8( uint8_t const * buffer, uint8_t * out ) {
    *out = buffer[0];
    return sizeof(uint8_t);
  }

  static
  uint32_t
  buffer_to_uint32( uint8_t const * buffer, uint32_t * out ) {
    uint32_t tmp0 = buffer[0];
    uint32_t tmp1 = buffer[1];
    uint32_t tmp2 = buffer[2];
    uint32_t tmp3 = buffer[3];
    *out = tmp0|(tmp1<<8)|(tmp2<<16)|(tmp3<<24);
    return sizeof(int32_t);
  }

  static
  uint32_t
  buffer_to_int32( uint8_t const * buffer, int32_t * out ) {
    return buffer_to_uint32( buffer, reinterpret_cast<uint32_t*>(out) );
  }

  static
  uint32_t
  buffer_to_uint64( uint8_t const * buffer, uint64_t * out ) {
    uint64_t tmp0 = buffer[0];
    uint64_t tmp1 = buffer[1];
    uint64_t tmp2 = buffer[2];
    uint64_t tmp3 = buffer[3];
    uint64_t tmp4 = buffer[4];
    uint64_t tmp5 = buffer[5];
    uint64_t tmp6 = buffer[6];
    uint64_t tmp7 = buffer[7];
    *out = tmp0|(tmp1<<8)|(tmp2<<16)|(tmp3<<24)|(tmp4<<32)|(tmp5<<40)|(tmp6<<48)|(tmp7<<56);
    return sizeof(uint64_t);
  }

  static
  uint32_t
  buffer_to_int64( uint8_t const * buffer, int64_t * out ) {
    return buffer_to_uint64( buffer, reinterpret_cast<uint64_t*>(out) );
  }

  static
  uint32_t
  buffer_to_double( uint8_t const * buffer, double * out ) {
    union { double f; uint64_t i; } tmp;
    buffer_to_uint64( buffer, &tmp.i );
    *out = tmp.f;
    return sizeof(double);
  }

  int32_t
  GenericContainer::mem_size() const {
    int32_t header_size = sizeof(int32_t);
    int32_t ptr_size    = 8;
    int32_t res         = 0;
    switch (m_data_type) {
    case GC_type::NOTYPE:      res = int32_t(header_size); break;
    case GC_type::BOOL:        res = int32_t(header_size+1); break;
    case GC_type::INTEGER:     res = int32_t(header_size+sizeof(int_type)); break;
    case GC_type::LONG:        res = int32_t(header_size+sizeof(long_type)); break;
    case GC_type::REAL:        res = int32_t(header_size+sizeof(real_type)); break;
    case GC_type::POINTER:     res = int32_t(header_size+ptr_size); break;
    case GC_type::STRING:      res = int32_t(header_size+m_data.s->length()+5); break;
    case GC_type::COMPLEX:     res = int32_t(header_size+sizeof(complex_type)); break;
    case GC_type::VEC_POINTER: res = int32_t(header_size+sizeof(int32_t)+ptr_size*m_data.v_p->size()); break;
    case GC_type::VEC_BOOL:    res = int32_t(header_size+sizeof(int32_t)+m_data.v_b->size()); break;
    case GC_type::VEC_INTEGER: res = int32_t(header_size+sizeof(int32_t)+sizeof(int_type)*m_data.v_i->size()); break;
    case GC_type::VEC_LONG:    res = int32_t(header_size+sizeof(int32_t)+sizeof(long_type)*m_data.v_l->size()); break;
    case GC_type::VEC_REAL:    res = int32_t(header_size+sizeof(int32_t)+sizeof(real_type)*m_data.v_r->size()); break;
    case GC_type::VEC_COMPLEX: res = int32_t(header_size+sizeof(int32_t)+sizeof(complex_type)*m_data.v_c->size()); break;

    case GC_type::MAT_INTEGER: res = int32_t(header_size+2*sizeof(int32_t)+sizeof(int_type)*m_data.m_i->size()); break;
    case GC_type::MAT_LONG:    res = int32_t(header_size+2*sizeof(int32_t)+sizeof(long_type)*m_data.m_l->size()); break;
    case GC_type::MAT_REAL:    res = int32_t(header_size+2*sizeof(int32_t)+sizeof(real_type)*m_data.m_r->size()); break;
    case GC_type::MAT_COMPLEX: res = int32_t(header_size+2*sizeof(int32_t)+sizeof(complex_type)*m_data.m_c->size()); break;
    case GC_type::VEC_STRING:
      res = int32_t(header_size+sizeof(int32_t));
      for ( auto & s : *m_data.v_s ) res += int32_t(sizeof(int32_t)+s.length()+1);
      break;
    case GC_type::VECTOR:
      res = int32_t(header_size+sizeof(int32_t));
      for ( auto & S : *m_data.v ) res += S.mem_size();
      break;
    case GC_type::MAP:
      res = int32_t(header_size+sizeof(int32_t));
      for ( auto & S : *m_data.m )
        res += int32_t(sizeof(int32_t)+S.first.length()+1+S.second.mem_size());
      break;
    }
    return res;
  }

  int32_t
  GenericContainer::serialize( std::vector<uint8_t> & buffer ) const {
    int32_t sz = int32_t(buffer.size());
    return this->serialize( sz, buffer.data() );
  }

  int32_t
  GenericContainer::serialize( int32_t buffer_dim, uint8_t * buffer ) const {
    //int_type ptr_size = 8;
    int32_t sz, nb, nbyte;

    nbyte = nb = int32_to_buffer( static_cast<int32_t>(m_data_type), buffer );
    buffer += nb;

    switch (m_data_type) {
    case GC_type::NOTYPE:
      break;
    case GC_type::BOOL:
      nb = int8_to_buffer( m_data.b ? 1 : 0, buffer );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::INTEGER:
      nb = int32_to_buffer( m_data.i, buffer );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::LONG:
      nb = int64_to_buffer( m_data.l, buffer );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::REAL:
      nb = double_to_buffer( m_data.r, buffer );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::POINTER:
      break;
    case GC_type::STRING:
      sz = int32_t(m_data.s->length()+1);
      nb = int32_to_buffer( sz, buffer );
      buffer += nb; nbyte += nb;
      memcpy( buffer, &m_data.s->front(), sz );
      buffer += sz; nbyte += sz;
      break;
    case GC_type::COMPLEX:
      nb = double_to_buffer( m_data.c->real(), buffer );
      buffer += nb; nbyte += nb;
      nb = double_to_buffer( m_data.c->imag(), buffer );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::VEC_POINTER:
      break;
    case GC_type::VEC_BOOL:
      nb = int32_to_buffer( int32_t(m_data.v_b->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto b : *m_data.v_b ) {
        int8_to_buffer( b ? 1 : 0, buffer );
        ++nb; ++buffer; ++nbyte;
      }
      break;
    case GC_type::VEC_INTEGER:
      nb = int32_to_buffer( int32_t(m_data.v_i->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & i : *m_data.v_i ) {
        sz = int32_to_buffer( i, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_LONG:
      nb = int32_to_buffer( int32_t(m_data.v_l->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & i : *m_data.v_l ) {
        sz = int64_to_buffer( i, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_REAL:
      nb = int32_to_buffer( int32_t(m_data.v_r->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & r : *m_data.v_r ) {
        sz = double_to_buffer( r, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_COMPLEX:
      nb = int32_to_buffer( int32_t(m_data.v_c->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & c : *m_data.v_c ) {
        sz = double_to_buffer( c.real(), buffer );
        nb += sz; buffer += sz; nbyte += sz;
        sz = double_to_buffer( c.imag(), buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_INTEGER:
      nb = int32_to_buffer( m_data.m_i->num_rows(), buffer );
      buffer += nb; nbyte += nb;
      nb = int32_to_buffer( m_data.m_i->num_cols(), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & i : *m_data.m_i ) {
        sz = int32_to_buffer( i, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_LONG:
      nb = int32_to_buffer( m_data.m_l->num_rows(), buffer );
      buffer += nb; nbyte += nb;
      nb = int32_to_buffer( m_data.m_l->num_cols(), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & i : *m_data.m_l ) {
        sz = int64_to_buffer( i, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_REAL:
      nb = int32_to_buffer( m_data.m_r->num_rows(), buffer );
      buffer += nb; nbyte += nb;
      nb = int32_to_buffer( m_data.m_r->num_cols(), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & r : *m_data.m_r ) {
        sz = double_to_buffer( r, buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_COMPLEX:
      nb = int32_to_buffer( m_data.m_c->num_rows(), buffer );
      buffer += nb; nbyte += nb;
      nb = int32_to_buffer( m_data.m_c->num_cols(), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & c : *m_data.m_c ) {
        sz = double_to_buffer( c.real(), buffer );
        nb += sz; buffer += sz; nbyte += sz;
        sz = double_to_buffer( c.imag(), buffer );
        nb += sz; buffer += sz; nbyte += sz;
      }
    case GC_type::VEC_STRING:
      nb = int32_to_buffer( int32_t(m_data.v_s->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & s : *m_data.v_s ) {
        sz = int32_t(s.length()+1);
        nb = int32_to_buffer( sz, buffer );
        buffer += nb; nbyte += nb;
        memcpy( buffer, &s.front(), sz );
        buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VECTOR:
      nb = int32_to_buffer( int32_t(m_data.v->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & S : *m_data.v ) {
        sz = S.serialize( buffer_dim-nbyte, buffer );
        buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAP:
      nb = int32_to_buffer( int32_t(m_data.m->size()), buffer );
      buffer += nb; nbyte += nb;
      for ( auto & S : *m_data.m ) {
        sz = int32_t(S.first.length()+1);
        nb = int32_to_buffer( sz, buffer );
        buffer += nb; nbyte += nb;
        memcpy( buffer, &S.first.front(), sz );
        buffer += sz; nbyte += sz;
        nb = S.second.serialize( buffer_dim-nbyte, buffer );
        buffer += nb; nbyte += nb;
      }
      break;
    }
    GC_ASSERT(
      nbyte <= buffer_dim,
      "GenericContainer::serialize, buffer overflow"
    );
    return nbyte;
  }

  int32_t
  GenericContainer::de_serialize( std::vector<uint8_t> const & buffer ) {
    int32_t sz = int32_t(buffer.size());
    return this->de_serialize( sz, buffer.data() );
  }

  int32_t
  GenericContainer::de_serialize( int32_t buffer_dim, uint8_t const * buffer ) {
    //int_type ptr_size = 8;
    int32_t sz, nb, nbyte;
    double  bf;

    this->clear();

    int32_t nr, nc, i32;
    nbyte = nb = buffer_to_int32( buffer, &i32 );
    buffer += nb;

    switch (TypeAllowed(i32)) {
    case GC_type::NOTYPE:
      m_data_type = GC_type::NOTYPE;
      break;
    case GC_type::BOOL:
      m_data_type = GC_type::BOOL;
      uint8_t b;
      nb = buffer_to_uint8( buffer, &b );
      buffer += nb; nbyte += nb;
      m_data.b = b > 0;
      break;
    case GC_type::INTEGER:
      m_data_type = GC_type::INTEGER;
      nb = buffer_to_int32( buffer, &m_data.i );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::LONG:
      m_data_type = GC_type::LONG;
      nb = buffer_to_int64( buffer, &m_data.l );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::REAL:
      m_data_type = GC_type::REAL;
      nb = buffer_to_double( buffer, &m_data.r );
      buffer += nb; nbyte += nb;
      break;
    case GC_type::POINTER:
      break;
    case GC_type::STRING:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_string();
      *m_data.s = reinterpret_cast<char const*>(buffer);
      buffer += i32; nbyte += i32;
      break;
    case GC_type::COMPLEX:
      allocate_complex();
      nb = buffer_to_double( buffer, &bf ); m_data.c->real(bf);
      buffer += nb; nbyte += nb;
      nb = buffer_to_double( buffer, &bf ); m_data.c->imag(bf);
      buffer += nb; nbyte += nb;
      break;
    case GC_type::VEC_POINTER:
      break;
    case GC_type::VEC_BOOL:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_bool(0);
      m_data.v_b->reserve(i32);
      for ( int32_t i = 0; i < i32; ++i ) {
        uint8_t i8;
        buffer_to_uint8( buffer, &i8 );
        ++nb; ++buffer; ++nbyte;
        m_data.v_b->push_back( i8 > 0 );
      }
      break;
    case GC_type::VEC_INTEGER:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_int(i32);
      for ( auto & i : *m_data.v_i ) {
        sz = buffer_to_int32( buffer, &i );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_LONG:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_long(i32);
      for ( auto & i : *m_data.v_l ) {
        sz = buffer_to_int64( buffer, &i );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_REAL:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_real(i32);
      for ( auto & r : *m_data.v_r ) {
        sz = buffer_to_double( buffer, &r );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::VEC_COMPLEX:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_complex(i32);
      for ( auto & c : *m_data.v_c ) {
        sz = buffer_to_double( buffer, &bf ); c.real(bf);
        nb += sz; buffer += sz; nbyte += sz;
        sz = buffer_to_double( buffer, &bf ); c.imag(bf);
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_INTEGER:
      nb = buffer_to_int32( buffer, &nr );
      buffer += nb; nbyte += nb;
      nb = buffer_to_int32( buffer, &nc );
      buffer += nb; nbyte += nb;
      allocate_mat_int( nr, nc );
      for ( auto & i : *m_data.m_i ) {
        sz = buffer_to_int32( buffer, &i );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_LONG:
      nb = buffer_to_int32( buffer, &nr );
      buffer += nb; nbyte += nb;
      nb = buffer_to_int32( buffer, &nc );
      buffer += nb; nbyte += nb;
      allocate_mat_long( nr, nc );
      for ( auto & i : *m_data.m_l ) {
        sz = buffer_to_int64( buffer, &i );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_REAL:
      nb = buffer_to_int32( buffer, &nr );
      buffer += nb; nbyte += nb;
      nb = buffer_to_int32( buffer, &nc );
      buffer += nb; nbyte += nb;
      allocate_mat_real( nr, nc );
      for ( auto & r : *m_data.m_r ) {
        sz = buffer_to_double( buffer, &r );
        nb += sz; buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAT_COMPLEX:
      nb = buffer_to_int32( buffer, &nr );
      buffer += nb; nbyte += nb;
      nb = buffer_to_int32( buffer, &nc );
      buffer += nb; nbyte += nb;
      allocate_mat_complex( nr, nc );
      for ( auto & c : *m_data.m_c ) {
        sz = buffer_to_double( buffer, &bf ); c.real(bf);
        nb += sz; buffer += sz; nbyte += sz;
        sz = buffer_to_double( buffer, &bf ); c.imag(bf);
        nb += sz; buffer += sz; nbyte += sz;
      }
    case GC_type::VEC_STRING:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vec_string(i32);
      for ( auto & s : *m_data.v_s ) {
        nb = buffer_to_int32( buffer, &i32 );
        buffer += nb; nbyte += nb;
        s = reinterpret_cast<char const*>(buffer);
        buffer += i32; nbyte += i32;
      }
      break;
    case GC_type::VECTOR:
      nb = buffer_to_int32( buffer, &i32 );
      buffer += nb; nbyte += nb;
      allocate_vector(i32);
      for ( auto & S : *m_data.v ) {
        sz = S.de_serialize( buffer_dim-nbyte, buffer );
        buffer += sz; nbyte += sz;
      }
      break;
    case GC_type::MAP:
      nb = buffer_to_int32( buffer, &i32 ); nr = i32;
      buffer += nb; nbyte += nb;
      allocate_map();
      for ( int32_t i = 0; i < nr; ++i ) {
        nb = buffer_to_int32( buffer, &i32 );
        buffer += nb; nbyte += nb;
        string_type key = reinterpret_cast<char const*>(buffer);
        buffer += i32; nbyte += i32;
        //std::cout << "key[" << i << "] = " << key << '\n';
        sz = (*m_data.m)[key].de_serialize( buffer_dim-nbyte, buffer );
        buffer += sz; nbyte += sz;
      }
      break;
    }
    GC_ASSERT(
      nbyte <= buffer_dim,
      "GenericContainer::serialize, buffer overflow"
    );
    return nbyte;
  }

}

//
// eof: GenericContainerSerialize.cc
//
