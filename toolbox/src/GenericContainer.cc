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
  #pragma warning(disable : 4661)
#endif

#include "GenericContainer/GenericContainer.hh"
#include <iomanip>
#include <cmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
#else
  #error This library needs at least a C++11 compliant compiler
#endif

#include <regex>
#include <fstream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define CHECK_RESIZE(pV,I) if ( pV->size() <= (I) ) pV->resize((I)+1)

using std::fpclassify;
using GC_namespace::real_type;

//static
//inline
//bool isZero( real_type x )
//{ return FP_ZERO == fpclassify(x); }

static bool isZero0( real_type const x )
{ int const c = fpclassify(x); return FP_ZERO == c || FP_SUBNORMAL == c; }

static bool isInteger( real_type const x )
{ using std::round; return isZero0( x-round(x) ); }

static bool isUnsigned( real_type const x )
{ return isInteger(x) && x >= 0; }

#endif

namespace GC_namespace {

  //!
  //! precision used in printing number
  //!
  unsigned stream_number_precision{12};

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  string
  to_string( complex_type const & v ) {
    ostringstream data;
    data.precision(stream_number_precision);
    data << v.real();
    if ( v.imag() > 0 ) data << '+' << v.imag() << 'i';
    if ( v.imag() < 0 ) data << '-' << -v.imag() << 'i';
    return data.str();
  }

  void
  string_escape( ostream_type & stream, string const & s ) {
    stream << '"';
    for ( auto const c : s ) {
      if      ( c == '"'  ) { stream << "\\\""; }
      else if ( c == '\n' ) { stream << "\\n"; }
      else if ( c == '\r' ) { stream << "\\r"; }
      else if ( c == '\t' ) { stream << "\\t"; }
      else if ( c == '\v' ) { stream << "\\v"; }
      else if ( c == '\b' ) { stream << "\\b"; }
      else if ( c == '\a' ) { stream << "\\a"; }
      else if ( c == '\\' ) { stream << "\\\\"; }
      else                    stream << c;
    }
    stream << '"';
  }

  template <typename TYPE>
  ostream_type &
  operator << ( ostream_type & s, vector<TYPE> const & v ) {
    s << '[';
    for ( TYPE const & vi : v ) s << ' ' << vi;
    s << " ]";
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, vec_bool_type const & v ) {
    s << '[';
    for ( bool const vi : v ) s << (vi?" true":" false");
    s << " ]";
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, vec_complex_type const & v ) {
    s << '[';
    for ( complex_type const & vi : v ) s << ' ' << to_string(vi);
    s << " ]";
    return s;
  }

  template ostream_type & operator << ( ostream_type & s, vec_int_type const & v );
  template ostream_type & operator << ( ostream_type & s, vec_long_type const & v );
  template ostream_type & operator << ( ostream_type & s, vec_real_type const & v );
  //template ostream_type & operator << ( ostream_type & s, vec_complex_type const & v );

  #endif

  template <typename TYPE>
  TYPE const &
  mat_type<TYPE>::operator () ( unsigned i, unsigned j ) const {
    try {
      return this->at(i+j*m_num_rows);
    } catch ( std::exception const & exc ) {
      GC_DO_ERROR( "mat_type::operator() (" << i << ", " << j << "): " << exc.what() << '\n' );
    } catch ( ... ) {
      GC_DO_ERROR( "mat_type::operator() (" << i << ", " << j << "): unknown error\n" );
    }
  }

  template <typename TYPE>
  TYPE &
  mat_type<TYPE>::operator () ( unsigned i, unsigned j ) {
    try {
      return this->at(i+j*m_num_rows);
    } catch ( std::exception const & exc ) {
      GC_DO_ERROR( "mat_type::operator() (" << i << ", " << j << "): " << exc.what() << '\n' );
    } catch ( ... ) {
      GC_DO_ERROR( "mat_type::operator() (" << i << ", " << j << "): unknown error\n" );
    }
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::get_column( unsigned nc, std::vector<TYPE> & C ) const {
    GC_ASSERT(
      nc < m_num_cols,
      "mat_type::get_column(" << nc << ",C) column index out of range max = " << m_num_cols-1
    );
    C.clear();
    C.reserve(m_num_rows);
    for ( unsigned i{0}; i < m_num_rows; ++i )
      C.push_back( (*this)(i,nc) );
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::get_column( unsigned nc, TYPE * C ) const {
    GC_ASSERT(
      nc < m_num_cols,
      "mat_type::get_column(" << nc << ",C) column index out of range max = " << m_num_cols-1
    );
    for ( unsigned i{0}; i < m_num_rows; ++i )
      *C++ = (*this)(i,nc);
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::get_row( unsigned nr, std::vector<TYPE> & R ) const {
    GC_ASSERT(
      nr < m_num_rows,
      "mat_type::get_row(" << nr << ",C) row index out of range max = " << m_num_rows-1
    );
    R.clear();
    R.reserve(m_num_cols);
    for ( unsigned j{0}; j < m_num_cols; ++j )
      R.push_back( (*this)(nr,j) );
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::get_row( unsigned nr, TYPE * R ) const {
    GC_ASSERT(
      nr < m_num_rows,
      "mat_type::get_row(" << nr << ",C) row index out of range max = " << m_num_rows-1
    );
    for ( unsigned j{0}; j < m_num_cols; ++j )
      *R++ = (*this)(nr,j);
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::info( ostream_type & stream ) const {
    stream
      << "Matrix of floating point number of size "
      << m_num_rows << " x " << m_num_cols
      << '\n';
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename TYPE>
  ostream_type &
  operator << ( ostream_type & s, mat_type<TYPE> const & mat ) {
    if ( mat.num_rows() > 0 && mat.num_cols() ) {
      for ( unsigned i{0}; i < mat.num_rows(); ++i ) {
        s << std::setw(8) << mat(i,0);
        for ( unsigned j{1}; j < mat.num_cols(); ++j )
          s << " " << std::setw(8) << mat(i,j);
        s << '\n';
      }
    } else {
      s << mat.num_rows() << " by " << mat.num_cols() << " matrix\n";
    }
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, mat_complex_type const & m ) {
    if ( m.num_rows() > 0 && m.num_cols() ) {
      for ( unsigned i{0}; i < m.num_rows(); ++i ) {
        s << std::setw(8) << m(i,0);
        for ( unsigned j{1}; j < m.num_cols(); ++j )
          s << " " << std::setw(12) << to_string(m(i,j));
        s << '\n';
      }
    } else {
      s << m.num_rows() << " by " << m.num_cols() << " matrix\n";
    }
    return s;
  }

  template ostream_type & operator << ( ostream_type & s, mat_type<int_type>     const & m );
  template ostream_type & operator << ( ostream_type & s, mat_type<long_type>    const & m );
  template ostream_type & operator << ( ostream_type & s, mat_type<real_type>    const & m );
  //template ostream_type & operator << ( ostream_type & s, mat_type<complex_type> const & m );

  class Pcre_for_GC {

    std::regex  reCompiled;
    std::smatch reMatches;

  public:

    Pcre_for_GC()
    : reCompiled(R"(^\s*\d+\s*(##?)(-|=|~|_|)\s*(.*)$)")
    { }

    ~Pcre_for_GC() = default;

    int
    exec( string_view const str_in, string_type matches[4] ) {
      if ( string const str(str_in); std::regex_match( str, reMatches, reCompiled ) ) {
        for ( size_t i{0}; i < reMatches.size(); ++i ) matches[i] = reMatches[i].str();
        return static_cast<int>(reMatches.size());
      }
      return 0;
    }
  };

  static Pcre_for_GC pcre_for_GC;

  #endif

  string_view
  to_string( GC_type const s ) {
    switch ( s ) {
      case GC_type::NOTYPE:      return "NOTYPE";
      case GC_type::POINTER:     return "pointer";
      case GC_type::BOOL:        return "bool_type";
      case GC_type::INTEGER:     return "int_type";
      case GC_type::LONG:        return "long_type";
      case GC_type::REAL:        return "real_type";
      case GC_type::COMPLEX:     return "complex_type";
      case GC_type::STRING:      return "string_type";
      case GC_type::VEC_POINTER: return "vec_pointer_type";
      case GC_type::VEC_BOOL:    return "vec_bool_type";
      case GC_type::VEC_INTEGER: return "vec_int_type";
      case GC_type::VEC_LONG:    return "vec_long_type";
      case GC_type::VEC_REAL:    return "vec_real_type";
      case GC_type::VEC_COMPLEX: return "vec_complex_type";
      case GC_type::VEC_STRING:  return "vec_string_type";
      case GC_type::MAT_INTEGER: return "mat_int_type";
      case GC_type::MAT_LONG:    return "mat_long_type";
      case GC_type::MAT_REAL:    return "mat_real_type";
      case GC_type::MAT_COMPLEX: return "mat_complex_type";
      case GC_type::VECTOR:      return "vector_type";
      case GC_type::MAP:         return "map_type";
    }
    return "";
  }

  #ifdef GENERIC_CONTAINER_ON_WINDOWS
  bool
  GenericContainer::simple_data() const {
    return m_data_type <= GC_type::STRING;
  }
  bool
  GenericContainer::simple_vec_data() const {
    return m_data_type < GC_type::VEC_STRING;
  }
  #endif

  void
  GenericContainer::get_keys( vec_string_type & keys ) const {
    keys.clear();
    if ( GC_type::MAP == m_data_type ) {
      keys.reserve( m_data.m->size() );
      for ( const auto &[fst, snd] : *m_data.m ) keys.emplace_back( fst );
    }
  }

  string
  GenericContainer::get_keys() const {
    string res;
    if ( GC_type::MAP == m_data_type ) {
      for ( const auto &[fst, snd] : *m_data.m ) { res += fst; res += ", "; }
      if ( !m_data.m->empty() ) { res.pop_back(); res.pop_back(); }
    }
    return res;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_bool_type const & a ) {
    set_vec_bool( static_cast<unsigned>(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_b->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_int_type const & a ) {
    set_vec_int( static_cast<unsigned>(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_i->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_long_type const & a ) {
    set_vec_long( a.size() );
    std::copy( a.begin(), a.end(), m_data.v_l->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_real_type const & a ) {
    set_vec_real( a.size() );
    std::copy( a.begin(), a.end(), m_data.v_r->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_complex_type const & a ) {
    set_vec_complex( a.size() );
    std::copy( a.begin(), a.end(), m_data.v_c->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_string_type const & a ) {
    set_vec_string( a.size() );
    std::copy( a.begin(), a.end(), m_data.v_s->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( mat_int_type const & a ) {
    set_mat_int( a );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( mat_long_type const & a ) {
    set_mat_long( a );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( mat_real_type const & a ) {
    set_mat_real( a );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( mat_complex_type const & a ) {
    set_mat_complex( a );
    return *this;
  }

  template <typename TYPE>
  static
  string
  compare_vector( string_view const who, vector<TYPE> const * A, vector<TYPE> const * B ) {
    ostringstream data;
    // controllo valori
    if ( A->size() != B->size() ) {
      data << who << " size: " << A->size() << " <> " << B->size() << '\n';
    } else {
      auto it1 = A->begin();
      auto it2 = B->begin();
      int i{0};
      while ( it1 != A->end() && it2 != B->end() ) {
        if ( *it1 != *it2 ) {
          data << who << " at " << i << " values "
               << *it1 << " <> " << *it2 << '\n';
          break;
        }
        ++i; ++it1; ++it2;
      }
    }
    return data.str();
  }

  template <typename TYPE>
  static
  string
  compare_matrix( string_view const who, mat_type<TYPE> const * A, mat_type<TYPE> const * B ) {
    ostringstream data;
    if ( A->num_rows() == B->num_rows() && A->num_cols() == B->num_cols() ) {
      for ( unsigned i{0}; i < A->num_rows(); ++i ) {
        for ( unsigned j{0}; j < A->num_cols(); ++j ) {
          TYPE const & Aij = (*A)(i,j);
          TYPE const & Bij = (*B)(i,j);
          if ( Aij != Bij ) {
            data << who << " at (" << i << "," << j << ") values "
                 << Aij << " <> " << Bij << '\n';
            break;
          }
        }
      }
    } else {
      data << who << " size: "
           << A->num_rows() << " x " << A->num_cols() << " <> "
           << B->num_rows() << " x " << B->num_cols() << '\n';
    }
    return data.str();
  }

  string
  GenericContainer::compare_content( GenericContainer const & gc, string_view from ) const {
    ostringstream data;
    if ( m_data_type != gc.m_data_type ) {
      data << from
           << "different type: "
           << to_string(m_data_type) << " <> "
           << to_string(gc.m_data_type) << '\n';
    } else {
      string tmp;
      switch (m_data_type) {
      case GC_type::NOTYPE:
        break;
      case GC_type::BOOL:
        if ( m_data.b != gc.m_data.b )
          data << from << "boolean: different\n";
        break;
      case GC_type::INTEGER:
        if ( m_data.i != gc.m_data.i )
          data << from << "integer: " << m_data.i << " <> " << gc.m_data.i << '\n';
        break;
      case GC_type::LONG:
        if ( m_data.l != gc.m_data.l )
          data << from << "long: " << m_data.l << " <> " << gc.m_data.l << '\n';
        break;
      case GC_type::REAL:
        if ( m_data.r != gc.m_data.r )
          data << from << "real: " << m_data.r << " <> " << gc.m_data.r << '\n';
        break;
      case GC_type::POINTER:
        if ( m_data.p != gc.m_data.p )
          data << from << "pointer: 0x" << std::hex << m_data.p << " <> 0x" <<  std::hex << gc.m_data.p << '\n';
        break;
      case GC_type::STRING:
        if ( *m_data.s != *gc.m_data.s )
          data << from << "string: \"" << *m_data.s << "\" <> \"" << *gc.m_data.s << "\"\n";
        break;
      case GC_type::COMPLEX:
        if ( *m_data.c != *gc.m_data.c )
          data << from << "complex: " << *m_data.c << " <> " << *gc.m_data.c << '\n';
        break;
      case GC_type::VEC_POINTER:
        tmp = compare_vector( "vector of pointer", m_data.v_p, gc.m_data.v_p );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_BOOL:
        tmp = compare_vector( "vector of boolean", m_data.v_b, gc.m_data.v_b );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_INTEGER:
        tmp = compare_vector( "vector of integer", m_data.v_i, gc.m_data.v_i );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_LONG:
        tmp = compare_vector( "vector of long", m_data.v_l, gc.m_data.v_l );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_REAL:
        tmp = compare_vector( "vector of double", m_data.v_r, gc.m_data.v_r );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_COMPLEX:
        tmp = compare_vector( "vector of complex", m_data.v_c, gc.m_data.v_c );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::MAT_INTEGER:
        tmp = compare_matrix( "mat of integer", m_data.m_i, gc.m_data.m_i );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::MAT_LONG:
        tmp = compare_matrix( "mat of long", m_data.m_l, gc.m_data.m_l );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::MAT_REAL:
        tmp = compare_matrix( "mat of double", m_data.m_r, gc.m_data.m_r );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::MAT_COMPLEX:
        tmp = compare_matrix( "mat of complex", m_data.m_c, gc.m_data.m_c );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VEC_STRING:
        tmp = compare_vector( "vector of string", m_data.v_s, gc.m_data.v_s );
        if ( !tmp.empty() ) data << from << tmp;
        break;
      case GC_type::VECTOR:
        if ( m_data.v->size() == gc.m_data.v->size() ) {
          // controllo contenutp
          auto it1 = m_data.v->begin();
          auto it2 = gc.m_data.v->begin();
          unsigned i{0};
          while ( it1 != m_data.v->end() ) {
            if ( string const res{ it1->compare_content( *it2, "> " ) }; !res.empty() ) {
              data << from << "position: " << i << '\n' << res;
              break;
            }
            ++i; ++it1; ++it2;
          }
        } else {
          data << from
               << "vector of GC size do not match: "
               << m_data.v->size() << " <> " << gc.m_data.v->size() << '\n';
        }
        break;
      case GC_type::MAP:
        if ( m_data.m->size() == gc.m_data.m->size() ) {
          // controllo le chiavi
          auto it1 = m_data.m->begin();
          auto it2 = gc.m_data.m->begin();
          while ( it1 != m_data.m->end() ) {
            if ( it1->first == it2->first ) {
              if ( string const res{ it1->second.compare_content( it2->second, "> " ) }; !res.empty() ) {
                data << from << "key: '" << it1->first << "'\n" << res;
                break;
              }
            } else {
              data << from
                   << "map of GC keys do not match: "
                   << it1->first << " <> " << it2->first << '\n';
              break;
            }
            ++it1; ++it2;
          }
        } else {
          data << from
               << "map of GC size do not match: " << m_data.m->size()
               << " <> " << gc.m_data.m->size() << '\n';
        }
        break;
      }
    }
    return data.str();
  }

  void
  GenericContainer::clear() {
    switch (m_data_type) {
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::POINTER:
      // removed annoying warning. To be re-thinked...
      //GC_WARNING( _data.p == nullptr, "find a pointer not deallocated!" );
      break;
    case GC_type::STRING:      delete m_data.s; break;
    case GC_type::COMPLEX:     delete m_data.c; break;

    case GC_type::VEC_POINTER: delete m_data.v_p; break;
    case GC_type::VEC_BOOL:    delete m_data.v_b; break;
    case GC_type::VEC_INTEGER: delete m_data.v_i; break;
    case GC_type::VEC_LONG:    delete m_data.v_l; break;
    case GC_type::VEC_REAL:    delete m_data.v_r; break;
    case GC_type::VEC_COMPLEX: delete m_data.v_c; break;

    case GC_type::MAT_INTEGER: delete m_data.m_i; break;
    case GC_type::MAT_LONG:    delete m_data.m_l; break;
    case GC_type::MAT_REAL:    delete m_data.m_r; break;
    case GC_type::MAT_COMPLEX: delete m_data.m_c; break;
    case GC_type::VEC_STRING:  delete m_data.v_s; break;

    case GC_type::VECTOR:
      for ( auto & it : *m_data.v ) it.clear();
      delete m_data.v;
      break;
    case GC_type::MAP:
      for ( auto &[fst, snd] : *m_data.m ) snd.clear();
      delete m_data.m;
      break;
    }
    m_data_type = GC_type::NOTYPE;
  }

  void
  GenericContainer::erase( string_view const name ) const {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      "GenericContainer::erase('" << name << "') bad data type\nexpect: " <<
      to_string(GC_type::POINTER) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    m_data.m->erase(string(name));
  }

  unsigned
  GenericContainer::get_num_elements() const {
    switch (m_data_type) {
    case GC_type::POINTER:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::COMPLEX:
    case GC_type::STRING:      return 1;

    case GC_type::VEC_POINTER: return static_cast<unsigned>(m_data.v_p->size());
    case GC_type::VEC_BOOL:    return static_cast<unsigned>(m_data.v_b->size());
    case GC_type::VEC_INTEGER: return static_cast<unsigned>(m_data.v_i->size());
    case GC_type::VEC_LONG:    return static_cast<unsigned>(m_data.v_l->size());
    case GC_type::VEC_REAL:    return static_cast<unsigned>(m_data.v_r->size());
    case GC_type::VEC_COMPLEX: return static_cast<unsigned>(m_data.v_c->size());
    case GC_type::VEC_STRING:  return static_cast<unsigned>(m_data.v_s->size());

    case GC_type::MAT_INTEGER: return static_cast<unsigned>(m_data.m_i->size());
    case GC_type::MAT_LONG:    return static_cast<unsigned>(m_data.m_l->size());
    case GC_type::MAT_REAL:    return static_cast<unsigned>(m_data.m_r->size());
    case GC_type::MAT_COMPLEX: return static_cast<unsigned>(m_data.m_c->size());

    case GC_type::VECTOR:      return static_cast<unsigned>(m_data.v->size());
    case GC_type::MAP:         return static_cast<unsigned>(m_data.m->size());
    case GC_type::NOTYPE:      return 0;
    }
    return 0;
  }

  unsigned
  GenericContainer::num_rows() const {
    switch (m_data_type) {
    case GC_type::POINTER:
    case GC_type::BOOL:
    case GC_type::INTEGER:
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
    case GC_type::VECTOR:      return 1;
    case GC_type::MAT_INTEGER: return m_data.m_i->num_rows();
    case GC_type::MAT_LONG:    return m_data.m_l->num_rows();
    case GC_type::MAT_REAL:    return m_data.m_r->num_rows();
    case GC_type::MAT_COMPLEX: return m_data.m_c->num_rows();
    case GC_type::MAP:         return 1;
    case GC_type::NOTYPE:      return 0;
    }
    return 0;
  }

  unsigned
  GenericContainer::num_cols() const {
    switch (m_data_type) {
    case GC_type::POINTER:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::COMPLEX:
    case GC_type::STRING:      return 1;
    case GC_type::VEC_POINTER: return static_cast<unsigned>(m_data.v_p->size());
    case GC_type::VEC_BOOL:    return static_cast<unsigned>(m_data.v_b->size());
    case GC_type::VEC_INTEGER: return static_cast<unsigned>(m_data.v_i->size());
    case GC_type::VEC_LONG:    return static_cast<unsigned>(m_data.v_l->size());
    case GC_type::VEC_REAL:    return static_cast<unsigned>(m_data.v_r->size());
    case GC_type::VEC_COMPLEX: return static_cast<unsigned>(m_data.v_c->size());
    case GC_type::VEC_STRING:  return static_cast<unsigned>(m_data.v_s->size());

    case GC_type::MAT_INTEGER: return m_data.m_i->num_cols();
    case GC_type::MAT_LONG:    return m_data.m_l->num_cols();
    case GC_type::MAT_REAL:    return m_data.m_r->num_cols();
    case GC_type::MAT_COMPLEX: return m_data.m_c->num_cols();

    case GC_type::VECTOR:      return static_cast<unsigned>(m_data.v->size());
    case GC_type::MAP:         return static_cast<unsigned>(m_data.m->size());
    case GC_type::NOTYPE:      return 0;
    }
    return 0;
  }

  int
  GenericContainer::ck( TypeAllowed const tp ) const {
    if ( tp == m_data_type     ) return 0; // ok
    if ( tp == GC_type::NOTYPE ) return 1; //
    return 2;
  }

  void
  GenericContainer::ck( string_view const where, TypeAllowed const tp ) const {
    GC_ASSERT(
      tp == m_data_type,
      where
        << " bad data type, expect: " << to_string(tp)
        << " but data stored is of type: " << to_string(m_data_type)
    )
  }

  void
  GenericContainer::ck_or_set( string_view const where, TypeAllowed const tp ) {
    if ( m_data_type == GC_type::NOTYPE ) m_data_type = tp;
    else                                  ck(where,tp);
  }

  /*
   //      _    _ _                 _
   //     / \  | | | ___   ___ __ _| |_ ___
   //    / _ \ | | |/ _ \ / __/ _` | __/ _ \
   //   / ___ \| | | (_) | (_| (_| | ||  __/
   //  /_/   \_\_|_|\___/ \___\__,_|\__\___|
   */
  void
  GenericContainer::allocate_string() {
    if ( m_data_type != GC_type::STRING ) {
      clear();
      m_data_type = GC_type::STRING;
      m_data.s    = new string_type;
    }
  }

  void
  GenericContainer::allocate_complex() {
    if ( m_data_type != GC_type::COMPLEX ) {
      clear();
      m_data_type = GC_type::COMPLEX;
      m_data.c    = new complex_type;
    }
  }

  void
  GenericContainer::allocate_vec_pointer( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_POINTER ) {
      clear();
      m_data_type = GC_type::VEC_POINTER;
      m_data.v_p  = new vec_pointer_type();
    }
    if ( sz > 0 ) m_data.v_p->resize( sz );
  }

  GenericContainer &
  GenericContainer::free_pointer() {
    GC_ASSERT(
      GC_type::POINTER == m_data_type ||
      GC_type::NOTYPE  == m_data_type,
      " free_pointer() bad data type\n"
      "expect: " << to_string(GC_type::POINTER) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    m_data.p = nullptr;
    m_data_type = GC_type::NOTYPE;
    return *this;
  }

  void
  GenericContainer::allocate_vec_bool( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_BOOL ) {
      clear();
      m_data_type = GC_type::VEC_BOOL;
      m_data.v_b  = new vec_bool_type();
    }
    if ( sz > 0 ) m_data.v_b->resize( sz );
  }

  void
  GenericContainer::allocate_vec_int( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_INTEGER ) {
      clear();
      m_data_type = GC_type::VEC_INTEGER;
      m_data.v_i  = new vec_int_type();
    }
    if ( sz > 0 ) m_data.v_i->resize( sz );
  }

  void
  GenericContainer::allocate_vec_long( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_LONG ) {
      clear();
      m_data_type = GC_type::VEC_LONG;
      m_data.v_l  = new vec_long_type();
    }
    if ( sz > 0 ) m_data.v_l->resize( sz );
  }

  void
  GenericContainer::allocate_vec_real( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_REAL ) {
      clear();
      m_data_type = GC_type::VEC_REAL;
      m_data.v_r  = new vec_real_type();
    }
    if ( sz > 0 ) m_data.v_r->resize( sz );
  }

  void
  GenericContainer::allocate_vec_complex( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_COMPLEX ) {
      clear();
      m_data_type = GC_type::VEC_COMPLEX;
      m_data.v_c  = new vec_complex_type();
    }
    if ( sz > 0 ) m_data.v_c->resize( sz );
  }

  void
  GenericContainer::allocate_mat_int( unsigned const nr, unsigned const nc ) {
    if ( m_data_type != GC_type::MAT_INTEGER ) {
      clear();
      m_data_type = GC_type::MAT_INTEGER;
      m_data.m_i  = new mat_int_type( nr, nc );
    } else {
      m_data.m_i->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_long( unsigned const nr, unsigned const nc ) {
    if ( m_data_type != GC_type::MAT_LONG ) {
      clear();
      m_data_type = GC_type::MAT_LONG;
      m_data.m_l  = new mat_long_type( nr, nc );
    } else {
      m_data.m_l->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_real( unsigned const nr, unsigned const nc ) {
    if ( m_data_type != GC_type::MAT_REAL ) {
      clear();
      m_data_type = GC_type::MAT_REAL;
      m_data.m_r  = new mat_real_type( nr, nc );
    } else {
      m_data.m_r->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_complex( unsigned const nr, unsigned const nc ) {
    if ( m_data_type != GC_type::MAT_COMPLEX ) {
      clear();
      m_data_type = GC_type::MAT_COMPLEX;
      m_data.m_c  = new mat_complex_type( nr, nc );
    } else {
      m_data.m_c->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_vec_string( unsigned const sz ) {
    if ( m_data_type != GC_type::VEC_STRING ) {
      clear();
      m_data_type = GC_type::VEC_STRING;
      m_data.v_s  = new vec_string_type();
    }
    if ( sz > 0 ) m_data.v_s->resize( sz );
  }

  void
  GenericContainer::allocate_vector( unsigned const sz ) {
    if ( m_data_type != GC_type::VECTOR ) {
      clear();
      m_data_type = GC_type::VECTOR;
      m_data.v    = new vector_type();
    }
    if ( sz > 0 ) m_data.v->resize( sz );
  }

  void
  GenericContainer::allocate_map() {
    if ( m_data_type != GC_type::MAP ) {
      clear();
      m_data_type = GC_type::MAP;
      m_data.m    = new map_type();
    }
  }

  /*
  //   ____       _
  //  / ___|  ___| |_
  //  \___ \ / _ \ __|
  //   ___) |  __/ |_
  //  |____/ \___|\__|
  */

  pointer_type &
  GenericContainer::set_pointer( pointer_type const value ) {
    clear();
    m_data_type = GC_type::POINTER;
    return (m_data.p = value);
  }

  bool_type &
  GenericContainer::set_bool( bool_type const value ) {
    clear();
    m_data_type = GC_type::BOOL;
    return (m_data.b = value);
  }

  int_type &
  GenericContainer::set_int( int_type const value ) {
    clear();
    m_data_type = GC_type::INTEGER;
    return (m_data.i = value);
  }

  long_type &
  GenericContainer::set_long( long_type const value ) {
    clear();
    m_data_type = GC_type::LONG;
    return (m_data.l = value);
  }

  real_type &
  GenericContainer::set_real( real_type const value ) {
    clear();
    m_data_type = GC_type::REAL;
    return (m_data.r = value);
  }

  complex_type &
  GenericContainer::set_complex( complex_type const & value ) {
    clear();
    m_data_type = GC_type::COMPLEX;
    m_data.c    = new complex_type;
    return (*m_data.c=value);
  }

  complex_type &
  GenericContainer::set_complex( real_type const re, real_type const im ) {
    clear();
    m_data_type = GC_type::COMPLEX;
    m_data.c    = new complex_type(re,im);
    return *m_data.c;
  }

  string_type &
  GenericContainer::set_string( string_view const value ) {
    allocate_string();
    return (*m_data.s = value);
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( unsigned const sz ) {
    allocate_vec_pointer( sz );
    return *m_data.v_p;
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( vec_pointer_type const & v ) {
    allocate_vec_pointer( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_p->begin() );
    return *m_data.v_p;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( unsigned const sz ) {
    allocate_vec_bool( sz ); return *m_data.v_b;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( vec_bool_type const & v ) {
    allocate_vec_bool( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_b->begin() );
    return *m_data.v_b;
  }

  vec_int_type &
  GenericContainer::set_vec_int( unsigned const sz ) {
    allocate_vec_int( sz );
    return *m_data.v_i;
  }

  vec_int_type &
  GenericContainer::set_vec_int( vec_int_type const & v ) {
    allocate_vec_int( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_i->begin() );
    return *m_data.v_i;
  }

  vec_long_type &
  GenericContainer::set_vec_long( unsigned const sz ) {
    allocate_vec_long( sz );
    return *m_data.v_l;
  }

  vec_long_type &
  GenericContainer::set_vec_long( vec_long_type const & v ) {
    allocate_vec_long( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_l->begin() );
    return *m_data.v_l;
  }

  vec_real_type &
  GenericContainer::set_vec_real( unsigned const sz ) {
    allocate_vec_real( sz );
    return *m_data.v_r;
  }

  vec_real_type &
  GenericContainer::set_vec_real( vec_real_type const & v ) {
    allocate_vec_real( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_r->begin() );
    return *m_data.v_r;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( unsigned const sz ) {
    allocate_vec_complex( sz );
    return *m_data.v_c;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( vec_complex_type const & v ) {
    allocate_vec_complex( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_c->begin() );
    return *m_data.v_c;
  }

  mat_int_type &
  GenericContainer::set_mat_int( unsigned const nr, unsigned const nc ) {
    allocate_mat_int( nr, nc );
    return *m_data.m_i;
  }

  mat_int_type &
  GenericContainer::set_mat_int( mat_int_type const & m ) {
    allocate_mat_int( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), m_data.m_i->begin() );
    return *m_data.m_i;
  }

  mat_long_type &
  GenericContainer::set_mat_long( unsigned const nr, unsigned const nc ) {
    allocate_mat_long( nr, nc );
    return *m_data.m_l;
  }

  mat_long_type &
  GenericContainer::set_mat_long( mat_long_type const & m ) {
    allocate_mat_long( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), m_data.m_l->begin() );
    return *m_data.m_l;
  }

  mat_real_type &
  GenericContainer::set_mat_real( unsigned const nr, unsigned const nc ) {
    allocate_mat_real( nr, nc );
    return *m_data.m_r;
  }

  mat_real_type &
  GenericContainer::set_mat_real( mat_real_type const & m ) {
    allocate_mat_real( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), m_data.m_r->begin() );
    return *m_data.m_r;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( unsigned const nr, unsigned const nc ) {
    allocate_mat_complex( nr, nc );
    return *m_data.m_c;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( mat_complex_type const & m ) {
    allocate_mat_complex( m.num_rows(), m.num_cols() );
    std::copy( m.begin(), m.end(), m_data.m_c->begin() );
    return *m_data.m_c;
  }

  vec_string_type &
  GenericContainer::set_vec_string( unsigned const sz ) {
    allocate_vec_string( sz );
    return *m_data.v_s;
  }

  vec_string_type &
  GenericContainer::set_vec_string( vec_string_type const & v ) {
    allocate_vec_string( static_cast<unsigned>( v.size() ) );
    std::copy( v.begin(), v.end(), m_data.v_s->begin() );
    return *m_data.v_s;
  }

  vector_type &
  GenericContainer::set_vector( unsigned const sz ) {
    allocate_vector( sz );
    return *m_data.v;
  }

  map_type &
  GenericContainer::set_map() {
    allocate_map();
    return *m_data.m;
  }

  /*
  //   ____            _
  //  |  _ \ _   _ ___| |__
  //  | |_) | | | / __| '_ \
  //  |  __/| |_| \__ \ | | |
  //  |_|    \__,_|___/_| |_|
  */
  void
  GenericContainer::push_bool( bool const b ) const {
    if ( m_data_type == GC_type::VEC_BOOL ) {
      m_data.v_b->push_back( b );
    } else if ( m_data_type == GC_type::VEC_INTEGER ) {
      m_data.v_i->emplace_back( b ? 1 : 0 );
    } else if ( m_data_type == GC_type::VEC_LONG ) {
      m_data.v_l->emplace_back( b ? 1 : 0 );
    } else if ( m_data_type == GC_type::VEC_REAL ) {
      m_data.v_r->emplace_back( b ? 1 : 0 );
    } else if ( m_data_type == GC_type::VEC_COMPLEX ) {
      m_data.v_c->emplace_back( b ? 1 : 0, 0 );
    } else if ( m_data_type == GC_type::VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_bool( b );
    } else {
      GC_DO_ERROR( "push_bool, bad data stored: " << get_type_name() )
    }
  }

  void
  GenericContainer::push_int( int_type const i ) {
    if ( m_data_type == GC_type::VEC_INTEGER ) {
      m_data.v_i->emplace_back( i );
    } else if ( m_data_type == GC_type::VEC_LONG ) {
      m_data.v_l->emplace_back( static_cast<long_type>(i) );
    } else if ( m_data_type == GC_type::VEC_REAL ) {
      m_data.v_r->emplace_back( i );
    } else if ( m_data_type == GC_type::VEC_COMPLEX ) {
      m_data.v_c->emplace_back( static_cast<real_type>(i), 0);
    } else if ( m_data_type == GC_type::VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_int( i );
    } else {
      if ( m_data_type != GC_type::VEC_INTEGER ) promote_to_vec_int();
      m_data.v_i->emplace_back( i );
    }
  }

  void
  GenericContainer::push_long( long_type const l ) {
    if ( m_data_type == GC_type::VEC_LONG ) {
      m_data.v_l->emplace_back( l );
    } else if ( m_data_type == GC_type::VEC_REAL ) {
      m_data.v_r->emplace_back( static_cast<real_type>(l) );
    } else if ( m_data_type == GC_type::VEC_COMPLEX ) {
      m_data.v_c->emplace_back( static_cast<real_type>(l), 0 );
    } else if ( m_data_type == GC_type::VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_long( l );
    } else {
      if ( m_data_type != GC_type::VEC_LONG ) promote_to_vec_long();
      m_data.v_l->emplace_back( l );
    }
  }

  void
  GenericContainer::push_real( real_type const r ) {
    if ( m_data_type == GC_type::VEC_REAL ) {
      m_data.v_r->emplace_back( r );
    } else if ( m_data_type == GC_type::VEC_COMPLEX ) {
      m_data.v_c->emplace_back( r, 0 );
    } else if ( m_data_type == GC_type::VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_real( r );
    } else {
      if ( m_data_type != GC_type::VEC_REAL ) promote_to_vec_real();
      m_data.v_r->emplace_back( r );
    }
  }

  void
  GenericContainer::push_complex( complex_type & c ) {
    if ( m_data_type == GC_type::VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_complex( c );
    } else {
      if ( m_data_type != GC_type::VEC_COMPLEX ) promote_to_vec_complex();
      m_data.v_c->emplace_back( c );
    }
  }

  void
  GenericContainer::push_complex( real_type const re, real_type const im ) {
    complex_type tmp( re, im );
    push_complex( tmp );
  }

  void
  GenericContainer::push_string( string_view const s ) {
    if ( m_data_type != GC_type::VEC_STRING ) promote_to_vector();
    if ( m_data_type == GC_type::VEC_STRING ) {
      m_data.v_s->emplace_back( s );
    } else {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_string( s );
    }
  }

  /*
  //    ____      _
  //   / ___| ___| |_
  //  | |  _ / _ \ __|
  //  | |_| |  __/ |_
  //   \____|\___|\__|
  */

  void *
  GenericContainer::get_pvoid( string_view const where ) const {
    ck(where,GC_type::POINTER);
    return m_data.p;
  }

  void **
  GenericContainer::get_ppvoid( string_view const where ) const {
    ck(where,GC_type::POINTER);
    return const_cast<void **>(&m_data.p);
  }

  int_type const *
  GenericContainer::get_int_pointer() const {
    switch (m_data_type) {
    case GC_type::INTEGER:
      return &m_data.i;
    case GC_type::VEC_INTEGER:
      return m_data.v_i->data();
    case GC_type::MAT_INTEGER:
      return m_data.m_i->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `int_type const*'"
      )
    }
    return nullptr;
  }

  int_type *
  GenericContainer::get_int_pointer() {
    switch (m_data_type) {
    case GC_type::INTEGER:
      return &m_data.i;
    case GC_type::VEC_INTEGER:
      return m_data.v_i->data();
    case GC_type::MAT_INTEGER:
      return m_data.m_i->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `int_type*'"
      )
    }
    return nullptr;
  }

  long_type const *
  GenericContainer::get_long_pointer() const {
    switch (m_data_type) {
    case GC_type::LONG:
      return &m_data.l;
    case GC_type::VEC_LONG:
      return m_data.v_l->data();
    case GC_type::MAT_LONG:
      return m_data.m_l->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_REAL:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_long_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `long_type const*'"
      )
    }
    return nullptr;
  }

  long_type *
  GenericContainer::get_long_pointer() {
    switch (m_data_type) {
    case GC_type::LONG:
      return &m_data.l;
    case GC_type::VEC_LONG:
      return m_data.v_l->data();
    case GC_type::MAT_LONG:
      return m_data.m_l->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_REAL:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_REAL:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_long_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `long_type*'"
      )
    }
    return nullptr;
  }

  real_type const *
  GenericContainer::get_real_pointer() const {
    switch (m_data_type) {
    case GC_type::REAL:
      return &m_data.r;
    case GC_type::VEC_REAL:
      return m_data.v_r->data();
    case GC_type::MAT_REAL:
      return m_data.m_r->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_real_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `real_type cont *'"
      )
    }
    return nullptr;
  }

  real_type *
  GenericContainer::get_real_pointer() {
    switch (m_data_type) {
    case GC_type::REAL:
      return &m_data.r;
    case GC_type::VEC_REAL:
      return m_data.v_r->data();
    case GC_type::MAT_REAL:
      return m_data.m_r->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::COMPLEX:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_COMPLEX:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_real_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `real_type*'"
      )
    }
    return nullptr;
  }

  complex_type const *
  GenericContainer::get_complex_pointer() const {
    switch (m_data_type) {
    case GC_type::COMPLEX:
      return m_data.c;
    case GC_type::VEC_COMPLEX:
      return m_data.v_c->data();
    case GC_type::MAT_COMPLEX:
      return m_data.m_c->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `complex_type const*'"
      )
    }
    return nullptr;
  }

  complex_type *
  GenericContainer::get_complex_pointer() {
    switch (m_data_type) {
    case GC_type::COMPLEX:
      return m_data.c;
    case GC_type::VEC_COMPLEX:
      return m_data.v_c->data();
    case GC_type::MAT_COMPLEX:
      return m_data.m_c->data();
    case GC_type::NOTYPE:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::POINTER:
    case GC_type::STRING:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_BOOL:
    case GC_type::VEC_INTEGER:
    case GC_type::VEC_LONG:
    case GC_type::VEC_REAL:
    case GC_type::MAT_INTEGER:
    case GC_type::MAT_LONG:
    case GC_type::MAT_REAL:
    case GC_type::VEC_STRING:
    case GC_type::VECTOR:
    case GC_type::MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << to_string(m_data_type) <<
        "' cannot be referred as `complex_type const*'"
      )
    }
    return nullptr;
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <>
  void
  GenericContainer::get_value( uint_type & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        where << " in get_value(...) negative `integer` value '" << m_data.i
              << "' cannot be converted into `uint_type'"
      )
      v = static_cast<unsigned>(m_data.i);
      break;
    case GC_type::LONG:
      GC_ASSERT(
        m_data.l >= 0,
        where << " in get_value(...) negative `long` value '" << m_data.l
              << "' cannot be converted into `uint_type'"
      )
      v = static_cast<unsigned>(m_data.l);
      break;
    case GC_type::REAL:
      GC_ASSERT(
        m_data.r >= 0 && isUnsigned(m_data.r),
        where << " in get_value(...) negative or fractional `real` value '" << m_data.r
              << "' cannot be converted into `uint_type'"
      )
      v = static_cast<unsigned>(m_data.r);
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()) && isUnsigned(m_data.c->real()),
        where << " in get_value(...) `complex` value = "
              << to_string(*m_data.c) << " cannot be converted into `uint_type'"
      )
      v = static_cast<unsigned>(m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << "\nbad data type: `" << to_string(m_data_type)
              << "' cannot be converted into `uint_type'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( int_type & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      v = m_data.i;
      break;
    case GC_type::LONG:
      v = static_cast<int>(m_data.l);
      break;
    case GC_type::REAL:
      GC_ASSERT(
        isInteger(m_data.r),
        where << " in get_value(...) fractional `real` value '" << m_data.r
              << "' cannot be converted into `int_type'"
      )
      v = static_cast<int>(m_data.r);
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()) && isInteger(m_data.c->real()),
        where << " in get_value(...) `complex` value = "
              << to_string(*m_data.c) << " cannot be converted into `int_type'"
      )
      v = static_cast<int>(m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " in get_value(...) bad data type: `" << to_string(m_data_type)
              << "' cannot be converted into `int_type'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( ulong_type & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        where << " in get_value(...) negative `integer` value '" << m_data.i
              << "' cannot be converted into `ulong_type'"
      )
      v = static_cast<ulong_type>(m_data.i);
      break;
    case GC_type::LONG:
      GC_ASSERT(
        m_data.l >= 0,
        where << " in get_value(...) negative `long` value '" << m_data.l
              << "' cannot be converted into `ulong_type'"
      )
      v = static_cast<ulong_type>(m_data.l);
      break;
    case GC_type::REAL:
      GC_ASSERT(
        m_data.r >= 0 && isUnsigned(m_data.r),
        where << " in get_value(...) negative or fractional `real` value '" << m_data.r
              << "' cannot be converted into `ulong_type'"
      )
      v = static_cast<ulong_type>(m_data.r);
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()) && isUnsigned(m_data.c->real()),
        where << " in get_value(...) `complex` value "
              << to_string(*m_data.c) << " cannot be converted into `ulong_type'"
      )
      v = static_cast<ulong_type>(m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " in get_value(...) bad data type: `"
              << to_string(m_data_type)
              << "' cannot be converted into `ulong_type'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( long_type & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      v = static_cast<long>(m_data.i);
      break;
    case GC_type::LONG:
      v = static_cast<long>(m_data.l);
      break;
    case GC_type::REAL:
      GC_ASSERT(
        isInteger(m_data.r),
        where << " in get_value(...) fractional `real` value '" << m_data.r
              << "' cannot be converted into `long_type'"
      )
      v = static_cast<long>(m_data.r);
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()) && isInteger(m_data.c->real()),
        where << " in get_value(...) `complex` value = "
              << to_string(*m_data.c) << " cannot be converted into `long_type'"
      )
      v = static_cast<long>(m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " in get_value(...) bad data type: `" << to_string(m_data_type)
              << "' cannot be converted into `long_type'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( float & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = static_cast<float>(m_data.b ? 1 : 0);
      break;
    case GC_type::INTEGER:
      v = static_cast<float>(m_data.i);
      break;
    case GC_type::LONG:
      v = static_cast<float>(m_data.l);
      break;
    case GC_type::REAL:
      v = static_cast<float>(m_data.r);
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()),
        where << " in get_value(...) `complex` value = "
              << to_string(*m_data.c) << " cannot be converted into `float'"
      )
      v = static_cast<float>(m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " in get_value(...) bad data type: `" << to_string(m_data_type)
              << "' cannot be converted into `float'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( double & v, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:
      v = m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      v = static_cast<double>(m_data.i);
      break;
    case GC_type::LONG:
      v = static_cast<double>(m_data.l);
      break;
    case GC_type::REAL:
      v = m_data.r;
      break;
    case GC_type::COMPLEX:
      GC_ASSERT(
        isZero0(m_data.c->imag()),
        where << " in get_value(...) `complex` value = "
              << to_string(*m_data.c) << " cannot be converted into `double'"
      )
      v = m_data.c->real();
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " in get_value(...) bad data type: `" << to_string(m_data_type)
              << "' cannot be converted into `double'"
      )
    }
  }

  #endif

  real_type
  GenericContainer::get_number( string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:    return m_data.b ? 1 : 0;
    case GC_type::INTEGER: return static_cast<real_type>(m_data.i);
    case GC_type::LONG:    return static_cast<real_type>(m_data.l);
    case GC_type::REAL:    return m_data.r;
    case GC_type::NOTYPE:
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
    case GC_type::MAP:
      GC_DO_ERROR(
        where << " get_number() type: " <<
        to_string(m_data_type) << " cannot be converted to double.\n"
      );
    }
    return 0;
  }

  complex_type
  GenericContainer::get_complex_number( string_view const where ) const {
    switch (m_data_type) {
    case GC_type::BOOL:    return {m_data.b?static_cast<real_type>(1):0,0};
    case GC_type::INTEGER: return {static_cast<real_type>(m_data.i),0};
    case GC_type::LONG:    return {static_cast<real_type>(m_data.l),0};
    case GC_type::REAL:    return {m_data.r,0};
    case GC_type::COMPLEX: return *m_data.c;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      GC_DO_ERROR( where << " get_number() type: " << to_string(m_data_type) << " cannot be converted to complex." );
    }
    return 0;
  }

  void
  GenericContainer::get_complex_number( real_type & re, real_type & im ) const {
    complex_type const tmp{ get_complex_number() };
    re = tmp.real();
    im = tmp.imag();
  }

  real_type
  GenericContainer::get_number_at( unsigned const i, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::VEC_BOOL:    return (*m_data.v_b)[i] ? 1 : 0;
    case GC_type::VEC_INTEGER: return static_cast<real_type>( (*m_data.v_i)[i] );
    case GC_type::VEC_LONG:    return static_cast<real_type>( (*m_data.v_l)[i] );
    case GC_type::VEC_REAL:    return (*m_data.v_r)[i];
    case GC_type::MAT_INTEGER: return static_cast<real_type>( (*m_data.m_i)[i] );
    case GC_type::MAT_LONG:    return static_cast<real_type>( (*m_data.m_l)[i] );
    case GC_type::MAT_REAL:    return (*m_data.m_r)[i];
    case GC_type::VECTOR:      return (*m_data.v)[i].get_number();
    case GC_type::NOTYPE:
    case GC_type::POINTER:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::COMPLEX:
    case GC_type::STRING:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_COMPLEX:
    case GC_type::VEC_STRING:
    case GC_type::MAT_COMPLEX:
    case GC_type::MAP:
      GC_DO_ERROR(
        where << "get_number_at( " << i << " ) type: " <<
        to_string(m_data_type) << " cannot be converted to double.\n"
      );
    }
    return 0;
  }

  complex_type
  GenericContainer::get_complex_number_at( unsigned const i, string_view const where ) const {
    switch (m_data_type) {
    case GC_type::VEC_BOOL:    return {static_cast<real_type>( (*m_data.v_b)[i]?1:0 ),0};
    case GC_type::VEC_INTEGER: return {static_cast<real_type>( (*m_data.v_i)[i] ),0};
    case GC_type::VEC_LONG:    return {static_cast<real_type>( (*m_data.v_l)[i] ),0};
    case GC_type::VEC_REAL:    return {(*m_data.v_r)[i],0};
    case GC_type::VEC_COMPLEX: return (*m_data.v_c)[i];
    case GC_type::MAT_INTEGER: return {static_cast<real_type>( (*m_data.m_i)[i] ),0};
    case GC_type::MAT_LONG:    return {static_cast<real_type>( (*m_data.m_l)[i] ),0};
    case GC_type::MAT_REAL:    return {(*m_data.m_r)[i],0};
    case GC_type::MAT_COMPLEX: return (*m_data.m_c)[i];
    case GC_type::VECTOR:      return (*m_data.v)[i].get_complex_number();
    case GC_type::NOTYPE:
    case GC_type::POINTER:
    case GC_type::BOOL:
    case GC_type::INTEGER:
    case GC_type::LONG:
    case GC_type::REAL:
    case GC_type::COMPLEX:
    case GC_type::STRING:
    case GC_type::VEC_POINTER:
    case GC_type::VEC_STRING:
    case GC_type::MAP:
      GC_DO_ERROR(
        where << "get_complex_number_at( " << i << " ) type: " <<
        to_string(m_data_type) << " cannot be converted to complex.\n"
      );
    }
    return 0;
  }

  void
  GenericContainer::get_complex_number_at( unsigned const i, real_type & re, real_type & im, string_view const where ) const {
    complex_type const tmp{ get_complex_number_at( i, where ) };
    re = tmp.real();
    im = tmp.imag();
  }

  bool_type &
  GenericContainer::get_bool( string_view const where ) {
    ck_or_set(where,GC_type::BOOL);
    return m_data.b;
  }

  bool_type const &
  GenericContainer::get_bool( string_view const where ) const {
    ck(where,GC_type::BOOL);
    return m_data.b;
  }

  int_type &
  GenericContainer::get_int( string_view const where ) {
    ck_or_set(where,GC_type::INTEGER);
    return m_data.i;
  }

  int_type const &
  GenericContainer::get_int( string_view const where ) const {
    ck(where,GC_type::INTEGER);
    return m_data.i;
  }

  long_type &
  GenericContainer::get_long( string_view const where ) {
    ck_or_set(where,GC_type::LONG);
    return m_data.l;
  }

  long_type const &
  GenericContainer::get_long( string_view const where ) const {
    ck(where,GC_type::LONG);
    return m_data.l;
  }

  int_type
  GenericContainer::get_as_int( string_view const where ) const {
    int_type res;
    this->get_value<int_type>( res, where );
    return res;
  }

  uint_type
  GenericContainer::get_as_uint( string_view const where ) const {
    uint_type res;
    this->get_value<uint_type>( res, where );
    return res;
  }

  long_type
  GenericContainer::get_as_long( string_view const where ) const {
    long_type res;
    this->get_value<long_type>( res, where );
    return res;
  }

  ulong_type
  GenericContainer::get_as_ulong( string_view const where ) const {
    ulong_type res;
    this->get_value<ulong_type>( res, where );
    return res;
  }

  real_type &
  GenericContainer::get_real( string_view const where ) {
    ck_or_set(where,GC_type::REAL);
    return m_data.r;
  }

  real_type const &
  GenericContainer::get_real( string_view const where ) const {
    ck(where,GC_type::REAL);
    return m_data.r;
  }

  complex_type &
  GenericContainer::get_complex( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) {
      clear();
      m_data_type = GC_type::COMPLEX;
      m_data.c    = new complex_type;
    } else {
      ck(where,GC_type::COMPLEX);
    }
    return *m_data.c;
  }

  complex_type const &
  GenericContainer::get_complex( string_view const where ) const {
    ck(where,GC_type::COMPLEX);
    return *m_data.c;
  }

  string_type &
  GenericContainer::get_string( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) {
      clear();
      m_data_type = GC_type::STRING;
      m_data.s    = new string_type("");
    } else {
      ck(where,GC_type::STRING);
    }
    return *m_data.s;
  }

  string const &
  GenericContainer::get_string( string_view const where ) const {
    ck(where,GC_type::STRING);
    return *m_data.s;
  }

  vector_type &
  GenericContainer::get_vector( string_view const where ) {
    ck(where,GC_type::VECTOR);
    return *m_data.v;
  }

  vector_type const &
  GenericContainer::get_vector( string_view const where ) const {
    ck(where,GC_type::VECTOR);
    return *m_data.v;
  }

  vec_pointer_type &
  GenericContainer::get_vec_pointer( string_view const where ) {
    ck(where,GC_type::VEC_POINTER);
    return *m_data.v_p;
  }

  vec_pointer_type const &
  GenericContainer::get_vec_pointer( string_view const where ) const {
    ck(where,GC_type::VEC_POINTER);
    return *m_data.v_p;
  }

  vec_bool_type &
  GenericContainer::get_vec_bool( string_view const where ) {
    ck(where,GC_type::VEC_BOOL);
    return *m_data.v_b;
  }

  vec_bool_type const &
  GenericContainer::get_vec_bool( string_view const where ) const {
    ck(where,GC_type::VEC_BOOL);
    return *m_data.v_b;
  }

  vec_int_type &
  GenericContainer::get_vec_int( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE   ) set_vec_int();
    if ( m_data_type == GC_type::VEC_BOOL ) promote_to_vec_int();
    ck(where,GC_type::VEC_INTEGER);
    return *m_data.v_i;
  }

  vec_int_type const &
  GenericContainer::get_vec_int( string_view const where ) const {
    ck(where,GC_type::VEC_INTEGER);
    return *m_data.v_i;
  }

  vec_long_type &
  GenericContainer::get_vec_long( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_vec_long();
    if ( m_data_type == GC_type::VEC_BOOL ||
         m_data_type == GC_type::VEC_INTEGER ) promote_to_vec_long();
    ck(where,GC_type::VEC_LONG);
    return *m_data.v_l;
  }

  vec_long_type const &
  GenericContainer::get_vec_long( string_view const where ) const {
    ck(where,GC_type::VEC_LONG);
    return *m_data.v_l;
  }

  vec_real_type &
  GenericContainer::get_vec_real( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_vec_real();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG ) promote_to_vec_real();
    ck(where,GC_type::VEC_REAL);
    return *m_data.v_r;
  }

  vec_real_type const &
  GenericContainer::get_vec_real( string_view const where ) const {
    ck(where,GC_type::VEC_REAL);
    return *m_data.v_r;
  }

  vec_complex_type &
  GenericContainer::get_vec_complex( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_vec_complex();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG    ||
         m_data_type == GC_type::VEC_REAL ) promote_to_vec_complex();
    ck(where,GC_type::VEC_COMPLEX);
    return *m_data.v_c;
  }

  vec_complex_type const &
  GenericContainer::get_vec_complex( string_view const where ) const {
    ck(where,GC_type::VEC_COMPLEX);
    return *m_data.v_c;
  }

  mat_int_type &
  GenericContainer::get_mat_int( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_mat_int();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG    ||
         m_data_type == GC_type::VEC_REAL ) promote_to_mat_int();
    ck(where,GC_type::MAT_INTEGER);
    return *m_data.m_i;
  }

  mat_int_type const &
  GenericContainer::get_mat_int( string_view const where ) const {
    ck(where,GC_type::MAT_INTEGER);
    return *m_data.m_i;
  }

  mat_long_type &
  GenericContainer::get_mat_long( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_mat_long();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG    ||
         m_data_type == GC_type::VEC_REAL    ||
         m_data_type == GC_type::MAT_INTEGER ) promote_to_mat_long();
    ck(where,GC_type::MAT_LONG);
    return *m_data.m_l;
  }

  mat_long_type const &
  GenericContainer::get_mat_long( string_view const where ) const {
    ck(where,GC_type::MAT_LONG);
    return *m_data.m_l;
  }

  mat_real_type &
  GenericContainer::get_mat_real( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_mat_real();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG    ||
         m_data_type == GC_type::VEC_REAL    ||
         m_data_type == GC_type::MAT_INTEGER ||
         m_data_type == GC_type::MAT_LONG ) promote_to_mat_real();
    ck(where,GC_type::MAT_REAL);
    return *m_data.m_r;
  }

  mat_real_type const &
  GenericContainer::get_mat_real( string_view const where ) const {
    ck(where,GC_type::MAT_REAL);
    return *m_data.m_r;
  }

  mat_complex_type &
  GenericContainer::get_mat_complex( string_view const where ) {
    if ( m_data_type == GC_type::NOTYPE ) set_mat_complex();
    if ( m_data_type == GC_type::VEC_BOOL    ||
         m_data_type == GC_type::VEC_INTEGER ||
         m_data_type == GC_type::VEC_LONG    ||
         m_data_type == GC_type::VEC_REAL    ||
         m_data_type == GC_type::MAT_REAL    ||
         m_data_type == GC_type::VEC_COMPLEX ) promote_to_mat_complex();
    ck(where,GC_type::MAT_COMPLEX);
    return *m_data.m_c;
  }

  mat_complex_type const &
  GenericContainer::get_mat_complex( string_view const where ) const {
    ck(where,GC_type::MAT_COMPLEX);
    return *m_data.m_c;
  }

  vec_string_type &
  GenericContainer::get_vec_string( string_view const where ) {
    ck(where,GC_type::VEC_STRING);
    return *m_data.v_s;
  }

  vec_string_type const &
  GenericContainer::get_vec_string( string_view const where ) const {
    ck(where,GC_type::VEC_STRING);
    return *m_data.v_s;
  }

  map_type &
  GenericContainer::get_map( string_view const where ) {
    ck(where,GC_type::MAP);
    return *m_data.m;
  }

  map_type const &
  GenericContainer::get_map( string_view const where ) const {
    ck(where,GC_type::MAP);
    return *m_data.m;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  bool_type
  GenericContainer::get_map_bool(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(string(key)).get_bool( where );
  }

  bool_type
  GenericContainer::get_map_bool( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_bool(msg);
  }

  bool_type
  GenericContainer::get_map_bool(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_bool( where );
  }

  int_type
  GenericContainer::get_map_int(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_as_int( where );
  }

  int_type
  GenericContainer::get_map_int( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_int(msg);
  }

  int_type
  GenericContainer::get_map_int(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_as_int( where );
  }

  real_type
  GenericContainer::get_map_number(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_number( where );
  }

  real_type
  GenericContainer::get_map_number( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_number(msg);
  }

  real_type
  GenericContainer::get_map_number(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_number( where );
  }

  string const &
  GenericContainer::get_map_string(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_string( where );
  }

  string const &
  GenericContainer::get_map_string( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_string(msg);
  }

  string const &
  GenericContainer::get_map_string(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_string( where );
  }

  vec_real_type const &
  GenericContainer::get_map_vec_real(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_vec_real( where );
  }

  vec_real_type const &
  GenericContainer::get_map_vec_real( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_vec_real(msg);
  }

  vec_real_type const &
  GenericContainer::get_map_vec_real(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_vec_real( where );
  }

  vec_complex_type const &
  GenericContainer::get_map_vec_complex(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_vec_complex( where );
  }

  vec_complex_type const &
  GenericContainer::get_map_vec_complex( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_vec_complex(msg);
  }

  vec_complex_type const &
  GenericContainer::get_map_vec_complex(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_vec_complex( where );
  }

  vec_string_type const &
  GenericContainer::get_map_vec_string(
    string_view const key,
    string_view const where
  ) const {
    GC_ASSERT( this->exists(key), where << " key: `" << key << "` is missing" );
    return this->m_data.m->at(key.data()).get_vec_string( where );
  }

  vec_string_type const &
  GenericContainer::get_map_vec_string( std::initializer_list<string> const args ) const {
    auto p{ this };
    string msg;
    for ( string const & key : args ) {
      msg += " -> `";
      msg += key;
      msg += '`';
      GC_ASSERT( p->exists(key), "in the map sequence: " << msg << " last key is missing" );
      p = &(*p)(key);
    }
    return p->get_vec_string(msg);
  }

  vec_string_type const &
  GenericContainer::get_map_vec_string(
    vec_string_type const & keys,
    string_view     const   where
  ) const {
    string const who{ must_exists( keys, where ) };
    return this->m_data.m->at(who).get_vec_string( where );
  }

  // --------------------------------------------------------------

  bool
  GenericContainer::exists( string_view const s ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(s.data()) };
    return iv != m_data.m->end();
  }

  bool
  GenericContainer::exists( vec_string_type const & vs ) const {
    if (m_data_type != GC_type::MAP) return false;
    return std::any_of(vs.begin(), vs.end(),
      [this](string_type const & s) { return m_data.m->find(s.data()) != m_data.m->end(); }
    );
  }
  string
  GenericContainer::must_exists( vec_string_type const & vs, string_view const where ) const {
    GC_ASSERT(
      m_data_type == GC_type::MAP,
      where
        << " bad data type, expect: " << to_string(GC_type::MAP)
        << " but data stored is of type: " << to_string(m_data_type)
    )
    for ( string_type const & s : vs ) {
      if ( auto iv{ m_data.m->find(s) }; iv != m_data.m->end() ) return s;
    }
    GC_DO_ERROR( where << " cant find keys: " << vs )
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists( string_view const field, bool & value ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    if ( iv->second.m_data_type != GC_type::BOOL ) return false;
    value = iv->second.m_data.b;
    return true;
  }

  bool
  GenericContainer::get_if_exists(
    vec_string_type const & fields,
    bool                  & value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    map_type & m = *m_data.m;
    for ( string_type const & field : fields ) {
      if ( auto iv{ m.find(field) }; iv != m.end() && iv->second.m_data_type != GC_type::BOOL ) {
        value = iv->second.m_data.b;
        return true;
      }
    }
    return false;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    int_type &        value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      value = iv->second.m_data.i;
      break;
    case GC_type::LONG:
      value = static_cast<int_type>(iv->second.m_data.l);
      break;
    case GC_type::REAL:
      if ( !isInteger(iv->second.m_data.r) ) return false;
      value = static_cast<int_type>(iv->second.m_data.r);
      break;
    case GC_type::COMPLEX:
      if ( ! ( isInteger(iv->second.m_data.c->real()) &&
               isZero0(iv->second.m_data.c->imag()) ) ) return false;
      value = static_cast<int_type>(iv->second.m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    uint_type &       value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      if ( iv->second.m_data.i < 0 ) return false;
      value = static_cast<uint_type>(iv->second.m_data.i);
      break;
    case GC_type::LONG:
      if ( iv->second.m_data.l < 0 ) return false;
      value = static_cast<uint_type>(iv->second.m_data.l);
      break;
    case GC_type::REAL:
      if ( ! isUnsigned(iv->second.m_data.r) ) return false;
      value = static_cast<uint_type>(iv->second.m_data.r);
      break;
    case GC_type::COMPLEX:
      if ( ! ( isUnsigned(iv->second.m_data.c->real()) &&
               isZero0(iv->second.m_data.c->imag()) ) ) return false;
      value = static_cast<uint_type>(iv->second.m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    long_type &       value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      value = static_cast<long_type>(iv->second.m_data.i);
      break;
    case GC_type::LONG:
      value = iv->second.m_data.l;
      break;
    case GC_type::REAL:
      if ( ! isInteger(iv->second.m_data.r) ) return false;
      value = static_cast<long_type>(iv->second.m_data.r);
      break;
    case GC_type::COMPLEX:
      if ( ! ( isInteger(iv->second.m_data.c->real()) &&
               isZero0(iv->second.m_data.c->imag()) ) ) return false;
      value = static_cast<long_type>(iv->second.m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    ulong_type &      value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_type::INTEGER:
      if ( iv->second.m_data.i < 0 ) return false;
      value = static_cast<ulong_type>(iv->second.m_data.i);
      break;
    case GC_type::LONG:
      if ( iv->second.m_data.l < 0 ) return false;
      value = static_cast<ulong_type>(iv->second.m_data.l);
      break;
    case GC_type::REAL:
      if ( ! isUnsigned(iv->second.m_data.r) ) return false;
      value = static_cast<ulong_type>(iv->second.m_data.r);
      break;
    case GC_type::COMPLEX:
      if ( ! ( isUnsigned(iv->second.m_data.c->real()) &&
               isZero0(iv->second.m_data.c->imag()) ) ) return false;
      value = static_cast<ulong_type>(iv->second.m_data.c->real());
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    real_type &       value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = static_cast<real_type>(iv->second.m_data.b ? 1 : 0);
      break;
    case GC_type::INTEGER:
      value = static_cast<real_type>(iv->second.m_data.i);
      break;
    case GC_type::LONG:
      value = static_cast<real_type>(iv->second.m_data.l);
      break;
    case GC_type::REAL:
      value = iv->second.m_data.r;
      break;
    case GC_type::COMPLEX:
      if ( ! isZero0(iv->second.m_data.c->imag()) ) return false;
      value = iv->second.m_data.c->real();
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    complex_type &    value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_type::BOOL:
      value = complex_type(iv->second.m_data.b?1:0,0);
      break;
    case GC_type::INTEGER:
      value = complex_type(static_cast<real_type>(iv->second.m_data.i),0);
      break;
    case GC_type::LONG:
      value = complex_type(static_cast<real_type>(iv->second.m_data.l),0);
      break;
    case GC_type::REAL:
      value = complex_type(iv->second.m_data.r,0);
      break;
    case GC_type::COMPLEX:
      value = *iv->second.m_data.c;
      break;
    case GC_type::NOTYPE:
    case GC_type::POINTER:
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
    case GC_type::MAP:
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------

  bool
  GenericContainer::get_if_exists(
    string_view const field,
    string_type &     value
  ) const {
    if ( m_data_type != GC_type::MAP ) return false;
    auto const iv{ m_data.m->find(field.data()) };
    if ( iv == m_data.m->end() ) return false;
    if ( iv->second.m_data_type != GC_type::STRING ) return false;
    value = *iv->second.m_data.s;
    return true;
  }

  // --------------------------------------------------------------
  bool_type
  GenericContainer::get_bool_at( unsigned const i ) {
    if ( m_data_type == GC_type::NOTYPE   ) set_vec_bool();
    if ( m_data_type == GC_type::VEC_BOOL ) {
      CHECK_RESIZE(m_data.v_b,i); // correct type, check size
      return (*m_data.v_b)[i];
    }
    if ( m_data_type != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i].set_bool(false);
  }

  bool_type
  GenericContainer::get_bool_at( unsigned const i, string_view const where ) const {
    ck(where,GC_type::VEC_BOOL);
    GC_ASSERT(
      i < m_data.v_b->size(),
      where << " get_bool_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_b)[i];
  }

  int_type &
  GenericContainer::get_int_at( unsigned const i ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_vec_int();
    else if ( m_data_type == GC_type::BOOL    ||
              m_data_type == GC_type::INTEGER ||
              m_data_type == GC_type::VEC_BOOL ) promote_to_vec_int();
    if ( m_data_type == GC_type::VEC_INTEGER ) {
      CHECK_RESIZE(m_data.v_i,i); // correct type, check size
      return (*m_data.v_i)[i];
    }
    if ( m_data_type != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i].set_int(0);
  }

  int_type const &
  GenericContainer::get_int_at( unsigned const i, string_view const where ) const {
    ck(where,GC_type::VEC_INTEGER);
    GC_ASSERT(
      i < m_data.v_i->size(),
      where << " get_int_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_i)[i];
  }

  int_type &
  GenericContainer::get_int_at( unsigned const i, unsigned const j ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ) promote_to_mat_int();
    GC_ASSERT(
      GC_type::MAT_INTEGER == m_data_type,
      "get_int_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << to_string(GC_type::MAT_INTEGER) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    return (*m_data.m_i)(i,j);
  }

  int_type const &
  GenericContainer::get_int_at( unsigned const i, unsigned const j, string_view const where ) const  {
    ck(where,GC_type::MAT_INTEGER);
    GC_ASSERT(
      i < m_data.m_i->num_rows() && j < m_data.m_i->num_cols(),
      where << " get_int_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_i)(i,j);
  }

  long_type &
  GenericContainer::get_long_at( unsigned const i ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_vec_long();
    else if ( m_data_type == GC_type::BOOL     ||
              m_data_type == GC_type::INTEGER  ||
              m_data_type == GC_type::LONG     ||
              m_data_type == GC_type::VEC_BOOL ||
              m_data_type == GC_type::VEC_INTEGER ) promote_to_vec_long();
    if ( m_data_type == GC_type::VEC_LONG ) {
      CHECK_RESIZE(m_data.v_l,i); // correct type, check size
      return (*m_data.v_l)[i];
    }
    if ( m_data_type != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i].set_long(0);
  }

  long_type const &
  GenericContainer::get_long_at( unsigned const i, string_view const where ) const {
    ck(where,GC_type::VEC_LONG);
    GC_ASSERT(
      i < m_data.v_l->size(),
      where << " get_long_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_l)[i];
  }

  long_type &
  GenericContainer::get_long_at( unsigned const i, unsigned const j ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ||
              m_data_type == GC_type::VEC_LONG ) promote_to_mat_long();
    GC_ASSERT(
      GC_type::MAT_LONG == m_data_type,
      "get_long_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << to_string(GC_type::MAT_LONG) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    return (*m_data.m_l)(i,j);
  }

  long_type const &
  GenericContainer::get_long_at( unsigned const i, unsigned const j, string_view const where ) const  {
    ck(where,GC_type::MAT_LONG);
    GC_ASSERT(
      i < m_data.m_l->num_rows() && j < m_data.m_l->num_cols(),
      where << " get_long_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_l)(i,j);
  }

  real_type &
  GenericContainer::get_real_at( unsigned const i ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_vec_real();
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ||
              m_data_type == GC_type::VEC_LONG ) promote_to_vec_real();
    if ( m_data_type == GC_type::VEC_REAL ) {
      CHECK_RESIZE(m_data.v_r,i); // correct type, check size
      return (*m_data.v_r)[i];
    }
    if ( m_data_type != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i].set_real(0);
  }

  real_type const &
  GenericContainer::get_real_at( unsigned const i, string_view const where ) const  {
    ck(where,GC_type::VEC_REAL);
    GC_ASSERT(
      i < m_data.v_r->size(),
      where << " get_real_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_r)[i];
  }

  real_type &
  GenericContainer::get_real_at( unsigned const i, unsigned const j ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ||
              m_data_type == GC_type::VEC_LONG ||
              m_data_type == GC_type::VEC_REAL ) promote_to_mat_real();
    GC_ASSERT(
      GC_type::MAT_REAL == m_data_type,
      "get_real_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << to_string(GC_type::MAT_REAL) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    return (*m_data.m_r)(i,j);
  }

  real_type const &
  GenericContainer::get_real_at( unsigned const i, unsigned const j, string_view const where ) const  {
    ck(where,GC_type::MAT_REAL);
    GC_ASSERT(
      i < m_data.v_r->size(),
      where << " get_real_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_r)(i,j);
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned const i ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_vec_complex();
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ||
              m_data_type == GC_type::VEC_LONG    ||
              m_data_type == GC_type::VEC_REAL ) promote_to_vec_complex();
    if ( m_data_type == GC_type::VEC_COMPLEX ) {
      CHECK_RESIZE(m_data.v_c,i); // correct type, check size
      return (*m_data.v_c)[i];
    }
    if ( m_data_type != GC_type::VECTOR ) promote_to_vector();
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i].set_complex(0,0);
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned const i, string_view const where ) const  {
    ck(where,GC_type::VEC_COMPLEX);
    GC_ASSERT(
      i < m_data.v_c->size(),
      where << " get_complex_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_c)[i];
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned const i, unsigned const j ) {
    if      ( m_data_type == GC_type::NOTYPE ) set_mat_complex(i,j);
    else if ( m_data_type == GC_type::VEC_BOOL    ||
              m_data_type == GC_type::VEC_INTEGER ||
              m_data_type == GC_type::VEC_LONG    ||
              m_data_type == GC_type::VEC_REAL    ||
              m_data_type == GC_type::VEC_COMPLEX ||
              m_data_type == GC_type::MAT_REAL ) promote_to_mat_complex();
    GC_ASSERT(
      GC_type::MAT_COMPLEX == m_data_type,
      "get_complex_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << to_string(GC_type::MAT_COMPLEX) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    return (*m_data.m_c)(i,j);
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned const i, unsigned const j, string_view const where ) const  {
    ck(where,GC_type::MAT_COMPLEX);
    return (*m_data.m_c)(i,j);
  }

  string_type &
  GenericContainer::get_string_at( unsigned const i ) {
    if ( m_data_type == GC_type::NOTYPE ) set_vec_string();
    if ( m_data_type == GC_type::VEC_STRING ) {
      CHECK_RESIZE(m_data.v_s,i);
      return (*m_data.v_s)[i];
    }
    promote_to_vector();
    return (*this)[i].set_string("");
  }

  string const &
  GenericContainer::get_string_at( unsigned const i, string_view const where ) const {
    ck(where,GC_type::VEC_STRING);
    GC_ASSERT(
      i < m_data.v_s->size(),
      where << " get_string_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_s)[i];
  }

  GenericContainer &
  GenericContainer::get_gc_at( unsigned const i )
  { return (*this)[i]; }

  GenericContainer const &
  GenericContainer::get_gc_at( unsigned const i, string_view const where ) const {
    return (*this)(i,where);
  }

  /*
  //   _        __
  //  (_)_ __  / _| ___
  //  | | '_ \| |_ / _ \
  //  | | | | |  _| (_) |
  //  |_|_| |_|_|  \___/
  */

  GenericContainer const &
  GenericContainer::info( ostream_type & stream ) const {
    switch ( m_data_type ) {
    case GC_type::NOTYPE:
      stream << "GenericContainer: No data stored\n";
      break;
    case GC_type::POINTER:
      stream << "Generic pointer: " << m_data.p << '\n';
      break;
    case GC_type::BOOL:
      stream << "Boolean: " << (m_data.b?"true":"false") << '\n';
      break;
    case GC_type::INTEGER:
      stream << "Integer: " << m_data.i << '\n';
      break;
    case GC_type::LONG:
      stream << "Long: " << m_data.l << '\n';
      break;
    case GC_type::REAL:
      stream << "Floating Point: " << m_data.r << '\n';
      break;
    case GC_type::COMPLEX:
      stream << "Complex Floating Point: " << to_string(*m_data.c) << '\n';
      break;
    case GC_type::STRING:
      stream << "String: " << m_data.s->c_str() << '\n';
      break;
    case GC_type::VEC_POINTER:
      stream << "Vector of generic pointer of size " << m_data.v_p->size() << '\n';
      break;
    case GC_type::VEC_BOOL:
      stream << "Vector of boolean of size " << m_data.v_b->size() << '\n';
      break;
    case GC_type::VEC_INTEGER:
      stream << "Vector of integer of size " << m_data.v_i->size() << '\n';
      break;
    case GC_type::VEC_LONG:
      stream << "Vector of long integer of size " << m_data.v_l->size() << '\n';
      break;
    case GC_type::VEC_REAL:
      stream << "Vector of floating point number of size " << m_data.v_r->size() << '\n';
      break;
    case GC_type::VEC_COMPLEX:
      stream << "Vector of complex floating point number of size " << m_data.v_c->size() << '\n';
      break;
    case GC_type::VEC_STRING:
      stream << "Vector of string of size " << m_data.v_s->size() << '\n';
      break;
    case GC_type::MAT_INTEGER:
      m_data.m_i->info(stream);
      break;
    case GC_type::MAT_LONG:
      m_data.m_l->info(stream);
      break;
    case GC_type::MAT_REAL:
      m_data.m_r->info(stream);
      break;
    case GC_type::MAT_COMPLEX:
      m_data.m_c->info(stream);
      break;
    case GC_type::VECTOR:
      stream << "Vector of generic data type of size " << m_data.v->size() << '\n';
      break;
    case GC_type::MAP:
      stream << "Map\n";
      break;
    //default:
    //  stream << "Type N. " << _data_type << " not recognized\n";
    //  break;
    }
    return *this;
  }

  /*
  //                              _               __ __
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __  | _|_ |
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | | | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | | | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | | | |
  //        |_|                                  |__|__|
  */

  GenericContainer &
  GenericContainer::operator [] ( unsigned const i ) {
    switch ( ck( GC_type::VECTOR ) ) {
      case 0: break; // data present
      default: set_vector(); // data must be allocated;
    }
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i];
  }

  GenericContainer const &
  GenericContainer::operator [] ( unsigned const i ) const {
    GC_ASSERT(
      GC_type::VECTOR == m_data_type,
      "operator [] integer argument = " << i <<
      "\nexpect: " << to_string(GC_type::VECTOR) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    GC_ASSERT(
      i < m_data.v->size(),
      "operator [] const, index " << i << " out of range"
    )
    return (*m_data.v)[i];
  }

  /*
  //                              _                ____
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __   / /\ \
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | |  | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | |  | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | |  | |
  //        |_|                                   \_\/_/
  //
  */

  GenericContainer &
  GenericContainer::operator () ( unsigned const i, string_view const where ) {
    ck(where,GC_type::VECTOR);
    GC_ASSERT(
      i < m_data.v->size(),
      where << " operator () const, index " << i << " out of range"
    )
    return (*m_data.v)[i];
  }

  GenericContainer const &
  GenericContainer::operator () ( unsigned const i, string_view const where ) const {
    ck(where,GC_type::VECTOR);
    GC_ASSERT(
      i < m_data.v->size(),
      where << " operator () const, index " << i << " out of range"
    )
    return (*m_data.v)[i];
  }

  GenericContainer &
  GenericContainer::operator () ( string_view s, string_view const where ) {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      where << " operator (), with string argument ``" << s << "''"
      "\nexpect: " << to_string(GC_type::MAP) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    auto iv{ m_data.m->find(s.data()) };
    if ( iv == m_data.m->end() ) {
      GC_DO_ERROR(
        where << " operator(): Cannot find key '" << s <<
        "'!\npossibile keys: " << get_keys()
      )
    }
    return iv->second;
  }

  GenericContainer const &
  GenericContainer::operator () ( string_view s, string_view const where ) const {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      where << "\noperator() const, with string argument ``" << s << "''"
      "\nexpect: " << to_string(GC_type::MAP) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    auto iv{ m_data.m->find(s.data()) };
    if ( iv == m_data.m->end() ) {
      GC_DO_ERROR(
        where << "\noperator() const: Cannot find key '" << s <<
        "'!\npossibile keys: " << get_keys()
      )
    }
    return iv->second;
  }

  // ------------

  GenericContainer &
  GenericContainer::operator () ( vec_string_type const & vs, string_view const where ) {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      where << " operator (), with vector of string argument\n"
      "expect: " << to_string(GC_type::MAP) <<
      " but data stored is of type: " << to_string(m_data_type)
    )
    map_type & m = *m_data.m;
    for ( string_type const & s : vs ) {
      if ( auto iv{ m.find(s) }; iv != m.end() ) return iv->second;
    }
    GC_DO_ERROR( where << " operator(): Cannot find the key!" );
    //return *this;
  }

  GenericContainer const &
  GenericContainer::operator () ( vec_string_type const & vs, string_view const where ) const {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      where << " operator (), with vector of string argument\n"
      "expect: " << to_string(GC_type::MAP) <<
      " but data stored is of type: " << to_string(m_data_type)
    )
    map_type const & m = *m_data.m;
    for ( string_type const & s : vs ) {
      if ( auto iv{ m.find(s) }; iv != m.end() ) return iv->second;
    }
    GC_DO_ERROR( where << "\noperator(): Cannot find the key!" );
    //return *this;
  }

  /*
  //                              _               __ __
  //    ___  _ __   ___ _ __ __ _| |_ ___  _ __  | _|_ |
  //   / _ \| '_ \ / _ \ '__/ _` | __/ _ \| '__| | | | |
  //  | (_) | |_) |  __/ | | (_| | || (_) | |    | | | |
  //   \___/| .__/ \___|_|  \__,_|\__\___/|_|    | | | |
  //        |_|                                  |__|__|
  */

  GenericContainer &
  GenericContainer::operator [] ( string_view const s ) {
    if ( ck( GC_type::MAP ) != 0 ) set_map(); // if not data present allocate!
    return (*m_data.m)[s.data()];
  }

  GenericContainer const &
  GenericContainer::operator [] ( string_view const s ) const {
    GC_ASSERT(
      GC_type::MAP == m_data_type,
      "operator [] string argument ``" << s << "''"
      "\nexpect: " << to_string(GC_type::MAP) <<
      "\nbut data stored is of type: " << to_string(m_data_type)
    )
    return (*m_data.m)[s.data()];
  }

  /*
  //              _       _
  //   _ __  _ __(_)_ __ | |_
  //  | '_ \| '__| | '_ \| __|
  //  | |_) | |  | | | | | |_
  //  | .__/|_|  |_|_| |_|\__|
  //  |_|
  */

  void
  GenericContainer::dump(
    ostream_type &    stream,
    string_view const prefix,
    string_view const indent
  ) const {

    switch (m_data_type) {
    case GC_type::NOTYPE:
      stream << prefix << "null\n";
      break;
    case GC_type::POINTER:
      stream << prefix << std::hex << std::showbase << reinterpret_cast<uintptr_t>(m_data.p) << '\n';
      break;
    case GC_type::BOOL:
      stream << prefix << (m_data.b?"true":"false") << '\n';
      break;
    case GC_type::INTEGER:
      stream << prefix << m_data.i << '\n';
      break;
    case GC_type::LONG:
      stream << prefix << m_data.l << '\n';
      break;
    case GC_type::REAL:
      stream << prefix << m_data.r << '\n';
      break;
    case GC_type::COMPLEX:
      stream << prefix << to_string(*m_data.c) << '\n';
      break;
    case GC_type::STRING:
      stream << prefix << "\"" << *m_data.s << "\"\n";
      break;
    case GC_type::VEC_POINTER:
      { vec_pointer_type const & v{*m_data.v_p};
        for ( vec_pointer_type::size_type i{0}; i < v.size(); ++i )
          stream << prefix << "vec_pointer(" << i << "): "
                 << std::hex << std::showbase << reinterpret_cast<uintptr_t>(v[i]) << '\n';
      }
      break;
    case GC_type::VEC_BOOL:
      { vec_bool_type const & v{*m_data.v_b};
        stream << prefix << v << '\n'; }
      break;
    case GC_type::VEC_INTEGER:
      { vec_int_type const & v{*m_data.v_i};
        stream << prefix << v << '\n'; }
      break;
    case GC_type::VEC_LONG:
      { vec_long_type const & v{*m_data.v_l};
        stream << prefix << v << '\n'; }
      break;
    case GC_type::VEC_REAL:
      { vec_real_type const & v{*m_data.v_r};
        stream << prefix << v << '\n'; }
      break;
    case GC_type::VEC_COMPLEX:
      { vec_complex_type const & v{*m_data.v_c};
        stream << prefix << v << '\n'; }
      break;
    case GC_type::MAT_INTEGER:
      { mat_int_type const & m{*m_data.m_i};
        stream << m; }
      break;
    case GC_type::MAT_LONG:
      { mat_long_type const & m{*m_data.m_l};
        stream << m; }
      break;
    case GC_type::MAT_REAL:
      { mat_real_type const & m{*m_data.m_r};
        stream << m; }
      break;
    case GC_type::MAT_COMPLEX:
      { mat_complex_type const & m{*m_data.m_c};
        stream << m; }
      break;
    case GC_type::VEC_STRING:
      { vec_string_type const & v{*m_data.v_s};
        for ( vec_string_type::size_type i{0}; i < v.size(); ++i )
          stream << prefix << i << ": \"" << v[i] << "\"\n";
      }
      break;

    case GC_type::VECTOR:
      { vector_type const & v{*m_data.v};
        string prefix2{prefix}; prefix2 += indent;
        for ( vector_type::size_type i{0}; i < v.size(); ++i ) {
          GenericContainer const & vi{v[i]};
          if ( vi.simple_data() ||
               ( vi.simple_vec_data() && vi.get_num_elements() <= 10 ) ) {
            stream << prefix << i << ": ";
            vi.dump(stream,"");
          } else {
            stream << prefix << i << ":\n";
            vi.dump(stream,prefix2);
          }
        }
      }
      break;
    case GC_type::MAP:
      { map_type const & m{*m_data.m};
        string prefix2{prefix}; prefix2 += indent;
        for ( const auto &[fst, snd] : m ) {
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          string_type matches[4];
          if ( int const pcreExecRet{pcre_for_GC.exec( fst, matches )}; pcreExecRet == 4 ) {
            string_type header = matches[3]; // header
            // found formatting
            if ( snd.simple_data() ) {
              stream << prefix << header << ": ";
              snd.dump(stream,"");
            } else {
              if ( matches[1].length() > 1 ) stream << '\n'; // double ## --> add nel line
              stream << prefix << header;
              if ( matches[2].length() > 0 ) {
                stream << '\n' << prefix;
                char const fmt{ matches[2][0] }; // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt; // underline header
              } else {
                stream << ':';
              }
              stream << '\n';
              snd.dump(stream,prefix2);
            }
          } else {
            string_type header = pcreExecRet == 3 ? matches[3] : fst;
            if ( snd.simple_data() ) {
              stream << prefix << header << ": ";
              snd.dump(stream,"");
            } else {
              stream << prefix << header << ":\n";
              snd.dump(stream,prefix2);
            }
          }
        }
      }
      break;

    //default:
    //  GC_DO_ERROR( "Error, print(...) unknown type!\n");
    //  break;
    }
  }


  void
  GenericContainer::print_content_types(
    ostream_type &    stream,
    string_view const prefix,
    string_view const indent
  ) const {

    switch (m_data_type) {
    case GC_type::NOTYPE:
      stream << prefix << "Empty!\n";
      break;
    case GC_type::POINTER:
      stream << prefix << "(*void)\n";
      break;
    case GC_type::BOOL:
      stream << prefix << "bool\n";
      break;
    case GC_type::INTEGER:
      stream << prefix << "int\n";
      break;
    case GC_type::LONG:
      stream << prefix << "long int\n";
      break;
    case GC_type::REAL:
      stream << prefix << "double\n";
      break;
    case GC_type::COMPLEX:
      stream << prefix << "complex\n";
      break;
    case GC_type::STRING:
      stream << prefix << "string\n";
      break;
    case GC_type::VEC_POINTER:
      { vec_pointer_type const & v{*m_data.v_p};
        stream << "vector of pointer[" << v.size() << "]\n"; }
      break;
    case GC_type::VEC_BOOL:
      { vec_bool_type const & v{*m_data.v_b};
        stream << "vector of bool[" << v.size() << "]\n"; }
      break;
    case GC_type::VEC_INTEGER:
      { vec_int_type const & v{*m_data.v_i};
        stream << "vector of int[" << v.size() << "]\n"; }
      break;
    case GC_type::VEC_LONG:
      { vec_long_type const & v{*m_data.v_l};
        stream << "vector of long[" << v.size() << "]\n"; }
      break;
    case GC_type::VEC_REAL:
      { vec_real_type const & v{*m_data.v_r};
        stream << "vector of double[" << v.size() << "]\n"; }
      break;
    case GC_type::VEC_COMPLEX:
      { vec_complex_type const & v{*m_data.v_c};
        stream << "vector of complex[" << v.size() << "]\n"; }
      break;
    case GC_type::MAT_INTEGER:
      { mat_int_type const & m{*m_data.m_i};
        stream << "matrix of int[" << m.num_rows() << "," << m.num_cols() << "]\n"; }
      break;
    case GC_type::MAT_LONG:
      { mat_long_type const & m{*m_data.m_l};
        stream << "matrix of long[" << m.num_rows() << "," << m.num_cols() << "]\n"; }
      break;
    case GC_type::MAT_REAL:
      { mat_real_type const & m{*m_data.m_r};
        stream << "matrix of double[" << m.num_rows() << "," << m.num_cols() << "]\n"; }
      break;
    case GC_type::MAT_COMPLEX:
      { mat_complex_type const & m{*m_data.m_c};
        stream << "matrix of complex[" << m.num_rows() << "," << m.num_cols() << "]\n"; }
      break;
    case GC_type::VEC_STRING:
      { vec_string_type const & v{*m_data.v_s};
        stream << "vector of string[" << v.size() << "]\n"; }
      break;
    case GC_type::VECTOR:
      { vector_type const & v{*m_data.v};
        string prefix2{prefix}; prefix2 += indent;
        for ( vector_type::size_type i{0}; i < v.size(); ++i ) {
          GenericContainer const & vi{v[i]};
          if ( vi.simple_data() || vi.simple_vec_data()) {
            stream << prefix << i << ": ";
            vi.print_content_types(stream,"");
          } else {
            stream << prefix << i << ":\n";
            vi.print_content_types(stream,prefix2,indent);
          }
        }
      }
      break;
    case GC_type::MAP:
      { map_type const & m{*m_data.m};
        string prefix2{prefix}; prefix2 += indent;
        for ( const auto &[fst, snd] : m ) {
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          string_type matches[4];
          if ( int const pcreExecRet{pcre_for_GC.exec( fst, matches )}; pcreExecRet == 4 ) {
            string_type header = matches[3]; // header
            // found formatting
            if ( snd.simple_data() || snd.simple_vec_data() ) {
              stream << prefix << header << ": ";
              snd.print_content_types(stream,"");
            } else {
              if ( matches[1].length() > 1 ) stream << '\n'; // double ## --> add nel line
              stream << prefix << header;
              if ( !matches[2].empty() ) {
                stream << '\n' << prefix;
                char const fmt{ matches[2][0] }; // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt; // underline header
              } else {
                stream << ':';
              }
              stream << '\n';
              snd.print_content_types(stream,prefix2,indent);
            }
          } else {
            string_type header = pcreExecRet == 3 ? matches[3] : fst;
            if ( snd.simple_data() || snd.simple_vec_data() ) {
              stream << prefix << header << ": ";
              snd.print_content_types(stream,"");
            } else {
              stream << prefix << header << ":\n";
              snd.print_content_types(stream,prefix2,indent);
            }
          }
        }
      }
      break;

    //default:
    //  GC_DO_ERROR( "Error, print(...) unknown type!\n");
    //  break;
    }
  }

  /*
  //   _
  //  | |_ ___      __ _  ___
  //  | __/ _ \    / _` |/ __|
  //  | || (_) |  | (_| | (__
  //   \__\___/____\__, |\___|
  //         |_____|___/
  */

  void
  GenericContainer::to_gc( GenericContainer & gc ) const {
    switch (m_data_type) {
    case GC_type::NOTYPE:      gc.clear();                break;
    case GC_type::BOOL:        gc = m_data.b;    break;
    case GC_type::INTEGER:     gc = m_data.i;    break;
    case GC_type::LONG:        gc = m_data.l;    break;
    case GC_type::REAL:        gc = m_data.r;    break;
    case GC_type::COMPLEX:     gc = *m_data.c;   break;
    case GC_type::STRING:      gc = *m_data.s;   break;
    case GC_type::VEC_BOOL:    gc = *m_data.v_b; break;
    case GC_type::VEC_INTEGER: gc = *m_data.v_i; break;
    case GC_type::VEC_LONG:    gc = *m_data.v_l; break;
    case GC_type::VEC_REAL:    gc = *m_data.v_r; break;
    case GC_type::VEC_STRING:  gc = *m_data.v_s; break;

    case GC_type::VECTOR:
      gc.set_vector();
      { vector_type const & v{*m_data.v};
        vector_type       & vv = gc.set_vector( static_cast<unsigned>(v.size()) );
        for ( vector_type::size_type i{0}; i < v.size(); ++i )
          v[i].to_gc(vv[i]);
      }
      break;
    case GC_type::MAP:
      gc.set_map();
      { map_type const & m{*m_data.m};
        for ( const auto &[fst, snd] : m )
          snd.to_gc(gc[fst]);
      }
      break;
    case GC_type::MAT_INTEGER: gc = *m_data.m_i; break;
    case GC_type::MAT_LONG:    gc = *m_data.m_l; break;
    case GC_type::MAT_REAL:    gc = *m_data.m_r; break;
    case GC_type::VEC_COMPLEX: gc = *m_data.v_c; break;
    case GC_type::MAT_COMPLEX: gc = *m_data.m_c; break;
    case GC_type::POINTER:     gc = this->get_pointer<void *>(); break;
    case GC_type::VEC_POINTER:
      { vec_pointer_type const & v  = this->get_vec_pointer();
        vec_pointer_type       & vv = gc.set_vec_pointer(v.size());
        for ( vec_pointer_type::size_type i{0}; i < v.size(); ++i )
          vv[i] = v[i];
      }
      break;
    }
  }

  /*
  //    __
  //   / _|_ __ ___  _ __ ___       __ _  ___
  //  | |_| '__/ _ \| '_ ` _ \     / _` |/ __|
  //  |  _| | | (_) | | | | | |   | (_| | (__
  //  |_| |_|  \___/|_| |_| |_|____\__, |\___|
  //                         |_____|___/
  */

  void
  GenericContainer::from_gc( GenericContainer const & gc ) {
    switch ( gc.get_type()) {
    case GC_type::NOTYPE:
      this->clear();
      break;
    case GC_type::BOOL:
      this->set_bool(gc.m_data.b);
      break;
    case GC_type::INTEGER:
      this->set_int(gc.m_data.i);
      break;
    case GC_type::LONG:
      this->set_long(gc.m_data.l);
      break;
    case GC_type::REAL:
      this->set_real(gc.m_data.r);
      break;
    case GC_type::COMPLEX:
      this->set_complex(*gc.m_data.c);
      break;
    case GC_type::STRING:
      this->set_string(*gc.m_data.s);
      break;
    case GC_type::VEC_BOOL:
      this->set_vec_bool(*gc.m_data.v_b);
      break;
    case GC_type::VEC_INTEGER:
      this->set_vec_int(*gc.m_data.v_i);
      break;
    case GC_type::VEC_LONG:
      this->set_vec_long(*gc.m_data.v_l);
      break;
    case GC_type::VEC_REAL:
      this->set_vec_real(*gc.m_data.v_r);
      break;
    case GC_type::VEC_STRING:
      this->set_vec_string(*gc.m_data.v_s);
      break;

    case GC_type::VECTOR:
      this->set_vector();
      { vector_type const & v{*gc.m_data.v};
        vector_type       & vv = this->set_vector( static_cast<unsigned>(v.size()) );
        for ( vector_type::size_type i{0}; i < v.size(); ++i )
          vv[i].from_gc(v[i]);
      }
      break;
    case GC_type::MAP:
      this->set_map();
      { map_type const & m{*gc.m_data.m};
        for ( const auto &[fst, snd] : m )
          (*this)[fst].from_gc(snd);
      }
      break;
    case GC_type::MAT_INTEGER:
      this->set_mat_int(*gc.m_data.m_i);
      break;
    case GC_type::MAT_LONG:
      this->set_mat_long(*gc.m_data.m_l);
      break;
    case GC_type::MAT_REAL:
      this->set_mat_real(*gc.m_data.m_r);
      break;
    case GC_type::VEC_COMPLEX:
      this->set_vec_complex(*gc.m_data.v_c);
      break;
    case GC_type::MAT_COMPLEX:
      this->set_mat_complex(*gc.m_data.m_c);
      break;
    case GC_type::POINTER:
      this->set_pointer( gc.get_pointer<void *>() );
      break;
    case GC_type::VEC_POINTER:
      { vec_pointer_type const & v  = gc.get_vec_pointer();
        vec_pointer_type       & vv = this->set_vec_pointer(v.size());
        for ( vec_pointer_type::size_type i{0}; i < v.size(); ++i )
          vv[i] = v[i];
      }
      break;
    }
  }

  //
  // -------------------------------------------------------------------
  //

  void
  GenericContainer::merge( GenericContainer const & gc, string_view const where ) {
    if ( gc.get_type() == GC_type::NOTYPE ) return;
    GC_ASSERT(
      gc.get_type() == GC_type::MAP,
      where
        << " in merge data expected to be of type: " << to_string(GC_type::MAP)
        << " but data stored is of type: " << gc.get_type_name()
    )
    if ( m_data_type == GC_type::NOTYPE ) this->set_map();
    ck( where, GC_type::MAP );
    { map_type const & m{gc.get_map()};
      for ( const auto &[fst, snd] : m )
        (*this)[fst].from_gc(snd);
    }
  }
  
  bool
  GenericContainer::from_file( string_view file_name ) {
    std::ifstream file(file_name.data());
    if ( !file.is_open() ) return false;

    // Utility to check file extension
    auto ends_with = []( string_view str, string_view suffix ) -> bool {
      return str.size() >= suffix.size() &&
             str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    bool success{ false };
    if ( ends_with( file_name, ".yaml" ) ||
         ends_with( file_name, ".yml"  ) ) {
      success = this->from_yaml(file);
    } else if ( ends_with(file_name, ".json" ) ) {
      success = this->from_json(file);
    } else if (ends_with(file_name, ".toml")) {
      success = this->from_toml(file);
    }
    file.close();
    return success;
  }

  void
  GenericContainer::exception( string_view const where ) {
    throw std::runtime_error(where.data());
  }

  // instantate classes
  template class mat_type<int_type>;
  template class mat_type<long_type>;
  template class mat_type<real_type>;
  template class mat_type<complex_type>;

}

//
// eof: GenericContainer.cc
//
