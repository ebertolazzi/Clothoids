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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: GenericContainer.cc
//

#include "GenericContainer.hh"
#include <iomanip>
#include <cmath>
#include <ctgmath>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#if __cplusplus >= 201103L &&                             \
    (!defined(__GLIBCXX__) || (__cplusplus >= 201402L) || \
        (defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) || \
         defined(_GLIBCXX_REGEX_STATE_LIMIT)           || \
             (defined(_GLIBCXX_RELEASE)                && \
             _GLIBCXX_RELEASE > 4)))
#define HAVE_WORKING_REGEX 1
#endif

#ifdef HAVE_WORKING_REGEX
  #include <regex>
#endif


#include <fstream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define CHECK_RESIZE(pV,I) if ( pV->size() <= (I) ) pV->resize((I)+1)

using std::fpclassify;
using GC_namespace::real_type;

static
inline
bool isZero( real_type x )
{ return FP_ZERO == fpclassify(x); }

static
inline
bool isInteger32( real_type x )
{ return isZero( x-static_cast<int32_t>(floor(x)) ); }

static
inline
bool isUnsigned32( real_type x )
{ return isInteger32(x) && x >= 0; }

static
inline
bool isInteger64( real_type x )
{ return isZero( x-static_cast<int64_t>(floor(x)) ); }

static
inline
bool isUnsigned64( real_type x )
{ return isInteger64(x) && x >= 0; }

#endif

namespace GC_namespace {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <typename TYPE>
  ostream_type &
  operator << ( ostream_type & s, std::vector<TYPE> const & v ) {
    if ( v.size() > 0 ) {
      s << "[ " << v[0];
      for ( unsigned i = 1; i < v.size(); ++i ) s << " " << v[i];
      s << " ]";
    } else {
      s << "[]";
    }
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, vec_bool_type const & v ) {
    if ( v.size() > 0 ) {
      s << (v[0]?"[ true":"[ false");
      for ( unsigned i = 1; i < v.size(); ++i )
        s << (v[i]?" true":" false");
      s << " ]";
    } else {
      s << "[]";
    }
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, vec_complex_type const & v ) {
    if ( v.size() > 0 ) {
      s << "[ (" << v[0].real() << ", " << v[0].imag() << " )";
      for ( unsigned i = 1; i < v.size(); ++i )
        s << " (" << v[i].real() << ", " << v[i].imag() << " )";
      s << " ]";
    } else {
      s << "[]";
    }
    return s;
  }

  template ostream_type & operator << ( ostream_type & s, vec_int_type const & v );
  template ostream_type & operator << ( ostream_type & s, vec_long_type const & v );
  template ostream_type & operator << ( ostream_type & s, vec_real_type const & v );
  #endif

  template <typename TYPE>
  TYPE const &
  mat_type<TYPE>::operator () ( unsigned i, unsigned j ) const {
    GC_ASSERT(
      i < m_numRows && j < m_numCols,
      "mat_type::operator() (" << i << ", " << j <<
      ") index out of range: [0," << m_numRows <<
      ") x [0," << m_numCols << ")\n"
    )
    return (*this)[std::size_t(i+j*m_numRows)];
  }

  template <typename TYPE>
  TYPE &
  mat_type<TYPE>::operator () ( unsigned i, unsigned j ) {
    GC_ASSERT(
      i < m_numRows && j < m_numCols,
      "mat_type::operator() (" << i << ", " << j <<
      ") index out of range: [0," << m_numRows <<
      ") x [0," << m_numCols << ")\n"
    )
    return (*this)[std::size_t(i+j*m_numRows)];
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::getColumn( unsigned nc, std::vector<TYPE> & C ) const {
    GC_ASSERT(
      nc < m_numCols,
      "mat_type::getColumn(" << nc <<
      ",C) column index out of range max = " << m_numCols-1
    )
    C.clear();
    C.reserve(m_numRows);
    for ( unsigned i = 0; i < m_numRows; ++i )
      C.push_back( (*this)(i,nc) );
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::getColumn( unsigned nc, TYPE * C ) const {
    GC_ASSERT(
      nc < m_numCols,
      "mat_type::getColumn(" << nc <<
      ",C) column index out of range max = " << m_numCols-1
    )
    for ( unsigned i = 0; i < m_numRows; ++i )
      *C++ = (*this)(i,nc);
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::getRow( unsigned nr, std::vector<TYPE> & R ) const {
    GC_ASSERT(
      nr < m_numRows,
      "mat_type::getRow(" << nr <<
      ",C) row index out of range max = " << m_numRows-1
    );
    R.clear();
    R.reserve(m_numCols);
    for ( unsigned j = 0; j < m_numCols; ++j )
      R.push_back( (*this)(nr,j) );
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::getRow( unsigned nr, TYPE * R ) const {
    GC_ASSERT(
      nr < m_numRows,
      "mat_type::getRow(" << nr <<
      ",C) row index out of range max = " << m_numRows-1
    )
    for ( unsigned j = 0; j < m_numCols; ++j )
      *R++ = (*this)(nr,j);
  }

  template <typename TYPE>
  void
  mat_type<TYPE>::info( ostream_type & stream ) const {
    stream
      << "Matrix of floating point number of size "
      << m_numRows << " x " << m_numCols
      << '\n';
  }

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename TYPE>
  ostream_type &
  operator << ( ostream_type & s, mat_type<TYPE> const & m ) {
    for ( unsigned i = 0; i < m.numRows(); ++i ) {
      s << std::setw(8) << m(i,0);
      for ( unsigned j = 1; j < m.numCols(); ++j )
        s << " " << std::setw(8) << m(i,j);
      s << '\n';
    }
    return s;
  }

  template <>
  ostream_type &
  operator << ( ostream_type & s, mat_complex_type const & m ) {
    for ( unsigned i = 0; i < m.numRows(); ++i ) {
      for ( unsigned j = 0; j < m.numCols(); ++j )
        s << " (" << std::setw(8) << m(i,j).real() << ", "
                  << std::setw(8) << m(i,j).imag() << " )";
      s << '\n';
    }
    return s;
  }

  template ostream_type & operator << ( ostream_type & s, mat_type<int_type> const & m );
  template ostream_type & operator << ( ostream_type & s, mat_type<long_type> const & m );
  template ostream_type & operator << ( ostream_type & s, mat_type<real_type> const & m );
  #endif

  #ifdef HAVE_WORKING_REGEX
  class Pcre_for_GC {

  private:

    std::regex  reCompiled;
    std::smatch reMatches;

  public:

    Pcre_for_GC()
    : reCompiled("^\\s*\\d+\\s*(##?)(-|=|~|_|)\\s*(.*)$")
    { }

    ~Pcre_for_GC()
    { }

    int
    exec( std::string const & str, std::string matches[4] ) {
      if ( std::regex_match( str, reMatches, reCompiled ) ) {
        for ( std::size_t i = 0; i < reMatches.size(); ++i )
          matches[i] = reMatches[i].str();
        return int(reMatches.size());
      } else {
        return 0;
      }
    }
  };

  static Pcre_for_GC pcre_for_GC;

  #endif

  static char const *typeName[] = {

    "NOTYPE",

    "pointer",
    "bool_type",
    "int_type",
    "long_type",
    "real_type",
    "complex_type",
    "string_type",

    "vec_pointer_type",
    "vec_bool_type",
    "vec_int_type",
    "vec_long_type",
    "vec_real_type",
    "vec_complex_type",
    "vec_string_type",

    "mat_int_type",
    "mat_long_type",
    "mat_real_type",
    "mat_complex_type",

    "vector_type",
    "map_type"
  };

  // costruttore
  GenericContainer::GenericContainer()
  : m_data_type(GC_NOTYPE)
  {
  }

  #ifdef GENERIC_CONTAINER_ON_WINDOWS
  bool
  GenericContainer::simple_data() const {
    return m_data_type <= GC_STRING;
  }
  bool
  GenericContainer::simple_vec_data() const {
    return m_data_type <= GC_VEC_STRING;
  }
  #endif

  GenericContainer &
  GenericContainer::operator = ( vec_bool_type const & a ) {
    set_vec_bool( unsigned(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_b->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_int_type const & a ) {
    set_vec_int( unsigned(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_i->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_long_type const & a ) {
    set_vec_long( unsigned(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_l->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_real_type const & a ) {
    set_vec_real( unsigned(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_r->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_complex_type const & a ) {
    set_vec_complex( unsigned(a.size()) );
    std::copy( a.begin(), a.end(), m_data.v_c->begin() );
    return *this;
  }

  GenericContainer &
  GenericContainer::operator = ( vec_string_type const & a ) {
    set_vec_string( unsigned(a.size()) );
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

  // distruttore
  void
  GenericContainer::clear() {
    switch (m_data_type) {
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
      // removed annoying warning. To be re-thinked...
      //GC_WARNING( _data.p == nullptr, "find a pointer not deallocated!" );
      break;
    case GC_STRING:      delete m_data.s; break;
    case GC_COMPLEX:     delete m_data.c; break;

    case GC_VEC_POINTER: delete m_data.v_p; break;
    case GC_VEC_BOOL:    delete m_data.v_b; break;
    case GC_VEC_INTEGER: delete m_data.v_i; break;
    case GC_VEC_LONG:    delete m_data.v_l; break;
    case GC_VEC_REAL:    delete m_data.v_r; break;
    case GC_VEC_COMPLEX: delete m_data.v_c; break;

    case GC_MAT_INTEGER: delete m_data.m_i; break;
    case GC_MAT_LONG:    delete m_data.m_l; break;
    case GC_MAT_REAL:    delete m_data.m_r; break;
    case GC_MAT_COMPLEX: delete m_data.m_c; break;
    case GC_VEC_STRING:  delete m_data.v_s; break;

    case GC_VECTOR:
      {
        vector_type::iterator it = m_data.v->begin();
        for ( ; it != m_data.v->end(); ++it ) it->clear();
        delete m_data.v;
      }
      break;
    case GC_MAP:
      {
        map_type::iterator it = m_data.m->begin();
        for ( ; it != m_data.m->end(); ++it ) it->second.clear();
        delete m_data.m;
      }
      break;
    }
    m_data_type = GC_NOTYPE;
  }

  // distruttore
  void
  GenericContainer::erase( char const * name ) {
    GC_ASSERT(
      GC_MAP == m_data_type,
      "erase('" << name << "') bad data type\nexpect: " <<
      typeName[GC_POINTER] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    m_data.m->erase(name);
  }

  char const *
  GenericContainer::get_type_name() const {
    return typeName[m_data_type];
  }

  unsigned
  GenericContainer::get_num_elements() const {
    switch (m_data_type) {
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:      return 1;

    case GC_VEC_POINTER: return unsigned(m_data.v_p->size());
    case GC_VEC_BOOL:    return unsigned(m_data.v_b->size());
    case GC_VEC_INTEGER: return unsigned(m_data.v_i->size());
    case GC_VEC_LONG:    return unsigned(m_data.v_l->size());
    case GC_VEC_REAL:    return unsigned(m_data.v_r->size());
    case GC_VEC_COMPLEX: return unsigned(m_data.v_c->size());
    case GC_VEC_STRING:  return unsigned(m_data.v_s->size());

    case GC_MAT_INTEGER: return unsigned(m_data.m_i->size());
    case GC_MAT_LONG:    return unsigned(m_data.m_l->size());
    case GC_MAT_REAL:    return unsigned(m_data.m_r->size());
    case GC_MAT_COMPLEX: return unsigned(m_data.m_c->size());

    case GC_VECTOR:      return unsigned(m_data.v->size());
    case GC_MAP:         return unsigned(m_data.m->size());
    case GC_NOTYPE:      return 0;
    }
    return 0;
  }

  unsigned
  GenericContainer::get_numRows() const {
    switch (m_data_type) {
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:      return 1;
    case GC_MAT_INTEGER: return unsigned(m_data.m_i->numRows());
    case GC_MAT_LONG:    return unsigned(m_data.m_l->numRows());
    case GC_MAT_REAL:    return unsigned(m_data.m_r->numRows());
    case GC_MAT_COMPLEX: return unsigned(m_data.m_c->numRows());
    case GC_MAP:         return 1;
    case GC_NOTYPE:      return 0;
    }
    return 0;
  }

  unsigned
  GenericContainer::get_numCols() const {
    switch (m_data_type) {
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:      return 1;
    case GC_VEC_POINTER: return unsigned(m_data.v_p->size());
    case GC_VEC_BOOL:    return unsigned(m_data.v_b->size());
    case GC_VEC_INTEGER: return unsigned(m_data.v_i->size());
    case GC_VEC_LONG:    return unsigned(m_data.v_l->size());
    case GC_VEC_REAL:    return unsigned(m_data.v_r->size());
    case GC_VEC_COMPLEX: return unsigned(m_data.v_c->size());
    case GC_VEC_STRING:  return unsigned(m_data.v_s->size());

    case GC_MAT_INTEGER: return unsigned(m_data.m_i->numCols());
    case GC_MAT_LONG:    return unsigned(m_data.m_l->numCols());
    case GC_MAT_REAL:    return unsigned(m_data.m_r->numCols());
    case GC_MAT_COMPLEX: return unsigned(m_data.m_c->numCols());

    case GC_VECTOR:      return unsigned(m_data.v->size());
    case GC_MAP:         return unsigned(m_data.m->size());
    case GC_NOTYPE:      return 0;
    }
    return 0;
  }

  void
  GenericContainer::load( GenericContainer const & gc ) {
    this->clear();
    switch (gc.m_data_type) {
    case GC_NOTYPE:      break;
    case GC_POINTER:     this->set_pointer(gc.m_data.p);  break;
    case GC_BOOL:        this->set_bool(gc.m_data.b);     break;
    case GC_INTEGER:     this->set_int(gc.m_data.i);      break;
    case GC_LONG:        this->set_long(gc.m_data.l);     break;
    case GC_REAL:        this->set_real(gc.m_data.r);     break;
    case GC_COMPLEX:     this->set_complex(*gc.m_data.c); break;
    case GC_STRING:      this->set_string(*gc.m_data.s);  break;

    case GC_VEC_POINTER: this->set_vec_pointer(*gc.m_data.v_p); break;
    case GC_VEC_BOOL:    this->set_vec_bool(*gc.m_data.v_b);    break;
    case GC_VEC_INTEGER: this->set_vec_int(*gc.m_data.v_i);     break;
    case GC_VEC_LONG:    this->set_vec_long(*gc.m_data.v_l);    break;
    case GC_VEC_REAL:    this->set_vec_real(*gc.m_data.v_r);    break;
    case GC_VEC_COMPLEX: this->set_vec_complex(*gc.m_data.v_c); break;
    case GC_VEC_STRING:  this->set_vec_string(*gc.m_data.v_s);  break;

    case GC_MAT_INTEGER: this->set_mat_int(*gc.m_data.m_i);     break;
    case GC_MAT_LONG:    this->set_mat_long(*gc.m_data.m_l);    break;
    case GC_MAT_REAL:    this->set_mat_real(*gc.m_data.m_r);    break;
    case GC_MAT_COMPLEX: this->set_mat_complex(*gc.m_data.m_c); break;

    case GC_VECTOR:
      { unsigned N = unsigned(gc.m_data.v->size());
        allocate_vector( N );
        std::copy( gc.m_data.v->begin(),
                   gc.m_data.v->end(),
                   this->m_data.v->begin() );
      }
      break;
    case GC_MAP:
      { allocate_map();
        // this->_data.m->insert( gc._data.m->begin(), gc._data.m->end() ); !!!! DO NOT WORK ON CLANG
        for ( map_type::iterator it = gc.m_data.m->begin();
              it != gc.m_data.m->end(); ++it )
          (*this->m_data.m)[it->first] = it->second;
      }
      break;
    }
  }

  int
  GenericContainer::ck( TypeAllowed tp ) const {
    if ( tp == m_data_type ) return 0; // ok
    if ( tp == GC_NOTYPE   ) return 1; //
    return 2;
  }

  void
  GenericContainer::ck( char const * msg_in, TypeAllowed tp ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    GC_ASSERT(
      tp == m_data_type,
      msg
        << "bad data type, expect: " << typeName[tp]
        << " but data stored is of type: " << typeName[m_data_type]
    )
  }

  void
  GenericContainer::ck_or_set( char const * msg, TypeAllowed tp) {
    if ( m_data_type == GC_NOTYPE ) m_data_type = tp;
    else                            ck(msg,tp);
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
    if ( m_data_type != GC_STRING ) {
      clear();
      m_data_type = GC_STRING;
      m_data.s    = new string_type;
    }
  }

  void
  GenericContainer::allocate_complex() {
    if ( m_data_type != GC_COMPLEX ) {
      clear();
      m_data_type = GC_COMPLEX;
      m_data.c    = new complex_type;
    }
  }

  void
  GenericContainer::allocate_vec_pointer( unsigned sz ) {
    if ( m_data_type != GC_VEC_POINTER ) {
      clear();
      m_data_type = GC_VEC_POINTER;
      m_data.v_p  = new vec_pointer_type();
    }
    if ( sz > 0 ) m_data.v_p->resize( sz );
  }

  GenericContainer &
  GenericContainer::free_pointer() {
    GC_ASSERT(
      GC_POINTER == m_data_type || GC_NOTYPE == m_data_type,
      "free_pointer() bad data type\nexpect: " << typeName[GC_POINTER] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    m_data.p = nullptr;
    m_data_type = GC_NOTYPE;
    return *this;
  }

  void
  GenericContainer::allocate_vec_bool( unsigned sz ) {
    if ( m_data_type != GC_VEC_BOOL ) {
      clear();
      m_data_type = GC_VEC_BOOL;
      m_data.v_b  = new vec_bool_type();
    }
    if ( sz > 0 ) m_data.v_b->resize( sz );
  }

  void
  GenericContainer::allocate_vec_int( unsigned sz ) {
    if ( m_data_type != GC_VEC_INTEGER ) {
      clear();
      m_data_type = GC_VEC_INTEGER;
      m_data.v_i  = new vec_int_type();
    }
    if ( sz > 0 ) m_data.v_i->resize( sz );
  }

  void
  GenericContainer::allocate_vec_long( unsigned sz ) {
    if ( m_data_type != GC_VEC_LONG ) {
      clear();
      m_data_type = GC_VEC_LONG;
      m_data.v_l  = new vec_long_type();
    }
    if ( sz > 0 ) m_data.v_l->resize( sz );
  }

  void
  GenericContainer::allocate_vec_real( unsigned sz ) {
    if ( m_data_type != GC_VEC_REAL ) {
      clear();
      m_data_type = GC_VEC_REAL;
      m_data.v_r  = new vec_real_type();
    }
    if ( sz > 0 ) m_data.v_r->resize( sz );
  }

  void
  GenericContainer::allocate_vec_complex( unsigned sz ) {
    if ( m_data_type != GC_VEC_COMPLEX ) {
      clear();
      m_data_type = GC_VEC_COMPLEX;
      m_data.v_c  = new vec_complex_type();
    }
    if ( sz > 0 ) m_data.v_c->resize( sz );
  }

  void
  GenericContainer::allocate_mat_int( unsigned nr, unsigned nc ) {
    if ( m_data_type != GC_MAT_INTEGER ) {
      clear();
      m_data_type = GC_MAT_INTEGER;
      m_data.m_i  = new mat_int_type( nr, nc );
    } else {
      m_data.m_i->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_long( unsigned nr, unsigned nc ) {
    if ( m_data_type != GC_MAT_LONG ) {
      clear();
      m_data_type = GC_MAT_LONG;
      m_data.m_l  = new mat_long_type( nr, nc );
    } else {
      m_data.m_l->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_real( unsigned nr, unsigned nc ) {
    if ( m_data_type != GC_MAT_REAL ) {
      clear();
      m_data_type = GC_MAT_REAL;
      m_data.m_r  = new mat_real_type( nr, nc );
    } else {
      m_data.m_r->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_mat_complex( unsigned nr, unsigned nc ) {
    if ( m_data_type != GC_MAT_COMPLEX ) {
      clear();
      m_data_type = GC_MAT_COMPLEX;
      m_data.m_c  = new mat_complex_type( nr, nc );
    } else {
      m_data.m_c->resize( nr, nc );
    }
  }

  void
  GenericContainer::allocate_vec_string( unsigned sz ) {
    if ( m_data_type != GC_VEC_STRING ) {
      clear();
      m_data_type = GC_VEC_STRING;
      m_data.v_s  = new vec_string_type();
    }
    if ( sz > 0 ) m_data.v_s->resize( sz );
  }

  void
  GenericContainer::allocate_vector( unsigned sz ) {
    if ( m_data_type != GC_VECTOR ) {
      clear();
      m_data_type = GC_VECTOR;
      m_data.v    = new vector_type();
    }
    if ( sz > 0 ) m_data.v->resize( sz );
  }

  void
  GenericContainer::allocate_map() {
    if ( m_data_type != GC_MAP ) {
      clear();
      m_data_type = GC_MAP;
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
  GenericContainer::set_pointer( pointer_type value ) {
    clear();
    m_data_type = GC_POINTER;
    return (m_data.p = value);
  }

  bool_type &
  GenericContainer::set_bool( bool_type value ) {
    clear();
    m_data_type = GC_BOOL;
    return (m_data.b = value);
  }

  int_type &
  GenericContainer::set_int( int_type value ) {
    clear();
    m_data_type = GC_INTEGER;
    return (m_data.i = value);
  }

  long_type &
  GenericContainer::set_long( long_type value ) {
    clear();
    m_data_type = GC_LONG;
    return (m_data.l = value);
  }

  real_type &
  GenericContainer::set_real( real_type value ) {
    clear();
    m_data_type = GC_REAL;
    return (m_data.r = value);
  }

  complex_type &
  GenericContainer::set_complex( complex_type & value ) {
    clear();
    m_data_type = GC_COMPLEX;
    m_data.c    = new complex_type;
    return (*m_data.c=value);
  }

  complex_type &
  GenericContainer::set_complex( real_type re, real_type im ) {
    clear();
    m_data_type = GC_COMPLEX;
    m_data.c    = new complex_type(re,im);
    return *m_data.c;
  }

  string_type &
  GenericContainer::set_string( string_type const & value ) {
    allocate_string();
    return (*m_data.s = value);
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( unsigned sz ) {
    allocate_vec_pointer( sz );
    return *m_data.v_p;
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( vec_pointer_type const & v ) {
    allocate_vec_pointer( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_p->begin() );
    return *m_data.v_p;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( unsigned sz ) {
    allocate_vec_bool( sz ); return *m_data.v_b;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( vec_bool_type const & v ) {
    allocate_vec_bool( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_b->begin() );
    return *m_data.v_b;
  }

  vec_int_type &
  GenericContainer::set_vec_int( unsigned sz ) {
    allocate_vec_int( sz );
    return *m_data.v_i;
  }

  vec_int_type &
  GenericContainer::set_vec_int( vec_int_type const & v ) {
    allocate_vec_int( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_i->begin() );
    return *m_data.v_i;
  }

  vec_long_type &
  GenericContainer::set_vec_long( unsigned sz ) {
    allocate_vec_long( sz );
    return *m_data.v_l;
  }

  vec_long_type &
  GenericContainer::set_vec_long( vec_long_type const & v ) {
    allocate_vec_long( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_l->begin() );
    return *m_data.v_l;
  }

  vec_real_type &
  GenericContainer::set_vec_real( unsigned sz ) {
    allocate_vec_real( sz );
    return *m_data.v_r;
  }

  vec_real_type &
  GenericContainer::set_vec_real( vec_real_type const & v ) {
    allocate_vec_real( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_r->begin() );
    return *m_data.v_r;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( unsigned sz ) {
    allocate_vec_complex( sz );
    return *m_data.v_c;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( vec_complex_type const & v ) {
    allocate_vec_complex( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_c->begin() );
    return *m_data.v_c;
  }

  mat_int_type &
  GenericContainer::set_mat_int( unsigned nr, unsigned nc ) {
    allocate_mat_int( nr, nc );
    return *m_data.m_i;
  }

  mat_int_type &
  GenericContainer::set_mat_int( mat_int_type const & m ) {
    allocate_mat_int( m.numRows(), m.numCols() );
    std::copy( m.begin(), m.end(), m_data.m_i->begin() );
    return *m_data.m_i;
  }

  mat_long_type &
  GenericContainer::set_mat_long( unsigned nr, unsigned nc ) {
    allocate_mat_long( nr, nc );
    return *m_data.m_l;
  }

  mat_long_type &
  GenericContainer::set_mat_long( mat_long_type const & m ) {
    allocate_mat_long( m.numRows(), m.numCols() );
    std::copy( m.begin(), m.end(), m_data.m_l->begin() );
    return *m_data.m_l;
  }

  mat_real_type &
  GenericContainer::set_mat_real( unsigned nr, unsigned nc ) {
    allocate_mat_real( nr, nc );
    return *m_data.m_r;
  }

  mat_real_type &
  GenericContainer::set_mat_real( mat_real_type const & m ) {
    allocate_mat_real( m.numRows(), m.numCols() );
    std::copy( m.begin(), m.end(), m_data.m_r->begin() );
    return *m_data.m_r;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( unsigned nr, unsigned nc ) {
    allocate_mat_complex( nr, nc );
    return *m_data.m_c;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( mat_complex_type const & m ) {
    allocate_mat_complex( m.numRows(), m.numCols() );
    std::copy( m.begin(), m.end(), m_data.m_c->begin() );
    return *m_data.m_c;
  }

  vec_string_type &
  GenericContainer::set_vec_string( unsigned sz ) {
    allocate_vec_string( sz );
    return *m_data.v_s;
  }

  vec_string_type &
  GenericContainer::set_vec_string( vec_string_type const & v ) {
    allocate_vec_string( unsigned(v.size()) );
    std::copy( v.begin(), v.end(), m_data.v_s->begin() );
    return *m_data.v_s;
  }

  vector_type &
  GenericContainer::set_vector( unsigned sz ) {
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
  GenericContainer::push_bool( bool val ) {
    if ( m_data_type == GC_VEC_BOOL ) {
      m_data.v_b->push_back( val );
    } else if ( m_data_type == GC_VEC_INTEGER ) {
      m_data.v_i->push_back( val ? 1 : 0 );
    } else if ( m_data_type == GC_VEC_LONG ) {
      m_data.v_l->push_back( val ? 1 : 0 );
    } else if ( m_data_type == GC_VEC_REAL ) {
      m_data.v_r->push_back( val ? 1 : 0 );
    } else if ( m_data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val ? 1 : 0, 0 );
      m_data.v_c->push_back( tmp );
    } else if ( m_data_type == GC_VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_bool( val );
    } else {
      GC_DO_ERROR( "push_bool, bad data stored: " << get_type_name() )
    }
  }

  void
  GenericContainer::push_int( int_type val ) {
    if ( m_data_type == GC_VEC_INTEGER ) {
      m_data.v_i->push_back( val );
    } else if ( m_data_type == GC_VEC_LONG ) {
      m_data.v_l->push_back( long_type(val) );
    } else if ( m_data_type == GC_VEC_REAL ) {
      m_data.v_r->push_back( val );
    } else if ( m_data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val, 0 );
      m_data.v_c->push_back( tmp );
    } else if ( m_data_type == GC_VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_int( val );
    } else {
      if ( m_data_type != GC_VEC_INTEGER ) promote_to_vec_int();
      m_data.v_i->push_back( val );
    }
  }

  void
  GenericContainer::push_long( long_type val ) {
    if ( m_data_type == GC_VEC_LONG ) {
      m_data.v_l->push_back( val );
    } else if ( m_data_type == GC_VEC_REAL ) {
      m_data.v_r->push_back( real_type(val) );
    } else if ( m_data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( real_type(val), 0 );
      m_data.v_c->push_back( tmp );
    } else if ( m_data_type == GC_VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_long( val );
    } else {
      if ( m_data_type != GC_VEC_LONG ) promote_to_vec_long();
      m_data.v_l->push_back( val );
    }
  }

  void
  GenericContainer::push_real( real_type val ) {
    if ( m_data_type == GC_VEC_REAL ) {
      m_data.v_r->push_back( val );
    } else if ( m_data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val, 0 );
      m_data.v_c->push_back( tmp );
    } else if ( m_data_type == GC_VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_real( val );
    } else {
      if ( m_data_type != GC_VEC_REAL ) promote_to_vec_real();
      m_data.v_r->push_back( val );
    }
  }

  void
  GenericContainer::push_complex( complex_type & val ) {
    if ( m_data_type == GC_VECTOR ) {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_complex( val );
    } else {
      if ( m_data_type != GC_VEC_COMPLEX ) promote_to_vec_complex();
      m_data.v_c->push_back( val );
    }
  }

  void
  GenericContainer::push_complex( real_type re, real_type im ) {
    complex_type tmp( re, im );
    push_complex( tmp );
  }

  void
  GenericContainer::push_string( string_type const & val ) {
    if ( m_data_type != GC_VEC_STRING ) promote_to_vector();
    if ( m_data_type == GC_VEC_STRING ) {
      m_data.v_s->push_back( val );
    } else {
      m_data.v->resize(m_data.v->size()+1);
      m_data.v->back().set_string( val );
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
  GenericContainer::get_pvoid( char const * msg ) const {
    ck(msg,GC_POINTER);
    return m_data.p;
  }

  void **
  GenericContainer::get_ppvoid( char const * msg ) const {
    ck(msg,GC_POINTER);
    return const_cast<void **>(&m_data.p);
  }

  int_type const *
  GenericContainer::get_int_pointer() const {
    switch (m_data_type) {
    case GC_INTEGER:
      return &m_data.i;
    case GC_VEC_INTEGER:
      return m_data.v_i->data();
    case GC_MAT_INTEGER:
      return m_data.m_i->data();
    case GC_NOTYPE:
    case GC_BOOL:
    //case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    //case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    //case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `int_type const*'"
      )
    }
    return nullptr;
  }

  int_type *
  GenericContainer::get_int_pointer() {
    switch (m_data_type) {
    case GC_INTEGER:
      return &m_data.i;
    case GC_VEC_INTEGER:
      return m_data.v_i->data();
    case GC_MAT_INTEGER:
      return m_data.m_i->data();
    case GC_NOTYPE:
    case GC_BOOL:
    //case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    //case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    //case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `int_type*'"
      )
    }
    return nullptr;
  }

  long_type const *
  GenericContainer::get_long_pointer() const {
    switch (m_data_type) {
    case GC_LONG:
      return &m_data.l;
    case GC_VEC_LONG:
      return m_data.v_l->data();
    case GC_MAT_LONG:
      return m_data.m_l->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    //case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    //case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    //case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_long_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `long_type const*'"
      )
    }
    return nullptr;
  }

  long_type *
  GenericContainer::get_long_pointer() {
    switch (m_data_type) {
    case GC_LONG:
      return &m_data.l;
    case GC_VEC_LONG:
      return m_data.v_l->data();
    case GC_MAT_LONG:
      return m_data.m_l->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    //case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    //case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    //case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_long_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `long_type*'"
      )
    }
    return nullptr;
  }

  real_type const *
  GenericContainer::get_real_pointer() const {
    switch (m_data_type) {
    case GC_REAL:
      return &m_data.r;
    case GC_VEC_REAL:
      return m_data.v_r->data();
    case GC_MAT_REAL:
      return m_data.m_r->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    //case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    //case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    //case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_real_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `real_type cont *'"
      )
    }
    return nullptr;
  }

  real_type *
  GenericContainer::get_real_pointer() {
    switch (m_data_type) {
    case GC_REAL:
      return &m_data.r;
    case GC_VEC_REAL:
      return m_data.v_r->data();
    case GC_MAT_REAL:
      return m_data.m_r->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    //case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    //case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    //case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_real_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `real_type*'"
      )
    }
    return nullptr;
  }

  complex_type const *
  GenericContainer::get_complex_pointer() const {
    switch (m_data_type) {
    case GC_COMPLEX:
      return m_data.c;
    case GC_VEC_COMPLEX:
      return m_data.v_c->data();
    case GC_MAT_COMPLEX:
      return m_data.m_c->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    //case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    //case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    //case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `complex_type const*'"
      )
    }
    return nullptr;
  }

  complex_type *
  GenericContainer::get_complex_pointer() {
    switch (m_data_type) {
    case GC_COMPLEX:
      return m_data.c;
    case GC_VEC_COMPLEX:
      return m_data.v_c->data();
    case GC_MAT_COMPLEX:
      return m_data.m_c->data();
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    //case GC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    //case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    //case GC_MAT_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        "get_int_pointer, bad data type: `" << typeName[m_data_type] <<
        "' cannot be referred as `complex_type const*'"
      )
    }
    return nullptr;
  }

  template <>
  void
  GenericContainer::get_value( unsigned & v, char const * msg ) const {
    switch (m_data_type) {
    case GC_BOOL:
      v = m_data.b?1:0;
      break;
    case GC_INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        msg << "value '" << m_data.i << "' cannot be converted in `unsigned'"
      )
      v = unsigned(m_data.i);
      break;
    case GC_LONG:
      GC_ASSERT(
        m_data.l >= 0,
        msg << "value '" << m_data.l << "' cannot be converted in `unsigned'"
      )
      v = unsigned(m_data.l);
      break;
    case GC_REAL:
      GC_ASSERT(
        m_data.r >= 0 && isUnsigned32(m_data.r),
        msg << "value '" << m_data.r << "' cannot be converted in `unsigned'"
      )
      v = unsigned(m_data.r);
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()) && isUnsigned32(m_data.c->real()),
        msg << "value = (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") cannot be converted in `unsigned'"
      )
      v = unsigned(m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `unsigned'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( int & v, char const * msg ) const {
    switch (m_data_type) {
    case GC_BOOL:
      v = m_data.b?1:0;
      break;
    case GC_INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        msg << "value '" << m_data.i << "' cannot be converted in `int'"
      )
      v = m_data.i;
      break;
    case GC_LONG:
      v = int(m_data.l);
      break;
    case GC_REAL:
      GC_ASSERT(
        isInteger32(m_data.r),
        msg << "value '" << m_data.r << "' cannot be converted in `int'"
      )
      v = int(m_data.r);
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()) && isInteger32(m_data.c->real()),
        msg << "value = (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") cannot be converted in `int'"
      )
      v = int(m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `unsigned'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( unsigned long & v, char const * msg ) const {
    typedef unsigned long ul;
    switch (m_data_type) {
    case GC_BOOL:
      v = m_data.b?1:0;
      break;
    case GC_INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        msg << "value '" << m_data.i << "' cannot be converted in `unsigned long'"
      )
      v = ul(m_data.i);
      break;
    case GC_LONG:
      GC_ASSERT(
        m_data.l >= 0,
        msg << "value '" << m_data.l << "' cannot be converted in `unsigned long'"
      )
      v = ul(m_data.l);
      break;
    case GC_REAL:
      GC_ASSERT(
        m_data.r >= 0 && isUnsigned64(m_data.r),
        msg << "value '" << m_data.r << "' cannot be converted in `unsigned long'"
      )
      v = ul(m_data.r);
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()) && isUnsigned64(m_data.c->real()),
        msg << "value (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") cannot be converted in `unsigned long'"
      )
      v = ul(m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `unsigned long'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( long & v, char const * msg ) const {
    switch (m_data_type) {
    case GC_BOOL:
      v = m_data.b?1:0;
      break;
    case GC_INTEGER:
      v = long(m_data.i);
      break;
    case GC_LONG:
      v = long(m_data.l);
      break;
    case GC_REAL:
      GC_ASSERT(
        isInteger64(m_data.r),
        msg << "value '" << m_data.r << "' cannot be converted in `long'"
      )
      v = long(m_data.r);
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()) && isInteger64(m_data.c->real()),
        msg << "value (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") cannot be converted in `long'"
      )
      v = long(m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `long'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( float & v, char const * msg ) const {
    switch (m_data_type) {
    case GC_BOOL:
      v = float(m_data.b?1:0);
      break;
    case GC_INTEGER:
      v = float(m_data.i);
      break;
    case GC_LONG:
      v = float(m_data.l);
      break;
    case GC_REAL:
      v = float(m_data.r);
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()),
        msg << "value (" << m_data.c->real() << " " << m_data.c->imag() <<
        ") cannot be converted in `float'"
      )
      v = float(m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `float'"
      )
    }
  }

  template <>
  void
  GenericContainer::get_value( double & v, char const * msg ) const {
    switch (m_data_type) {
    case GC_BOOL:
      v = m_data.b?1:0;
      break;
    case GC_INTEGER:
      v = double(m_data.i);
      break;
    case GC_LONG:
      v = double(m_data.l);
      break;
    case GC_REAL:
      v = m_data.r;
      break;
    case GC_COMPLEX:
      GC_ASSERT(
        isZero(m_data.c->imag()),
        msg << "value (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") cannot be converted in `double'"
      )
      v = m_data.c->real();
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `double'"
      )
    }
  }

  real_type
  GenericContainer::get_number() const {
    switch (m_data_type) {
    case GC_BOOL:    return real_type(m_data.b?1:0);
    case GC_INTEGER: return real_type(m_data.i);
    case GC_LONG:    return real_type(m_data.l);
    case GC_REAL:    return m_data.r;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    break;
    }
    return 0;
  }

  complex_type
  GenericContainer::get_complex_number() const {
    switch (m_data_type) {
    case GC_BOOL:    return complex_type(m_data.b?1:0);
    case GC_INTEGER: return complex_type(real_type(m_data.i),0);
    case GC_LONG:    return complex_type(real_type(m_data.l),0);
    case GC_REAL:    return complex_type(m_data.r,0);
    case GC_COMPLEX: return *m_data.c;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    break;
    }
    return 0;
  }

  void
  GenericContainer::get_complex_number( real_type & re, real_type & im ) const {
    complex_type tmp = get_complex_number();
    re = tmp.real();
    im = tmp.imag();
  }

  real_type
  GenericContainer::get_number_at( unsigned i ) const {
    switch (m_data_type) {
    case GC_VEC_BOOL:    return (*m_data.v_b)[i] ? 1 : 0;
    case GC_VEC_INTEGER: return real_type((*m_data.v_i)[i]);
    case GC_VEC_LONG:    return real_type((*m_data.v_l)[i]);
    case GC_VEC_REAL:    return real_type((*m_data.v_r)[i]);
    case GC_MAT_INTEGER: return real_type((*m_data.m_i)[i]);
    case GC_MAT_LONG:    return real_type((*m_data.m_l)[i]);
    case GC_MAT_REAL:    return (*m_data.m_r)[i];
    case GC_VECTOR:      return (*m_data.v)[i].get_number();
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_COMPLEX:
    case GC_MAP:
    break;
    }
    return 0;
  }

  complex_type
  GenericContainer::get_complex_number_at( unsigned i ) const {
    switch (m_data_type) {
    case GC_VEC_BOOL:    return complex_type(real_type( (*m_data.v_b)[i]?1:0 ),0);
    case GC_VEC_INTEGER: return complex_type(real_type( (*m_data.v_i)[i] ),0);
    case GC_VEC_LONG:    return complex_type(real_type( (*m_data.v_l)[i] ),0);
    case GC_VEC_REAL:    return complex_type((*m_data.v_r)[i],0);
    case GC_VEC_COMPLEX: return (*m_data.v_c)[i];
    case GC_MAT_INTEGER: return complex_type(real_type( (*m_data.m_i)[i] ),0);
    case GC_MAT_LONG:    return complex_type(real_type( (*m_data.m_l)[i] ),0);
    case GC_MAT_REAL:    return complex_type((*m_data.m_r)[i],0);
    case GC_MAT_COMPLEX: return (*m_data.m_c)[i];
    case GC_VECTOR:      return (*m_data.v)[i].get_complex_number();
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    break;
    }
    return 0;
  }

  void
  GenericContainer::get_complex_number_at( unsigned i, real_type & re, real_type & im ) const {
    complex_type tmp = get_complex_number_at(i);
    re = tmp.real();
    im = tmp.imag();
  }

  bool_type &
  GenericContainer::get_bool( char const * msg ) {
    ck_or_set(msg,GC_BOOL);
    return m_data.b;
  }

  bool_type const &
  GenericContainer::get_bool( char const * msg ) const {
    ck(msg,GC_BOOL);
    return m_data.b;
  }

  int_type &
  GenericContainer::get_int( char const * msg ) {
    ck_or_set(msg,GC_INTEGER);
    return m_data.i;
  }

  int_type const &
  GenericContainer::get_int( char const * msg ) const {
    ck(msg,GC_INTEGER);
    return m_data.i;
  }

  long_type &
  GenericContainer::get_long( char const * msg ) {
    ck_or_set(msg,GC_LONG);
    return m_data.l;
  }

  long_type const &
  GenericContainer::get_long( char const * msg ) const {
    ck(msg,GC_LONG);
    return m_data.l;
  }

  int_type
  GenericContainer::get_as_int( char const * msg_in ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    switch (m_data_type) {
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "bad data type: " << typeName[m_data_type] <<
        " can't be converted to int"
      )
    case GC_BOOL:
      return m_data.b ? 1 : 0;
    case GC_INTEGER:
      return m_data.i;
    case GC_LONG:
      GC_ASSERT(
        (std::abs(m_data.l)>>31) == 0,
        msg << "type: " << typeName[m_data_type] <<
        " value = " << m_data.l << " can't be converted to int"
      )
      return int_type(m_data.l);
    case GC_REAL:
      GC_ASSERT(
        isInteger32(m_data.r),
        msg << "type: " << typeName[m_data_type] <<
        " value = " << m_data.r << " can't be converted to int"
      )
      return int_type(m_data.r);
    case GC_COMPLEX:
      GC_ASSERT(
        isInteger32(m_data.c->real()) && isZero(m_data.c->imag()),
        msg << "type: " << typeName[m_data_type] <<
        " value = (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") can't be converted to int"
      )
      return int_type(m_data.c->real());
    }
    return 0;
  }

  uint_type
  GenericContainer::get_as_uint( char const * msg_in ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    switch (m_data_type) {
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "bad data type: " << typeName[m_data_type] <<
        " can't be converted to uint"
      )
    case GC_BOOL:
      return m_data.b ? 1 : 0;
    case GC_INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        msg << "type: " << typeName[m_data_type] <<
        " value = " << m_data.i << " can't be converted to uint"
      )
      return uint_type(m_data.i);
    case GC_LONG:
      GC_ASSERT(
        (std::abs(m_data.l)>>32) == 0 && m_data.l >= 0,
        msg << "type: " << typeName[m_data_type] <<
        " value = " << m_data.l << " can't be converted to uint"
      )
      return uint_type(m_data.l);
    case GC_REAL:
      GC_ASSERT(
        isUnsigned32(m_data.r),
        msg << "type: " << typeName[m_data_type] <<
        " value = " << m_data.r << " can't be converted to uint"
      )
      return uint_type(m_data.r);
    case GC_COMPLEX:
      GC_ASSERT(
        isUnsigned32(m_data.c->real()) && isZero(m_data.c->imag()),
        msg << "type: " << typeName[m_data_type] <<
        "v alue = (" << m_data.c->real() << "," <<  m_data.c->imag() <<
        ") can't be converted to uint"
      )
      return uint_type(m_data.c->real());
    }
    return 0;
  }

  long_type
  GenericContainer::get_as_long( char const * msg_in ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    switch (m_data_type) {
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "bad data type: " << typeName[m_data_type] <<
        " can't be converted to long"
      )
    case GC_BOOL:
      return m_data.b ? 1 : 0;
    case GC_INTEGER:
      return long_type(m_data.i);
    case GC_LONG:
      return m_data.l;
    case GC_REAL:
      GC_ASSERT(
        isInteger64(m_data.r),
        msg << "value = " << m_data.r << " can't be converted to long"
      )
      return long_type(m_data.r);
    case GC_COMPLEX:
      GC_ASSERT(
        isInteger64(m_data.c->real()) && isZero(m_data.c->imag()),
        msg << "value = (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") can't be converted to long"
      )
      return long_type(m_data.c->real());
    }
    return 0;
  }

  ulong_type
  GenericContainer::get_as_ulong( char const * msg_in ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    switch (m_data_type) {
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "bad data type: " << typeName[m_data_type] <<
        " can't be converted to ulong"
      )
    case GC_BOOL:
      return m_data.b ? 1 : 0;
    case GC_INTEGER:
      GC_ASSERT(
        m_data.i >= 0,
        msg << "value = " << m_data.i << " can't be converted to ulong"
      )
      return ulong_type(m_data.i);
    case GC_LONG:
      GC_ASSERT(
        m_data.l >= 0,
        msg << "value = " << m_data.l << " can't be converted to ulong"
      )
      return ulong_type(m_data.l);
    case GC_REAL:
      GC_ASSERT(
        isUnsigned64(m_data.r),
        msg << "value = " << m_data.r << " can't be converted to ulong"
      )
      return ulong_type(m_data.r);
    case GC_COMPLEX:
      GC_ASSERT(
        isUnsigned64(m_data.c->real()) && isZero(m_data.c->imag()),
        msg << "value = (" << m_data.c->real() << "," << m_data.c->imag() <<
        ") can't be converted to ulong"
      )
      return ulong_type(m_data.c->real());
    }
    return 0;
  }


  real_type &
  GenericContainer::get_real( char const * msg ) {
    ck_or_set(msg,GC_REAL);
    return m_data.r;
  }

  real_type const &
  GenericContainer::get_real( char const * msg ) const {
    ck(msg,GC_REAL);
    return m_data.r;
  }

  complex_type &
  GenericContainer::get_complex( char const * msg ) {
    ck_or_set(msg,GC_COMPLEX);
    return *m_data.c;
  }

  complex_type const &
  GenericContainer::get_complex( char const * msg ) const {
    ck(msg,GC_COMPLEX);
    return *m_data.c;
  }

  string_type &
  GenericContainer::get_string( char const * msg ) {
    ck_or_set(msg,GC_STRING);
    return *m_data.s;
  }

  string_type const &
  GenericContainer::get_string( char const * msg ) const {
    ck(msg,GC_STRING);
    return *m_data.s;
  }

  vector_type &
  GenericContainer::get_vector( char const * msg ) {
    ck(msg,GC_VECTOR);
    return *m_data.v;
  }

  vector_type const &
  GenericContainer::get_vector( char const * msg ) const {
    ck(msg,GC_VECTOR);
    return *m_data.v;
  }

  vec_pointer_type &
  GenericContainer::get_vec_pointer( char const * msg ) {
    ck(msg,GC_VEC_POINTER);
    return *m_data.v_p;
  }

  vec_pointer_type const &
  GenericContainer::get_vec_pointer( char const * msg ) const {
    ck(msg,GC_VEC_POINTER);
    return *m_data.v_p;
  }

  vec_bool_type &
  GenericContainer::get_vec_bool( char const * msg ) {
    ck(msg,GC_VEC_BOOL);
    return *m_data.v_b;
  }

  vec_bool_type const &
  GenericContainer::get_vec_bool( char const * msg ) const {
    ck(msg,GC_VEC_BOOL);
    return *m_data.v_b;
  }

  vec_int_type &
  GenericContainer::get_vec_int( char const * msg ) {
    if ( m_data_type == GC_NOTYPE   ) set_vec_int();
    if ( m_data_type == GC_VEC_BOOL ) promote_to_vec_int();
    ck(msg,GC_VEC_INTEGER);
    return *m_data.v_i;
  }

  vec_int_type const &
  GenericContainer::get_vec_int( char const * msg ) const {
    ck(msg,GC_VEC_INTEGER);
    return *m_data.v_i;
  }

  vec_long_type &
  GenericContainer::get_vec_long( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_vec_long();
    if ( m_data_type == GC_VEC_BOOL || m_data_type == GC_VEC_INTEGER ) promote_to_vec_long();
    ck(msg,GC_VEC_LONG);
    return *m_data.v_l;
  }

  vec_long_type const &
  GenericContainer::get_vec_long( char const * msg ) const {
    ck(msg,GC_VEC_LONG);
    return *m_data.v_l;
  }

  vec_real_type &
  GenericContainer::get_vec_real( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_vec_real();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG ) promote_to_vec_real();
    ck(msg,GC_VEC_REAL);
    return *m_data.v_r;
  }

  vec_real_type const &
  GenericContainer::get_vec_real( char const * msg ) const {
    ck(msg,GC_VEC_REAL);
    return *m_data.v_r;
  }

  vec_complex_type &
  GenericContainer::get_vec_complex( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_vec_complex();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG    ||
         m_data_type == GC_VEC_REAL ) promote_to_vec_complex();
    ck(msg,GC_VEC_COMPLEX);
    return *m_data.v_c;
  }

  vec_complex_type const &
  GenericContainer::get_vec_complex( char const * msg ) const {
    ck(msg,GC_VEC_COMPLEX);
    return *m_data.v_c;
  }

  mat_int_type &
  GenericContainer::get_mat_int( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_mat_int();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG    ||
         m_data_type == GC_VEC_REAL ) promote_to_mat_int();
    ck(msg,GC_MAT_INTEGER);
    return *m_data.m_i;
  }

  mat_int_type const &
  GenericContainer::get_mat_int( char const * msg ) const {
    ck(msg,GC_MAT_INTEGER);
    return *m_data.m_i;
  }

  mat_long_type &
  GenericContainer::get_mat_long( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_mat_long();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG    ||
         m_data_type == GC_VEC_REAL    ||
         m_data_type == GC_MAT_INTEGER ) promote_to_mat_long();
    ck(msg,GC_MAT_LONG);
    return *m_data.m_l;
  }

  mat_long_type const &
  GenericContainer::get_mat_long( char const * msg ) const {
    ck(msg,GC_MAT_LONG);
    return *m_data.m_l;
  }

  mat_real_type &
  GenericContainer::get_mat_real( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_mat_real();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG    ||
         m_data_type == GC_VEC_REAL    ||
         m_data_type == GC_MAT_INTEGER ||
         m_data_type == GC_MAT_LONG ) promote_to_mat_real();
    ck(msg,GC_MAT_REAL);
    return *m_data.m_r;
  }

  mat_real_type const &
  GenericContainer::get_mat_real( char const * msg ) const {
    ck(msg,GC_MAT_REAL);
    return *m_data.m_r;
  }

  mat_complex_type &
  GenericContainer::get_mat_complex( char const * msg ) {
    if ( m_data_type == GC_NOTYPE ) set_mat_complex();
    if ( m_data_type == GC_VEC_BOOL    ||
         m_data_type == GC_VEC_INTEGER ||
         m_data_type == GC_VEC_LONG    ||
         m_data_type == GC_VEC_REAL    ||
         m_data_type == GC_MAT_REAL    ||
         m_data_type == GC_VEC_COMPLEX ) promote_to_mat_complex();
    ck(msg,GC_MAT_COMPLEX);
    return *m_data.m_c;
  }

  mat_complex_type const &
  GenericContainer::get_mat_complex( char const * msg ) const {
    ck(msg,GC_MAT_COMPLEX);
    return *m_data.m_c;
  }

  vec_string_type &
  GenericContainer::get_vec_string( char const * msg ) {
    ck(msg,GC_VEC_STRING);
    return *m_data.v_s;
  }

  vec_string_type const &
  GenericContainer::get_vec_string( char const * msg ) const {
    ck(msg,GC_VEC_STRING);
    return *m_data.v_s;
  }

  map_type &
  GenericContainer::get_map( char const * msg ) {
    ck(msg,GC_MAP);
    return *m_data.m;
  }

  map_type const &
  GenericContainer::get_map( char const * msg ) const {
    ck(msg,GC_MAP);
    return *m_data.m;
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_int( vec_int_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    long_type    lval;
    real_type    rval;
    complex_type cval;
    int_type     val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = m_data.b ? 1 : 0;
        break;
      case GC_INTEGER:
        val = m_data.i;
        break;
      case GC_LONG:
        lval = m_data.l;
        GC_ASSERT(
          int_type(lval) == lval,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `integer'"
        );
        val = int_type(lval);
        break;
      case GC_REAL:
        rval = m_data.r;
        GC_ASSERT(
          isInteger32(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `integer'"
        )
        val = int_type(rval);
        break;
      case GC_VEC_BOOL:
        val = (*m_data.v_b)[i] ? 1 : 0;
        break;
      case GC_VEC_INTEGER:
        val = (*m_data.v_i)[i];
        break;
      case GC_VEC_LONG:
        lval = (*m_data.v_l)[i];
        GC_ASSERT(
          int_type(lval) == lval,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `integer'"
        )
        val = int_type(lval);
        break;
      case GC_VEC_REAL:
        rval = (*m_data.v_r)[i];
        GC_ASSERT(
          isInteger32(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `integer'"
        )
        val = int_type(rval);
        break;
      case GC_COMPLEX:
        cval = *m_data.c;
        GC_ASSERT(
          isZero(cval.imag()) && isInteger32(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `integer'"
        )
        val = int_type(cval.real());
        break;
      case GC_VEC_COMPLEX:
        cval = (*m_data.v_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isInteger32(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `integer'"
        )
        val = int_type(cval.real());
        break;
      case GC_MAT_INTEGER:
        val = (*m_data.m_i)[i];
        break;
      case GC_MAT_LONG:
        lval = (*m_data.m_l)[i];
        GC_ASSERT(
          int_type(lval) == lval,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `integer'"
        );
        val = int_type(lval);
        break;
      case GC_MAT_REAL:
        rval = (*m_data.m_r)[i];
        GC_ASSERT(
          isInteger32(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `integer'"
        )
        val = int_type(rval);
        break;
      case GC_MAT_COMPLEX:
        cval = (*m_data.m_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isInteger32(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `integer'"
        )
        val = int_type(cval.real());
        break;
      case GC_VECTOR:
        val = (*this)(i).get_as_int("GenericContainer::copyto_vec_int ");
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_int_type'"
        )
      }
      v.push_back(val);
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_uint( vec_uint_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    uint_type    val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = m_data.b ? 1 : 0;
        break;
      case GC_INTEGER:
        ival = m_data.i;
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned integer'"
        );
        val = uint_type(ival);
        break;
      case GC_LONG:
        lval = m_data.l;
        GC_ASSERT(
          int_type(lval) == lval && lval >= 0,
          msg << "copyto_vec_uint: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(lval);
        break;
      case GC_REAL:
        rval = m_data.r;
        GC_ASSERT(
          isUnsigned32(rval),
          msg << "copyto_vec_uint: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(rval);
        break;
      case GC_VEC_BOOL:
        val = (*m_data.v_b)[i] ? 1 : 0;
        break;
      case GC_VEC_INTEGER:
        ival = (*m_data.v_i)[i];
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned integer'"
        );
        val = uint_type(ival);
        break;
      case GC_VEC_LONG:
        lval = (*m_data.v_l)[i];
        GC_ASSERT(
          int_type(lval) == lval && lval >= 0,
          msg << "copyto_vec_uint: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(lval);
        break;
      case GC_VEC_REAL:
        rval = (*m_data.v_r)[i];
        GC_ASSERT(
          isUnsigned32(rval),
          msg << "copyto_vec_uint: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(rval);
        break;
      case GC_COMPLEX:
        cval = *m_data.c;
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned32(cval.real()),
          msg << "copyto_vec_uint: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned integer'"
        )
        val = uint_type(cval.real());
        break;
      case GC_VEC_COMPLEX:
        cval = (*m_data.v_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned32(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned integer'"
        );
        val = uint_type(cval.real());
        break;
      case GC_MAT_INTEGER:
        ival = (*m_data.m_i)[i];
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(ival);
        break;
      case GC_MAT_LONG:
        lval = (*m_data.m_l)[i];
        GC_ASSERT(
          int_type(lval) == lval && lval >= 0,
          msg << "copyto_vec_uint: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(lval);
        break;
      case GC_MAT_REAL:
        rval = (*m_data.m_r)[i];
        GC_ASSERT(
          isUnsigned32(rval),
          msg << "copyto_vec_uint: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned integer'"
        )
        val = uint_type(rval);
        break;
      case GC_MAT_COMPLEX:
        cval = (*m_data.m_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned32(cval.real()),
          msg << "copyto_vec_uint: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned integer'"
        )
        val = uint_type(cval.real());
        break;
      case GC_VECTOR:
        val = (*this)(i).get_as_uint("GenericContainer::copyto_vec_uint ");
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_uint_type'"
        )
      }
      v.push_back(val);
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_long( vec_long_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    real_type    rval;
    complex_type cval;
    long_type    val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = m_data.b ? 1 : 0;
        break;
      case GC_INTEGER:
        val = long_type(m_data.i);
        break;
      case GC_LONG:
        val = m_data.l;
        break;
      case GC_REAL:
        rval = m_data.r;
        GC_ASSERT(
          isInteger64(rval),
          msg << "copyto_vec_long: v[" << i << "] = " << rval <<
          " cannot be converted to `long'"
        )
        val = long_type(rval);
        break;
      case GC_VEC_BOOL:
        val = (*m_data.v_b)[i] ? 1 : 0;
        break;
      case GC_VEC_INTEGER:
        val = long_type((*m_data.v_i)[i]);
        break;
      case GC_VEC_LONG:
        val = (*m_data.v_l)[i];
        break;
      case GC_VEC_REAL:
        rval = (*m_data.v_r)[i];
        GC_ASSERT(
          isInteger64(rval),
          msg << "copyto_vec_long: v[" << i << "] = " << rval <<
          " cannot be converted to `long'"
        )
        val = long_type(rval);
        break;
      case GC_COMPLEX:
        cval = *m_data.c;
        GC_ASSERT(
          isZero(cval.imag()) && isInteger64(cval.real()),
          msg << "copyto_vec_long: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `long'"
        )
        val = long_type(cval.real());
        break;
      case GC_VEC_COMPLEX:
        cval = (*m_data.v_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isInteger64(cval.real()),
          msg << "copyto_vec_long: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `long'"
        )
        val = long_type(cval.real());
        break;
      case GC_MAT_INTEGER:
        val = long_type((*m_data.m_i)[i]);
        break;
      case GC_MAT_LONG:
        val = (*m_data.m_l)[i];
        break;
      case GC_MAT_REAL:
        rval = (*m_data.m_r)[i];
        GC_ASSERT(
          isInteger64(rval),
          msg << "copyto_vec_long: v[" << i << "] = " << rval <<
          " cannot be converted to `long'"
        )
        val = long_type(rval);
        break;
      case GC_MAT_COMPLEX:
        cval = (*m_data.m_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isInteger64(cval.real()),
          msg << "copyto_vec_long: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `long'"
        )
        val = long_type(cval.real());
        break;
      case GC_VECTOR:
        val = (*this)(i).get_as_long("GenericContainer::copyto_vec_long ");
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_long_type'"
        )
      }
      v.push_back(val);
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_ulong( vec_ulong_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    int_type     ival;
    long_type    lval;
    real_type    rval;
    complex_type cval;
    ulong_type   val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = m_data.b ? 1 : 0;
        break;
      case GC_INTEGER:
        ival = m_data.i;
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned long'"
        );
        val = ulong_type(ival);
        break;
      case GC_LONG:
        lval = m_data.l;
        GC_ASSERT(
          lval >= 0,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(lval);
        break;
      case GC_REAL:
        rval = m_data.r;
        GC_ASSERT(
          isUnsigned64(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(rval);
        break;
      case GC_VEC_BOOL:
        val = (*m_data.v_b)[i] ? 1 : 0;
        break;
      case GC_VEC_INTEGER:
        ival = (*m_data.v_i)[i];
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(ival);
        break;
      case GC_VEC_LONG:
        lval = (*m_data.v_l)[i];
        GC_ASSERT(
          lval >= 0,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned long'"
        );
        val = ulong_type(lval);
        break;
      case GC_VEC_REAL:
        rval = (*m_data.v_r)[i];
        GC_ASSERT(
          isUnsigned64(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(rval);
        break;
      case GC_COMPLEX:
        cval = *m_data.c;
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned64(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned long'"
        )
        val = ulong_type(cval.real());
        break;
      case GC_VEC_COMPLEX:
        cval = (*m_data.v_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned64(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned long'"
        );
        val = ulong_type(cval.real());
        break;
      case GC_MAT_INTEGER:
        ival = (*m_data.m_i)[i];
        GC_ASSERT(
          ival >= 0,
          msg << "copyto_vec_uint: value = " << ival <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(ival);
        break;
      case GC_MAT_LONG:
        lval = (*m_data.m_l)[i];
        GC_ASSERT(
          lval >= 0,
          msg << "copyto_vec_int: v[" << i << "] = " << lval <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(lval);
        break;
      case GC_MAT_REAL:
        rval = (*m_data.m_r)[i];
        GC_ASSERT(
          isUnsigned64(rval),
          msg << "copyto_vec_int: v[" << i << "] = " << rval <<
          " cannot be converted to `unsigned long'"
        )
        val = ulong_type(rval);
        break;
      case GC_MAT_COMPLEX:
        cval = (*m_data.m_c)[i];
        GC_ASSERT(
          isZero(cval.imag()) && isUnsigned64(cval.real()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `unsigned long'"
        )
        val = ulong_type(cval.real());
        break;
      case GC_VECTOR:
        val = (*this)(i).get_as_ulong("GenericContainer::copyto_vec_ulong ");
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_ulong_type'"
        )
      }
      v.push_back(val);
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_real( vec_real_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    complex_type cval;
    real_type    val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = m_data.b ? 1 : 0;
        break;
      case GC_INTEGER:
        val = real_type(m_data.i);
        break;
      case GC_LONG:
        val = real_type(m_data.l);
        break;
      case GC_REAL:
        val = m_data.r;
        break;
      case GC_VEC_BOOL:
        val = (*m_data.v_b)[i] ? 1 : 0;
        break;
      case GC_VEC_INTEGER:
        val = real_type((*m_data.v_i)[i]);
        break;
      case GC_VEC_LONG:
        val = real_type((*m_data.v_l)[i]);
        break;
      case GC_VEC_REAL:
        val = (*m_data.v_r)[i];
        break;
      case GC_COMPLEX:
        cval = *m_data.c;
        GC_ASSERT(
          isZero(cval.imag()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `real'"
        )
        val = cval.real();
        break;
      case GC_VEC_COMPLEX:
        cval = (*m_data.v_c)[i];
        GC_ASSERT(
          isZero(cval.imag()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `real'"
        )
        val = cval.real();
        break;
      case GC_MAT_INTEGER:
        val = real_type((*m_data.m_i)[i]);
        break;
      case GC_MAT_LONG:
        val = real_type((*m_data.m_l)[i]);
        break;
      case GC_MAT_REAL:
        val = (*m_data.m_r)[i];
        break;
      case GC_MAT_COMPLEX:
        cval = (*m_data.m_c)[i];
        GC_ASSERT(
          isZero(cval.imag()),
          msg << "copyto_vec_int: v[" << i << "] = (" <<
          cval.real() << "," << cval.imag() <<
          ") cannot be converted to `real'"
        )
        val = cval.real();
        break;
      case GC_VECTOR:
        val = (*this)(i).get_real();
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_real_type'"
        )
      }
      v.push_back(val);
    }
  }


  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void
  GenericContainer::copyto_vec_complex( vec_complex_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    v.reserve(ne);
    complex_type val = 0;
    v.reserve(ne);
    for ( unsigned i = 0; i < ne; ++i ) {
      switch (m_data_type) {
      case GC_BOOL:
        val = real_type(m_data.b ? 1 : 0);
        break;
      case GC_INTEGER:
        val = complex_type(real_type( m_data.i ),0);
        break;
      case GC_LONG:
        val = complex_type(real_type( m_data.l ),0);
        break;
      case GC_REAL:
        val = complex_type(m_data.r,0);
        break;
      case GC_VEC_BOOL:
        val = complex_type( real_type( (*m_data.v_b)[i] ? 1 : 0), 0 );
        break;
      case GC_VEC_INTEGER:
        val = complex_type( real_type( (*m_data.v_i)[i] ),0);
        break;
      case GC_VEC_LONG:
        val = complex_type( real_type( (*m_data.v_l)[i] ),0);
        break;
      case GC_VEC_REAL:
        val = complex_type((*m_data.v_r)[i],0);
        break;
      case GC_COMPLEX:
        val = *m_data.c;
        break;
      case GC_VEC_COMPLEX:
        val = (*m_data.v_c)[i];
        break;
      case GC_MAT_INTEGER:
        val = complex_type( real_type( (*m_data.m_i)[i] ),0);
        break;
      case GC_MAT_LONG:
        val = complex_type( real_type( (*m_data.m_l)[i] ),0);
        break;
      case GC_MAT_REAL:
        val = complex_type((*m_data.m_r)[i],0);
        break;
      case GC_MAT_COMPLEX:
        val = (*m_data.m_c)[i];
        break;
      case GC_VECTOR:
        val = (*this)(i).get_complex();
        break;
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_STRING:
      case GC_VEC_POINTER:
      case GC_VEC_STRING:
      case GC_MAP:
        GC_DO_ERROR(
          msg << "\nbad data type: `" << typeName[m_data_type] <<
          "' cannot be converted in `vec_complex_type'"
        )
      }
      v.push_back(val);
    }
  }

  void
  GenericContainer::copyto_vec_string( vec_string_type & v, char const * msg ) const {
    v.clear();
    unsigned ne = get_num_elements();
    switch (m_data_type) {
    case GC_STRING:
      v.reserve(ne);
      v.push_back( *m_data.s );
      break;
    case GC_VEC_STRING:
      v.resize(ne);
      std::copy( m_data.v_s->begin(), m_data.v_s->end(), v.begin() );
      break;
    case GC_VECTOR:
      v.reserve(ne);
      for ( unsigned i = 0; i < ne; ++i ) {
        GenericContainer const & gc = get_gc_at(i,msg);
        v.push_back( gc.get_string(msg) );
      }
      break;
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_MAT_INTEGER:
    case GC_VEC_LONG:
    case GC_MAT_LONG:
    case GC_VEC_REAL:
    case GC_MAT_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
    case GC_POINTER:
    case GC_VEC_POINTER:
    case GC_MAP:
      GC_DO_ERROR(
        msg << "\nbad data type: `" << typeName[m_data_type] <<
        "' cannot be converted in `vec_string_type'"
      )
    }
  }

  // --------------------------------------------------------------

  bool
  GenericContainer::exists( std::string const & s ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(s);
    return iv != (*m_data.m).end();
  }

  bool
  GenericContainer::get_if_exists( char const * field, bool & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    if ( iv->second.m_data_type != GC_BOOL ) return false;
    value = iv->second.m_data.b;
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, int_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_INTEGER:
      value = iv->second.m_data.i;
      break;
    case GC_LONG:
      value = int_type(iv->second.m_data.l);
      break;
    case GC_REAL:
      if ( !isInteger32(iv->second.m_data.r) ) return false;
      value = int_type(iv->second.m_data.r);
      break;
    case GC_COMPLEX:
      if ( ! ( isInteger32(iv->second.m_data.c->real()) &&
               isZero(iv->second.m_data.c->imag()) ) ) return false;
      value = int_type(iv->second.m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, uint_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_INTEGER:
      if ( iv->second.m_data.i < 0 ) return false;
      value = uint_type(iv->second.m_data.i);
      break;
    case GC_LONG:
      if ( iv->second.m_data.l < 0 ) return false;
      value = uint_type(iv->second.m_data.l);
      break;
    case GC_REAL:
      if ( ! isUnsigned32(iv->second.m_data.r) ) return false;
      value = uint_type(iv->second.m_data.r);
      break;
    case GC_COMPLEX:
      if ( ! ( isUnsigned32(iv->second.m_data.c->real()) &&
               isZero(iv->second.m_data.c->imag()) ) ) return false;
      value = uint_type(iv->second.m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, long_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_INTEGER:
      value = long_type(iv->second.m_data.i);
      break;
    case GC_LONG:
      value = iv->second.m_data.l;
      break;
    case GC_REAL:
      if ( ! isInteger64(iv->second.m_data.r) ) return false;
      value = long_type(iv->second.m_data.r);
      break;
    case GC_COMPLEX:
      if ( ! ( isInteger64(iv->second.m_data.c->real()) &&
               isZero(iv->second.m_data.c->imag()) ) ) return false;
      value = long_type(iv->second.m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, ulong_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = iv->second.m_data.b?1:0;
      break;
    case GC_INTEGER:
      if ( iv->second.m_data.i < 0 ) return false;
      value = ulong_type(iv->second.m_data.i);
      break;
    case GC_LONG:
      if ( iv->second.m_data.l < 0 ) return false;
      value = ulong_type(iv->second.m_data.l);
      break;
    case GC_REAL:
      if ( ! isUnsigned64(iv->second.m_data.r) ) return false;
      value = ulong_type(iv->second.m_data.r);
      break;
    case GC_COMPLEX:
      if ( ! ( isUnsigned64(iv->second.m_data.c->real()) &&
               isZero(iv->second.m_data.c->imag()) ) ) return false;
      value = ulong_type(iv->second.m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, real_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = real_type(iv->second.m_data.b?1:0);
      break;
    case GC_INTEGER:
      value = real_type(iv->second.m_data.i);
      break;
    case GC_LONG:
      value = real_type(iv->second.m_data.l);
      break;
    case GC_REAL:
      value = iv->second.m_data.r;
      break;
    case GC_COMPLEX:
      if ( ! isZero(iv->second.m_data.c->imag()) ) return false;
      value = real_type(iv->second.m_data.c->real());
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, complex_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    switch (iv->second.m_data_type) {
    case GC_BOOL:
      value = complex_type(iv->second.m_data.b?1:0,0);
      break;
    case GC_INTEGER:
      value = complex_type(real_type( iv->second.m_data.i),0);
      break;
    case GC_LONG:
      value = complex_type(real_type( iv->second.m_data.l),0);
      break;
    case GC_REAL:
      value = complex_type(iv->second.m_data.r,0);
      break;
    case GC_COMPLEX:
      value = *iv->second.m_data.c;
      break;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false;
    }
    return true;
  }

  bool
  GenericContainer::get_if_exists( char const * field, string_type & value ) const {
    if ( m_data_type != GC_MAP ) return false;
    map_type::iterator iv = (*m_data.m).find(field);
    if ( iv == (*m_data.m).end() ) return false;
    if ( iv->second.m_data_type != GC_STRING ) return false;
    value = *iv->second.m_data.s;
    return true;
  }

  // --------------------------------------------------------------
  bool_type
  GenericContainer::get_bool_at( unsigned i ) {
    if ( m_data_type == GC_NOTYPE   ) set_vec_bool();
    if ( m_data_type == GC_VEC_BOOL ) {
      CHECK_RESIZE(m_data.v_b,i); // correct type, check size
      return (*m_data.v_b)[i];
    } else {
      if ( m_data_type != GC_VECTOR ) promote_to_vector();
      CHECK_RESIZE(m_data.v,i);
      return (*m_data.v)[i].set_bool(false);
    }
  }

  bool_type
  GenericContainer::get_bool_at( unsigned i, char const * msg ) const {
    ck(msg,GC_VEC_BOOL);
    GC_ASSERT(
      i < m_data.v_b->size(),
      msg << "\nget_bool_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_b)[i];
  }

  int_type &
  GenericContainer::get_int_at( unsigned i ) {
    if      ( m_data_type == GC_NOTYPE ) set_vec_int();
    else if ( m_data_type == GC_BOOL    ||
              m_data_type == GC_INTEGER ||
              m_data_type == GC_VEC_BOOL ) promote_to_vec_int();
    if ( m_data_type == GC_VEC_INTEGER ) {
      CHECK_RESIZE(m_data.v_i,i); // correct type, check size
      return (*m_data.v_i)[i];
    } else {
      if ( m_data_type != GC_VECTOR ) promote_to_vector();
      CHECK_RESIZE(m_data.v,i);
      return (*m_data.v)[i].set_int(0);
    }
  }

  int_type const &
  GenericContainer::get_int_at( unsigned i, char const * msg ) const {
    ck(msg,GC_VEC_INTEGER);
    GC_ASSERT(
      i < m_data.v_i->size(),
      msg << "\nget_int_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_i)[i];
  }

  int_type &
  GenericContainer::get_int_at( unsigned i, unsigned j ) {
    if      ( m_data_type == GC_NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ) promote_to_mat_int();
    GC_ASSERT(
      GC_MAT_INTEGER == m_data_type,
      "get_int_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << typeName[GC_MAT_INTEGER] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    return (*m_data.m_i)(i,j);
  }

  int_type const &
  GenericContainer::get_int_at( unsigned i, unsigned j, char const * msg ) const  {
    ck(msg,GC_MAT_INTEGER);
    GC_ASSERT(
      i < m_data.m_i->numRows() && j < m_data.m_i->numCols(),
      msg << "\bget_int_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_i)(i,j);
  }

  long_type &
  GenericContainer::get_long_at( unsigned i ) {
    if      ( m_data_type == GC_NOTYPE ) set_vec_long();
    else if ( m_data_type == GC_BOOL     ||
              m_data_type == GC_INTEGER  ||
              m_data_type == GC_LONG     ||
              m_data_type == GC_VEC_BOOL ||
              m_data_type == GC_VEC_INTEGER ) promote_to_vec_long();
    if ( m_data_type == GC_VEC_LONG ) {
      CHECK_RESIZE(m_data.v_l,i); // correct type, check size
      return (*m_data.v_l)[i];
    } else {
      if ( m_data_type != GC_VECTOR ) promote_to_vector();
      CHECK_RESIZE(m_data.v,i);
      return (*m_data.v)[i].set_long(0);
    }
  }

  long_type const &
  GenericContainer::get_long_at( unsigned i, char const * msg ) const {
    ck(msg,GC_VEC_LONG);
    GC_ASSERT(
      i < m_data.v_l->size(),
      msg << "\nget_long_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_l)[i];
  }

  long_type &
  GenericContainer::get_long_at( unsigned i, unsigned j ) {
    if      ( m_data_type == GC_NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ||
              m_data_type == GC_VEC_LONG ) promote_to_mat_long();
    GC_ASSERT(
      GC_MAT_LONG == m_data_type,
      "get_long_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << typeName[GC_MAT_LONG] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    return (*m_data.m_l)(i,j);
  }

  long_type const &
  GenericContainer::get_long_at( unsigned i, unsigned j, char const * msg ) const  {
    ck(msg,GC_MAT_LONG);
    GC_ASSERT(
      i < m_data.m_l->numRows() && j < m_data.m_l->numCols(),
      msg << "\bget_long_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_l)(i,j);
  }

  real_type &
  GenericContainer::get_real_at( unsigned i ) {
    if      ( m_data_type == GC_NOTYPE ) set_vec_real();
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ||
              m_data_type == GC_VEC_LONG ) promote_to_vec_real();
    if ( m_data_type == GC_VEC_REAL ) {
      CHECK_RESIZE(m_data.v_r,i); // correct type, check size
      return (*m_data.v_r)[i];
    } else {
      if ( m_data_type != GC_VECTOR ) promote_to_vector();
      CHECK_RESIZE(m_data.v,i);
      return (*m_data.v)[i].set_real(0);
    }
  }

  real_type const &
  GenericContainer::get_real_at( unsigned i, char const * msg ) const  {
    ck(msg,GC_VEC_REAL);
    GC_ASSERT(
      i < m_data.v_r->size(),
      msg << "\nget_real_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_r)[i];
  }

  real_type &
  GenericContainer::get_real_at( unsigned i, unsigned j ) {
    if      ( m_data_type == GC_NOTYPE ) set_mat_real(i,j);
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ||
              m_data_type == GC_VEC_LONG ||
              m_data_type == GC_VEC_REAL ) promote_to_mat_real();
    GC_ASSERT(
      GC_MAT_REAL == m_data_type,
      "get_real_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << typeName[GC_MAT_REAL] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    return (*m_data.m_r)(i,j);
  }

  real_type const &
  GenericContainer::get_real_at( unsigned i, unsigned j, char const * msg ) const  {
    ck(msg,GC_MAT_REAL);
    GC_ASSERT(
      i < m_data.v_r->size(),
      msg << "\bget_real_at( " << i << ", " << j << " ) const, out of range"
    )
    return (*m_data.m_r)(i,j);
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned i ) {
    if      ( m_data_type == GC_NOTYPE ) set_vec_complex();
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ||
              m_data_type == GC_VEC_LONG    ||
              m_data_type == GC_VEC_REAL ) promote_to_vec_complex();
    if ( m_data_type == GC_VEC_COMPLEX ) {
      CHECK_RESIZE(m_data.v_c,i); // correct type, check size
      return (*m_data.v_c)[i];
    } else {
      if ( m_data_type != GC_VECTOR ) promote_to_vector();
      CHECK_RESIZE(m_data.v,i);
      return (*m_data.v)[i].set_complex(0,0);
    }
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned i, char const * msg ) const  {
    ck(msg,GC_VEC_COMPLEX);
    GC_ASSERT(
      i < m_data.v_c->size(),
      msg << "\nget_complex_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_c)[i];
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned i, unsigned j ) {
    if      ( m_data_type == GC_NOTYPE ) set_mat_complex(i,j);
    else if ( m_data_type == GC_VEC_BOOL    ||
              m_data_type == GC_VEC_INTEGER ||
              m_data_type == GC_VEC_LONG    ||
              m_data_type == GC_VEC_REAL    ||
              m_data_type == GC_VEC_COMPLEX ||
              m_data_type == GC_MAT_REAL ) promote_to_mat_complex();
    GC_ASSERT(
      GC_MAT_COMPLEX == m_data_type,
      "get_complex_at( " << i << ", " << j << " ) bad data type" <<
      "\nexpect: " << typeName[GC_MAT_COMPLEX] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    return (*m_data.m_c)(i,j);
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned i, unsigned j, char const * msg ) const  {
    ck(msg,GC_MAT_COMPLEX);
    return (*m_data.m_c)(i,j);
  }

  string_type &
  GenericContainer::get_string_at( unsigned i ) {
    if ( m_data_type == GC_NOTYPE ) set_vec_string();
    if ( m_data_type == GC_VEC_STRING ) {
      CHECK_RESIZE(m_data.v_s,i);
      return (*m_data.v_s)[i];
    } else {
      promote_to_vector();
      return (*this)[i].set_string("");
    }
  }

  string_type const &
  GenericContainer::get_string_at( unsigned i, char const * msg ) const {
    ck(msg,GC_VEC_STRING);
    GC_ASSERT(
      i < m_data.v_s->size(),
      msg << "\nget_string_at( " << i << " ) const, out of range"
    )
    return (*m_data.v_s)[i];
  }

  GenericContainer &
  GenericContainer::get_gc_at( unsigned i )
  { return (*this)[i]; }

  GenericContainer const &
  GenericContainer::get_gc_at( unsigned i, char const * msg ) const {
    return (*this)(i,msg);
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
    case GC_NOTYPE:
      stream << "GenericContainer: No data stored\n";
      break;
    case GC_POINTER:
      stream << "Generic pointer: " << m_data.p << '\n';
      break;
    case GC_BOOL:
      stream << "Boolean: " << (m_data.b?"true":"false") << '\n';
      break;
    case GC_INTEGER:
      stream << "Integer: " << m_data.i << '\n';
      break;
    case GC_LONG:
      stream << "Long: " << m_data.l << '\n';
      break;
    case GC_REAL:
      stream << "Floating Point: " << m_data.r << '\n';
      break;
    case GC_COMPLEX:
      stream << "Complex Floating Point: [" << m_data.c->real() << ", " << m_data.c->imag() << " ]\n";
      break;
    case GC_STRING:
      stream << "String: " << m_data.s->c_str() << '\n';
      break;
    case GC_VEC_POINTER:
      stream << "Vector of generic pointer of size " << m_data.v_p->size() << '\n';
      break;
    case GC_VEC_BOOL:
      stream << "Vector of boolean of size " << m_data.v_b->size() << '\n';
      break;
    case GC_VEC_INTEGER:
      stream << "Vector of integer of size " << m_data.v_i->size() << '\n';
      break;
    case GC_VEC_LONG:
      stream << "Vector of long integer of size " << m_data.v_l->size() << '\n';
      break;
    case GC_VEC_REAL:
      stream << "Vector of floating point number of size " << m_data.v_r->size() << '\n';
      break;
    case GC_VEC_COMPLEX:
      stream << "Vector of complex floating point number of size " << m_data.v_c->size() << '\n';
      break;
    case GC_VEC_STRING:
      stream << "Vector of string of size " << m_data.v_s->size() << '\n';
      break;
    case GC_MAT_INTEGER:
      m_data.m_i->info(stream);
      break;
    case GC_MAT_LONG:
      m_data.m_l->info(stream);
      break;
    case GC_MAT_REAL:
      m_data.m_r->info(stream);
      break;
    case GC_MAT_COMPLEX:
      m_data.m_c->info(stream);
      break;
    case GC_VECTOR:
      stream << "Vector of generic data type of size " << m_data.v->size() << '\n';
      break;
    case GC_MAP:
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
  GenericContainer::operator [] ( unsigned i ) {
    switch ( ck( GC_VECTOR ) ) {
      case 0: break; // data present
      default: set_vector(); // data must be allocated;
    }
    CHECK_RESIZE(m_data.v,i);
    return (*m_data.v)[i];
  }

  GenericContainer const &
  GenericContainer::operator [] ( unsigned i ) const {
    GC_ASSERT(
      GC_VECTOR == m_data_type,
      "operator [] integer argument = " << i <<
      "\nexpect: " << typeName[GC_VECTOR] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
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
  GenericContainer::operator () ( unsigned i, char const * msg ) {
    ck(msg,GC_VECTOR);
    GC_ASSERT(
      i < m_data.v->size(),
      msg << "\noperator () const, index " << i << " out of range"
    )
    return (*m_data.v)[i];
  }

  GenericContainer const &
  GenericContainer::operator () ( unsigned i, char const * msg_in ) const {
    ck(msg_in,GC_VECTOR);
    char const * msg = msg_in == nullptr ? "" : msg_in;
    GC_ASSERT(
      i < m_data.v->size(),
      msg << "\noperator () const, index " << i << " out of range"
    )
    return (*m_data.v)[i];
  }

  GenericContainer &
  GenericContainer::operator () ( std::string const & s, char const * msg_in ) {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    GC_ASSERT(
      GC_MAP == m_data_type,
      "operator (), with string argument ``" << s.c_str() << "''"
      "\nexpect: " << typeName[GC_MAP] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    map_type::iterator iv = m_data.m->find(s);
    GC_ASSERT(
      iv != m_data.m->end(),
      msg << "\nGenericContainer::operator(): Cannot find key '" << s.c_str() << "'!"
    )
    return iv->second;
  }

  GenericContainer const &
  GenericContainer::operator () ( std::string const & s, char const * msg_in ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in;
    GC_ASSERT(
      GC_MAP == m_data_type,
      "operator () const, with string argument ``" << s.c_str() << "''"
      "\nexpect: " << typeName[GC_MAP] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    map_type::const_iterator iv = m_data.m->find(s);
    GC_ASSERT(
      iv != m_data.m->end(),
      msg << "\nGenericContainer::operator() const: Cannot find key '" << s.c_str() << "'!"
    )
    return iv->second;
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
  GenericContainer::operator [] ( std::string const & s ) {
    if ( ck( GC_MAP ) != 0 ) set_map(); // if not data present allocate!
    return (*m_data.m)[s];
  }

  GenericContainer const &
  GenericContainer::operator [] ( std::string const & s ) const {
    GC_ASSERT(
      GC_MAP == m_data_type,
      "operator [] string argument ``" << s.c_str() << "''"
      "\nexpect: " << typeName[GC_MAP] <<
      "\nbut data stored is of type: " << typeName[m_data_type]
    )
    return (*m_data.m)[s];
  }

  /*
  //   ____                            _
  //  |  _ \ _ __ ___  _ __ ___   ___ | |_ ___
  //  | |_) | '__/ _ \| '_ ` _ \ / _ \| __/ _ \
  //  |  __/| | | (_) | | | | | | (_) | ||  __/
  //  |_|   |_|  \___/|_| |_| |_|\___/ \__\___|
  */

  GenericContainer const &
  GenericContainer::promote_to_int() {
    switch (m_data_type) {
    case GC_NOTYPE:
      set_int(0);
      break;
    case GC_BOOL:
      set_int(m_data.b?1:0);
      break;
    case GC_INTEGER:
      break;
    case GC_POINTER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_int() cannot promote " << get_type_name() << " to int"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_long() {
    switch (m_data_type) {
    case GC_NOTYPE:
      set_long(0);
      break;
    case GC_BOOL:
      set_long(m_data.b?1:0);
      break;
    case GC_INTEGER:
      set_long(m_data.i);
      break;
    case GC_LONG:
      break;
    case GC_POINTER:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_LONG:
    case GC_MAT_INTEGER:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_long() cannot promote " << get_type_name() << " to long"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_real() {
    switch (m_data_type) {
    case GC_NOTYPE:
      set_real(0);
      break;
    case GC_BOOL:
      set_real(m_data.b?1:0);
      break;
    case GC_INTEGER:
      set_real(real_type(m_data.i));
      break;
    case GC_LONG:
      set_real(real_type(m_data.l));
      break;
    case GC_REAL:
      break;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_real() cannot promote " << get_type_name() << " to real"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_complex() {
    switch (m_data_type) {
    case GC_NOTYPE:
      set_complex(0,0);
      break;
    case GC_BOOL:
      set_complex(m_data.b?1:0,0);
      break;
    case GC_INTEGER:
      set_complex(real_type(m_data.i),0);
      break;
    case GC_LONG:
      set_complex(real_type(m_data.l),0);
      break;
    case GC_REAL:
      set_complex(real_type(m_data.r),0);
      break;
    case GC_COMPLEX:
      break;
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
      promote_to_vec_complex();
      break;
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_complex_type() cannot promote " << get_type_name() <<
        " to real type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_vec_int() {
    switch (m_data_type) {
    case GC_NOTYPE:
    { set_vec_int(1); get_int_at(0) = 0; }
      break;
    case GC_BOOL:
    { int_type tmp = m_data.b?1:0;
      set_vec_int(1);
      get_int_at(0) = tmp; }
      break;
    case GC_INTEGER:
    { int_type tmp = m_data.i;
      set_vec_int(1);
      get_int_at(0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_vec_int(unsigned(v_b->size()));
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.v_i)[i] = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      break;
    case GC_POINTER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_LONG:
    case GC_VEC_POINTER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_vec_int() cannot promote " << get_type_name() <<
        " to vec_int_type"
      )
    }
    return *this;

  }

  GenericContainer const &
  GenericContainer::promote_to_vec_long() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_vec_long(1); get_long_at(0) = 0; }
      break;
    case GC_BOOL:
      { long_type tmp = m_data.b?1:0;
        set_vec_long(1);
        get_long_at(0) = tmp; }
      break;
    case GC_INTEGER:
      { long_type tmp = long_type(m_data.i);
        set_vec_long(1);
        get_long_at(0) = tmp; }
      break;
    case GC_LONG:
      { long_type tmp = m_data.l;
        set_vec_long(1);
        get_long_at(0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_vec_long(unsigned(v_b->size()));
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.v_l)[i] = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_vec_long(unsigned(v_i->size()));
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.v_l)[i] = int_type( (*v_i)[i] );
        delete v_i;
      }
      break;
    case GC_VEC_LONG: // nothing to do
      break;
    case GC_POINTER:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_vec_long() cannot promote " << get_type_name() <<
        " to vec_long_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_vec_real() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_vec_real(1); get_real_at(0) = 0; }
      break;
    case GC_BOOL:
      { real_type tmp = m_data.b?1:0;
        set_vec_real(1);
        get_real_at(0) = tmp; }
      break;
    case GC_INTEGER:
      { real_type tmp = real_type( m_data.i );
        set_vec_real(1);
        get_real_at(0) = tmp; }
      break;
    case GC_LONG:
      { real_type tmp = real_type( m_data.l );
        set_vec_real(1);
        get_real_at(0) = tmp; }
      break;
    case GC_REAL:
      { real_type tmp = m_data.r;
        set_vec_real(1);
        get_real_at(0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b; // salva puntatore
        m_data_type = GC_NOTYPE;
        set_vec_real(unsigned(v_b->size()));
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.v_r)[i] = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i; // salva puntatore
        m_data_type = GC_NOTYPE;
        set_vec_real(unsigned(v_i->size()));
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.v_r)[i] = real_type( (*v_i)[i] );
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l; // salva puntatore
        m_data_type = GC_NOTYPE;
        set_vec_real(unsigned(v_l->size()));
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.v_r)[i] = real_type( (*v_l)[i] );
        delete v_l;
      }
      break;
    case GC_VEC_REAL:
      break;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_vec_real() cannot promote " << get_type_name() <<
        " vec_real_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_vec_complex() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_vec_complex(1); get_complex_at(0) = 0; }
      break;
    case GC_BOOL:
      { real_type tmp = m_data.b?1:0;
        set_vec_complex(1);
        get_complex_at(0) = tmp; }
      break;
    case GC_INTEGER:
      { real_type tmp = real_type( m_data.i );
        set_vec_complex(1);
        get_complex_at(0) = tmp; }
      break;
    case GC_LONG:
      { real_type tmp = real_type( m_data.l );
        set_vec_complex(1);
        get_complex_at(0) = tmp; }
      break;
    case GC_REAL:
      { real_type tmp = m_data.r;
        set_vec_complex(1);
        get_complex_at(0) = tmp; }
      break;
    case GC_COMPLEX:
      { complex_type tmp = *m_data.c;
        set_vec_complex(1);
        get_complex_at(0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_vec_complex(unsigned(v_b->size()));
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.v_c)[i] = complex_type( (*v_b)[i] ? 1: 0, 0 );
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_vec_complex(unsigned(v_i->size()));
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.v_c)[i] = complex_type( real_type( (*v_i)[i] ) , 0 );
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l;
        m_data_type = GC_NOTYPE;
        set_vec_complex(unsigned(v_l->size()));
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.v_c)[i] = complex_type( real_type( (*v_l)[i] ), 0 );
        delete v_l;
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type * v_r = m_data.v_r;
        m_data_type = GC_NOTYPE;
        set_vec_complex(unsigned(v_r->size()));
        for ( unsigned i = 0; i < v_r->size(); ++i )
          (*m_data.v_c)[i] = complex_type( (*v_r)[i], 0 );
        delete v_r;
      }
      break;
    case GC_VEC_COMPLEX:
      break;
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_vec_real() cannot promote " << get_type_name() <<
        " to vec_complex_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_mat_int() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_mat_int(1,1); get_int_at(0,0) = 0; }
      break;
    case GC_BOOL:
      { int_type tmp = m_data.b?1:0;
        set_mat_int(1,1);
        get_int_at(0,0) = tmp; }
      break;
    case GC_INTEGER:
      { int_type tmp = m_data.i;
        set_mat_int(1,1);
        get_int_at(0,0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_mat_int(unsigned(v_b->size()),1);
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.m_r)(i,0) = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_mat_int(unsigned(v_i->size()),1);
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.m_r)(i,0) = real_type( (*v_i)[i] );
        delete v_i;
      }
      break;
    case GC_MAT_INTEGER:
      break;
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_REAL:
    case GC_POINTER:
    case GC_STRING:
    case GC_LONG:
    case GC_COMPLEX:
    case GC_VEC_LONG:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_mat_int() cannot promote " << get_type_name() <<
        " to mat_int_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_mat_long() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_mat_long(1,1); get_long_at(0,0) = 0; }
      break;
    case GC_BOOL:
      { long_type tmp = m_data.b?1:0;
        set_mat_long(1,1);
        get_long_at(0,0) = tmp; }
      break;
    case GC_INTEGER:
      { long_type tmp = m_data.i;
        set_mat_long(1,1);
        get_long_at(0,0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_mat_real(unsigned(v_b->size()),1);
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.m_r)(i,0) = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_mat_long(unsigned(v_i->size()),1);
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.m_l)(i,0) = long_type((*v_i)[i]);
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l;
        m_data_type = GC_NOTYPE;
        set_mat_long(unsigned(v_l->size()),1);
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.m_l)(i,0) = (*v_l)[i];
        delete v_l;
      }
      break;
    case GC_MAT_INTEGER:
      { mat_int_type * m_i = m_data.m_i;
        m_data_type = GC_NOTYPE;
        set_mat_long(m_i->numRows(),m_i->numCols());
        for ( unsigned i = 0; i < m_i->size(); ++i )
          (*m_data.m_l)[i] = long_type((*m_i)[i]);
        delete m_i;
      }
      break;
    case GC_MAT_LONG:
      break;
    case GC_POINTER:
    case GC_STRING:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_mat_long() cannot promote " << get_type_name() <<
        " to mat_long_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_mat_real() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_mat_real(1,1); get_real_at(0,0) = 0; }
      break;
    case GC_BOOL:
      { real_type tmp = m_data.b?1:0;
        set_mat_real(1,1);
        get_real_at(0,0) = tmp; }
      break;
    case GC_INTEGER:
      { real_type tmp = m_data.i;
        set_mat_real(1,1);
        get_real_at(0,0) = tmp; }
      break;
    case GC_REAL:
      { real_type tmp = m_data.r;
        set_mat_real(1,1);
        get_real_at(0,0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_mat_real(unsigned(v_b->size()),1);
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.m_r)(i,0) = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_mat_real(unsigned(v_i->size()),1);
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.m_r)(i,0) = real_type( (*v_i)[i] );
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l;
        m_data_type = GC_NOTYPE;
        set_mat_real(unsigned(v_l->size()),1);
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.m_r)(i,0) = real_type( (*v_l)[i] );
        delete v_l;
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type * v_r = m_data.v_r;
        m_data_type = GC_NOTYPE;
        set_mat_real(unsigned(v_r->size()),1);
        for ( unsigned i = 0; i < v_r->size(); ++i )
          (*m_data.m_r)(i,0) = (*v_r)[i];
        delete v_r;
      }
      break;
    case GC_MAT_INTEGER:
      { mat_int_type * m_i = m_data.m_i;
        m_data_type = GC_NOTYPE;
        set_mat_real(m_i->numRows(),m_i->numCols());
        for ( unsigned i = 0; i < m_i->size(); ++i )
          (*m_data.m_r)[i] = real_type( (*m_i)[i] );
        delete m_i;
      }
      break;
    case GC_MAT_LONG:
      { mat_long_type * m_l = m_data.m_l;
        m_data_type = GC_NOTYPE;
        set_mat_real(m_l->numRows(),m_l->numCols());
        for ( unsigned i = 0; i < m_l->size(); ++i )
          (*m_data.m_r)[i] = real_type( (*m_l)[i] );
        delete m_l;
      }
      break;
    case GC_MAT_REAL:
      break;
    case GC_POINTER:
    case GC_STRING:
    case GC_LONG:
    case GC_COMPLEX:
    case GC_VEC_COMPLEX:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_mat_real() cannot promote " << get_type_name() <<
        " to mat_real_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_mat_complex() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_mat_complex(1,1); get_complex_at(0,0) = 0; }
      break;
    case GC_BOOL:
      { real_type tmp = m_data.b?1:0;
        set_mat_complex(1,1);
        get_complex_at(0,0) = tmp; }
      break;
    case GC_INTEGER:
      { real_type tmp = real_type( m_data.i );
        set_mat_complex(1,1);
        get_complex_at(0,0) = tmp; }
      break;
    case GC_LONG:
      { real_type tmp = real_type( m_data.l );
        set_mat_complex(1,1);
        get_complex_at(0,0) = tmp; }
      break;
    case GC_REAL:
      { real_type tmp = m_data.r;
        set_mat_complex(1,1);
        get_complex_at(0,0) = tmp; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_mat_complex(unsigned(v_b->size()),1);
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.m_r)(i,0) = ((*v_b)[i]?1:0);
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_mat_complex(unsigned(v_i->size()),1);
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.m_r)(i,0) = real_type( (*v_i)[i] );
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l;
       m_data_type = GC_NOTYPE;
        set_mat_complex(unsigned(v_l->size()),1);
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.m_r)(i,0) = real_type( (*v_l)[i] );
        delete v_l;
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type * v_r = m_data.v_r;
        m_data_type = GC_NOTYPE;
        set_mat_complex(unsigned(v_r->size()),1);
        for ( unsigned i = 0; i < v_r->size(); ++i )
          (*m_data.m_r)(i,0) = (*v_r)[i];
        delete v_r;
      }
      break;
    case GC_MAT_INTEGER:
      { mat_int_type * m_i = m_data.m_i;
        m_data_type = GC_NOTYPE;
        set_mat_complex(m_i->numRows(),m_i->numCols());
        for ( unsigned i = 0; i < m_i->size(); ++i )
          (*m_data.m_c)[i] = complex_type( real_type( (*m_i)[i] ), 0 );
        delete m_i;
      }
      break;
    case GC_MAT_LONG:
      { mat_long_type * m_l = m_data.m_l;
        m_data_type = GC_NOTYPE;
        set_mat_complex(m_l->numRows(),m_l->numCols());
        for ( unsigned i = 0; i < m_l->size(); ++i )
          (*m_data.m_c)[i] = complex_type( real_type( (*m_l)[i] ), 0 );
        delete m_l;
      }
      break;
    case GC_MAT_REAL:
      { mat_real_type * m_r = m_data.m_r;
        m_data_type = GC_NOTYPE;
        set_mat_complex(m_r->numRows(),m_r->numCols());
        for ( unsigned i = 0; i < m_r->size(); ++i )
          (*m_data.m_c)[i] = complex_type((*m_r)[i],0);
        delete m_r;
      }
      break;
    case GC_MAT_COMPLEX:
      break;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_mat_real() cannot promote " << get_type_name() <<
        " to mat_complex_type"
      )
    }
    return *this;
  }

  GenericContainer const &
  GenericContainer::promote_to_vector() {
    switch (m_data_type) {
    case GC_NOTYPE:
      { set_vector(1); (*this)[0].clear(); } // set data to no type
      break;
    case GC_POINTER:
      { set_vector(1); (*this)[0] = m_data.p; }
      break;
    case GC_BOOL:
      { set_vector(1); (*this)[0] = m_data.b; }
      break;
    case GC_INTEGER:
      { set_vector(1); (*this)[0] = m_data.i; }
      break;
    case GC_LONG:
      { set_vector(1); (*this)[0] = m_data.l; }
      break;
    case GC_REAL:
      { set_vector(1); (*this)[0] = m_data.r; }
      break;
    case GC_COMPLEX:
      { set_vector(1); (*this)[0] = *m_data.c; }
      break;
    case GC_STRING:
      { set_vector(1); (*this)[0] = *m_data.s; }
      break;
    case GC_VEC_POINTER:
      { vec_pointer_type * v_p = m_data.v_p;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_p->size()));
        for ( unsigned i = 0; i < v_p->size(); ++i )
          (*m_data.v)[i] = (*v_p)[i];
        delete v_p;
      }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = m_data.v_b;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_b->size()));
        for ( unsigned i = 0; i < v_b->size(); ++i )
          (*m_data.v)[i] = (*v_b)[i];
        delete v_b;
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = m_data.v_i;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_i->size()));
        for ( unsigned i = 0; i < v_i->size(); ++i )
          (*m_data.v)[i] = (*v_i)[i];
        delete v_i;
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type * v_l = m_data.v_l;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_l->size()));
        for ( unsigned i = 0; i < v_l->size(); ++i )
          (*m_data.v)[i] = (*v_l)[i];
        delete v_l;
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type * v_r = m_data.v_r;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_r->size()));
        for ( unsigned i = 0; i < v_r->size(); ++i )
          (*m_data.v)[i] = (*v_r)[i];
        delete v_r;
      }
      break;
    case GC_VEC_COMPLEX:
      { vec_complex_type * v_c = m_data.v_c;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_c->size()));
        for ( unsigned i = 0; i < v_c->size(); ++i )
          (*m_data.v)[i] = (*v_c)[i];
        delete v_c;
      }
      break;
    case GC_VEC_STRING:
      { vec_string_type * v_s = m_data.v_s;
        m_data_type = GC_NOTYPE;
        set_vector(unsigned(v_s->size()));
        for ( unsigned i = 0; i < v_s->size(); ++i )
          (*m_data.v)[i] = (*v_s)[i];
        delete v_s;
      }
      break;
    case GC_VECTOR:
      break;
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_MAP:
      GC_DO_ERROR(
        ":promote_to_vector() cannot promote " << get_type_name() <<
        " to vector_type"
      )
    }
    return *this;
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
    ostream_type      & stream,
    std::string const & prefix,
    std::string const & indent
  ) const {

    switch (m_data_type) {
    case GC_NOTYPE:
      stream << prefix.c_str() << "Empty!\n";
      break;
    case GC_POINTER:
      stream << prefix.c_str() << this->get_pvoid() << '\n';
      break;
    case GC_BOOL:
      stream << prefix.c_str() << (this->get_bool()?"true":"false") << '\n';
      break;
    case GC_INTEGER:
      stream << prefix.c_str() << this->get_int() << '\n';
      break;
    case GC_LONG:
      stream << prefix.c_str() << this->get_long() << '\n';
      break;
    case GC_REAL:
      stream << prefix.c_str() << this->get_real() << '\n';
      break;
    case GC_COMPLEX:
      stream << prefix.c_str() << "( " << this->get_complex().real()
             << ", " << this->get_complex().imag() << " )\n";
      break;
    case GC_STRING:
      stream << prefix.c_str() << "\"" << this->get_string().c_str() << "\"\n";
      break;
    case GC_VEC_POINTER:
      { vec_pointer_type const & v = this->get_vec_pointer();
        for ( vec_pointer_type::size_type i = 0; i < v.size(); ++i )
          stream << prefix.c_str() << "vec_pointer(" << i << "): " << v[i] << '\n';
      }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this->get_vec_bool();
        stream << prefix.c_str() << v << '\n'; }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type const & v = this->get_vec_int();
        stream << prefix.c_str() << v << '\n'; }
      break;
    case GC_VEC_LONG:
      { vec_long_type const & v = this->get_vec_long();
        stream << prefix.c_str() << v << '\n'; }
      break;
    case GC_VEC_REAL:
      { vec_real_type const & v = this->get_vec_real();
        stream << prefix.c_str() << v << '\n'; }
      break;
    case GC_VEC_COMPLEX:
      { vec_complex_type const & v = this->get_vec_complex();
        stream << prefix.c_str() << v << '\n'; }
      break;
    case GC_MAT_INTEGER:
      { mat_int_type const & m = this->get_mat_int();
        stream << m; }
      break;
    case GC_MAT_LONG:
      { mat_long_type const & m = this->get_mat_long();
        stream << m; }
      break;
    case GC_MAT_REAL:
      { mat_real_type const & m = this->get_mat_real();
        stream << m; }
      break;
    case GC_MAT_COMPLEX:
      { mat_complex_type const & m = this->get_mat_complex();
        stream << m; }
      break;
    case GC_VEC_STRING:
      { vec_string_type const & v = this->get_vec_string();
        stream << '\n';
        for ( vec_string_type::size_type i = 0; i < v.size(); ++i )
          stream << (prefix+indent).c_str() << i << ": \"" << v[i].c_str() << "\"\n";
      }
      break;

    case GC_VECTOR:
      { vector_type const & v = this->get_vector();
        for ( vector_type::size_type i = 0; i < v.size(); ++i ) {
          GenericContainer const & vi = v[i];
          if ( vi.simple_data() ||
               ( vi.simple_vec_data() && vi.get_num_elements() <= 10 ) ) {
            stream << prefix.c_str() << i << ": ";
            vi.dump(stream,"");
          } else {
            stream << prefix.c_str() << i << ":\n";
            vi.dump(stream,prefix+indent);
          }
        }
      }
      break;
    case GC_MAP:
      { map_type const & m = this->get_map();
        for ( map_type::const_iterator im = m.begin(); im != m.end(); ++im ) {
          #ifndef HAVE_WORKING_REGEX
          stream << prefix.c_str() << im->first << ":\n";
          im->second.dump(stream,prefix+indent);
          #else
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          std::string matches[4];
          int pcreExecRet = pcre_for_GC.exec( im->first, matches );
          if ( pcreExecRet == 4 ) {
            std::string header = matches[3]; // header
            // found formatting
            if ( im->second.simple_data() ) {
              stream << prefix.c_str() << header << ": ";
              im->second.dump(stream,"");
            } else {
              if ( matches[1].length() > 1 ) stream << '\n'; // double ## --> add nel line
              stream << prefix.c_str() << header.c_str();
              if ( matches[2].length() > 0 ) {
                stream << '\n' << prefix;
                char fmt = matches[2][0]; // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt; // underline header
              } else {
                stream << ':';
              }
              stream << '\n';
              im->second.dump(stream,prefix+indent);
            }
          } else {
            std::string header = pcreExecRet == 3 ? matches[3] : im->first;
            if ( im->second.simple_data() ) {
              stream << prefix.c_str() << header.c_str() << ": ";
              im->second.dump(stream,"");
            } else {
              stream << prefix.c_str() << header.c_str() << ":\n";
              im->second.dump(stream,prefix+indent);
            }
          }
          #endif
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
    ostream_type      & stream,
    std::string const & prefix,
    std::string const & indent
  ) const {

    switch (m_data_type) {
    case GC_NOTYPE:
      stream << prefix.c_str() << "Empty!\n";
      break;
    case GC_POINTER:
      stream << prefix.c_str() << "(*void)\n";
      break;
    case GC_BOOL:
      stream << prefix.c_str() << "bool\n";
      break;
    case GC_INTEGER:
      stream << prefix.c_str() << "int\n";
      break;
    case GC_LONG:
      stream << prefix.c_str() << "long int\n";
      break;
    case GC_REAL:
      stream << prefix.c_str() << "double\n";
      break;
    case GC_COMPLEX:
      stream << prefix.c_str() << "complex\n";
      break;
    case GC_STRING:
      stream << prefix.c_str() << "string\n";
      break;
    case GC_VEC_POINTER:
      { vec_pointer_type const & v = this->get_vec_pointer();
        stream << "vector of pointer[" << v.size() << "]\n"; }
      break;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this->get_vec_bool();
        stream << "vector of bool[" << v.size() << "]\n"; }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type const & v = this->get_vec_int();
        stream << "vector of int[" << v.size() << "]\n"; }
      break;
    case GC_VEC_LONG:
      { vec_long_type const & v = this->get_vec_long();
        stream << "vector of long[" << v.size() << "]\n"; }
      break;
    case GC_VEC_REAL:
      { vec_real_type const & v = this->get_vec_real();
        stream << "vector of double[" << v.size() << "]\n"; }
      break;
    case GC_VEC_COMPLEX:
      { vec_complex_type const & v = this->get_vec_complex();
        stream << "vector of complex[" << v.size() << "]\n"; }
      break;
    case GC_MAT_INTEGER:
      { mat_int_type const & m = this->get_mat_int();
        stream << "matrix of int[" << m.numRows() << "," << m.numCols() << "]\n"; }
      break;
    case GC_MAT_LONG:
      { mat_long_type const & m = this->get_mat_long();
        stream << "matrix of long[" << m.numRows() << "," << m.numCols() << "]\n"; }
      break;
    case GC_MAT_REAL:
      { mat_real_type const & m = this->get_mat_real();
        stream << "matrix of double[" << m.numRows() << "," << m.numCols() << "]\n"; }
      break;
    case GC_MAT_COMPLEX:
      { mat_complex_type const & m = this->get_mat_complex();
        stream << "matrix of complex[" << m.numRows() << "," << m.numCols() << "]\n"; }
      break;
    case GC_VEC_STRING:
      { vec_string_type const & v = this->get_vec_string();
        stream << "vector of string[" << v.size() << "]\n"; }
      break;

    case GC_VECTOR:
      { vector_type const & v = this->get_vector();
        for ( vector_type::size_type i = 0; i < v.size(); ++i ) {
          GenericContainer const & vi = v[i];
          if ( vi.simple_data() || vi.simple_vec_data()) {
            stream << prefix.c_str() << i << ": ";
            vi.print_content_types(stream,"");
          } else {
            stream << prefix.c_str() << i << ":\n";
            vi.print_content_types(stream,prefix+indent,indent);
          }
        }
      }
      break;
    case GC_MAP:
      { map_type const & m = this->get_map();
        for ( map_type::const_iterator im = m.begin(); im != m.end(); ++im ) {
          #ifndef HAVE_WORKING_REGEX
          stream << prefix.c_str() << im->first << ":\n";
          im->second.print_content_types(stream,prefix+indent,indent);
          #else
          // check formatting using pcre
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          std::string matches[4];
          int pcreExecRet = pcre_for_GC.exec( im->first, matches );
          if ( pcreExecRet == 4 ) {
            std::string header = matches[3]; // header
            // found formatting
            if ( im->second.simple_data() || im->second.simple_vec_data() ) {
              stream << prefix.c_str() << header << ": ";
              im->second.print_content_types(stream,"");
            } else {
              if ( matches[1].length() > 1 ) stream << '\n'; // double ## --> add nel line
              stream << prefix.c_str() << header.c_str();
              if ( matches[2].length() > 0 ) {
                stream << '\n' << prefix;
                char fmt = matches[2][0]; // underline char
                std::size_t m3 = header.length();
                while ( m3-- > 0 ) stream << fmt; // underline header
              } else {
                stream << ':';
              }
              stream << '\n';
              im->second.print_content_types(stream,prefix+indent,indent);
            }
          } else {
            std::string header = pcreExecRet == 3 ? matches[3] : im->first;
            if ( im->second.simple_data() || im->second.simple_vec_data() ) {
              stream << prefix.c_str() << header.c_str() << ": ";
              im->second.print_content_types(stream,"");
            } else {
              stream << prefix.c_str() << header.c_str() << ":\n";
              im->second.print_content_types(stream,prefix+indent,indent);
            }
          }
          #endif
        }
      }
      break;

    //default:
    //  GC_DO_ERROR( "Error, print(...) unknown type!\n");
    //  break;
    }
  }

  /*
  //   _                                 _
  //  | |_ ___     _   _  __ _ _ __ ___ | |
  //  | __/ _ \   | | | |/ _` | '_ ` _ \| |
  //  | || (_) |  | |_| | (_| | | | | | | |
  //   \__\___/____\__, |\__,_|_| |_| |_|_|
  //         |_____|___/
  */

  void
  GenericContainer::to_yaml(
    ostream_type      & stream,
    std::string const & prefix
  ) const {
    switch (m_data_type) {
    case GC_NOTYPE:
      stream << "Empty!\n";
      break;
    case GC_BOOL:
      stream << (this->get_bool()?"true":"false") << '\n';
      break;
    case GC_INTEGER:
      stream << this->get_int() << '\n';
      break;
    case GC_LONG:
      stream << this->get_long() << '\n';
      break;
    case GC_REAL:
      stream << this->get_real() << '\n';
      break;
    case GC_COMPLEX:
      stream << this->get_complex().real() << ' '
             << this->get_complex().imag() << '\n';
      break;
    case GC_STRING:
      stream << "'" << this->get_string().c_str() << "'\n";
      break;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this->get_vec_bool();
        stream << "[ " << (v[0]?"true":"false");
        for ( vec_bool_type::size_type i = 1; i < v.size(); ++i )
          stream << ", " << (v[i]?"true":"false");
        stream << " ]\n";
      }
      break;
    case GC_VEC_INTEGER:
      { vec_int_type const & v = this->get_vec_int();
        stream << "[ " << v[0];
        for ( vec_int_type::size_type i = 1; i < v.size(); ++i )
          stream << ", " << v[i];
        stream << " ]\n";
      }
      break;
    case GC_VEC_LONG:
      { vec_long_type const & v = this->get_vec_long();
        stream << "[ " << v[0];
        for ( vec_long_type::size_type i = 1; i < v.size(); ++i )
          stream << ", " << v[i];
        stream << " ]\n";
      }
      break;
    case GC_VEC_REAL:
      { vec_real_type const & v = this->get_vec_real();
        stream << "[ " << v[0];
        for ( vec_real_type::size_type i = 1; i < v.size(); ++i )
          stream << ", " << v[i];
        stream << " ]\n";
      }
      break;
    case GC_VEC_STRING:
      { vec_string_type const & v = this->get_vec_string();
        stream << "[ '" << v[0].c_str() << "'";
        for ( vec_string_type::size_type i = 1; i < v.size(); ++i )
          stream << ", '" << v[i].c_str() << "'";
        stream << " ]\n";
      }
      break;

    case GC_VECTOR:
      { vector_type const & v = this->get_vector();
        stream << '\n';
        for ( vector_type::size_type i = 0; i < v.size(); ++i ) {
          stream << prefix.c_str() << "- ";
          v[i].to_yaml(stream,prefix+"  ");
        }
      }
      break;
    case GC_MAP:
      { map_type const & m = this->get_map();
        stream << '\n';
        for ( map_type::const_iterator im = m.begin(); im != m.end(); ++im ) {
          stream << prefix.c_str() << im->first.c_str() << ": ";
          im->second.to_yaml(stream,prefix+"  ");
        }
      }
      break;
    case GC_MAT_INTEGER:
    case GC_MAT_LONG:
    case GC_MAT_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
    case GC_POINTER:
    case GC_VEC_POINTER:
      { /* DA FARE */ }
      break;
    }
  }

  GenericContainer &
  GenericContainer::readFormattedData(
    char const * fname,
    char const * commentChars,
    char const * delimiters
  ) {
    std::ifstream file( fname );
    GC_ASSERT(
      file.good(),
      "readFormattedData, failed to open file: ``" << fname << "''"
    )
    return readFormattedData( file, commentChars, delimiters );
  }

  void
  GenericContainer::exception( char const * msg ) {
    throw std::runtime_error(msg);
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
