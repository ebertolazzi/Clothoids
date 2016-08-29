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

// use pcre for pattern matching
#ifdef GENERIC_CONTAINER_USE_CXX11
  #include <regex>
#endif

#define CHECK_RESIZE(pV,I) if ( pV->size() <= (I) ) pV->resize((I)+1)

namespace GenericContainerNamespace {

  using std::fpclassify;

  static
  inline
  bool isZero( real_type x )
  { return FP_ZERO == fpclassify(x) ; }

  static
  inline
  bool isInteger( real_type x )
  { return isZero(x-floor(x)) ; }

  static
  inline
  bool isUnsigned( real_type x )
  { return isInteger(x) && x >= 0 ; }

  std::ostream &
  operator << ( std::ostream & s, vec_bool_type const & v ) {
    if ( v.size() > 0 ) {
      s << (v[0]?"[ true":"[ false") ;
      for ( unsigned i = 1 ; i < v.size() ; ++i )
        s << (v[i]?" true":" false") ;
      s << " ]" ;
    } else {
      s << "[]" ;
    }
    return s ;
  }

  std::ostream &
  operator << ( std::ostream & s, vec_int_type const & v ) {
    if ( v.size() > 0 ) {
      s << "[ " << v[0] ;
      for ( unsigned i = 1 ; i < v.size() ; ++i ) s << " " << v[i] ;
      s << " ]" ;
    } else {
      s << "[]" ;
    }
    return s ;
  }

  std::ostream &
  operator << ( std::ostream & s, vec_real_type const & v ) {
    if ( v.size() > 0 ) {
      s << "[ " << v[0] ;
      for ( unsigned i = 1 ; i < v.size() ; ++i ) s << " " << v[i] ;
      s << " ]" ;
    } else {
      s << "[]" ;
    }
    return s ;
  }

  std::ostream &
  operator << ( std::ostream & s, vec_complex_type const & v ) {
    if ( v.size() > 0 ) {
      s << "[ (" << v[0].real() << ", " << v[0].imag() << " )"  ;
      for ( unsigned i = 1 ; i < v.size() ; ++i )
        s << " (" << v[i].real() << ", " << v[i].imag() << " )"  ;
      s << " ]" ;
    } else {
      s << "[]" ;
    }
    return s ;
  }

  real_type const &
  mat_real_type::operator () ( unsigned i, unsigned j ) const {
    GC_ASSERT( i < _numRows && j < _numCols,
               "mat_real_type::operator() (" << i << ", " << j <<
               ") index out of range: [0," << _numRows <<
               ") x [0," << _numCols << ")\n" ) ;
    return (*this)[std::size_t(i+j*_numRows)] ;
  }

  real_type &
  mat_real_type::operator () ( unsigned i, unsigned j ) {
    GC_ASSERT( i < _numRows && j < _numCols,
               "mat_real_type::operator() (" << i << ", " << j <<
               ") index out of range: [0," << _numRows <<
               ") x [0," << _numCols << ")\n" ) ;
    return (*this)[std::size_t(i+j*_numRows)] ;
  }

  void
  mat_real_type::getColumn( unsigned nc, vec_real_type & C ) const {
    GC_ASSERT( nc < _numCols,
               "mat_real_type::getColumn(" << nc <<
               ",C) column index out of range max = " << _numCols-1 ) ;
    C.clear() ;
    C.reserve(_numRows) ;
    for ( unsigned i = 0 ; i < _numRows ; ++i )
      C.push_back( (*this)(i,nc) ) ;
  }

  void
  mat_real_type::getRow( unsigned nr, vec_real_type & R ) const {
    GC_ASSERT( nr < _numRows,
               "mat_real_type::getRow(" << nr <<
               ",C) row index out of range max = " << _numRows-1 ) ;
    R.clear() ;
    R.reserve(_numCols) ;
    for ( unsigned j = 0 ; j < _numCols ; ++j )
      R.push_back( (*this)(nr,j) ) ;
  }

  void
  mat_real_type::info( std::basic_ostream<char> & stream ) const {
    stream << "Matrix of floating point number of size "
           << _numRows << " x " << _numCols
           << '\n' ;
  }

  complex_type const &
  mat_complex_type::operator () ( unsigned i, unsigned j ) const {
    GC_ASSERT( i < _numRows && j < _numCols,
               "mat_real_type::operator() (" << i << ", " << j <<
               ") index out of range: [0," << _numRows <<
               ") x [0," << _numCols << ")\n" ) ;
    return (*this)[std::size_t(i+j*_numRows)] ;
  }

  complex_type &
  mat_complex_type::operator () ( unsigned i, unsigned j ) {
    GC_ASSERT( i < _numRows && j < _numCols,
               "mat_real_type::operator() (" << i << ", " << j <<
               ") index out of range: [0," << _numRows <<
               ") x [0," << _numCols << ")\n" ) ;
    return (*this)[std::size_t(i+j*_numRows)] ;
  }

  void
  mat_complex_type::getColumn( unsigned nc, vec_complex_type & C ) const {
    GC_ASSERT( nc < _numCols,
               "mat_real_type::getColumn(" << nc <<
               ",C) column index out of range max = " << _numCols-1 ) ;
    C.clear() ;
    C.reserve(_numRows) ;
    for ( unsigned i = 0 ; i < _numRows ; ++i )
      C.push_back( (*this)(i,nc) ) ;
  }

  void
  mat_complex_type::getRow( unsigned nr, vec_complex_type & R ) const {
    GC_ASSERT( nr < _numRows,
               "mat_real_type::getRow(" << nr <<
               ",C) row index out of range max = " << _numRows-1 ) ;
    R.clear() ;
    R.reserve(_numCols) ;
    for ( unsigned j = 0 ; j < _numCols ; ++j )
      R.push_back( (*this)(nr,j) ) ;
  }

  void
  mat_complex_type::info( std::basic_ostream<char> & stream ) const {
    stream << "Matrix of complex floating point number of size "
           << _numRows << " x " << _numCols
           << '\n' ;
  }

  std::ostream &
  operator << ( std::ostream & s, mat_real_type const & m ) {
    for ( unsigned i = 0 ; i < m.numRows() ; ++i ) {
      s << std::setw(8) << m(i,0) ;
      for ( unsigned j = 1 ; j < m.numCols() ; ++j )
        s << " " << std::setw(8) << m(i,j) ;
      s << '\n' ;
    }
    return s ;
  }

  std::ostream &
  operator << ( std::ostream & s, mat_complex_type const & m ) {
    for ( unsigned i = 0 ; i < m.numRows() ; ++i ) {
      s << std::setw(8) << m(i,0) ;
      for ( unsigned j = 1 ; j < m.numCols() ; ++j )
        s << " (" << std::setw(8) << m(i,j).real() << ", " << std::setw(8) << m(i,j).imag() << " )" ;
      s << '\n' ;
    }
    return s ;
  }

  #ifdef GENERIC_CONTAINER_USE_CXX11

  class GENERIC_CONTAINER_API_DLL Pcre_for_GC {

  private:
  
    std::regex  reCompiled ;
    std::smatch reMatches ;

  public:

    Pcre_for_GC()
    : reCompiled("^\\s*\\d+\\s*(##?)(-|=|~|_|)\\s*(.*)$")
    { }

    int
    exec( std::string const & str, std::string matches[4] ) {
      if ( std::regex_match( str, reMatches, reCompiled ) ) {
        for ( std::size_t i = 0 ; i < reMatches.size() ; ++i )
          matches[i] = reMatches[i].str() ;
        return int(reMatches.size()) ;
      } else {
        return 0 ;
      }
    }
  } ;

  static Pcre_for_GC pcre_for_GC ;

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
    "vec_real_type",
    "vec_complex_type",
    "vec_string_type",

    "mat_real_type",
    "mat_complex_type",

    "vector_type",
    "map_type"
  } ;

  // costruttore
  GenericContainer::GenericContainer()
  : _data_type(GC_NOTYPE)
  {
  }

  #ifdef GENERIC_CONTAINER_ON_WINDOWS
  bool
  GenericContainer::simple_data() const {
    return _data_type <= GC_STRING ;
  }
  bool
  GenericContainer::simple_vec_data() const {
    return _data_type <= GC_VEC_STRING ;
  }
  #endif

  // distruttore
  void
  GenericContainer::clear() {
    switch (_data_type) {
    case GC_NOTYPE:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_POINTER:
      // removed annoying warning. To be re-thinked...
      //GC_WARNING( _data.p == nullptr, "find a pointer not deallocated!" ) ;
      break ;
    case GC_STRING:      delete _data.s   ; break ;
    case GC_COMPLEX:     delete _data.c   ; break ;

    case GC_VEC_POINTER: delete _data.v_p ; break ;
    case GC_VEC_BOOL:    delete _data.v_b ; break ;
    case GC_VEC_INTEGER: delete _data.v_i ; break ;
    case GC_VEC_REAL:    delete _data.v_r ; break ;
    case GC_VEC_COMPLEX: delete _data.v_c ; break ;
    case GC_MAT_REAL:    delete _data.m_r ; break ;
    case GC_MAT_COMPLEX: delete _data.m_c ; break ;
    case GC_VEC_STRING:  delete _data.v_s ; break ;

    case GC_VECTOR:      delete _data.v   ; break ;
    case GC_MAP:         delete _data.m   ; break ;
    }
    _data_type = GC_NOTYPE ;
  }

  //! Return a string representing the type of data stored
  char const *
  GenericContainer::get_type_name() const {
    return typeName[_data_type] ;
  }

  unsigned
  GenericContainer::get_num_elements() const {
    switch (_data_type) {
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:      return 1 ;
    case GC_VEC_POINTER: return unsigned(_data.v_p->size()) ;
    case GC_VEC_BOOL:    return unsigned(_data.v_b->size()) ;
    case GC_VEC_INTEGER: return unsigned(_data.v_i->size()) ;
    case GC_VEC_REAL:    return unsigned(_data.v_r->size()) ;
    case GC_VEC_COMPLEX: return unsigned(_data.v_c->size()) ;
    case GC_MAT_REAL:    return unsigned(_data.m_r->size()) ;
    case GC_MAT_COMPLEX: return unsigned(_data.m_c->size()) ;
    case GC_VEC_STRING:  return unsigned(_data.v_s->size()) ;
    case GC_VECTOR:      return unsigned(_data.v->size()) ;
    case GC_MAP:         return unsigned(_data.m->size()) ;
    case GC_NOTYPE:      return 0 ;
    }
    return 0 ;
  }

  unsigned
  GenericContainer::get_numRows() const {
    switch (_data_type) {
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
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:      return 1 ;
    case GC_MAT_REAL:    return unsigned(_data.m_r->numRows()) ;
    case GC_MAT_COMPLEX: return unsigned(_data.m_c->numRows()) ;
    case GC_MAP:         return 1 ;
    case GC_NOTYPE:      return 0 ;
    }
    return 0 ;
  }

  unsigned
  GenericContainer::get_numCols() const {
    switch (_data_type) {
    case GC_POINTER:
    case GC_BOOL:
    case GC_INTEGER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:      return 1 ;
    case GC_VEC_POINTER: return unsigned(_data.v_p->size()) ;
    case GC_VEC_BOOL:    return unsigned(_data.v_b->size()) ;
    case GC_VEC_INTEGER: return unsigned(_data.v_i->size()) ;
    case GC_VEC_REAL:    return unsigned(_data.v_r->size()) ;
    case GC_VEC_COMPLEX: return unsigned(_data.v_c->size()) ;
    case GC_MAT_REAL:    return unsigned(_data.m_r->numCols()) ;
    case GC_MAT_COMPLEX: return unsigned(_data.m_c->numCols()) ;
    case GC_VEC_STRING:  return unsigned(_data.v_s->size()) ;
    case GC_VECTOR:      return unsigned(_data.v->size()) ;
    case GC_MAP:         return unsigned(_data.m->size()) ;
    case GC_NOTYPE:      return 0 ;
    }
    return 0 ;
  }

  //! Assign a generic container `a` to the generic container.
  void
  GenericContainer::load( GenericContainer const & gc ) {
    this -> clear() ;
    switch (gc._data_type) {
    case GC_NOTYPE:      break ;
    case GC_POINTER:     this -> set_pointer(gc._data.p)  ; break ;
    case GC_BOOL:        this -> set_bool(gc._data.b)     ; break ;
    case GC_INTEGER:     this -> set_int(gc._data.i)      ; break ;
    case GC_LONG:        this -> set_long(gc._data.l)      ; break ;
    case GC_REAL:        this -> set_real(gc._data.r)     ; break ;
    case GC_COMPLEX:     this -> set_complex(*gc._data.c) ; break ;
    case GC_STRING:      this -> set_string(*gc._data.s)  ; break ;

    case GC_VEC_POINTER: this -> set_vec_pointer(*gc._data.v_p) ; break ;
    case GC_VEC_BOOL:    this -> set_vec_bool(*gc._data.v_b)    ; break ;
    case GC_VEC_INTEGER: this -> set_vec_int(*gc._data.v_i)     ; break ;
    case GC_VEC_REAL:    this -> set_vec_real(*gc._data.v_r)    ; break ;
    case GC_VEC_COMPLEX: this -> set_vec_complex(*gc._data.v_c) ; break ;
    case GC_VEC_STRING:  this -> set_vec_string(*gc._data.v_s)  ; break ;

    case GC_MAT_REAL:    this -> set_mat_real(*gc._data.m_r)    ; break ;
    case GC_MAT_COMPLEX: this -> set_mat_complex(*gc._data.m_c) ; break ;

    case GC_VECTOR:
      { unsigned N = unsigned(gc._data.v->size()) ;
        allocate_vector( N ) ;
        std::copy( gc._data.v->begin(),
                   gc._data.v->end(),
                   this->_data.v->begin() ) ;
      }
      break ;
    case GC_MAP:
      { allocate_map() ;
        // this->_data.m->insert( gc._data.m->begin(), gc._data.m->end() ) ; !!!! DO NOT WORK ON CLANG
        for ( map_type::iterator it = gc._data.m->begin() ;
              it != gc._data.m->end() ; ++it )
          (*this->_data.m)[it->first] = it->second ;
      }
      break ;
    }
  }

  int
  GenericContainer::ck( TypeAllowed tp ) const {
    if ( tp == _data_type ) return 0 ; // ok
    if ( tp == GC_NOTYPE  ) return 1 ; //
    return 2 ;
  }

  void
  GenericContainer::ck( char const msg_in[], TypeAllowed tp ) const {
    char const * msg = msg_in == nullptr ? "" : msg_in ;
    GC_ASSERT( tp == _data_type,
               msg <<
               "bad data type, expect: " << typeName[tp] <<
               " but data stored is of type: " << typeName[_data_type] ) ;
  }

  void
  GenericContainer::ck_or_set( char const msg[], TypeAllowed tp) {
    if ( _data_type == GC_NOTYPE ) _data_type = tp ;
    else                           ck(msg,tp) ;
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
    if ( _data_type != GC_STRING ) {
      clear() ;
      _data_type = GC_STRING ;
      _data.s    = new string_type ;
    }
  }

  void
  GenericContainer::allocate_complex() {
    if ( _data_type != GC_COMPLEX ) {
      clear() ;
      _data_type = GC_COMPLEX ;
      _data.c    = new complex_type ;
    }
  }

  void
  GenericContainer::allocate_vec_pointer( unsigned sz ) {
    if ( _data_type != GC_VEC_POINTER ) {
      clear() ;
      _data_type = GC_VEC_POINTER ;
      _data.v_p  = new vec_pointer_type() ;
    }
    if ( sz > 0 ) _data.v_p -> resize( sz ) ;
  }

  GenericContainer &
  GenericContainer::free_pointer() {
    GC_ASSERT( GC_POINTER == _data_type || GC_NOTYPE == _data_type,
               "free_pointer() bad data type\nexpect: " << typeName[GC_POINTER] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    _data.p = nullptr ;
    _data_type = GC_NOTYPE ;
    return *this ;
  }

  void
  GenericContainer::allocate_vec_bool( unsigned sz ) {
    if ( _data_type != GC_VEC_BOOL ) {
      clear() ;
      _data_type = GC_VEC_BOOL ;
      _data.v_b  = new vec_bool_type() ;
    }
    if ( sz > 0 ) _data.v_b -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_vec_int( unsigned sz ) {
    if ( _data_type != GC_VEC_INTEGER ) {
      clear() ;
      _data_type = GC_VEC_INTEGER ;
      _data.v_i  = new vec_int_type() ;
    }
    if ( sz > 0 ) _data.v_i -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_vec_real( unsigned sz ) {
    if ( _data_type != GC_VEC_REAL ) {
      clear() ;
      _data_type = GC_VEC_REAL ;
      _data.v_r  = new vec_real_type() ;
    }
    if ( sz > 0 ) _data.v_r -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_vec_complex( unsigned sz ) {
    if ( _data_type != GC_VEC_COMPLEX ) {
      clear() ;
      _data_type = GC_VEC_COMPLEX ;
      _data.v_c  = new vec_complex_type() ;
    }
    if ( sz > 0 ) _data.v_c -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_mat_real( unsigned nr, unsigned nc ) {
    if ( _data_type != GC_MAT_REAL ) {
      clear() ;
      _data_type = GC_MAT_REAL ;
      _data.m_r  = new mat_real_type( nr, nc ) ;
    } else {
      _data.m_r -> resize( nr, nc ) ;
    }
  }

  void
  GenericContainer::allocate_mat_complex( unsigned nr, unsigned nc ) {
    if ( _data_type != GC_MAT_COMPLEX ) {
      clear() ;
      _data_type = GC_MAT_COMPLEX ;
      _data.m_c  = new mat_complex_type( nr, nc ) ;
    } else {
      _data.m_c -> resize( nr, nc ) ;
    }
  }

  void
  GenericContainer::allocate_vec_string( unsigned sz ) {
    if ( _data_type != GC_VEC_STRING ) {
      clear() ;
      _data_type = GC_VEC_STRING ;
      _data.v_s  = new vec_string_type() ;
    }
    if ( sz > 0 ) _data.v_s -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_vector( unsigned sz ) {
    if ( _data_type != GC_VECTOR ) {
      clear() ;
      _data_type = GC_VECTOR ;
      _data.v    = new vector_type() ;
    }
    if ( sz > 0 ) _data.v -> resize( sz ) ;
  }

  void
  GenericContainer::allocate_map() {
    if ( _data_type != GC_MAP ) {
      clear() ;
      _data_type = GC_MAP ;
      _data.m    = new map_type() ;
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
    clear() ;
    _data_type = GC_POINTER ;
    return (_data.p = value) ;
  }

  bool_type &
  GenericContainer::set_bool( bool_type value ) {
    clear() ;
    _data_type = GC_BOOL ;
    return (_data.b = value) ;
  }

  int_type &
  GenericContainer::set_int( int_type value ) {
    clear() ;
    _data_type = GC_INTEGER ;
    return (_data.i = value) ;
  }

  long_type &
  GenericContainer::set_long( long_type value ) {
    clear() ;
    _data_type = GC_LONG ;
    return (_data.l = value) ;
  }

  real_type &
  GenericContainer::set_real( real_type value ) {
    clear() ;
    _data_type = GC_REAL ;
    return (_data.r = value) ;
  }

  complex_type &
  GenericContainer::set_complex( complex_type & value ) {
    clear() ;
    _data_type = GC_COMPLEX ;
    _data.c    = new complex_type ;
    return (*_data.c=value) ;
  }

  complex_type &
  GenericContainer::set_complex( real_type re, real_type im ) {
    clear() ;
    _data_type = GC_COMPLEX ;
    _data.c    = new complex_type(re,im) ;
    return *_data.c ;
  }

  string_type &
  GenericContainer::set_string( string_type const & value ) {
    allocate_string() ;
    return (*_data.s = value) ;
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( unsigned sz ) {
    allocate_vec_pointer( sz ) ;
    return *_data.v_p ;
  }

  vec_pointer_type &
  GenericContainer::set_vec_pointer( vec_pointer_type const & v ) {
    allocate_vec_pointer( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_p->begin() ) ;
    return *_data.v_p ;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( unsigned sz ) {
    allocate_vec_bool( sz ) ; return *_data.v_b ;
  }

  vec_bool_type &
  GenericContainer::set_vec_bool( vec_bool_type const & v ) {
    allocate_vec_bool( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_b->begin() ) ;
    return *_data.v_b ;
  }

  vec_int_type &
  GenericContainer::set_vec_int( unsigned sz ) {
    allocate_vec_int( sz ) ;
    return *_data.v_i ;
  }

  vec_int_type &
  GenericContainer::set_vec_int( vec_int_type const & v ) {
    allocate_vec_int( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_i->begin() ) ;
    return *_data.v_i ;
  }

  vec_real_type &
  GenericContainer::set_vec_real( unsigned sz ) {
    allocate_vec_real( sz ) ;
    return *_data.v_r ;
  }

  vec_real_type &
  GenericContainer::set_vec_real( vec_real_type const & v ) {
    allocate_vec_real( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_r->begin() ) ;
    return *_data.v_r ;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( unsigned sz ) {
    allocate_vec_complex( sz ) ;
    return *_data.v_c ;
  }

  vec_complex_type &
  GenericContainer::set_vec_complex( vec_complex_type const & v ) {
    allocate_vec_complex( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_c->begin() ) ;
    return *_data.v_c ;
  }

  mat_real_type &
  GenericContainer::set_mat_real( unsigned nr, unsigned nc ) {
    allocate_mat_real( nr, nc ) ;
    return *_data.m_r ;
  }

  mat_real_type &
  GenericContainer::set_mat_real( mat_real_type const & m ) {
    allocate_mat_real( m.numRows(), m.numCols() ) ;
    std::copy( m.begin(), m.end(), _data.m_r->begin() ) ;
    return *_data.m_r ;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( unsigned nr, unsigned nc ) {
    allocate_mat_complex( nr, nc ) ;
    return *_data.m_c ;
  }

  mat_complex_type &
  GenericContainer::set_mat_complex( mat_complex_type const & m ) {
    allocate_mat_complex( m.numRows(), m.numCols() ) ;
    std::copy( m.begin(), m.end(), _data.m_c->begin() ) ;
    return *_data.m_c ;
  }

  vec_string_type &
  GenericContainer::set_vec_string( unsigned sz ) {
    allocate_vec_string( sz ) ;
    return *_data.v_s ;
  }

  vec_string_type &
  GenericContainer::set_vec_string( vec_string_type const & v ) {
    allocate_vec_string( unsigned(v.size()) ) ;
    std::copy( v.begin(), v.end(), _data.v_s->begin() ) ;
    return *_data.v_s ;
  }

  vector_type &
  GenericContainer::set_vector( unsigned sz ) {
    allocate_vector( sz ) ;
    return *_data.v ;
  }

  map_type &
  GenericContainer::set_map() {
    allocate_map() ;
    return *_data.m ;
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
    if ( _data_type == GC_VEC_BOOL ) {
      _data.v_b->push_back( val ) ;
    } else if ( _data_type == GC_VEC_INTEGER ) {
      _data.v_i->push_back( val ? 1 : 0 ) ;
    } else if ( _data_type == GC_VEC_REAL ) {
      _data.v_r->push_back( val ? 1 : 0 ) ;
    } else if ( _data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val ? 1 : 0, 0 ) ;
      _data.v_c->push_back( tmp ) ;
    } else if ( _data_type == GC_VECTOR   ) {
      _data.v->resize(_data.v->size()+1) ;
      _data.v->back().set_bool( val ) ;
    } else {
      GC_ASSERT( false, "push_bool, bad data stored: " << get_type_name() ) ;
    }
  }

  void
  GenericContainer::push_int( int_type val ) {
    if ( _data_type == GC_VEC_INTEGER ) {
      _data.v_i->push_back( val ) ;
    } else if ( _data_type == GC_VEC_REAL ) {
      _data.v_r->push_back( val ) ;
    } else if ( _data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val, 0 ) ;
      _data.v_c->push_back( tmp ) ;
    } else if ( _data_type == GC_VECTOR ) {
      _data.v->resize(_data.v->size()+1) ;
      _data.v->back().set_int( val ) ;
    } else {
      if ( _data_type != GC_VEC_INTEGER ) promote_to_vec_int() ;
      _data.v_i->push_back( val ) ;
    }
  }

  void
  GenericContainer::push_real( real_type val ) {
    if ( _data_type == GC_VEC_REAL ) {
      _data.v_r->push_back( val ) ;
    } else if ( _data_type == GC_VEC_COMPLEX ) {
      complex_type tmp( val, 0 ) ;
      _data.v_c->push_back( tmp ) ;
    } else if ( _data_type == GC_VECTOR ) {
      _data.v->resize(_data.v->size()+1) ;
      _data.v->back().set_real( val ) ;
    } else {
      if ( _data_type != GC_VEC_REAL ) promote_to_vec_real() ;
      _data.v_r->push_back( val ) ;
    }
  }

  void
  GenericContainer::push_complex( complex_type & val ) {
    if ( _data_type == GC_VECTOR ) {
      _data.v->resize(_data.v->size()+1) ;
      _data.v->back().set_complex( val ) ;
    } else {
      if ( _data_type != GC_VEC_COMPLEX ) promote_to_vec_complex() ;
      _data.v_c->push_back( val ) ;
    }
  }

  void
  GenericContainer::push_complex( real_type re, real_type im ) {
    complex_type tmp( re, im ) ;
    push_complex( tmp ) ;
  }

  void
  GenericContainer::push_string( string_type const & val ) {
    if ( _data_type != GC_VEC_STRING ) promote_to_vector() ;
    if ( _data_type == GC_VEC_STRING ) {
      _data.v_s->push_back( val ) ;
    } else {
      _data.v->resize(_data.v->size()+1) ;
      _data.v->back().set_string( val ) ;
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
  GenericContainer::get_pvoid( char const msg[] ) const {
    ck(msg,GC_POINTER) ;
    return _data.p ;
  }

  void **
  GenericContainer::get_ppvoid( char const msg[] ) const {
    ck(msg,GC_POINTER) ;
    return const_cast<void **>(&_data.p) ;
  }

  template <>
  void
  GenericContainer::get_value( unsigned & v, char const msg[] ) const {
    switch (_data_type) {
    case GC_BOOL:
      v = _data.b?1:0 ;
      break ;
    case GC_INTEGER:
      GC_ASSERT( _data.i >= 0,
                 msg <<
                 "value " << _data.i <<
                 "' cannot be converted in `unsigned'" ) ;
      v = unsigned(_data.i) ;
      break ;
    case GC_LONG:
      GC_ASSERT( _data.l >= 0,
                 msg <<
                 "value " << _data.l <<
                 "' cannot be converted in `unsigned'" ) ;
      v = unsigned(_data.l) ;
      break ;
    case GC_REAL:
      GC_ASSERT( _data.r >= 0 && isUnsigned(_data.r),
                 msg <<
                 "value " << _data.r <<
                 "' cannot be converted in `unsigned'" ) ;
      v = unsigned(_data.r) ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()) && isUnsigned(_data.c->real()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `unsigned'" ) ;
      v = unsigned(_data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `unsigned'" ) ;
    }
  }

  template <>
  void
  GenericContainer::get_value( int & v, char const msg[] ) const {
    switch (_data_type) {
    case GC_BOOL:
      v = _data.b?1:0 ;
      break ;
    case GC_INTEGER:
      GC_ASSERT( _data.i >= 0,
                 msg <<
                 "value " << _data.i <<
                 "' cannot be converted in `int'" ) ;
      v = _data.i ;
      break ;
    case GC_LONG:
      v = int(_data.l) ;
      break ;
    case GC_REAL:
      GC_ASSERT( isInteger(_data.r),
                 msg <<
                 "value " << _data.r <<
                 "' cannot be converted in `int'" ) ;
      v = int(_data.r) ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()) && isInteger(_data.c->real()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `int'" ) ;
      v = int(_data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `unsigned'" ) ;
    }
  }

  template <>
  void
  GenericContainer::get_value( unsigned long & v, char const msg[] ) const {
    typedef unsigned long ul ;
    switch (_data_type) {
    case GC_BOOL:
      v = _data.b?1:0 ;
      break ;
    case GC_INTEGER:
      GC_ASSERT( _data.i >= 0,
                 msg <<
                 "value " << _data.i <<
                 "' cannot be converted in `unsigned long'" ) ;
      v = ul(_data.i) ;
      break ;
    case GC_LONG:
      GC_ASSERT( _data.l >= 0,
                 msg <<
                 "value " << _data.l <<
                 "' cannot be converted in `unsigned long'" ) ;
      v = ul(_data.l) ;
      break ;
    case GC_REAL:
      GC_ASSERT( _data.r >= 0 && isUnsigned(_data.r),
                 msg <<
                 "value " << _data.r <<
                 "' cannot be converted in `unsigned long'" ) ;
      v = ul(_data.r) ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()) && isUnsigned(_data.c->real()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `unsigned long'" ) ;
      v = ul(_data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `unsigned long'" ) ;
    }
  }

  template <>
  void
  GenericContainer::get_value( long & v, char const msg[] ) const {
    switch (_data_type) {
    case GC_BOOL:
      v = _data.b?1:0 ;
      break ;
    case GC_INTEGER:
      v = long(_data.i) ;
      break ;
    case GC_LONG:
      v = _data.l ;
      break ;
    case GC_REAL:
      GC_ASSERT( isInteger(_data.r),
                 msg <<
                 "value " << _data.r <<
                 "' cannot be converted in `long'" ) ;
      v = long(_data.r) ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()) && isInteger(_data.c->real()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `long'" ) ;
      v = long(_data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `long'" ) ;
    }
  }

  template <>
  void
  GenericContainer::get_value( float & v, char const msg[] ) const {
    switch (_data_type) {
    case GC_BOOL:
      v = float(_data.b?1:0) ;
      break ;
    case GC_INTEGER:
      v = float(_data.i) ;
      break ;
    case GC_LONG:
      v = float(_data.l) ;
      break ;
    case GC_REAL:
      v = float(_data.r) ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `float'" ) ;
      v = float(_data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `float'" ) ;
    }
  }

  template <>
  void
  GenericContainer::get_value( double & v, char const msg[] ) const {
    switch (_data_type) {
    case GC_BOOL:
      v = _data.b?1:0 ;
      break ;
    case GC_INTEGER:
      v = double(_data.i) ;
      break ;
    case GC_LONG:
      v = double(_data.l) ;
      break ;
    case GC_REAL:
      v = _data.r ;
      break ;
    case GC_COMPLEX:
      GC_ASSERT( isZero(_data.c->imag()),
                 msg <<
                 "value " << *_data.c <<
                 "' cannot be converted in `double'" ) ;
      v = _data.c->real() ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `double'" ) ;
    }
  }

  //! If data is boolean, integer or floating point return number, otherwise return `0`.
  real_type
  GenericContainer::get_number() const {
    switch (_data_type) {
    case GC_BOOL:    return real_type(_data.b?1:0) ;
    case GC_INTEGER: return real_type(_data.i) ;
    case GC_LONG:    return real_type(_data.l) ;
    case GC_REAL:    return _data.r ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    break ;
    }
    return 0 ;
  }

  //! If data is boolean, integer or floating point return number, otherwise return `0`.
  complex_type
  GenericContainer::get_complex_number() const {
    switch (_data_type) {
    case GC_BOOL:    return complex_type(_data.b?1:0) ;
    case GC_INTEGER: return complex_type(_data.i,0) ;
    case GC_LONG:    return complex_type(_data.l,0) ;
    case GC_REAL:    return complex_type(_data.r,0) ;
    case GC_COMPLEX: return *_data.c ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
    break ;
    }
    return 0 ;
  }

  void
  GenericContainer::get_complex_number( real_type & re, real_type & im ) const {
    complex_type tmp = get_complex_number() ;
    re = tmp.real() ;
    im = tmp.imag() ;
  }

  real_type
  GenericContainer::get_number_at( unsigned i ) const {
    switch (_data_type) {
    case GC_VEC_BOOL:    return (*_data.v_b)[i] ? 1 : 0 ;
    case GC_VEC_INTEGER: return real_type((*_data.v_i)[i]) ;
    case GC_VEC_REAL:    return real_type((*_data.v_r)[i]) ;
    case GC_MAT_REAL:    return real_type((*_data.m_r)[i]) ;
    case GC_VECTOR:      return (*_data.v)[i].get_number() ;
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
    break ;
    }
    return 0 ;
  }

  complex_type
  GenericContainer::get_complex_number_at( unsigned i ) const {
    switch (_data_type) {
    case GC_VEC_BOOL:    return complex_type((*_data.v_b)[i]) ;
    case GC_VEC_INTEGER: return complex_type((*_data.v_i)[i],0) ;
    case GC_VEC_REAL:    return complex_type((*_data.v_r)[i],0) ;
    case GC_VEC_COMPLEX: return (*_data.v_c)[i] ;
    case GC_MAT_REAL:    return complex_type((*_data.m_r)[i],0) ;
    case GC_MAT_COMPLEX: return (*_data.m_c)[i] ;
    case GC_VECTOR:      return (*_data.v)[i].get_complex_number() ;
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
    return 0 ;
  }

  void
  GenericContainer::get_complex_number_at( unsigned i, real_type & re, real_type & im ) const {
    complex_type tmp = get_complex_number_at(i) ;
    re = tmp.real() ;
    im = tmp.imag() ;
  }

  bool_type &
  GenericContainer::get_bool( char const msg[] ) {
    ck_or_set(msg,GC_BOOL) ;
    return _data.b ;
  }

  bool_type const &
  GenericContainer::get_bool( char const msg[] ) const {
    ck(msg,GC_BOOL) ;
    return _data.b ;
  }

  int_type &
  GenericContainer::get_int( char const msg[] ) {
    ck_or_set(msg,GC_INTEGER) ;
    return _data.i ;
  }

  int_type const &
  GenericContainer::get_int( char const msg[] ) const {
    ck(msg,GC_INTEGER) ;
    return _data.i ;
  }

  long_type &
  GenericContainer::get_long( char const msg[] ) {
    ck_or_set(msg,GC_LONG) ;
    return _data.l ;
  }

  long_type const &
  GenericContainer::get_long( char const msg[] ) const {
    ck(msg,GC_LONG) ;
    return _data.l ;
  }

  real_type &
  GenericContainer::get_real( char const msg[] ) {
    ck_or_set(msg,GC_REAL) ;
    return _data.r ;
  }

  real_type const &
  GenericContainer::get_real( char const msg[] ) const {
    ck(msg,GC_REAL) ;
    return _data.r ;
  }

  complex_type &
  GenericContainer::get_complex( char const msg[] ) {
    ck_or_set(msg,GC_COMPLEX) ;
    return *_data.c ;
  }

  complex_type const &
  GenericContainer::get_complex( char const msg[] ) const {
    ck(msg,GC_REAL) ;
    return *_data.c ;
  }

  string_type &
  GenericContainer::get_string( char const msg[] ) {
    ck_or_set(msg,GC_STRING) ;
    return *_data.s ;
  }

  string_type const &
  GenericContainer::get_string( char const msg[] ) const {
    ck(msg,GC_STRING) ;
    return *_data.s ;
  }

  vector_type &
  GenericContainer::get_vector( char const msg[] ) {
    ck(msg,GC_VECTOR) ;
    return *_data.v ;
  }

  vector_type const &
  GenericContainer::get_vector( char const msg[] ) const {
    ck(msg,GC_VECTOR) ;
    return *_data.v ;
  }

  vec_pointer_type &
  GenericContainer::get_vec_pointer( char const msg[] ) {
    ck(msg,GC_VEC_POINTER) ;
    return *_data.v_p ;
  }

  vec_pointer_type const &
  GenericContainer::get_vec_pointer( char const msg[] ) const {
    ck(msg,GC_VEC_POINTER) ;
    return *_data.v_p ;
  }

  vec_bool_type &
  GenericContainer::get_vec_bool( char const msg[] ) {
    ck(msg,GC_VEC_BOOL) ;
    return *_data.v_b ;
  }

  vec_bool_type const &
  GenericContainer::get_vec_bool( char const msg[] ) const {
    ck(msg,GC_VEC_BOOL) ;
    return *_data.v_b ;
  }

  vec_int_type &
  GenericContainer::get_vec_int( char const msg[] ) {
    if ( _data_type == GC_NOTYPE   ) set_vec_int() ;
    if ( _data_type == GC_VEC_BOOL ) promote_to_vec_int() ;
    ck(msg,GC_VEC_INTEGER) ;
    return *_data.v_i ;
  }

  vec_int_type const &
  GenericContainer::get_vec_int( char const msg[] ) const {
    ck(msg,GC_VEC_INTEGER) ;
    return *_data.v_i ;
  }

  vec_real_type &
  GenericContainer::get_vec_real( char const msg[] ) {
    if ( _data_type == GC_NOTYPE   ) set_vec_real() ;
    if ( _data_type == GC_VEC_BOOL ||
         _data_type == GC_VEC_INTEGER ) promote_to_vec_real() ;
    ck(msg,GC_VEC_REAL) ;
    return *_data.v_r ;
  }

  vec_real_type const &
  GenericContainer::get_vec_real( char const msg[] ) const {
    ck(msg,GC_VEC_REAL) ;
    return *_data.v_r ;
  }

  vec_complex_type &
  GenericContainer::get_vec_complex( char const msg[] ) {
    if ( _data_type == GC_NOTYPE ) set_vec_complex() ;
    if ( _data_type == GC_VEC_BOOL    ||
         _data_type == GC_VEC_INTEGER ||
         _data_type == GC_VEC_REAL ) promote_to_vec_complex() ;
    ck(msg,GC_VEC_COMPLEX) ;
    return *_data.v_c ;
  }

  vec_complex_type const &
  GenericContainer::get_vec_complex( char const msg[] ) const {
    ck(msg,GC_VEC_COMPLEX) ;
    return *_data.v_c ;
  }

  mat_real_type &
  GenericContainer::get_mat_real( char const msg[] ) {
    if ( _data_type == GC_NOTYPE ) set_mat_real() ;
    if ( _data_type == GC_VEC_BOOL    ||
         _data_type == GC_VEC_INTEGER ||
         _data_type == GC_VEC_REAL ) promote_to_mat_real() ;
    ck(msg,GC_MAT_REAL) ;
    return *_data.m_r ;
  }

  mat_real_type const &
  GenericContainer::get_mat_real( char const msg[] ) const {
    ck(msg,GC_MAT_REAL) ;
    return *_data.m_r ;
  }

  mat_complex_type &
  GenericContainer::get_mat_complex( char const msg[] ) {
    if ( _data_type == GC_NOTYPE ) set_mat_complex() ;
    if ( _data_type == GC_VEC_BOOL    ||
         _data_type == GC_VEC_INTEGER ||
         _data_type == GC_VEC_REAL    ||
         _data_type == GC_MAT_REAL    ||
         _data_type == GC_VEC_COMPLEX ) promote_to_mat_complex() ;
    ck(msg,GC_MAT_COMPLEX) ;
    return *_data.m_c ;
  }

  mat_complex_type const &
  GenericContainer::get_mat_complex( char const msg[] ) const {
    ck(msg,GC_MAT_COMPLEX) ;
    return *_data.m_c ;
  }

  vec_string_type &
  GenericContainer::get_vec_string( char const msg[] ) {
    ck(msg,GC_VEC_STRING) ;
    return *_data.v_s ;
  }

  vec_string_type const &
  GenericContainer::get_vec_string( char const msg[] ) const {
    ck(msg,GC_VEC_STRING) ;
    return *_data.v_s ;
  }

  map_type &
  GenericContainer::get_map( char const msg[] ) {
    ck(msg,GC_MAP) ;
    return *_data.m ;
  }

  map_type const &
  GenericContainer::get_map( char const msg[] ) const {
    ck(msg,GC_MAP) ;
    return *_data.m ;
  }

  // --------------------------------------------------------------
  void
  GenericContainer::copyto_vec_int( vec_int_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.resize(ne) ;
      copy( _data.v_b->begin(), _data.v_b->end(), v.begin() ) ;
      break ;
    case GC_VEC_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = (*_data.v_r)[i] ;
        GC_ASSERT( isInteger(val),
                   msg << "copyto_vec_int: v[" << i << "] = " << val <<
                   " cannot be converted to `integer'" ) ;
        v.push_back( int_type(val) ) ;
      }
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type val = get_complex_number_at(i) ;
        GC_ASSERT( isZero(val.imag()) && isInteger(val.real()),
                   msg << "copyto_vec_int: v[" << i << "] = " << val <<
                   " cannot be converted to `integer'" ) ;
        v.push_back( int_type(val.real()) ) ;
      }
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = get_number_at(i) ;
        GC_ASSERT( isInteger(val),
                   msg << "copyto_vec_int: v[" << i << "] = " << val <<
                   " cannot be converted to `integer'" ) ;
        v.push_back( int_type(val) ) ;
      }
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      v.push_back( _data.i ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      v.push_back( int_type(_data.l) ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      GC_ASSERT( isInteger(_data.r),
                msg << "copyto_vec_int: value " << _data.r <<
                " cannot be converted to `integer'" ) ;
      v.push_back( int_type(_data.r) ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      GC_ASSERT( isZero(_data.c->imag()) && isInteger(_data.c->real()),
                 msg << "copyto_vec_int: value " << *_data.c <<
                 " cannot be converted to `integer'" ) ;
      v.push_back( int_type(_data.c->real()) ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_int_type'" ) ;
    }
  }

  void
  GenericContainer::copyto_vec_uint( vec_uint_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        int_type val = (*_data.v_i)[i] ;
        GC_ASSERT( val >= 0,
                   msg << "copyto_vec_uint: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned integer'" ) ;
        v.push_back( uint_type(val) ) ;
      }
      break ;
    case GC_VEC_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = (*_data.v_r)[i] ;
        GC_ASSERT( isUnsigned(val),
                   msg << "copyto_vec_uint: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned integer'" ) ;
        v.push_back( uint_type(val) ) ;
      }
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type val = get_complex_number_at(i) ;
        GC_ASSERT( isZero(val.imag()) && isUnsigned(val.real()),
                   msg << "copyto_vec_uint: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned integer'" ) ;
        v.push_back( uint_type(val.real()) ) ;
      }
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = get_number_at(i) ;
        GC_ASSERT( isUnsigned(val),
                   msg << "copyto_vec_uint: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned integer'" ) ;
        v.push_back( uint_type(val) ) ;
      }
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      GC_ASSERT( _data.i >= 0,
                 msg << "copyto_vec_uint: value = " << _data.i <<
                 " cannot be converted to `unsigned integer'" ) ;
      v.push_back( uint_type(_data.i) ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      GC_ASSERT( _data.l >= 0,
                 msg << "copyto_vec_uint: value = " << _data.l <<
                 " cannot be converted to `unsigned integer'" ) ;
      v.push_back( uint_type(_data.l) ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      GC_ASSERT( isUnsigned(_data.r),
                msg << "copyto_vec_uint: value = " << _data.r <<
                " cannot be converted to `unsigned integer'" ) ;
      v.push_back( uint_type(_data.r) ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      GC_ASSERT( isZero(_data.c->imag()) && isUnsigned(_data.c->real()),
                 msg << "copyto_vec_uint: value " << *_data.c <<
                 " cannot be converted to `unsigned integer'" ) ;
      v.push_back( uint_type(_data.c->real()) ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_uint_type'" ) ;
    }
  }

  void
  GenericContainer::copyto_vec_long( vec_long_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        int_type val = (*_data.v_i)[i] ;
        v.push_back( long_type(val) ) ;
      }
      break ;
    case GC_VEC_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = (*_data.v_r)[i] ;
        GC_ASSERT( isInteger(val),
                   msg << "copyto_vec_long: v[" << i << "] = " << val <<
                   " cannot be converted to `long'" ) ;
        v.push_back( long_type(val) ) ;
      }
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type val = get_complex_number_at(i) ;
        GC_ASSERT( isZero(val.imag()) && isInteger(val.real()),
                   msg << "copyto_vec_long: v[" << i << "] = " << val <<
                   " cannot be converted to `long'" ) ;
        v.push_back( long_type(val.real()) ) ;
      }
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = get_number_at(i) ;
        GC_ASSERT( isInteger(val),
                   msg << "copyto_vec_long: v[" << i << "] = " << val <<
                   " cannot be converted to `long'" ) ;
        v.push_back( long_type(val) ) ;
      }
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      v.push_back( long_type(_data.i) ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      v.push_back( _data.l ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      GC_ASSERT( isInteger(_data.r),
                msg << "copyto_vec_long: value " << _data.r <<
                " cannot be converted to `long'" ) ;
      v.push_back( long_type(_data.r) ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      GC_ASSERT( isZero(_data.c->imag()) && isInteger(_data.c->real()),
                 msg << "copyto_vec_long: value " << *_data.c <<
                 " cannot be converted to `long'" ) ;
      v.push_back( long_type(_data.c->real()) ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_long_type'" ) ;
    }

  }

  void
  GenericContainer::copyto_vec_ulong( vec_ulong_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        int_type val = (*_data.v_i)[i] ;
        GC_ASSERT( val >= 0,
                   msg << "copyto_vec_ulong: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned long'" ) ;
        v.push_back( ulong_type(val) ) ;
      }
      break ;
    case GC_VEC_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = (*_data.v_r)[i] ;
        GC_ASSERT( isInteger(val),
                   msg << "copyto_vec_ulong: v[" << i << " = " << val <<
                   " cannot be converted to `unsigned long'" ) ;
        v.push_back( ulong_type(val) ) ;
      }
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type val = get_complex_number_at(i) ;
        GC_ASSERT( isZero(val.imag()) && isUnsigned(val.real()),
                   msg << "copyto_vec_ulong: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned long'" ) ;
        v.push_back( ulong_type(val.real()) ) ;
      }
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        real_type val = get_number_at(i) ;
        GC_ASSERT( isUnsigned(val),
                   msg << "copyto_vec_ulong: v[" << i << "] = " << val <<
                   " cannot be converted to `unsigned long'" ) ;
        v.push_back( ulong_type(val) ) ;
      }
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      GC_ASSERT( _data.i >= 0,
                 msg << "copyto_vec_ulong: value = " << _data.i <<
                 " cannot be converted to `unsigned long'" ) ;
      v.push_back( ulong_type(_data.i) ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      GC_ASSERT( _data.l >= 0,
                 msg << "copyto_vec_ulong: value = " << _data.l <<
                 " cannot be converted to `unsigned long'" ) ;
      v.push_back( ulong_type(_data.l) ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      GC_ASSERT( isUnsigned(_data.r),
                msg << "copyto_vec_ulong: value = " << _data.r <<
                " cannot be converted to `unsigned long'" ) ;
      v.push_back( ulong_type(_data.r) ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      GC_ASSERT( isZero(_data.c->imag()) && isUnsigned(_data.c->real()),
                 msg << "copyto_vec_ulong: value " << *_data.c <<
                 " cannot be converted to `unsigned long'" ) ;
      v.push_back( ulong_type(_data.c->real()) ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_ulong_type'" ) ;
    }

  }

  void
  GenericContainer::copyto_vec_real( vec_real_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( real_type( (*_data.v_i)[i] ) ) ;
      break ;
    case GC_VEC_REAL:
      v.resize(ne) ;
      copy( _data.v_r->begin(), _data.v_r->end(), v.begin() ) ;
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type val = get_complex_number_at(i) ;
        GC_ASSERT( isZero(val.imag()),
                   msg << "copyto_vec_real: v[" << i << "] = " << val <<
                   " cannot be converted to `real'" ) ;
        v.push_back( val.real() ) ;
      }
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i )
        v.push_back( get_number_at(i) ) ;
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      v.push_back( real_type(_data.i) ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      v.push_back( real_type(_data.l) ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      v.push_back( _data.r ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      GC_ASSERT( isZero(_data.c->imag()),
                 msg << "copyto_vec_real: value " << *_data.c <<
                 " cannot be converted to `real'" ) ;
      v.push_back( _data.c->real() ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_real_type'" ) ;
    }
  }

  void
  GenericContainer::copyto_vec_complex( vec_complex_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_VEC_BOOL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) v.push_back( (*_data.v_b)[i] ? 1 : 0 ) ;
      break ;
    case GC_VEC_INTEGER:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type tmp(real_type( (*_data.v_i)[i] ),0) ;
        v.push_back( tmp ) ;
      }
      break ;
    case GC_VEC_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type tmp((*_data.v_r)[i],0) ;
        v.push_back( tmp ) ;
      }
      break ;
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i )
        v.push_back( get_complex_number_at(i) ) ;
      break ;
    case GC_VECTOR:
    case GC_MAT_REAL:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        complex_type tmp(get_number_at(i),0) ;
        v.push_back( tmp ) ;
      }
      break;
    case GC_BOOL:
      v.reserve(1) ;
      v.push_back( _data.b ? 1 : 0 ) ;
      break ;
    case GC_INTEGER:
      v.reserve(1) ;
      v.push_back( complex_type( real_type(_data.i), 0 ) ) ;
      break ;
    case GC_LONG:
      v.reserve(1) ;
      v.push_back( complex_type( real_type(_data.l), 0 ) ) ;
      break ;
    case GC_REAL:
      v.reserve(1) ;
      v.push_back( complex_type( _data.r, 0 ) ) ;
      break ;
    case GC_COMPLEX:
      v.reserve(1) ;
      v.push_back( *_data.c ) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_complex_type'" ) ;
    }
  }

  void
  GenericContainer::copyto_vec_string( vec_string_type & v, char const msg[] ) const {
    v.clear() ;
    unsigned ne = get_num_elements() ;
    switch (_data_type) {
    case GC_STRING:
      v.reserve(ne) ;
      v.push_back( *_data.s ) ;
      break;
    case GC_VEC_STRING:
      v.resize(ne) ;
      std::copy( _data.v_s->begin(), _data.v_s->end(), v.begin() ) ;
      break;
    case GC_VECTOR:
      v.reserve(ne) ;
      for ( unsigned i = 0 ; i < ne ; ++i ) {
        GenericContainer const & gc = get_gc_at(i,msg) ;
        v.push_back( gc.get_string(msg) ) ;
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
    case GC_VEC_REAL:
    case GC_MAT_REAL:
    case GC_VEC_COMPLEX:
    case GC_MAT_COMPLEX:
    case GC_POINTER:
    case GC_VEC_POINTER:
    case GC_MAP:
    GC_DO_ERROR( msg <<
                 "\nbad data type: `" << typeName[_data_type] <<
                 "' cannot be converted in `vec_string_type'" ) ;
    }
  }

  // --------------------------------------------------------------

  bool
  GenericContainer::exists( std::string const & s ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(s) ;
    return iv != (*_data.m).end() ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], bool & value ) const {
    if ( _data_type != GC_MAP ) return false;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    if ( iv->second._data_type != GC_BOOL ) return false ;
    value = iv->second._data.b ;
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], int_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = iv->second._data.b?1:0 ;
      break ;
    case GC_INTEGER:
      value = iv->second._data.i ;
      break ;
    case GC_LONG:
      value = int_type(iv->second._data.l) ;
      break ;
    case GC_REAL:
      if ( isInteger(iv->second._data.r) ) value = int_type(iv->second._data.r) ;
      break ;
    case GC_COMPLEX:
      if ( isInteger(iv->second._data.c->real()) && isZero(iv->second._data.c->imag()) )
        value = int_type(iv->second._data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], uint_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = iv->second._data.b?1:0 ;
      break ;
    case GC_INTEGER:
      if ( isUnsigned(iv->second._data.i) ) value = uint_type(iv->second._data.i) ;
      break ;
    case GC_LONG:
      if ( isUnsigned(iv->second._data.l) ) value = uint_type(iv->second._data.l) ;
      break ;
    case GC_REAL:
      if ( isUnsigned(iv->second._data.r) ) value = uint_type(iv->second._data.r) ;
      break ;
    case GC_COMPLEX:
      if ( isUnsigned(iv->second._data.c->real()) && isZero(iv->second._data.c->imag()) )
        value = uint_type(iv->second._data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], long_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = iv->second._data.b?1:0 ;
      break ;
    case GC_INTEGER:
      value = long_type(iv->second._data.i) ;
      break ;
    case GC_LONG:
      value = iv->second._data.l ;
      break ;
    case GC_REAL:
      if ( isInteger(iv->second._data.r) ) value = long_type(iv->second._data.r) ;
      break ;
    case GC_COMPLEX:
      if ( isInteger(iv->second._data.c->real()) && isZero(iv->second._data.c->imag()) )
        value = long_type(iv->second._data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], ulong_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = iv->second._data.b?1:0 ;
      break ;
    case GC_INTEGER:
      if ( isUnsigned(iv->second._data.i) ) value = ulong_type(iv->second._data.i) ;
      break ;
    case GC_LONG:
      if ( isUnsigned(iv->second._data.l) ) value = ulong_type(iv->second._data.l) ;
      break ;
    case GC_REAL:
      if ( isUnsigned(iv->second._data.r) ) value = ulong_type(iv->second._data.r) ;
      break ;
    case GC_COMPLEX:
      if ( isUnsigned(iv->second._data.c->real()) && isZero(iv->second._data.c->imag()) )
        value = ulong_type(iv->second._data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], real_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = real_type(iv->second._data.b?1:0) ;
      break ;
    case GC_INTEGER:
      value = real_type(iv->second._data.i) ;
      break ;
    case GC_LONG:
      value = real_type(iv->second._data.l) ;
      break ;
    case GC_REAL:
      value = iv->second._data.r ;
      break ;
    case GC_COMPLEX:
      if ( isZero(iv->second._data.c->imag()) )
        value = real_type(iv->second._data.c->real()) ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], complex_type & value ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    switch (iv->second._data_type) {
    case GC_BOOL:
      value = complex_type(iv->second._data.b?1:0,0) ;
      break ;
    case GC_INTEGER:
      value = complex_type(iv->second._data.i,0) ;
      break ;
    case GC_LONG:
      value = complex_type(iv->second._data.l,0) ;
      break ;
    case GC_REAL:
      value = complex_type(iv->second._data.r,0) ;
      break ;
    case GC_COMPLEX:
      value = *iv->second._data.c ;
      break ;
    case GC_NOTYPE:
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      return false ;
    }
    return true ;
  }

  bool
  GenericContainer::get_if_exists( char const field[], string_type & value ) const {
    if ( _data_type != GC_MAP ) return false;
    map_type::iterator iv = (*_data.m).find(field) ;
    if ( iv == (*_data.m).end() ) return false ;
    if ( iv->second._data_type != GC_STRING ) return false ;
    value = *iv->second._data.s ;
    return true ;
  }

  // --------------------------------------------------------------
  bool_type
  GenericContainer::get_bool_at( unsigned i ) {
    if ( _data_type == GC_NOTYPE   ) set_vec_bool() ;
    if ( _data_type == GC_VEC_BOOL ) {
      CHECK_RESIZE(_data.v_b,i) ; // correct type, check size
      return (*_data.v_b)[i] ;
    } else {
      if ( _data_type != GC_VECTOR ) promote_to_vector() ;
      CHECK_RESIZE(_data.v,i) ;
      return (*_data.v)[i].set_bool(false) ;
    }
  }

  bool_type
  GenericContainer::get_bool_at( unsigned i, char const msg[] ) const {
    ck(msg,GC_VEC_BOOL) ;
    GC_ASSERT( i < _data.v_b->size(), msg << "\nget_bool_at( " << i << " ) const, out of range" ) ;
    return (*_data.v_b)[i] ;
  }

  int_type &
  GenericContainer::get_int_at( unsigned i ) {
    if      ( _data_type == GC_NOTYPE ) set_vec_int() ;
    else if ( _data_type == GC_VEC_BOOL ) promote_to_vec_int() ;
    if ( _data_type == GC_VEC_INTEGER ) {
      CHECK_RESIZE(_data.v_i,i) ; // correct type, check size
      return (*_data.v_i)[i] ;
    } else {
      if ( _data_type != GC_VECTOR ) promote_to_vector() ;
      CHECK_RESIZE(_data.v,i) ;
      return (*_data.v)[i].set_int(0) ;
    }
  }

  int_type const &
  GenericContainer::get_int_at( unsigned i, char const msg[] ) const {
    ck(msg,GC_VEC_INTEGER) ;
    GC_ASSERT( i < _data.v_i->size(), msg << "\nget_int_at( " << i << " ) const, out of range" ) ;
    return (*_data.v_i)[i] ;
  }

  real_type &
  GenericContainer::get_real_at( unsigned i ) {
    if      ( _data_type == GC_NOTYPE ) set_vec_real() ;
    else if ( _data_type == GC_VEC_BOOL || _data_type == GC_VEC_INTEGER ) promote_to_vec_real() ;
    if ( _data_type == GC_VEC_REAL ) {
      CHECK_RESIZE(_data.v_r,i) ; // correct type, check size
      return (*_data.v_r)[i] ;
    } else {
      if ( _data_type != GC_VECTOR ) promote_to_vector() ;
      CHECK_RESIZE(_data.v,i) ;
      return (*_data.v)[i].set_real(0) ;
    }
  }

  real_type const &
  GenericContainer::get_real_at( unsigned i, char const msg[] ) const  {
    ck(msg,GC_VEC_REAL) ;
    GC_ASSERT( i < _data.v_r->size(), msg << "\nget_real_at( " << i << " ) const, out of range" ) ;
    return (*_data.v_r)[i] ;
  }

  real_type &
  GenericContainer::get_real_at( unsigned i, unsigned j ) {
    if      ( _data_type == GC_NOTYPE ) set_mat_real(i,j) ;
    else if ( _data_type == GC_VEC_BOOL    ||
              _data_type == GC_VEC_INTEGER ||
              _data_type == GC_VEC_REAL ) promote_to_mat_real() ;
    GC_ASSERT( GC_MAT_REAL == _data_type,
               "get_real_at( " << i << ", " << j << " ) bad data type" <<
               "\nexpect: " << typeName[GC_MAT_REAL] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    return (*_data.m_r)(i,j) ;
  }

  real_type const &
  GenericContainer::get_real_at( unsigned i, unsigned j, char const msg[] ) const  {
    ck(msg,GC_MAT_REAL) ;
    GC_ASSERT( i < _data.v_r->size(),
               msg << "\bget_real_at( " << i << ", " << j << " ) const, out of range" ) ;
    return (*_data.m_r)(i,j) ;
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned i ) {
    if      ( _data_type == GC_NOTYPE ) set_vec_complex() ;
    else if ( _data_type == GC_VEC_BOOL    ||
              _data_type == GC_VEC_INTEGER ||
              _data_type == GC_VEC_REAL ) promote_to_vec_complex() ;
    if ( _data_type == GC_VEC_COMPLEX ) {
      CHECK_RESIZE(_data.v_c,i) ; // correct type, check size
      return (*_data.v_c)[i] ;
    } else {
      if ( _data_type != GC_VECTOR ) promote_to_vector() ;
      CHECK_RESIZE(_data.v,i) ;
      return (*_data.v)[i].set_complex(0,0) ;
    }
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned i, char const msg[] ) const  {
    ck(msg,GC_VEC_COMPLEX) ;
    GC_ASSERT( i < _data.v_c->size(),
               msg << "\nget_complex_at( " << i << " ) const, out of range" ) ;
    return (*_data.v_c)[i] ;
  }

  complex_type &
  GenericContainer::get_complex_at( unsigned i, unsigned j ) {
    if      ( _data_type == GC_NOTYPE ) set_mat_complex(i,j) ;
    else if ( _data_type == GC_VEC_BOOL    ||
              _data_type == GC_VEC_INTEGER ||
              _data_type == GC_VEC_REAL    ||
              _data_type == GC_VEC_COMPLEX ||
              _data_type == GC_MAT_REAL ) promote_to_mat_complex() ;
    GC_ASSERT( GC_MAT_COMPLEX == _data_type,
               "get_complex_at( " << i << ", " << j << " ) bad data type" <<
               "\nexpect: " << typeName[GC_MAT_COMPLEX] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    return (*_data.m_c)(i,j) ;
  }

  complex_type const &
  GenericContainer::get_complex_at( unsigned i, unsigned j, char const msg[] ) const  {
    ck(msg,GC_MAT_COMPLEX) ;
    return (*_data.m_c)(i,j) ;
  }

  string_type &
  GenericContainer::get_string_at( unsigned i ) {
    if ( _data_type == GC_NOTYPE ) set_vec_string() ;
    if ( _data_type == GC_VEC_STRING ) {
      CHECK_RESIZE(_data.v_s,i) ;
      return (*_data.v_s)[i] ;
    } else {
      promote_to_vector() ;
      return (*this)[i].set_string("") ;
    }
  }

  string_type const &
  GenericContainer::get_string_at( unsigned i, char const msg[] ) const {
    ck(msg,GC_VEC_STRING) ;
    GC_ASSERT( i < _data.v_s->size(),
               msg << "\nget_string_at( " << i << " ) const, out of range" ) ;
    return (*_data.v_s)[i] ;
  }
  
  GenericContainer &
  GenericContainer::get_gc_at( unsigned i )
  { return (*this)[i] ; }
  
  GenericContainer const &
  GenericContainer::get_gc_at( unsigned i, char const msg[] ) const {
    return (*this)(i,msg) ;
  }

  /*
  //   _        __
  //  (_)_ __  / _| ___
  //  | | '_ \| |_ / _ \
  //  | | | | |  _| (_) |
  //  |_|_| |_|_|  \___/
  */

  //! Print to stream the kind of data stored
  GenericContainer const &
  GenericContainer::info( std::basic_ostream<char> & stream ) const {
    switch ( _data_type ) {
    case GC_NOTYPE:
      stream << "GenericContainer: No data stored\n" ;
      break ;
    case GC_POINTER:
      stream << "Generic pointer: " << _data.p << '\n' ;
      break ;
    case GC_BOOL:
      stream << "Boolean: " << (_data.b?"true":"false") << '\n' ;
      break ;
    case GC_INTEGER:
      stream << "Integer: " << _data.i << '\n' ;
      break ;
    case GC_LONG:
      stream << "Long: " << _data.l << '\n' ;
      break ;
    case GC_REAL:
      stream << "Floating Point: " << _data.r << '\n' ;
      break ;
    case GC_COMPLEX:
      stream << "Complex Floating Point: [" << _data.c->real() << ", " << _data.c->imag() << " ]\n" ;
      break ;
    case GC_STRING:
      stream << "String: " << *_data.s << '\n' ;
      break ;
    case GC_VEC_POINTER:
      stream << "Vector of generic pointer of size " << _data.v_p->size() << '\n' ;
      break ;
    case GC_VEC_BOOL:
      stream << "Vector of boolean of size " << _data.v_b->size() << '\n' ;
      break ;
    case GC_VEC_INTEGER:
      stream << "Vector of integer of size " << _data.v_i->size() << '\n' ;
      break ;
    case GC_VEC_REAL:
      stream << "Vector of floating point number of size " << _data.v_r->size() << '\n' ;
      break ;
    case GC_VEC_COMPLEX:
      stream << "Vector of complex floating point number of size " << _data.v_c->size() << '\n' ;
      break ;
    case GC_VEC_STRING:
      stream << "Vector of string of size " << _data.v_s->size() << '\n' ;
      break ;
    case GC_MAT_REAL:
      _data.m_r->info(stream) ;
      break ;
    case GC_MAT_COMPLEX:
      _data.m_c->info(stream) ;
      break ;
    case GC_VECTOR:
      stream << "Vector of generic data type of size " << _data.v->size() << '\n' ;
      break ;
    case GC_MAP:
      stream << "Map\n" ;
      break ;
    //default:
    //  stream << "Type N. " << _data_type << " not recognized\n" ;
    //  break ;
    }
    return *this ;
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
      case 0: break ; // data present
      default: set_vector() ; // data must be allocated ;
    }
    CHECK_RESIZE(_data.v,i) ;
    return (*_data.v)[i] ;
  }

  GenericContainer const &
  GenericContainer::operator [] ( unsigned i ) const {
    GC_ASSERT( GC_VECTOR == _data_type,
               "operator [] integer argument = " << i <<
               "\nexpect: " << typeName[GC_VECTOR] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    GC_ASSERT( i < _data.v->size(), "operator [] const, index " << i << " out of range" ) ;
    return (*_data.v)[i] ;
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
  GenericContainer::operator () ( unsigned i, char const msg[] ) {
    ck(msg,GC_VECTOR) ;
    GC_ASSERT( i < _data.v->size(),
               msg << "\noperator () const, index " << i << " out of range" ) ;
    return (*_data.v)[i] ;
  }

  GenericContainer const &
  GenericContainer::operator () ( unsigned i, char const msg_in[] ) const {
    ck(msg_in,GC_VECTOR) ;
    char const * msg = msg_in == nullptr ? "" : msg_in ;
    GC_ASSERT( i < _data.v->size(),
               msg << "\noperator () const, index " << i << " out of range" ) ;
    return (*_data.v)[i] ;
  }

  GenericContainer &
  GenericContainer::operator () ( std::string const & s, char const msg_in[] ) {
    ck(msg_in,GC_MAP) ;
    map_type::iterator iv = (*_data.m) . find(s) ;
    char const * msg = msg_in == nullptr ? "" : msg_in ;
    GC_ASSERT( iv != (*_data.m) . end(),
               msg << "\noperator(): Cannot find key '" << s << "'!" ) ;
    return iv -> second ;
  }

  GenericContainer const &
  GenericContainer::operator () ( std::string const & s, char const msg_in[] ) const {
    ck(msg_in,GC_MAP) ;
    map_type::const_iterator iv = (*_data.m) . find(s) ;
    char const * msg = msg_in == nullptr ? "" : msg_in ;
    GC_ASSERT( iv != (*_data.m) . end(),
               msg << "\noperator(): Cannot find key '" << s << "'!" ) ;
    return iv -> second ;
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
    if ( ck( GC_MAP ) != 0 ) set_map() ; // if not data present allocate!
    return (*_data.m)[s] ;
  }

  GenericContainer const &
  GenericContainer::operator [] ( std::string const & s ) const {
    GC_ASSERT( GC_MAP == _data_type,
               "operator [] string argument ``" << s << "''"
               "\nexpect: " << typeName[GC_MAP] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    return (*_data.m)[s] ;
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
    switch (_data_type) {
    case GC_NOTYPE:
      set_int(0) ;
      break ;
    case GC_BOOL:
      set_int(_data.b?1:0) ;
      break ;
    case GC_INTEGER:
      break ;
    case GC_POINTER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_int() cannot promote " << get_type_name() << " to int") ;
    }
    return *this ;
  }

  GenericContainer const &
  GenericContainer::promote_to_long() {
    switch (_data_type) {
    case GC_NOTYPE:
      set_long(0) ;
      break ;
    case GC_BOOL:
      set_long(_data.b?1:0) ;
      break ;
    case GC_INTEGER:
      set_long(_data.i) ;
      break ;
    case GC_LONG:
      break ;
    case GC_POINTER:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_long() cannot promote " << get_type_name() << " to long") ;
    }
    return *this ;
  }

  GenericContainer const &
  GenericContainer::promote_to_real() {
    switch (_data_type) {
    case GC_NOTYPE:
      set_real(0) ;
      break ;
    case GC_BOOL:
      set_real(_data.b?1:0) ;
      break ;
    case GC_INTEGER:
      set_real(_data.i) ;
      break ;
    case GC_LONG:
      set_real(_data.l) ;
      break ;
    case GC_REAL:
      break ;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_BOOL:
    case GC_VEC_INTEGER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_real() cannot promote " << get_type_name() << " to real") ;
    }
    return *this ;
  }

  GenericContainer const &
  GenericContainer::promote_to_vec_int() {
    switch (_data_type) {
    case GC_NOTYPE:
    { set_vec_int(1) ; get_int_at(0) = 0 ; }
      break ;
    case GC_BOOL:
    { int_type tmp = _data.b?1:0 ; set_vec_int(1) ; get_int_at(0) = tmp ; }
      break ;
    case GC_INTEGER:
    { int_type tmp = _data.i ; set_vec_int(1) ; get_int_at(0) = tmp ; }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_vec_int(unsigned(v_b->size())) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.v_i)[i] = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      break ;
    case GC_POINTER:
    case GC_LONG:
    case GC_REAL:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_REAL:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_vec_int() cannot promote " << get_type_name() << " to vec_int_type") ;
    }
    return *this ;

  }

  //! If data contains vector of booleans or integer it is promoted to a vector of real.
  GenericContainer const &
  GenericContainer::promote_to_vec_real() {
    switch (_data_type) {
    case GC_NOTYPE:
      { set_vec_real(1) ; get_real_at(0) = 0 ; }
      break ;
    case GC_BOOL:
      { real_type tmp = _data.b?1:0 ; set_vec_real(1) ; get_real_at(0) = tmp ; }
      break ;
    case GC_INTEGER:
      { real_type tmp = _data.i ; set_vec_real(1) ; get_real_at(0) = tmp ; }
      break ;
    case GC_LONG:
      { real_type tmp = _data.l ; set_vec_real(1) ; get_real_at(0) = tmp ; }
      break ;
    case GC_REAL:
      { real_type tmp = _data.r ; set_vec_real(1) ; get_real_at(0) = tmp ; }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_vec_real(unsigned(v_b->size())) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.v_r)[i] = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = _data.v_i ;
        _data_type = GC_NOTYPE ;
        set_vec_real(unsigned(v_i->size())) ;
        for ( unsigned i = 0 ; i < v_i->size() ; ++i ) (*_data.v_r)[i] = (*v_i)[i] ;
        delete v_i ;
      }
      break ;
    case GC_VEC_REAL:
      break ;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_vec_real() cannot promote " << get_type_name() << " vec_real_type") ;
    }
    return *this ;
  }

  //! If data contains vector of booleans or integer it is promoted to a vector of real.
  GenericContainer const &
  GenericContainer::promote_to_vec_complex() {
    switch (_data_type) {
    case GC_NOTYPE:
      { set_vec_complex(1) ; get_complex_at(0) = 0 ; }
      break ;
    case GC_BOOL:
      { real_type tmp = _data.b?1:0 ; set_vec_complex(1) ; get_complex_at(0) = tmp ; }
      break ;
    case GC_INTEGER:
      { real_type tmp = _data.i ; set_vec_complex(1) ; get_complex_at(0) = tmp ; }
      break ;
    case GC_LONG:
      { real_type tmp = _data.l ; set_vec_complex(1) ; get_complex_at(0) = tmp ; }
      break ;
    case GC_REAL:
      { real_type tmp = _data.r ; set_vec_complex(1) ; get_complex_at(0) = tmp ; }
      break ;
    case GC_COMPLEX:
      { complex_type tmp = *_data.c ; set_vec_complex(1) ; get_complex_at(0) = tmp ; }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_vec_complex(unsigned(v_b->size())) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.v_c)[i] = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = _data.v_i ;
        _data_type = GC_NOTYPE ;
        set_vec_complex(unsigned(v_i->size())) ;
        for ( unsigned i = 0 ; i < v_i->size() ; ++i ) (*_data.v_c)[i] = (*v_i)[i] ;
        delete v_i ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type * v_r = _data.v_r ;
        _data_type = GC_NOTYPE ;
        set_vec_complex(unsigned(v_r->size())) ;
        for ( unsigned i = 0 ; i < v_r->size() ; ++i ) (*_data.v_c)[i] = (*v_r)[i] ;
        delete v_r ;
      }
      break ;
    case GC_VEC_COMPLEX:
      break ;
    case GC_POINTER:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_STRING:
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_vec_real() cannot promote " << get_type_name() << " to vec_complex_type") ;
    }
    return *this ;
  }

  //! If data contains vector of booleans, integer or real it is promoted to a vector of real.
  GenericContainer const &
  GenericContainer::promote_to_mat_real() {
    switch (_data_type) {
    case GC_NOTYPE:
      { set_mat_real(1,1) ; get_real_at(0,0) = 0 ; }
      break ;
    case GC_BOOL:
      { real_type tmp = _data.b?1:0 ; set_mat_real(1,1) ; get_real_at(0,0) = tmp ; }
      break ;
    case GC_INTEGER:
      { real_type tmp = _data.i ; set_mat_real(1,1) ; get_real_at(0,0) = tmp ; }
      break ;
    case GC_REAL:
      { real_type tmp = _data.r ; set_mat_real(1,1) ; get_real_at(0,0) = tmp ; }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_mat_real(unsigned(v_b->size()),1) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.m_r)(i,0) = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = _data.v_i ;
        _data_type = GC_NOTYPE ;
        set_mat_real(unsigned(v_i->size()),1) ;
        for ( unsigned i = 0 ; i < v_i->size() ; ++i ) (*_data.m_r)(i,0) = (*v_i)[i] ;
        delete v_i ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type * v_r = _data.v_r ;
        _data_type = GC_NOTYPE ;
        set_mat_real(unsigned(v_r->size()),1) ;
        for ( unsigned i = 0 ; i < v_r->size() ; ++i ) (*_data.m_r)(i,0) = (*v_r)[i] ;
        delete v_r ;
      }
      break ;
    case GC_MAT_REAL:
      break ;
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
      GC_DO_ERROR( ":promote_to_mat_real() cannot promote " << get_type_name() << " to mat_real_type") ;
    }
    return *this ;
  }

  //! If data contains vector of booleans, integer or real it is promoted to a vector of real.
  GenericContainer const &
  GenericContainer::promote_to_mat_complex() {
    switch (_data_type) {
    case GC_NOTYPE:
      { set_mat_complex(1,1) ; get_complex_at(0,0) = 0 ; }
      break ;
    case GC_BOOL:
      { real_type tmp = _data.b?1:0 ; set_mat_complex(1,1) ; get_complex_at(0,0) = tmp ; }
      break ;
    case GC_INTEGER:
      { real_type tmp = _data.i ; set_mat_complex(1,1) ; get_complex_at(0,0) = tmp ; }
      break ;
    case GC_LONG:
      { real_type tmp = _data.l ; set_mat_complex(1,1) ; get_complex_at(0,0) = tmp ; }
      break ;
    case GC_REAL:
      { real_type tmp = _data.r ; set_mat_complex(1,1) ; get_complex_at(0,0) = tmp ; }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_mat_complex(unsigned(v_b->size()),1) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.m_r)(i,0) = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = _data.v_i ;
        _data_type = GC_NOTYPE ;
        set_mat_complex(unsigned(v_i->size()),1) ;
        for ( unsigned i = 0 ; i < v_i->size() ; ++i ) (*_data.m_r)(i,0) = (*v_i)[i] ;
        delete v_i ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type * v_r = _data.v_r ;
        _data_type = GC_NOTYPE ;
        set_mat_complex(unsigned(v_r->size()),1) ;
        for ( unsigned i = 0 ; i < v_r->size() ; ++i ) (*_data.m_r)(i,0) = (*v_r)[i] ;
        delete v_r ;
      }
      break ;
    case GC_MAT_REAL:
      { mat_real_type * m_r = _data.m_r ;
        _data_type = GC_NOTYPE ;
        set_mat_complex(m_r->numRows(),m_r->numCols()) ;
        for ( unsigned i = 0 ; i < m_r->numRows() ; ++i )
          for ( unsigned j = 0 ; j < m_r->numCols() ; ++j )
            (*_data.m_c)(i,j) = (*m_r)(i,j) ;
        delete m_r ;
      }
      break ;
    case GC_MAT_COMPLEX:
      break ;
    case GC_POINTER:
    case GC_COMPLEX:
    case GC_STRING:
    case GC_VEC_POINTER:
    case GC_VEC_COMPLEX:
    case GC_VEC_STRING:
    case GC_VECTOR:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_mat_real() cannot promote " << get_type_name() << " to mat_complex_type") ;
    }
    return *this ;
  }


  //! If data contains vector of someting it is promoted to a vector of `GenericContainer`.
  GenericContainer const &
  GenericContainer::promote_to_vector() {
    switch (_data_type) {
    case GC_NOTYPE:
      { set_vector(1) ; (*this)[0].initialize() ; } // set data to no type
      break ;
    case GC_POINTER:
      { set_vector(1) ; (*this)[0] = _data.p ; }
      break ;
    case GC_BOOL:
      { set_vector(1) ; (*this)[0] = _data.b ; }
      break ;
    case GC_INTEGER:
      { set_vector(1) ; (*this)[0] = _data.i ; }
      break ;
    case GC_LONG:
      { set_vector(1) ; (*this)[0] = _data.l ; }
      break ;
    case GC_REAL:
      { set_vector(1) ; (*this)[0] = _data.r ; }
      break ;
    case GC_COMPLEX:
      { set_vector(1) ; (*this)[0] = *_data.c ; }
      break ;
    case GC_STRING:
      { set_vector(1) ; (*this)[0] = *_data.s ; }
      break ;
    case GC_VEC_POINTER:
      { vec_pointer_type * v_p = _data.v_p ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_p->size())) ;
        for ( unsigned i = 0 ; i < v_p->size() ; ++i ) (*_data.v)[i] = (*v_p)[i] ;
        delete v_p ;
      }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type * v_b = _data.v_b ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_b->size())) ;
        for ( unsigned i = 0 ; i < v_b->size() ; ++i ) (*_data.v)[i] = (*v_b)[i] ;
        delete v_b ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type * v_i = _data.v_i ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_i->size())) ;
        for ( unsigned i = 0 ; i < v_i->size() ; ++i ) (*_data.v)[i] = (*v_i)[i] ;
        delete v_i ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type * v_r = _data.v_r ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_r->size())) ;
        for ( unsigned i = 0 ; i < v_r->size() ; ++i ) (*_data.v)[i] = (*v_r)[i] ;
        delete v_r ;
      }
      break ;
    case GC_VEC_COMPLEX:
      { vec_complex_type * v_c = _data.v_c ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_c->size())) ;
        for ( unsigned i = 0 ; i < v_c->size() ; ++i ) (*_data.v)[i] = (*v_c)[i] ;
        delete v_c ;
      }
      break ;
    case GC_VEC_STRING:
      { vec_string_type * v_s = _data.v_s ;
        _data_type = GC_NOTYPE ;
        set_vector(unsigned(v_s->size())) ;
        for ( unsigned i = 0 ; i < v_s->size() ; ++i ) (*_data.v)[i] = (*v_s)[i] ;
        delete v_s ;
      }
      break ;
    case GC_VECTOR:
      break ;
    case GC_MAT_REAL:
    case GC_MAT_COMPLEX:
    case GC_MAP:
      GC_DO_ERROR( ":promote_to_vector() cannot promote " << get_type_name() << " to vector_type") ;
    }
    return *this ;
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
  GenericContainer::print( std::basic_ostream<char> & stream,
                           std::string const        & prefix,
                           std::string const        & indent ) const {

    switch (_data_type) {
    case GC_NOTYPE:
      stream << prefix << "Empty!\n" ;
      break ;
    case GC_POINTER:
      stream << prefix << this -> get_pvoid() << '\n' ;
      break ;
    case GC_BOOL:
      stream << prefix << (this -> get_bool()?"true":"false") << '\n' ;
      break ;
    case GC_INTEGER:
      stream << prefix << this -> get_int() << '\n' ;
      break ;
    case GC_LONG:
      stream << prefix << this -> get_long() << '\n' ;
      break ;
    case GC_REAL:
      stream << prefix << this -> get_real() << '\n' ;
      break ;
    case GC_COMPLEX:
      stream << prefix << "( " << this -> get_complex().real() << ", " << this -> get_complex().imag() << " )\n" ;
      break ;
    case GC_STRING:
      stream << prefix << "\"" << this -> get_string() << "\"\n" ;
      break ;
    case GC_VEC_POINTER:
      { vec_pointer_type const & v = this -> get_vec_pointer() ;
        for ( vec_pointer_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_pointer(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this -> get_vec_bool() ;
        stream << prefix << v << '\n' ; }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type const & v = this -> get_vec_int() ;
        stream << prefix << v << '\n' ; }
      break ;
    case GC_VEC_REAL:
      { vec_real_type const & v = this -> get_vec_real() ;
        stream << prefix << v << '\n' ; }
      break ;
    case GC_VEC_COMPLEX:
      { vec_complex_type const & v = this -> get_vec_complex() ;
        stream << prefix << v << '\n' ; }
      break ;
    case GC_MAT_REAL:
      { mat_real_type const & m = this -> get_mat_real() ;
        stream << m ; }
      break ;
    case GC_MAT_COMPLEX:
      { mat_complex_type const & m = this -> get_mat_complex() ;
        stream << m ; }
      break ;
    case GC_VEC_STRING:
      { vec_string_type const & v = this -> get_vec_string() ;
        stream << '\n';
        for ( vec_string_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix+indent << i << ": \"" << v[i] << "\"\n" ;
      }
      break ;

    case GC_VECTOR:
      { vector_type const & v = this -> get_vector() ;
        for ( vector_type::size_type i = 0 ; i < v.size() ; ++i ) {
          GenericContainer const & vi = v[i] ;
          if ( vi.simple_data() ||
               ( vi.simple_vec_data() && vi.get_num_elements() <= 10 ) ) {
            stream << prefix << i << ": " ;
            vi.print(stream,"") ;
          } else {
            stream << prefix << i << ":\n" ;
            vi.print(stream,prefix+indent) ;
          }
        }
      }
      break ;
    case GC_MAP:
      { map_type const & m = this -> get_map() ;
        for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
          // check formatting using pcre
          #ifndef GENERIC_CONTAINER_USE_CXX11
          if ( im->second.simple_data() ||
               ( im->second.simple_vec_data() && im->second.get_num_elements() <= 10 ) ) {
            stream << prefix << im->first << ": " ;
            im->second.print(stream,"") ;
          } else {
            stream << prefix << im->first << ":\n" ;
            im->second.print(stream,prefix+indent) ;
          }
          #else
          // num+"@"+"underline character"
          // Try to find the regex in aLineToMatch, and report results.
          std::string matches[4] ;
          int pcreExecRet = pcre_for_GC.exec( im->first, matches ) ;
          if ( pcreExecRet == 4 ) {
            std::string header = matches[3] ; // header
            // found formatting
            if ( im->second.simple_data() ) {
              stream << prefix << header << ": " ;
              im->second.print(stream,"") ;
            } else {
              if ( matches[1].length() > 1 ) stream << '\n' ; // double ## --> add nel line
              stream << prefix << header ;
              if ( matches[2].length() > 0 ) {
                stream << '\n' << prefix ;
                char fmt = matches[2][0] ; // underline char
                std::size_t m3 = header.length() ;
                while ( m3-- > 0 ) stream << fmt ; // underline header
              } else {
                stream << ':' ;
              }
              stream << '\n' ;
              im->second.print(stream,prefix+indent) ;
            }
          } else {
            std::string header = pcreExecRet == 3 ? matches[3] : im->first ;
            if ( im->second.simple_data() ) {
              stream << prefix << header << ": " ;
              im->second.print(stream,"") ;
            } else {
              stream << prefix << header << ":\n" ;
              im->second.print(stream,prefix+indent) ;
            }
          }
          #endif
        }
      }
      break ;

    //default:
    //  GC_DO_ERROR( "Error, print(...) unknown type!\n") ;
    //  break ;
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
  GenericContainer::to_yaml( std::basic_ostream<char> & stream, std::string const & prefix ) const {
    switch (_data_type) {
    case GC_NOTYPE:
      stream << "Empty!\n" ;
      break ;
    case GC_BOOL:
      stream << (this -> get_bool()?"true":"false") << '\n' ;
      break ;
    case GC_INTEGER:
      stream << this -> get_int() << '\n' ;
      break ;
    case GC_LONG:
      stream << this -> get_long() << '\n' ;
      break ;
    case GC_REAL:
      stream << this -> get_real() << '\n' ;
      break ;
    case GC_COMPLEX:
      stream << this -> get_complex().real() << ' ' 
             << this -> get_complex().imag() << '\n' ;
      break ;
    case GC_STRING:
      stream << "'" << this -> get_string() << "'\n" ;
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this -> get_vec_bool() ;
        stream << "[ " << (v[0]?"true":"false") ;
        for ( vec_bool_type::size_type i = 1 ; i < v.size() ; ++i )
          stream << ", " << (v[i]?"true":"false") ;
        stream << " ]\n" ;
      }
      break ;
    case GC_VEC_INTEGER:
      { vec_int_type const & v = this -> get_vec_int() ;
        stream << "[ " << v[0] ;
        for ( vec_int_type::size_type i = 1 ; i < v.size() ; ++i )
          stream << ", " << v[i] ;
        stream << " ]\n" ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type const & v = this -> get_vec_real() ;
        stream << "[ " << v[0] ;
        for ( vec_real_type::size_type i = 1 ; i < v.size() ; ++i )
          stream << ", " << v[i] ;
        stream << " ]\n" ;
      }
      break ;
    case GC_VEC_STRING:
      { vec_string_type const & v = this -> get_vec_string() ;
        stream << "[ '" << v[0] << "'" ;
        for ( vec_string_type::size_type i = 1 ; i < v.size() ; ++i )
          stream << ", '" << v[i] << "'" ;
        stream << " ]\n" ;
      }
      break ;

    case GC_VECTOR:
      { vector_type const & v = this -> get_vector() ;
        stream << '\n' ;
        for ( vector_type::size_type i = 0 ; i < v.size() ; ++i ) {
          stream << prefix << "- " ;
          v[i].to_yaml(stream,prefix+"  ") ;
        }
      }
      break ;
    case GC_MAP:
      { map_type const & m = this -> get_map() ;
        stream << '\n' ;
        for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
          stream << prefix << im->first << ": " ;
          im->second.to_yaml(stream,prefix+"  ") ;
        }
      }
      break ;
    case GC_MAT_REAL:
      case GC_VEC_COMPLEX:
      case GC_MAT_COMPLEX:
      case GC_POINTER:
      case GC_VEC_POINTER:
      { /* DA FARE */ }
      break ;
    }
  }

  void
  GenericContainer::exception( char const msg[] ) {
    throw std::runtime_error(msg) ;
  }
}

//
// eof: GenericContainer.cc
//
