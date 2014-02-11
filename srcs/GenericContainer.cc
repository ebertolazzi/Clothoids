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

// use pcre for pattern matching
#ifndef GENERIC_CONTAINER_NO_PCRE
  #include <pcre.h>
#endif

#define CHECK_RESIZE(pV,I) if ( pV->size() <= (I) ) pV->resize((I)+1)

namespace GC {

  static char const *typeName[] = {
    "NOTYPE",
    "pointer",
    "bool_type",
    "int_type",
    "real_type",
    "string_type",
    "vec_pointer_type",
    "vec_bool_type",
    "vec_int_type",
    "vec_real_type",
    "vec_string_type",
    "vector_type",
    "map_type"
  } ;

  #ifndef GENERIC_CONTAINER_NO_PCRE
  void * GenericContainer::_reCompiled = NULL ;
  void * GenericContainer::_pcreExtra  = NULL ;
  #endif

  // costruttore
  GenericContainer::GenericContainer()
  : _data_type(GC::GC_NOTYPE)
  {
    // compilo regex
    #ifndef GENERIC_CONTAINER_NO_PCRE
	static bool DoInit = true ; // initialize only once!
	if ( DoInit ) {
      pcre       *& reCompiled = *reinterpret_cast<pcre**>(&_reCompiled)      ; // pcre *
      pcre_extra *& pcreExtra  = *reinterpret_cast<pcre_extra**>(&_pcreExtra) ; // pcre_extra *

      reCompiled = pcre_compile("^\\s*\\d+\\s*(##?)(-|=|~|_|)\\s*(.*)$",
                                0,
                                &pcreErrorStr,
                                &pcreErrorOffset,
                                NULL);
      // pcre_compile returns NULL on error, and sets pcreErrorOffset & pcreErrorStr
      GC_ASSERT( reCompiled != nullptr,
                 "GenericContainer: Could not compile regex for print GenericContainer\n" ) ;
      // Optimize the regex
      pcreExtra = pcre_study(reCompiled, 0, &pcreErrorStr);
      GC_ASSERT( pcreExtra != nullptr,
                "GenericContainer: Could not optimize regex for print GenericContainer\n" ) ;
	  DoInit = false ;
	}
    #endif
  }

  // distruttore
  void
  GenericContainer::clear() {
    switch (_data_type) {
      case GC_POINTER:
        // removed annoying warning. To be re-thinked...
        //GC_WARNING( _data.p == nullptr, "find a pointer not deallocated!" ) ;
        break ;
      case GC_STRING:      delete _data.s   ; break ;

      case GC_VEC_POINTER: delete _data.v_p ; break ;
      case GC_VEC_BOOL:    delete _data.v_b ; break ;
      case GC_VEC_INT:     delete _data.v_i ; break ;
      case GC_VEC_REAL:    delete _data.v_r ; break ;
      case GC_VEC_STRING:  delete _data.v_s ; break ;

      case GC_VECTOR:      delete _data.v   ; break ;
      case GC_MAP:         delete _data.m   ; break ;
      default:
        break ;
    }
    _data_type = GC::GC_NOTYPE ;
  }

  //! Return a string representing the type of data stored
  char const *
  GenericContainer::get_type_name() const {
    return typeName[_data_type] ;
  }

  //! Assign a generic container `a` to the generic container.
  void
  GenericContainer::load( GenericContainer const & gc ) {
    this -> clear() ;
    switch (gc._data_type) {
      //case GC_NOTYPE:      this -> clear()                 ; break ;
      case GC_POINTER:     this -> set_pointer(gc._data.p) ; break ;
      case GC_BOOL:        this -> set_bool(gc._data.b)    ; break ;
      case GC_INT:         this -> set_int(gc._data.i)     ; break ;
      case GC_REAL:        this -> set_real(gc._data.r)    ; break ;
      case GC_STRING:      this -> set_string(*gc._data.s) ; break ;

      case GC_VEC_POINTER: this -> set_vec_pointer(*gc._data.v_p) ; break ;
      case GC_VEC_BOOL:    this -> set_vec_bool(*gc._data.v_b)    ; break ;
      case GC_VEC_INT:     this -> set_vec_int(*gc._data.v_i)     ; break ;
      case GC_VEC_REAL:    this -> set_vec_real(*gc._data.v_r)    ; break ;
      case GC_VEC_STRING:  this -> set_vec_string(*gc._data.v_s)  ; break ;

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
      default:
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
  GenericContainer::ck(char const who[], TypeAllowed tp) const {
    GC_ASSERT( tp == _data_type,
               who <<
               " bad data type\nexpect: " << typeName[tp] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
  }

  void
  GenericContainer::ck_or_set(char const who[], TypeAllowed tp) {
    if ( _data_type == GC_NOTYPE ) {
      _data_type = tp ;
    } else {
      GC_ASSERT( tp == _data_type,
                 who <<
                 " bad data type\nexpect: " << typeName[tp] <<
                 "\nbut data stored is of type: " << typeName[_data_type] ) ;
    }
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
    if ( _data_type != GC::GC_STRING ) {
      clear() ;
      _data_type = GC::GC_STRING ;
      _data.s    = new string_type ;
    }
  }

  void
  GenericContainer::allocate_vec_pointer( unsigned sz ) {
    if ( _data_type != GC::GC_VEC_POINTER ) {
      clear() ;
      _data_type = GC::GC_VEC_POINTER ;
      _data.v_p  = new vec_pointer_type() ;
      if ( sz > 0 ) _data.v_p -> resize( sz ) ;
    }
  }

  GenericContainer &
  GenericContainer::free_pointer() {
    GC_ASSERT( GC_POINTER == _data_type || GC_NOTYPE == _data_type,
               "free_pointer() bad data type\nexpect: " << typeName[GC_POINTER] <<
               "\nbut data stored is of type: " << typeName[_data_type] ) ;
    _data.p = nullptr ;
    _data_type = GC::GC_NOTYPE ;
    return *this ;
  }


  void
  GenericContainer::allocate_vec_bool( unsigned sz ) {
    if ( _data_type != GC::GC_VEC_BOOL ) {
      clear() ;
      _data_type = GC::GC_VEC_BOOL ;
      _data.v_b  = new vec_bool_type() ;
      if ( sz > 0 ) _data.v_b -> resize( sz ) ;
    }
  }

  void
  GenericContainer::allocate_vec_int( unsigned sz ) {
    if ( _data_type != GC::GC_VEC_INT ) {
      clear() ;
      _data_type = GC::GC_VEC_INT ;
      _data.v_i  = new vec_int_type() ;
      if ( sz > 0 ) _data.v_i -> resize( sz ) ;
    }
  }

  void
  GenericContainer::allocate_vec_real( unsigned sz ) {
    if ( _data_type != GC::GC_VEC_REAL ) {
      clear() ;
      _data_type = GC::GC_VEC_REAL ;
      _data.v_r  = new vec_real_type() ;
      if ( sz > 0 ) _data.v_r -> resize( sz ) ;
    }
  }

  void
  GenericContainer::allocate_vec_string( unsigned sz ) {
    if ( _data_type != GC::GC_VEC_STRING ) {
      clear() ;
      _data_type = GC::GC_VEC_STRING ;
      _data.v_s  = new vec_string_type() ;
      if ( sz > 0 ) _data.v_s -> resize( sz ) ;
    }
  }

  void
  GenericContainer::allocate_vector( unsigned sz ) {
    if ( _data_type != GC::GC_VECTOR ) {
      clear() ;
      _data_type = GC::GC_VECTOR ;
      _data.v    = new vector_type() ;
      if ( sz > 0 ) _data.v -> resize( sz ) ;
    }
  }

  void
  GenericContainer::allocate_map() {
    if ( _data_type != GC::GC_MAP ) {
      clear() ;
      _data_type = GC::GC_MAP ;
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
    _data_type = GC::GC_POINTER ;
    return (_data.p = value) ;
  }

  bool_type &
  GenericContainer::set_bool( bool_type value ) {
    clear() ;
    _data_type = GC::GC_BOOL ;
    return (_data.b = value) ;
  }

  int_type &
  GenericContainer::set_int( int_type value ) {
    clear() ;
    _data_type = GC::GC_INT ;
    return (_data.i = value) ;
  }

  real_type &
  GenericContainer::set_real( real_type value ) {
    clear() ;
    _data_type = GC::GC_REAL ;
    return (_data.r = value) ;
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
  //    ____      _
  //   / ___| ___| |_
  //  | |  _ / _ \ __|
  //  | |_| |  __/ |_
  //   \____|\___|\__|
  */

#if defined(_WIN32) || defined(_WIN64)
  void *
  GenericContainer::get_pvoid() const {
    ck("get_pvoid",GC_POINTER) ;
    return _data.p ;
  }

  void **
  GenericContainer::get_ppvoid() const {
    ck("get_ppvoid",GC_POINTER) ;
    return (void **)&_data.p ;
  }
#endif

  //! If data is boolean, integer or floating point return number, otherwise return `0`.
  real_type
  GenericContainer::get_number() const {
    switch (_data_type) {
      case GC_BOOL: return (_data.b?1:0) ;
      case GC_INT:  return _data.i ;
      case GC_REAL: return _data.r ;
      default:
        break ;
    }
    return 0 ;
  }

  real_type
  GenericContainer::get_number( unsigned i ) const {
    switch (_data_type) {
      case GC_VEC_BOOL: return (*_data.v_b)[i] ;
      case GC_VEC_INT:  return (*_data.v_i)[i] ;
      case GC_VEC_REAL: return (*_data.v_r)[i] ;
      case GC_VECTOR:   return (*_data.v)[i].get_number() ;
      default: break ; // to quiet warnings
    }
    return 0 ;
  }

  bool_type &
  GenericContainer::get_bool() {
    ck_or_set("get_bool",GC_BOOL) ;
    return _data.b ;
  }

  bool_type const &
  GenericContainer::get_bool() const {
    ck("get_bool",GC_BOOL) ;
    return _data.b ;
  }

  int_type &
  GenericContainer::get_int() {
    ck_or_set("get_int",GC_INT) ;
    return _data.i ;
  }

  int_type const &
  GenericContainer::get_int() const {
    ck("get_int",GC_INT) ;
    return _data.i ;
  }

  real_type &
  GenericContainer::get_real() {
    ck_or_set("get_real",GC_REAL) ;
    return _data.r ;
  }

  real_type const &
  GenericContainer::get_real() const {
    ck("get_real",GC_REAL) ;
    return _data.r ;
  }

  string_type &
  GenericContainer::get_string() {
    ck_or_set("get_string",GC_STRING) ;
    return *_data.s ;
  }

  string_type const &
  GenericContainer::get_string() const {
    ck("get_string",GC_STRING) ;
    return *_data.s ;
  }

  vector_type &
  GenericContainer::get_vector() {
    ck("get_vector",GC_VECTOR) ;
    return *_data.v ;
  }

  vector_type const &
  GenericContainer::get_vector() const {
    ck("get_vector",GC_VECTOR) ;
    return *_data.v ;
  }

  vec_pointer_type &
  GenericContainer::get_vec_pointer() {
    ck("get_vec_pointer",GC_VEC_POINTER) ;
    return *_data.v_p ;
  }

  vec_pointer_type const &
  GenericContainer::get_vec_pointer() const {
    ck("get_vec_pointer",GC_VEC_POINTER) ;
    return *_data.v_p ;
  }

  vec_bool_type &
  GenericContainer::get_vec_bool() {
    ck("get_vec_bool",GC_VEC_BOOL) ;
    return *_data.v_b ;
  }

  vec_bool_type const &
  GenericContainer::get_vec_bool() const {
    ck("get_vec_bool",GC_VEC_BOOL) ;
    return *_data.v_b ;
  }

  vec_int_type &
  GenericContainer::get_vec_int() {
    if ( _data_type == GC_NOTYPE   ) set_vec_int() ;
    if ( _data_type == GC_VEC_BOOL ) promote_to_vec_int() ;
    ck("get_vec_int",GC_VEC_INT) ;
    return *_data.v_i ;
  }

  vec_int_type const &
  GenericContainer::get_vec_int() const {
    ck("get_vec_int",GC_VEC_INT) ;
    return *_data.v_i ;
  }

  vec_real_type &
  GenericContainer::get_vec_real() {
    if ( _data_type == GC_NOTYPE   ) set_vec_int() ;
    if ( _data_type == GC_VEC_BOOL || _data_type == GC_VEC_INT ) promote_to_vec_real() ;
    ck("get_vec_real",GC_VEC_REAL) ;
    return *_data.v_r ;
  }

  vec_real_type const &
  GenericContainer::get_vec_real() const {
    ck("get_vec_real",GC_VEC_REAL) ;
    return *_data.v_r ;
  }

  vec_string_type &
  GenericContainer::get_vec_string() {
    ck("get_vec_string",GC_VEC_STRING) ;
    return *_data.v_s ;
  }

  vec_string_type const &
  GenericContainer::get_vec_string() const {
    ck("get_vec_string",GC_VEC_STRING) ;
    return *_data.v_s ;
  }

  map_type &
  GenericContainer::get_map() {
    ck("get_map",GC_MAP) ;
    return *_data.m ;
  }

  map_type const &
  GenericContainer::get_map() const {
    ck("get_map",GC_MAP) ;
    return *_data.m ;
  }

  bool
  GenericContainer::exists( std::string const & s ) const {
    if ( _data_type != GC_MAP ) return false ;
    map_type::iterator iv = (*_data.m).find(s) ;
    return iv != (*_data.m).end() ;
  }

  // --------------------------------------------------------------
  bool_type
  GenericContainer::get_bool( unsigned i ) {
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
  GenericContainer::get_bool( unsigned i ) const {
    ck("get_bool",GC_VEC_BOOL) ;
    GC_ASSERT( i < _data.v_b->size(), "get_bool( " << i << " ) const, out of range" ) ;
    return (*_data.v_b)[i] ;
  }

  int_type &
  GenericContainer::get_int( unsigned i ) {
    if      ( _data_type == GC_NOTYPE ) set_vec_int() ;
    else if ( _data_type == GC_VEC_BOOL ) promote_to_vec_int() ;
    if ( _data_type == GC_VEC_INT ) {
      CHECK_RESIZE(_data.v_i,i) ; // correct type, check size
      return (*_data.v_i)[i] ;
    } else {
      if ( _data_type != GC_VECTOR ) promote_to_vector() ;
      CHECK_RESIZE(_data.v,i) ;
      return (*_data.v)[i].set_int(0) ;
    }
  }

  int_type const &
  GenericContainer::get_int( unsigned i ) const {
    ck("get_int",GC_VEC_INT) ;
    GC_ASSERT( i < _data.v_i->size(), "get_int( " << i << " ) const, out of range" ) ;
    return (*_data.v_i)[i] ;
  }

  real_type &
  GenericContainer::get_real( unsigned i ) {
    if      ( _data_type == GC_NOTYPE ) set_vec_real() ;
    else if ( _data_type == GC_VEC_BOOL || _data_type == GC_VEC_INT ) promote_to_vec_real() ;
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
  GenericContainer::get_real( unsigned i ) const  {
    ck("get_real",GC_VEC_REAL) ;
    GC_ASSERT( i < _data.v_r->size(), "get_real( " << i << " ) const, out of range" ) ;
    return (*_data.v_r)[i] ;
  }

  string_type &
  GenericContainer::get_string( unsigned i ) {
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
  GenericContainer::get_string( unsigned i ) const {
    ck("get_string",GC_VEC_STRING) ;
    GC_ASSERT( i < _data.v_s->size(), "get_string( " << i << " ) const, out of range" ) ;
    return (*_data.v_s)[i] ;
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
  GenericContainer::info( std::ostream & stream ) const {
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
      case GC_INT:
        stream << "Integer: " << _data.i << '\n' ;
        break ;
      case GC_REAL:
        stream << "Floating Point: " << _data.r << '\n' ;
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
      case GC_VEC_INT:
        stream << "Vector of integer of size " << _data.v_i->size() << '\n' ;
        break ;
      case GC_VEC_REAL:
        stream << "Vector of floating point number of size " << _data.v_r->size() << '\n' ;
        break ;
      case GC_VEC_STRING:
        stream << "Vector of string of size " << _data.v_s->size() << '\n' ;
        break ;
      case GC_VECTOR:
        stream << "Vector of generic data type of size " << _data.v->size() << '\n' ;
        break ;
      case GC_MAP:
        stream << "Map\n" ;
        break ;
      default:
        stream << "Type N. " << _data_type << " not recognized\n" ;
        break ;
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
    ck("operator []",GC_VECTOR) ;
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
  GenericContainer::operator () ( unsigned i ) {
    ck("operator ()",GC_VECTOR) ;
    GC_ASSERT( i < _data.v->size(), "operator (), index " << i << " out of range" ) ;
    return (*_data.v)[i] ;
  }

  GenericContainer const &
  GenericContainer::operator () ( unsigned i ) const {
    ck("operator ()",GC_VECTOR) ;
    GC_ASSERT( i < _data.v->size(), "operator () const, index " << i << " out of range" ) ;
    return (*_data.v)[i] ;
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
    ck("operator []",GC_MAP) ;
    return (*_data.m)[s] ;
  }

  GenericContainer &
  GenericContainer::operator () ( std::string const & s ) {
    ck("operator ()",GC_MAP) ;
    map_type::iterator iv = (*_data.m) . find(s) ;
    GC_ASSERT( iv != (*_data.m) . end(), "operator(): Cannot find key '" << s << "'!" ) ;
    return iv -> second ;
  }

  GenericContainer const &
  GenericContainer::operator () ( std::string const & s ) const {
    ck("operator ()",GC_MAP) ;
    map_type::const_iterator iv = (*_data.m) . find(s) ;
    GC_ASSERT( iv != (*_data.m) . end(), "operator(): Cannot find key '" << s << "'!" ) ;
    return iv -> second ;
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
      case GC_INT:
        break ;
      default:
        GC_ASSERT( false, ":promote_to_int() cannot promote " << get_type_name() << " to int") ;
        break ;
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
      case GC_INT:
        set_real(_data.i) ;
        break ;
      case GC_REAL:
        break ;
      default:
        GC_ASSERT( false, ":promote_to_real() cannot promote " << get_type_name() << " to real") ;
        break ;
    }
    return *this ;
  }

  GenericContainer const &
  GenericContainer::promote_to_vec_int() {
    switch (_data_type) {
      case GC_NOTYPE:
      { set_vec_int(1) ; get_int(0) = 0 ; }
        break ;
      case GC_BOOL:
      { int_type tmp = _data.b?1:0 ; set_vec_int(1) ; get_int(0) = tmp ; }
        break ;
      case GC_INT:
      { int_type tmp = _data.i ; set_vec_int(1) ; get_int(0) = tmp ; }
        break ;
      case GC_VEC_BOOL:
      { unsigned sz = unsigned(_data.v_b->size()) ;
        vec_int_type tmp(sz) ;
        for ( unsigned i = 0 ; i < sz ; ++i ) tmp[i] = (*_data.v_b)[i] ? 1 : 0 ;
        set_vec_int(sz) ;
        for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v_i)[i] = tmp[i] ;
      }
        break ;
      case GC_VEC_INT:
        break ;
      default:
        GC_ASSERT( false, ":promote_to_vec_int() cannot promote " << get_type_name() << " to real") ;
        break ;
    }
    return *this ;

  }

  //! If data contains vector of booleans or integer it is promoted to a vector of real.
  GenericContainer const &
  GenericContainer::promote_to_vec_real() {
    switch (_data_type) {
      case GC_NOTYPE:
        { set_vec_real(1) ; get_real(0) = 0 ; }
        break ;
      case GC_BOOL:
        { real_type tmp = _data.b?1:0 ; set_vec_real(1) ; get_real(0) = tmp ; }
        break ;
      case GC_INT:
        { real_type tmp = _data.i ; set_vec_real(1) ; get_real(0) = tmp ; }
        break ;
      case GC_REAL:
        { real_type tmp = _data.r ; set_vec_real(1) ; get_real(0) = tmp ; }
        break ;
      case GC_VEC_BOOL:
        { unsigned sz = unsigned(_data.v_b->size()) ;
          vec_real_type tmp(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) tmp[i] = (*_data.v_b)[i] ? 1 : 0 ;
          set_vec_real(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v_r)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_INT:
        { unsigned sz = unsigned(_data.v_i->size()) ;
          vec_real_type tmp(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) tmp[i] = (*_data.v_i)[i] ;
          set_vec_real(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v_r)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_REAL:
        break ;
      default:
        GC_ASSERT( false, ":promote_to_vec_real() cannot promote " << get_type_name() << " to real") ;
        break ;
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
      case GC_INT:
        { set_vector(1) ; (*this)[0] = _data.i ; }
        break ;
      case GC_REAL:
        { set_vector(1) ; (*this)[0] = _data.r ; }
        break ;
      case GC_STRING:
        { set_vector(1) ; (*this)[0] = *_data.s ; }
        break ;
      case GC_VEC_POINTER:
        { vec_pointer_type tmp(_data.v_p->begin(), _data.v_p->end()); // range-based constructor
          unsigned sz = unsigned(tmp.size()) ;
          set_vector(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_BOOL:
        { vec_bool_type tmp(_data.v_b->begin(), _data.v_b->end()); // range-based constructor
          unsigned sz = unsigned(tmp.size()) ;
          set_vector(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_INT:
        { vec_int_type tmp(_data.v_i->begin(), _data.v_i->end()); // range-based constructor
          unsigned sz = unsigned(tmp.size()) ;
          set_vector(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_REAL:
        { vec_real_type tmp(_data.v_r->begin(), _data.v_r->end()); // range-based constructor
          unsigned sz = unsigned(tmp.size()) ;
          set_vector(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v)[i] = tmp[i] ;
        }
        break ;
      case GC_VEC_STRING:
        { vec_string_type tmp(_data.v_s->begin(), _data.v_s->end()); // range-based   constructor
          unsigned sz = unsigned(tmp.size()) ;
          set_vector(sz) ;
          for ( unsigned i = 0 ; i < sz ; ++i ) (*_data.v)[i] = tmp[i] ;
        }
        break ;
      case GC_VECTOR:
        break ;
      default:
        GC_ASSERT( false, ":promote_to_vector() cannot promote " << get_type_name() << " to real") ;
        break ;
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
  GenericContainer::print( std::ostream      & stream,
                           std::string const & prefix,
                           std::string const & indent ) const {

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
      case GC_INT:
        stream << prefix << this -> get_int() << '\n' ;
        break ;
      case GC_REAL:
        stream << prefix << this -> get_real() << '\n' ;
        break ;
      case GC_STRING:
        stream << prefix << "\"" << this -> get_string() << "\"\n" ;
        break ;
      case GC_VEC_POINTER:
        { vec_pointer_type const & v = this -> get_vec_pointer() ;
          for ( vec_pointer_type::size_type i = 0 ; i < v.size() ; ++i )
            stream << prefix << "vec_pointer(" << i << "): " << (unsigned long)v[i] << '\n' ;
        }
        break ;
      case GC_VEC_BOOL:
        { vec_bool_type const & v = this -> get_vec_bool() ;
          stream << prefix << "[" ;
          for ( vec_bool_type::size_type i = 0 ; i < v.size() ; ++i )
            stream << (v[i]?" true":" false")  ;
          stream << " ]\n" ;
        }
        break ;
      case GC_VEC_INT:
        { vec_int_type const & v = this -> get_vec_int() ;
          stream << prefix << "[" ;
          for ( vec_int_type::size_type i = 0 ; i < v.size() ; ++i ) stream << " " << v[i] ;
          stream << " ]\n" ;
        }
        break ;
      case GC_VEC_REAL:
        { vec_real_type const & v = this -> get_vec_real() ;
          stream << prefix << "[" ;
          for ( vec_real_type::size_type i = 0 ; i < v.size() ; ++i ) stream << " " << v[i] ;
          stream << " ]\n" ;
        }
        break ;
      case GC_VEC_STRING:
        { vec_string_type const & v = this -> get_vec_string() ;
          for ( vec_string_type::size_type i = 0 ; i < v.size() ; ++i )
            stream << prefix << i << ": \"" << v[i] << "\"\n" ;
        }
        break ;

      case GC_VECTOR:
        { vector_type const & v = this -> get_vector() ;
          for ( vector_type::size_type i = 0 ; i < v.size() ; ++i ) {
            if ( v[i].simple_data() ) {
              stream << prefix << i << ": " ;
              v[i].print(stream,"") ;
            } else {
              stream << prefix << i << ":\n" ;
              v[i].print(stream,prefix+indent) ;
            }
          }
        }
        break ;
      case GC_MAP:
        { map_type const & m = this -> get_map() ;
          for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
            // check formatting using pcre
            #ifdef GENERIC_CONTAINER_NO_PCRE
            if ( im->second.simple_data() ) {
              stream << prefix << im->first << ": " ;
              im->second.print(stream,"") ;
            } else {
              stream << prefix << im->first << ":\n" ;
              im->second.print(stream,prefix+indent) ;
            }
            #else
            pcre       *const& reCompiled = *reinterpret_cast<pcre*const*>(&_reCompiled)      ;
            pcre_extra *const& pcreExtra  = *reinterpret_cast<pcre_extra*const*>(&_pcreExtra) ;
            // num+"@"+"underline character"
            // Try to find the regex in aLineToMatch, and report results.
            int imatch[30];
            int pcreExecRet = pcre_exec(reCompiled,
                                        pcreExtra,
                                        im->first.c_str(),
                                        int(im->first.length()), // length of string
                                        0,                       // Start looking at this point
                                        0,                       // OPTIONS
                                        imatch,
                                        30);                     // Length of subStrVec

            if ( pcreExecRet == 4 ) {
              // extract match
              int m1 = imatch[3]-imatch[2] ; // # or ##
              int m2 = imatch[5]-imatch[4] ; // -,= etc
              int m3 = imatch[7]-imatch[6] ; // # or ##
              std::string header = im->first.substr((std::string::size_type)imatch[6],
                                                    (std::string::size_type)imatch[7]) ; // header
              // found formatting
              if ( im->second.simple_data() ) {
                stream << prefix << header << ": " ;
                im->second.print(stream,"") ;
              } else {
                if ( m1 > 1 ) stream << '\n' ; // double ## --> add nel line
                stream << prefix << header ;
                if ( m2 > 0 ) {
                  stream << '\n' << prefix ;
                  char fmt = im->first[(std::string::size_type)imatch[4]] ; // underline char
                  while ( m3-- > 0 ) stream << fmt ; // underline header
                } else {
                  stream << ':' ;
                }
                stream << '\n' ;
                im->second.print(stream,prefix+indent) ;
              }
            } else {
              std::string header = pcreExecRet == 3 ?
                                   im->first.substr((std::string::size_type)imatch[4],
                                                    (std::string::size_type)imatch[5]) :
                                   im->first ;
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

      default:
        GC_ASSERT(false,"Error, print(...) unknown type!\n") ;
        break ;
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
  GenericContainer::to_yaml( std::ostream & stream, std::string const & prefix ) const {
    switch (_data_type) {
        
      case GC_NOTYPE:
        stream << "Empty!\n" ;
        break ;
      case GC_BOOL:
        stream << (this -> get_bool()?"true":"false") << '\n' ;
        break ;
      case GC_INT:
        stream << this -> get_int() << '\n' ;
        break ;
      case GC_REAL:
        stream << this -> get_real() << '\n' ;
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
      case GC_VEC_INT:
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
        
      default:
        GC_ASSERT( false, "Error, print(...) unknown type!\n" ) ;
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
