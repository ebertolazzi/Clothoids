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

// costruttore
GenericContainer::GenericContainer()
: data_type(GenericContainer::GC_NOTYPE)
{}

// distruttore
void
GenericContainer::clear() {
  switch (data_type) {
    case GC_POINTER:
      GC_WARNING( data.p == nullptr, "In deleting GenericContainer find a pointer not deallocated!") ;
      break ;
    case GC_STRING:      delete data.s   ; break ;

    case GC_VEC_POINTER: delete data.v_p ; break ;
    case GC_VEC_BOOL:    delete data.v_b ; break ;
    case GC_VEC_INT:     delete data.v_i ; break ;
    case GC_VEC_REAL:    delete data.v_r ; break ;
    case GC_VEC_STRING:  delete data.v_s ; break ;

    case GC_VECTOR:      delete data.v   ; break ;
    case GC_MAP:         delete data.m   ; break ;
    default:
    break ;
  }
  data_type = GenericContainer::GC_NOTYPE ;
}

//! Return a string representing the type of data stored
std::string
GenericContainer::get_type_name() const {
  return typeName[data_type] ;
}

//! Assign a generic container `a` to the generic container.
GenericContainer const &
GenericContainer::operator = ( GenericContainer const & gc ) {
  this -> clear() ;
  switch (gc.data_type) {
    case GC_NOTYPE:      this -> clear()                ; break ;
    case GC_POINTER:     this -> set_pointer(gc.data.p) ; break ;
    case GC_BOOL:        this -> set_bool(gc.data.b)    ; break ;
    case GC_INT:         this -> set_int(gc.data.i)     ; break ;
    case GC_REAL:        this -> set_real(gc.data.r)    ; break ;
    case GC_STRING:      this -> set_string(*gc.data.s) ; break ;
      
    case GC_VEC_POINTER: this -> set_vec_pointer(*gc.data.v_p) ; break ;
    case GC_VEC_BOOL:    this -> set_vec_bool(*gc.data.v_b)    ; break ;
    case GC_VEC_INT:     this -> set_vec_int(*gc.data.v_i)     ; break ;
    case GC_VEC_REAL:    this -> set_vec_real(*gc.data.v_r)    ; break ;
    case GC_VEC_STRING:  this -> set_vec_string(*gc.data.v_s)  ; break ;
      
    case GC_VECTOR:
      { unsigned N = unsigned(gc.data.v->size()) ;
        allocate_vector( N ) ;
        std::copy( gc.data.v->begin(),
                   gc.data.v->end(),
                   this->data.v->begin() ) ;
      }
    break ;
    case GC_MAP:
      { allocate_map() ;
        this->data.m->insert( gc.data.m->begin(), gc.data.m->end() ) ;
      }
    break ;
    default:
    break ;
  }
  return * this ;
}

int
GenericContainer::ck( TypeAllowed tp ) const {
  if ( tp == data_type ) return 0 ; // ok
  if ( tp == GC_NOTYPE ) return 1 ; //
  return 2 ;
}

void
GenericContainer::ck(char const who[], TypeAllowed tp) const {
  GC_ASSERT( tp == data_type,
             who <<
             " bad data type\nexpect: " << typeName[tp] <<
             "\nbut data stored is of type: " << typeName[data_type] ) ;
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
  if ( data_type != GenericContainer::GC_STRING ) {
    clear() ;
    data_type = GenericContainer::GC_STRING ;
    data.s    = new string_type ;
  }
}

void
GenericContainer::allocate_vec_pointer( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VEC_POINTER ) {
    clear() ;
    data_type = GenericContainer::GC_VEC_POINTER ;
    data.v_p  = new vec_pointer_type() ;
    if ( sz > 0 ) data.v_p -> resize( sz ) ;
  }
}

GenericContainer &
GenericContainer::free_pointer() {
  GC_ASSERT( GC_POINTER == data_type || GC_NOTYPE == data_type,
             "free_pointer() bad data type\nexpect: " << typeName[GC_POINTER] <<
             "\nbut data stored is of type: " << typeName[data_type] ) ;
  data.p = nullptr ;
  data_type = GenericContainer::GC_NOTYPE ;
  return *this ;
}


void
GenericContainer::allocate_vec_bool( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VEC_BOOL ) {
    clear() ;
    data_type = GenericContainer::GC_VEC_BOOL ;
    data.v_b  = new vec_bool_type() ;
    if ( sz > 0 ) data.v_b -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_int( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VEC_INT ) {
    clear() ;
    data_type = GenericContainer::GC_VEC_INT ;
    data.v_i  = new vec_int_type() ;
    if ( sz > 0 ) data.v_i -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_real( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VEC_REAL ) {
    clear() ;
    data_type = GenericContainer::GC_VEC_REAL ;
    data.v_r  = new vec_real_type() ;
    if ( sz > 0 ) data.v_r -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_string( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VEC_STRING ) {
    clear() ;
    data_type = GenericContainer::GC_VEC_STRING ;
    data.v_s  = new vec_string_type() ;
    if ( sz > 0 ) data.v_s -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vector( unsigned sz ) {
  if ( data_type != GenericContainer::GC_VECTOR ) {
    clear() ;
    data_type = GenericContainer::GC_VECTOR ;
    data.v    = new vector_type() ;
    if ( sz > 0 ) data.v -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_map() {
  if ( data_type != GenericContainer::GC_MAP ) {
    clear() ;
    data_type = GenericContainer::GC_MAP ;
    data.m    = new map_type() ;
  }
}

/*
//   ____       _
//  / ___|  ___| |_
//  \___ \ / _ \ __|
//   ___) |  __/ |_
//  |____/ \___|\__|
*/

GenericContainer::pointer_type &
GenericContainer::set_pointer( pointer_type value ) {
  clear() ;
  data_type = GenericContainer::GC_POINTER ;
  return (data.p = value) ;
}

GenericContainer::bool_type &
GenericContainer::set_bool( bool_type value ) {
  clear() ;
  data_type = GenericContainer::GC_BOOL ;
  return (data.b = value) ;
}

GenericContainer::int_type &
GenericContainer::set_int( int_type value ) {
  clear() ;
  data_type = GenericContainer::GC_INT ;
  return (data.i = value) ;
}

GenericContainer::real_type &
GenericContainer::set_real( real_type value ) {
  clear() ;
  data_type = GenericContainer::GC_REAL ;
  return (data.r = value) ;
}

GenericContainer::string_type &
GenericContainer::set_string( string_type const & value ) {
  allocate_string() ;
  return (*data.s = value) ;
}

GenericContainer::vec_pointer_type &
GenericContainer::set_vec_pointer( unsigned sz ) {
  allocate_vec_pointer( sz ) ;
  return *data.v_p ;
}

GenericContainer::vec_pointer_type &
GenericContainer::set_vec_pointer( vec_pointer_type const & v ) {
  allocate_vec_pointer( unsigned(v.size()) ) ;
  std::copy( v.begin(), v.end(), data.v_p->begin() ) ;
  return *data.v_p ;
}

GenericContainer::vec_bool_type &
GenericContainer::set_vec_bool( unsigned sz ) {
  allocate_vec_bool( sz ) ; return *data.v_b ;
}

GenericContainer::vec_bool_type &
GenericContainer::set_vec_bool( vec_bool_type const & v ) {
  allocate_vec_bool( unsigned(v.size()) ) ;
  std::copy( v.begin(), v.end(), data.v_b->begin() ) ;
  return *data.v_b ;
}

GenericContainer::vec_int_type &
GenericContainer::set_vec_int( unsigned sz ) {
  allocate_vec_int( sz ) ;
  return *data.v_i ;
}

GenericContainer::vec_int_type &
GenericContainer::set_vec_int( vec_int_type const & v ) {
  allocate_vec_int( unsigned(v.size()) ) ;
  std::copy( v.begin(), v.end(), data.v_i->begin() ) ;
  return *data.v_i ;
}

GenericContainer::vec_real_type &
GenericContainer::set_vec_real( unsigned sz ) {
  allocate_vec_real( sz ) ;
  return *data.v_r ;
}

GenericContainer::vec_real_type &
GenericContainer::set_vec_real( vec_real_type const & v ) {
  allocate_vec_real( unsigned(v.size()) ) ;
  std::copy( v.begin(), v.end(), data.v_r->begin() ) ;
  return *data.v_r ;
}

GenericContainer::vec_string_type &
GenericContainer::set_vec_string( unsigned sz ) {
  allocate_vec_string( sz ) ;
  return *data.v_s ;
}

GenericContainer::vec_string_type &
GenericContainer::set_vec_string( vec_string_type const & v ) {
  allocate_vec_string( unsigned(v.size()) ) ;
  std::copy( v.begin(), v.end(), data.v_s->begin() ) ;
  return *data.v_s ;
}

GenericContainer::vector_type &
GenericContainer::set_vector( unsigned sz ) {
  allocate_vector( sz ) ;
  return *data.v ;
}

GenericContainer::map_type &
GenericContainer::set_map() {
  allocate_map() ;
  return *data.m ;
}

/*
//    ____      _
//   / ___| ___| |_
//  | |  _ / _ \ __|
//  | |_| |  __/ |_
//   \____|\___|\__|
*/

//! If data is boolean, integer or floating point return number, otherwise return `0`.
GenericContainer::real_type
GenericContainer::get_number() const {
  switch (data_type) {
    case GC_BOOL: return (data.b?1:0) ;
    case GC_INT:  return data.i ;
    case GC_REAL: return data.r ;
    default:
      break ;
  }
  return 0 ;
}

GenericContainer::real_type
GenericContainer::get_number( unsigned i ) const {
  switch (data_type) {
    case GC_VEC_BOOL: return (*data.v_b)[i] ;
    case GC_VEC_INT:  return (*data.v_i)[i] ;
    case GC_VEC_REAL: return (*data.v_r)[i] ;
    case GC_VECTOR:   return (*data.v)[i].get_number() ;
    default: break ; // to quiet warnings
  }
  return 0 ;
}

GenericContainer::bool_type &
GenericContainer::get_bool() {
  ck("get_bool",GC_BOOL) ;
  return data.b ;
}

GenericContainer::bool_type const &
GenericContainer::get_bool() const {
  ck("get_bool",GC_BOOL) ;
  return data.b ;
}

GenericContainer::int_type &
GenericContainer::get_int() {
  ck("get_int",GC_INT) ;
  return data.i ;
}

GenericContainer::int_type const &
GenericContainer::get_int() const {
  ck("get_int",GC_INT) ;
  return data.i ;
}

GenericContainer::real_type &
GenericContainer::get_real() {
  ck("get_real",GC_REAL) ;
  return data.r ;
}

GenericContainer::real_type const &
GenericContainer::get_real() const {
  ck("get_real",GC_REAL) ;
  return data.r ;
}

GenericContainer::string_type &
GenericContainer::get_string() {
  ck("get_string",GC_STRING) ;
  return *data.s ;
}

GenericContainer::string_type const &
GenericContainer::get_string() const {
  ck("get_string",GC_STRING) ;
  return *data.s ;
}

GenericContainer::vector_type &
GenericContainer::get_vector() {
  ck("get_vector",GC_VECTOR) ;
  return *data.v ;
}

GenericContainer::vector_type const &
GenericContainer::get_vector() const {
  ck("get_vector",GC_VECTOR) ;
  return *data.v ;
}

GenericContainer::vec_pointer_type &
GenericContainer::get_vec_pointer() {
  ck("get_vec_pointer",GC_VEC_POINTER) ;
  return *data.v_p ;
}

GenericContainer::vec_pointer_type const &
GenericContainer::get_vec_pointer() const {
  ck("get_vec_pointer",GC_VEC_POINTER) ;
  return *data.v_p ;
}

GenericContainer::vec_bool_type &
GenericContainer::get_vec_bool() {
  ck("get_vec_bool",GC_VEC_BOOL) ;
  return *data.v_b ;
}

GenericContainer::vec_bool_type const &
GenericContainer::get_vec_bool() const {
  ck("get_vec_bool",GC_VEC_BOOL) ;
  return *data.v_b ;
}

GenericContainer::vec_int_type &
GenericContainer::get_vec_int() {
  ck("get_vec_int",GC_VEC_INT) ;
  return *data.v_i ;
}

GenericContainer::vec_int_type const &
GenericContainer::get_vec_int() const {
  ck("get_vec_int",GC_VEC_INT) ;
  return *data.v_i ;
}

GenericContainer::vec_real_type &
GenericContainer::get_vec_real() {
  ck("get_vec_real",GC_VEC_REAL) ;
  return *data.v_r ;
}

GenericContainer::vec_real_type const &
GenericContainer::get_vec_real() const {
  ck("get_vec_real",GC_VEC_REAL) ;
  return *data.v_r ;
}

GenericContainer::vec_string_type &
GenericContainer::get_vec_string() {
  ck("get_vec_string",GC_VEC_STRING) ;
  return *data.v_s ;
}

GenericContainer::vec_string_type const &
GenericContainer::get_vec_string() const {
  ck("get_vec_string",GC_VEC_STRING) ;
  return *data.v_s ;
}

GenericContainer::map_type &
GenericContainer::get_map() {
  ck("get_map",GC_MAP) ;
  return *data.m ;
}

GenericContainer::map_type const &
GenericContainer::get_map() const {
  ck("get_map",GC_MAP) ;
  return *data.m ;
}

bool
GenericContainer::exists( std::string const & s ) const {
  ck("exists",GC_MAP) ;
  map_type::iterator iv = (*data.m).find(s) ;
  return iv != (*data.m).end() ;
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
  switch ( data_type ) {
    case GC_NOTYPE:
      stream << "GenericContainer: No data stored\n" ;
      break ;
    case GC_POINTER:
      stream << "Generic pointer: " << std::hex << (unsigned long)(data.p) << '\n' ;
      break ;
    case GC_BOOL:
      stream << "Boolean: " << (data.b?"true":"false") << '\n' ;
      break ;
    case GC_INT:
      stream << "Integer: " << data.i << '\n' ;
      break ;
    case GC_REAL:
      stream << "Floating Point: " << data.r << '\n' ;
      break ;
    case GC_STRING:
      stream << "String: " << *data.s << '\n' ;
      break ;
    case GC_VEC_POINTER:
      stream << "Vector of generic pointer of size " << data.v_p->size() << '\n' ;
      break ;
    case GC_VEC_BOOL:
      stream << "Vector of boolean of size " << data.v_b->size() << '\n' ;
      break ;
    case GC_VEC_INT:
      stream << "Vector of integer of size " << data.v_i->size() << '\n' ;
      break ;
    case GC_VEC_REAL:
      stream << "Vector of floating point number of size " << data.v_r->size() << '\n' ;
      break ;
    case GC_VEC_STRING:
      stream << "Vector of string of size " << data.v_s->size() << '\n' ;
      break ;
    case GC_VECTOR:
      stream << "Vector of generic data type of size " << data.v->size() << '\n' ;
      break ;
    case GC_MAP:
      stream << "Map\n" ;
      break ;
    default:
      stream << "Type N. " << data_type << " not recognized\n" ;
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
  ck("operator []",GC_VECTOR) ;
  GC_ASSERT( i < data.v->size(),
             "operator [ " << i << " ] out of range\n"<<
             "Vector size: " << data.v->size() ) ;
  return (*data.v)[i] ;
}

GenericContainer const &
GenericContainer::operator [] ( unsigned i ) const {
  ck("operator []",GC_VECTOR) ;
  GC_ASSERT( i < data.v->size(),
             "operator [ " << i << " ] out of range\n"<<
             "Vector size: " << data.v->size() ) ;
  return (*data.v)[i] ;
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
  return (*data.v)[i] ;
}

GenericContainer const &
GenericContainer::operator () ( unsigned i ) const {
  ck("operator ()",GC_VECTOR) ;
  return (*data.v)[i] ;
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
  switch ( ck( GC_MAP ) ) {
    case 0: break ; // data present
    default: set_map() ; // data must be allocated ;
  }
  if ( ck( GC_MAP ) != 0 ) set_map() ; // if not data present allocate!
  return (*data.m)[s] ;
}

GenericContainer const &
GenericContainer::operator [] ( std::string const & s ) const {
  ck("operator []",GC_MAP) ;
  return (*data.m)[s] ;
}

GenericContainer const &
GenericContainer::operator () ( std::string const & s ) const {
  ck("operator ()",GC_MAP) ;
  map_type::const_iterator iv = (*data.m) . find(s) ;
  GC_ASSERT( iv != (*data.m) . end(), "operator['" << s << "'] Cant find !" ) ;
  return iv -> second ;
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
GenericContainer::print( std::ostream & stream, std::string const & prefix ) const {
  switch (data_type) {

    case GC_NOTYPE:
      stream << prefix << "Empty!\n" ;
      break ;
    case GC_POINTER:
      stream << prefix << "pointer: " << (unsigned long)get_pointer<void*>() << '\n' ;
      break ;
    case GC_BOOL:
      stream << prefix << "bool: " << (this -> get_bool()?"true":"false") << '\n' ;
      break ;
    case GC_INT:
      stream << prefix << "int: " << this -> get_int() << '\n' ;
      break ;
    case GC_REAL:
      stream << prefix << "real: " << this -> get_real() << '\n' ;
      break ;
    case GC_STRING:
      stream << prefix << "string: " << this -> get_string() << '\n' ;
      break ;

    case GC_VEC_POINTER:
      { vec_pointer_type const & v = this -> get_vec_pointer() ;
        for ( vec_pointer_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_pointer(" << i << "): " << (unsigned long)v[i] << '\n' ;
      }
      break ;
    case GC_VEC_BOOL:
      { vec_bool_type const & v = this -> get_vec_bool() ;
        for ( vec_bool_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_bool(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case GC_VEC_INT:
      { vec_int_type const & v = this -> get_vec_int() ;
        for ( vec_int_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_int(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case GC_VEC_REAL:
      { vec_real_type const & v = this -> get_vec_real() ;
        for ( vec_real_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_real(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case GC_VEC_STRING:
      { vec_string_type const & v = this -> get_vec_string() ;
        for ( vec_string_type::size_type i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_string(" << i << "): " << v[i] << "\n" ;
      }
      break ;

    case GC_VECTOR:
      { vector_type const & v = this -> get_vector() ;
        for ( vector_type::size_type i = 0 ; i < v.size() ; ++i ) {
          stream << prefix << "vec(" << i << "):\n" ;
          v[i].print(stream,prefix+"   ") ;
        }
      }
      break ;
    case GC_MAP:
      { map_type const & m = this -> get_map() ;
        for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
          stream << prefix << "map['" << im->first << "']:\n" ;
          im->second.print(stream,prefix+"   ") ;
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
  switch (data_type) {
      
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

//
// eof: GenericContainer.cc
//
