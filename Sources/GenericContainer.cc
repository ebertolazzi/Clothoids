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

///
/// file: GenericContainer.cc
///

#include "GenericContainer.hh"
#include <iomanip>

static char const *typeName[] = {
  "NOTYPE",
  "POINTER",
  "BOOL",
  "INT",
  "REAL",
  "STRING",
  "VEC_POINTER",
  "VEC_BOOL",
  "VEC_INT",
  "VEC_REAL",
  "VEC_STRING",
  "VECTOR",
  "MAP"
} ;


// costruttore
GenericContainer::GenericContainer()
: data_type(GenericContainer::NOTYPE)
{} ;


// distruttore
void
GenericContainer::clear() {
  switch (data_type) {
    case STRING:      delete data.s   ; break ;

    case VEC_POINTER: delete data.v_p ; break ;
    case VEC_BOOL:    delete data.v_b ; break ;
    case VEC_INT:     delete data.v_i ; break ;
    case VEC_REAL:    delete data.v_r ; break ;
    case VEC_STRING:  delete data.v_s ; break ;

    case VECTOR:      delete data.v   ; break ;
    case MAP:         delete data.m   ; break ;
    default:
    break ;
  }
  data_type = GenericContainer::NOTYPE ;
}

//! Assign a generic container `a` to the generic container.
GenericContainer const &
GenericContainer::operator = ( GenericContainer const & gc ) {
  this -> clear() ;
  switch (gc.data_type) {
    case NOTYPE:      this -> clear()                ; break ;
    case POINTER:     this -> set_pointer(gc.data.p) ; break ;
    case BOOL:        this -> set_bool(gc.data.b)    ; break ;
    case INT:         this -> set_int(gc.data.i)     ; break ;
    case REAL:        this -> set_real(gc.data.r)    ; break ;
    case STRING:      this -> set_string(*gc.data.s) ; break ;
      
    case VEC_POINTER: this -> set_vec_pointer(*gc.data.v_p) ; break ;
    case VEC_BOOL:    this -> set_vec_bool(*gc.data.v_b)    ; break ;
    case VEC_INT:     this -> set_vec_int(*gc.data.v_i)     ; break ;
    case VEC_REAL:    this -> set_vec_real(*gc.data.v_r)    ; break ;
    case VEC_STRING:  this -> set_vec_string(*gc.data.v_s)  ; break ;
      
    case VECTOR:
    { unsigned N = unsigned(gc.data.v->size()) ;
      allocate_vector( N ) ;
      std::copy( gc.data.v->begin(),
                 gc.data.v->end(),
                 this->data.v->begin() ) ;
    }
      break ;
    case MAP:
    { allocate_map() ;
      this->data.m->insert( gc.data.m->begin(), gc.data.m->end() ) ;
    }
    default:
      break ;
  }
  return * this ;
}

void
GenericContainer::ck(char const who[], TypeAllowed tp) const {
  ASSERT( tp == data_type,
          "GenericContainer::" << who <<
          " bad data type\nexpect: " << typeName[tp] <<
          "\nbut data stored is of type: " << typeName[data_type] ) ;
}

void
GenericContainer::allocate_string() {
  if ( data_type != GenericContainer::STRING ) {
    clear() ;
    data_type = GenericContainer::STRING ;
    data.s    = new string_type ;
  }
}

void
GenericContainer::allocate_vec_pointer( unsigned sz ) {
  if ( data_type != GenericContainer::VEC_POINTER ) {
    clear() ;
    data_type = GenericContainer::VEC_POINTER ;
    data.v_p  = new vec_pointer_type() ;
    if ( sz > 0 ) data.v_p -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_bool( unsigned sz ) {
  if ( data_type != GenericContainer::VEC_BOOL ) {
    clear() ;
    data_type = GenericContainer::VEC_BOOL ;
    data.v_b  = new vec_bool_type() ;
    if ( sz > 0 ) data.v_b -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_int( unsigned sz ) {
  if ( data_type != GenericContainer::VEC_INT ) {
    clear() ;
    data_type = GenericContainer::VEC_INT ;
    data.v_i  = new vec_int_type() ;
    if ( sz > 0 ) data.v_i -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_real( unsigned sz ) {
  if ( data_type != GenericContainer::VEC_REAL ) {
    clear() ;
    data_type = GenericContainer::VEC_REAL ;
    data.v_r  = new vec_real_type() ;
    if ( sz > 0 ) data.v_r -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vec_string( unsigned sz ) {
  if ( data_type != GenericContainer::VEC_STRING ) {
    clear() ;
    data_type = GenericContainer::VEC_STRING ;
    data.v_s  = new vec_string_type() ;
    if ( sz > 0 ) data.v_s -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_vector( unsigned sz ) {
  if ( data_type != GenericContainer::VECTOR ) {
    clear() ;
    data_type = GenericContainer::VECTOR ;
    data.v    = new vector_type() ;
    if ( sz > 0 ) data.v -> resize( sz ) ;
  }
}

void
GenericContainer::allocate_map() {
  if ( data_type != GenericContainer::MAP ) {
    clear() ;
    data_type = GenericContainer::MAP ;
    data.m    = new map_type() ;
  }
}

GenericContainer::pointer_type &
GenericContainer::set_pointer( pointer_type value ) {
  clear() ;
  data_type = GenericContainer::POINTER ;
  return (data.p = value) ;
}

GenericContainer::bool_type &
GenericContainer::set_bool( bool_type value ) {
  clear() ;
  data_type = GenericContainer::BOOL ;
  return (data.b = value) ;
}

GenericContainer::int_type &
GenericContainer::set_int( int_type value ) {
  clear() ;
  data_type = GenericContainer::INT ;
  return (data.i = value) ;
}

GenericContainer::real_type &
GenericContainer::set_real( real_type value ) {
  clear() ;
  data_type = GenericContainer::REAL ;
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

//! If data is boolean, integer or floating point return number, otherwise return `0`.
GenericContainer::real_type
GenericContainer::get_number() const {
  switch (data_type) {
    case BOOL: return (data.b?1:0) ;
    case INT:  return data.i ;
    case REAL: return data.r ;
    default:
      break ;
  }
  return 0 ;
}

GenericContainer::pointer_type &
GenericContainer::get_pointer() {
  ck("get_pointer",POINTER) ;
  return data.p ;
}

GenericContainer::pointer_type const &
GenericContainer::get_pointer() const {
  ck("get_pointer",POINTER) ;
  return data.p ;
}

GenericContainer::bool_type &
GenericContainer::get_bool() {
  ck("get_bool",BOOL) ;
  return data.b ;
}

GenericContainer::bool_type const &
GenericContainer::get_bool() const {
  ck("get_bool",BOOL) ;
  return data.b ;
}

GenericContainer::int_type &
GenericContainer::get_int() {
  ck("get_int",INT) ;
  return data.i ;
}

GenericContainer::int_type const &
GenericContainer::get_int() const {
  ck("get_int",INT) ;
  return data.i ;
}

GenericContainer::real_type &
GenericContainer::get_real() {
  ck("get_real",REAL) ;
  return data.r ;
}

GenericContainer::real_type const &
GenericContainer::get_real() const {
  ck("get_real",REAL) ;
  return data.r ;
}

GenericContainer::string_type &
GenericContainer::get_string() {
  ck("get_string",STRING) ;
  return *data.s ;
}

GenericContainer::string_type const &
GenericContainer::get_string() const {
  ck("get_string",STRING) ;
  return *data.s ;
}

GenericContainer::vector_type &
GenericContainer::get_vector() {
  ck("get_vector",VECTOR) ;
  return *data.v ;
}

GenericContainer::vector_type const &
GenericContainer::get_vector() const {
  ck("get_vector",VECTOR) ;
  return *data.v ;
}

GenericContainer::vec_pointer_type &
GenericContainer::get_vec_pointer() {
  ck("get_vec_pointer",VEC_POINTER) ;
  return *data.v_p ;
}

GenericContainer::vec_pointer_type const &
GenericContainer::get_vec_pointer() const {
  ck("get_vec_pointer",VEC_POINTER) ;
  return *data.v_p ;
}

GenericContainer::vec_bool_type &
GenericContainer::get_vec_bool() {
  ck("get_vec_bool",VEC_BOOL) ;
  return *data.v_b ;
}

GenericContainer::vec_bool_type const &
GenericContainer::get_vec_bool() const {
  ck("get_vec_bool",VEC_BOOL) ;
  return *data.v_b ;
}

GenericContainer::vec_int_type &
GenericContainer::get_vec_int() {
  ck("get_vec_int",VEC_INT) ;
  return *data.v_i ;
}

GenericContainer::vec_int_type const &
GenericContainer::get_vec_int() const {
  ck("get_vec_int",VEC_INT) ;
  return *data.v_i ;
}

GenericContainer::vec_real_type &
GenericContainer::get_vec_real() {
  ck("get_vec_real",VEC_REAL) ;
  return *data.v_r ;
}

GenericContainer::vec_real_type const &
GenericContainer::get_vec_real() const {
  ck("get_vec_real",VEC_REAL) ;
  return *data.v_r ;
}

GenericContainer::vec_string_type &
GenericContainer::get_vec_string() {
  ck("get_vec_string",VEC_STRING) ;
  return *data.v_s ;
}

GenericContainer::vec_string_type const &
GenericContainer::get_vec_string() const {
  ck("get_vec_string",VEC_STRING) ;
  return *data.v_s ;
}

GenericContainer::map_type &
GenericContainer::get_map() {
  ck("get_map",MAP) ;
  return *data.m ;
}

GenericContainer::map_type const &
GenericContainer::get_map() const {
  ck("get_map",MAP) ;
  return *data.m ;
}

bool
GenericContainer::exists( std::string const & s ) const {
  ck("exists",MAP) ;
  map_type::iterator iv = (*data.m).find(s) ;
  return iv != (*data.m).end() ;
}

//! Print to stream the kind of data stored
GenericContainer const &
GenericContainer::info( std::ostream & stream ) const {
  switch ( data_type ) {
    case NOTYPE:
      stream << "No data stored\n" ;
      break ;
    case POINTER:
      stream << "Generic pointer: " << std::hex << (unsigned long)(data.p) << '\n' ;
      break ;
    case BOOL:
      stream << "Boolean: " << (data.b?"true":"false") << '\n' ;
      break ;
    case INT:
      stream << "Integer: " << data.i << '\n' ;
      break ;
    case REAL:
      stream << "Floating Point: " << data.r << '\n' ;
      break ;
    case STRING:
      stream << "String: " << *data.s << '\n' ;
      break ;
    case VEC_POINTER:
      stream << "Vector of generic pointer of size " << data.v_p->size() << '\n' ;
      break ;
    case VEC_BOOL:
      stream << "Vector of boolean of size " << data.v_b->size() << '\n' ;
      break ;
    case VEC_INT:
      stream << "Vector of integer of size " << data.v_i->size() << '\n' ;
      break ;
    case VEC_REAL:
      stream << "Vector of floating point number of size " << data.v_r->size() << '\n' ;
      break ;
    case VEC_STRING:
      stream << "Vector of string of size " << data.v_s->size() << '\n' ;
      break ;
    case VECTOR:
      stream << "Vector of generic data type of size " << data.v->size() << '\n' ;
      break ;
    case MAP:
      stream << "Map\n" ;
      break ;
  }
  return *this ;
}

// --------------------------------------------------------

GenericContainer &
GenericContainer::operator [] ( unsigned i ) {
  ck("operator []",VECTOR) ;
  ASSERT( i < data.v->size(),
          "GenericContainer::operator [ " << i << " ] out of range\n"<<
          "Vector size: " << data.v->size() ) ;
  return (*data.v)[i] ;
}

GenericContainer const &
GenericContainer::operator [] ( unsigned i ) const {
  ck("operator []",VECTOR) ;
  ASSERT( i < data.v->size(),
         "GenericContainer::operator [ " << i << " ] out of range\n"<<
         "Vector size: " << data.v->size() ) ;
  return (*data.v)[i] ;
}

GenericContainer &
GenericContainer::operator () ( unsigned i ) {
  ck("operator ()",VECTOR) ;
  return (*data.v)[i] ;
}

GenericContainer const &
GenericContainer::operator () ( unsigned i ) const {
  ck("operator ()",VECTOR) ;
  return (*data.v)[i] ;
}

// --------------------------------------------------------

GenericContainer &
GenericContainer::operator [] ( std::string const & s ) {
  ck("operator []",MAP) ;
  return (*data.m)[s] ;
}

GenericContainer const &
GenericContainer::operator () ( std::string const & s ) const {
  ck("operator []",MAP) ;
  map_type::const_iterator iv = (*data.m) . find(s) ;
  ASSERT( iv != (*data.m) . end(), "GenericContainer, operator['" << s << "'] Cant find !" ) ;
  return iv -> second ;
}

// --------------------------------------------------------

void
GenericContainer::print( std::ostream & stream, std::string const & prefix ) const {
  switch (data_type) {

    case NOTYPE:
      stream << prefix << "Empty!\n" ;
      break ;
    case POINTER:
      stream << prefix << "pointer: " << (unsigned long)this -> get_pointer() << '\n' ;
      break ;
    case BOOL:
      stream << prefix << "bool: " << (this -> get_bool()?"true":"false") << '\n' ;
      break ;
    case INT:
      stream << prefix << "int: " << this -> get_int() << '\n' ;
      break ;
    case REAL:
      stream << prefix << "real: " << this -> get_real() << '\n' ;
      break ;
    case STRING:
      stream << prefix << "string: " << this -> get_string() << '\n' ;
      break ;

    case VEC_POINTER:
      { vec_pointer_type const & v = this -> get_vec_pointer() ;
        for ( int i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_pointer(" << i << "): " << (unsigned long)v[i] << '\n' ;
      }
      break ;
    case VEC_BOOL:
      { vec_bool_type const & v = this -> get_vec_bool() ;
        for ( int i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_bool(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case VEC_INT:
      { vec_int_type const & v = this -> get_vec_int() ;
        for ( int i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_int(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case VEC_REAL:
      { vec_real_type const & v = this -> get_vec_real() ;
        for ( int i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_real(" << i << "): " << v[i] << '\n' ;
      }
      break ;
    case VEC_STRING:
      { vec_string_type const & v = this -> get_vec_string() ;
        for ( int i = 0 ; i < v.size() ; ++i )
          stream << prefix << "vec_string(" << i << "): " << v[i] << "\n" ;
      }
      break ;

    case VECTOR:
      { vector_type const & v = this -> get_vector() ;
        for ( int i = 0 ; i < v.size() ; ++i ) {
          stream << prefix << "vec(" << i << "):\n" ;
          v[i].print(stream,prefix+"   ") ;
        }
      }
      break ;
    case MAP:
      { map_type const & m = this -> get_map() ;
        for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
          stream << prefix << "map['" << im->first << "']:\n" ;
          im->second.print(stream,prefix+"   ") ;
        }
      }
      break ;

    default:
      ASSERT(false,"Error, GenericContainer::print(...) unknown type!\n") ;
      break ;
  }
}

void
GenericContainer::to_yaml( std::ostream & stream, std::string const & prefix ) const {
  switch (data_type) {
      
    case NOTYPE:
      stream << "Empty!\n" ;
      break ;
    case BOOL:
      stream << (this -> get_bool()?"true":"false") << '\n' ;
      break ;
    case INT:
      stream << this -> get_int() << '\n' ;
      break ;
    case REAL:
      stream << this -> get_real() << '\n' ;
      break ;
    case STRING:
      stream << "'" << this -> get_string() << "'\n" ;
      break ;
      
    case VEC_BOOL:
      { vec_bool_type const & v = this -> get_vec_bool() ;
        stream << "[ " << (v[0]?"true":"false") ;
        for ( int i = 1 ; i < v.size() ; ++i )
          stream << ", " << (v[i]?"true":"false") ;
        stream << " ]\n" ;
      }
      break ;
    case VEC_INT:
      { vec_int_type const & v = this -> get_vec_int() ;
        stream << "[ " << v[0] ;
        for ( int i = 1 ; i < v.size() ; ++i )
          stream << ", " << v[i] ;
        stream << " ]\n" ;
      }
      break ;
    case VEC_REAL:
      { vec_real_type const & v = this -> get_vec_real() ;
        stream << "[ " << v[0] ;
        for ( int i = 1 ; i < v.size() ; ++i )
          stream << ", " << v[i] ;
        stream << " ]\n" ;
      }
      break ;
    case VEC_STRING:
      { vec_string_type const & v = this -> get_vec_string() ;
        stream << "[ '" << v[0] << "'" ;
        for ( int i = 1 ; i < v.size() ; ++i )
          stream << ", '" << v[i] << "'" ;
        stream << " ]\n" ;
      }
      break ;

    case VECTOR:
      { vector_type const & v = this -> get_vector() ;
        stream << '\n' ;
        for ( int i = 0 ; i < v.size() ; ++i ) {
          stream << prefix << "- " ;
          v[i].to_yaml(stream,prefix+"  ") ;
        }
      }
      break ;
    case MAP:
      { map_type const & m = this -> get_map() ;
        stream << '\n' ;
        for ( map_type::const_iterator im = m.begin() ; im != m.end() ; ++im ) {
          stream << prefix << im->first << ": " ;
          im->second.to_yaml(stream,prefix+"  ") ;
        }
      }
      break ;
      
    default:
      ASSERT( false, "Error, GenericContainer::print(...) unknown type!\n" ) ;
      break ;
  }
}

///
/// eof: GenericContainer.cc
///
