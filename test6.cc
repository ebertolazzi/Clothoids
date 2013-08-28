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

#include "GenericContainer.hh"
#include "GenericContainerCinterface.h"

using namespace std ;

#define CK( A ) ok = A ; ASSERT( ok == GENERIC_CONTAINER_OK, "Error = " << ok )

int
main() {

  int ok ;

  cout << "\n\n\n"
       << "***********************\n"
       << "       test N.6        \n"
       << "***********************\n\n" ;
  
  try {
  
    CK( GC_select( "generic_container" ) ) ;
    
    CK( GC_set_vector(10) ) ;

    CK( GC_set_vector_position(0) ) ;
    CK( GC_set_int( 1 ) ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(1) ) ;
    CK( GC_set_empty_vector_of_real() ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(2) ) ;
    CK( GC_set_map() ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(3) ) ;
    CK( GC_set_empty_vector_of_string() ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(4) ) ;
    CK( GC_set_real( 1.3 ) ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(5) ) ;
    CK( GC_set_string( "pippo" ) ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(6) ) ;
    CK( GC_set_map() ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(7) ) ;
    CK( GC_set_empty_vector() ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(8) ) ;
    CK( GC_set_bool(true) ) ;
    CK( GC_pop_head() ) ; // return to vector
    
    // return to element 1 of the vector
    CK( GC_set_vector_position(1) ) ;
      CK( GC_push_real(1.1) ) ;
      CK( GC_push_real(2.2) ) ;
      CK( GC_push_real(3.3) ) ;
      CK( GC_push_real(4.4) ) ;
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(2) ) ;

      CK( GC_set_map_position("pippo") ) ;
      CK( GC_set_int(13) ) ;
      CK( GC_pop_head() ) ;

      CK( GC_set_map_position("pluto") ) ;
      CK( GC_set_int(1) ) ;
      CK( GC_pop_head() ) ;

      CK( GC_set_map_position("paperino") ) ;
      CK( GC_set_int(3) ) ;
      CK( GC_pop_head() ) ;

      CK( GC_set_map_position("aaa") ) ;
      CK( GC_set_string("stringa") ) ;
      CK( GC_pop_head() ) ;

    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(3) ) ;
      CK( GC_push_string("string1") ) ;
      CK( GC_push_string("string2") ) ;
      CK( GC_push_string("string3") ) ;
      CK( GC_push_string("string4") ) ;
    CK( GC_pop_head() ) ; // return to vector
    
    CK( GC_set_vector_position(6) ) ;
    
    CK( GC_set_map_position("aaa") ) ;
    CK( GC_set_int(123) ) ;
    CK( GC_pop_head() ) ;
    
    CK( GC_set_map_position("bbb") ) ;
    CK( GC_set_real(3.4) ) ;
    CK( GC_pop_head() ) ;
    
    CK( GC_set_map_position("vector") ) ;
      CK( GC_set_empty_vector_of_int() ) ;
      CK( GC_push_int(12) ) ;
      CK( GC_push_int(10) ) ;
      CK( GC_push_int(1) ) ;
    CK( GC_pop_head() ) ;
    
    CK( GC_pop_head() ) ; // return to vector

    CK( GC_set_vector_position(7) ) ;
      CK( GC_set_empty_vector() ) ;
      CK( GC_push_int(123) ) ;
      CK( GC_push_real(3.14) ) ;
      CK( GC_push_string("nonna papera") ) ;

    //CK( GC_pop_head() ) ;

    GC_print() ;
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n'  ;
  }
  catch (...) {
    cout << "Unknonwn error\n" ;
  }

  cout << "ALL DONE!\n\n\n\n" ;
}
