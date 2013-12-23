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

/*!
 \example example4.cc
 
 Complex usage of `GenericContainer` class.
 The class is used to store a generic map, i.e.
 a map whose elements are of type `GenericContainer`
 or basic types.
 
 */

#include "GenericContainer.hh"

using namespace std ;
using namespace GC ;

int
main() {

  cout << "\n\n\n"
       << "***********************\n"
       << "      example N.4      \n"
       << "***********************\n\n" ;

  // Using complex data
  
  try {
    
    GenericContainer gc ;
    cout << "gc: " ; gc.info(cout) ;
    
    gc . set_map() ;
    cout << "gc: " ; gc.info(cout) ;
    
    GC::map_type & m = gc . set_map() ;
    cout << "gc: " ; gc.info(cout) ;

    // access using map and vector like syntax 
    m["a"] = 1 ;
    m["b"] = 1.2 ;
    m["c"] = true ;
    m["d"] = "pippo" ;
    
    // or using overloading
    gc["e"] = 1 ;
    gc["f"] = 1.2 ;
    gc["g"] = true ;
    gc["h"] = "pippo" ;
    gc["pointer"] = &gc ;
    cout << "Print contents of gc:\n" ;
    gc.print(cout) ;

    // issue an error!
    gc[0] = 1.2 ;
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n'  ;
  }
  catch (...) {
    cout << "Unknonwn error\n" ;
  }
  cout << "ALL DONE!\n\n\n\n" ;
}
