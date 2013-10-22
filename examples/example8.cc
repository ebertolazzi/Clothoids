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
 \example example7.cc
 
 Example of I/O utilities. `GenericContainer` is used 
 to read and store a table from a file and to print
 a formatted table on a stream.
 */

#include "GenericContainerLuaInterface.hh"
#include <iostream>
#include <fstream>

using namespace std ;

int
main() {

  cout << "\n\n\n"
       << "***********************\n"
       << "      example N.8      \n"
       << "***********************\n\n" ;
  
  try {
    GenericContainer gc ;
    LuaInterpreter   lua ;

    lua.load("test.lua") ;
    lua.global_to_GC("DATA",gc) ;
    gc.print(cout) ;
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n'  ;
  }
  catch (...) {
    cout << "Unknonwn error\n" ;
  }
  
  cout << "ALL DONE!\n\n\n\n" ;
}
