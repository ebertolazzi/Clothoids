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
 \example example11.cc

 */

#include "GenericContainer.hh"
#include <iostream>
#include <fstream>

using namespace std;
using namespace GC;

int
main() {

  cout
    << "\n\n\n"
    << "***********************\n"
    << "      example N.11     \n"
    << "***********************\n\n";

  try {
    GC::GenericContainer gc, gc_res;
    gc.set_mat_real(2,2);
    gc.get_real_at(1,1) = 2;
    gc.get_real_at(0,1) = 3;
    cout << "Result:\n";
    gc.dump(cout);
    gc.info(cout);
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n';
  }
  catch (...) {
    cout << "Unknonwn error\n";
  }
  
  cout << "ALL DONE!\n\n\n\n";
}
