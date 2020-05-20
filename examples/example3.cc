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
 \example example3.cc

 Complex usage of `GenericContainer` class.
 The class is used to store a generic vector, i.e.
 a vector whose elements are of type `GenericContainer`
 or basic types.
 An error is forced to show an error message.

 */

#include "GenericContainer.hh"

using namespace std;
using namespace GC;

int
main() {

  cout
    << "\n\n\n"
    << "***********************\n"
    << "      example N.3      \n"
    << "***********************\n\n";

  // Using complex data

  try {

    GenericContainer gc;
    cout << "gc: "; gc.info(cout);

    gc . set_vector();
    cout << "gc: "; gc.info(cout);

    GC::vector_type & v = gc.get_vector();
    v.resize(10);
    cout << "gc: "; gc.info(cout);

    // access using vector
    v[0] = 1;
    v[1] = 1.2;
    v[2] = true;
    v[3] = "pippo";

    // or using overloading
    gc[4] = 1;
    gc[5] = 1.2;
    gc[6] = true;
    gc[7] = "pippo";
    cout << "Print contents of gc:\n";
    gc.dump(cout);

    // issue an error!
    gc[15] = 1.2;
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n';
  }
  catch (...) {
    cout << "Unknonwn error\n";
  }
  cout << "ALL DONE!\n\n\n\n";
}
