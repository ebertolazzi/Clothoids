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
 \example example1.cc

 Simple usage of `GenericContainer` class.
 The class is used to store one element
 of the basic available types.

 */

#include "GenericContainer.hh"

using namespace std;
using namespace GC;

int
main() {
  cout
    << "\n\n\n"
    << "***********************\n"
    << "     example N.1       \n"
    << "***********************\n\n";
  // Simple example using simple data
  GenericContainer gc1 = 1;
  GenericContainer gc2 = 1.2;
  GenericContainer gc3 = true;
  GenericContainer gc4 = "pippo";
  GenericContainer gc5;

  cout << "GenericContainer simple usage\n";
  cout << "gc1: "; gc1.info(cout);
  cout << "gc2: "; gc2.info(cout);
  cout << "gc3: "; gc3.info(cout);
  cout << "gc4: "; gc4.info(cout);
  cout << "gc5: "; gc5.info(cout);

  cout << "ALL DONE!\n\n\n\n";
}
