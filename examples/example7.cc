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

#include "GenericContainer.hh"
#include <iostream>
#include <fstream>

using namespace std;

int
main() {

  cout
    << "\n\n\n"
    << "***********************\n"
    << "      example N.7      \n"
    << "***********************\n\n";

  try {
    GC::GenericContainer gc;
    ifstream file("data_example.txt");
    if ( file.fail() ) throw std::runtime_error("file to open file");
    gc.readFormattedData( file, "#", "\t " );
    gc.dump(cout);
    cout << "\n\nData Read:\n";
    gc.writeFormattedData( cout, '\t' );
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n';
  }
  catch (...) {
    cout << "Unknonwn error\n";
  }

  cout << "ALL DONE!\n\n\n\n";
}
