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
 \example example10.cc

 Example of use of Lua interface. `GenericContainer` is used
 to call a lua function and read the results.

 */

#include "GenericContainerLuaInterface.hh"
#include <iostream>
#include <fstream>

using namespace std;
using namespace GC;

int
main() {

  cout
    << "\n\n\n"
    << "***********************\n"
    << "      example N.10     \n"
    << "***********************\n\n";

  try {
    LuaInterpreter lua;
    lua.do_file("test_call.lua");

    GC::GenericContainer gc, gc_res;
    gc["function"] = "pippo";
    GC::GenericContainer & vec = gc["args"];
    vec[0] = 12;
    vec[1] = 13;
    vec[2] = "aaa";
    GC::GenericContainer & map = vec[3];
    map["nonna"] = "papera";
    GC::GenericContainer & vec1 = map["abc123"];
    vec1[0] = 12.3;
    vec1[1] = "a string";
    vec1[2] = 1;
    //dump(cout);
    lua.call( gc, gc_res );
    cout << "Result:\n";
    gc_res.dump(cout);
  }
  catch ( std::exception & exc ) {
    cout << exc.what() << '\n';
  }
  catch (...) {
    cout << "Unknonwn error\n";
  }

  cout << "ALL DONE!\n\n\n\n";
}
