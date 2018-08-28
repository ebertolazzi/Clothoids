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
// file: GenericContainer_doc.hh
//

/*! 
\mainpage  Generic container class
\author    Enrico Bertolazzi (enrico.bertolazzi@unitn.it), homepage: http://www.ing.unitn.it/~bertolaz
\version   1.0.5
\date      2013
\copyright GNU Public License.

\details
 
This library available at
 
 - https://github.com/ebertolazzi/GenericContainer

implement `GenericContainer` a class which permit to store eterogeneous data:

- pointer
- boolean
- integer
- floating point
- complex floating point
- string
- vector of pointer
- vector of boolean
- vector of integer
- vector of floating point
- vector of complex floating point
- vector of string

in addition to this data type the following two container are added

- vector of `GenericContainer`
- map of `GenericContainer`

this permits to build complex recursive data.
The main usage of the class is in interchange data with
scripting language like `Ruby`, `Lua`, `MATLAB`.


The usage is simple, for example it
can be used as an associative array with eterogenous data

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc;
gc["one"]  = 1;             // store integer
gc["two"]  = true;          // store a boolean
gc["3"]    = 1.4;           // store floating point number
gc["four"] = "pippo";       // store a string
gc["five"].set_vec_int(10); // store a vector of integer of 10 elements
~~~~~~~~~~~~~

and to retrieve elements

~~~~~~~~~~~~~{.cc}
cout << gc["one"].get_int()     << '\n';
cout << gc["two"].get_bool()    << '\n';
cout << gc["3"].get_real()      << '\n';
cout << gc["four"].get_string() << '\n';
GC::vec_int_type & v = gc["five"].get_vec_int();
cout << v[1] << '\n';
~~~~~~~~~~~~~

For more complex examples and recursive data see examples files
in the distribution.

============================

\section sec1 Initialization
 
Getting an instance of `GenericContainer`
 
~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc; // initialize empty container
~~~~~~~~~~~~~

if can be initialized to a boolean
~~~~~~~~~~~~~{.cc}
gc.set_bool(true);
gc.set_bool(false);
~~~~~~~~~~~~~

to an integer

~~~~~~~~~~~~~{.cc}
gc.set_int(123);
~~~~~~~~~~~~~

to a floating point number

~~~~~~~~~~~~~{.cc}
gc.set_real(1.23);
gc.set_real(3);
~~~~~~~~~~~~~

to a string

~~~~~~~~~~~~~{.cc}
gc.set_string("a C string");
string s = "a C++ sring";
gc.set_string(s);
~~~~~~~~~~~~~

to a pointer

~~~~~~~~~~~~~{.cc}
gc.set_pointer(&cout);
~~~~~~~~~~~~~

to a vector of boolean, integer or floating points

~~~~~~~~~~~~~{.cc}
gc.set_vec_bool(10);  // a vector of 10 booleans
GC::vec_bool_type bv; // initialize an empty vector of booleans
bv.push_bach(true); bv.push_bach(false);
gc.set_vec_bool(bv);  // a vector of 2 booleans copy of bv

gc.set_vec_int(10);   // a vector of 10 integers
GC::vec_int_type iv;  // initialize an empty vector of integers
iv.push_back(1); iv.push_back(2); iv.push_back(-1);
gc.set_vec_int(iv);   // a vector of 3 integers copy of iv

gc.set_vec_real(10);  // a vector of 10 floating point numbers
GC::vec_real_type rv; // initialize an empty vector of integers
rv.push_back(1.4); rv.push_back(2.1); rv.push_back(-1);
gc.set_vec_int(rv);   // a vector of 3 floating point copy of rv
~~~~~~~~~~~~~

to a vector of strings or pointers
 
~~~~~~~~~~~~~{.cc}
gc.set_vec_string(10);  // a vector of 10 strings
GC::vec_string_type sv; // initialize an empty vector of booleans
sv.push_bach("pippo"); sv.push_bach("pluto");
gc.set_vec_string(sv);  // a vector of 2 string copy of sv
 
gc.set_vec_pointer(10);  // a vector of 10 pointers
GC::vec_pointer_type pv; // initialize an empty vector of pointers
pv.push_back(&cout); pv.push_back(&cin);
gc.set_vec_pointer(pv);  // a vector of 2 pointers copy of pv
~~~~~~~~~~~~~
 
To build complex aggregate data a generic vector and generic
map data type are available:
 
~~~~~~~~~~~~~{.cc}
gc.set_vector(10); // a generic vector of 10 elements
gc.set_map();      // an empty generic map of
~~~~~~~~~~~~~

How to access to the data stored in `GenericContainer` objects
are discussed in section \ref sec3

============================

\section sec2 Implicit type initialization

A generic container can be initialized empty or to a specific value

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc1;          // initialize empty container
GC::GenericContainer gc2(1);       // store an integer
GC::GenericContainer gc3(1.2);     // store a floating point
GC::GenericContainer gc4("pippo"); // store a string
GC::GenericContainer gc5(true);    // store a bool
GC::GenericContainer gc6(&cout);   // store a pointer
GC::GenericContainer gc7(gc6);     // store a copy of gc6, a pointer
GC::GenericContainer gc8(gc1);     // store a copy of gc1, no data
~~~~~~~~~~~~~

getting information 

~~~~~~~~~~~~~{.cc}
gc1.info(cout); // print the type stored in the `GenericContainer`
gc2.info(cout);
gc3.info(cout);
gc4.info(cout);
gc5.info(cout);
gc6.info(cout);
gc7.info(cout);
gc8.info(cout);
~~~~~~~~~~~~~

result in

~~~~~~~~~~~~~{.cc}
GenericContainer: No data stored
Integer: 1
Floating Point: 1.2
String: pippo
Boolean: true
Generic pointer: 7fff74272f48
Generic pointer: 7fff74272f48
GenericContainer: No data stored
~~~~~~~~~~~~~

Initialization with operator =
------------------------------
 
A generic container can be initialized using `operator =`

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc, gc1;
gc.info(cout);
gc = 1;       gc.info(cout);
gc = 1.2;     gc.info(cout);
gc = "pippo"; gc.info(cout);
gc = true;    gc.info(cout);
gc = &cout;   gc.info(cout);
gc = gc1;     gc.info(cout);
~~~~~~~~~~~~~

the output is:

~~~~~~~~~~~~~{.cc}
GC::GenericContainer: No data stored
Integer: 1
Floating Point: 1.2
String: pippo
Boolean: true
Generic pointer: 7fff74272f48
~~~~~~~~~~~~~

============================
 
\section sec3 Accessing data stored in vector
 
To retrieve the data stored in a `GenericContainer` you can use the 
following methods:
 
~~~~~~~~~~~~~{.cc}
bool   b = gc.get_bool();          // to access a boolean
int    i = gc.get_int();           // to access an integer
double r = gc.get_real();          // to access a floating point number
string s = gc.get_string();        // to access a string
int * p  = gc.get_pointer<int*>(); // to access a pointer as pointer to integer
~~~~~~~~~~~~~

if you request to access, for example, an integer with `gc.get_int()`
and the container store, for example, a `string` a run time error
is issued.

The access to generic vector can be done in 3 way

accessing by using references (alias)
-------------------------------------

~~~~~~~~~~~~~{.cc}
GC::vec_bool_type & bv = gc.get_bool_vec(); // make a reference of the vector of booleans
bv[0] = true; // Access the elements [read/write]
bv[1] = false;

GC::vec_int_type & iv = gc.get_int_vec(); // make a reference of the vector of integers
iv[0] = 1; // Access the elements [read/write]
iv[1] = 4;

GC::vec_real_type & rv = gc.get_real_vec(); // make a reference of the vector of floating point numbers
rv[0] = 1; // Access the elements [read/write]
rv[1] = 4.5;

GC::vec_string_type & sv = gc.get_string_vec(); // make a reference of the vector of strings
sv[0] = "pippo"; // Access the elements [read/write]
sv[1] = "pluto";

GC::vec_pointer_type & pv = gc.get_pointer_vec(); // make a reference of the vector of pointers
pv[0] = &cout; // Access the elements [read/write]
pv[1] = &cin;
~~~~~~~~~~~~~

elements can be generic vector or generic maps

~~~~~~~~~~~~~{.cc}
GC::vector_type & gv = gc.get_vector(); // make a reference of the generic vector
gv[0] = 1;   // access first element of generic vector
gv[1] = 1.3; // access second element of generic vector

GC::map_type & m = gc.get_map(); // make a reference of the generic map
m["pippo"] = 1; // access element "pippo" of the generic map
m["pluto"] = 4; // access element "pluto" of the generic map
~~~~~~~~~~~~~

accessing directly the i-th element
-------------------------------------

~~~~~~~~~~~~~{.cc}
gc.get_bool(i)           = true; // Access the i-th element of vector of booleans
gc.get_int(i)            = 123;  // Access the i-th element of vector of integers
gc.get_real(i)           = 1.23; // Access the i-th element of vector of floating point numbers
gc.get_string(i)         = "pippo"; // Access the i-th element of vector of strings
gc.get_pointer<void*>(i) = &cout;   // Access the i-th element of vector of pointers
gc.get_pointer<void*>(i) = &cout;   // Access the i-th element of vector of pointers
~~~~~~~~~~~~~

from a generic vector

~~~~~~~~~~~~~{.cc}
GC::GenericContainer & c = gc.get_gc(i); // make a reference of a `GenericContainer` at i-th position
c.get_bool() = true; // if the element is a boolean set it
c.set_bool(true);    // equivalent way
c = true;            // equivalent way
~~~~~~~~~~~~~

accessing directly the i-th element using operator [] and ()
------------------------------------------------------------

~~~~~~~~~~~~~{.cc}
gc[i] = true;    // Access the i-th element of vector of booleans
gc[i] = 123;     // Access the i-th element of vector of integers
gc[i] = 1.23;    // Access the i-th element of vector of floating point numbers
gc[i] = "pippo"; // Access the i-th element of vector of strings
gc[i] = &cout;   // Access the i-th element of vector of pointers
~~~~~~~~~~~~~

from a generic vector

~~~~~~~~~~~~~{.cc}
gc[i].get_bool() = true; // Access the i-th element and set it
gc[i].set_bool(true);    // equivalent way
gc[i] = true;            // equivalent way
~~~~~~~~~~~~~

operator () do the same work. The difference is that operator []
rewrite the object with a new type if assigned with a different type.
For example if gc store a generic vector:

~~~~~~~~~~~~~{.cc}
gc[i] = true;    // set to boolean
gc[i] = "pippo"; // change type to string
~~~~~~~~~~~~~

while operator () cannot change type of the object nor initialize it:

~~~~~~~~~~~~~{.cc}
gc[i] = true;    // set to boolean
gc(i) = "pippo"; // run time error cannot change allocation type
~~~~~~~~~~~~~

============================
 
\section sec4 Accessing data stored in map

Map are associative array indexed with strings.
To define a map you can initialize in many ways:

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc; // empty object
gc.set_map();
GC::map_type & m = gc.get_map(); // get an alias of the map data
m["pippo"] = 1; // access element "pippo" of the generic map
m["pluto"] = 4; // access element "pluto" of the generic map
// equivalent way
gc["pippo"] = 1; // access element "pippo" of the generic map
gc["pluto"] = 4; // access element "pluto" of the generic map
~~~~~~~~~~~~~

operator [] can initialize a map

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc = 1; // create `GenericContainer` which store integer 1
gc["pippo"] = 1; // gc is reallocated as a map and store 1 at index "pippo"
gc["pluto"] = 4; // access element "pluto" of the generic map
~~~~~~~~~~~~~

============================
 
\section sec5 Build complex data structures

For more complex examples and recursive data see examples files
in the distribution.

*/

//
// eof: GenericContainer_doc.hh
//
