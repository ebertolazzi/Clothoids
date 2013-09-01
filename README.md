GenericContainer
================

`GenericContainer` is a C++ class with permit to store eterogeneous data:
  
  - pointer
  - boolean
  - integer
  - floating point
  - string
  - vector of pointer
  - vector of boolean
  - vector of integer
  - vector of floating point
  - vector of string

in addition to this data type the following two container are added

  - vector of `GenericContainer`
  - map of `GenericContainer`

this permits to build complex recursive data.
The main usage of the class is in interchange data with
scripting language like `Ruby`, `Lua`, `MATLAB`.

Simple Usage
------------

The usage is simple, for example it
can be used as an associative array with eterogenous data

~~~~~~~~~~~~~
GenericContainer gc ;
gc["one"]  = 1       ; // store integer
gc["two"]  = true    ; // store a boolean
gc["3"]    = 1.4     ; // store floating point number
gc["four"] = "pippo" ; // store a string
gc["five"].set_vec_int(10) ; // store a vector of integer of 10 elements
~~~~~~~~~~~~~

and to retrieve elements

~~~~~~~~~~~~~
cout << gc["one"].get_int()     << '\n' ;
cout << gc["two"].get_bool()    << '\n' ;
cout << gc["3"].get_real()      << '\n' ;
cout << gc["four"].get_string() << '\n' ;
GenericContainer::vec_int_type & v = gc["five"].get_vec_int();
cout << v[1] << '\n' ;
~~~~~~~~~~~~~

For more complex emxamples and recursive data see example test files
in the distribution.


Compile and tests
-----------------

Edit makefile file to match compiler of your OS and do:

  make

To run the test

  make run

To generate documentation (using DOXYGEN: http://www.stack.nl/~dimitri/doxygen/index.html)

make doc

DOXYGEN documentation
---------------------
Available at: http://www.ing.unitn.it/~bertolaz/4-software/genericContainer/index.html

* * *

Enrico Bertolazzi<br>
Dipartimento di Ingegneria Industriale<br>
Universita` degli Studi di Trento<br>
email: enrico.bertolazzi@unitn.it
