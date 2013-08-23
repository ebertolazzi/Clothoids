GenericContainer is a C++ class with permit to store eterogeneous data:
  
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

------------------

To compile the test do:

make

To run the test

make run

To generate documentation

make doc

-----------------

Enrico Bertolazzi
Dipartimento di Ingegneria Industriale
Universita` degli Studi di Trento
email: enrico.bertolazzi@unitn.it
