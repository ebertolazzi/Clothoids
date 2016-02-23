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
// file: GenericContainer.hh
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
 - https://bitbucket.org/ebertolazzi/genericcontainer

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
GC::GenericContainer gc ;
gc["one"]  = 1       ; // store integer
gc["two"]  = true    ; // store a boolean
gc["3"]    = 1.4     ; // store floating point number
gc["four"] = "pippo" ; // store a string
gc["five"].set_vec_int(10) ; // store a vector of integer of 10 elements
~~~~~~~~~~~~~

and to retrieve elements

~~~~~~~~~~~~~{.cc}
cout << gc["one"].get_int()     << '\n' ;
cout << gc["two"].get_bool()    << '\n' ;
cout << gc["3"].get_real()      << '\n' ;
cout << gc["four"].get_string() << '\n' ;
GC::vec_int_type & v = gc["five"].get_vec_int();
cout << v[1] << '\n' ;
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
gc.set_bool(true) ;
gc.set_bool(false) ;
~~~~~~~~~~~~~

to an integer

~~~~~~~~~~~~~{.cc}
gc.set_int(123) ;
~~~~~~~~~~~~~

to a floating point number

~~~~~~~~~~~~~{.cc}
gc.set_real(1.23) ;
gc.set_real(3) ;
~~~~~~~~~~~~~

to a string

~~~~~~~~~~~~~{.cc}
gc.set_string("a C string") ;
string s = "a C++ sring";
gc.set_string(s) ;
~~~~~~~~~~~~~

to a pointer

~~~~~~~~~~~~~{.cc}
gc.set_pointer(&cout) ;
~~~~~~~~~~~~~

to a vector of boolean, integer or floating points

~~~~~~~~~~~~~{.cc}
gc.set_vec_bool(10)  ; // a vector of 10 booleans
GC::vec_bool_type bv ; // initialize an empty vector of booleans
bv.push_bach(true) ; bv.push_bach(false) ;
gc.set_vec_bool(bv)  ; // a vector of 2 booleans copy of bv

gc.set_vec_int(10)  ; // a vector of 10 integers
GC::vec_int_type iv ; // initialize an empty vector of integers
iv.push_back(1) ; iv.push_back(2) ; iv.push_back(-1) ;
gc.set_vec_int(iv)  ; // a vector of 3 integers copy of iv

gc.set_vec_real(10)  ; // a vector of 10 floating point numbers
GC::vec_real_type rv ; // initialize an empty vector of integers
rv.push_back(1.4) ; rv.push_back(2.1) ; rv.push_back(-1) ;
gc.set_vec_int(rv)   ; // a vector of 3 floating point copy of rv
~~~~~~~~~~~~~

to a vector of strings or pointers
 
~~~~~~~~~~~~~{.cc}
gc.set_vec_string(10)  ; // a vector of 10 strings
GC::vec_string_type sv ; // initialize an empty vector of booleans
sv.push_bach("pippo") ; sv.push_bach("pluto") ;
gc.set_vec_string(sv)  ; // a vector of 2 string copy of sv
 
gc.set_vec_pointer(10)  ; // a vector of 10 pointers
GC::vec_pointer_type pv ; // initialize an empty vector of pointers
pv.push_back(&cout) ; pv.push_back(&cin) ;
gc.set_vec_pointer(pv)  ; // a vector of 2 pointers copy of pv
~~~~~~~~~~~~~
 
To build complex aggregate data a generic vector and generic
map data type are available:
 
~~~~~~~~~~~~~{.cc}
gc.set_vector(10) ; // a generic vector of 10 elements
gc.set_map()      ; // an empty generic map of 
~~~~~~~~~~~~~

How to access to the data stored in `GenericContainer` objects
are discussed in section \ref sec3

============================

\section sec2 Implicit type initialization

A generic container can be initialized empty or to a specific value

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc1          ; // initialize empty container
GC::GenericContainer gc2(1)       ; // store an integer
GC::GenericContainer gc3(1.2)     ; // store a floating point
GC::GenericContainer gc4("pippo") ; // store a string
GC::GenericContainer gc5(true)    ; // store a bool
GC::GenericContainer gc6(&cout)   ; // store a pointer
GC::GenericContainer gc7(gc6)     ; // store a copy of gc6, a pointer
GC::GenericContainer gc8(gc1)     ; // store a copy of gc1, no data
~~~~~~~~~~~~~

getting information 

~~~~~~~~~~~~~{.cc}
gc1.info(cout) ; // print the type stored in the `GenericContainer`
gc2.info(cout) ;
gc3.info(cout) ;
gc4.info(cout) ;
gc5.info(cout) ;
gc6.info(cout) ;
gc7.info(cout) ;
gc8.info(cout) ;
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
GC::GenericContainer gc, gc1 ;
gc.info(cout) ;
gc = 1       ; gc.info(cout) ;
gc = 1.2     ; gc.info(cout) ;
gc = "pippo" ; gc.info(cout) ;
gc = true    ; gc.info(cout) ;
gc = &cout   ; gc.info(cout) ;
gc = gc1     ; gc.info(cout) ;
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
bool   b = gc.get_bool()    ; // to access a boolean
int    i = gc.get_int()     ; // to access an integer
double r = gc.get_real()    ; // to access a floating point number
string s = gc.get_string()  ; // to access a string
int * p  = gc.get_pointer<int*>() ; // to access a pointer as pointer to integer
~~~~~~~~~~~~~

if you request to access, for example, an integer with `gc.get_int()`
and the container store, for example, a `string` a run time error
is issued.

The access to generic vector can be done in 3 way

accessing by using references (alias)
-------------------------------------

~~~~~~~~~~~~~{.cc}
GC::vec_bool_type & bv = gc.get_bool_vec() ; // make a reference of the vector of booleans
bv[0] = true  ; // Access the elements [read/write]
bv[1] = false ;

GC::vec_int_type & iv = gc.get_int_vec() ; // make a reference of the vector of integers
iv[0] = 1 ; // Access the elements [read/write]
iv[1] = 4 ;

GC::vec_real_type & rv = gc.get_real_vec() ; // make a reference of the vector of floating point numbers
rv[0] = 1 ; // Access the elements [read/write]
rv[1] = 4.5 ;

GC::vec_string_type & sv = gc.get_string_vec() ; // make a reference of the vector of strings
sv[0] = "pippo" ; // Access the elements [read/write]
sv[1] = "pluto" ;

GC::vec_pointer_type & pv = gc.get_pointer_vec() ; // make a reference of the vector of pointers
pv[0] = &cout ; // Access the elements [read/write]
pv[1] = &cin ;
~~~~~~~~~~~~~

elements can be generic vector or generic maps

~~~~~~~~~~~~~{.cc}
GC::vector_type & gv = gc.get_vector() ; // make a reference of the generic vector
gv[0] = 1   ; // access first element of generic vector
gv[1] = 1.3 ; // access second element of generic vector

GC::map_type & m = gc.get_map() ; // make a reference of the generic map
m["pippo"] = 1 ; // access element "pippo" of the generic map
m["pluto"] = 4 ; // access element "pluto" of the generic map
~~~~~~~~~~~~~

accessing directly the i-th element
-------------------------------------

~~~~~~~~~~~~~{.cc}
gc.get_bool(i)           = true    ; // Access the i-th element of vector of booleans
gc.get_int(i)            = 123     ; // Access the i-th element of vector of integers
gc.get_real(i)           = 1.23    ; // Access the i-th element of vector of floating point numbers
gc.get_string(i)         = "pippo" ; // Access the i-th element of vector of strings
gc.get_pointer<void*>(i) = &cout   ; // Access the i-th element of vector of pointers
gc.get_pointer<void*>(i) = &cout   ; // Access the i-th element of vector of pointers
~~~~~~~~~~~~~

from a generic vector

~~~~~~~~~~~~~{.cc}
GC::GenericContainer & c = gc.get_gc(i) ; // make a reference of a `GenericContainer` at i-th position
c.get_bool() = true ; // if the element is a boolean set it
c.set_bool(true)    ; // equivalent way
c = true            ; // equivalent way
~~~~~~~~~~~~~

accessing directly the i-th element using operator [] and ()
------------------------------------------------------------

~~~~~~~~~~~~~{.cc}
gc[i] = true    ; // Access the i-th element of vector of booleans
gc[i] = 123     ; // Access the i-th element of vector of integers
gc[i] = 1.23    ; // Access the i-th element of vector of floating point numbers
gc[i] = "pippo" ; // Access the i-th element of vector of strings
gc[i] = &cout   ; // Access the i-th element of vector of pointers
~~~~~~~~~~~~~

from a generic vector

~~~~~~~~~~~~~{.cc}
gc[i].get_bool() = true ; // Access the i-th element and set it
gc[i].set_bool(true)    ; // equivalent way
gc[i] = true            ; // equivalent way
~~~~~~~~~~~~~

operator () do the same work. The difference is that operator []
rewrite the object with a new type if assigned with a different type.
For example if gc store a generic vector:

~~~~~~~~~~~~~{.cc}
gc[i] = true    ; // set to boolean
gc[i] = "pippo" ; // change type to string
~~~~~~~~~~~~~

while operator () cannot change type of the object nor initialize it:

~~~~~~~~~~~~~{.cc}
gc[i] = true    ; // set to boolean
gc(i) = "pippo" ; // run time error cannot change allocation type
~~~~~~~~~~~~~

============================
 
\section sec4 Accessing data stored in map

Map are associative array indexed with strings.
To define a map you can initialize in many ways:

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc ; // empty object
gc.set_map() ;
GC::map_type & m = gc.get_map() ; // get an alias of the map data
m["pippo"] = 1 ; // access element "pippo" of the generic map
m["pluto"] = 4 ; // access element "pluto" of the generic map
// equivalent way
gc["pippo"] = 1 ; // access element "pippo" of the generic map
gc["pluto"] = 4 ; // access element "pluto" of the generic map
~~~~~~~~~~~~~

operator [] can initialize a map

~~~~~~~~~~~~~{.cc}
GC::GenericContainer gc = 1 ; // create `GenericContainer` which store integer 1
gc["pippo"] = 1 ; // gc is reallocated as a map and store 1 at index "pippo"
gc["pluto"] = 4 ; // access element "pluto" of the generic map
~~~~~~~~~~~~~

============================
 
\section sec5 Build complex data structures

For more complex examples and recursive data see examples files
in the distribution.

*/

#ifndef GENERIC_CONTAINER_HH
#define GENERIC_CONTAINER_HH

#include <iostream>
#include <string>
#include <complex>
#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <stdexcept>

#if __cplusplus > 199711L
  #define GENERIC_CONTAINER_USE_CXX11
#endif

// if C++ < C++11 define nullptr
#ifndef GENERIC_CONTAINER_USE_CXX11
  #include <cstdlib>
  #ifndef nullptr
    #include <cstddef>
    #define nullptr NULL
  #endif
#endif

#ifndef GC_ASSERT
  #define GC_ASSERT(COND,MSG)                            \
    if ( !(COND) ) {                                     \
      std::ostringstream ost ;                           \
      ost << "in GenericContainer: " << MSG << '\n' ;    \
      GenericContainer::exception( ost.str().c_str() ) ; \
    }
#endif

#ifndef GC_WARNING
  #define GC_WARNING(COND,MSG)                                 \
    if ( !(COND) ) {                                           \
      std::cout << "On line: " << __LINE__                     \
                << " file: " << __FILE__                       \
                << " in GenericContainer\nWARNING: " << MSG << '\n' ;  \
    }
#endif

#if defined(_WIN32) || defined(_WIN64)
  #define GENERIC_CONTAINER_ON_WINDOWS
#endif

#ifndef GENERIC_CONTAINER_API_DLL
  #ifdef GENERIC_CONTAINER_ON_WINDOWS
    #ifdef GENERIC_CONTAINER_EXPORT
      #define GENERIC_CONTAINER_API_DLL __declspec(dllexport)
    #elif defined(GENERIC_CONTAINER_IMPORT)
     #define GENERIC_CONTAINER_API_DLL __declspec(dllimport)
    #else
      #define GENERIC_CONTAINER_API_DLL
    #endif
  #else
    #define GENERIC_CONTAINER_API_DLL
  #endif
#endif

namespace GenericContainerNamepace {

  class GenericContainer ;

  typedef void*                   pointer_type ; //!< generic pointer type
  typedef bool                    bool_type    ; //!< boolean type data
  typedef int                     int_type     ; //!< integer type data
  typedef long                    long_type    ; //!< long integer type data
  typedef double                  real_type    ; //!< floating point type data
  typedef std::complex<real_type> complex_type ; //!< complex floating point type data
  typedef std::string             string_type  ; //!< string type data

  typedef std::vector<pointer_type> vec_pointer_type ; //!< vector of generic pointer
  typedef std::vector<bool_type>    vec_bool_type    ; //!< vector of boolean
  typedef std::vector<int_type>     vec_int_type     ; //!< vector of integer
  typedef std::vector<real_type>    vec_real_type    ; //!< vector of floating point
  typedef std::vector<complex_type> vec_complex_type ; //!< vector of complex floating point
  typedef std::vector<string_type>  vec_string_type  ; //!< vector of strings

  typedef std::vector<GenericContainer>          vector_type ; //!< vector of `GenericContainer`
  typedef std::map<string_type,GenericContainer> map_type    ; //!< associative array of `GenericContainer`

  // ---------------------------------------------------------------------------

  class mat_real_type : public vec_real_type {
    unsigned _numRows ;
    unsigned _numCols ;
  public:

    mat_real_type()
    : _numRows(0)
    , _numCols(0)
    {}

    mat_real_type( unsigned nr, unsigned nc )
    : _numRows(nr)
    , _numCols(nc)
    { this->resize(size_type(nr*nc)) ; }

    unsigned numRows() const { return _numRows ; }
    unsigned numCols() const { return _numCols ; }

    real_type const & operator () ( unsigned i, unsigned j ) const ;
    real_type       & operator () ( unsigned i, unsigned j ) ;

  } ;

  // ---------------------------------------------------------------------------

  class mat_complex_type : public vec_complex_type {
    unsigned _numRows ;
    unsigned _numCols ;
  public:

    mat_complex_type()
    : _numRows(0)
    , _numCols(0)
    {}

    mat_complex_type( unsigned nr, unsigned nc )
    : _numRows(nr)
    , _numCols(nc)
    { this->resize(size_type(nr*nc)) ; }
    
    unsigned numRows() const { return _numRows ; }
    unsigned numCols() const { return _numCols ; }

    complex_type const & operator () ( unsigned i, unsigned j ) const ;
    complex_type       & operator () ( unsigned i, unsigned j ) ;

  } ;

  // ---------------------------------------------------------------------------

  std::ostream & operator << ( std::ostream & s, vec_pointer_type const & v ) ;
  std::ostream & operator << ( std::ostream & s, vec_bool_type const & v ) ;
  std::ostream & operator << ( std::ostream & s, vec_int_type const & v ) ;
  std::ostream & operator << ( std::ostream & s, vec_real_type const & ) ;
  std::ostream & operator << ( std::ostream & s, vec_complex_type const & ) ;
  std::ostream & operator << ( std::ostream & s, mat_real_type const & ) ;
  std::ostream & operator << ( std::ostream & s, mat_complex_type const & ) ;

  // ---------------------------------------------------------------------------

  //! Type allowed for the `GenericContainer`
  enum TypeAllowed {
    // simple type
    GC_NOTYPE=0,
    GC_POINTER,
    GC_BOOL,
    GC_INTEGER,
    GC_LONG,
    GC_REAL,
    GC_COMPLEX,
    GC_STRING,

    // vector type
    GC_VEC_POINTER,
    GC_VEC_BOOL,
    GC_VEC_INTEGER,
    GC_VEC_REAL,
    GC_VEC_COMPLEX,
    GC_VEC_STRING,

    // matrix type
    GC_MAT_REAL,
    GC_MAT_COMPLEX,

    // complex type
    GC_VECTOR,
    GC_MAP
  } ;

  /*!
    \brief `GenericContainer` is a class which permit to store eterogeneous data:

    - pointer
    - boolean
    - integer
    - long integer
    - floating point
    - complex floating point
    - string
    - vector of pointer
    - vector of boolean
    - vector of integer
    - vector of floating point
    - vector of complex floating point
    - matrix of floating point
    - matrix of complex floating point
    - vector of string

    in addition to this data type the following two container are added

    - vector of `GenericContainer`
    - map of `GenericContainer`

   */
  class GENERIC_CONTAINER_API_DLL GenericContainer {
  private:

    //! Data is stored in a union
    typedef union {
      pointer_type     p ;
      bool_type        b ;
      int_type         i ;
      long_type        l ;
      real_type        r ;
      complex_type     * c ;
      string_type      * s ;

      vec_pointer_type * v_p ;
      vec_bool_type    * v_b ;
      vec_int_type     * v_i ;
      vec_real_type    * v_r ;
      vec_complex_type * v_c ;
      mat_real_type    * m_r ;
      mat_complex_type * m_c ;
      vec_string_type  * v_s ;

      vector_type      * v ;
      map_type         * m ;

    } DataStorage ;

    DataStorage _data      ; //!< The data stored in the class instance
    TypeAllowed _data_type ; //!< The kind of data stored

    void allocate_string() ;
    void allocate_complex() ;

    void allocate_vec_pointer( unsigned sz ) ;
    void allocate_vec_bool( unsigned sz ) ;
    void allocate_vec_int( unsigned sz ) ;
    void allocate_vec_real( unsigned sz ) ;
    void allocate_vec_complex( unsigned sz ) ;
    void allocate_mat_real( unsigned nr, unsigned nc ) ;
    void allocate_mat_complex( unsigned nr, unsigned nc ) ;
    void allocate_vec_string( unsigned sz ) ;

    void allocate_vector( unsigned sz ) ;
    void allocate_map() ;

    void ck(char const [],TypeAllowed) const ;
    int  ck(TypeAllowed) const ;
    void ck_or_set(char const [], TypeAllowed) ;

    void initialize() { _data_type = GC_NOTYPE ; }

    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    bool simple_data()     const ;
    bool simple_vec_data() const ;
    #else
    bool simple_data()     const { return _data_type <= GC_STRING     ; }
    bool simple_vec_data() const { return _data_type <= GC_VEC_STRING ; }
    #endif

  public:

    //! build an instance of `GenericContainer` with empty data
    GenericContainer() ;

    //! destroy the instance of `GenericContainer`
    ~GenericContainer() { clear() ; }

    //! free memory of the data stored in `GenericContainer`, data type become `NOTYPE`
    void clear() ;
    
    //! \name Initialize simple data
    //@{
    //! Set data to `pointer_type` initialize and return a reference to the data
    pointer_type & set_pointer( pointer_type value );

    //! Free pointer to memory pointed, set data to `NO TYPE` initialize and return a reference to the data
    GenericContainer & free_pointer() ;

    //! Set data to `bool_type` initialize and return a reference to the data
    bool_type & set_bool( bool_type value ) ;

    //! Set data to `int_type` initialize and return a reference to the data
    int_type & set_int( int_type value ) ;

    //! Set data to `int_type` initialize and return a reference to the data
    long_type & set_long( long_type value ) ;

    //! Set data to `real_type` initialize and return a reference to the data
    real_type & set_real( real_type value ) ;

    //! Set data to `complex_type` initialize and return a reference to the data
    complex_type & set_complex( complex_type & value ) ;

    //! Set data to `complex_type` initialize and return a reference to the data
    complex_type & set_complex( real_type r, real_type i ) ;

    //! Set data to `string_type`, allocate and initialize. Return a reference to the data
    string_type & set_string( string_type const & value ) ;
    //@}

    //! \name Initialize vector data
    //@{
    /*! \brief 
        Set data to `vec_pointer_type`, allocate and initialize.
        Return a reference to vector of pointer.
        If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_pointer_type & set_vec_pointer( unsigned sz = 0 ) ;

    /*! \brief
     Set data to `vec_pointer_type`, allocate and initialize.
     Return a reference to vector of pointer.
     Copy the data from vector `v`.
     */
    vec_pointer_type & set_vec_pointer( vec_pointer_type const & v ) ;

    /*! \brief
        Set data to `vec_bool_type`, allocate and initialize.
        Return a reference to vector of booleans.
        If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_bool_type & set_vec_bool( unsigned sz = 0 ) ;

    /*! \brief
     Set data to `vec_bool_type`, allocate and initialize.
     Return a reference to vector of `bool`.
     Copy the data from vector `v`.
     */
    vec_bool_type & set_vec_bool( vec_bool_type const & v ) ;

    /*! \brief
        Set data to `vec_int_type`, allocate and initialize.
        Return a reference to vector of integers.
        If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_int_type & set_vec_int( unsigned sz = 0 ) ;

    /*! \brief
     Set data to `vec_int_type`, allocate and initialize.
     Return a reference to vector of integer.
     Copy the data from vector `v`.
     */
    vec_int_type & set_vec_int( vec_int_type const & v ) ;

    /*! \brief
        Set data to `vec_real_type`, allocate and initialize.
        Return a reference to vector of floating point numbers.
     If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_real_type & set_vec_real( unsigned sz = 0 ) ;

    /*! \brief
     Set data to `vec_real_type`, allocate and initialize.
     Return a reference to vector of floating point number.
     Copy the data from vector `v`.
     */
    vec_real_type & set_vec_real( vec_real_type const & v ) ;


    /*! \brief
        Set data to `vec_complex_type`, allocate and initialize.
        Return a reference to vector of complex floating point numbers.
     If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_complex_type & set_vec_complex( unsigned sz = 0 ) ;

    /*! \brief
     Set data to `vec_complex_type`, allocate and initialize.
     Return a reference to vector of complex floating point number.
     Copy the data from vector `v`.
     */
    vec_complex_type & set_vec_complex( vec_complex_type const & v ) ;

    /*! \brief
        Set data to `mat_real_type`, allocate and initialize.
        Return a reference to a matrix of floating point numbers.
     If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
     */
    mat_real_type & set_mat_real( unsigned nr = 0, unsigned nc = 0 ) ;

    /*! \brief
     Set data to `mat_real_type`, allocate and initialize.
     Return a reference to a matrix of floating point number.
     Copy the data from matrix `m`.
     */
    mat_real_type & set_mat_real( mat_real_type const & m ) ;

    /*! \brief
        Set data to `mat_complex_type`, allocate and initialize.
        Return a reference to a matrix of complex floating point numbers.
     If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
     */
    mat_complex_type & set_mat_complex( unsigned nr = 0, unsigned nc = 0 ) ;

    /*! \brief
     Set data to `mat_complex_type`, allocate and initialize.
     Return a reference to a matrix of floating point number.
     Copy the data from matrix `m`.
     */
    mat_complex_type & set_mat_complex( mat_complex_type const & m ) ;

    /*! \brief
        Set data to `vec_string_type`, allocate and initialize.
        Return a reference to vector of strings.
        If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vec_string_type & set_vec_string( unsigned sz = 0 ) ;
    
    /*! \brief
     Set data to `vec_string_type`, allocate and initialize.
     Return a reference to vector of strings.
     Copy the data from vector `v`.
     */
    vec_string_type & set_vec_string( vec_string_type const & v ) ;

    //! push boolean if data is `vec_bool_type` or `vector_type'
    void push_bool( bool ) ;

    //! push integer if data is `vec_int_type` or `vector_type'
    void push_int( int_type ) ;

    //! push integer if data is `vec_int_type` or `vector_type'
    void push_long( long_type ) ;

    //! push real if data is `vec_real_type` or `vector_type'
    void push_real( real_type ) ;

    //! push complex if data is `vec_complex_type` or `vector_type'
    void push_complex( complex_type & ) ;

    //! push complex if data is `vec_complex_type` or `vector_type'
    void push_complex( real_type re, real_type im ) ;

    //! push complex if data is `vec_string_type` or `vector_type'
    void push_string( string_type const & ) ;

    //@}

    //! \name Initialize generic data
    //@{
    /*! \brief 
        Set data to `vector_type`, allocate an empty generic vector and return a reference to it.
        If `sz` > 0 then the vector is allocated to size `sz`.
     */
    vector_type & set_vector( unsigned sz = 0 ) ;
    
    //! Set data to `map_type`, allocate an empty generic map and return a reference to it.
    map_type & set_map() ;
    //@}
    
    //! \name Access to a single element
    //@{
    
    //! Return an integer representing the type of data stored
    /*!
       Integer to data type map
       ------------------------
   
         -   No data stored (return 0)
         1.  `pointer_type`
         2.  `bool_type`
         3.  `int_type`
         4.  `long_type`
         5.  `real_type`
         6.  `complex_type`
         7.  `string_data`
         8.  `vec_pointer_type`
         9.  `vec_bool_type`
         10. `vec_int_type`
         11. `vec_real_type`
         12. `vec_complex_type`
         13. `vec_string_type`
         14. `mat_real_type`
         15. `mat_complex_type`
         16. `vector_type`
         17. `map_type`

    */
    TypeAllowed get_type() const { return _data_type ; }
    
    //! Return a string pointer representing the type of data stored
    char const * get_type_name() const ;

    //! Print to stream the kind of data stored
    GenericContainer const & info( std::basic_ostream<char> & stream ) const ;

    /*! Return the number of the elements of the first level of the generic container
    //  1 for single element, the size of vector or map, or 0
    */
    unsigned get_num_elements() const ;

    unsigned get_numRows() const ;
    unsigned get_numCols() const ;

    //! If data is boolean, integer or floating point return number, otherwise return `0`.
    real_type get_number() const ;
    complex_type get_complex_number() const ;
    void get_complex_number( real_type & re, real_type & im ) const ;

    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    void * get_pvoid() const ;
    void ** get_ppvoid() const ;

    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_POINTER) ; return *(T*)(get_ppvoid()) ; }

    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_POINTER) ; return static_cast<T>(get_pvoid()) ; }

    #else
    void * get_pvoid() const
    { ck("get_pointer",GC_POINTER) ; return _data.p ; }

    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_POINTER) ; return *(T*)(&(_data.p)) ; }

    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_POINTER) ; return static_cast<T>(_data.p) ; }
    #endif

    bool_type       & get_bool( char const msg[] = "get_bool()" ) ;
    bool_type const & get_bool( char const msg[] = "get_bool()" ) const ;
    //!< Return the stored boolean (if fails issue an error).

    int_type       & get_int( char const msg[] = "get_int()" ) ;
    int_type const & get_int( char const msg[] = "get_int()" ) const ;
    //!< Return the stored integer (if fails issue an error).

    long_type       & get_long( char const msg[] = "get_long()") ;
    long_type const & get_long( char const msg[] = "get_long()") const ;
    //!< Return the stored long integer (if fails issue an error).

    real_type       & get_real( char const msg[] = "get_real()" ) ;
    real_type const & get_real( char const msg[] = "get_real()" ) const ;
    //!< Return the stored floating point (if fails issue an error).

    complex_type       & get_complex( char const msg[] = "get_complex()" ) ;
    complex_type const & get_complex( char const msg[] = "get_complex()" ) const ;
    //!< Return the stored complex floating point (if fails issue an error).

    string_type       & get_string( char const msg[] = "get_string()" ) ;
    string_type const & get_string( char const msg[] = "get_string()" ) const ;
    //!< Return the stored string (if fails issue an error).
    //@}

    //! \name Access to vector type data
    //@{
    vector_type       & get_vector( char const msg[] = "get_vector()" ) ;
    vector_type const & get_vector( char const msg[] = "get_vector()" ) const ;
    //!< Return reference to a generic vector (if fails issue an error).

    vec_pointer_type       & get_vec_pointer( char const msg[] = "get_vec_pointer()" ) ;
    vec_pointer_type const & get_vec_pointer( char const msg[] = "get_vec_pointer()" ) const ;
    //!< Return reference to a vector of pointer (if fails issue an error).
    
    vec_bool_type       & get_vec_bool( char const msg[] = "get_vec_bool()" ) ;
    vec_bool_type const & get_vec_bool( char const msg[] = "get_vec_bool()" ) const ;
    //!< Return reference to a vector of booleans (if fails issue an error).

    vec_int_type       & get_vec_int( char const msg[] = "get_vec_int()" ) ;
    vec_int_type const & get_vec_int( char const msg[] = "get_vec_int()" ) const ;
    //!< Return reference to a vector of integers (if fails issue an error).

    vec_real_type       & get_vec_real( char const msg[] = "get_vec_real()" ) ;
    vec_real_type const & get_vec_real( char const msg[] = "get_vec_real()" ) const ;
    //!< Return reference to a vector of floating point number (if fails issue an error).

    vec_complex_type       & get_vec_complex( char const msg[] = "get_vec_complex()" ) ;
    vec_complex_type const & get_vec_complex( char const msg[] = "get_vec_complex()" ) const ;
    //!< Return reference to a vector of complex floating point number (if fails issue an error).

    mat_real_type       & get_mat_real( char const msg[] = "get_mat_real()" ) ;
    mat_real_type const & get_mat_real( char const msg[] = "get_mat_real()" ) const ;
    //!< Return reference to a matrix of floating point number (if fails issue an error).

    mat_complex_type       & get_mat_complex( char const msg[] = "get_mat_complex()" ) ;
    mat_complex_type const & get_mat_complex( char const msg[] = "get_mat_complex()" ) const ;
    //!< Return reference to a matrix of complex floating point number (if fails issue an error).

    vec_string_type       & get_vec_string( char const msg[] = "get_vec_string()" ) ;
    vec_string_type const & get_vec_string( char const msg[] = "get_vec_string()" ) const ;
    //!< Return reference to a vector of strings (if fails issue an error).
    //@}

    //! \name Access to element of vector type data
    //@{
    //! If `i`-th element of the vector is boolean, integer or floating point then return number, otherwise return `0`.
    real_type get_number_at( unsigned i ) const ;
    complex_type get_complex_number_at( unsigned i ) const ;
    void get_complex_number_at( unsigned i, real_type & re, real_type & im ) const ;

    template <typename T>
    T& get_pointer_at( unsigned i )
    { return (*this)[i].get_pointer<T>() ; }

    template <typename T>
    T get_pointer_at( unsigned i ) const
    { return (*this)[i].get_pointer<T>() ; }
    //!< Return `i`-th generic pointer (if fails issue an error).

    bool_type get_bool_at( unsigned i ) ;
    bool_type get_bool_at( unsigned i ) const ;
    //!< Return `i`-th boolean (if fails issue an error).

    int_type       & get_int_at( unsigned i ) ;
    int_type const & get_int_at( unsigned i, char const msg[] = "get_int_at(...)" ) const ;
    //!< Return `i`-th integer (if fails issue an error).

    real_type       & get_real_at( unsigned i ) ;
    real_type const & get_real_at( unsigned i ) const ;
    //!< Return `i`-th floating point number (if fails issue an error).

    complex_type       & get_complex_at( unsigned i ) ;
    complex_type const & get_complex_at( unsigned i ) const ;
    //!< Return `i`-th complex floating point number (if fails issue an error).

    real_type       & get_real_at( unsigned i, unsigned j ) ;
    real_type const & get_real_at( unsigned i, unsigned j ) const ;
    //!< Return `i`-th floating point number (if fails issue an error).

    complex_type       & get_complex_at( unsigned i, unsigned j ) ;
    complex_type const & get_complex_at( unsigned i, unsigned j ) const ;
    //!< Return `i`-th complex floating point number (if fails issue an error).

    string_type       & get_string_at( unsigned i ) ;
    string_type const & get_string_at( unsigned i ) const ;
    //!< Return `i`-th string (if fails issue an error).

    GenericContainer       & get_gc_at( unsigned i )       { return (*this)[i] ; }
    GenericContainer const & get_gc_at( unsigned i ) const { return (*this)[i] ; }
    //!< Return `i`-th `GenericContainer` of a generic vector (if fails issue an error).
    //@}

    //! \name Access to map type element
    //@{
    map_type       & get_map( char const msg[] = "get_map()" ) ;
    map_type const & get_map( char const msg[] = "get_map()" ) const ;
    //!< Return the reference of the stored map or issue an error.

    //! Check if string `s` is a key of the stored map (if fails issue an error).
    bool exists( std::string const & s ) const ;
    //@}

    //! \name Access using operators
    //@{
    
    GenericContainer & operator [] ( unsigned i ) ;
    
    /*! \brief
        Overload of the `()` operator to access the `i`-th element of a stored generic vector.
        Run time bound check.
     */
    GenericContainer const & operator [] ( unsigned i ) const ;
    
    GenericContainer & operator () ( unsigned i ) ;
    
    /*! \brief
     Overload of the `()` operator to access the `i`-th element of a stored generic vector.
     NO run time bound check.
     */
    GenericContainer const & operator () ( unsigned i ) const ;

    GenericContainer & operator [] ( std::string const & s ) ;

    /*! \brief
       Overload of the `[]` operator to access the `s`-th element of a stored generic map. 
       If the element do not exists it is created.
    */
    GenericContainer const & operator [] ( std::string const & s ) const ;

    GenericContainer & operator () ( std::string const & s ) ;

    /*! \brief
       Overload of the `()` operator to access the `i`-th element of a stored generic map. 
       If the element do not exists an error is issued.
     */
    GenericContainer const & operator () ( std::string const & s ) const ;

    //@}
    
    //! \name Initialize data using operators
    //! The `=` operator is overloaded to initialize the `GenericContainer` on its left side.
    //@{

    //! Assign a boolean to the generic container.
    GenericContainer & operator = ( bool a )
    { this -> set_bool(a) ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( unsigned a )
    { this -> set_int(int_type(a)) ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( int_type a )
    { this -> set_int(a) ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( unsigned long a )
    { this -> set_long(long_type(a)) ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( long_type a )
    { this -> set_long(a) ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( float a )
    { this -> set_real(a) ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( double a )
    { this -> set_real(a) ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( std::complex<float> & a )
    { this -> set_complex(a.real(),a.imag()) ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( std::complex<double> & a )
    { this -> set_complex(a.real(),a.imag()) ; return * this ; }

    //! Assign a string to the generic container.
    GenericContainer & operator = ( char const a[] )
    { this -> set_string(a) ; return * this ; }

    //! Assign a string to the generic container.
    GenericContainer & operator = ( std::string const & a )
    { this -> set_string(a) ; return * this ; }

    //! Assign a pointer to the generic container.
    GenericContainer & operator = ( void * a )
    { this -> set_pointer(a) ; return * this ; }

    //! Assign a generic container `a` to the generic container.
    GenericContainer const & operator = ( GenericContainer const & a )
    { this -> load( a ) ; return * this ; }

    //! Copy a generic container `a` to the generic container.
    void load( GenericContainer const & a ) ;
    //@}

    //! \name Promotion to a ``bigger'' data
    //@{
    //! If data contains a boolean it is promoted to an integer.
    GenericContainer const & promote_to_int() ;

    //! If data contains a boolean it is promoted to an integer.
    GenericContainer const & promote_to_long() ;

    //! If data contains a boolean or an integer it is promoted to a real.
    GenericContainer const & promote_to_real() ;

    //! If data contains a boolean or an integer or floating point it is promoted to a complex floating point.
    GenericContainer const & promote_to_complex() ;

    //! If data contains vector of booleans it is promoted to a vector of integer.
    GenericContainer const & promote_to_vec_int() ;

    //! If data contains vector of booleans or integer it is promoted to a vector of real.
    GenericContainer const & promote_to_vec_real() ;

    //! If data contains vector of booleans or integer or real it is promoted to a vector of complex.
    GenericContainer const & promote_to_vec_complex() ;

    //! If data contains vector of booleans, integer or real it is promoted to a matrix of real.
    GenericContainer const & promote_to_mat_real() ;

    //! If data contains vector of booleans, integer or real or complex or matrix of real it is promoted to a matrix of complex.
    GenericContainer const & promote_to_mat_complex() ;

    //! If data contains vector of someting it is promoted to a vector of `GenericContainer`.
    GenericContainer const & promote_to_vector() ;
    //@}

    //! \name Initialize data by overloading constructor
    //@{
    
    //! Construct a generic container storing a boolean
    GenericContainer( bool a ) { initialize() ; this->operator=(a) ; }
    
    //! Construct a generic container storing an integer
    GenericContainer( unsigned a ) { initialize() ; *this = a ; }
    
    //! Construct a generic container storing an integer
    GenericContainer( int a ) { initialize() ; this->operator=(a) ; }
    
    //! Construct a generic container storing a floating point number
    GenericContainer( float a ) { initialize() ; this->operator=(a) ; }

    //! Construct a generic container storing a floating point number
    GenericContainer( double a ) { initialize() ; this->operator=(a) ; }

    //! Construct a generic container storing a complex floating point number
    GenericContainer( std::complex<float> & a ) { initialize() ; this->operator=(a) ; }

    //! Construct a generic container storing a complex floating point number
    GenericContainer( std::complex<double> & a ) { initialize() ; this->operator=(a) ; }
    
    //! Construct a generic container storing a string
    GenericContainer( char const a[] ) { initialize() ; this->operator=(a) ; }

    //! Construct a generic container storing a string
    GenericContainer( std::string const & a ) { initialize() ; this->operator=(a) ; }
    
    //! Construct a generic container storing a pointer
    GenericContainer( void * a ) { initialize() ; this->operator=(a) ; }

    //! Construct a generic container copying container `gc`
    GenericContainer( GenericContainer const & gc ) { initialize() ; this->operator=(gc) ; }
    //@}

    //! \name I/O for `GenericContainer` objects
    //@{

    //! print the contents of the object in a human readable way
    void print( std::basic_ostream<char> &,
                std::string const & prefix = "",
                std::string const & indent = "    " ) const ;

    //! print the contents of the object in yaml syntax
    void to_yaml( std::basic_ostream<char> &, std::string const & prefix = "" ) const ;

    /*! 
      \brief write `GenericContainer` as regular formatted data
       
      Write the contents of the `GenericContainer` object to stream.
     
      `GenericContainer` must be a map which contains the fields:
       
      - "headers" this element must be a `vec_string_type` which contains
                  the strings of the headers of the columns of the data

      - "data"    this element must be a `vector_type` which contais the
                  vectors which are the columns of the data to be saved.
                  Each column can be of type
                  
                  1. `vec_bool_type`
                  2. `vec_int_type`
                  3. `vec_real_type`
     
                  all the vector must have the same size.

       \param stream     stream to write the output
       \param delimiter  desired delimiter (optional). Default is tab.
     */

    GenericContainer const &
    writeFormattedData( std::basic_ostream<char> & stream, char const delimiter = '\t' ) const ;

    /*! 
      \brief read regular formatted data from `stream` to `GenericContainer`.
     
      After successful read `GenericContainer` will be a map which contains the fields:
       
      - "headers"  a `vec_string_type` which contains
                   the strings of the headers of the columns of the data

      - "data"     a `vector_type` which contais the vectors which are the
                   columns of the data readed of type `vec_real_type`.

       \param stream       stream to write the output
       \param commentChars lines beginnig with one of this chars are treated as comments. 
                           Default are `#` and `%`
       \param delimiters   caracters used as delimiter for headers
     */
    GenericContainer &
    readFormattedData( std::istream & stream,
                       char const commentChars[] = "#%",
                       char const delimiters[] = " \t" ) ;

    static void exception( char const msg[] ) ;

  } ;

  // -------------------------------------------------------
  // support functions
  void
  writeTable( vec_string_type const    & headers,
              vector_type     const    & data,
              std::basic_ostream<char> & stream,
              char const delimiter = '\t' ) ;

  void
  writeTableFormatted( vec_string_type const    & headers,
                       vector_type     const    & data,
                       std::basic_ostream<char> & stream ) ;

  // -------------------------------------------------------
  class GENERIC_CONTAINER_API_DLL GenericContainerExplorer {
  private:
    enum {
      GENERIC_CONTAINER_OK        = 0,
      GENERIC_CONTAINER_BAD_TYPE  = 1,
      GENERIC_CONTAINER_NO_DATA   = 2,
      GENERIC_CONTAINER_NOT_EMPTY = 3,
      GENERIC_CONTAINER_BAD_HEAD  = 4
    } ;

    GenericContainer              data ;
    std::deque<GenericContainer*> head ;

    map_type         * ptr_map ;
    map_type::iterator map_iterator ;

  public:
  
    GenericContainerExplorer() { head.push_back(&data) ; }
    ~GenericContainerExplorer() {}

    GenericContainer *
    top() {
      GC_ASSERT( head.size() > 0, "GenericContainerExplorer::top() empty stack!" ) ;
      GC_ASSERT( head.back() != nullptr, "GenericContainerExplorer::top() bad top pointer!" ) ;
      return head.back() ;
    }

    GenericContainer const *
    top() const {
      GC_ASSERT( head.size() > 0, "GenericContainerExplorer::top() empty stack!" ) ;
      GC_ASSERT( head.back() != nullptr, "GenericContainerExplorer::top() bad top pointer!" ) ;
      return head.back() ;
    }

    void       * mem_ptr()       { return &data ; }
    void const * mem_ptr() const { return &data ; }

    int
    check( int data_type ) const {
      if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD ;
      if ( data_type == head.back() -> get_type() ) {
        if ( head.back() == nullptr ) return GENERIC_CONTAINER_NO_DATA ;
        return GENERIC_CONTAINER_OK ;
      } else {
        return GENERIC_CONTAINER_BAD_TYPE ;
      }
    }

    int
    check_no_data( int data_type ) const {
      if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD ;
      if ( GC_NOTYPE == head.back() -> get_type() ||
           data_type == head.back() -> get_type() ) return GENERIC_CONTAINER_OK ;
      return GENERIC_CONTAINER_NOT_EMPTY ;
    }

    int
    pop() {
      if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA ;
      head.pop_back() ;
      return GENERIC_CONTAINER_OK ;
    }

    int
    push( GenericContainer * gc ) {
      head.push_back( gc ) ;
      return GENERIC_CONTAINER_OK ;
    }

    int
    push_vector_position( unsigned pos ) {
      int ok = check( GC_VECTOR ) ;
      if ( ok == GENERIC_CONTAINER_OK  ) {
        GenericContainer::GenericContainer * gc = &(*head.back())[pos] ;
        head.push_back( gc ) ;
      }
      return ok ;
    }

    int
    push_map_position( char const pos[] ) {
      int ok = check( GC_MAP ) ;
      if ( ok == GENERIC_CONTAINER_OK  ) {
        GenericContainer::GenericContainer * gc = &(*head.back())[pos] ;
        head.push_back( gc ) ;
      }
      return ok ;
    }

    int
    init_map_key() {
      int ok = check( GC_MAP ) ;
      if ( ok == GENERIC_CONTAINER_OK ) {
        ptr_map = &head.back()->get_map() ;
        map_iterator = ptr_map->begin() ;
      }
      return ok ;
    }

    char const *
    next_map_key() {
      if ( map_iterator != ptr_map->end() )
        return map_iterator++ -> first.c_str() ;
      else
        return nullptr ;
    }
  
    int
    reset() {
      if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA ;
      while ( head.size() > 1 ) head.pop_back() ;
      return GENERIC_CONTAINER_OK ;
    }

  } ;
}

// do not define alias GC if use X11
#ifndef XlibSpecificationRelease
namespace GC = GenericContainerNamepace ;
#endif

#endif

//
// eof: GenericContainer.hh
//
