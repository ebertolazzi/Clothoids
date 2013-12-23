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
#include <map>
#include <vector>

// if C++ < C++11 define nullptr
#if __cplusplus <= 199711L
  #include <cstdlib>
  #ifndef nullptr
    #include <cstddef>
    #define nullptr NULL
  #endif
#endif

#ifndef GC_ASSERT
  #include <sstream>
  #include <stdexcept>
  #define GC_ASSERT(COND,MSG)                            \
    if ( !(COND) ) {                                     \
      std::ostringstream ost ;                           \
      ost << "in GenericContainer: " << MSG << '\n' ;    \
      GenericContainer::exception( ost.str().c_str() ) ; \
    }
#endif

#ifndef GC_WARNING
  #include <sstream>
  #include <stdexcept>
  #define GC_WARNING(COND,MSG)                                 \
    if ( !(COND) ) {                                           \
      std::cout << "On line: " << __LINE__                     \
                << " file: " << __FILE__                       \
                << " in GenericContainer\nWARNING: " << MSG << '\n' ;  \
    }
#endif

namespace GC {

  class GenericContainer ;

  typedef void*       pointer_type ; //!< generic pointer type
  typedef bool        bool_type    ; //!< boolean type data
  typedef int         int_type     ; //!< integer type data
  typedef double      real_type    ; //!< floating point type data
  typedef std::string string_type  ; //!< string type data

  typedef std::vector<pointer_type> vec_pointer_type ; //!< vector of generic pointer
  typedef std::vector<bool_type>    vec_bool_type    ; //!< vector of boolean
  typedef std::vector<int_type>     vec_int_type     ; //!< vector of integer
  typedef std::vector<real_type>    vec_real_type    ; //!< vector of floating point
  typedef std::vector<string_type>  vec_string_type  ; //!< vector of strings

  typedef std::vector<GenericContainer>          vector_type ; //!< vector of `GenericContainer`
  typedef std::map<string_type,GenericContainer> map_type    ; //!< associative array of `GenericContainer`

  //! Type allowed for the `GenericContainer`
  enum TypeAllowed {
    GC_NOTYPE=0,
    GC_POINTER,
    GC_BOOL,
    GC_INT,
    GC_REAL,
    GC_STRING,
    GC_VEC_POINTER,
    GC_VEC_BOOL,
    GC_VEC_INT,
    GC_VEC_REAL,
    GC_VEC_STRING,
    GC_VECTOR,
    GC_MAP
  } ;

  /*!
    \brief `GenericContainer` is a class which permit to store eterogeneous data:

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

   */
  class GenericContainer {

  public:

    #ifndef GENERIC_CONTAINER_NO_PCRE
    static void * _reCompiled      ; // pcre *
    static void * _pcreExtra       ; // pcre_extra *
    char const  * pcreErrorStr     ; // const char *
    int           pcreErrorOffset  ;
    #endif

  private:

    //! Data is stored in a union
    typedef union {
      pointer_type     p ;
      bool_type        b ;
      int_type         i ;
      real_type        r ;
      string_type      * s ;

      vec_pointer_type * v_p ;
      vec_bool_type    * v_b ;
      vec_int_type     * v_i ;
      vec_real_type    * v_r ;
      vec_string_type  * v_s ;

      vector_type      * v ;
      map_type         * m ;

    } DataStorage ;

    DataStorage _data      ; //!< The data stored in the class instance
    TypeAllowed _data_type ; //!< The kind of data stored

    void allocate_string() ;

    void allocate_vec_pointer( unsigned sz ) ;
    void allocate_vec_bool( unsigned sz ) ;
    void allocate_vec_int( unsigned sz ) ;
    void allocate_vec_real( unsigned sz ) ;
    void allocate_vec_string( unsigned sz ) ;

    void allocate_vector( unsigned sz ) ;
    void allocate_map() ;

    void ck(char const [],TypeAllowed) const ;
    int  ck(TypeAllowed) const ;
    void ck_or_set(char const [], TypeAllowed) ;

    void initialize() { _data_type = GC_NOTYPE ; }

    bool simple_data() const {
      return _data_type != GC_VECTOR &&
             _data_type != GC_MAP    &&
             _data_type != GC_VEC_STRING ;
    }

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
    
    //! Set data to `real_type` initialize and return a reference to the data
    real_type & set_real( real_type value ) ;

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
         4.  `real_type`
         5.  `string_data`
         6.  `vec_pointer_type`
         7.  `vec_bool_type`
         8.  `vec_int_type`
         9.  `vec_real_type`
         10. `vec_string_type`
         11. `vector_type`
         12. `map_type`

    */
    TypeAllowed get_type() const { return _data_type ; }
    
    //! Return a string pointer representing the type of data stored
    char const * get_type_name() const ;

    //! Print to stream the kind of data stored
    GenericContainer const & info( std::ostream & stream ) const ;

    //! If data is boolean, integer or floating point return number, otherwise return `0`.
    real_type get_number() const ;
    
    template <typename T>
    T& get_pointer() { ck("get_pointer",GC_POINTER) ; return *(T*)(&(_data.p)) ; }

    template <typename T>
    T get_pointer() const { ck("get_pointer",GC_POINTER) ; return static_cast<T>(_data.p) ; }
    //!< Return the stored generic pointer (if fails issue an error).

    bool_type       & get_bool() ;
    bool_type const & get_bool() const ;
    //!< Return the stored boolean (if fails issue an error).

    int_type       &  get_int() ;
    int_type const &  get_int() const ;
    //!< Return the stored integer (if fails issue an error).

    real_type       & get_real() ;
    real_type const & get_real() const ;
    //!< Return the stored floating point (if fails issue an error).

    string_type       & get_string() ;
    string_type const & get_string() const ;
    //!< Return the stored string (if fails issue an error).
    //@}

    //! \name Access to vector type data
    //@{
    vector_type       & get_vector() ;
    vector_type const & get_vector() const ;
    //!< Return reference to a generic vector (if fails issue an error).

    vec_pointer_type       & get_vec_pointer() ;
    vec_pointer_type const & get_vec_pointer() const ;
    //!< Return reference to a vector of pointer (if fails issue an error).
    
    vec_bool_type       & get_vec_bool() ;
    vec_bool_type const & get_vec_bool() const ;
    //!< Return reference to a vector of booleans (if fails issue an error).

    vec_int_type       & get_vec_int() ;
    vec_int_type const & get_vec_int() const ;
    //!< Return reference to a vector of integers (if fails issue an error).

    vec_real_type       & get_vec_real() ;
    vec_real_type const & get_vec_real() const ;
    //!< Return reference to a vector of floating point number (if fails issue an error).

    vec_string_type       & get_vec_string() ;
    vec_string_type const & get_vec_string() const ;
    //!< Return reference to a vector of strings (if fails issue an error).
    //@}

    //! \name Access to element of vector type data
    //@{
    //! If `i`-th element of the vector is boolean, integer or floating point then return number, otherwise return `0`.
    real_type get_number( unsigned i ) const ;

    template <typename T>
    T& get_pointer( unsigned i ) { return (*this)[i].get_pointer<T>() ; }

    template <typename T>
    T get_pointer( unsigned i ) const { return (*this)[i].get_pointer<T>() ; }
    //!< Return `i`-th generic pointer (if fails issue an error).

    bool_type get_bool( unsigned i ) ;
    bool_type get_bool( unsigned i ) const ;
    //!< Return `i`-th boolean (if fails issue an error).
    
    int_type       & get_int( unsigned i ) ;
    int_type const & get_int( unsigned i ) const ;
    //!< Return `i`-th integer (if fails issue an error).

    real_type       & get_real( unsigned i ) ;
    real_type const & get_real( unsigned i ) const ;
    //!< Return `i`-th floating point number (if fails issue an error).

    string_type       & get_string( unsigned i ) ;
    string_type const & get_string( unsigned i ) const ;
    //!< Return `i`-th string (if fails issue an error).

    GenericContainer       & get_gc( unsigned i )       { return (*this)[i] ; }
    GenericContainer const & get_gc( unsigned i ) const { return (*this)[i] ; }
    //!< Return `i`-th `GenericContainer` of a generic vector (if fails issue an error).
    //@}

    //! \name Access to map type element
    //@{
    map_type       & get_map() ;
    map_type const & get_map() const ;
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
    { this -> set_int(int(a))  ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( unsigned long a )
    { this -> set_int(int(a))  ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( int a )
    { this -> set_int(a)  ; return * this ; }

    //! Assign an integer to the generic container.
    GenericContainer & operator = ( long a )
    { this -> set_int(int(a))  ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( float a )
    { this -> set_real(a) ; return * this ; }

    //! Assign a floating point number to the generic container.
    GenericContainer & operator = ( double a )
    { this -> set_real(a) ; return * this ; }

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

    //! If data contains a boolean or an integer it is promoted to a real.
    GenericContainer const & promote_to_real() ;

    //! If data contains vector of booleans it is promoted to a vector of integer.
    GenericContainer const & promote_to_vec_int() ;

    //! If data contains vector of booleans or integer it is promoted to a vector of real.
    GenericContainer const & promote_to_vec_real() ;

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
    void print( std::ostream &,
                std::string const & prefix = "",
                std::string const & indent = "    " ) const ;

    //! print the contents of the object in yaml syntax
    void to_yaml( std::ostream &, std::string const & prefix = "" ) const ;

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
    writeFormattedData( std::ostream & stream, char const delimiter = '\t' ) const ;

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

}

#endif

//
// eof: GenericContainer.hh
//
