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
\version   1.0.2
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
GenericContainer gc ;
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
GenericContainer::vec_int_type & v = gc["five"].get_vec_int();
cout << v[1] << '\n' ;
~~~~~~~~~~~~~

*/

#ifndef GENERIC_CONTAINER_HH
#define GENERIC_CONTAINER_HH

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

// if C++ < C++11 define nullptr
#if __cplusplus <= 199711L
  #include <cstddef>
  #define nullptr NULL
#endif

#ifndef ASSERT
  #if defined(__MINGW32__) || defined(__MINGW64__)
    #include <cstdlib>
  #endif
  #include <stdexcept>
  #define ASSERT(COND,MSG)                  \
    if ( !(COND) ) {                        \
      std::ostringstream ost ;              \
      ost << "On line: " << __LINE__        \
          << " file: " << __FILE__          \
          << '\n' << MSG << '\n' ;          \
      throw std::runtime_error(ost.str()) ; \
    }
#endif

/*!
  \brief

  `GenericContainer` is a class which permit to store eterogeneous data:

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

  enum TypeAllowed { GC_NOTYPE=0,
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
                     GC_MAP } ;

private:

  //! Data is stored in a union
  typedef union {
    pointer_type      p ;
    bool_type         b ;
    int_type          i ;
    real_type         r ;
    string_type       * s ;

    vec_pointer_type  * v_p ;
    vec_bool_type     * v_b ;
    vec_int_type      * v_i ;
    vec_real_type     * v_r ;
    vec_string_type   * v_s ;

    vector_type       * v ;
    map_type          * m ;

  } DataStorage ;

  DataStorage data      ; //!< The data stored in the class instance
  TypeAllowed data_type ; //!< The kind of data stored

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
  TypeAllowed get_type() const { return data_type ; }
  
  //! Return a string representing the type of data stored
  std::string get_type_name() const ;

  //! Print to stream the kind of data stored
  GenericContainer const & info( std::ostream & stream ) const ;

  //! If data is boolean, integer or floating point return number, otherwise return `0`.
  real_type get_number() const ;
  
  template <typename T>
  T * & get_pointer() { ck("get_pointer",GC_POINTER) ; return (T*)data.p ; }

  template <typename T>
  T const * get_pointer() const { ck("get_pointer",GC_POINTER) ; return (T const*)data.p ; }
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
  T * get_pointer( unsigned i ) { return (*this)[i].get_pointer<T>() ; }

  template <typename T>
  T const * get_pointer( unsigned i ) const
  { return (*this)[i].get_pointer<T>() ; }
  //!< Return `i`-th generic pointer (if fails issue an error).
  
  bool_type       & get_bool( unsigned i )       { return (*this)[i].get_bool() ; }
  bool_type const & get_bool( unsigned i ) const { return (*this)[i].get_bool() ; }
  //!< Return `i`-th boolean (if fails issue an error).
  
  int_type       & get_int( unsigned i )       { return (*this)[i].get_int() ; }
  int_type const & get_int( unsigned i ) const { return (*this)[i].get_int() ; }
  //!< Return `i`-th integer (if fails issue an error).

  real_type       & get_real( unsigned i )       { return (*this)[i].get_real() ; }
  real_type const & get_real( unsigned i ) const { return (*this)[i].get_real() ; }
  //!< Return `i`-th floating point number (if fails issue an error).

  string_type       & get_string( unsigned i )       { return (*this)[i].get_string() ; }
  string_type const & get_string( unsigned i ) const { return (*this)[i].get_string() ; }
  //!< Return `i`-th string (if fails issue an error).
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
  { this -> set_int(a)  ; return * this ; }

  //! Assign an integer to the generic container.
  GenericContainer & operator = ( int a )
  { this -> set_int(a)  ; return * this ; }

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

  //! Assign a generic container `a` to the generic container.
  GenericContainer const & operator = ( GenericContainer const & a ) ;
  //@}

  //! \name Initialize data by overloading constructor
  //@{
  
  //! Construct a generic container storing a boolean
  GenericContainer( bool a ) { *this = a ; }
  
  //! Construct a generic container storing an integer
  GenericContainer( unsigned a ) { *this = a ; }
  
  //! Construct a generic container storing an integer
  GenericContainer( int a ) { *this = a ; }
  
  //! Construct a generic container storing a floating point number
  GenericContainer( float a ) { *this = a ; }
  
  //! Construct a generic container storing a floating point number
  GenericContainer( double a ) { *this = a ; }
  
  //! Construct a generic container storing a string
  GenericContainer( char const a[] ) { *this = a ; }

  //! Construct a generic container storing a string
  GenericContainer( std::string const & a ) { *this = a ; }

  //! Construct a generic container copying container `gc`
  GenericContainer( GenericContainer const & gc ) { *this = gc ; }
  //@}

  //! \name I/O for `GenericContainer` objects
  //@{

  //! print the contents of the object in a human readable way
  void print( std::ostream &, std::string const & prefix = "" ) const ;

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

} ;

#endif

//
// eof: GenericContainer.hh
//
