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

#ifndef GENERIC_CONTAINER_HH
#define GENERIC_CONTAINER_HH

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wpoison-system-directories"
#endif

#include <iostream>
#include <string>
#include <complex>
#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <stdexcept>

#include "GenericContainerConfig.hh"

#ifndef GC_DO_ERROR
  #define GC_DO_ERROR(MSG) {                          \
    std::ostringstream ost;                           \
    ost << "in GenericContainer: " << MSG << '\n';    \
    GenericContainer::exception( ost.str().c_str() ); \
  }
#endif

#ifndef GC_ASSERT
  #define GC_ASSERT(COND,MSG) if ( !(COND) ) GC_DO_ERROR(MSG)
#endif

#ifndef GC_WARNING
  #define GC_WARNING(COND,MSG)                                        \
    if ( !(COND) ) {                                                  \
      std::cout << "On line: " << __LINE__                            \
                << " file: " << __FILE__                              \
                << " in GenericContainer\nWARNING: " << MSG << '\n';  \
    }
#endif

#ifdef __GNUC__
  #define GC_NO_RETURN __attribute__ ((noreturn))
#else
  #define GC_NO_RETURN
#endif

namespace GC_namespace {

  typedef std::basic_ostream<char> ostream_type;

  #ifndef _MSC_VER
  using std::int32_t;
  using std::int64_t;
  using std::uint32_t;
  using std::uint64_t;
  #endif

  class GenericContainer;

  typedef void*                   pointer_type; //!< generic pointer type
  typedef bool                    bool_type;    //!< boolean type data
  typedef int32_t                 int_type;     //!< integer type data
  typedef int64_t                 long_type;    //!< long integer type data
  typedef double                  real_type;    //!< floating point type data
  typedef std::complex<real_type> complex_type; //!< complex floating point type data
  typedef std::string             string_type;  //!< string type data

  typedef std::vector<pointer_type> vec_pointer_type; //!< vector of generic pointer
  typedef std::vector<bool_type>    vec_bool_type;    //!< vector of boolean
  typedef std::vector<int_type>     vec_int_type;     //!< vector of integer
  typedef std::vector<long_type>    vec_long_type;    //!< vector of long integer
  typedef std::vector<real_type>    vec_real_type;    //!< vector of floating point
  typedef std::vector<complex_type> vec_complex_type; //!< vector of complex floating point
  typedef std::vector<string_type>  vec_string_type;  //!< vector of strings

  typedef std::vector<GenericContainer>          vector_type; //!< vector of `GenericContainer`
  typedef std::map<string_type,GenericContainer> map_type;    //!< associative array of `GenericContainer`

  // ---------------------------------------------------------------------------

  typedef uint32_t uint_type;  //!< integer type data
  typedef uint64_t ulong_type; //!< long integer type data

  typedef std::vector<uint_type>  vec_uint_type;  //!< vector of unsigned integer
  typedef std::vector<ulong_type> vec_ulong_type; //!< vector of unsigned integer

  // ---------------------------------------------------------------------------
  template <typename TYPE>
  class mat_type : public std::vector<TYPE> {
    unsigned m_numRows;
    unsigned m_numCols;
    typedef typename std::vector<TYPE>::size_type size_type;
  public:

    mat_type()
    : m_numRows(0)
    , m_numCols(0)
    {}

    //!
    //! Create a matrix of size `nr x nc`
    //!
    mat_type( unsigned nr, unsigned nc )
    : m_numRows(nr)
    , m_numCols(nc)
    { std::vector<TYPE>::resize(size_type(nr*nc)); }

    //!
    //! Resize the matrix to size `nr x nc`
    //!
    void
    resize( unsigned nr, unsigned nc ) {
      m_numRows = nr;
      m_numCols = nc;
      std::vector<TYPE>::resize(size_type(nr*nc));
    }

    //!
    //! Get the `nc`-th column of the matrix and store to vector `C`.
    //!
    void getColumn( unsigned nc, std::vector<TYPE> & C ) const;

    //!
    //! Get the `nr`-th row of the matrix and store to vector `R`.
    //!
    void getRow( unsigned nr, std::vector<TYPE> & R ) const;

    //!
    //! Get the `nc`-th column of the matrix and store to memory pointed by `C`.
    //!
    void getColumn( unsigned nc, TYPE * C ) const;

    //!
    //! Get the `nr`-th row of the matrix and store to memory pointed by `R`.
    //!
    void getRow( unsigned nr, TYPE * R ) const;

    //!
    //! Get the number of rows of the matrix
    //!
    unsigned numRows() const { return m_numRows; }

    //!
    //! Get the number of columns of the matrix
    //!
    unsigned numCols() const { return m_numCols; }

    //!
    //! Access to the element `A(i,j)` og the matrix
    //!
    TYPE const & operator () ( unsigned i, unsigned j ) const;

    //!
    //! Access to the element `A(i,j)` og the matrix
    //!
    TYPE       & operator () ( unsigned i, unsigned j );

    void info( ostream_type & stream ) const;

    std::string
    info() const {
      std::ostringstream ostr;
      this->info(ostr);
      return ostr.str();
    }

    //!
    //! Access the memory block used by the matrix object
    //!
    TYPE       * data()       { return &std::vector<TYPE>::front(); }

    //!
    //! Access the memory block used by the matrix object
    //!
    TYPE const * data() const { return &std::vector<TYPE>::front(); }
  };

  // ---------------------------------------------------------------------------

  #ifndef GENERIC_CONTAINER_ON_WINDOWS
  extern template class mat_type<int_type>;
  extern template class mat_type<long_type>;
  extern template class mat_type<real_type>;
  extern template class mat_type<complex_type>;
  #endif

  typedef mat_type<int_type>     mat_int_type;
  typedef mat_type<long_type>    mat_long_type;
  typedef mat_type<real_type>    mat_real_type;
  typedef mat_type<complex_type> mat_complex_type;

  // ---------------------------------------------------------------------------
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename TYPE>
  ostream_type & operator << ( ostream_type & s, std::vector<TYPE> const & v );

  template <typename TYPE>
  ostream_type & operator << ( ostream_type & s, mat_type<TYPE> const & );
  #endif
  // ---------------------------------------------------------------------------

  //!
  //! Type allowed for the `GenericContainer`
  //!
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
    GC_VEC_LONG,
    GC_VEC_REAL,
    GC_VEC_COMPLEX,
    GC_VEC_STRING,

    // matrix type
    GC_MAT_INTEGER,
    GC_MAT_LONG,
    GC_MAT_REAL,
    GC_MAT_COMPLEX,

    // complex type
    GC_VECTOR,
    GC_MAP
  };

  //!
  //! `GenericContainer` is a class which permit to store eterogeneous data:
  //!
  //! - pointer
  //! - boolean
  //! - integer
  //! - long integer
  //! - floating point
  //! - complex floating point
  //! - string
  //! - vector of pointer
  //! - vector of boolean
  //! - vector of integer
  //! - vector of floating point
  //! - vector of complex floating point
  //! - matrix of floating point
  //! - matrix of complex floating point
  //! - vector of string
  //!
  //! in addition to this data type the following two container are added
  //!
  //! - vector of `GenericContainer`
  //! - map of `GenericContainer`
  //!
  class GenericContainer {
  public:
    // import type
    typedef GC_namespace::pointer_type     pointer_type;
    typedef GC_namespace::bool_type        bool_type;
    typedef GC_namespace::int_type         int_type;
    typedef GC_namespace::uint_type        uint_type;
    typedef GC_namespace::long_type        long_type;
    typedef GC_namespace::ulong_type       ulong_type;
    typedef GC_namespace::real_type        real_type;
    typedef GC_namespace::complex_type     complex_type;
    typedef GC_namespace::string_type      string_type;
    typedef GC_namespace::vec_pointer_type vec_pointer_type;
    typedef GC_namespace::vec_bool_type    vec_bool_type;
    typedef GC_namespace::vec_int_type     vec_int_type;
    typedef GC_namespace::vec_uint_type    vec_uint_type;
    typedef GC_namespace::vec_long_type    vec_long_type;
    typedef GC_namespace::vec_ulong_type   vec_ulong_type;
    typedef GC_namespace::vec_real_type    vec_real_type;
    typedef GC_namespace::vec_complex_type vec_complex_type;
    typedef GC_namespace::vec_string_type  vec_string_type;
    typedef GC_namespace::vector_type      vector_type;
    typedef GC_namespace::map_type         map_type;
    typedef GC_namespace::mat_int_type     mat_int_type;
    typedef GC_namespace::mat_long_type    mat_long_type;
    typedef GC_namespace::mat_real_type    mat_real_type;
    typedef GC_namespace::mat_complex_type mat_complex_type;

  private:

    //!
    //! Data is stored in a union
    //!
    typedef union {
      pointer_type     p;
      bool_type        b;
      int_type         i;
      long_type        l;
      real_type        r;
      complex_type     * c;
      string_type      * s;

      vec_pointer_type * v_p;
      vec_bool_type    * v_b;
      vec_int_type     * v_i;
      vec_long_type    * v_l;
      vec_real_type    * v_r;
      vec_complex_type * v_c;
      vec_string_type  * v_s;

      mat_int_type     * m_i;
      mat_long_type    * m_l;
      mat_real_type    * m_r;
      mat_complex_type * m_c;

      vector_type      * v;
      map_type         * m;

    } DataStorage;

    DataStorage m_data;      //!< The data stored in the class instance
    TypeAllowed m_data_type; //!< The kind of data stored

    void allocate_string();
    void allocate_complex();

    void allocate_vec_pointer( unsigned sz );
    void allocate_vec_bool( unsigned sz );
    void allocate_vec_int( unsigned sz );
    void allocate_vec_long( unsigned sz );
    void allocate_vec_real( unsigned sz );
    void allocate_vec_complex( unsigned sz );
    void allocate_mat_int( unsigned nr, unsigned nc );
    void allocate_mat_long( unsigned nr, unsigned nc );
    void allocate_mat_real( unsigned nr, unsigned nc );
    void allocate_mat_complex( unsigned nr, unsigned nc );
    void allocate_vec_string( unsigned sz );

    void allocate_vector( unsigned sz );
    void allocate_map();

    void ck( char const *, TypeAllowed ) const;
    int  ck( TypeAllowed ) const;
    void ck_or_set( char const *, TypeAllowed );

    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    bool simple_data()     const;
    bool simple_vec_data() const;
    #else
    bool simple_data()     const { return m_data_type <= GC_STRING; }
    bool simple_vec_data() const { return m_data_type <= GC_VEC_STRING; }
    #endif

  public:

    //!
    //! build an instance of `GenericContainer` with empty data
    //!
    GenericContainer();

    //!
    //! destroy the instance of `GenericContainer`
    //!
    ~GenericContainer() { clear(); }

    //!
    //! free memory of the data stored in `GenericContainer`,
    //! data type become `NOTYPE`
    //!
    void clear();

    //!
    //! free memory of the data stored in `GenericContainer`,
    //! data type must be 'MAP'
    //!
    void erase( char const * name );

    //! \name Initialize simple data
    ///@{
    //!
    //! Set data to `pointer_type` initialize and
    //! return a reference to the data
    //!
    pointer_type & set_pointer( pointer_type value );

    //!
    //! Free pointer to memory pointed, set data to `NO TYPE`
    //! initialize and return a reference to the data
    //!
    GenericContainer & free_pointer();

    //!
    //! Set data to `bool_type` initialize and return
    //! a reference to the data
    //!
    bool_type & set_bool( bool_type value );

    //!
    //! Set data to `int_type` initialize and return
    //! a reference to the data
    //!
    int_type & set_int( int_type value );

    //!
    //! Set data to `int_type` initialize and
    //! return a reference to the data
    //!
    long_type & set_long( long_type value );

    //!
    //! Set data to `real_type` initialize and
    //! return a reference to the data
    //!
    real_type & set_real( real_type value );

    //!
    //! Set data to `complex_type` initialize and
    //! return a reference to the data
    //!
    complex_type & set_complex( complex_type & value );

    //!
    //! Set data to `complex_type` initialize and
    //! return a reference to the data
    //!
    complex_type & set_complex( real_type r, real_type i );

    //!
    //! Set data to `string_type`, allocate and initialize.
    //! Return a reference to the data
    //!
    string_type & set_string( string_type const & value );
    ///@}

    //! \name Initialize vector data
    ///@{
    //!
    //! Set data to `vec_pointer_type`, allocate and initialize.
    //!
    //! Return a reference to vector of pointer.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of pointers that will be allocated
    //! \return the reference of the internally allocated vector of pointers
    //!
    vec_pointer_type & set_vec_pointer( unsigned sz = 0 );

    //!
    //! Set data to `vec_pointer_type`, allocate and initialize.
    //!
    //! Return a reference to vector of pointer.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of pointers  used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of pointers
    //!
    vec_pointer_type & set_vec_pointer( vec_pointer_type const & v );

    //!
    //! Set data to `vec_bool_type`, allocate and initialize.
    //!
    //! Return a reference to vector of booleans.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of bools that will be allocated
    //! \return the reference of the internally allocated vector of bools
    //!
    vec_bool_type & set_vec_bool( unsigned sz = 0 );

    //! 
    //! Set data to `vec_bool_type`, allocate and initialize.
    //!
    //! Return a reference to vector of `bool`.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of pointers  used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of pointers
    //!
    vec_bool_type & set_vec_bool( vec_bool_type const & v );

    //!
    //! Set data to `vec_int_type`, allocate and initialize.
    //!
    //! Return a reference to vector of integers.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of integers that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_int_type & set_vec_int( unsigned sz = 0 );

    //!
    //! Set data to `vec_int_type`, allocate and initialize.
    //!
    //! Return a reference to vector of integer.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of integer used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_int_type & set_vec_int( vec_int_type const & v );

    //!
    //! Set data to `vec_int_type`, allocate and initialize.
    //!
    //! Return a reference to vector of integers.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of integers that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_long_type & set_vec_long( unsigned sz = 0 );

    //!
    //! Set data to `vec_int_type`, allocate and initialize.
    //!
    //! Return a reference to vector of integer.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of long integers used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of long integers
    //!
    vec_long_type & set_vec_long( vec_long_type const & v );

    //!
    //! Set data to `vec_real_type`, allocate and initialize.
    //!
    //! Return a reference to vector of `real_type` numbers.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of integers that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_real_type & set_vec_real( unsigned sz = 0 );

    //!
    //! Set data to `vec_real_type`, allocate and initialize.
    //!
    //! Return a reference to vector of `real_type` number.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of reals  used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of reals
    //!
    vec_real_type & set_vec_real( vec_real_type const & v );

    //!
    //! Set data to `vec_complex_type`, allocate and initialize.
    //! Return a reference to vector of complex numbers.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of integers that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_complex_type & set_vec_complex( unsigned sz = 0 );

    //! 
    //! Set data to `vec_complex_type`, allocate and initialize.
    //!
    //! Return a reference to vector of complex number.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of complex used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of complex
    //!
    vec_complex_type & set_vec_complex( vec_complex_type const & v );

    //!
    //! Set data to `vec_string_type`, allocate and initialize.
    //!
    //! Return a reference to vector of strings.
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz dimension of the vector of integers that will be allocated
    //! \return the reference of the internally allocated vector of integers
    //!
    vec_string_type & set_vec_string( unsigned sz = 0 );

    //!
    //! Set data to `vec_string_type`, allocate and initialize.
    //!
    //! Return a reference to vector of strings.
    //! Copy the data from vector `v`.
    //!
    //! \param[in] v  vector of strings used to initialize that will be allocated
    //! \return the reference of the internally allocated vector of strings
    //!
    vec_string_type & set_vec_string( vec_string_type const & v );

    //!
    //! Set data to `mat_int_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
    //!
    //! \param[in] nr number of rows
    //! \param[in] nc number of columns
    //! \return the reference of the internally allocated matrix of integers
    //!
    mat_int_type & set_mat_int( unsigned nr = 0, unsigned nc = 0 );

    //!
    //! Set data to `mat_int_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! Copy the data from matrix `m`.
    //!
    //! \param[in] m matrix of integers used to initialize data
    //! \return the reference of the internally allocated matrix of integers
    //!
    mat_int_type & set_mat_int( mat_int_type const & m );

    //!
    //! Set data to `mat_long_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
    //!
    //! \param[in] nr number of rows
    //! \param[in] nc number of columns
    //! \return the reference of the internally allocated matrix of long integers
    //!
    mat_long_type & set_mat_long( unsigned nr = 0, unsigned nc = 0 );

    //!
    //! Set data to `mat_long_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! Copy the data from matrix `m`.
    //!
    //! \param[in] m matrix of long used to initialize data
    //! \return the reference of the internally allocated matrix of long
    //!
    mat_long_type & set_mat_long( mat_long_type const & m );

    //!
    //! Set data to `mat_real_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
    //!
    //! \param[in] nr number of rows
    //! \param[in] nc number of columns
    //! \return the reference of the internally allocated matrix of reals
    //!
    mat_real_type & set_mat_real( unsigned nr = 0, unsigned nc = 0 );

    //!
    //! Set data to `mat_real_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `real_type`.
    //! Copy the data from matrix `m`.
    //!
    //! \param[in] m matrix of reals used to initialize data
    //! \return the reference of the internally allocated matrix of reals
    //!
    mat_real_type & set_mat_real( mat_real_type const & m );

    //!
    //! Set data to `mat_complex_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `complex_type`.
    //! If `nr` > 0 and `nc` > 0 then the matrix is allocated to size `nr` x `nc`.
    //!
    //! \param[in] nr number of rows
    //! \param[in] nc number of columns
    //! \return the reference of the internally allocated matrix of complex
    //!
    mat_complex_type & set_mat_complex( unsigned nr = 0, unsigned nc = 0 );

    //!
    //! Set data to `mat_complex_type`, allocate and initialize.
    //!
    //! Return a reference to a matrix of `complex_type`.
    //! Copy the data from matrix `m`.
    //!
    //! \param[in] m matrix of integers used to initialize data
    //! \return the reference of the internally allocated matrix of integers
    //!
    mat_complex_type & set_mat_complex( mat_complex_type const & m );

    //!
    //! push boolean if data is `vec_bool_type` or `vector_type`
    //!
    void push_bool( bool );

    //!
    //! push integer if data is `vec_int_type` or `vector_type`
    //!
    void push_int( int_type );

    //!
    //! push integer if data is `vec_int_type` or `vector_type`
    //!
    void push_long( long_type );

    //!
    //! push real if data is `vec_real_type` or `vector_type`
    //!
    void push_real( real_type );

    //!
    //! push complex if data is `vec_complex_type` or `vector_type`
    //!
    void push_complex( complex_type & );

    //!
    //! push complex if data is `vec_complex_type` or `vector_type`
    //!
    void push_complex( real_type re, real_type im );

    //!
    //! push complex if data is `vec_string_type` or `vector_type`
    //!
    void push_string( string_type const & );

    ///@}

    //! \name Initialize generic data
    ///@{
    //!
    //! Set data to `vector_type`, allocate an empty generic
    //! vector and return a reference to it.
    //!
    //! If `sz` > 0 then the vector is allocated to size `sz`.
    //!
    //! \param[in] sz the size of the allocated vector
    //!
    vector_type & set_vector( unsigned sz = 0 );

    //!
    //! Set data to `map_type`, allocate an empty generic map
    //! and return a reference to it.
    //!
    map_type & set_map();
    ///@}

    //! \name Access to a single element
    ///@{

    //! Return an integer representing the type of data stored
    //!
    //! Integer to data type map
    //! ------------------------
    //!
    //! -   No data stored (return 0)
    //! 1.  `pointer_type`
    //! 2.  `bool_type`
    //! 3.  `int_type`
    //! 4.  `long_type`
    //! 5.  `real_type`
    //! 6.  `complex_type`
    //! 7.  `string_data`
    //! 8.  `vec_pointer_type`
    //! 9.  `vec_bool_type`
    //! 10. `vec_int_type`
    //! 11. `vec_long_type`
    //! 12. `vec_real_type`
    //! 13. `vec_complex_type`
    //! 14. `vec_string_type`
    //! 15. `mat_int_type`
    //! 16. `mat_long_type`
    //! 17. `mat_real_type`
    //! 18. `mat_complex_type`
    //! 19. `vector_type`
    //! 20. `map_type`
    //!
    //! \return the type of the internally stored data 
    //!
    TypeAllowed get_type() const { return m_data_type; }

    //!
    //! Return a string pointer representing the type of data stored
    //!
    char const * get_type_name() const;

    //!
    //! Print to stream the kind of data stored
    //!
    GenericContainer const & info( ostream_type & stream ) const;

    //!
    //! Print to string the kind of data stored
    //!
    std::string
    info() const {
      std::ostringstream sstr;
      this->info(sstr);
      return sstr.str();
    }

    //!
    //! Return the number of the elements of the first level of the generic container
    //! 1 for single element, the size of vector or map, or 0
    //!
    unsigned get_num_elements() const;

    //!
    //! Return the number of rows ot the internally stored matrix
    //!
    unsigned get_numRows() const;

    //!
    //! Return the number of columns ot the internally stored matrix
    //!
    unsigned get_numCols() const;

    //!
    //! If data is boolean, integer or `real_type`
    //! return number, otherwise return `0`.
    //!
    real_type get_number() const;

    //!
    //! If data is boolean, integer, `real_type` or
    //! `complex_type` return number, otherwise return `0`.
    //!
    complex_type get_complex_number() const;

    //!
    //! If data is boolean, integer, `real_type` or
    //! `complex_type` return number, otherwise return `0`.
    //!
    void get_complex_number( real_type & re, real_type & im ) const;

    //!
    //! Return the stored data as a pointer
    //!
    void * get_pvoid( char const * msg = nullptr ) const;

    //!
    //! Return the stored data as a double pointer
    //!
    void ** get_ppvoid( char const * msg = nullptr ) const;

    //!
    //! Return the stored data as a pointer to const integer
    //!
    int_type const * get_int_pointer() const;

    //!
    //! Return the stored data as a pointer to integer
    //!
    int_type * get_int_pointer();

    //!
    //! Return the stored data as a pointer to const long
    //!
    long_type const * get_long_pointer() const;

    //!
    //! Return the stored data as a pointer to long
    //!
    long_type * get_long_pointer();

    //!
    //! Return the stored data as a pointer to const `real_type`
    //!
    real_type const * get_real_pointer() const;

    //!
    //! Return the stored data as a pointer to `real_type`
    //!
    real_type * get_real_pointer();

    //!
    //! Return the stored data as a pointer to const `complex_type`
    //!
    complex_type const * get_complex_pointer() const;

    //!
    //! Return the stored data as a pointer to `complex_type`
    //!
    complex_type * get_complex_pointer();

    //!
    //! Get the stored value
    //!
    //!
    //! \param[out] v    a copy of the stored value
    //! \param[in]  msg  message of error in case of failure
    //!
    template <typename T>
    void
    get_value( T & v, char const * msg = "" ) const;

    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_POINTER); return *reinterpret_cast<T*>(get_ppvoid()); }

    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_POINTER); return reinterpret_cast<T>(get_pvoid()); }
    #else
    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_POINTER); return *reinterpret_cast<T*>(&(m_data.p)); }

    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_POINTER); return reinterpret_cast<T>(m_data.p); }
    #endif

    //!
    //! Get the stored value in the map as boolean
    //!
    //! \param[in] key key of the map to be selected
    //! \return the boolean stored in the container
    //!
    bool_type
    get_map_bool( char const * key ) const {
      bool_type ret = false;
      if ( exists(key) ) ret = (*this)(key).get_bool();
      return ret;
    }

    //!
    //! Get the stored value as a boolean
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the boolean stored in the container
    //!
    bool_type & get_bool( char const * msg = nullptr );

    //!
    //! Get the stored value as a const boolean
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the boolean stored in the container
    //!
    bool_type const & get_bool( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as an integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the boolean stored in the container
    //!
    int_type & get_int( char const * msg = nullptr );

    //!
    //! Get the stored value as a const integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the integer stored in the container
    //!
    int_type const & get_int( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a long integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the long integer stored in the container
    //!
    long_type & get_long( char const * msg = nullptr );

    //!
    //! Get the stored value as a long integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the long integer stored in the container
    //!
    long_type const & get_long( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as an integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as an integer
    //!
    int_type get_as_int( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as an unsigned
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as an unsigned
    //!
    uint_type get_as_uint( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a long integer
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a long integer
    //!
    long_type get_as_long( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as an unsigned long
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as an unsigned long
    //!
    ulong_type get_as_ulong( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as `real_type`
    //!
    real_type & get_real( char const * msg = nullptr );

    //!
    //! Get the stored value as a const `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const `real_type`
    //!
    real_type const & get_real( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a `complex_type`
    //!
    complex_type & get_complex( char const * msg = nullptr );

    //!
    //! Get the stored value as a const `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const `complex_type`
    //!
    complex_type const & get_complex( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a string
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a string
    //!
    string_type & get_string( char const * msg = nullptr );

    //!
    //! Get the stored value as a const string
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const string
    //!
    string_type const & get_string( char const * msg = nullptr ) const;
    ///@}

    //! \name Access to vector type data
    ///@{
    //!
    //! Get the stored value as a vector of `GenericoContainer`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of `GenericoContainer`
    //!
    vector_type & get_vector( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of `GenericoContainer`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of `GenericoContainer`
    //!
    vector_type const & get_vector( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of pointers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of pointers
    //!
    vec_pointer_type & get_vec_pointer( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of pointers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of pointers
    //!
    vec_pointer_type const & get_vec_pointer( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of booleans
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of booleans
    //!
    vec_bool_type & get_vec_bool( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of booleans
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of booleans
    //!
    vec_bool_type const & get_vec_bool( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of integers
    //!
    vec_int_type & get_vec_int( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of integers
    //!
    vec_int_type const & get_vec_int( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of long integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of long integers
    //!
    vec_long_type & get_vec_long( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of long integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of long integers
    //!
    vec_long_type const & get_vec_long( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of `real_type`
    //!
    vec_real_type & get_vec_real( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of `real_type`
    //!
    vec_real_type const & get_vec_real( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of `complex_type`
    //!
    vec_complex_type & get_vec_complex( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of `complex_type`
    //!
    vec_complex_type const & get_vec_complex( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a matrix of integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a matrix of integers
    //!
    mat_int_type & get_mat_int( char const * msg = nullptr );

    //!
    //! Get the stored value as a const matrix of integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const matrix of integers
    //!
    mat_int_type const & get_mat_int( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a matrix of long integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a matrix of long integers
    //!
    mat_long_type & get_mat_long( char const * msg = nullptr );

    //!
    //! Get the stored value as a const matrix of long integers
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const matrix of long integers
    //!
    mat_long_type const & get_mat_long( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a matrix of `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a matrix of `real_type`
    //!
    mat_real_type & get_mat_real( char const * msg = nullptr );

    //!
    //! Get the stored value as a const matrix of `real_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const matrix of `real_type`
    //!
    mat_real_type const & get_mat_real( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a matrix of `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a matrix of `complex_type`
    //!
    mat_complex_type & get_mat_complex( char const * msg = nullptr );

    //!
    //! Get the stored value as a const matrix of `complex_type`
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const matrix of `complex_type`
    //!
    mat_complex_type const & get_mat_complex( char const * msg = nullptr ) const;

    //!
    //! Get the stored value as a vector of string
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a vector of string
    //!
    vec_string_type & get_vec_string( char const * msg = nullptr );

    //!
    //! Get the stored value as a const vector of string
    //!
    //! \param[in] msg message of error in case of failure
    //! \return the data stored in the container as a const vector of string
    //!
    vec_string_type const & get_vec_string( char const * msg = nullptr ) const;

    ///@}

    //! \name Access to vector type data and convert
    ///@{

    //!
    //! Copy internal data a vector of integers
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_int( vec_int_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of unsigned integer
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_uint( vec_uint_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of long integer
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_long( vec_long_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of unsigned long integer
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_ulong( vec_ulong_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of `real_type`
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_real( vec_real_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of `complex_type`
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_complex( vec_complex_type & v, char const * msg = "" ) const;

    //!
    //! Copy internal data a vector of strings
    //!
    //! \param[out] v   vector to store the data
    //! \param[in]  msg message of error in case of failure
    //!
    void copyto_vec_string( vec_string_type & v, char const * msg = "" ) const;

    ///@}

    //! \name Access to element of vector type data

    ///@{
    //!
    //! If `i`-th element of the vector is boolean, 
    //! integer or floating point then return number, otherwise return `0`.
    //!
    //! \param[in] i the position of the element in the vector
    //! \return the value as `real_type`
    //!
    real_type get_number_at( unsigned i ) const;

    //!
    //! If `i`-th element of the vector is convertible to
    //! compex return number, otherwise return `0`.
    //!
    //! \param[in] i the position of the element in the vector
    //! \return the value as `complex_type`
    //!
    complex_type get_complex_number_at( unsigned i ) const;

    //!
    //! If `i`-th element of the vector is convertible to
    //! compex return number, otherwise return `0`.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[out] re real part of the number
    //! \param[out] im imaginary part of the number
    //!
    void get_complex_number_at( unsigned i, real_type & re, real_type & im ) const;

    //!
    //! Get the `i`-th pointer of the vector of pointers.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the reference to the pointer
    //!
    template <typename T>
    T& get_pointer_at( unsigned i )
    { return (*this)[i].get_pointer<T>(); }

    //!
    //! Get the `i`-th pointer of the vector of pointers.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored pointer
    //!
    template <typename T>
    T get_pointer_at( unsigned i ) const
    { return (*this)[i].get_pointer<T>(); }
    //!< Return `i`-th generic pointer (if fails issue an error).

    //!
    //! Get the `i`-th boolean of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored value
    //!
    bool_type get_bool_at( unsigned i );

    //!
    //! Get the `i`-th boolean of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    bool_type get_bool_at( unsigned i, char const * msg ) const;

    //!
    //! Get the `i`-th integer of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored value
    //!
    int_type & get_int_at( unsigned i );

    //!
    //! Get the `i`-th const integer of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    int_type const & get_int_at( unsigned i, char const * msg ) const;

    //!
    //! Get the `i`-th long integer of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored value
    //!
    long_type & get_long_at( unsigned i );

    //!
    //! Get the `i`-th const long integer of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    long_type const & get_long_at( unsigned i, char const * msg ) const;

    //!
    //! Get the `i`-th `real_type` of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored value
    //!
    real_type & get_real_at( unsigned i );

    //!
    //! Get the `i`-th const `real_type` of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    real_type const & get_real_at( unsigned i, char const * msg ) const;
 
    //!
    //! Get the `i`-th `complex_type` of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \return the stored value
    //!
    complex_type & get_complex_at( unsigned i );

    //!
    //! Get the `i`-th const `complex_type` of the stored data.
    //!
    //! \param[in]  i  the position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    complex_type const & get_complex_at( unsigned i, char const * msg ) const;

    //!
    //! Get the `i`-th integer of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \return the stored value
    //!
    int_type & get_int_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const integer of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    int_type const & get_int_at( unsigned i, unsigned j, char const * msg ) const;

    //!
    //! Get the `i`-th long integer of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \return the stored value
    //!
    long_type & get_long_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const long integer of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    long_type const & get_long_at( unsigned i, unsigned j, char const * msg ) const;

    //!
    //! Get the `i`-th `real_type` of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \return the stored value
    //!
    real_type & get_real_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const `real_type` of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    real_type const & get_real_at( unsigned i, unsigned j, char const * msg ) const;

    //!
    //! Get the `i`-th `complex_type` of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \return the stored value
    //!
    complex_type & get_complex_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const `complex_type` of the stored data.
    //!
    //! \param[in]  i  row position of the element in the matrix
    //! \param[in]  j  colun position of the element in the matrix
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    complex_type const & get_complex_at( unsigned i, unsigned j, char const * msg ) const;

    //!
    //! Get the `i`-th string of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \return the stored value
    //!
    string_type & get_string_at( unsigned i );

    //!
    //! Get the `i`-th const string of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    string_type const & get_string_at( unsigned i, char const * msg ) const;

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \return the stored value
    //!
    GenericContainer & get_gc_at( unsigned i );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    GenericContainer const & get_gc_at( unsigned i, char const * msg ) const;

    ///@}

    //! \name Access to map type element
    ///@{

    //!
    //! Get the stored data as a map of `GenericContainer`.
    //!
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    map_type & get_map( char const * msg = nullptr );

    //!
    //! Get the stored data as a const map of `GenericContainer`.
    //!
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    map_type const & get_map( char const * msg = nullptr ) const;

    ///@}

    //! \name Access using operators
    ///@{

    //!
    //! Get the `i`-th `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \return the stored value
    //!
    GenericContainer & operator [] ( unsigned i );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  position of the element in the vector
    //! \return the stored value
    //!
    GenericContainer const & operator [] ( unsigned i ) const;

    //!
    //! Get the `i`-th `GenericContainer` of the stored data.
    //!
    //! \param[in]  s  key string of the element in the map
    //! \return the stored value
    //!
    GenericContainer & operator [] ( std::string const & s );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  s  key string of the element in the map
    //! \return the stored value
    //!
    GenericContainer const & operator [] ( std::string const & s ) const;

    //!
    //! Get the `i`-th `GenericContainer` of the stored data.
    //!
    //! \param[in]  i   position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    GenericContainer & operator () ( unsigned i, char const * msg = nullptr );

    //!
    //! Get the `i`-th `GenericContainer` of the stored data.
    //!
    //! \param[in]  i   position of the element in the vector
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    GenericContainer const & operator () ( unsigned i, char const * msg = nullptr ) const;

    //!
    //! Get a `GenericContainer` in the stored data.
    //!
    //! \param[in]  s   key string of the element in the map
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    GenericContainer & operator () ( std::string const & s, char const * msg = nullptr );

    //!
    //! Get a const `GenericContainer` in the stored data.
    //!
    //! \param[in]  s   key string of the element in the map
    //! \param[in]  msg message of error in case of failure
    //! \return the stored value
    //!
    GenericContainer const & operator () ( std::string const & s, char const * msg = nullptr ) const;
    ///@}

    //! \name Initialize data using set command
    ///@{

    //! 
    //! Assign a boolean to the generic container.
    //! 
    //! \param[in] a boolean to be stored
    //! 
    void set( bool const & a ) { this->set_bool(a); }

    //! 
    //! Assign an integer to the generic container.
    //! 
    //! \param[in] a integer to be stored
    //! 
    void set( uint_type const & a ) { this->set_int(int_type(a)); }

    //! 
    //! Assign an integer to the generic container.
    //! 
    //! \param[in] a integer to be stored
    //! 
    void set( int_type const & a ) { this->set_int(a); }

    //! 
    //! Assign an unsigned integer to the generic container.
    //! 
    //! \param[in] a unsigned integer to be stored
    //! 
    void set( ulong_type const & a ) { this->set_long(long_type(a)); }

    //! 
    //! Assign a long integer to the generic container.
    //! 
    //! \param[in] a long integer to be stored
    //! 
    void set( long_type const & a ) { this->set_long(a); }

    //! 
    //! Assign a float to the generic container.
    //! 
    //! \param[in] a float to be stored
    //! 
    void set( float const & a ) { this->set_real(real_type(a)); }

    //! 
    //! Assign a double to the generic container.
    //! 
    //! \param[in] a double to be stored
    //! 
    void set( double const & a ) { this->set_real(real_type(a)); }

    //! 
    //! Assign a complex of float to the generic container.
    //! 
    //! \param[in] a complex of float to be stored
    //! 
    void set( std::complex<float> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); }

    //! 
    //! Assign a complex of double to the generic container.
    //! 
    //! \param[in] a complex of double to be stored
    //! 
    void set( std::complex<double> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); }

    //! 
    //! Assign a string to the generic container.
    //! 
    //! \param[in] a string to be stored
    //! 
    void set( char const * a ) { this->set_string(a); }

    //! 
    //! Assign a string to the generic container.
    //! 
    //! \param[in] a string to be stored
    //! 
    void set( std::string const & a ) { this->set_string(a); }

    //! 
    //! Assign a pointer to the generic container.
    //! 
    //! \param[in] a pointer to be stored
    //! 
    void set( pointer_type a ) { this->set_pointer(a); }

    ///@}

    //! \name Initialize data using operators
    //!
    //! The `=` operator is overloaded to initialize the 
    //! `GenericContainer` on its left side.
    //!
    ///@{

    //! 
    //! Assign a boolean to the generic container.
    //! 
    //! \param[in] a boolean to be stored
    //! 
    GenericContainer & operator = ( bool const & a )
    { this->set_bool(a); return * this; }

    //! 
    //! Assign an integer to the generic container.
    //! 
    //! \param[in] a integer to be stored
    //! 
    GenericContainer & operator = ( uint_type const & a )
    { this->set_int(int_type(a)); return * this; }

    //! 
    //! Assign an integer to the generic container.
    //! 
    //! \param[in] a integer to be stored
    //! 
    GenericContainer & operator = ( int_type const & a )
    { this->set_int(a); return * this; }

    //! 
    //! Assign an unsigned long to the generic container.
    //! 
    //! \param[in] a unsigned long to be stored
    //! 
    GenericContainer & operator = ( ulong_type const & a )
    { this->set_long(long_type(a)); return * this; }

    //! 
    //! Assign a long integer to the generic container.
    //! 
    //! \param[in] a long integer to be stored
    //! 
    GenericContainer & operator = ( long_type const & a )
    { this->set_long(a); return * this; }

    //! 
    //! Assign a float to the generic container.
    //! 
    //! \param[in] a float to be stored
    //! 
    GenericContainer & operator = ( float const & a )
    { this->set_real(real_type(a)); return * this; }

    //! 
    //! Assign a double to the generic container.
    //! 
    //! \param[in] a double to be stored
    //! 
    GenericContainer & operator = ( double const & a )
    { this->set_real(real_type(a)); return * this; }

    //! 
    //! Assign a complex of float to the generic container.
    //! 
    //! \param[in] a complex of float to be stored
    //! 
    GenericContainer & operator = ( std::complex<float> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); return * this; }

    //! 
    //! Assign a complex of double to the generic container.
    //! 
    //! \param[in] a complex of double to be stored
    //! 
    GenericContainer & operator = ( std::complex<double> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); return * this; }

    //! 
    //! Assign a vector of bool to the generic container.
    //! 
    //! \param[in] a vector of bool to be stored
    //! 
    GenericContainer & operator = ( vec_bool_type const & a );

    //! 
    //! Assign a vector of integer to the generic container.
    //! 
    //! \param[in] a vector of integer to be stored
    //! 
    GenericContainer & operator = ( vec_int_type const & a );

    //! 
    //! Assign a vector of long integer to the generic container.
    //! 
    //! \param[in] a vector of long integer to be stored
    //! 
    GenericContainer & operator = ( vec_long_type const & a );

    //! 
    //! Assign a vector of `real_type` to the generic container.
    //! 
    //! \param[in] a vector of `real_type` to be stored
    //! 
    GenericContainer & operator = ( vec_real_type const & a );

    //! 
    //! Assign a vector of `complex_type` to the generic container.
    //! 
    //! \param[in] a vector of `complex_type` to be stored
    //! 
    GenericContainer & operator = ( vec_complex_type const & a );

    //! 
    //! Assign a vector of string to the generic container.
    //! 
    //! \param[in] a vector of string to be stored
    //! 
    GenericContainer & operator = ( vec_string_type const & a );

    //! 
    //! Assign a matrix of integer to the generic container.
    //! 
    //! \param[in] a matrix of integer to be stored
    //! 
    GenericContainer & operator = ( mat_int_type const & a );

    //! 
    //! Assign a matrix of long integer to the generic container.
    //! 
    //! \param[in] a matrix of long integer to be stored
    //! 
    GenericContainer & operator = ( mat_long_type const & a );

    //! 
    //! Assign a matrix of `real_type` to the generic container.
    //! 
    //! \param[in] a matrix of `real_type` to be stored
    //! 
    GenericContainer & operator = ( mat_real_type const & a );

    //! 
    //! Assign a matrix of `complex_type` to the generic container.
    //! 
    //! \param[in] a matrix of `complex_type` to be stored
    //! 
    GenericContainer & operator = ( mat_complex_type const & a );

    //! 
    //! Assign a string to the generic container.
    //! 
    //! \param[in] a string to be stored
    //! 
    GenericContainer & operator = ( char const * a )
    { this->set_string(a); return * this; }

    //! 
    //! Assign a string to the generic container.
    //! 
    //! \param[in] a string to be stored
    //! 
    GenericContainer & operator = ( std::string const & a )
    { this->set_string(a); return * this; }

    //! 
    //! Assign a pointer to the generic container.
    //! 
    //! \param[in] a pointer to be stored
    //! 
    GenericContainer & operator = ( pointer_type a )
    { this->set_pointer(a); return * this; }

    //! 
    //! Assign a `GenericContainer` to the generic container (deep copy).
    //! 
    //! \param[in] a `GenericContainer` to be stored
    //! 
    GenericContainer const & operator = ( GenericContainer const & a )
    { this->load( a ); return * this; }

    //! 
    //! Load a `GenericContainer` to the generic container (deep copy).
    //! 
    //! \param[in] a `GenericContainer` to be stored
    //! 
    void load( GenericContainer const & a );
    ///@}

    //! \name Promotion to a ``bigger'' data
    ///@{
    //!
    //! If data contains a boolean it is promoted to an integer.
    //!
    GenericContainer const & promote_to_int();

    //!
    //! If data contains a boolean it is promoted to an integer.
    //!
    GenericContainer const & promote_to_long();

    //!
    //! If data contains a boolean or an integer it is promoted to a real.
    //!
    GenericContainer const & promote_to_real();

    //!
    //! If data contains a boolean or an integer or floating point
    //! it is promoted to a complex floating point.
    //!
    GenericContainer const & promote_to_complex();

    //!
    //! If data contains vector of booleans it is promoted to a vector of integer.
    //!
    GenericContainer const & promote_to_vec_int();

    //!
    //! If data contains vector of booleans it is promoted to a vector of integer.
    //!
    GenericContainer const & promote_to_vec_long();

    //!
    //! If data contains vector of booleans or integer
    //! it is promoted to a vector of real.
    //!
    GenericContainer const & promote_to_vec_real();

    //!
    //! If data contains vector of booleans or integer
    //! or real it is promoted to a vector of complex.
    //!
    GenericContainer const & promote_to_vec_complex();

    //!
    //! If data contains vector of booleans, integer or
    //! real it is promoted to a matrix of integer.
    //!
    GenericContainer const & promote_to_mat_int();

    //!
    //! If data contains vector of booleans, integer or
    //! real it is promoted to a matrix of long integer.
    //!
    GenericContainer const & promote_to_mat_long();

    //!
    //! If data contains vector of booleans, integer or
    //! real it is promoted to a matrix of real.
    //!
    GenericContainer const & promote_to_mat_real();

    //!
    //! If data contains vector of booleans, integer or
    //! real or complex or matrix of real it is promoted to a matrix of complex.
    //!
    GenericContainer const & promote_to_mat_complex();

    //!
    //! If data contains vector of someting it is promoted
    //! to a vector of `GenericContainer`.
    //!
    GenericContainer const & promote_to_vector();

    ///@}

    //! \name Initialize data by overloading constructor
    ///@{

    //! 
    //! Construct a generic container storing a boolean
    //! \param[in] a initializer data
    //! 
    GenericContainer( bool const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //! 
    GenericContainer( uint_type const & a )
    : m_data_type(GC_NOTYPE) { *this = a; }

    //! 
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //! 
    GenericContainer( int_type const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //! 
    GenericContainer( ulong_type const & a )
    : m_data_type(GC_NOTYPE) { *this = a; }

    //! 
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //! 
    GenericContainer( long_type const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a floating point number
    //! \param[in] a initializer data
    //! 
    GenericContainer( float const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a floating point number
    //! \param[in] a initializer data
    //! 
    GenericContainer( double const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a complex floating point number
    //! \param[in] a initializer data
    //! 
    GenericContainer( std::complex<float> const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a complex floating point number
    //! \param[in] a initializer data
    //! 
    GenericContainer( std::complex<double> const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a string
    //! \param[in] a initializer data
    //! 
    GenericContainer( char const * a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a string
    //! \param[in] a initializer data
    //! 
    GenericContainer( std::string const & a )
    : m_data_type(GC_NOTYPE) { this->operator=(a); }

    //! 
    //! Construct a generic container storing a pointer
    //! \param[in] a initializer data
    //! 
    GenericContainer( pointer_type a )
    : m_data_type(GC_NOTYPE) { this->set_pointer(a); }

    //! 
    //! Construct a generic container copying container `gc`
    //! \param[in] gc initializer data
    //! 
    GenericContainer( GenericContainer const & gc )
    : m_data_type(GC_NOTYPE) { this->load(gc); }

    ///@}

    //! \name Utilities methods

    //!
    //! Check if string `s` is a key of the stored map (if fails issue an error).
    //! \param[in] s key to be checked
    //!
    bool exists( std::string const & s ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, bool & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, int_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, uint_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, long_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, ulong_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, real_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, complex_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( char const * field, string_type & value ) const;

    ///@}

    //! \name I/O for GenericContainer objects
    ///@{

    //!
    //! Print the contents of the object in a human readable way (without data)
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    void
    print_content_types(
      ostream_type      & stream,
      std::string const & prefix = "",
      std::string const & indent = "    "
    ) const;

    //!
    //! Dump the contents of the object in a human readable way
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    void
    dump(
      ostream_type      & stream,
      std::string const & prefix = "",
      std::string const & indent = "    "
    ) const;

    //!
    //! Dump the contents of the object in a human readable way (aliad of dump)
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    void
    print(
      ostream_type      & stream,
      std::string const & prefix = "",
      std::string const & indent = "    "
    ) const {
      this->dump( stream, prefix, indent );
    }

    //!
    //! Dump the contents of the object in a human readable way (aliad of dump)
    //!
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    std::string
    print(
      std::string const & prefix = "",
      std::string const & indent = "    "
    ) const {
      std::ostringstream ostr;
      this->print(ostr,prefix,indent);
      return ostr.str();
    }

    //!
    //! Print the contents of the object in yaml syntax
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //!
    void to_yaml( ostream_type & stream, std::string const & prefix = "" ) const;

    //!
    //! Write `GenericContainer` as regular formatted data to stream.
    //!
    //! `GenericContainer` must be a map which contains the fields:
    //!
    //! - "headers" this element must be a `vec_string_type` which contains
    //!             the strings of the headers of the columns of the data
    //!
    //! - "data"    this element must be a `vector_type` which contais the
    //!             vectors which are the columns of the data to be saved.
    //!             Each column can be of type
    //!
    //!             1. `vec_bool_type`
    //!             2. `vec_int_type`
    //!             3. `vec_real_type`
    //!
    //!             all the vector must have the same size.
    //!
    //! \param stream     stream to write the output
    //! \param delimiter  desired delimiter (optional). Default is tab.
    //!
    GenericContainer const &
    writeFormattedData( ostream_type & stream, char const delimiter = '\t' ) const;

    //!
    //! Read regular formatted data from `stream` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map
    //! which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `vector_type` which contais the vectors which are the
    //!              columns of the data readed of type `vec_real_type`.
    //!
    //! \param stream       stream to write the output
    //! \param commentChars lines beginnig with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   caracters used as delimiter for headers
    //!
    GenericContainer &
    readFormattedData(
      std::istream & stream,
      char const * commentChars = "#%",
      char const * delimiters = " \t"
    );

    //!
    //! Read regular formatted data from file `fname` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `vector_type` which contais the vectors which are the
    //!              columns of the data readed of type `vec_real_type`.
    //!
    //! \param fname        file name to be read
    //! \param commentChars lines beginnig with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   caracters used as delimiter for headers
    //!
    GenericContainer &
    readFormattedData(
      char const * fname,
      char const * commentChars = "#%",
      char const * delimiters   = " \t"
    );

    //!
    //! Do an exception
    //!
    //! \param[in] msg message of the exception
    //!
    static
    void
    exception( char const * msg ) GC_NO_RETURN;

  };

  // -------------------------------------------------------
  // support functions
  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      vector of `GenericContainer` with the columns of the table
  //! \param[in] stream    output stream
  //! \param[in] delimiter delimiter character between columns
  //!
  void
  writeTable(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  );

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //! \param[in] delimiter delimiter character between columns
  //!
  void
  writeTable(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  );

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //!
  void
  writeTableFormatted(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream
  );

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //!
  void
  writeTableFormatted(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream
  );

  // -------------------------------------------------------
  class GenericContainerExplorer {
  private:

    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    enum {
      GENERIC_CONTAINER_OK        = 0,
      GENERIC_CONTAINER_BAD_TYPE  = 1,
      GENERIC_CONTAINER_NO_DATA   = 2,
      GENERIC_CONTAINER_NOT_EMPTY = 3,
      GENERIC_CONTAINER_BAD_HEAD  = 4
    };
    #endif

    GenericContainer              data;
    std::deque<GenericContainer*> head;

    map_type         * ptr_map;
    map_type::iterator map_iterator;

  public:

    GenericContainerExplorer() { head.push_back(&data); }
    ~GenericContainerExplorer() {}

    GenericContainer *
    top() {
      GC_ASSERT(
        head.size() > 0,
        "GenericContainerExplorer::top() empty stack!"
      )
      GC_ASSERT(
        head.back() != nullptr,
        "GenericContainerExplorer::top() bad top pointer!"
      )
      return head.back();
    }

    GenericContainer const *
    top() const {
      GC_ASSERT(
        head.size() > 0,
        "GenericContainerExplorer::top() empty stack!"
      )
      GC_ASSERT(
        head.back() != nullptr,
        "GenericContainerExplorer::top() bad top pointer!"
      )
      return head.back();
    }

    void       * mem_ptr()       { return &data; }
    void const * mem_ptr() const { return &data; }

    int
    check( int data_type ) const {
      if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD;
      if ( data_type == head.back()->get_type() ) {
        if ( head.back() == nullptr ) return GENERIC_CONTAINER_NO_DATA;
        return GENERIC_CONTAINER_OK;
      } else {
        return GENERIC_CONTAINER_BAD_TYPE;
      }
    }

    int
    check_no_data( int data_type ) const {
      if ( head.empty() ) return GENERIC_CONTAINER_BAD_HEAD;
      if ( GC_NOTYPE == head.back()->get_type() ||
           data_type == head.back()->get_type() ) return GENERIC_CONTAINER_OK;
      return GENERIC_CONTAINER_NOT_EMPTY;
    }

    int
    pop() {
      if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA;
      head.pop_back();
      return GENERIC_CONTAINER_OK;
    }

    int
    push( GenericContainer * gc ) {
      head.push_back( gc );
      return GENERIC_CONTAINER_OK;
    }

    int
    push_vector_position( unsigned pos ) {
      int ok = check( GC_VECTOR );
      if ( ok == GENERIC_CONTAINER_OK  ) {
        GenericContainer * gc = &((*head.back())[pos]);
        head.push_back( gc );
      }
      return ok;
    }

    int
    push_map_position( char const * pos ) {
      int ok = check( GC_MAP );
      if ( ok == GENERIC_CONTAINER_OK  ) {
        GenericContainer * gc = &((*head.back())[pos]);
        head.push_back( gc );
      }
      return ok;
    }

    int
    init_map_key() {
      int ok = check( GC_MAP );
      if ( ok == GENERIC_CONTAINER_OK ) {
        ptr_map = &head.back()->get_map();
        map_iterator = ptr_map->begin();
      }
      return ok;
    }

    char const *
    next_map_key() {
      if ( map_iterator != ptr_map->end() )
        return map_iterator++->first.c_str();
      else
        return nullptr;
    }

    int
    reset() {
      if ( head.empty() ) return GENERIC_CONTAINER_NO_DATA;
      while ( head.size() > 1 ) head.pop_back();
      return GENERIC_CONTAINER_OK;
    }

  };
}

// do not define alias GC if use X11
#ifndef XlibSpecificationRelease
namespace GC = GC_namespace;
#endif

// for backward compatibility
namespace GenericContainerNamespace = GC_namespace;

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

//
// eof: GenericContainer.hh
//
