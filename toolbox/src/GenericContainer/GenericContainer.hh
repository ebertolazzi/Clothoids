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

#pragma once

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
#include <string_view>
#include <complex>
#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "GenericContainerConfig.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef GC_DO_ERROR
  #define GC_DO_ERROR(MSG) {                          \
    ostringstream ost;                                \
    ost << "in GenericContainer: " << MSG << '\n';    \
    GenericContainer::exception( ost.str().data() ); \
  }
#endif

#ifndef GC_ASSERT
  #define GC_ASSERT(COND,MSG) if ( !(COND) ) GC_DO_ERROR(MSG)
#endif

#ifndef GC_WARNING
  #define GC_WARNING(COND,MSG)                                   \
    if ( !(COND) ) {                                             \
      cout << "On line: " << __LINE__                            \
           << " file: " << __FILE__                              \
           << " in GenericContainer\nWARNING: " << MSG << '\n';  \
    }
#endif

#ifdef __GNUC__
  #define GC_NO_RETURN __attribute__ ((noreturn))
#else
  #define GC_NO_RETURN
#endif

#endif

//!
//! Namespace for the Generic Container
//!
namespace GC_namespace {

  using std::cin;
  using std::cout;
  using std::endl;
  using std::string;
  using std::string_view;
  using std::complex;
  using std::vector;
  using std::map;
  using std::deque;
  using std::istringstream;
  using std::ostringstream;
  using std::runtime_error;
  using std::exception;

  extern unsigned stream_number_precision;

  //!
  //! \brief Alias for a character-based output stream.
  //!
  //! This alias represents a `basic_ostream` specialized for `char` types,
  //! which is typically used for standard output operations like `cout` or
  //! file output streams.
  //!
  using ostream_type = std::basic_ostream<char>;

  //!
  //! \brief Alias for a character-based input stream.
  //!
  //! This alias represents a `basic_istream` specialized for `char` types,
  //! which is typically used for standard input operations like `cin` or
  //! file input streams.
  //!
  using istream_type = std::basic_istream<char>;

  #ifndef DOXYGEN_SHOULD_SKIP_THIS

  #if defined(GENERIC_CONTAINER_ON_WINDOWS) && defined(GENERIC_CONTAINER_USE_WINDOWS_TYPES)
  #else
  using std::uint8_t;
  using std::int32_t;
  using std::int64_t;
  using std::uint32_t;
  using std::uint64_t;
  #endif

  class GenericContainer;

  typedef void* pointer_type;

  using string           = string;
  using bool_type        = bool;
  using int_type         = int32_t;
  using long_type        = int64_t;
  using uint_type        = uint32_t;
  using ulong_type       = uint64_t;
  using real_type        = double;
  using complex_type     = complex<real_type>;
  using string_type      = string;
  using vec_pointer_type = vector<pointer_type>;
  using vec_bool_type    = vector<bool_type>;
  using vec_int_type     = vector<int_type>;
  using vec_long_type    = vector<long_type>;
  using vec_real_type    = vector<real_type>;
  using vec_complex_type = vector<complex_type>;
  using vec_string_type  = vector<string_type>;
  using vector_type      = vector<GenericContainer>;
  using map_type         = map<string_type,GenericContainer>;
  using vec_uint_type    = vector<uint_type>;
  using vec_ulong_type   = vector<ulong_type>;

  #endif

  // ---------------------------------------------------------------------------
  //!
  //! \brief Generic matrix storage type.
  //!
  //! This template class defines a matrix type that extends vector<TYPE>
  //! to store and manipulate a 2D matrix of elements of type `TYPE`.
  //! The matrix is stored internally as a 1D vector in row-major order.
  //!
  //! @tparam TYPE The type of elements stored in the matrix.
  //!
  template <typename TYPE>
  class mat_type : public vector<TYPE> {
    unsigned m_num_rows{0};  //!< Number of rows in the matrix.
    unsigned m_num_cols{0};  //!< Number of columns in the matrix.
    typedef typename vector<TYPE>::size_type size_type;
  public:

    //! Default constructor that creates an empty matrix.
    mat_type() = default;

    //!
    //! Constructs a matrix with given number of rows and columns.
    //!
    //! \param nr Number of rows.
    //! \param nc Number of columns.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 4);  // Creates a 3x4 matrix
    //! \endcode
    //!
    mat_type( unsigned nr, unsigned nc )
    : m_num_rows(nr)
    , m_num_cols(nc)
    { vector<TYPE>::resize(size_type(nr*nc)); }

    //!
    //! Resizes the matrix to the specified dimensions.
    //!
    //! \param nr New number of rows.
    //! \param nc New number of columns.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(2, 3);
    //! matrix.resize(4, 5);  // Resizes the matrix to 4x5
    //! \endcode
    //!
    void
    resize( unsigned nr, unsigned nc ) {
      m_num_rows = nr;
      m_num_cols = nc;
      vector<TYPE>::resize(size_type(nr*nc));
    }

    //!
    //! Copies the specified column of the matrix to a vector.
    //!
    //! \param nc The index of the column to be copied (0-based).
    //! \param C The vector that will be filled with the column elements.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! vector<int> column;
    //! matrix.get_column(1, column);  // Copies the second column into `column`
    //! \endcode
    //!
    void get_column( unsigned nc, vector<TYPE> & C ) const;

    //!
    //! \deprecated
    //!
    void
    getColumn( unsigned nc, vector<TYPE> & C ) const
    { this->get_column( nc, C ); }

    //!
    //! Copies the specified row of the matrix to a vector.
    //!
    //! \param nr The index of the row to be copied (0-based).
    //! \param R The vector that will be filled with the row elements.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! vector<int> row;
    //! matrix.get_row(0, row);  // Copies the first row into `row`
    //! \endcode
    //!
    void get_row( unsigned nr, vector<TYPE> & R ) const;

    //!
    //! \deprecated
    //!
    void
    getRow( unsigned nr, vector<TYPE> & R ) const
    { this->get_row( nr, R ); }

    //!
    //! Copies the specified column of the matrix to the given memory.
    //!
    //! \param nc The index of the column to be copied (0-based).
    //! \param C Pointer to memory where the column elements will be stored.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! int column[3];
    //! matrix.get_column(1, column);  // Copies the second column into `column`
    //! \endcode
    //!
    void get_column( unsigned nc, TYPE * C ) const;

    //!
    //! \deprecated
    //!
    void
    getColumn( unsigned nc, TYPE * C ) const
    { get_column( nc, C ); }

    //!
    //! Copies the specified row of the matrix to the given memory.
    //!
    //! \param nr The index of the row to be copied (0-based).
    //! \param R Pointer to memory where the row elements will be stored.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! int row[3];
    //! matrix.get_row(0, row);  // Copies the first row into `row`
    //! \endcode
    //!
    void get_row( unsigned nr, TYPE * R ) const;

    //!
    //! \deprecated
    //!
    void
    getRow( unsigned nr, TYPE * R ) const
    { this->get_row( nr, R ); }

    //!
    //! Returns the number of rows in the matrix.
    //!
    //! \return The number of rows.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! unsigned rows = matrix.num_rows();  // Returns 3
    //! \endcode
    //!
    [[nodiscard]] unsigned num_rows() const { return m_num_rows; }

    //!
    //! \deprecated
    //!
    [[nodiscard]] unsigned numRows() const { return m_num_rows; }

    //!
    //! Returns the number of columns in the matrix.
    //!
    //! \return The number of columns.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! unsigned cols = matrix.num_cols();  // Returns 3
    //! \endcode
    //!
    [[nodiscard]] unsigned num_cols() const { return m_num_cols; }

    //!
    //! \deprecated
    //!
    [[nodiscard]] unsigned numCols() const { return m_num_cols; }

    //!
    //! Provides constant access to the element at position (i, j).
    //!
    //! \param i Row index (0-based).
    //! \param j Column index (0-based).
    //! \return A constant reference to the element at (i, j).
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! const int &value = matrix(1, 1);  // Accesses element at (1,1)
    //! \endcode
    //!
    TYPE const & operator () ( unsigned i, unsigned j ) const;

    //!
    //! Provides mutable access to the element at position (i, j).
    //!
    //! \param i Row index (0-based).
    //! \param j Column index (0-based).
    //! \return A reference to the element at (i, j).
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! matrix(1, 1) = 42;  // Sets the element at (1,1) to 42
    //! \endcode
    //!
    TYPE & operator () ( unsigned i, unsigned j );

    //!
    //! Prints matrix information (dimensions and content) to the given output stream.
    //!
    //! \param stream The output stream where matrix information will be printed.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! matrix.info(cout);  // Prints matrix info to standard output
    //! \endcode
    //!
    void info( ostream_type & stream ) const;

    //!
    //! Returns a string containing matrix information (dimensions and content).
    //!
    //! \return A string with matrix information.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! string info = matrix.info();  // Returns matrix info as string
    //! \endcode
    //!
    string_type
    info() const {
      ostringstream ostr;
      this->info(ostr);
      return ostr.str();
    }

    //!
    //! Returns a pointer to the underlying data block of the matrix.
    //!
    //! \return A pointer to the first element of the matrix data array.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! int *dataPtr = matrix.data();  // Pointer to matrix data
    //! \endcode
    //!
    TYPE * data() { return &vector<TYPE>::front(); }

    //!
    //! Returns a constant pointer to the underlying data block of the matrix.
    //!
    //! \return A constant pointer to the first element of the matrix data array.
    //!
    //! ### Example
    //! \code
    //! mat_type<int> matrix(3, 3);
    //! const int *dataPtr = matrix.data();  // Constant pointer to matrix data
    //! \endcode
    //!
    TYPE const * data() const { return &vector<TYPE>::front(); }
  };

  // ---------------------------------------------------------------------------
  #ifndef GENERIC_CONTAINER_ON_WINDOWS
  extern template class mat_type<int_type>;
  extern template class mat_type<long_type>;
  extern template class mat_type<real_type>;
  extern template class mat_type<complex_type>;
  #endif

  using mat_int_type     = mat_type<int_type>;
  using mat_long_type    = mat_type<long_type>;
  using mat_real_type    = mat_type<real_type>;
  using mat_complex_type = mat_type<complex_type>;

  string to_string( complex_type const & v );

  //!
  //! \brief Overload of the `operator<<` for printing a vector of elements of a generic type.
  //!
  //! This function allows the printing of a `vector` containing elements of type `TYPE`.
  //! The elements are typically printed in a comma-separated format inside square brackets.
  //!
  //! \tparam TYPE The type of the elements contained in the vector.
  //! \param s The output stream to write to (e.g., `cout` or a `ostringstream`).
  //! \param v The `vector` object to print.
  //! \return The modified output stream after writing the vector elements.
  //!
  template <typename TYPE>
  ostream_type & operator << ( ostream_type & s, vector<TYPE> const & v );

  //!
  //! \brief Overload of the `operator<<` for printing a matrix of elements of a generic type.
  //!
  //! This function allows the printing of a `mat_type` object, which represents a matrix
  //! containing elements of type `TYPE`. The matrix is typically printed in a row-by-row format.
  //!
  //! \tparam TYPE The type of the elements contained in the matrix.
  //! \param  s    The output stream to write to (e.g., `cout` or a `ostringstream`).
  //! \param  mat  The `mat_type` object to print, representing the matrix.
  //! \return The modified output stream after writing the matrix elements.
  //!
  template <typename TYPE>
  ostream_type & operator << ( ostream_type & s, mat_type<TYPE> const & mat );

  // ---------------------------------------------------------------------------

  //!
  //! \brief Enum class representing types allowed for the `GenericContainer`.
  //!
  //! This enum class defines the types that are allowed to be used in the `GenericContainer`.
  //! The types are categorized as simple types, vector types, matrix types, and complex types.
  //!
  //! ### Example
  //! \code
  //! TypeAllowed type = TypeAllowed::INTEGER;  // Defines an integer type for GenericContainer
  //! \endcode
  //!
  using TypeAllowed = enum class GC_type : int_type {
    NOTYPE,      //!< No type assigned
    POINTER,     //!< Pointer type
    BOOL,        //!< Boolean type
    INTEGER,     //!< Integer type
    LONG,        //!< Long integer type
    REAL,        //!< Real number (floating-point) type
    COMPLEX,     //!< Complex number type
    STRING,      //!< String type

    // Vector types
    VEC_POINTER, //!< Vector of pointers
    VEC_BOOL,    //!< Vector of booleans
    VEC_INTEGER, //!< Vector of integers
    VEC_LONG,    //!< Vector of long integers
    VEC_REAL,    //!< Vector of real numbers
    VEC_COMPLEX, //!< Vector of complex numbers
    VEC_STRING,  //!< Vector of strings

    // Matrix types
    MAT_INTEGER, //!< Matrix of integers
    MAT_LONG,    //!< Matrix of long integers
    MAT_REAL,    //!< Matrix of real numbers
    MAT_COMPLEX, //!< Matrix of complex numbers

    // Complex types
    VECTOR,      //!< Vector container
    MAP          //!< Map (key-value container)
  };

  //!
  //! \brief Converts the `GC_type` enum value to a string representation.
  //!
  //! This function takes a `GC_type` enum value and returns a corresponding
  //! string representation for easier debugging and logging.
  //!
  //! \param s The `GC_type` enum value to convert.
  //! \return A constant C-string representing the type.
  //!
  //! ### Example
  //! \code
  //! GC_type type = GC_type::INTEGER;
  //! const char* typeStr = to_string(type);  // Returns "INTEGER"
  //! cout << "Type: " << typeStr << endl;  // Output: Type: INTEGER
  //! \endcode
  //!
  string_view to_string(GC_type s);

  //!
  //! \brief The `GenericContainer` class provides a flexible container for storing heterogeneous data types.
  //!
  //! This class allows storage of various data types, including primitive types (e.g., integer, floating-point),
  //! complex types, and containers like vectors and maps. It supports recursive structures where elements of the
  //! container can themselves be other `GenericContainer` objects.
  //!
  //! ### Supported Data Types:
  //! - **Primitive types:**
  //!   - Pointer
  //!   - Boolean
  //!   - Integer
  //!   - Long integer
  //!   - Floating point (real numbers)
  //!   - Complex floating point
  //!   - String
  //!
  //! - **Vector types:**
  //!   - Vector of pointers
  //!   - Vector of booleans
  //!   - Vector of integers
  //!   - Vector of long integers
  //!   - Vector of floating-point numbers
  //!   - Vector of complex numbers
  //!   - Vector of strings
  //!
  //! - **Matrix types:**
  //!   - Matrix of integers
  //!   - Matrix of long integers
  //!   - Matrix of floating-point numbers
  //!   - Matrix of complex numbers
  //!
  //! - **Complex container types:**
  //!   - Vector of `GenericContainer` (recursive container)
  //!   - Map of `GenericContainer` (key-value structure, recursive container)
  //!
  //! These capabilities make `GenericContainer` a highly versatile container for managing mixed-type data in C++.
  //!
  //! ### Example:
  //! \code
  //! GenericContainer gc;
  //! gc.set_integer(42);  // Store an integer
  //! gc.set_string("Hello, World!");  // Store a string
  //!
  //! GenericContainer vec_gc;
  //! vec_gc.set_vector(3);  // Store a vector of GenericContainers
  //! vec_gc[0].set_real(3.14);  // Set first element as a floating-point number
  //! \endcode
  //!
  class GenericContainer {
  public:
    // Import type aliases from namespace `GC_namespace`
    using pointer_type     = GC_namespace::pointer_type;       //!< Alias for pointer type
    using bool_type        = GC_namespace::bool_type;          //!< Alias for boolean type
    using int_type         = GC_namespace::int_type;           //!< Alias for integer type
    using uint_type        = GC_namespace::uint_type;          //!< Alias for unsigned integer type
    using long_type        = GC_namespace::long_type;          //!< Alias for long integer type
    using ulong_type       = GC_namespace::ulong_type;         //!< Alias for unsigned long integer type
    using real_type        = GC_namespace::real_type;          //!< Alias for real (floating point) type
    using complex_type     = GC_namespace::complex_type;       //!< Alias for complex number type
    using string_type      = GC_namespace::string_type;        //!< Alias for string type
    using vec_pointer_type = GC_namespace::vec_pointer_type;   //!< Alias for vector of pointers type
    using vec_bool_type    = GC_namespace::vec_bool_type;      //!< Alias for vector of booleans type
    using vec_int_type     = GC_namespace::vec_int_type;       //!< Alias for vector of integers type
    using vec_uint_type    = GC_namespace::vec_uint_type;      //!< Alias for vector of unsigned integers type
    using vec_long_type    = GC_namespace::vec_long_type;      //!< Alias for vector of long integers type
    using vec_ulong_type   = GC_namespace::vec_ulong_type;     //!< Alias for vector of unsigned long integers type
    using vec_real_type    = GC_namespace::vec_real_type;      //!< Alias for vector of real numbers type
    using vec_complex_type = GC_namespace::vec_complex_type;   //!< Alias for vector of complex numbers type
    using vec_string_type  = GC_namespace::vec_string_type;    //!< Alias for vector of strings type
    using vector_type      = GC_namespace::vector_type;        //!< Alias for vector of `GenericContainer` type
    using map_type         = GC_namespace::map_type;           //!< Alias for map of `GenericContainer` type
    using mat_int_type     = GC_namespace::mat_int_type;       //!< Alias for matrix of integers type
    using mat_long_type    = GC_namespace::mat_long_type;      //!< Alias for matrix of long integers type
    using mat_real_type    = GC_namespace::mat_real_type;      //!< Alias for matrix of real numbers type
    using mat_complex_type = GC_namespace::mat_complex_type;   //!< Alias for matrix of complex numbers type

  private:

    //!
    //! \brief Union for internal data storage.
    //!
    //! This union holds the actual data for the container in its various possible forms.
    //! Depending on the current data type stored in the container, different members of the union will be used.
    //!
    using DataStorage = union {
      pointer_type     p;    ///< Pointer data
      bool_type        b;    ///< Boolean data
      int_type         i;    ///< Integer data
      long_type        l;    ///< Long integer data
      real_type        r;    ///< Floating point (real) data
      complex_type     * c;  ///< Pointer to complex number data
      string_type      * s;  ///< Pointer to string data

      vec_pointer_type * v_p;   ///< Pointer to vector of pointers
      vec_bool_type    * v_b;   ///< Pointer to vector of booleans
      vec_int_type     * v_i;   ///< Pointer to vector of integers
      vec_long_type    * v_l;   ///< Pointer to vector of long integers
      vec_real_type    * v_r;   ///< Pointer to vector of real numbers
      vec_complex_type * v_c;   ///< Pointer to vector of complex numbers
      vec_string_type  * v_s;   ///< Pointer to vector of strings

      mat_int_type     * m_i;   ///< Pointer to matrix of integers
      mat_long_type    * m_l;   ///< Pointer to matrix of long integers
      mat_real_type    * m_r;   ///< Pointer to matrix of real numbers
      mat_complex_type * m_c;   ///< Pointer to matrix of complex numbers

      vector_type      * v;     ///< Pointer to vector of `GenericContainer`
      map_type         * m;     ///< Pointer to map of `GenericContainer`
    };

    DataStorage m_data;                       //!< The actual data stored in the container.
    TypeAllowed m_data_type{GC_type::NOTYPE}; //!< The type of data currently stored.

    //! \brief Allocates memory for a string.
    void allocate_string();

    //! \brief Allocates memory for a complex number.
    void allocate_complex();

    //! \brief Allocates memory for a vector of pointers of size `sz`.
    void allocate_vec_pointer(unsigned sz);

    //! \brief Allocates memory for a vector of booleans of size `sz`.
    void allocate_vec_bool(unsigned sz);

    //! \brief Allocates memory for a vector of integers of size `sz`.
    void allocate_vec_int(unsigned sz);

    //! \brief Allocates memory for a vector of long integers of size `sz`.
    void allocate_vec_long(unsigned sz);

    //! \brief Allocates memory for a vector of real numbers of size `sz`.
    void allocate_vec_real(unsigned sz);

    //! \brief Allocates memory for a vector of complex numbers of size `sz`.
    void allocate_vec_complex(unsigned sz);

    //! \brief Allocates memory for a matrix of integers of size `nr` x `nc`.
    void allocate_mat_int(unsigned nr, unsigned nc);

    //! \brief Allocates memory for a matrix of long integers of size `nr` x `nc`.
    void allocate_mat_long(unsigned nr, unsigned nc);

    //! \brief Allocates memory for a matrix of real numbers of size `nr` x `nc`.
    void allocate_mat_real(unsigned nr, unsigned nc);

    //! \brief Allocates memory for a matrix of complex numbers of size `nr` x `nc`.
    void allocate_mat_complex(unsigned nr, unsigned nc);

    //! \brief Allocates memory for a vector of strings of size `sz`.
    void allocate_vec_string(unsigned sz);

    //! \brief Allocates memory for a vector of `GenericContainer` of size `sz`.
    void allocate_vector(unsigned sz);

    //! \brief Allocates memory for a map of `GenericContainer`.
    void allocate_map();

    //! \brief Checks the type of data stored, throws error if type mismatch.
    void ck( string_view, TypeAllowed) const;

    //! \brief Checks the type of data stored, returns an error code for type mismatch.
    [[nodiscard]] int ck(TypeAllowed) const;

    //! \brief Checks or sets the type of data stored.
    void ck_or_set(string_view, TypeAllowed);

    //! \brief Returns true if the data type is a simple type (e.g., primitive).
    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    bool simple_data() const;
    #else
    [[nodiscard]] bool simple_data() const { return m_data_type <= GC_type::STRING; }
    #endif

    //! \brief Returns true if the data type is a simple vector type.
    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    bool simple_vec_data() const;
    #else
    [[nodiscard]] bool simple_vec_data() const { return m_data_type < GC_type::VEC_STRING; }
    #endif

  public:

    //!
    //! \brief Constructs a `GenericContainer` with an initial empty state.
    //!
    //! This constructor initializes a `GenericContainer` object with no data. The internal
    //! data type is set to `NOTYPE`, indicating that no data is currently stored in the container.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! // gc is now an empty container with no type assigned
    //! \endcode
    //!
    GenericContainer() : m_data_type(GC_type::NOTYPE) {}

    //!
    //! \brief Destroys the `GenericContainer` and releases any allocated resources.
    //!
    //! This destructor ensures that any dynamically allocated memory or resources associated with
    //! the data stored in the `GenericContainer` are properly freed. It automatically calls the
    //! `clear()` method to ensure the object is fully cleaned up before destruction.
    //!
    //! ### Example:
    //! \code
    //! {
    //!   GenericContainer gc;
    //!   // Perform operations with gc...
    //! }
    //! // When gc goes out of scope, the destructor is called automatically, freeing resources
    //! \endcode
    //!
    ~GenericContainer() { clear(); }

    //!
    //! \brief Clears the content of the `GenericContainer`, resetting it to an empty state.
    //!
    //! This method frees any memory or resources associated with the data currently stored in
    //! the container. After calling `clear()`, the container's data type is set back to `NOTYPE`.
    //!
    //! This function is particularly useful when the container is no longer needed or when you
    //! want to reuse the same `GenericContainer` object for different data.
    //!
    //! ### Example:
    //! \code
    //! gc.set_integer(42);  // Store an integer
    //! gc.clear();          // Clear the container, now it is empty
    //! \endcode
    //!
    void clear();

    //!
    //! \brief Removes an item from the map stored in the `GenericContainer` by its key.
    //!
    //! This method deletes an entry in the map contained within the `GenericContainer`,
    //! identified by the provided key `name`. The `GenericContainer` must currently store
    //! a map (`MAP` data type), otherwise, an error is thrown.
    //!
    //! \param name The key of the element to be removed from the map.
    //!
    //! \throws runtime_error If the current data type is not `MAP`.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_map();  // Initialize as a map
    //! gc["key1"].set_integer(10);  // Add an entry with key "key1"
    //! gc.erase("key1");  // Remove the entry with key "key1"
    //! \endcode
    //!
    void erase( string_view name ) const;

    //!
    //! \name Methods for Initializing Simple Data Types
    //!
    //! This group of methods allows setting simple data types (e.g., pointer, boolean, integer, etc.) within a `GenericContainer`.
    //! Each method initializes the container to the specified type, sets the value, and returns a reference to the stored data.
    //! @{

    //!
    //! \brief Set the data type to `pointer_type` and assign a value.
    //!
    //! This method initializes the `GenericContainer` to hold a pointer, stores the provided `value`,
    //! and returns a reference to the internal pointer storage.
    //!
    //! \param value The pointer value to store in the container.
    //! \return A reference to the stored pointer value.
    //!
    //! ### Example:
    //! \code
    //! int x = 42;
    //! GenericContainer gc;
    //! gc.set_pointer(&x);  // Store the pointer to x
    //! \endcode
    //!
    pointer_type & set_pointer( pointer_type value );

    //!
    //! \brief Free the pointer and reset the container to an empty state.
    //!
    //! This method deallocates any memory or resources associated with the current pointer,
    //! resets the container's type to `NOTYPE`, and returns a reference to the `GenericContainer` instance itself.
    //!
    //! \return A reference to the `GenericContainer` after the pointer has been freed.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_pointer(malloc(100));  // Allocate memory and store the pointer
    //! gc.free_pointer();            // Free the allocated memory and reset container
    //! \endcode
    //!
    GenericContainer & free_pointer();

    //!
    //! \brief Extracts the keys of the map stored in the `GenericContainer` into a vector of strings.
    //!
    //! If the container holds a map (`MAP` type), this method populates the provided `keys` vector
    //! with the keys from the map. If the container is not a map, the method throws an exception.
    //!
    //! \param[out] keys A vector of strings that will be filled with the map's keys.
    //! \throws runtime_error if the container does not hold a map.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_map();
    //! gc["key1"].set_int(1);
    //! gc["key2"].set_real(3.14);
    //! vector<string> keys;
    //! gc.get_keys(keys);  // keys will contain "key1" and "key2"
    //! \endcode
    //!
    void get_keys( vec_string_type & keys ) const;

    //!
    //! \brief Extracts the keys of the map stored in the `GenericContainer` as a comma-separated string.
    //!
    //! If the container holds a map (`MAP` type), this method returns the map's keys as a single string,
    //! where each key is separated by a comma. If the container is not a map, an exception is thrown.
    //!
    //! \return A string containing the keys of the map, separated by ", ".
    //! \throws runtime_error if the container does not hold a map.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_map();
    //! gc["key1"].set_int(1);
    //! gc["key2"].set_real(3.14);
    //! string keys = gc.get_keys();  // "key1, key2"
    //! \endcode
    //!
    string get_keys() const;

    //!
    //! \brief Set the data type to `bool_type` and assign a boolean value.
    //!
    //! This method initializes the `GenericContainer` to hold a boolean value, stores the provided value,
    //! and returns a reference to the stored boolean.
    //!
    //! \param value The boolean value to store.
    //! \return A reference to the stored boolean value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_bool(true);  // Store a boolean value
    //! \endcode
    //!
    bool_type & set_bool( bool_type value );

    //!
    //! \brief Set the data type to `int_type` and assign an integer value.
    //!
    //! This method initializes the `GenericContainer` to hold an integer value, stores the provided value,
    //! and returns a reference to the stored integer.
    //!
    //! \param value The integer value to store.
    //! \return A reference to the stored integer value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_int(42);  // Store an integer value
    //! \endcode
    //!
    int_type & set_int( int_type value );

    //!
    //! \brief Set the data type to `long_type` and assign a long integer value.
    //!
    //! This method initializes the `GenericContainer` to hold a long integer value, stores the provided value,
    //! and returns a reference to the stored long integer.
    //!
    //! \param value The long integer value to store.
    //! \return A reference to the stored long integer value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_long(123456789L);  // Store a long integer value
    //! \endcode
    //!
    long_type & set_long( long_type value );

    //!
    //! \brief Set the data type to `real_type` and assign a floating-point value.
    //!
    //! This method initializes the `GenericContainer` to hold a floating-point value, stores the provided value,
    //! and returns a reference to the stored real value.
    //!
    //! \param value The floating-point value to store.
    //! \return A reference to the stored floating-point value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_real(3.14);  // Store a floating-point value
    //! \endcode
    //!
    real_type & set_real( real_type value );

    //!
    //! \brief Set the data type to `complex_type` and assign a complex value.
    //!
    //! This method initializes the `GenericContainer` to hold a complex number, stores the provided value,
    //! and returns a reference to the stored complex value.
    //!
    //! \param value The complex number to store.
    //! \return A reference to the stored complex value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! complex<double> c(1.0, 2.0);
    //! gc.set_complex(c);  // Store a complex number
    //! \endcode
    //!
    complex_type & set_complex( complex_type const & value );

    //!
    //! \brief Set the data type to `complex_type` and assign a complex value from real and imaginary parts.
    //!
    //! This method initializes the `GenericContainer` to hold a complex number. It takes two real numbers
    //! (representing the real and imaginary components), stores the complex number, and returns a reference
    //! to the stored complex value.
    //!
    //! \param r The real part of the complex number.
    //! \param i The imaginary part of the complex number.
    //! \return A reference to the stored complex value.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_complex(1.0, 2.0);  // Store a complex number (1.0 + 2.0i)
    //! \endcode
    //!
    complex_type & set_complex( real_type re, real_type im );

    //!
    //! \brief Set the data type to `string_type`, allocate memory, and assign a string value.
    //!
    //! This method initializes the `GenericContainer` to hold a string, allocates the necessary memory,
    //! stores the provided string value, and returns a reference to the stored string.
    //!
    //! \param value The string to store.
    //! \return A reference to the stored string.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_string("Hello, World!");  // Store a string
    //! \endcode
    //!
    string_type & set_string( string_view value );
    ///@}





    //!
    //! \name Methods for Initializing Vector and Matrix Data
    //!
    //! This section includes methods for setting, allocating, and initializing vector and matrix data types
    //! within a `GenericContainer`. Each method allows for the creation of vectors or matrices, either by specifying
    //! the size or by copying from an existing vector or matrix.
    //! @{

    //! \brief Set the data to `vec_pointer_type`, allocate and initialize.
    //!
    //! This method allocates memory for a vector of pointers. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of pointers to allocate (default is 0).
    //! \return A reference to the internally allocated vector of pointers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_pointer(10);  // Allocates a vector of 10 pointers
    //! \endcode
    //!
    vec_pointer_type & set_vec_pointer( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_pointer_type& set_vec_pointer(T sz)
    { return set_vec_pointer(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_pointer_type` by copying from another vector.
    //!
    //! This method initializes the vector of pointers by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of pointers used for initialization.
    //! \return A reference to the internally allocated vector of pointers.
    //!
    //! ### Example:
    //! \code
    //! vec_pointer_type v = {ptr1, ptr2, ptr3};
    //! GenericContainer gc;
    //! gc.set_vec_pointer(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_pointer_type & set_vec_pointer( vec_pointer_type const & v );

    //! \brief Set the data to `vec_bool_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of boolean values. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of booleans to allocate (default is 0).
    //! \return A reference to the internally allocated vector of booleans.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_bool(5);  // Allocates a vector of 5 booleans
    //! \endcode
    //!
    vec_bool_type & set_vec_bool( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_bool_type& set_vec_bool(T sz)
    { return set_vec_bool(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_bool_type` by copying from another vector.
    //!
    //! This method initializes the vector of booleans by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of booleans used for initialization.
    //! \return A reference to the internally allocated vector of booleans.
    //!
    //! ### Example:
    //! \code
    //! vec_bool_type v = {true, false, true};
    //! GenericContainer gc;
    //! gc.set_vec_bool(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_bool_type & set_vec_bool( vec_bool_type const & v );

    //! \brief Set the data to `vec_int_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of integers. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of integers to allocate (default is 0).
    //! \return A reference to the internally allocated vector of integers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_int(10);  // Allocates a vector of 10 integers
    //! \endcode
    //!
    vec_int_type & set_vec_int( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_int_type& set_vec_int(T sz)
    { return set_vec_int(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_int_type` by copying from another vector.
    //!
    //! This method initializes the vector of integers by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of integers used for initialization.
    //! \return A reference to the internally allocated vector of integers.
    //!
    //! ### Example:
    //! \code
    //! vec_int_type v{1, 2, 3};
    //! GenericContainer gc;
    //! gc.set_vec_int(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_int_type & set_vec_int( vec_int_type const & v );

    //! \brief Set the data to `vec_long_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of long integers. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of long integers to allocate (default is 0).
    //! \return A reference to the internally allocated vector of long integers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_long(10);  // Allocates a vector of 10 long integers
    //! \endcode
    //!
    vec_long_type & set_vec_long( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_long_type& set_vec_long(T sz)
    { return set_vec_long(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_long_type` by copying from another vector.
    //!
    //! This method initializes the vector of long integers by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of long integers used for initialization.
    //! \return A reference to the internally allocated vector of long integers.
    //!
    //! ### Example:
    //! \code
    //! vec_long_type v{100000L, 200000L};
    //! GenericContainer gc;
    //! gc.set_vec_long(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_long_type & set_vec_long( vec_long_type const & v );

    //! \brief Set the data to `vec_real_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of `real_type` numbers. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of real numbers to allocate (default is 0).
    //! \return A reference to the internally allocated vector of real numbers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_real(5);  // Allocates a vector of 5 real numbers
    //! \endcode
    //!
    vec_real_type & set_vec_real( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_real_type& set_vec_real(T sz)
    { return set_vec_real(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_real_type` by copying from another vector.
    //!
    //! This method initializes the vector of `real_type` numbers by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of `real_type` numbers used for initialization.
    //! \return A reference to the internally allocated vector of `real_type` numbers.
    //!
    //! ### Example:
    //! \code
    //! vec_real_type v = {1.5, 2.5, 3.5};
    //! GenericContainer gc;
    //! gc.set_vec_real(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_real_type & set_vec_real( vec_real_type const & v );

    //! \brief Set the data to `vec_complex_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of complex numbers. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of complex numbers to allocate (default is 0).
    //! \return A reference to the internally allocated vector of complex numbers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_complex(10);  // Allocates a vector of 10 complex numbers
    //! \endcode
    //!
    vec_complex_type & set_vec_complex( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_complex_type& set_vec_complex(T sz)
    { return set_vec_complex(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_complex_type` by copying from another vector.
    //!
    //! This method initializes the vector of complex numbers by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of complex numbers used for initialization.
    //! \return A reference to the internally allocated vector of complex numbers.
    //!
    //! ### Example:
    //! \code
    //! vec_complex_type v{ {1.0, 2.0}, {3.0, 4.0} };
    //! GenericContainer gc;
    //! gc.set_vec_complex(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_complex_type & set_vec_complex( vec_complex_type const & v );

    //! \brief Set the data to `vec_string_type`, allocate and initialize.
    //!
    //! Allocates memory for a vector of strings. If the specified size `sz` is greater than zero,
    //! the vector will be allocated with that size.
    //!
    //! \param[in] sz The size of the vector of strings to allocate (default is 0).
    //! \return A reference to the internally allocated vector of strings.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_vec_string(5);  // Allocates a vector of 5 strings
    //! \endcode
    //!
    vec_string_type & set_vec_string( unsigned sz = 0 );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    vec_string_type& set_vec_string(T sz)
    { return set_vec_string(static_cast<unsigned>(sz)); }

    //! \brief Set the data to `vec_string_type` by copying from another vector.
    //!
    //! This method initializes the vector of strings by copying data from the provided vector `v`.
    //!
    //! \param[in] v The vector of strings used for initialization.
    //! \return A reference to the internally allocated vector of strings.
    //!
    //! ### Example:
    //! \code
    //! vec_string_type v{"Hello", "World"};
    //! GenericContainer gc;
    //! gc.set_vec_string(v);  // Copies data from vector v
    //! \endcode
    //!
    vec_string_type & set_vec_string( vec_string_type const & v );

    //! \brief Set the data to `mat_int_type`, allocate and initialize.
    //!
    //! Allocates memory for a matrix of integers. If the specified number of rows `nr` and columns `nc`
    //! are greater than zero, the matrix will be allocated with size `nr` x `nc`.
    //!
    //! \param[in] nr The number of rows in the matrix (default is 0).
    //! \param[in] nc The number of columns in the matrix (default is 0).
    //! \return A reference to the internally allocated matrix of integers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_mat_int(3, 4);  // Allocates a 3x4 matrix of integers
    //! \endcode
    //!
    mat_int_type & set_mat_int( unsigned nr = 0, unsigned nc = 0 );

    //! \brief Set the data to `mat_int_type` by copying from another matrix.
    //!
    //! This method initializes the matrix of integers by copying data from the provided matrix `m`.
    //!
    //! \param[in] m The matrix of integers used for initialization.
    //! \return A reference to the internally allocated matrix of integers.
    //!
    //! ### Example:
    //! \code
    //! mat_int_type m{{1, 2}, {3, 4}};
    //! GenericContainer gc;
    //! gc.set_mat_int(m);  // Copies data from matrix m
    //! \endcode
    //!
    mat_int_type & set_mat_int( mat_int_type const & m );

    //! \brief Set the data to `mat_long_type`, allocate and initialize.
    //!
    //! Allocates memory for a matrix of long integers. If the specified number of rows `nr` and columns `nc`
    //! are greater than zero, the matrix will be allocated with size `nr` x `nc`.
    //!
    //! \param[in] nr The number of rows in the matrix (default is 0).
    //! \param[in] nc The number of columns in the matrix (default is 0).
    //! \return A reference to the internally allocated matrix of long integers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_mat_long(3, 4);  // Allocates a 3x4 matrix of long integers
    //! \endcode
    //!
    mat_long_type & set_mat_long( unsigned nr = 0, unsigned nc = 0 );

    //! \brief Set the data to `mat_long_type` by copying from another matrix.
    //!
    //! This method initializes the matrix of long integers by copying data from the provided matrix `m`.
    //!
    //! \param[in] m The matrix of long integers used for initialization.
    //! \return A reference to the internally allocated matrix of long integers.
    //!
    //! ### Example:
    //! \code
    //! mat_long_type m = {{100000L, 200000L}, {300000L, 400000L}};
    //! GenericContainer gc;
    //! gc.set_mat_long(m);  // Copies data from matrix m
    //! \endcode
    //!
    mat_long_type & set_mat_long( mat_long_type const & m );

    //! \brief Set the data to `mat_real_type`, allocate and initialize.
    //!
    //! Allocates memory for a matrix of `real_type` numbers. If the specified number of rows `nr` and columns `nc`
    //! are greater than zero, the matrix will be allocated with size `nr` x `nc`.
    //!
    //! \param[in] nr The number of rows in the matrix (default is 0).
    //! \param[in] nc The number of columns in the matrix (default is 0).
    //! \return A reference to the internally allocated matrix of `real_type` numbers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_mat_real(3, 4);  // Allocates a 3x4 matrix of real numbers
    //! \endcode
    //!
    mat_real_type & set_mat_real( unsigned nr = 0, unsigned nc = 0 );

    //! \brief Set the data to `mat_real_type` by copying from another matrix.
    //!
    //! This method initializes the matrix of `real_type` numbers by copying data from the provided matrix `m`.
    //!
    //! \param[in] m The matrix of `real_type` numbers used for initialization.
    //! \return A reference to the internally allocated matrix of `real_type` numbers.
    //!
    //! ### Example:
    //! \code
    //! mat_real_type m = {{1.1, 2.2}, {3.3, 4.4}};
    //! GenericContainer gc;
    //! gc.set_mat_real(m);  // Copies data from matrix m
    //! \endcode
    //!
    mat_real_type & set_mat_real( mat_real_type const & m );

    //! \brief Set the data to `mat_complex_type`, allocate and initialize.
    //!
    //! Allocates memory for a matrix of `complex_type` numbers. If the specified number of rows `nr` and columns `nc`
    //! are greater than zero, the matrix will be allocated with size `nr` x `nc`.
    //!
    //! \param[in] nr The number of rows in the matrix (default is 0).
    //! \param[in] nc The number of columns in the matrix (default is 0).
    //! \return A reference to the internally allocated matrix of complex numbers.
    //!
    //! ### Example:
    //! \code
    //! GenericContainer gc;
    //! gc.set_mat_complex(3, 4);  // Allocates a 3x4 matrix of complex numbers
    //! \endcode
    //!
    mat_complex_type & set_mat_complex( unsigned nr = 0, unsigned nc = 0 );

    //! \brief Set the data to `mat_complex_type` by copying from another matrix.
    //!
    //! This method initializes the matrix of complex numbers by copying data from the provided matrix `m`.
    //!
    //! \param[in] m The matrix of complex numbers used for initialization.
    //! \return A reference to the internally allocated matrix of complex numbers.
    //!
    //! ### Example:
    //! \code
    //! mat_complex_type m = {{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}};
    //! GenericContainer gc;
    //! gc.set_mat_complex(m);  // Copies data from matrix m
    //! \endcode
    //!
    mat_complex_type & set_mat_complex( mat_complex_type const & m );

    //! \brief Push a boolean value into the vector or matrix.
    //!
    //! This method adds a boolean value to `vec_bool_type` or a general `vector_type`.
    //!
    //! \param[in] b value The boolean value to push.
    //!
    void push_bool( bool b ) const;

    //! \brief Push an integer value into the vector or matrix.
    //!
    //! This method adds an integer value to `vec_int_type` or a general `vector_type`.
    //!
    //! \param[in] i value The integer value to push.
    //!
    void push_int( int_type i );

    //! \brief Push a long integer value into the vector or matrix.
    //!
    //! This method adds a long integer value to `vec_long_type` or a general `vector_type`.
    //!
    //! \param[in] l value The long integer value to push.
    //!
    void push_long( long_type l );

    //! \brief Push a real number into the vector or matrix.
    //!
    //! This method adds a real number to `vec_real_type` or a general `vector_type`.
    //!
    //! \param[in] r value The real number to push.
    //!
    void push_real( real_type r );

    //! \brief Push a complex number into the vector or matrix using an existing complex object.
    //!
    //! This method adds a complex number to `vec_complex_type` or a general `vector_type`.
    //!
    //! \param[in] c value The complex number to push.
    //!
    void push_complex( complex_type & c );

    //! \brief Push a complex number into the vector or matrix using real and imaginary parts.
    //!
    //! This method adds a complex number to `vec_complex_type` or a general `vector_type`.
    //!
    //! \param[in] re The real part of the complex number.
    //! \param[in] im The imaginary part of the complex number.
    //!
    void push_complex( real_type re, real_type im );

    //! \brief Push a string value into the vector or matrix.
    //!
    //! This method adds a string to `vec_string_type` or a general `vector_type`.
    //!
    //! \param[in] s value The string value to push.
    //!
    void push_string( string_view s );

    ///@}






    //!
    //! \name Initialize generic data.
    //!
    ///@{

    //! \brief Initializes a generic vector.
    //!
    //! This function sets the data type to `vector_type`, allocates an empty vector, and returns a reference to it.
    //! If the specified size `sz` is greater than 0, the vector is allocated with that size.
    //!
    //! \param[in] sz The size of the allocated vector (default is 0).
    //! \return A reference to the initialized vector.
    //!
    //! \code
    //! GenericContainer container;
    //! auto & vec = container.set_vector(5); // Initializes a vector of size 5
    //! \endcode
    vector_type & set_vector( unsigned sz = 0 );

    //! \brief Initializes a generic map.
    //!
    //! This function sets the data type to `map_type`, allocates an empty map, and returns a reference to it.
    //!
    //! \return A reference to the initialized map.
    //!
    //! \code
    //! GenericContainer container;
    //! auto & my_map = container.set_map(); // Initializes a map
    //! \endcode
    map_type & set_map();

    ///@}






    //!
    //! \name Access to a single element.
    //!
    ///@{

    //!
    //! \brief Return an integer representing the type of data stored.
    //!
    //! This function returns an integer corresponding to the type of data
    //! that is currently stored in the container. The integer is mapped
    //! to specific data types as follows:
    //!
    //! - 0: No data stored
    //! - 1: pointer_type
    //! - 2: bool_type
    //! - 3: int_type
    //! - 4: long_type
    //! - 5: real_type
    //! - 6: complex_type
    //! - 7: string_data
    //! - 8: vec_pointer_type
    //! - 9: vec_bool_type
    //! - 10: vec_int_type
    //! - 11: vec_long_type
    //! - 12: vec_real_type
    //! - 13: vec_complex_type
    //! - 14: vec_string_type
    //! - 15: mat_int_type
    //! - 16: mat_long_type
    //! - 17: mat_real_type
    //! - 18: mat_complex_type
    //! - 19: vector_type
    //! - 20: map_type
    //!
    //! \return The type of the internally stored data as an integer.
    //!
    TypeAllowed get_type() const { return m_data_type; }

    //!
    //! \brief Return a string representing the type of data stored.
    //!
    //! This function returns a pointer to a string that describes the
    //! type of data currently held by the container. This is helpful
    //! for debugging and logging purposes.
    //!
    //! \return A pointer to a string representation of the data type.
    //!
    string_view get_type_name() const { return to_string(get_type()); }

    //!
    //! \brief Print information about the kind of data stored to a stream.
    //!
    //! This method outputs information about the data type and its
    //! properties to the provided output stream.
    //!
    //! \param[out] stream The output stream to write the information to.
    //! \return A reference to the current GenericContainer object.
    //!
    GenericContainer const & info( ostream_type & stream ) const;

    //!
    //! \brief Print information about the kind of data stored as a string.
    //!
    //! This method returns a string containing information about the
    //! data type and its properties.
    //!
    //! \return A string representation of the data information.
    //!
    string_type
    info() const {
      ostringstream sstr;
      this->info(sstr);
      return sstr.str();
    }

    //!
    //! \brief Return the number of elements in the first level of the generic container.
    //!
    //! - For scalar elements (e.g., boolean, integer, double), it returns 1.
    //! - For vector types, it returns the size of the vector.
    //! - For maps, it returns the number of keys.
    //!
    //! \return The number of elements in the first level of the container.
    //!
    unsigned get_num_elements() const;

    //!
    //! \brief Return the number of rows in the internally stored matrix.
    //!
    //! \return The number of rows in the matrix.
    //!
    unsigned num_rows() const;
    //!
    //! \deprecated
    //!
    unsigned get_numRows() const { return this->num_rows(); }

    //!
    //! \brief Return the number of columns in the internally stored matrix.
    //!
    //! \return The number of columns in the matrix.
    //!
    unsigned num_cols() const;
    //!
    //! \deprecated
    //!
    unsigned get_numCols() const { return this->num_cols(); }

    //!
    //! \brief Get a stored numeric value if the data is boolean, integer, or real type.
    //!
    //! \param[in] where Optional parameter to provide context for error messages.
    //! \return The numeric value, or `0` if the data is of an unsupported type.
    //!
    real_type get_number( string_view const where = "" ) const;

    //!
    //! \brief Get a stored complex number if the data is boolean, integer, real, or complex type.
    //!
    //! \param[in] where Optional parameter to provide context for error messages.
    //! \return The complex number, or `0` if the data is of an unsupported type.
    //!
    complex_type get_complex_number( string_view const where = "" ) const;

    //!
    //! \brief Get the real and imaginary parts of a stored complex number.
    //!
    //! This function extracts the real and imaginary parts of the complex number
    //! and stores them in the provided references.
    //!
    //! \param[out] re Reference to store the real part of the complex number.
    //! \param[out] im Reference to store the imaginary part of the complex number.
    //!
    void get_complex_number( real_type & re, real_type & im ) const;

    //!
    //! \brief Return the stored data as a generic pointer.
    //!
    //! \param[in] where Optional parameter to provide context for error messages.
    //! \return A void pointer to the stored data.
    //!
    void * get_pvoid( string_view const where = "" ) const;

    //!
    //! \brief Return the stored data as a double pointer.
    //!
    //! \param[in] where Optional parameter to provide context for error messages.
    //! \return A double pointer to the stored data.
    //!
    void ** get_ppvoid( string_view const where = "" ) const;

    //!
    //! Return the stored data as a pointer to const integer
    //!
    int_type const * get_int_pointer() const;

    //!
    //! \brief Return the stored data as a pointer to const integer.
    //!
    //! \return A pointer to const integer data.
    //!
    int_type * get_int_pointer();

    //!
    //! \brief Return the stored data as a pointer to const long.
    //!
    //! \return A pointer to const long data.
    //!
    long_type const * get_long_pointer() const;

    //!
    //! \brief Return the stored data as a pointer to long.
    //!
    //! \return A pointer to long data.
    //!
    long_type * get_long_pointer();

    //!
    //! \brief Return the stored data as a pointer to const real_type.
    //!
    //! \return A pointer to const real_type data.
    //!
    real_type const * get_real_pointer() const;

    //!
    //! \brief Return the stored data as a pointer to real_type.
    //!
    //! \return A pointer to real_type data.
    //!
    real_type * get_real_pointer();

    //!
    //! \brief Return the stored data as a pointer to const complex_type.
    //!
    //! \return A pointer to const complex_type data.
    //!
    complex_type const * get_complex_pointer() const;

    //!
    //! \brief Return the stored data as a pointer to complex_type.
    //!
    //! \return A pointer to complex_type data.
    //!
    complex_type * get_complex_pointer();

    //!
    //! \brief Get the stored value.
    //!
    //! This template function retrieves the stored value and assigns it to the
    //! provided reference variable.
    //!
    //! \param[out] v The reference to store the copied value.
    //! \param[in] where Optional parameter to provide context for error messages.
    //!
    template <typename T>
    void
    get_value( T & v, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a pointer.
    //!
    //! This function retrieves the stored pointer value.
    //!
    #ifdef GENERIC_CONTAINER_ON_WINDOWS
    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_type::POINTER); return *reinterpret_cast<T*>(get_ppvoid()); }

    //!
    //! Get the stored value as a pointer
    //!
    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_type::POINTER); return reinterpret_cast<T>(get_pvoid()); }
    #else
    template <typename T>
    T& get_pointer()
    { ck("get_pointer",GC_type::POINTER); return *reinterpret_cast<T*>(&(m_data.p)); }

    //!
    //! \brief Get the stored value as a pointer (const version).
    //!
    //! This function retrieves the stored pointer value as a const.
    //!
    template <typename T>
    T get_pointer() const
    { ck("get_pointer",GC_type::POINTER); return reinterpret_cast<T>(m_data.p); }
    #endif

    //!
    //! \brief Get the stored value in the map as boolean.
    //!
    //! This function retrieves a boolean value from the map using the specified key.
    //!
    //! \param[in] key The key of the map to be selected.
    //! \param[in] where Optional parameter to provide context for error messages.
    //! \return The boolean value stored in the container.
    //!
    bool_type get_map_bool( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as boolean.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the boolean stored in the container
    //!
    bool_type get_map_bool( std::initializer_list<string> args ) const;

    //!
    //! Get the stored value in the map as boolean.
    //! The key is searched between the strings in `keys`
    //! the first one found return the boolena value.
    //!
    //! \param[in] keys  keys of the map to be selected
    //! \param[in] where position added to the error message
    //! \return the boolean stored in the container
    //!
    bool_type get_map_bool( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as an integer.
    //!
    //! This function retrieves an integer value from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The integer value stored in the container.
    //!
    int_type get_map_int( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as an integer.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the integer stored in the container
    //!
    int_type get_map_int( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as an integer from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding integer value.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The integer value stored in the container for the first found key.
    //!
    int_type get_map_int( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as a real number.
    //!
    //! This function retrieves a real number from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The real number stored in the container.
    //!
    real_type get_map_number( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as a real number.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the real number stored in the container
    //!
    real_type get_map_number( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as a real number from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding real number value.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The real number stored in the container for the first found key.
    //!
    real_type get_map_number( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as a string.
    //!
    //! This function retrieves a string from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the string stored in the container.
    //!
    string_view get_map_string( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as a string.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the string stored in the container
    //!
    string_view get_map_string( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as a string from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding string value.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the string stored in the container for the first found key.
    //!
    string_view get_map_string( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as a vector of real numbers.
    //!
    //! This function retrieves a vector of real numbers from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of real numbers stored in the container.
    //!
    vec_real_type const & get_map_vec_real( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as a  vector of real numbers.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the vector of real numbers stored in the container
    //!
    vec_real_type const & get_map_vec_real( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as a vector of real numbers from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding vector of real numbers.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of real numbers stored in the container for the first found key.
    //!
    vec_real_type const & get_map_vec_real( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as a vector of complex numbers.
    //!
    //! This function retrieves a vector of complex numbers from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of complex numbers stored in the container.
    //!
    vec_complex_type const & get_map_vec_complex( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as a  vector of complex numbers.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the vector of complex numbers stored in the container
    //!
    vec_complex_type const & get_map_vec_complex( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as a vector of complex numbers from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding vector of complex numbers.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of complex numbers stored in the container for the first found key.
    //!
    vec_complex_type const & get_map_vec_complex( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value in the map as a vector of strings.
    //!
    //! This function retrieves a vector of strings from the map using the specified key.
    //!
    //! \param[in] key   The key of the map to be selected.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of strings stored in the container.
    //!
    vec_string_type const & get_map_vec_string( string_view const key, string_view const where = "" ) const;

    //!
    //! Get the stored value in the map as a  vector of strings.
    //!
    //! \param[in] args  keys lists sequence
    //! \return the vector of strings stored in the container
    //!
    vec_string_type const & get_map_vec_string( std::initializer_list<string> args ) const;

    //!
    //! \brief Get the stored value in the map as a vector of strings from a list of keys.
    //!
    //! This function searches for the first matching key among the provided strings in `keys`
    //! and returns the corresponding vector of strings.
    //!
    //! \param[in] keys  The list of keys to search in the map.
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the vector of strings stored in the container for the first found key.
    //!
    vec_string_type const & get_map_vec_string( vec_string_type const & keys, string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a boolean.
    //!
    //! This function retrieves a boolean value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the boolean stored in the container.
    //!
    bool_type & get_bool( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const boolean.
    //!
    //! This function retrieves a const boolean value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const boolean stored in the container.
    //!
    bool_type const & get_bool( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as an integer.
    //!
    //! This function retrieves an integer value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the integer stored in the container.
    //!
    int_type & get_int( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const integer.
    //!
    //! This function retrieves a const integer value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const integer stored in the container.
    //!
    int_type const & get_int( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a long integer.
    //!
    //! This function retrieves a long integer value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the long integer stored in the container.
    //!
    long_type & get_long( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const long integer.
    //!
    //! This function retrieves a const long integer value from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const long integer stored in the container.
    //!
    long_type const & get_long( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as an integer.
    //!
    //! This function retrieves the data stored in the container as an integer.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The data stored in the container as an integer.
    //!
    int_type get_as_int( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as an unsigned integer.
    //!
    //! This function retrieves the data stored in the container as an unsigned integer.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The data stored in the container as an unsigned integer.
    //!
    uint_type get_as_uint( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a long integer.
    //!
    //! This function retrieves the data stored in the container as a long integer.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The data stored in the container as a long integer.
    //!
    long_type get_as_long( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as an unsigned long integer.
    //!
    //! This function retrieves the data stored in the container as an unsigned long integer.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return The data stored in the container as an unsigned long integer.
    //!
    ulong_type get_as_ulong( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a real number.
    //!
    //! This function retrieves a real number from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the real number stored in the container.
    //!
    real_type & get_real( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const real number.
    //!
    //! This function retrieves a const real number from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const real number stored in the container.
    //!
    real_type const & get_real( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a complex number.
    //!
    //! This function retrieves a complex number from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the complex number stored in the container.
    //!
    complex_type & get_complex( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const complex number.
    //!
    //! This function retrieves a const complex number from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const complex number stored in the container.
    //!
    complex_type const & get_complex( string_view const where = "" ) const;

    //!
    //! \brief Get the stored value as a string.
    //!
    //! This function retrieves a string from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the string stored in the container.
    //!
    string_type & get_string( string_view const where = "" );

    //!
    //! \brief Get the stored value as a const string.
    //!
    //! This function retrieves a const string from the container.
    //!
    //! \param[in] where Optional context for error messages, indicating the position of the call.
    //! \return A reference to the const string stored in the container.
    //!
    string_view get_string( string_view const where = "" ) const;

    ///@}








    //!
    //! \name Access to vector type data.
    //!
    ///@{


    //!
    //! Get the stored value as a vector of `GenericoContainer`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of `GenericoContainer`
    //!
    //! \code
    //! GenericoContainer container;
    //! vector_type &data = container.get_vector("In function get_vector");
    //! // Use data...
    //! \endcode
    //!
    vector_type & get_vector( string_view const where = "" );

    //!
    //! Get the stored value as a const vector of `GenericoContainer`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of `GenericoContainer`
    //!
    //! \code
    //! const GenericoContainer container;
    //! vector_type const &data = container.get_vector("In function get_vector");
    //! // Use data...
    //! \endcode
    //!
    vector_type const & get_vector( string_view const where = "" ) const;

    //!
    //! Get the stored value as a vector of pointers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of pointers
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_pointer_type &ptr_data = container.get_vec_pointer("In function get_vec_pointer");
    //! // Use ptr_data...
    //! \endcode
    //!
    vec_pointer_type & get_vec_pointer( string_view const where = "" );

    //!
    //! Get the stored value as a const vector of pointers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of pointers
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_pointer_type const &ptr_data = container.get_vec_pointer("In function get_vec_pointer");
    //! // Use ptr_data...
    //! \endcode
    //!
    vec_pointer_type const & get_vec_pointer( string_view const where = "" ) const;

    //!
    //! Get the stored value as a vector of booleans
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of booleans
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_bool_type &bool_data = container.get_vec_bool("In function get_vec_bool");
    //! // Use bool_data...
    //! \endcode
    //!
    vec_bool_type & get_vec_bool( string_view const where = "" );

    //!
    //! Get the stored value as a const vector of booleans
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of booleans
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_bool_type const &bool_data = container.get_vec_bool("In function get_vec_bool");
    //! // Use bool_data...
    //! \endcode
    //!
    vec_bool_type const & get_vec_bool( string_view const where = "" ) const;

    //!
    //! Get the stored value as a vector of integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of integers
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_int_type &int_data = container.get_vec_int("In function get_vec_int");
    //! // Use int_data...
    //! \endcode
    //!
    vec_int_type & get_vec_int( string_view = "" );

    //!
    //! Get the stored value as a const vector of integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of integers
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_int_type const &int_data = container.get_vec_int("In function get_vec_int");
    //! // Use int_data...
    //! \endcode
    //!
    vec_int_type const & get_vec_int( string_view = "" ) const;

    //!
    //! Get the stored value as a vector of long integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of long integers
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_long_type &long_data = container.get_vec_long("In function get_vec_long");
    //! // Use long_data...
    //! \endcode
    //!
    vec_long_type & get_vec_long( string_view = "" );

    //!
    //! Get the stored value as a const vector of long integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of long integers
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_long_type const &long_data = container.get_vec_long("In function get_vec_long");
    //! // Use long_data...
    //! \endcode
    //!
    vec_long_type const & get_vec_long( string_view = "" ) const;

    //!
    //! Get the stored value as a vector of `real_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of `real_type`
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_real_type &real_data = container.get_vec_real("In function get_vec_real");
    //! // Use real_data...
    //! \endcode
    //!
    vec_real_type & get_vec_real( string_view = "" );

    //!
    //! Get the stored value as a const vector of `real_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of `real_type`
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_real_type const &real_data = container.get_vec_real("In function get_vec_real");
    //! // Use real_data...
    //! \endcode
    //!
    vec_real_type const & get_vec_real( string_view = "" ) const;

    //!
    //! Get the stored value as a vector of `complex_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of `complex_type`
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_complex_type &complex_data = container.get_vec_complex("In function get_vec_complex");
    //! // Use complex_data...
    //! \endcode
    //!
    vec_complex_type & get_vec_complex( string_view = "" );

    //!
    //! Get the stored value as a const vector of `complex_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of `complex_type`
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_complex_type const &complex_data = container.get_vec_complex("In function get_vec_complex");
    //! // Use complex_data...
    //! \endcode
    //!
    vec_complex_type const & get_vec_complex( string_view = "" ) const;

    //!
    //! Get the stored value as a matrix of integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a matrix of integers
    //!
    //! \code
    //! GenericoContainer container;
    //! mat_int_type &int_matrix = container.get_mat_int("In function get_mat_int");
    //! // Use int_matrix...
    //! \endcode
    //!
    mat_int_type & get_mat_int( string_view = "" );

    //!
    //! Get the stored value as a const matrix of integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const matrix of integers
    //!
    //! \code
    //! const GenericoContainer container;
    //! mat_int_type const &int_matrix = container.get_mat_int("In function get_mat_int");
    //! // Use int_matrix...
    //! \endcode
    //!
    mat_int_type const & get_mat_int( string_view = "" ) const;

    //!
    //! Get the stored value as a matrix of long integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a matrix of long integers
    //!
    //! \code
    //! GenericoContainer container;
    //! mat_long_type &long_matrix = container.get_mat_long("In function get_mat_long");
    //! // Use long_matrix...
    //! \endcode
    //!
    mat_long_type & get_mat_long( string_view = "" );

    //!
    //! Get the stored value as a const matrix of long integers
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const matrix of long integers
    //!
    //! \code
    //! const GenericoContainer container;
    //! mat_long_type const &long_matrix = container.get_mat_long("In function get_mat_long");
    //! // Use long_matrix...
    //! \endcode
    //!
    mat_long_type const & get_mat_long( string_view = "" ) const;

    //!
    //! Get the stored value as a matrix of `real_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a matrix of `real_type`
    //!
    //! \code
    //! GenericoContainer container;
    //! mat_real_type &real_matrix = container.get_mat_real("In function get_mat_real");
    //! // Use real_matrix...
    //! \endcode
    //!
    mat_real_type & get_mat_real( string_view = "" );

    //!
    //! Get the stored value as a const matrix of `real_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const matrix of `real_type`
    //!
    //! \code
    //! const GenericoContainer container;
    //! mat_real_type const &real_matrix = container.get_mat_real("In function get_mat_real");
    //! // Use real_matrix...
    //! \endcode
    //!
    mat_real_type const & get_mat_real( string_view = "" ) const;

    //!
    //! Get the stored value as a matrix of `complex_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a matrix of `complex_type`
    //!
    //! \code
    //! GenericoContainer container;
    //! mat_complex_type &complex_matrix = container.get_mat_complex("In function get_mat_complex");
    //! // Use complex_matrix...
    //! \endcode
    //!
    mat_complex_type & get_mat_complex( string_view = "" );

    //!
    //! Get the stored value as a const matrix of `complex_type`
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const matrix of `complex_type`
    //!
    //! \code
    //! const GenericoContainer container;
    //! mat_complex_type const &complex_matrix = container.get_mat_complex("In function get_mat_complex");
    //! // Use complex_matrix...
    //! \endcode
    //!
    mat_complex_type const & get_mat_complex( string_view = "" ) const;

    //!
    //! Get the stored value as a vector of strings
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a vector of strings
    //!
    //! \code
    //! GenericoContainer container;
    //! vec_string_type &string_data = container.get_vec_string("In function get_vec_string");
    //! // Use string_data...
    //! \endcode
    //!
    vec_string_type & get_vec_string( string_view = "" );

    //!
    //! Get the stored value as a const vector of strings
    //!
    //! \param[in]  where Position added to the error message
    //! \return The data stored in the container as a const vector of strings
    //!
    //! \code
    //! const GenericoContainer container;
    //! vec_string_type const &string_data = container.get_vec_string("In function get_vec_string");
    //! // Use string_data...
    //! \endcode
    //!
    vec_string_type const & get_vec_string( string_view = "" ) const;

    ///@}







    //!
    //! \name Access to vector type data and convert.
    //!
    ///@{

    //!
    //! Copy internal data to a vector of integers
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_int_type int_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_int(int_vector, "In function example_call");
    //! // Use int_vector...
    //! \endcode
    //!
    void copyto_vec_int( vec_int_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of unsigned integers
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_uint_type uint_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_uint(uint_vector, "In function example_call");
    //! // Use uint_vector...
    //! \endcode
    //!
    void copyto_vec_uint( vec_uint_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of long integers
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_long_type long_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_long(long_vector, "In function example_call");
    //! // Use long_vector...
    //! \endcode
    //!
    void copyto_vec_long( vec_long_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of unsigned long integers
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_ulong_type ulong_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_ulong(ulong_vector, "In function example_call");
    //! // Use ulong_vector...
    //! \endcode
    //!
    void copyto_vec_ulong( vec_ulong_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of `real_type`
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_real_type real_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_real(real_vector, "In function example_call");
    //! // Use real_vector...
    //! \endcode
    //!
    void copyto_vec_real( vec_real_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of `complex_type`
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_complex_type complex_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_complex(complex_vector, "In function example_call");
    //! // Use complex_vector...
    //! \endcode
    //!
    void copyto_vec_complex( vec_complex_type & v, string_view = "" ) const;

    //!
    //! Copy internal data to a vector of strings
    //!
    //! \param[out] v     Vector to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_string_type string_vector;
    //! GenericoContainer container;
    //! container.copyto_vec_string(string_vector, "In function example_call");
    //! // Use string_vector...
    //! \endcode
    //!
    void copyto_vec_string( vec_string_type & v, string_view = "" ) const;

    ///@}




    //!
    //! \name Access to matrix type data and convert.
    //!
    ///@{

    //!
    //! Copy internal data to a matrix of integers
    //!
    //! \param[out] m     Matrix to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! mat_int_type int_matrix;
    //! GenericoContainer container;
    //! container.copyto_mat_int(int_matrix, "In function example_call");
    //! // Use int_matrix...
    //! \endcode
    //!
    void copyto_mat_int( mat_int_type & m, string_view = "" ) const;

    //!
    //! Copy internal data to a matrix of long integers
    //!
    //! \param[out] m     Matrix to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! mat_long_type long_mat;
    //! GenericoContainer container;
    //! container.copyto_mat_long(long_mat, "In function example_call");
    //! // Use long_mat...
    //! \endcode
    //!
    void copyto_mat_long( mat_long_type & m, string_view = "" ) const;

    //!
    //! Copy internal data to a matrix of `real_type`
    //!
    //! \param[out] m     Matrix to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! mat_real_type real_mat;
    //! GenericoContainer container;
    //! container.copyto_mat_real(real_mat, "In function example_call");
    //! // Use real_mat...
    //! \endcode
    //!
    void copyto_mat_real( mat_real_type & m, string_view = "" ) const;

    //!
    //! Copy internal data to a matrix of `complex_type`
    //!
    //! \param[out] m     Matrix to store the data
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! vec_complex_type complex_mat;
    //! GenericoContainer container;
    //! container.copyto_vec_complex(complex_mat, "In function example_call");
    //! // Use complex_mat...
    //! \endcode
    //!
    void copyto_mat_complex( mat_complex_type & m, string_view = "" ) const;

    ///@}






    //!
    //! \name Access to element of vector type data
    //!
    ///@{


    //!
    //! If the `i`-th element of the vector is boolean,
    //! integer, or floating point then return the number, otherwise return `0`.
    //!
    //! \param[in] i     The position of the element in the vector
    //! \param[in] where Position added to the error message
    //! \return The value as `real_type`
    //!
    //! \code
    //! unsigned index = 2;
    //! real_type value = container.get_number_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    real_type get_number_at( unsigned i, string_view = "" ) const;

    //!
    //! If the `i`-th element of the vector is convertible to
    //! complex, return the number; otherwise return `0`.
    //!
    //! \param[in] i     The position of the element in the vector
    //! \param[in] where Position added to the error message
    //! \return The value as `complex_type`
    //!
    //! \code
    //! unsigned index = 2;
    //! complex_type value = container.get_complex_number_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    complex_type get_complex_number_at( unsigned i, string_view = "" ) const;

    //!
    //! If the `i`-th element of the vector is convertible to
    //! complex, return the number; otherwise return `0`.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \param[out] re Real part of the number
    //! \param[out] im Imaginary part of the number
    //! \param[in]  where Position added to the error message
    //!
    //! \code
    //! unsigned index = 2;
    //! real_type realPart, imagPart;
    //! container.get_complex_number_at(index, realPart, imagPart, "In function example_call");
    //! // Use realPart and imagPart...
    //! \endcode
    //!
    void get_complex_number_at( unsigned i, real_type & re, real_type & im, string_view = "" ) const;

    //!
    //! Get the `i`-th pointer of the vector of pointers.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return Reference to the pointer
    //!
    //! \code
    //! unsigned index = 1;
    //! auto& pointer = container.get_pointer_at<MyType>(index);
    //! // Use pointer...
    //! \endcode
    //!
    template <typename T>
    T& get_pointer_at( unsigned i )
    { return (*this)[i].get_pointer<T>(); }

    //!
    //! Get the `i`-th pointer of the vector of pointers.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored pointer
    //!
    //! \code
    //! unsigned index = 1;
    //! MyType pointer = container.get_pointer_at<MyType>(index);
    //! // Use pointer...
    //! \endcode
    //!
    template <typename T>
    T get_pointer_at( unsigned i ) const
    { return (*this)[i].get_pointer<T>(); }
    //!< Return `i`-th generic pointer (if fails issue an error).

    //!
    //! Get the `i`-th boolean of the stored data.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 3;
    //! bool_type value = container.get_bool_at(index);
    //! // Use value...
    //! \endcode
    //!
    bool_type get_bool_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    bool_type get_bool_at( T i )
    { return get_bool_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th boolean of the stored data.
    //!
    //! \param[in]  i     The position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 3;
    //! bool_type value = container.get_bool_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    bool_type get_bool_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th integer of the stored data.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 4;
    //! int_type &value = container.get_int_at(index);
    //! // Use value...
    //! \endcode
    //!
    int_type & get_int_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    int_type & get_int_at( T i )
    { return get_int_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const integer of the stored data.
    //!
    //! \param[in]  i     The position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 4;
    //! int_type const &value = container.get_int_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    int_type const & get_int_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th long integer of the stored data.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 5;
    //! long_type &value = container.get_long_at(index);
    //! // Use value...
    //! \endcode
    //!
    long_type & get_long_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    long_type & get_long_at( T i )
    { return get_long_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const long integer of the stored data.
    //!
    //! \param[in]  i     The position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 5;
    //! long_type const &value = container.get_long_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    long_type const & get_long_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th `real_type` of the stored data.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 6;
    //! real_type &value = container.get_real_at(index);
    //! // Use value...
    //! \endcode
    //!
    real_type & get_real_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    real_type & get_real_at( T i )
    { return get_real_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const `real_type` of the stored data.
    //!
    //! \param[in]  i     The position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 6;
    //! real_type const &value = container.get_real_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    real_type const & get_real_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th `complex_type` of the stored data.
    //!
    //! \param[in]  i  The position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 7;
    //! complex_type &value = container.get_complex_at(index);
    //! // Use value...
    //! \endcode
    //!
    complex_type & get_complex_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    complex_type & get_complex_at( T i )
    { return get_complex_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const `complex_type` of the stored data.
    //!
    //! \param[in]  i     The position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 7;
    //! complex_type const &value = container.get_complex_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    complex_type const & get_complex_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th integer of the stored data in a matrix.
    //!
    //! \param[in]  i  Row position of the element in the matrix
    //! \param[in]  j  Column position of the element in the matrix
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 0, col = 1;
    //! int_type &value = container.get_int_at(row, col);
    //! // Use value...
    //! \endcode
    //!
    int_type & get_int_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const integer of the stored data in a matrix.
    //!
    //! \param[in]  i     Row position of the element in the matrix
    //! \param[in]  j     Column position of the element in the matrix
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 0, col = 1;
    //! int_type const &value = container.get_int_at(row, col, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    int_type const & get_int_at( unsigned i, unsigned j, string_view const where ) const;

    //!
    //! Get the `i`-th long integer of the stored data in a matrix.
    //!
    //! \param[in]  i  Row position of the element in the matrix
    //! \param[in]  j  Column position of the element in the matrix
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 0, col = 2;
    //! long_type &value = container.get_long_at(row, col);
    //! // Use value...
    //! \endcode
    //!
    long_type & get_long_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const long integer of the stored data in a matrix.
    //!
    //! \param[in]  i     Row position of the element in the matrix
    //! \param[in]  j     Column position of the element in the matrix
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 0, col = 2;
    //! long_type const &value = container.get_long_at(row, col, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    long_type const & get_long_at( unsigned i, unsigned j, string_view const where ) const;

    //!
    //! Get the `i`-th `real_type` of the stored data in a matrix.
    //!
    //! \param[in]  i  Row position of the element in the matrix
    //! \param[in]  j  Column position of the element in the matrix
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 1, col = 1;
    //! real_type &value = container.get_real_at(row, col);
    //! // Use value...
    //! \endcode
    //!
    real_type & get_real_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const `real_type` of the stored data in a matrix.
    //!
    //! \param[in]  i     Row position of the element in the matrix
    //! \param[in]  j     Column position of the element in the matrix
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 1, col = 1;
    //! real_type const &value = container.get_real_at(row, col, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    real_type const & get_real_at( unsigned i, unsigned j, string_view const where ) const;

    //!
    //! Get the `i`-th `complex_type` of the stored data in a matrix.
    //!
    //! \param[in]  i  Row position of the element in the matrix
    //! \param[in]  j  Column position of the element in the matrix
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 1, col = 2;
    //! complex_type &value = container.get_complex_at(row, col);
    //! // Use value...
    //! \endcode
    //!
    complex_type & get_complex_at( unsigned i, unsigned j );

    //!
    //! Get the `i`-th const `complex_type` of the stored data in a matrix.
    //!
    //! \param[in]  i     Row position of the element in the matrix
    //! \param[in]  j     Column position of the element in the matrix
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned row = 1, col = 2;
    //! complex_type const &value = container.get_complex_at(row, col, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    complex_type const & get_complex_at( unsigned i, unsigned j, string_view const where ) const;

    //!
    //! Get the `i`-th string of the stored data.
    //!
    //! \param[in]  i  Position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 0;
    //! string_type &value = container.get_string_at(index);
    //! // Use value...
    //! \endcode
    //!
    string_type & get_string_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    string_type & get_string_at( T i )
    { return get_string_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const string of the stored data.
    //!
    //! \param[in]  i     Position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 0;
    //! string_view value = container.get_string_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    string_view get_string_at( unsigned i, string_view const where ) const;

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  Position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 0;
    //! GenericContainer &value = container.get_gc_at(index);
    //! // Use value...
    //! \endcode
    //!
    GenericContainer & get_gc_at( unsigned i );

    template <typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    GenericContainer & get_gc_at( T i )
    { return get_gc_at(static_cast<unsigned>(i)); }

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i     Position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! unsigned index = 0;
    //! GenericContainer const &value = container.get_gc_at(index, "In function example_call");
    //! // Use value...
    //! \endcode
    //!
    GenericContainer const & get_gc_at( unsigned i, string_view const where ) const;

    ///@}








    //!
    //! \name Access to map type element.
    //!
    ///@{

    //!
    //! Get the stored data as a map of `GenericContainer`.
    //!
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of get_map()
    //! try {
    //!     auto &myMap = container.get_map("In function example_call");
    //!     // Access elements in the map
    //!     auto value = myMap["key1"];
    //!     // Use value...
    //! } catch (const exception &e) {
    //!     cerr << "Error accessing map: " << e.what() << endl;
    //! }
    //! \endcode
    //!
    map_type & get_map( string_view = "" );

    //!
    //! Get the stored data as a const map of `GenericContainer`.
    //!
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of get_map() with const
    //! try {
    //!     const auto &myMap = container.get_map("In function example_call");
    //!     // Access elements in the const map
    //!     auto const value = myMap.at("key2");
    //!     // Use value...
    //! } catch (const exception &e) {
    //!     cerr << "Error accessing const map: " << e.what() << endl;
    //! }
    //! \endcode
    //!
    map_type const & get_map( string_view = "" ) const;

    ///@}



    //!
    //! \name Access using operators.
    //!
    ///@{

    //!
    //! Get the `i`-th `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  Position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator[] for non-const
    //! GenericContainer &containerItem = container[2];
    //! // Use containerItem...
    //! \endcode
    //!
    GenericContainer & operator [] ( unsigned i );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data.
    //!
    //! \param[in]  i  Position of the element in the vector
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator[] for const
    //! const GenericContainer &constContainerItem = container[2];
    //! // Use constContainerItem...
    //! \endcode
    //!
    GenericContainer const & operator [] ( unsigned i ) const;

    //!
    //! Get the `i`-th `GenericContainer` of the stored data using a string key.
    //!
    //! \param[in]  s  Key string of the element in the map
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator[] with string key
    //! GenericContainer &mapItem = container["myKey"];
    //! // Use mapItem...
    //! \endcode
    //!
    GenericContainer & operator [] ( string_view s );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data using a string key.
    //!
    //! \param[in]  s  Key string of the element in the map
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator[] with string key for const
    //! const GenericContainer &constMapItem = container["myKey"];
    //! // Use constMapItem...
    //! \endcode
    //!
    GenericContainer const & operator [] ( string_view s ) const;

    //!
    //! Get the `i`-th `GenericContainer` of the stored data with error message.
    //!
    //! \param[in]  i     Position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator() with error message
    //! try {
    //!     GenericContainer &item = container(i, "Accessing item at index 3");
    //!     // Use item...
    //! } catch (const exception &e) {
    //!     cerr << e.what() << endl;
    //! }
    //! \endcode
    //!
    GenericContainer & operator () ( unsigned i, string_view = "" );

    //!
    //! Get the `i`-th const `GenericContainer` of the stored data with error message.
    //!
    //! \param[in]  i     Position of the element in the vector
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator() with error message for const
    //! try {
    //!     const GenericContainer &constItem = container(i, "Accessing item at index 3");
    //!     // Use constItem...
    //! } catch (const exception &e) {
    //!     cerr << e.what() << endl;
    //! }
    //! \endcode
    //!
    GenericContainer const & operator () ( unsigned i, string_view = "" ) const;

    //!
    //! Get a `GenericContainer` in the stored data using a string key with error message.
    //!
    //! \param[in]  s     Key string of the element in the map
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator() with string key
    //! GenericContainer &item = container("myKey", "Accessing item with key 'myKey'");
    //! // Use item...
    //! \endcode
    //!
    GenericContainer & operator () ( string_view s, string_view = "" );

    //!
    //! Get a const `GenericContainer` in the stored data using a string key with error message.
    //!
    //! \param[in]  s     Key string of the element in the map
    //! \param[in]  where Position added to the error message
    //! \return The stored value
    //!
    //! \code
    //! // Example usage of operator() with string key for const
    //! const GenericContainer &constItem = container("myKey", "Accessing const item with key 'myKey'");
    //! // Use constItem...
    //! \endcode
    //!
    GenericContainer const & operator () ( string_view s, string_view = "" ) const;

    //!
    //! Get a `GenericContainer` in the stored data by searching for a matching key from a vector of keys.
    //!
    //! \param[in]  vs    Vector of keys strings
    //! \param[in]  where Position added to the error message
    //! \return The stored value of the first match
    //!
    //! \code
    //! // Example usage of operator() with vector of keys
    //! vec_string_type keys = {"key1", "key2", "key3"};
    //! GenericContainer &matchedItem = container(keys, "Searching for matching key");
    //! // Use matchedItem...
    //! \endcode
    //!
    GenericContainer & operator () ( vec_string_type const & vs, string_view = "" );

    //!
    //! Get a const `GenericContainer` in the stored data by searching for a matching key from a vector of keys.
    //!
    //! \param[in]  vs    Vector of keys string
    //! \param[in]  where Position added to the error message
    //! \return The stored value of the first match
    //!
    //! \code
    //! // Example usage of operator() with vector of keys for const
    //! const vec_string_type keys = {"key1", "key2", "key3"};
    //! const GenericContainer &constMatchedItem = container(keys, "Searching for matching key");
    //! // Use constMatchedItem...
    //! \endcode
    //!
    GenericContainer const & operator () ( vec_string_type const & vs, string_view = "" ) const;

    ///@}









    //!
    //! \name Initialize data using set command.
    //!
    ///@{

    //!
    //! Assign a boolean to the generic container.
    //!
    //! \param[in] a Boolean to be stored
    //!
    //! \code
    //! // Example usage of setting a boolean
    //! GenericContainer container;
    //! container.set(true); // Store true in the container
    //! \endcode
    //!
    void set( bool const & a ) { this->set_bool(a); }

    //!
    //! Assign an integer to the generic container.
    //!
    //! \param[in] a Integer to be stored
    //!
    //! \code
    //! // Example usage of setting an unsigned integer
    //! GenericContainer container;
    //! container.set(42u); // Store 42 as an unsigned integer in the container
    //! \endcode
    //!
    void set( uint_type const & a ) { this->set_int(int_type(a)); }

    //!
    //! Assign an integer to the generic container.
    //!
    //! \param[in] a Integer to be stored
    //!
    //! \code
    //! // Example usage of setting an integer
    //! GenericContainer container;
    //! container.set(-7); // Store -7 in the container
    //! \endcode
    //!
    void set( int_type const & a ) { this->set_int(a); }

    //!
    //! Assign an unsigned integer to the generic container.
    //!
    //! \param[in] a Unsigned integer to be stored
    //!
    //! \code
    //! // Example usage of setting an unsigned long integer
    //! GenericContainer container;
    //! container.set(100000ul); // Store 100000 as an unsigned long integer in the container
    //! \endcode
    //!
    void set( ulong_type const & a ) { this->set_long(long_type(a)); }

    //!
    //! Assign a long integer to the generic container.
    //!
    //! \param[in] a Long integer to be stored
    //!
    //! \code
    //! // Example usage of setting a long integer
    //! GenericContainer container;
    //! container.set(123456789L); // Store 123456789 as a long integer in the container
    //! \endcode
    //!
    void set( long_type const & a ) { this->set_long(a); }

    //!
    //! Assign a float to the generic container.
    //!
    //! \param[in] a Float to be stored
    //!
    //! \code
    //! // Example usage of setting a float
    //! GenericContainer container;
    //! container.set(3.14f); // Store 3.14 as a float in the container
    //! \endcode
    //!
    void set( float const & a ) { this->set_real(real_type(a)); }

    //!
    //! Assign a double to the generic container.
    //!
    //! \param[in] a Double to be stored
    //!
    //! \code
    //! // Example usage of setting a double
    //! GenericContainer container;
    //! container.set(2.7182818284); // Store 2.7182818284 as a double in the container
    //! \endcode
    //!
    void set( double const & a ) { this->set_real(real_type(a)); }

    //!
    //! Assign a complex float to the generic container.
    //!
    //! \param[in] a Complex float to be stored
    //!
    //! \code
    //! // Example usage of setting a complex float
    //! GenericContainer container;
    //! complex<float> complexFloat(1.0f, 2.0f);
    //! container.set(complexFloat); // Store complex number (1.0, 2.0) in the container
    //! \endcode
    //!
    void set( complex<float> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); }

    //!
    //! Assign a complex double to the generic container.
    //!
    //! \param[in] a Complex double to be stored
    //!
    //! \code
    //! // Example usage of setting a complex double
    //! GenericContainer container;
    //! complex<double> complexDouble(3.0, 4.0);
    //! container.set(complexDouble); // Store complex number (3.0, 4.0) in the container
    //! \endcode
    //!
    void set( complex<double> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); }

    //!
    //! Assign a string to the generic container.
    //!
    //! \param[in] a String to be stored
    //!
    //! \code
    //! // Example usage of setting a C-style string
    //! GenericContainer container;
    //! container.set("Hello, World!"); // Store "Hello, World!" in the container
    //! \endcode
    //!
    void set( char const * a ) { this->set_string(a); }

    //!
    //! Assign a string to the generic container.
    //!
    //! \param[in] a String to be stored
    //!
    //! \code
    //! // Example usage of setting a string
    //! GenericContainer container;
    //! container.set(string("Hello, C++!")); // Store "Hello, C++!" in the container
    //! \endcode
    //!
    void set( string_view a ) { this->set_string(a); }

    //!
    //! Assign a pointer to the generic container.
    //!
    //! \param[in] a Pointer to be stored
    //!
    //! \code
    //! // Example usage of setting a pointer
    //! GenericContainer container;
    //! int value = 42;
    //! container.set(&value); // Store the address of value in the container
    //! \endcode
    //!
    void set( pointer_type a ) { this->set_pointer(a); }

    ///@}





    //!
    //! \name Initialize data using operators.
    //!
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
    GenericContainer & operator = ( complex<float> const & a )
    { this->set_complex(real_type(a.real()),real_type(a.imag())); return * this; }

    //!
    //! Assign a complex of double to the generic container.
    //!
    //! \param[in] a complex of double to be stored
    //!
    GenericContainer & operator = ( complex<double> const & a )
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
    GenericContainer & operator = ( const char * a )
    { this->set_string(a); return * this; }

    //!
    //! Assign a string to the generic container.
    //!
    //! \param[in] a string to be stored
    //!
    GenericContainer & operator = ( string const & a )
    { this->set_string(a); return * this; }

    //!
    //! Assign a string to the generic container.
    //!
    //! \param[in] a string to be stored
    //!
    GenericContainer & operator = ( string_view a )
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
    { this->from_gc( a ); return * this; }

    //!
    //! Load a `GenericContainer` to the generic container (deep copy).
    //!
    //! \param[in] a `GenericContainer` to be stored
    //!
    void load( GenericContainer const & a ) { this->from_gc(a); }
    ///@}

    //!
    //! \name Promotion to a `bigger` data.
    //!
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

    //!
    //! \name Initialize data by overloading constructor.
    //!
    ///@{

    //!
    //! Construct a generic container storing a boolean
    //! \param[in] a initializer data
    //!
    GenericContainer( bool const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //!
    GenericContainer( uint_type const & a )
    : m_data_type(GC_type::NOTYPE) { *this = a; }

    //!
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //!
    GenericContainer( int_type const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //!
    GenericContainer( ulong_type const & a )
    : m_data_type(GC_type::NOTYPE) { *this = a; }

    //!
    //! Construct a generic container storing an integer
    //! \param[in] a initializer data
    //!
    GenericContainer( long_type const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a floating point number
    //! \param[in] a initializer data
    //!
    GenericContainer( float const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a floating point number
    //! \param[in] a initializer data
    //!
    GenericContainer( double const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a complex floating point number
    //! \param[in] a initializer data
    //!
    GenericContainer( complex<float> const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a complex floating point number
    //! \param[in] a initializer data
    //!
    GenericContainer( complex<double> const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a string or pointer
    //! \param[in] a initializer data
    //!
    GenericContainer( char const * a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a string or pointer
    //! \param[in] a initializer data
    //!
    GenericContainer( string const & a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a string or pointer
    //! \param[in] a initializer data
    //!
    GenericContainer( string_view a )
    : m_data_type(GC_type::NOTYPE) { this->operator=(a); }

    //!
    //! Construct a generic container storing a pointer
    //! \param[in] a initializer data
    //!
    GenericContainer( pointer_type a )
    : m_data_type(GC_type::NOTYPE) { this->set_pointer(a); }

    //!
    //! Construct a generic container copying container `gc`
    //! \param[in] gc initializer data
    //!
    GenericContainer( GenericContainer const & gc )
    : m_data_type(GC_type::NOTYPE) { this->from_gc(gc); }

    ///@}

    //!
    //! \name Utilities methods.
    //!
    ///@{

    //!
    //! Check if string `s` is a key of the stored map (if fails issue an error).
    //! \param[in] s key to be checked
    //!
    bool exists( string_view s ) const;

    //!
    //! Check if any string in `vs` is a key of the stored
    //! map (if fails issue an error).
    //! \param[in] vs vector of string with the keys to be checked
    //!
    bool exists( vec_string_type const & vs ) const;

    //!
    //! The data stored musty be a `map`.
    //! Search the
    //! \param[in] `vs` vector of string with the keys to be searched
    //! \param[in]  where position added to the error message
    //!
    string must_exists( vec_string_type const & vs, string_view const where ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, bool & value ) const;

    //!
    //! Check if vector of strings `fields` is a key of the
    //! stored map and extract value if exists
    //! \param[in] fields keys to be checked
    //! \param[in] value  value to be extracted
    //!
    bool get_if_exists( vec_string_type const & fields, bool & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, int_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, uint_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, long_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, ulong_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, real_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, complex_type & value ) const;

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    bool get_if_exists( string_view field, string_type & value ) const;

    //!
    //! Check if vector of strings `fields` is a key of the
    //! stored map and extract value if exists
    //! \param[in] fields keys to be checked
    //! \param[in] value  value to be extracted
    //!
    template <typename T>
    bool
    get_if_exists( vec_string_type const & fields, T & value ) const {
      for ( string_view field : fields ) {
        if ( get_if_exists( field, value ) ) return true;
      }
      return false;
    }

    //!
    //! Check if string `field` is a key of the stored map and extract value if exists
    //! \param[in] field key to be checked
    //! \param[in] value value to be extracted
    //!
    template <typename T>
    bool
    get_if_exists( char const field[], T & value ) const {
      string_type field_str(field);
      return get_if_exists( field_str, value );
    }

    ///@}

    //!
    //! \name I/O for GenericContainer objects.
    //!
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
      ostream_type & stream,
      string_view    prefix = "",
      string_view    indent = "    "
    ) const;

    //!
    //! Compare the contents of the object with `gc`
    //!
    //! \param[in] gc to object to compare
    //! \return a string with the first difference found
    //!
    string
    compare_content( GenericContainer const & gc, string_view from = "" ) const;

    //!
    //! Dump the contents of the object in a human readable way
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    void
    dump(
      ostream_type & stream,
      string_view    prefix = "",
      string_view    indent = "    "
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
      ostream_type & stream,
      string_view    prefix = "",
      string_view    indent = "    "
    ) const {
      this->dump( stream, prefix, indent );
    }

    //!
    //! Dump the contents of the object in a human readable way (aliad of dump)
    //!
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //! \param[in] indent indentation
    //!
    string_type
    print(
      string_view prefix = "",
      string_view indent = "    "
    ) const {
      ostringstream ostr;
      ostr.precision(stream_number_precision);
      this->print(ostr,prefix,indent);
      return ostr.str();
    }

    //!
    //! Copy the contents of the object into another object
    //!
    //! \param[out] gc output `GenericContainer`
    //!
    void to_gc( GenericContainer & gc ) const;

    //!
    //! Copy the contents of the object into another object
    //!
    //! \param[in] gc input `GenericContainer`
    //!
    void from_gc( GenericContainer const & gc );

    //!
    //! Merge two generic container
    //!
    //! \param[in] gc    input `GenericContainer`
    //! \param[in] where position added to the error message
    //!
    void merge( GenericContainer const & gc, string_view const where );

    //!
    //! Print the contents of the object in YAML syntax
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix string to be prepended to any field of the `GenericContainer`
    //!
    void to_yaml( ostream_type & stream, string_view prefix = "" ) const;

    //!
    //! Read the contents of stream in YAML syntax
    //!
    //! \param[in] stream input stream
    //! \return true if conversion successful
    //!
    bool from_yaml( istream_type & stream );

    //!
    //! Print the contents of the object in JSON syntax
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //!
    void to_json( ostream_type & stream, string_view prefix = "" ) const;

    //!
    //! Read the contents of stream in JSON syntax
    //!
    //! \param[in] stream input stream
    //! \return true if conversion successful
    //!
    bool from_json( istream_type & stream );
    bool from_json2( istream_type & stream );

    //!
    //! Print the contents of the object in TOML syntax
    //!
    //! \param[in] stream output stream
    //! \param[in] prefix strig to be prepended to any field of the `GenericContainer`
    //!
    bool to_toml( ostream_type & stream ) const;

    //!
    //! Read the contents of stream in JSON syntax
    //!
    //! \param[in] stream input stream
    //! \return true if conversion successful
    //!
    bool from_toml( istream_type & stream );

    //!
    //! Collapse heterogeneous vectors into a unified type.
    //! Attempts to collapse nested vectors into a matrix when possible.
    //!
    void collapse();

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
    write_formatted_data( ostream_type & stream, char const delimiter = '\t' ) const;

    //!
    //! \deprecated use `write_formatted_data`
    //!
    GenericContainer const &
    writeFormattedData( ostream_type & stream, char const delimiter = '\t' ) const {
      return this->write_formatted_data( stream, delimiter );
    }

    //!
    //! Read regular formatted data from `stream` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map
    //! which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `vector_type` which contains the vectors which are the
    //!              columns of the data red of type `vec_real_type`.
    //!
    //! \param stream       stream to write the output
    //! \param commentChars lines beginning with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   characters used as delimiter for headers
    //!
    GenericContainer &
    read_formatted_data(
      istream_type & stream,
      char const commentChars[] = "#%",
      char const delimiters[]   = " \t"
    );

    //!
    //! \deprecated use `read_formatted_data`
    //!
    GenericContainer &
    readFormattedData(
      istream_type & stream,
      char const commentChars[] = "#%",
      char const delimiters[]   = " \t"
    ) {
      return read_formatted_data( stream, commentChars, delimiters );
    }

    //!
    //! Read regular formatted data from file `fname` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `vector_type` which contains the vectors which are the
    //!              columns of the data red of type `vec_real_type`.
    //!
    //! \param fname        file name to be read
    //! \param commentChars lines beginning with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   characters used as delimiter for headers
    //!
    GenericContainer &
    read_formatted_data(
      char const fname[],
      char const commentChars[] = "#%",
      char const delimiters[]   = " \t"
    );

    //!
    //! \deprecated use `read_formatted_data`
    //!
    GenericContainer &
    readFormattedData(
      char const fname[],
      char const commentChars[] = "#%",
      char const delimiters[]   = " \t"
    ) {
      return read_formatted_data( fname, commentChars, delimiters );
    }

    //!
    //! Read regular formatted data from `stream` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map
    //! which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `map_type` which contains the vectors which are the
    //!              columns of the data red of type `vec_real_type`.
    //!
    //! \param stream       stream to write the output
    //! \param commentChars lines beginning with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   characters used as delimiter for headers
    //! \param ptr_pars     pointer to a `GenericContainer` which store poarameter parsed in the comment part of the file
    //!
    GenericContainer &
    read_formatted_data2(
      istream_type   & stream,
      char const       commentChars[] = "#%",
      char const       delimiters[]   = " \t",
      GenericContainer ptr_pars[]     = nullptr
    );

    //!
    //! \deprecated use `read_formatted_data2`
    //!
    GenericContainer &
    readFormattedData2(
      istream_type   & stream,
      char const       commentChars[] = "#%",
      char const       delimiters[]   = " \t",
      GenericContainer ptr_pars[]     = nullptr
    ) {
      return read_formatted_data2( stream, commentChars, delimiters, ptr_pars );
    }

    //!
    //! Read regular formatted data from file `fname` to `GenericContainer`.
    //!
    //! After successful read `GenericContainer` will be a map which contains the fields:
    //!
    //! - "headers"  a `vec_string_type` which contains
    //!              the strings of the headers of the columns of the data
    //!
    //! - "data"     a `map_type` which contains the vectors which are the
    //!              columns of the data red of type `vec_real_type`.
    //!
    //! \param fname        file name to be read
    //! \param commentChars lines beginning with one of this chars are treated as comments.
    //!                     Default are `#` and `%`
    //! \param delimiters   characters used as delimiter for headers
    //! \param ptr_pars     pointer to a `GenericContainer` which store the parameter parsed in the comment part of the file
    //!
    GenericContainer &
    read_formatted_data2(
      char const       fname[],
      char const       commentChars[] = "#%",
      char const       delimiters[]   = " \t",
      GenericContainer ptr_pars[]     = nullptr
    );

    //!
    //! \deprecated use `read_formatted_data2`
    //!
    GenericContainer &
    readFormattedData2(
      char const       fname[],
      char const       commentChars[] = "#%",
      char const       delimiters[]   = " \t",
      GenericContainer ptr_pars[]     = nullptr
    ) {
      return read_formatted_data2( fname, commentChars, delimiters, ptr_pars );
    }
    ///@}

    //!
    //! \name Serialization.
    //!
    ///@{

    //!
    //! Returns the size, in bytes, required to store the serialized version of the GenericContainer.
    //! This size represents the memory footprint of the container when serialized.
    //!
    int32_t mem_size() const;

    //!
    //! Serializes the GenericContainer into the provided buffer.
    //! The buffer must have enough space to store the serialized data (use \ref mem_size() to determine the required size).
    //!
    //! \param buffer_dim Size of the provided buffer, in bytes.
    //! \param buffer Pointer to the buffer where the serialized data will be written.
    //! \return The number of bytes written into the buffer, or -1 if an error occurs (e.g., insufficient buffer size).
    //!
    int32_t serialize( int32_t buffer_dim, uint8_t * buffer ) const;

    //!
    //! Serializes the GenericContainer into the provided vector buffer.
    //! The buffer will be resized if necessary to accommodate the serialized data.
    //!
    //! \param buffer A vector of bytes that will be filled with the serialized data.
    //! \return The number of bytes written into the buffer.
    //!
    int32_t serialize( vector<uint8_t> & buffer ) const;

    //!
    //! Deserializes the GenericContainer from the provided buffer.
    //! Reconstructs the container from its serialized version contained within the buffer.
    //!
    //! \param buffer_dim Size of the buffer containing the serialized data, in bytes.
    //! \param buffer Pointer to the buffer containing the serialized data.
    //! \return The number of bytes readed from the buffer.
    //!
    int32_t de_serialize( int32_t buffer_dim, uint8_t const * buffer );

    //!
    //! Deserializes the GenericContainer from the provided vector buffer.
    //! Reconstructs the container from its serialized version contained within the vector.
    //!
    //! \param buffer A vector of bytes containing the serialized data.
    //! \return  The number of bytes readed from the buffer.
    //!
    int32_t de_serialize( vector<uint8_t> const & buffer );

    ///@}

    //!
    //! Do an exception
    //!
    //! \param[in]  where position added to the error message
    //!
    static
    void
    exception( string_view const where ) GC_NO_RETURN;

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
  write_table(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  );

  //!
  //! \deprecated use `write_table`
  //!
  inline
  void
  writeTable(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  ) {
    write_table( headers, data, stream, delimiter );
  }

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //! \param[in] delimiter delimiter character between columns
  //!
  void
  write_table(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  );

  //!
  //! \deprecated use `write_table`
  //!
  inline
  void
  writeTable(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream,
    char delimiter = '\t'
  ) {
    write_table( headers, data, stream, delimiter );
  }

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //!
  void
  write_table_formatted(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream
  );

  //!
  //! \deprecated use `write_table_formatted`
  //!
  inline
  void
  writeTableFormatted(
    vec_string_type const & headers,
    vector_type     const & data,
    ostream_type          & stream
  ) {
    write_table_formatted( headers, data, stream );
  }

  //!
  //! Write data as a table
  //!
  //! \param[in] headers   vector of string with the header of the table
  //! \param[in] data      matrix of `real_type` with the columns of the table
  //! \param[in] stream    output stream
  //!
  void
  write_table_formatted(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream
  );

  //!
  //! \deprecated use `write_table_formatted`
  //!
  inline
  void
  writeTableFormatted(
    vec_string_type const & headers,
    mat_real_type   const & data,
    ostream_type          & stream
  ) {
    write_table_formatted( headers, data, stream );
  }

  //!
  //! Utrility to write sctring escaping non printable
  //!
  void string_escape( ostream_type & stream, string const & s );

}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  // do not define alias GC if use X11
  #ifndef XlibSpecificationRelease
  namespace GC = GC_namespace;
  #endif

  // for backward compatibility
  namespace GenericContainerNamespace = GC_namespace;

#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

//
// eof: GenericContainer.hh
//
