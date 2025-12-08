
/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2021                                                      |
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
 |      Università degli Studi di Trento                                    |
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <string_view>
#include <typeinfo>
#include <vector>

#include "mex.h"

#define CLASS_HANDLE_SIGNATURE 0xFF00F0A5

#define arg_in_0 prhs[0]
#define arg_in_1 prhs[1]
#define arg_in_2 prhs[2]
#define arg_in_3 prhs[3]
#define arg_in_4 prhs[4]
#define arg_in_5 prhs[5]
#define arg_in_6 prhs[6]
#define arg_in_7 prhs[7]
#define arg_in_8 prhs[8]
#define arg_in_9 prhs[9]
#define arg_in_10 prhs[10]
#define arg_in_11 prhs[11]
#define arg_in_12 prhs[12]
#define arg_in_13 prhs[13]
#define arg_in_14 prhs[14]
#define arg_in_15 prhs[15]
#define arg_in_16 prhs[16]
#define arg_in_17 prhs[17]
#define arg_in_18 prhs[18]
#define arg_in_19 prhs[19]

#define arg_out_0 plhs[0]
#define arg_out_1 plhs[1]
#define arg_out_2 plhs[2]
#define arg_out_3 plhs[3]
#define arg_out_4 plhs[4]
#define arg_out_5 plhs[5]
#define arg_out_6 plhs[6]
#define arg_out_7 plhs[7]
#define arg_out_8 plhs[8]
#define arg_out_9 plhs[9]
#define arg_out_10 plhs[10]
#define arg_out_11 plhs[11]
#define arg_out_12 plhs[12]
#define arg_out_13 plhs[13]
#define arg_out_14 plhs[14]
#define arg_out_15 plhs[15]
#define arg_out_16 plhs[16]
#define arg_out_17 plhs[17]
#define arg_out_18 plhs[18]
#define arg_out_19 plhs[19]

#define UTILS_MEX_ASSERT0( COND, MSG )                                                                                 \
  if ( !( COND ) ) Utils::mex_error_message( MSG )

#define UTILS_MEX_ASSERT( COND, FMT, ... ) UTILS_MEX_ASSERT0( COND, fmt::format( FMT, __VA_ARGS__ ) )

// -----------------------------------------------------------------------------

namespace Utils
{

  using std::string;
  using std::string_view;
  using std::vector;

  /*!
   * \addtogroup Mex
   * @{
   */

  //!
  //! \brief Sends an error message to MATLAB.
  //!
  //! \param msg The error message to display.
  //!
  inline void
  mex_error_message( string_view msg )
  {
    mexErrMsgTxt( msg.data() );
  }

  //!
  //! \brief Checks if the input argument is a scalar.
  //!
  //! \param arg The mxArray input argument.
  //! \param msg The error message to display if not scalar.
  //! \return True if arg is a scalar, false otherwise.
  //!
  inline bool
  mex_is_scalar( mxArray const * arg, string_view msg )
  {
    mwSize number_of_dimensions = mxGetNumberOfDimensions( arg );
    UTILS_MEX_ASSERT0( number_of_dimensions == 2, msg );
    mwSize const * dims{ mxGetDimensions( arg ) };
    return dims[0] == 1 && dims[1] == 1;
  }

  //!
  //! \brief Gets the scalar value from the input argument.
  //!
  //! \param arg The mxArray input argument.
  //! \param msg The error message to display if not scalar.
  //! \return The scalar value.
  //!
  inline double
  mex_get_scalar_value( mxArray const * arg, string_view msg )
  {
    mwSize number_of_dimensions = mxGetNumberOfDimensions( arg );
    UTILS_MEX_ASSERT0( number_of_dimensions == 2, msg );
    mwSize const * dims{ mxGetDimensions( arg ) };
    UTILS_MEX_ASSERT( dims[0] == 1 && dims[1] == 1, "{}, found {} x {} matrix\n", msg, dims[0], dims[1] );
    return mxGetScalar( arg );
  }

  //!
  //! \brief Gets a boolean value from the input argument.
  //!
  //! \param arg The mxArray input argument.
  //! \param msg The error message to display if not a logical scalar.
  //! \return The boolean value.
  //!
  inline bool
  mex_get_bool( mxArray const * arg, string_view msg )
  {
    UTILS_MEX_ASSERT0( mxIsLogicalScalar( arg ), msg );
    return mxIsLogicalScalarTrue( arg );
  }

  //!
  //! \brief Gets a 64-bit integer value from the input argument.
  //!
  //! \param arg The mxArray input argument.
  //! \param msg The error message to display if not a valid scalar.
  //! \return The 64-bit integer value.
  //!
  inline int64_t
  mex_get_int64( mxArray const * arg, string_view msg )
  {
    using std::floor;
    mwSize number_of_dimensions = mxGetNumberOfDimensions( arg );
    UTILS_MEX_ASSERT0( number_of_dimensions == 2, msg );
    mwSize const * dims = mxGetDimensions( arg );
    UTILS_MEX_ASSERT( dims[0] == 1 && dims[1] == 1, "{}, found {} x {} matrix\n", msg, dims[0], dims[1] );
    mxClassID category = mxGetClassID( arg );
    int64_t   res      = 0;
    void *    ptr      = mxGetData( arg );
    switch ( category )
    {
      case mxINT8_CLASS:
        res = *static_cast<int8_t *>( ptr );
        break;
      case mxUINT8_CLASS:
        res = *static_cast<uint8_t *>( ptr );
        break;
      case mxINT16_CLASS:
        res = *static_cast<int16_t *>( ptr );
        break;
      case mxUINT16_CLASS:
        res = *static_cast<uint16_t *>( ptr );
        break;
      case mxINT32_CLASS:
        res = *static_cast<int32_t *>( ptr );
        break;
      case mxUINT32_CLASS:
        res = *static_cast<uint32_t *>( ptr );
        break;
      case mxINT64_CLASS:
        res = *static_cast<int64_t *>( ptr );
        break;
      case mxUINT64_CLASS:
        res = *static_cast<uint64_t *>( ptr );
        break;
      case mxDOUBLE_CLASS:
      {
        double tmp{ *static_cast<double *>( ptr ) };
        UTILS_MEX_ASSERT( tmp == floor( tmp ), "{} expected int, found {}\n", msg, tmp );
        res = static_cast<int64_t>( tmp );
      }
      break;
      case mxSINGLE_CLASS:
      {
        float tmp{ *static_cast<float *>( ptr ) };
        UTILS_MEX_ASSERT( tmp == floor( tmp ), "{} expected int, found {}\n", msg, tmp );
        res = static_cast<int64_t>( tmp );
      }
      break;
      default:
        UTILS_MEX_ASSERT( false, "{} bad type scalar", msg );
        break;
    }
    return res;
  }

  //!
  //! \brief Gets a pointer to a vector from the input argument.
  //!
  //! \param arg The mxArray input argument.
  //! \param sz The size of the vector (output).
  //! \param msg The error message to display if not a valid vector.
  //! \return Pointer to the vector data.
  //!
  inline double const *
  mex_vector_pointer( mxArray const * arg, mwSize & sz, string_view msg )
  {
    mwSize number_of_dimensions = mxGetNumberOfDimensions( arg );
    UTILS_MEX_ASSERT0( number_of_dimensions == 2, msg );
    mwSize const * dims{ mxGetDimensions( arg ) };
    UTILS_MEX_ASSERT( dims[0] == 1 || dims[1] == 1 || dims[0] * dims[1] == 0,
                      "{}\nExpect (1 x n or n x 1 or empty) matrix, found {} x {}\n", msg, dims[0], dims[1] );
    sz = dims[0] * dims[1];
    return mxGetPr( arg );
  }

  //!
  //! \brief Gets a pointer to a matrix from the input argument.
  //!
  //! \param arg The mxArray input argument.
  //! \param nr The number of rows (output).
  //! \param nc The number of columns (output).
  //! \param msg The error message to display if not a valid matrix.
  //! \return Pointer to the matrix data.
  //!
  inline double const *
  mex_matrix_pointer( mxArray const * arg, mwSize & nr, mwSize & nc, string_view msg )
  {
    mwSize number_of_dimensions = mxGetNumberOfDimensions( arg );
    UTILS_MEX_ASSERT0( number_of_dimensions == 2, msg );
    mwSize const * dims{ mxGetDimensions( arg ) };
    nr = dims[0];
    nc = dims[1];
    return mxGetPr( arg );
  }

  // -----------------------------------------------------------------------------

  //!
  //! \brief Sets a scalar value in the output argument.
  //!
  //! \param arg Reference to the mxArray output argument.
  //! \param value The scalar value to set.
  //!
  inline void
  mex_set_scalar_value( mxArray *& arg, double value )
  {
    arg             = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
    *mxGetPr( arg ) = value;
  }

  //!
  //! \brief Sets a scalar integer value in the output argument.
  //!
  //! \param arg Reference to the mxArray output argument.
  //! \param value The integer value to set.
  //!
  inline void
  mex_set_scalar_int32( mxArray *& arg, int32_t value )
  {
    arg                                         = mxCreateNumericMatrix( 1, 1, mxINT32_CLASS, mxREAL );
    *static_cast<int32_t *>( mxGetData( arg ) ) = value;
  }

  //!
  //! \brief Sets a scalar 64-bit integer value in the output argument.
  //!
  //! \param arg Reference to the mxArray output argument.
  //! \param value The 64-bit integer value to set.
  //!
  inline void
  mex_set_scalar_int64( mxArray *& arg, int64_t value )
  {
    arg                                         = mxCreateNumericMatrix( 1, 1, mxINT64_CLASS, mxREAL );
    *static_cast<int64_t *>( mxGetData( arg ) ) = value;
  }

  //!
  //! \brief Sets a boolean value in the output argument.
  //!
  //! \param arg Reference to the mxArray output argument.
  //! \param value The boolean value to set.
  //!
  inline void
  mex_set_scalar_bool( mxArray *& arg, bool value )
  {
    arg = mxCreateLogicalScalar( value );
  }

  //!
  //! \brief Creates a numeric matrix of type int32 and returns a pointer to its
  //! data.
  //!
  //! This function allocates memory for a matrix of specified size (nrow x
  //! ncol) and initializes it as an int32 numeric matrix. The created matrix is
  //! assigned to the output argument `arg`, and a pointer to the data is
  //! returned.
  //!
  //! \param arg Reference to the mxArray pointer that will hold the created
  //! matrix.
  //! \param nrow Number of rows in the matrix.
  //! \param ncol Number of columns in the matrix.
  //! \return Pointer to the data of the created int32 matrix.
  //!
  inline int32_t *
  mex_create_matrix_int32( mxArray *& arg, mwSize nrow, mwSize ncol )
  {
    arg = mxCreateNumericMatrix( nrow, ncol, mxINT32_CLASS, mxREAL );
    return static_cast<int32_t *>( mxGetData( arg ) );
  }


  //!
  //! \brief Creates a numeric matrix of type int64 and returns a pointer to its
  //! data.
  //!
  //! This function allocates memory for a matrix of specified size (nrow x
  //! ncol) and initializes it as an int64 numeric matrix. The created matrix is
  //! assigned to the output argument `arg`, and a pointer to the data is
  //! returned.
  //!
  //! \param arg Reference to the mxArray pointer that will hold the created
  //! matrix.
  //! \param nrow Number of rows in the matrix.
  //! \param ncol Number of columns in the matrix.
  //! \return Pointer to the data of the created int64 matrix.
  //!
  inline int64_t *
  mex_create_matrix_int64( mxArray *& arg, mwSize nrow, mwSize ncol )
  {
    arg = mxCreateNumericMatrix( nrow, ncol, mxINT64_CLASS, mxREAL );
    return static_cast<int64_t *>( mxGetData( arg ) );
  }


  //!
  //! \brief Creates a numeric matrix of type double and returns a pointer to
  //! its data.
  //!
  //! This function allocates memory for a matrix of specified size (nrow x
  //! ncol) and initializes it as a double numeric matrix. The created matrix is
  //! assigned to the output argument `arg`, and a pointer to the data is
  //! returned.
  //!
  //! \param arg Reference to the mxArray pointer that will hold the created
  //! matrix.
  //! \param nrow Number of rows in the matrix.
  //! \param ncol Number of columns in the matrix.
  //! \return Pointer to the data of the created double matrix.
  //!
  inline double *
  mex_create_matrix_value( mxArray *& arg, mwSize nrow, mwSize ncol )
  {
    arg = mxCreateNumericMatrix( nrow, ncol, mxDOUBLE_CLASS, mxREAL );
    return mxGetPr( arg );
  }

  // -----------------------------------------------------------------------------
  //!
  //!  \brief Creates a sparse matrix in MATLAB format.
  //!
  //!  This function creates a sparse matrix using MATLAB's internal data
  //!  structures. The inputs are arrays representing the row indices, column
  //!  indices, and values of the non-zero elements of the sparse matrix.
  //!
  //!  \tparam R      Type of the values in the sparse matrix.
  //!  \tparam I      Type of the row and column indices.
  //!  \param  nnz    Number of non-zero elements in the sparse matrix.
  //!  \param  nrows  Number of rows in the sparse matrix.
  //!  \param  ncols  Number of columns in the sparse matrix.
  //!  \param  i_rows Array of row indices (0-based).
  //!  \param  j_cols Array of column indices (0-based).
  //!  \param  vals   Array of non-zero values.
  //!  \return int    Status code of the MATLAB function call.
  //!
  template <typename R, typename I>
  inline int
  mex_create_sparse_matrix( size_t    nnz,
                            size_t    nrows,
                            size_t    ncols,
                            I         i_rows[],
                            I         j_cols[],
                            R         vals[],
                            mxArray * arg_out[] )
  {
    mxArray * args[5];  // Array of arguments to be passed to MATLAB's sparse
                        // function.

    // Create MATLAB matrices to store row indices, column indices, and values.
    double * Irow{ mex_create_matrix_value( args[0], 1, nnz ) };
    double * Jcol{ mex_create_matrix_value( args[1], 1, nnz ) };
    double * VALS{ mex_create_matrix_value( args[2], 1, nnz ) };

    // Set the number of rows and columns in the sparse matrix.
    mex_set_scalar_value( args[3], nrows );
    mex_set_scalar_value( args[4], ncols );

    // Convert the row and column indices to 1-based indexing for MATLAB.
    std::transform( i_rows, i_rows + nnz, Irow, []( I val ) -> double { return double( val + 1 ); } );
    std::transform( j_cols, j_cols + nnz, Jcol, []( I val ) -> double { return double( val + 1 ); } );

    // Copy the values of the non-zero elements into the MATLAB array.
    std::copy_n( vals, nnz, VALS );

    // Call the MATLAB function 'sparse' to create the sparse matrix.
    return mexCallMATLAB( 1, arg_out, 5, args, "sparse" );
  }


  //!
  //! \brief Creates a MATLAB cell array and fills it with a vector of C++
  //! strings.
  //!
  //! \param arg Reference to the mxArray pointer that will hold the MATLAB cell
  //! array.
  //! \param str_vec Vector of C++ strings to be inserted into the cell array.
  //!

  inline void
  mex_create_string_cell_array( mxArray *& arg, vector<string> const & str_vec )
  {
    // Crea un cell array MATLAB con lo stesso numero di elementi del vettore di
    // stringhe
    arg = mxCreateCellMatrix( str_vec.size(), 1 );

    // Riempie il cell array con le stringhe C++
    for ( size_t i{ 0 }; i < str_vec.size(); ++i )
    {
      mxArray * str{ mxCreateString( str_vec[i].data() ) };  // Crea una stringa MATLAB dalla stringa C++
      mxSetCell( arg, i,
                 str );  // Imposta la stringa nella posizione corretta del cell array
    }
  }

  /*
  Class Handle by Oliver Woodford

  https://it.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class
  */

  //!
  //! \brief A class template that wraps and manages a C++ object for use in
  //! MATLAB.
  //!
  //! This handle safely encapsulates a pointer to a C++ object for use with MEX
  //! functions. It provides automatic memory management and type-safe access,
  //! and ensures validity through a unique signature and RTTI type checking.
  //!
  //! \tparam base The C++ class type being wrapped.
  //!
  template <typename base>
  class mex_class_handle
  {
    uint32_t m_signature{ CLASS_HANDLE_SIGNATURE };  ///< Signature used to verify handle validity
    base *   m_ptr{ nullptr };                       ///< Pointer to the managed C++ object..
    string   m_name;                                 ///< Name of the C++ class type (RTTI).

  public:
    /// Deleted copy assignment operator to prevent copying
    mex_class_handle & operator=( const mex_class_handle & ) = delete;

    /// Deleted default constructor to enforce explicit initialization
    mex_class_handle() = delete;

    //!
    //! \brief Constructs a handle for the given C++ object.
    //!
    //! \param ptr Pointer to the object to be managed.
    //!
    explicit mex_class_handle( base * ptr ) : m_ptr( ptr ), m_name( typeid( base ).name() ) {}

    //!
    //! \brief Destructor that deletes the managed object.
    //!
    ~mex_class_handle()
    {
      m_signature = 0;
      delete m_ptr;
      m_ptr = nullptr;
    }

    //!
    //! \brief Verifies the integrity and type of the handle.
    //!
    //! \return true if the handle is valid and matches the expected type.
    //!
    bool
    is_valid() const
    {
      return m_signature == CLASS_HANDLE_SIGNATURE && m_name == typeid( base ).name();
    }

    //!
    //! \brief Returns a pointer to the managed object.
    //!
    //! \return Pointer to the wrapped C++ object.
    //!
    base *
    ptr() const
    {
      return m_ptr;
    }
  };


  //!
  //! \brief Converts a C++ pointer into a MATLAB mxArray handle.
  //!
  //! Allocates a `mxUINT64_CLASS` array to store a pointer to a
  //! `mex_class_handle` wrapper. Also locks the MEX file to prevent unloading
  //! while the handle is in use.
  //!
  //! \tparam base The C++ class type being wrapped.
  //! \param ptr Pointer to the object to wrap.
  //! \return mxArray* containing the handle.
  //!
  template <typename base>
  inline mxArray *
  mex_convert_ptr_to_mx( base * ptr )
  {
    mexLock();  // Impedisce che MEX venga scaricato mentre il puntatore è
                // attivo

    // Alloca una matrice numerica di 2x1 uint64, anche se usi solo il primo
    // elemento
    mxArray * out{ mxCreateNumericMatrix( 2, 1, mxUINT64_CLASS, mxREAL ) };

    // Converte l'indirizzo del puntatore in uintptr_t e lo salva nell'array
    uintptr_t handle = reinterpret_cast<uintptr_t>( new mex_class_handle<base>( ptr ) );
    std::memcpy( mxGetData( out ), &handle, sizeof( uintptr_t ) );

    return out;
  }

  //!
  //! \brief Converts a MATLAB mxArray back to a mex_class_handle pointer.
  //!
  //! This function retrieves the handle from a MATLAB mxArray, checking
  //! that it is a valid uint64 scalar.
  //!
  //! \tparam base The type of the C++ class being wrapped.
  //! \param in The mxArray containing the handle.
  //! \return Pointer to the mex_class_handle.
  //! \throws std::runtime_error if the input is not a valid uint64 scalar.
  //!
  template <typename base>
  inline mex_class_handle<base> *
  mex_convert_mx_to_handle_ptr( const mxArray * in )
  {
    if ( mxGetNumberOfElements( in ) != 2 || mxGetClassID( in ) != mxUINT64_CLASS || mxIsComplex( in ) )
    {
      mexErrMsgTxt( "Expected a real 2x1 uint64 array" );
    }

    uintptr_t handle;
    std::memcpy( &handle, mxGetData( in ), sizeof( uintptr_t ) );

    auto * ptr = reinterpret_cast<mex_class_handle<base> *>( handle );
    if ( !ptr->is_valid() ) mexErrMsgTxt( "Invalid handle" );

    return ptr;
  }

  //!
  //! \brief Extracts the managed C++ object pointer from an mxArray.
  //!
  //! \tparam base The C++ class type being wrapped.
  //! \param in mxArray containing the handle.
  //! \return Pointer to the original C++ object.
  //!
  template <typename base>
  inline base *
  mex_convert_mx_to_ptr( mxArray const * in )
  {
    return mex_convert_mx_to_handle_ptr<base>( in )->ptr();
  }

  //!
  //! \brief Destroys the object wrapped by the mex_class_handle.
  //!
  //! This function deletes the handle and its managed object and
  //! sets the mxArray reference to nullptr.
  //!
  //! \tparam base The type of the C++ class being wrapped.
  //! \param in Reference to the mxArray containing the handle.
  //!
  template <typename base>
  inline void
  mex_destroy_object( mxArray const *& in )
  {
    if ( in != nullptr ) delete mex_convert_mx_to_handle_ptr<base>( in );
    in = nullptr;
    mexUnlock();
  }

  /*! @} */

}  // namespace Utils

//
// eof: mex_utils.hxx
//
