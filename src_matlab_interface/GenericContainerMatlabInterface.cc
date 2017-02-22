/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
\****************************************************************************/

#include "GenericContainerMatlabInterface.hh"
#include "matrix.h"

#include <iostream>
#include <cstdint>

namespace GenericContainerNamespace {

  using namespace std ;

  // ===========================================================================
  void
  mexPrint( GenericContainer const & gc ) {
    std::ostringstream ss ;
    gc.print(ss) ;
    mexPrintf("%s\n", ss.str().c_str()) ;
  }

  // ===========================================================================

  static
  void
  mx_to_vec_bool( mxArray const * mx, GenericContainer & gc ) {
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc = mxIsLogicalScalarTrue(mx) ;
    } else {
      mxLogical const * pr = mxGetLogicals(mx);
      mwSize total_num_of_elements = mxGetNumberOfElements(mx);
      vector_type & vec = gc.set_vector( total_num_of_elements ) ;
      for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
        vec[index] = *pr++ ;
    }
  }

  template <typename T>
  static
  void
  mx_to_vec_int( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims = mxGetDimensions(mx) ;
    T * pr = (T*)mxGetData(mx);
    if ( mxGetNumberOfElements(mx) == 1 ) {
      gc.set_int( int_type(*pr) ) ;
    } else {
      mwSize number_of_dimensions = mxGetNumberOfDimensions(mx) ;
      switch ( number_of_dimensions ) {
        case 1:
        { mwSize total_num_of_elements = mxGetNumberOfElements(mx);
          vec_int_type & vec = gc.set_vec_int( total_num_of_elements ) ;
          for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
            vec[index] = int_type(*pr++) ;
        }
        break ;
        case 2:
          if ( dims[0] == 1 ) {
            mwSize total_num_of_elements = mxGetNumberOfElements(mx);
            vec_int_type & vec = gc.set_vec_int( total_num_of_elements ) ;
            for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
              vec[index] = int_type(*pr++) ;
          } else {
            mat_real_type & mat = gc.set_mat_real( dims[0], dims[1] ) ;
            for ( mwSize j = 0 ; j < dims[1] ; ++j )
              for ( mwSize i = 0 ; i < dims[0] ; ++i )
                mat(i,j) = *pr++ ;
          }
        break ;
        default:
          mexPrintf("number_of_dimensions = %d\n", number_of_dimensions ) ;
        break ;
      }
    }
  }

  // ===========================================================================

  template <typename T>
  static
  void
  mx_to_vec_real( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims = mxGetDimensions(mx) ;
    T * pr = (T*)mxGetData(mx);
    if ( mxGetNumberOfElements(mx) == 1 ) {
      T val = *pr ;
      gc = val ;
    } else {
      mwSize number_of_dimensions = mxGetNumberOfDimensions(mx) ;
      switch ( number_of_dimensions ) {
        case 1:
        { mwSize total_num_of_elements = mxGetNumberOfElements(mx);
          vec_real_type & vec = gc.set_vec_real( total_num_of_elements ) ;
          for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
            vec[index] = *pr++ ;
        }
        break ;
        case 2:
          if ( dims[0] == 1 ) {
            mwSize total_num_of_elements = mxGetNumberOfElements(mx);
            vec_real_type & vec = gc.set_vec_real( total_num_of_elements ) ;
            for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
              vec[index] = *pr++ ;
          } else {
            mat_real_type & mat = gc.set_mat_real( dims[0], dims[1] ) ;
            mwSize total_num_of_elements = mxGetNumberOfElements(mx);
            for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
              mat[index] = *pr++ ;
          }
        break ;
        default:
          mexPrintf("number_of_dimensions = %d\n", number_of_dimensions ) ;
        break ;
      }
    }
  }

  // ===========================================================================

  template <typename T>
  static
  void
  mx_to_vec_complex( mxArray const * mx, GenericContainer & gc ) {
    mwSize const * dims = mxGetDimensions(mx) ;
    T * pr = (T*)mxGetData(mx);
    T * pi = (T*)mxGetImagData(mx);

    if ( mxGetNumberOfElements(mx) == 1 ) {
      complex_type c(*pr,*pi) ;
      gc = c ;
    } else {
      mwSize number_of_dimensions = mxGetNumberOfDimensions(mx) ;
      switch ( number_of_dimensions ) {
        case 1:
        { mwSize total_num_of_elements = mxGetNumberOfElements(mx);
          vec_complex_type & vec = gc.set_vec_complex( total_num_of_elements ) ;
          for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
            vec[index] = complex_type(*pr++,*pi++) ;
        }
        break ;
        case 2:
          if ( dims[0] == 1 ) {
            mwSize total_num_of_elements = mxGetNumberOfElements(mx);
            vec_complex_type & vec = gc.set_vec_complex( total_num_of_elements ) ;
            for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
              vec[index] = complex_type(*pr++,*pi++) ;
          } else {
            mat_complex_type & mat = gc.set_mat_complex( dims[0], dims[1] ) ;
            mwSize total_num_of_elements = mxGetNumberOfElements(mx);
            for ( mwSize index = 0 ; index < total_num_of_elements ; ++index )
              mat[index] = complex_type(*pr++,*pi++) ;
          }
        break ;
        default:
          mexPrintf("number_of_dimensions = %d\n", number_of_dimensions ) ;
        break ;
      }
    }
  }

  // ===========================================================================

  static
  void
  mx_to_vector( mxArray const * mx, GenericContainer & gc ) {
    unsigned ne = mxGetNumberOfElements(mx) ;
    gc.set_vector(ne) ;
    for ( unsigned index = 0 ; index < ne ; ++index ) {
      GenericContainer & gc1 = gc[index] ;
      mxArray const * cell = mxGetCell(mx, index);
      mxArray_to_GenericContainer( cell, gc1 ) ;
    }
  }

  // ===========================================================================

  static
  void
  mx_to_map( mxArray const * mx, GenericContainer & gc ) {
    gc.set_map() ;
    unsigned numFields = mxGetNumberOfFields(mx) ;
    for ( unsigned ifield = 0 ; ifield < numFields ; ++ifield ) {
      char const * field_name = mxGetFieldNameByNumber(mx,ifield) ;
      GenericContainer & gc1 = gc[field_name] ;
      unsigned nr = mxGetM(mx) ;
      unsigned nc = mxGetN(mx) ;
      if ( nc == 1 && nr == 1 ) {
        mxArray const * mxField = mxGetFieldByNumber(mx,0,ifield);
        mxArray_to_GenericContainer( mxField, gc1 ) ;
      } else {
        gc1.set_vector(nr*nc) ;
        for ( unsigned i = 0 ; i < nr*nc ; ++i ) {
          mxArray const * mxField = mxGetFieldByNumber(mx,i,ifield);
          mxArray_to_GenericContainer( mxField, gc1[i] ) ;
        }
      }
    }
  }

  // ===========================================================================

  void
  mxSparse_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) {

    size_t const * irs = mxGetIr(mx);
    size_t const * jcs = mxGetJc(mx);
    int            nc  = mxGetN(mx);
    int            nnz = jcs[nc];

    vec_int_type & jc = gc["jc"].set_vec_int( nc+1 ) ;
    vec_int_type & ir = gc["ir"].set_vec_int( nnz ) ;

    for ( unsigned i = 0 ; i <= nc ; ++i ) jc[i] = jcs[i] ;
    for ( unsigned i = 0 ; i < nnz ; ++i ) ir[i] = irs[i] ;

    double * sr = mxGetPr(mx);
    if ( mxIsComplex( mx ) ) {
      double * si = mxGetPi(mx);
      vec_complex_type & val = gc["values"].set_vec_complex( nnz ) ;
      for ( unsigned i = 0 ; i < nnz ; ++i ) val[i] = complex_type(sr[i],si[i]) ;
    } else {
      vec_real_type & val = gc["values"].set_vec_real( nnz ) ;
      for ( unsigned i = 0 ; i < nnz ; ++i ) val[i] = sr[i] ;
    }
  }

  // ===========================================================================

  void
  mxArray_to_GenericContainer( mxArray const * mx, GenericContainer & gc ) {
    if ( mx == nullptr ) return ;
    mxClassID category = mxGetClassID(mx);
    //mexPrintf("\n\n\n%s\n\n\n",mxGetClassName(mx)) ;
    if ( category == mxCELL_CLASS ) {
      mx_to_vector(mx,gc) ;
    } else if ( category == mxSTRUCT_CLASS ) {
      mx_to_map(mx,gc) ;
    } else if ( category == mxCHAR_CLASS ) {
      gc = mxArrayToString(mx) ;
    } else if ( mxIsSparse(mx) ) {
      mxSparse_to_GenericContainer( mx, gc ) ;
    } else {
      if ( mxIsComplex(mx) ) {
        switch (category)  {
          case mxLOGICAL_CLASS: mx_to_vec_bool(mx,gc);              break;
          case mxINT8_CLASS:    mx_to_vec_complex<int8_t>(mx,gc);   break;
          case mxUINT8_CLASS:   mx_to_vec_complex<uint8_t>(mx,gc);  break;
          case mxINT16_CLASS:   mx_to_vec_complex<int16_t>(mx,gc);  break;
          case mxUINT16_CLASS:  mx_to_vec_complex<uint16_t>(mx,gc); break;
          case mxINT32_CLASS:   mx_to_vec_complex<int32_t>(mx,gc);  break;
          case mxUINT32_CLASS:  mx_to_vec_complex<uint32_t>(mx,gc); break;
          case mxINT64_CLASS:   mx_to_vec_int<int64_t>(mx,gc);      break;
          case mxUINT64_CLASS:  mx_to_vec_int<uint64_t>(mx,gc);     break;
          case mxSINGLE_CLASS:  mx_to_vec_complex<float>(mx,gc);    break;
          case mxDOUBLE_CLASS:  mx_to_vec_complex<double>(mx,gc);   break;
          default:
            mexPrintf("Complex Class ID = %d not converted!\n", category ) ;
          break;
        }
      } else {
        switch (category)  {
          case mxLOGICAL_CLASS: mx_to_vec_bool(mx,gc);          break;
          case mxINT8_CLASS:    mx_to_vec_int<int8_t>(mx,gc);   break;
          case mxUINT8_CLASS:   mx_to_vec_int<uint8_t>(mx,gc);  break;
          case mxINT16_CLASS:   mx_to_vec_int<int16_t>(mx,gc);  break;
          case mxUINT16_CLASS:  mx_to_vec_int<uint16_t>(mx,gc); break;
          case mxINT32_CLASS:   mx_to_vec_int<int32_t>(mx,gc);  break;
          case mxUINT32_CLASS:  mx_to_vec_int<uint32_t>(mx,gc); break;
          case mxINT64_CLASS:   mx_to_vec_int<int64_t>(mx,gc);  break;
          case mxUINT64_CLASS:  mx_to_vec_int<uint64_t>(mx,gc); break;
          case mxSINGLE_CLASS:  mx_to_vec_real<float>(mx,gc);   break;
          case mxDOUBLE_CLASS:  mx_to_vec_real<double>(mx,gc);  break;
          default:
            mexPrintf("Class ID = %d not converted!\n", category ) ;
          break;
        }
      }
    }
  }

  // ===========================================================================

  void
  GenericContainer_to_mxArray( GenericContainer const & gc, mxArray * & mx ) {
    static char const where[] = "in GenericContainer_to_mxArray" ; 
    mwSize dims[2] = {1,1} ;
    switch ( gc.get_type() ) {
      case GC_NOTYPE:
      case GC_POINTER:
      case GC_VEC_POINTER:
        mx = mxCreateDoubleMatrix(0,0,mxREAL) ;
      break;
      case GC_BOOL:
        {
          mxLogical val = gc.get_bool() ? 1 : 0 ;
          mx = mxCreateLogicalScalar(val) ;
        }
      break;
      case GC_INTEGER:
        mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL) ;
        *(mwSize *)mxGetData(mx) = gc.get_int() ;
      break;
      case GC_LONG:
        mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL) ;
        *(mwSize *)mxGetData(mx) = gc.get_long() ;
      break;
      case GC_REAL:
        mx = mxCreateDoubleScalar(gc.get_real()) ;
      break;
      case GC_COMPLEX:
        {
          mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxCOMPLEX) ;
          GC::real_type re, im ;
          gc.get_complex_number(re,im) ;
          *mxGetPr(mx) = re ;
          *mxGetPi(mx) = im ;
        }
      break;
      case GC_STRING:
        mx = mxCreateString( gc.get_string().c_str() ) ;
      break;
      case GC_VEC_BOOL:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateNumericArray(2,dims,mxLOGICAL_CLASS,mxREAL) ;
          mxLogical * ptr = (mxLogical*)mxGetData(mx) ;
          for ( mwSize i = 0 ; i < dims[1] ; ++i ) ptr[i] = gc.get_bool_at(i,where) ;
        }
      break;
      case GC_VEC_INTEGER:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL) ;
          mwSize * ptr = (mwSize*)mxGetData(mx) ;
          for ( mwSize i = 0 ; i < dims[1] ; ++i ) ptr[i] = gc.get_int_at(i,where) ;
        }
      break;
      case GC_VEC_LONG:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL) ;
          mwSize * ptr = (mwSize*)mxGetData(mx) ;
          for ( mwSize i = 0 ; i < dims[1] ; ++i ) ptr[i] = gc.get_long_at(i,where) ;
        }
      break;
      case GC_VEC_REAL:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL) ;
          double * ptr = mxGetPr(mx) ;
          for ( mwSize i = 0 ; i < dims[1] ; ++i ) ptr[i] = gc.get_real_at(i,where) ;
        }
      break;
      case GC_VEC_COMPLEX:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxCOMPLEX) ;
          double * ptr = mxGetPr(mx) ;
          double * pti = mxGetPi(mx) ;
          for ( mwSize i = 0 ; i < dims[1] ; ++i ) gc.get_complex_number_at(i,ptr[i],pti[i]) ;
        }
      break;
      case GC_VEC_STRING:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateCellMatrix(dims[0], dims[1]) ;
          for( mwSize i = 0 ; i < dims[1] ; ++i )
            mxSetCell(mx,i,mxCreateString( gc.get_string_at(i,where).c_str()));
        }
      break;
      case GC_MAT_REAL:
        {
          dims[0] = gc.get_numRows() ;
          dims[1] = gc.get_numCols() ;
          mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL) ;
          double * ptr = mxGetPr(mx) ;
          mwSize k = 0 ;
          for ( mwSize j = 0 ; j < dims[1] ; ++j )
            for ( mwSize i = 0 ; i < dims[0] ; ++i )
              ptr[k++] = gc.get_real_at(i,j,where) ;
        }
      break;
      case GC_MAT_COMPLEX:
        {
          dims[0] = gc.get_numRows() ;
          dims[1] = gc.get_numCols() ;
          mx = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxCOMPLEX) ;
          double * ptr = mxGetPr(mx) ;
          double * pti = mxGetPi(mx) ;
          mwSize k = 0 ;
          for ( mwSize j = 0 ; j < dims[1] ; ++j ) {
            for ( mwSize i = 0 ; i < dims[0] ; ++i ) {
              GC::complex_type val = gc.get_complex_at(i,j,where) ;
              ptr[k] = val.real() ;
              pti[k] = val.imag() ;
              ++k ;
            }
          }
        }
      break;
      case GC_VECTOR:
        {
          dims[1] = gc.get_num_elements() ;
          mx = mxCreateCellMatrix(dims[0], dims[1]) ;
          for( mwSize i = 0 ; i < dims[1] ; ++i ) {
            mxArray * mxi = nullptr ;
            GenericContainer_to_mxArray( gc[i], mxi ) ;
            if ( mxi != nullptr ) mxSetCell( mx, i, mxi ) ;
          }
        }
      break;
      case GC_MAP:
        {
          GC::map_type const & mappa = gc.get_map() ;
          std::vector<char const *> fieldnames ;
          int nfield = mappa.size() ;
          fieldnames.reserve(nfield) ;
          for ( GC::map_type::const_iterator im = mappa.begin() ; im != mappa.end() ; ++im )
            fieldnames.push_back(im->first.c_str()) ;
          
          mx = mxCreateStructMatrix(1,1,nfield,&fieldnames.front());

          int ifield = 0  ;
          for ( GC::map_type::const_iterator im = mappa.begin() ; im != mappa.end() ; ++im, ++ifield ) {
            mxArray * mxi = nullptr ;
            GenericContainer_to_mxArray( im->second, mxi ) ;
            if ( mxi != nullptr ) mxSetFieldByNumber( mx, 0, ifield, mxi ) ;
          }
        }
      break;
    }
  }

}
