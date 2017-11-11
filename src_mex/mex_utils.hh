/****************************************************************************\
Copyright (c) 2015, Enrico Bertolazzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
\****************************************************************************/

#ifndef MEX_UTILS_HH
#define MEX_UTILS_HH

#include "mex.h"
#include <map>
#include <string>
#include <sstream>
#include <iostream>

#define arg_in_0 prhs[0]
#define arg_in_1 prhs[1]
#define arg_in_2 prhs[2]
#define arg_in_3 prhs[3]
#define arg_in_4 prhs[4]
#define arg_in_5 prhs[5]
#define arg_in_6 prhs[6]
#define arg_in_7 prhs[7]

#define arg_out_0 plhs[0]
#define arg_out_1 plhs[1]
#define arg_out_2 plhs[2]
#define arg_out_3 plhs[3]
#define arg_out_4 plhs[4]
#define arg_out_5 plhs[5]
#define arg_out_6 plhs[6]
#define arg_out_7 plhs[7]

#define MEX_ASSERT(COND,MSG)              \
  if ( !(COND) ) {                        \
    std::ostringstream ost ;              \
    ost << "Mex Error: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;     \
  }

// -----------------------------------------------------------------------------

static
inline
double
getScalarValue( mxArray const * arg, char const msg[] ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  MEX_ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  MEX_ASSERT( dims[0] == 1 && dims[1] == 1,
              msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  return mxGetScalar(arg) ;
}

static
inline
int64_t
getInt( mxArray const * arg, char const msg[] ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  MEX_ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  MEX_ASSERT( dims[0] == 1 && dims[1] == 1,
              msg << ", found " << dims[0] << " x " << dims[1] << " matrix" ) ;
  mxClassID category = mxGetClassID(arg);
  int64_t res = 0 ;
  void *ptr = mxGetData(arg) ;
  switch (category)  {
    case mxINT8_CLASS:   res = *static_cast<uint32_t*>(ptr); break;
    case mxUINT8_CLASS:  res = *static_cast<uint8_t*>(ptr);  break;
    case mxINT16_CLASS:  res = *static_cast<int16_t*>(ptr);  break;
    case mxUINT16_CLASS: res = *static_cast<uint16_t*>(ptr); break;
    case mxINT32_CLASS:  res = *static_cast<int32_t*>(ptr);  break;
    case mxUINT32_CLASS: res = *static_cast<uint32_t*>(ptr); break;
    case mxINT64_CLASS:  res = *static_cast<int64_t*>(ptr);  break;
    case mxUINT64_CLASS: res = *static_cast<uint64_t*>(ptr); break;
    default:
      MEX_ASSERT (false, msg << " bad type scalar" ) ;
    break;
  }
  return res ;
}

static
inline
double const *
getVectorPointer( mxArray const * arg, mwSize & sz, char const msg[] ) {
  mwSize number_of_dimensions = mxGetNumberOfDimensions(arg) ;
  MEX_ASSERT( number_of_dimensions == 2, msg ) ;
  mwSize const * dims = mxGetDimensions(arg) ;
  MEX_ASSERT( dims[0] == 1 || dims[1] == 1,
              msg << "\nExpect (1 x n or n x 1) matrix, found " <<
              dims[0] << " x " << dims[1] ) ;
  sz = dims[0]*dims[1] ;
  return mxGetPr(arg) ;
}

// -----------------------------------------------------------------------------

static
inline
void
setScalarValue( mxArray * & arg, double const & value ) {
  arg = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL) ;
  *mxGetPr(arg) = value;
}

static
inline
void
setScalarInt( mxArray * & arg, int64_t value ) {
  arg = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
  *static_cast<int64_t*>(mxGetData(arg)) = value ;
}

static
inline
int64_t *
createMatrixInt( mxArray * & arg, mwSize nrow, mwSize ncol ) {
  arg = mxCreateNumericMatrix( nrow, ncol, mxINT64_CLASS, mxREAL );
  return static_cast<int64_t*>(mxGetData(arg)) ;
}

static
inline
double *
createMatrixValue( mxArray * & arg, mwSize nrow, mwSize ncol ) {
  arg = mxCreateNumericMatrix( nrow, ncol, mxDOUBLE_CLASS, mxREAL );
  return static_cast<double*>(mxGetData(arg)) ;
}

// -----------------------------------------------------------------------------

/*
  Class Handle by Oliver Woodford

  https://it.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class
*/

#ifndef __MEX_CLASS_HANDLE_HH__
#define __MEX_CLASS_HANDLE_HH__
#include "mex.h"
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>

#define CLASS_HANDLE_SIGNATURE 0xFF00F0A5

template <typename base>
class class_handle {
  uint32_t    signature_m;
  base        *ptr_m;
  std::string name_m;
public:
  class_handle(base *ptr)
  : ptr_m(ptr)
  , name_m(typeid(base).name())
  { signature_m = CLASS_HANDLE_SIGNATURE; }

  ~class_handle()
  { signature_m = 0; delete ptr_m; }

  bool isValid()
  { return ((signature_m == CLASS_HANDLE_SIGNATURE) &&
            !strcmp(name_m.c_str(), typeid(base).name())); }

  base *ptr() { return ptr_m; }
};

template <typename base>
inline
mxArray *
convertPtr2Mat( base *ptr ) {
  mexLock();
  mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(new class_handle<base>(ptr));
  return out;
}

template <typename base>
inline
class_handle<base> *
convertMat2HandlePtr(const mxArray *in) {
  if ( mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
    mexErrMsgTxt("Input must be an uint64 scalar.");
  class_handle<base> *ptr = reinterpret_cast<class_handle<base> *>(*((uint64_t *)mxGetData(in)));
  if (!ptr->isValid())
    mexErrMsgTxt("Handle not valid.");
  return ptr;
}

template <typename base>
inline
base *
convertMat2Ptr(const mxArray *in) {
  return convertMat2HandlePtr<base>(in)->ptr();
}

template <typename base>
inline
void
destroyObject(const mxArray *in) {
  if ( in != nullptr ) delete convertMat2HandlePtr<base>(in);
  in = nullptr ;
  mexUnlock();
}

#endif // __CLASS_HANDLE_HPP__

#endif
