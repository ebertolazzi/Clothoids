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
