#if defined(__llvm__) || defined(__clang__)
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#pragma clang diagnostic ignored "-Wduplicate-enum"
#endif

#ifndef _CRT_SECURE_NO_WARNINGS
  #define _CRT_SECURE_NO_WARNINGS
#endif

#include "Utils_fmt.hh"
#include "Utils/3rd/os.cxx"
