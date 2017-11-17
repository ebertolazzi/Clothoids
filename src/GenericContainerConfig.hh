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
// file: GenericContainerConfig.hh
//

#ifndef GENERIC_CONTAINER_CONFIG_HH
#define GENERIC_CONTAINER_CONFIG_HH

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  #define GENERIC_CONTAINER_ON_WINDOWS
  //#pragma comment(lib, "kernel32.lib")
  //#pragma comment(lib, "user32.lib")
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus >= 201103L)
  #ifndef GENERIC_CONTAINER_DO_NOT_USE_CXX11
    #define GENERIC_CONTAINER_USE_CXX11
  #endif
#else
  #include <cstdlib>
  #ifndef nullptr
    #include <cstddef>
    #ifndef nullptr
      #define nullptr NULL
    #endif
  #endif
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

#ifdef GENERIC_CONTAINER_USE_CXX11
  #ifndef GENERIC_CONTAINER_USE_REGEX
    #define GENERIC_CONTAINER_USE_REGEX
  #endif
#endif

// Standard types
#ifdef GENERIC_CONTAINER_ON_WINDOWS
  // se in windows includo PRIMA windows.h per evitare conflitti!
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #ifdef _MSC_VER
    #include <stdint.h>
  #else
    typedef          __int8  int8_t   ;
    typedef          __int16 int16_t  ;
    typedef          __int32 int32_t  ;
    typedef          __int64 int64_t  ;
    typedef unsigned __int8  uint8_t  ;
    typedef unsigned __int16 uint16_t ;
    typedef unsigned __int32 uint32_t ;
    typedef unsigned __int64 uint64_t ;
  #endif
#else
  #ifdef GENERIC_CONTAINER_USE_CXX11
    #include <cstdint>
  #else
    #include <stdint.h>
  #endif
#endif

#endif

//
// eof: GenericContainerConfig.hh
//
