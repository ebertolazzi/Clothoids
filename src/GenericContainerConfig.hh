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

#endif

//
// eof: GenericContainerConfig.hh
//
