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

#pragma once

#ifndef GENERIC_CONTAINER_CONFIG_HH
#define GENERIC_CONTAINER_CONFIG_HH

// select computer architecture
#if defined(__APPLE__) && defined(__MACH__)
  // osx architecture
  #define GENERIC_CONTAINER_ON_OSX 1
#elif defined(__unix__)
  // linux architecture
  #define GENERIC_CONTAINER_ON_LINUX 1
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define GENERIC_CONTAINER_ON_WINDOWS
#else
  #warning "unsupported OS!"
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
#else
  #error "Lapack Wrapper must be compiled using C++ >= C++11"
#endif


// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus >= 201103L)
#else
  #error "must use a compiler >= c++11"
#endif

// Standard types
#ifdef GENERIC_CONTAINER_ON_WINDOWS
  #ifdef GENERIC_CONTAINER_USE_WINDOWS_TYPES
    typedef          __int8  int8_t;
    typedef          __int16 int16_t;
    typedef          __int32 int32_t;
    typedef          __int64 int64_t;
    typedef unsigned __int8  uint8_t;
    typedef unsigned __int16 uint16_t;
    typedef unsigned __int32 uint32_t;
    typedef unsigned __int64 uint64_t;
  #else
    #include <cstdint>
  #endif
#else
  #include <cstdint>
#endif

#endif

//
// eof: GenericContainerConfig.hh
//
