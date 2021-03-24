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

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus >= 201103L)
#else
  #error "must use a compiler >= c++11"
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
    typedef          __int8  int8_t;
    typedef          __int16 int16_t;
    typedef          __int32 int32_t;
    typedef          __int64 int64_t;
    typedef unsigned __int8  uint8_t;
    typedef unsigned __int16 uint16_t;
    typedef unsigned __int32 uint32_t;
    typedef unsigned __int64 uint64_t;
  #endif
#else
  #include <cstdint>
#endif

#endif

//
// eof: GenericContainerConfig.hh
//
