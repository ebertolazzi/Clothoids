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
 |      Università degli Studi di Trento                                    |
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
#if defined( __APPLE__ ) && defined( __MACH__ )
// osx architecture
#define GENERIC_CONTAINER_ON_OSX 1
#elif defined( __unix__ )
// linux architecture
#define GENERIC_CONTAINER_ON_LINUX 1
#elif defined( _WIN32 ) || defined( WIN32 ) || defined( _WIN64 ) || defined( WIN64 )
// windows architecture
#define GENERIC_CONTAINER_ON_WINDOWS
#else
#warning "unsupported OS!"
#endif

// require a C++20 compiler (MSVC reports __cplusplus honestly only with
// /Zc:__cplusplus, so accept _MSVC_LANG as well)
#if ( defined( _MSVC_LANG ) && _MSVC_LANG >= 202002L ) || ( defined( __cplusplus ) && __cplusplus >= 202002L )
#else
#error "GenericContainer must be compiled as C++20 or later"
#endif

// Standard types
#include <cstdint>

#endif

//
// eof: GenericContainerConfig.hh
//
