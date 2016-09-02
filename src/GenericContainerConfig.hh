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

// workaround for gcc < 4.9 that do not support regex
#if defined(__GNUC__) && !defined(__clang__)
  #if __GNUC__ > 5 || ( __GNUC__== 5 && __GNUC_MINOR__ > 8 )
  #else
    #ifndef DO_NOT_USE_CXX11
      #define DO_NOT_USE_CXX11
    #endif
  #endif
#endif

// comment to disable c++11 regex support
#define GENERIC_CONTAINER_USE_REGEX

#endif

//
// eof: GenericContainerConfig.hh
//
