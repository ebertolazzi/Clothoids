/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2011                                                      |
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
#pragma once

#ifndef UTILS_CPU_INFO_dot_HH
#define UTILS_CPU_INFO_dot_HH

#include <string>

namespace Utils {

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  using std::string;
  void cpuId( unsigned long CPUInfo[4], unsigned long l );
  #endif

  /*!
   * Get info of Intel CPU
   *
   *
   * \param[out] bMMX      booelan, support MMX instruction set
   * \param[out] bMMXplus  booelan, support MMXplus instruction set
   * \param[out] bSSE      booelan, support SSE instruction set
   * \param[out] bSSE2     booelan, support SSE2 instruction set
   * \param[out] bSSE3     booelan, support SSE3 instruction set
   * \param[out] bSSSE3    booelan, support SSEE3 instruction set
   * \param[out] bSSE41    booelan, support SSE41 instruction set
   * \param[out] bSSE42    booelan, support SSE42 instruction set
   * \param[out] bSSE4a    booelan, support SSE4a instruction set
   * \param[out] bSSE5     booelan, support SSE5 instruction set
   * \param[out] b3Dnow    booelan, support 3Dnow instruction set
   * \param[out] b3DnowExt booelan, support 3DnowExt instruction set
   */
  void
  info(
    bool & bMMX,
    bool & bMMXplus,
    bool & bSSE,
    bool & bSSE2,
    bool & bSSE3,
    bool & bSSSE3,
    bool & bSSE41,
    bool & bSSE42,
    bool & bSSE4a,
    bool & bSSE5,
    bool & b3Dnow,
    bool & b3DnowExt
  );

  bool has_MMX();      //!< check if CPU support MMX instruction set
  bool has_MMXplus();  //!< check if CPU support MMXplus instruction set
  bool has_SSE();      //!< check if CPU support SSE instruction set
  bool has_SSE2();     //!< check if CPU support SSE2 instruction set
  bool has_SSE3();     //!< check if CPU support SSE3 instruction set
  bool has_SSSE3();    //!< check if CPU support SSSE3 instruction set
  bool has_SSE41();    //!< check if CPU support SSE41 instruction set
  bool has_SSE42();    //!< check if CPU support SSE42 instruction set
  bool has_SSE4a();    //!< check if CPU support SSE4a instruction set
  bool has_SSE5();     //!< check if CPU support SSE5 instruction set
  bool has_3Dnow();    //!< check if CPU support 3Dnow instruction set
  bool has_3DnowExt(); //!< check if CPU support 3DnowExt instruction set

  //!
  //! Return a string describing the CPU.
  //!
  string cpuInfo();
}

#endif
