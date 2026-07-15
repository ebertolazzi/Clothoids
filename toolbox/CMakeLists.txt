#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  Copyright (C) 2022                                                   #
#                                                                       #
#      Enrico Bertolazzi                                                #
#      Dipartimento di Ingegneria Industriale                           #
#      Universita` degli Studi di Trento                                #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#-----------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.25)

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# Determina l'architettura nativa del sistema host - con fallback manuale
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  
  # Prima prova a rilevare l'architettura dal sistema
  if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    set(HOST_ARCH "arm64")
  elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
    # Potrebbe essere Rosetta - verifichiamo se il sistema supporta ARM
    execute_process(
      COMMAND uname -m
      OUTPUT_VARIABLE UNAME_M
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(UNAME_M MATCHES "arm64")
      set(HOST_ARCH "arm64")
      message(STATUS "Detected ARM64 via uname -m (running under Rosetta)")
    else()
      set(HOST_ARCH "x86_64")
    endif()
  else()
    set(HOST_ARCH "x86_64")  # fallback
  endif()
  
  # Verifica aggiuntiva: se esiste /usr/bin/arch, usalo per testare
  if(EXISTS /usr/bin/arch)
    execute_process(
      COMMAND /usr/bin/arch
      OUTPUT_VARIABLE ARCH_OUT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(ARCH_OUT MATCHES "arm64")
      set(HOST_ARCH "arm64")
      message(STATUS "Detected ARM64 via arch command")
    endif()
  endif()
  
  # Di default, compiliamo SOLO per l'architettura nativa
  set(UTILS_BUILD_FAT_DEFAULT OFF)
  
  # Se l'utente vuole esplicitamente FAT build, usiamo entrambe
  # altrimenti solo l'architettura nativa
  option(UTILS_BUILD_FAT "Build Apple universal/fat binaries" OFF)
  
  # Architetture da includere - di default solo quella nativa
  if(UTILS_BUILD_FAT)
    set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")
  else()
    set(UTILS_FAT_ARCHITECTURES "${HOST_ARCH}" CACHE STRING "Apple architectures to build" FORCE)
  endif()
  
  message(STATUS "Host architecture: ${HOST_ARCH}")
  message(STATUS "Target architectures: ${UTILS_FAT_ARCHITECTURES}")
  
  # Forza l'architettura in CMAKE_OSX_ARCHITECTURES
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
  
  message(STATUS "CMAKE_OSX_ARCHITECTURES set to: ${CMAKE_OSX_ARCHITECTURES}")
  
else()
  set(UTILS_BUILD_FAT_DEFAULT OFF)
  option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
  set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")
endif()

# ----------------------------------------------------------------------------
# Project version discovery
# ----------------------------------------------------------------------------
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(GIT_VERSION MATCHES "^v")
    string(SUBSTRING ${GIT_VERSION} 1 -1 GIT_VERSION)
  endif()
endif()
if(NOT GIT_VERSION)
  set(GIT_VERSION "0.0.0")
endif()

project(
  ClothoidsToolbox
  VERSION ${GIT_VERSION}
  DESCRIPTION "MATLAB MEX wrappers for Clothoids library"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
  LANGUAGES CXX C
)

# ----------------------------------------------------------------------------
# Build options
# ----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(CLOTHOIDS_TOOLBOX_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_TOOLBOX_STRICT_WARNINGS "Compile with narrowing/conversion warnings" OFF)

# ----------------------------------------------------------------------------
# Install prefix and directory layout
# ----------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)
set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "Install binary directory" FORCE)

set(CLOTHOIDS_LOCAL_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CLOTHOIDS_LOCAL_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib")
set(CLOTHOIDS_LOCAL_3RD_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/lib")

include(GNUInstallDirs)

# ----------------------------------------------------------------------------
# Find MATLAB - with multiple fallback strategies
# ----------------------------------------------------------------------------
# MEX wrappers are optional: if MATLAB can't be found, the Clothoids C++
# library still builds fine, only the MEX targets are skipped.
option(CLOTHOIDS_TOOLBOX_BUILD_MEX "Build MATLAB MEX wrappers (requires MATLAB)" ON)

if(CLOTHOIDS_TOOLBOX_BUILD_MEX)

  # A directory is only accepted as a MATLAB root if it really looks like one.
  # This prevents unrelated env vars (MATLAB_ROOT / MATLAB_HOME set to something
  # else entirely, e.g. a TeX distribution) from silently hijacking detection.
  function(_clothoids_is_valid_matlab_root _dir _out_var)
    if(_dir AND EXISTS "${_dir}/extern/include/matrix.h")
      set(${_out_var} TRUE PARENT_SCOPE)
    else()
      set(${_out_var} FALSE PARENT_SCOPE)
    endif()
  endfunction()

  # Candidate MATLAB "arch" subdirectories (bin/<arch>, extern/bin/<arch>) to try,
  # in priority order. On Apple Silicon, MATLAB R2023b+ ships native arm64
  # binaries under "maca64"; older/Intel installs (or MATLAB run under Rosetta)
  # use "maci64" - so both are tried, native first.
  if(APPLE)
    # Determina l'architettura del sistema per scegliere l'ordine delle architetture MATLAB
    # Usa HOST_ARCH per determinare l'architettura nativa
    if(HOST_ARCH MATCHES "arm64")
      set(_matlab_arch_dirs "maca64" "maci64")
    else()
      set(_matlab_arch_dirs "maci64" "maca64")
    endif()
  elseif(WIN32)
    set(_matlab_arch_dirs "win64")
  else()
    set(_matlab_arch_dirs "glnxa64")
  endif()

  # 1) An explicit -DMatlab_ROOT_DIR=... or MATLAB_ROOT / MATLAB_HOME env var
  #    wins, but only if it actually points at a MATLAB installation.
  set(_matlab_root_candidate "")
  if(DEFINED Matlab_ROOT_DIR)
    set(_matlab_root_candidate "${Matlab_ROOT_DIR}")
  elseif(DEFINED ENV{MATLAB_ROOT})
    set(_matlab_root_candidate "$ENV{MATLAB_ROOT}")
  elseif(DEFINED ENV{MATLAB_HOME})
    set(_matlab_root_candidate "$ENV{MATLAB_HOME}")
  endif()

  _clothoids_is_valid_matlab_root("${_matlab_root_candidate}" _matlab_root_ok)
  if(_matlab_root_candidate AND NOT _matlab_root_ok)
    message(WARNING "Ignoring '${_matlab_root_candidate}' as MATLAB root "
                     "(no extern/include/matrix.h found there) - searching common locations instead.")
    set(_matlab_root_candidate "")
    unset(Matlab_ROOT_DIR CACHE)
  endif()

  # 2) Otherwise scan common install locations. Globbing (instead of a fixed
  #    per-release list) means newer MATLAB releases are found automatically.
  if(NOT _matlab_root_candidate)
    if(APPLE)
      file(GLOB _matlab_glob_candidates "/Applications/MATLAB_R*.app" "/usr/local/MATLAB/R*")
    elseif(WIN32)
      file(GLOB _matlab_glob_candidates "C:/Program Files/MATLAB/R*")
    else()
      file(GLOB _matlab_glob_candidates "/usr/local/MATLAB/R*" "/opt/MATLAB/R*")
    endif()
    if(_matlab_glob_candidates)
      list(SORT _matlab_glob_candidates COMPARE NATURAL ORDER DESCENDING)
    endif()
    foreach(_cand ${_matlab_glob_candidates})
      _clothoids_is_valid_matlab_root("${_cand}" _cand_ok)
      if(_cand_ok)
        set(_matlab_root_candidate "${_cand}")
        break()
      endif()
    endforeach()
  endif()

  if(_matlab_root_candidate)
    set(Matlab_ROOT_DIR "${_matlab_root_candidate}" CACHE PATH "MATLAB root directory" FORCE)
  endif()

  # Try the standard FindMatlab module first
  find_package(Matlab QUIET)

  # Manual fallback: locate mx/mex directly, trying each possible arch subdir
  if(NOT Matlab_FOUND AND Matlab_ROOT_DIR)
    message(STATUS "MATLAB not found via standard FindMatlab, trying manual detection...")
    message(STATUS "Using Matlab_ROOT_DIR: ${Matlab_ROOT_DIR}")

    set(Matlab_INCLUDE_DIR "${Matlab_ROOT_DIR}/extern/include")

    foreach(_arch_dir ${_matlab_arch_dirs})
      set(_lib_dir "${Matlab_ROOT_DIR}/bin/${_arch_dir}")
      if(NOT EXISTS "${_lib_dir}")
        continue()
      endif()
      unset(MX_LIBRARY CACHE)
      unset(MEX_LIBRARY CACHE)
      find_library(MX_LIBRARY  NAMES mx  libmx  HINTS "${_lib_dir}" NO_DEFAULT_PATH)
      find_library(MEX_LIBRARY NAMES mex libmex HINTS "${_lib_dir}" NO_DEFAULT_PATH)
      if(MX_LIBRARY AND MEX_LIBRARY)
        set(Matlab_LIBRARY_DIR "${_lib_dir}")
        set(Matlab_MEX_ARCH "${_arch_dir}")
        break()
      endif()
    endforeach()

    find_program(Matlab_MAIN_PROGRAM NAMES matlab HINTS "${Matlab_ROOT_DIR}/bin" NO_DEFAULT_PATH)

    if(MX_LIBRARY AND MEX_LIBRARY)
      set(Matlab_FOUND TRUE)
      set(Matlab_INCLUDE_DIRS "${Matlab_INCLUDE_DIR}")
      set(Matlab_LIBRARIES "${MX_LIBRARY};${MEX_LIBRARY}")
      message(STATUS "MATLAB manually found:")
      message(STATUS "  arch:        ${Matlab_MEX_ARCH}")
      message(STATUS "  MX_LIBRARY:  ${MX_LIBRARY}")
      message(STATUS "  MEX_LIBRARY: ${MEX_LIBRARY}")
      message(STATUS "  INCLUDE:     ${Matlab_INCLUDE_DIR}")
    else()
      message(WARNING "Could not find MATLAB mx/mex libraries under any of: "
                       "${Matlab_ROOT_DIR}/bin/{${_matlab_arch_dirs}}")
    endif()
  elseif(Matlab_FOUND)
    # find_package(Matlab) succeeded on its own - figure out which arch
    # subdir it used, needed later for RPATH on the MEX targets.
    list(GET _matlab_arch_dirs 0 Matlab_MEX_ARCH)
    foreach(_arch_dir ${_matlab_arch_dirs})
      if(Matlab_LIBRARIES MATCHES "${_arch_dir}")
        set(Matlab_MEX_ARCH "${_arch_dir}")
        break()
      endif()
    endforeach()
  endif()

  if(NOT Matlab_FOUND)
    message(WARNING "
    ==================================================================
    MATLAB not found - MEX wrappers will NOT be built.

    Please specify the MATLAB root directory with:
      cmake .. -DMatlab_ROOT_DIR=/path/to/MATLAB
    Common paths:
      macOS:   /Applications/MATLAB_R2025b.app
      Windows: C:/Program Files/MATLAB/R2025b
      Linux:   /usr/local/MATLAB/R2025b
    Or set an environment variable:
      export MATLAB_ROOT=/path/to/MATLAB
      export MATLAB_HOME=/path/to/MATLAB

    To silence this warning and only build the Clothoids C++ library, re-run with:
      cmake .. -DCLOTHOIDS_TOOLBOX_BUILD_MEX=OFF
    ==================================================================
    ")
    set(CLOTHOIDS_TOOLBOX_BUILD_MEX OFF)
  endif()

endif() # CLOTHOIDS_TOOLBOX_BUILD_MEX

# ----------------------------------------------------------------------------
# Compiler and platform settings
# ----------------------------------------------------------------------------
if(UNIX AND NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
  add_compile_options(-fPIC)
endif()

# Clang-specific warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-poison-system-directories>"
    "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-switch-default>"
    "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-exit-time-destructors>"
    "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-missing-include-dirs>"
  )
  if(CLOTHOIDS_TOOLBOX_STRICT_WARNINGS)
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:C,CXX>:-Wconversion>"
      "$<$<COMPILE_LANGUAGE:C,CXX>:-Wsign-conversion>"
      "$<$<COMPILE_LANGUAGE:C,CXX>:-Wshadow>"
    )
  endif()
  add_compile_options(-Wno-deprecated-declarations)
  add_compile_options(-Wno-deprecated)
  add_compile_options(-Wno-covered-switch-default)
  add_compile_options(-Wno-unused-result)
  add_compile_options(-Wno-unused-member-function)
  add_compile_options(-Wno-reserved-identifier)
endif()

message(STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")

# ----------------------------------------------------------------------------
# Debug mode definitions
# ----------------------------------------------------------------------------
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "DEBUG MODE")
  add_compile_definitions(G2LIB_DEBUG)
  add_compile_definitions(DEBUG)
endif()

# ----------------------------------------------------------------------------
# Clothoids and dependencies - already in src/
# ----------------------------------------------------------------------------
set(CLOTHOIDS_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CLOTHOIDS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# ----------------------------------------------------------------------------
# Source discovery - Clothoids library
# ----------------------------------------------------------------------------
set(CLOTHOIDS_SOURCES)
file(GLOB_RECURSE S
  "${CLOTHOIDS_SRC_DIR}/*.cc"
)
foreach(F ${S})
  file(RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}")
  list(APPEND CLOTHOIDS_SOURCES ${RF})
endforeach()

# Exclude src_mex files from Clothoids library
list(FILTER CLOTHOIDS_SOURCES EXCLUDE REGEX "src_mex/")

# ----------------------------------------------------------------------------
# Public headers discovery - Clothoids
# ----------------------------------------------------------------------------
set(CLOTHOIDS_PUBLIC_HEADERS)
file(GLOB_RECURSE S
  "${CLOTHOIDS_SRC_DIR}/*.h"
  "${CLOTHOIDS_SRC_DIR}/*.hh"
  "${CLOTHOIDS_SRC_DIR}/*.hpp"
  "${CLOTHOIDS_SRC_DIR}/*.hxx"
)
foreach(F ${S})
  file(RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}")
  list(APPEND CLOTHOIDS_PUBLIC_HEADERS ${RF})
endforeach()

list(FILTER CLOTHOIDS_PUBLIC_HEADERS EXCLUDE REGEX "\\.DS")
list(FILTER CLOTHOIDS_PUBLIC_HEADERS EXCLUDE REGEX "src_mex/")

message(STATUS "Clothoids sources: ${CLOTHOIDS_SOURCES}")

# ----------------------------------------------------------------------------
# Build Clothoids library
# ----------------------------------------------------------------------------
add_library(Clothoids STATIC ${CLOTHOIDS_SOURCES})

target_include_directories(Clothoids
  PUBLIC
    "$<BUILD_INTERFACE:${CLOTHOIDS_SRC_DIR}>"
    "$<BUILD_INTERFACE:${CLOTHOIDS_SRC_DIR}/src>"
)

# Clothoids compile features and definitions
target_compile_features(Clothoids PUBLIC cxx_std_20)
target_compile_definitions(Clothoids PUBLIC FMT_HEADER_ONLY)

# Link flags for Clothoids
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
  target_compile_options(Clothoids PUBLIC -stdlib=libc++)
  target_link_options(Clothoids PUBLIC -stdlib=libc++)
endif()

# Compiler warnings for Clothoids
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
  target_compile_options(Clothoids PRIVATE -Wall -Wextra)
  if(CLOTHOIDS_TOOLBOX_STRICT_WARNINGS)
    target_compile_options(
      Clothoids PRIVATE
      -Wconversion -Wsign-conversion -Wshadow -Wdouble-promotion
    )
  endif()
endif()

if(MSVC)
  target_compile_options(Clothoids PUBLIC /utf-8 /bigobj)
endif()

# Output directory for Clothoids library
set_target_properties(Clothoids PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CLOTHOIDS_LOCAL_LIBRARY_DIR}"
  LIBRARY_OUTPUT_DIRECTORY "${CLOTHOIDS_LOCAL_LIBRARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${CLOTHOIDS_LOCAL_BIN_DIR}"
)

# Create alias for consistent linking
add_library(Clothoids::Clothoids ALIAS Clothoids)

# ----------------------------------------------------------------------------
# MEX file definitions
# ----------------------------------------------------------------------------
set(MEX_LIST
  ClothoidListMexWrapper
  FresnelCS
  XY_to_angle
  CircleArcMexWrapper
  BiarcMexWrapper
  BiarcListMexWrapper
  LineSegmentMexWrapper
  PolyLineMexWrapper
  Triangle2DMexWrapper
  ClothoidCurveMexWrapper
  ClothoidSplineG2MexWrapper
  TriTriOverlap
  DubinsMexWrapper
  Dubins3pMexWrapper
)

if(CLOTHOIDS_TOOLBOX_BUILD_MEX)

  # Create bin directory
  add_custom_command(
    TARGET Clothoids PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CLOTHOIDS_LOCAL_BIN_DIR}"
  )

  # --------------------------------------------------------------------------
  # MEX file generation
  # --------------------------------------------------------------------------
  foreach(MEX ${MEX_LIST})
    # Check if source file exists
    set(MEX_SRC "src_mex/mex_${MEX}.cc")
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${MEX_SRC}")
      message(WARNING "MEX source not found: ${MEX_SRC}, skipping...")
      continue()
    endif()

    matlab_add_mex(
      NAME mex_${MEX}
      SRC ${MEX_SRC}
      LINK_TO Clothoids::Clothoids
      OUTPUT_NAME ${MEX}
    )

    # Include directories for MEX
    target_include_directories(mex_${MEX}
      PRIVATE
        "${CLOTHOIDS_SRC_DIR}"
        "${CLOTHOIDS_SRC_DIR}/include"
        "${CLOTHOIDS_SRC_DIR}/src"
        "${CLOTHOIDS_SRC_DIR}/src/Utils"
        "${Matlab_INCLUDE_DIRS}"
    )

    # Platform-specific link flags. Matlab_MEX_ARCH is resolved once during
    # MATLAB detection (maca64/maci64 on Apple Silicon/Intel, glnxa64, win64).
    if(APPLE)
      set_target_properties(mex_${MEX} PROPERTIES
        LINK_FLAGS "-Wl,-undefined,dynamic_lookup"
        BUILD_RPATH "${Matlab_ROOT_DIR}/extern/bin/${Matlab_MEX_ARCH}"
        INSTALL_RPATH "${Matlab_ROOT_DIR}/extern/bin/${Matlab_MEX_ARCH}"
        MACOSX_RPATH ON
      )
      if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
        target_compile_options(mex_${MEX} PRIVATE -stdlib=libc++)
        target_link_options(mex_${MEX} PRIVATE -stdlib=libc++)
      endif()
    elseif(UNIX)
      # Linux
      set_target_properties(mex_${MEX} PROPERTIES
        BUILD_RPATH "${Matlab_ROOT_DIR}/bin/${Matlab_MEX_ARCH}"
        INSTALL_RPATH "${Matlab_ROOT_DIR}/bin/${Matlab_MEX_ARCH}"
      )
    elseif(WIN32)
      # Windows
      set_target_properties(mex_${MEX} PROPERTIES
        BUILD_RPATH "${Matlab_ROOT_DIR}/bin/${Matlab_MEX_ARCH}"
        INSTALL_RPATH "${Matlab_ROOT_DIR}/bin/${Matlab_MEX_ARCH}"
      )
    endif()

    # Copy MEX to bin directory
    add_custom_command(
      TARGET mex_${MEX}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:mex_${MEX}>"
        "${CLOTHOIDS_LOCAL_BIN_DIR}"
      COMMENT "Copying ${MEX} to bin directory"
    )
  endforeach()

else()
  message(STATUS "CLOTHOIDS_TOOLBOX_BUILD_MEX is OFF - skipping MEX wrapper targets")
endif()

# ----------------------------------------------------------------------------
# Configuration summary
# ----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "============================================================")
message(STATUS "ClothoidsToolbox configuration:")
message(STATUS "============================================================")
message(STATUS "  MATLAB root:     ${Matlab_ROOT_DIR}")
message(STATUS "  MATLAB include:  ${Matlab_INCLUDE_DIRS}")
message(STATUS "  Clothoids src:   ${CLOTHOIDS_SRC_DIR}")
message(STATUS "  Bin dir:         ${CLOTHOIDS_LOCAL_BIN_DIR}")
message(STATUS "  Lib dir:         ${CLOTHOIDS_LOCAL_LIBRARY_DIR}")
message(STATUS "  Build type:      ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler:        ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  MEX files:       ${MEX_LIST}")
message(STATUS "============================================================")
message(STATUS "")

# ----------------------------------------------------------------------------
# Install rules
# ----------------------------------------------------------------------------
if(CLOTHOIDS_TOOLBOX_INSTALL)
  install(
    TARGETS Clothoids
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  # Install MEX files
  install(
    DIRECTORY "${CLOTHOIDS_LOCAL_BIN_DIR}/"
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    FILES_MATCHING
      PATTERN "*.mex*"
      PATTERN "*.dll"
      PATTERN "*.so"
      PATTERN "*.dylib"
  )

  # Install public headers
  foreach(_header IN LISTS CLOTHOIDS_PUBLIC_HEADERS)
    get_filename_component(_dir "${_header}" DIRECTORY)
    install(
      FILES "${CMAKE_CURRENT_SOURCE_DIR}/${_header}"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${_dir}"
    )
  endforeach()

  # CPack configuration
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
  include(CPack)
endif()
