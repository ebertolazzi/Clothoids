#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0   date 9/3/2022                                         #
#                                                                       #
#  Copyright (C) 2022                                                   #
#                                                                       #
#      Enrico Bertolazzi                                                #
#      Dipartimento di Ingegneria Industriale                           #
#      Universita` degli Studi di Trento                                #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#-----------------------------------------------------------------------#

cmake_minimum_required( VERSION 3.14 )

project( Clothoids )

include(./cmake_utils/CMakeLists-common.txt)

project( ${PROJECT_NAME} VERSION ${UTILS_PROJECT_VERSION} )

if ( NOT DEFINED CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()

set(Matlab_ROOT_DIR "/Applications/MATLAB_R2025b.app")
find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY MEX_LIBRARY MAIN_PROGRAM )

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )

if ( (UNIX OR LINUX) AND NOT APPLE )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined")
  add_compile_options(-fPIC)
endif()

include(./cmake_utils/CMakeLists-cflags.txt)
include(./cmake_utils/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang" )
  add_compile_options( "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-poison-system-directories>" )
  add_compile_options( "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-switch-default>" )
  add_compile_options( "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-exit-time-destructors>" )
  add_compile_options( "$<$<COMPILE_LANGUAGE:C,CXX>:-Wno-missing-include-dirs>" )
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ____                                _     _ _
#  |  _ \ _ __ ___ _ __ ___  __ _ _   _(_)___(_) |_ ___
#  | |_) | '__/ _ \ '__/ _ \/ _` | | | | / __| | __/ _ \
#  |  __/| | |  __/ | |  __/ (_| | |_| | \__ \ | ||  __/
#  |_|   |_|  \___|_|  \___|\__, |\__,_|_|___/_|\__\___|
#                              |_|
# MUST be executed before the build of the header list

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#

set( SOURCES )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src_lib/*.cc
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h*
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX "\.DS")

set( PRECOMPILED_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hh
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PRECOMPILED_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PRECOMPILED_HEADERS EXCLUDE REGEX "\.DS")

# extra include directory for compilation
include_directories( src src_lib src/Utils )

message( STATUS "SOURCES = ${SOURCES}" )

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wno-deprecated-declarations)
  add_compile_options(-Wno-deprecated)
  add_compile_options(-Wno-covered-switch-default)
  add_compile_options(-Wno-unused-result)
  add_compile_options(-Wno-unused-member-function)
  add_compile_options(-Wno-reserved-identifier)
endif()


#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#

utils_final_messages()

add_library( Clothoids STATIC ${SOURCES} )

set(
  MEX_LIST
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

message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if ( CMAKE_BUILD_TYPE MATCHES Debug )
  message( STATUS "DEBUG MODE\n\n")
  add_compile_definitions(G2LIB_DEBUG)
  add_compile_definitions(DEBUG)
endif()

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/bin")

foreach(MEX ${MEX_LIST})

  matlab_add_mex(
    NAME mex_${MEX}
    SRC src_mex/mex_${MEX}.cc
    LINK_TO Clothoids Matlab::mx Matlab::mex
    OUTPUT_NAME ${MEX}
  )

  # Rimuove flag di linking inutili e usa dynamic lookup per macOS
  set_target_properties(
    mex_${MEX} PROPERTIES 
    LINK_FLAGS "-Wl,-undefined,dynamic_lookup"
    BUILD_RPATH "${Matlab_ROOT_DIR}/extern/bin/maca64"
    INSTALL_RPATH "${Matlab_ROOT_DIR}/extern/bin/maca64"
    MACOSX_RPATH ON
  )

  add_custom_command(
    TARGET mex_${MEX}
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:mex_${MEX}>"
    "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    COMMENT "Copying ${MEX} to bin directory"
  )

endforeach()
