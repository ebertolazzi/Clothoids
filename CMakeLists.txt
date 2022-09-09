############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

project( Clothoids )
include(./cmake/CMakeLists-common.txt)

project(
  ${PROJECT_NAME}
  VERSION ${EB_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
)

include(./cmake/CMakeLists-cflags.txt)
include(./cmake/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( EB_BUILD_SHARED )
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

set( SOURCES )
file( GLOB S ${SOURCE_DIR}/*.cc )
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB S
  ${SOURCE_DIR}/*.h*
  ${SOURCE_DIR}/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

# extra include directory for compilation
include_directories( ./lib3rd/include )

#       _                           _                 _
#    __| | ___ _ __   ___ _ __   __| | ___ _ __   ___(_) ___  ___
#   / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| |/ _ \/ __|
#  | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |  __/\__ \
#   \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|_|\___||___/
#             |_|
#
set(DEPEND_TARGETS)

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
set(UPDIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
find_library(
  UTILS Utils${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${UPDIR}/Utils/lib/lib
  NO_DEFAULT_PATH
)
if( UTILS )
  message(STATUS "Found Utils${ARTIFACTS_STATIC_SUFFIX}")
  ebertolazzi_copy_directory(
    "${UPDIR}/Utils/lib"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd"
  )
else()
  message(STATUS "NOT Found Utils${ARTIFACTS_STATIC_SUFFIX} use submodule")
  add_subdirectory(./submodules/Utils)
  include_directories(./submodules/Utils/src)
  set( UTILS ebertolazzi_Utils_Static )
  set( DEPEND_TARGETS ${DEPEND_TARGETS} ${UTILS} )
endif()

# DEPENDENCIES (TRY TO USE ALREADY COMPILED SUBMODULE IN LARGE PROJECT)
find_library(
  ROOTS quarticRootsFlocke${ARTIFACTS_STATIC_SUFFIX}
  HINTS ${UPDIR}/quarticRootsFlocke/lib/lib
  NO_DEFAULT_PATH
)
if( ROOTS )
  message(STATUS "Found quarticRootsFlocke${ARTIFACTS_STATIC_SUFFIX}")
  include_directories(${UPDIR}/quarticRootsFlocke/lib/include)
else()
  message(STATUS "NOT Found quarticRootsFlocke${ARTIFACTS_STATIC_SUFFIX}")
  add_subdirectory(./submodules/quarticRootsFlocke)
  include_directories(./submodules/quarticRootsFlocke/src)
  set( ROOTS ebertolazzi_quarticRootsFlocke_Static )
  set( DEPEND_TARGETS ${DEPEND_TARGETS} ${ROOTS} )
endif()

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
ebertolazzi_setup_target(
  ${PROJECT_NAME}
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
if ( EB_BUILD_SHARED )
  ebertolazzi_setup_target(
    ${PROJECT_NAME}
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
  target_link_libraries( ${EB_NAMESPACE}_${PROJECT_NAME} ${ROOTS} ${UTILS} )
endif()

if ( DEPEND_TARGETS )
  add_dependencies( ${EB_NAMESPACE}_${PROJECT_NAME}_Static ${DEPEND_TARGETS} )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( EB_ENABLE_TESTS )

  enable_testing()

  if(APPLE OR UNIX OR MINGW )
    find_package( ZLIB REQUIRED )
  else()
    if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      find_library( ZLIB_LIBRARIES zlibstaticd HINTS ./lib3rd/lib )
    else()
      find_library( ZLIB_LIBRARIES zlibstatic HINTS ./lib3rd/lib )
    endif()
    set( ZLIB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/zlib/include ")
  endif()

  set( LIBS ${EB_NAMESPACE}_${PROJECT_NAME}_Static ${UTILS} ${ROOTS} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )

  if( UNIX )
    if ( NOT APPLE )
      set( THREADS_PREFER_PTHREAD_FLAG ON )
      find_package( Threads REQUIRED )
      set( LIBS ${LIBS} Threads::Threads )
    endif()
  else()
    #libcmt.lib
    #find_library( WS2 ws2_32 REQUIRED )
    #find_library( IPH iphlpapi REQUIRED )
    #find_library( KER kernel32 REQUIRED )
    #find_library( IPP ippsmt REQUIRED )
    #find_library( WPTH pthread REQUIRED )
    #set( LIBS ${LIBS} ${WS2} ${IPH} ${KER} )
    set( LIBS ${LIBS} ws2_32 iphlpapi kernel32 )
    if (EB_BUILD_SHARED)
      target_link_libraries( ${EB_NAMESPACE}_${PROJECT_NAME} ${LIBS} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )
    endif()
  endif()

  set(
    EXELISTCPP
    testBiarc
    testDistance
    testDistanceCircuit
    testG2
    testG2plot
    testG2stat
    testG2stat2arc
    testG2statCLC
    testIntersect
    testIntersect2
    testIntersect3
    testPolyline
    testTriangle2D
  )

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    add_test(
      NAME "${S}"
      COMMAND ./bin/${S}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()

endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

include( ./cmake/CMakeLists-cpack.txt )

ebertolazzi_final_messages()
