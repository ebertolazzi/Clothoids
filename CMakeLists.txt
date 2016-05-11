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

cmake_minimum_required(VERSION 3.5)

SET( TARGET GenericContainer )

PROJECT( ${TARGET} CXX C )

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
SET( CMAKE_BUILD_TYPE                  Release )

IF( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC )
  SET( CMAKE_CXX_FLAGS "-Wall -Wcast-align" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-O1 -g -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "-fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "-fPIC -fno-common -O1 -g -DDEBUG " )
  ENDIF()
ELSEIF( MSVC )
ENDIF()
SET( CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS} )
SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./srcs/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS )
FILE( GLOB S ./srcs/*.hh )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

SET(HAVE_LUA false)
IF ( HAVE_LUA  )
  find_library(LUA_LIB NAMES lua lua53 HINTS "C:\\Mechatronix\\lib" )
  FILE( GLOB S ./srcs_lua_interface/*.cc )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    LIST( APPEND SOURCES ${RF} )
  ENDFOREACH (F ${S})
  FILE( GLOB S ./srcs_lua_interface/*.hh )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    LIST( APPEND HEADERS_WITH_LUA ${RF} )
  ENDFOREACH (F ${S})
  INCLUDE_DIRECTORIES( ./srcs ./srcs_lua_interface "C:\\MechatronixBuild\\include" )
ELSE()
  INCLUDE_DIRECTORIES( ./srcs )
ENDIF()

add_library( ${TARGET} ${SOURCES} )
ADD_EXECUTABLE( example1 examples/example1.cc ${HEADERS} )
ADD_EXECUTABLE( example2 examples/example2.cc ${HEADERS} )
ADD_EXECUTABLE( example3 examples/example3.cc ${HEADERS} )
ADD_EXECUTABLE( example4 examples/example4.cc ${HEADERS} )
ADD_EXECUTABLE( example5 examples/example5.cc ${HEADERS} )
ADD_EXECUTABLE( example6 examples/example6.c  ${HEADERS} ) # IT'S C!
ADD_EXECUTABLE( example7 examples/example7.cc ${HEADERS} )
IF ( HAVE_LUA )
  ADD_EXECUTABLE( example8  examples/example8.cc ${HEADERS} )
  ADD_EXECUTABLE( example9  examples/example9.cc ${HEADERS} )
  ADD_EXECUTABLE( example10 examples/example10.cc ${HEADERS} )
ENDIF()
ADD_EXECUTABLE( example10 examples/example11.cc ${HEADERS} )

INCLUDE_DIRECTORIES( srcs srcs_lua_interface )

TARGET_LINK_LIBRARIES( example1 ${TARGET} )
TARGET_LINK_LIBRARIES( example2 ${TARGET} )
TARGET_LINK_LIBRARIES( example3 ${TARGET} )
TARGET_LINK_LIBRARIES( example4 ${TARGET} )
TARGET_LINK_LIBRARIES( example5 ${TARGET} )
TARGET_LINK_LIBRARIES( example6 ${TARGET} )
IF ( HAVE_LUA )
  TARGET_LINK_LIBRARIES( example7 ${TARGET} ${LUA_LIB} )
  TARGET_LINK_LIBRARIES( example8 ${TARGET} ${LUA_LIB} )
  TARGET_LINK_LIBRARIES( example9 ${TARGET} ${LUA_LIB} )
ENDIF()
TARGET_LINK_LIBRARIES( example10 ${TARGET} )

MESSAGE( STATUS "Using ${SSE_FLAGS} extensions")
MESSAGE( STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}" )
MESSAGE( STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}" )
MESSAGE( STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}" )
MESSAGE( STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}" )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}" )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}" )
IF(APPLE)
  MESSAGE( STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}" )
  MESSAGE( STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}" )
  MESSAGE( STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}" )
ENDIF()
