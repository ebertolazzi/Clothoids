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
SET( CMAKE_VERBOSE_MAKEFILE OFF )
SET( CMAKE_INSTALL_MESSAGE NEVER )

IF ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  SET( DEBUG_SUFFIX "_debug" )
ELSE()
  SET( DEBUG_SUFFIX )
ENDIF()

IF ( BUILD_SHARED )
  IF ( APPLE )
    SET( TARGET GenericContainer_osx${DEBUG_SUFFIX} )
  ELSEIF( LINUX OR UNIX )
    SET( TARGET GenericContainer_linux${DEBUG_SUFFIX} )
  ELSE()
    SET( TARGET GenericContainer_win_${BITS}${DEBUG_SUFFIX} )
    SET( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )
  ENDIF()
ENDIF()

IF ( APPLE )
  SET( TARGETS GenericContainer_osx_static${DEBUG_SUFFIX} )
ELSEIF( LINUX OR UNIX )
  SET( TARGETS GenericContainer_linux_static${DEBUG_SUFFIX} )
ELSE()
  SET( TARGETS GenericContainer_win_${BITS}_static${DEBUG_SUFFIX} )
ENDIF()

PROJECT( GenericContainer CXX C )

INCLUDE( "CMakeLists-cflags.txt" )

MESSAGE( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )

SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

SET( SOURCES )
FILE( GLOB S ./src/*.cc )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND SOURCES ${RF} )
ENDFOREACH (F ${S})
SET( HEADERS )
FILE( GLOB S ./src/*.hh)
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./src/*.h)
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

INCLUDE( "./CMakeLists-lua.txt" )

IF ( LUA_FOUND )
  MESSAGE(STATUS "Adding Lua Support")
  FILE( GLOB S ./src_lua_interface/*.cc )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    LIST( APPEND SOURCES ${RF} )
  ENDFOREACH (F ${S})
  FILE( GLOB S ./src_lua_interface/*.hh )
  FOREACH (F ${S})
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    LIST( APPEND HEADERS ${RF} )
  ENDFOREACH (F ${S})
  INCLUDE_DIRECTORIES( src src_lua_interface ${LUA_INCLUDE_DIR} )
ELSE()
  MESSAGE(STATUS "Missing Lua Support")
  INCLUDE_DIRECTORIES( src )
ENDIF()

SET( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib )

IF ( BUILD_SHARED )
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  ADD_LIBRARY( ${TARGET}  SHARED ${SOURCES} ${HEADERS} )
  TARGET_LINK_LIBRARIES( ${TARGET} ${LUA_LIBRARIES} )
  INSTALL(
    TARGETS ${TARGET} ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ELSE()
  ADD_LIBRARY( ${TARGETS} STATIC ${SOURCES} ${HEADERS} )
  INSTALL(
    TARGETS ${TARGETS}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/dll
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
ENDIF()

SET_PROPERTY( TARGET ${TARGETS} PROPERTY POSITION_INDEPENDENT_CODE ON )

INSTALL( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )

IF( BUILD_EXECUTABLE )

  SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

  ADD_EXECUTABLE( example1 examples/example1.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example1 ${TARGETS} )

  ADD_EXECUTABLE( example2 examples/example2.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example2 ${TARGETS} )

  ADD_EXECUTABLE( example3 examples/example3.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example3 ${TARGETS} )

  ADD_EXECUTABLE( example4 examples/example4.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example4 ${TARGETS} )

  ADD_EXECUTABLE( example5 examples/example5.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example5 ${TARGETS} )

  ADD_EXECUTABLE( example6 examples/example6.c  ${HEADERS} ) # IT'S C!
  TARGET_LINK_LIBRARIES( example6 ${TARGETS} )

  ADD_EXECUTABLE( example7 examples/example7.cc ${HEADERS} )
  TARGET_LINK_LIBRARIES( example7 ${TARGETS} )

  IF ( LUA_FOUND )
    ADD_EXECUTABLE( example8 examples/example8.cc ${HEADERS} )
    TARGET_LINK_LIBRARIES( example8 ${TARGETS} ${LUA_LIBRARIES} )

    ADD_EXECUTABLE( example9 examples/example9.cc ${HEADERS} )
    TARGET_LINK_LIBRARIES( example9 ${TARGETS} ${LUA_LIBRARIES} )

    ADD_EXECUTABLE( example10 examples/example10.cc ${HEADERS} )
    TARGET_LINK_LIBRARIES( example10 ${TARGETS} ${LUA_LIBRARIES} )

    ADD_EXECUTABLE( example11 examples/example11.cc ${HEADERS} )
    TARGET_LINK_LIBRARIES( example11 ${TARGETS} ${LUA_LIBRARIES} )
  ENDIF()

ENDIF()

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
MESSAGE( STATUS "LUA_LIBRARIES               = ${LUA_LIBRARIES}" )
