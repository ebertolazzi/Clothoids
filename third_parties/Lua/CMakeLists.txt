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

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET( TARGET Lua )

PROJECT( ${TARGET} CXX C )

SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
SET( CMAKE_VERBOSE_MAKEFILE            TRUE )
SET( CMAKE_BUILD_TYPE                  Release )

IF( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC )
  SET( CMAKE_CXX_FLAGS "-Wall -Wcast-align" )
  IF ( APPLE )
    SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -O1 -g -DDEBUG " )
  ELSE()
    SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -fPIC -fno-common -O3 -g0 -funroll-loops " )
    SET( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -fPIC -fno-common -O1 -g -DDEBUG " )
  ENDIF()
ELSEIF( MSVC )
  SET( CMAKE_CXX_DEFS          "/D_WINDOWS /D_CRT_SECURE_NO_DEPRECATE /D_CRT_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC" )
  SET( CMAKE_CXX_COMMON        "/nologo /GS /W3 /WX- /Gm- /Gd /fp:precise /EHsc /FS" )
  SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_COMMON} /O2 /MD ${CMAKE_CXX_DEFS}" )
  SET( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_COMMON} /Od /Ob0 /MDd /Zi /RTC1 /D_DEBUG ${CMAKE_CXX_DEFS}" )
ENDIF()

SET( CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_CXX_FLAGS       ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
SET( CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG} )

# set 64 bit 
if (WIN32)
  set(WINXXBITS Win32)
endif(WIN32)
if (WIN64)
  set(WINXXBITS Win64)
endif(WIN64)

SET( SOURCES )
FILE( GLOB S ./src/*.c )
FOREACH (F ${S})
  IF ( ${F} MATCHES "wmain.c$" OR ${F} MATCHES "luac.c$"  ) 
  ELSE()
    FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
    LIST( APPEND SOURCES ${RF} )
  ENDIF()
ENDFOREACH (F ${S})

SET( HEADERS )
FILE( GLOB S ./src/*.h )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})
FILE( GLOB S ./src/*.hpp )
FOREACH (F ${S})
  FILE( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  LIST( APPEND HEADERS ${RF} )
ENDFOREACH (F ${S})

INCLUDE_DIRECTORIES( ./src )

ADD_LIBRARY( ${TARGET} ${SOURCES} )

INSTALL(TARGETS ${TARGET}
        RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib3rd
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib3rd
        ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib3rd )

INSTALL( FILES ${HEADERS} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib3rd/include/lua )

