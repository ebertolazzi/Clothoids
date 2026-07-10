############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Università degli Studi di Trento                                    #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.25)

# ----------------------------------------------------------------------------
# Estrai versione da git tag (se disponibile)
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
  Clothoids
  VERSION ${GIT_VERSION}
  DESCRIPTION "A C++ library for clothoid splines"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
  LANGUAGES CXX C
)

# ----------------------------------------------------------------------------
# Opzioni principali
# ----------------------------------------------------------------------------
# Default a Release se non specificato (per generatori a singola configurazione)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(BUILD_SHARED_LIBS "Build Clothoids as a shared library" OFF)
option(BUILD_TESTING "Build the Clothoids test executables" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_BUILD_BENCHMARKS "Build the micro-benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_STRICT_WARNINGS "Compile with narrowing/conversion warnings" OFF)

# Directory di installazione predefinita: ./lib
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)
set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "Install binary directory" FORCE)

# Local install layout. Keep generated libraries only in lib/lib; the
# project-root lib/ directory must not receive duplicate archives.
set(CLOTHOIDS_LOCAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/include")
set(CLOTHOIDS_LOCAL_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib")
set(CLOTHOIDS_LOCAL_RUNTIME_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/bin")
set(CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include")

include(GNUInstallDirs)
include(FetchContent)

# ----------------------------------------------------------------------------
# Helper per dipendenze: locale o FetchContent
# ----------------------------------------------------------------------------
function(clothoids_use_local_or_fetch NAME LOCAL_DIR)
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  string(TOUPPER "${NAME}" _name_upper)
  string(TOLOWER "${NAME}" _name_lower)

  if(FETCHCONTENT_BASE_DIR)
    set(_fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}")
  else()
    set(_fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps")
  endif()

  set(_dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake")
  set(_saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}")
  set(CMAKE_SKIP_INSTALL_RULES ON)

  if(EXISTS "${_local_abs}/CMakeLists.txt")
    message(STATUS "Using local ${NAME}: ${_local_abs}")
    set(FETCHCONTENT_SOURCE_DIR_${_name_upper} "${_local_abs}" CACHE PATH "Local ${NAME} checkout" FORCE)
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "local: ${_local_abs}" PARENT_SCOPE)
    set("${NAME}_SOURCE_DIR" "${_local_abs}" PARENT_SCOPE)
  else()
    message(STATUS "Fetching ${NAME}")
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "FetchContent" PARENT_SCOPE)
    if(DEFINED "${_name_lower}_SOURCE_DIR")
      set("${NAME}_SOURCE_DIR" "${${_name_lower}_SOURCE_DIR}" PARENT_SCOPE)
    endif()
  endif()

  file(MAKE_DIRECTORY "${_fetchcontent_base_dir}/${_name_lower}-build")
  if(NOT EXISTS "${_dependency_install_stub}")
    file(WRITE "${_dependency_install_stub}" "# Install disabled for dependency ${NAME} inside Clothoids.\n")
  endif()

  set(CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}")
endfunction()

# ----------------------------------------------------------------------------
# Percorsi locali per checkout dei dipendenti (override via cache)
# ----------------------------------------------------------------------------
set(CLOTHOIDS_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout")
set(CLOTHOIDS_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout")
set(CLOTHOIDS_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout")

# ----------------------------------------------------------------------------
# 1. GenericContainer
# ----------------------------------------------------------------------------
set(GENERIC_CONTAINER_INSTALL OFF CACHE BOOL "Disable GenericContainer install rules inside Clothoids" FORCE)
set(GENERIC_CONTAINER_BUILD_EXAMPLES OFF CACHE BOOL "Disable GenericContainer examples inside Clothoids" FORCE)

FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/ebertolazzi/GenericContainer.git
  GIT_TAG        GC2.0
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(GenericContainer "${CLOTHOIDS_GC_LOCAL_DIR}")

# Patch per Clang: forzare -stdlib=libc++ se necessario (come in Splines)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  foreach(_gc_target GenericContainer GenericContainerYaml)
    if(TARGET ${_gc_target})
      get_target_property(_gc_type ${_gc_target} TYPE)
      if(_gc_type STREQUAL "INTERFACE_LIBRARY")
        target_compile_options(${_gc_target} INTERFACE -stdlib=libc++)
        target_link_options(${_gc_target} INTERFACE -stdlib=libc++)
      else()
        target_compile_options(${_gc_target} PUBLIC -stdlib=libc++)
        target_link_options(${_gc_target} PUBLIC -stdlib=libc++)
      endif()
    endif()
  endforeach()
endif()

# GenericContainer install/examples are disabled before FetchContent_MakeAvailable.

# ----------------------------------------------------------------------------
# 2. UtilsLite
# ----------------------------------------------------------------------------
set(UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable UtilsLite install rules inside Clothoids" FORCE)
set(UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable UtilsLite test suite" FORCE)

FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        Utils2.0
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(UtilsLite "${CLOTHOIDS_UTILSLITE_LOCAL_DIR}")

# Tratta gli header di UtilsLite come SYSTEM per evitare warning
foreach(_utils_tgt utils_UtilsLite_Static utils_UtilsLite)
  if(TARGET ${_utils_tgt})
    get_target_property(_utils_inc ${_utils_tgt} INTERFACE_INCLUDE_DIRECTORIES)
    if(_utils_inc)
      set_target_properties(
        ${_utils_tgt} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_utils_inc}"
      )
    endif()
  endif()
endforeach()

# ----------------------------------------------------------------------------
# 3. quarticRootsFlocke
# ----------------------------------------------------------------------------
FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        1.2.0
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(quarticRootsFlocke "${CLOTHOIDS_QUARTIC_ROOTS_LOCAL_DIR}")

# Raccogli i target delle dipendenze
set(ROOTS utils_quarticRootsFlocke_Static)
set(UTILS utils::UtilsLite_Static)
set(GC GenericContainer::Yaml)

# ----------------------------------------------------------------------------
# Funzione per impostare le directory di output per un target
# ----------------------------------------------------------------------------
function(clothoids_set_target_output_dirs TARGET_NAME LIB_DIR)
  set_target_properties(${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# ----------------------------------------------------------------------------
# Funzione per impostare le directory di output per le dipendenze (alias)
# ----------------------------------------------------------------------------
function(clothoids_set_dependency_output_dirs LIB_DIR)
  foreach(tgt IN LISTS ARGN)
    get_target_property(real_tgt ${tgt} ALIASED_TARGET)
    if(NOT real_tgt)
      set(real_tgt ${tgt})
    endif()
    if(TARGET ${real_tgt})
      get_target_property(tgt_type ${real_tgt} TYPE)
      if(NOT tgt_type STREQUAL "INTERFACE_LIBRARY")
        set_target_properties(${real_tgt} PROPERTIES
          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        )
      endif()
    endif()
  endforeach()
endfunction()

# ----------------------------------------------------------------------------
# Raccolta delle fonti (come originale)
# ----------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB SOURCES "${SOURCE_DIR}/*.cc")
# Aggiungi anche eventuali sottodirectory se necessario (es. src/*/*.cc)
file(GLOB_RECURSE SOURCES_RECURSE "${SOURCE_DIR}/*.cc")
list(APPEND SOURCES ${SOURCES_RECURSE})
list(REMOVE_DUPLICATES SOURCES)

# Header pubblici da rendere disponibili dopo rake build.
# Alcuni header storici del progetto sono in include/, altri in src/.
# Li copiamo entrambi in lib/include preservando la struttura relativa.
set(CLOTHOIDS_PROJECT_HEADER_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${SOURCE_DIR}"
)

set(PUBLIC_HEADERS)
foreach(_header_dir IN LISTS CLOTHOIDS_PROJECT_HEADER_DIRS)
  if(EXISTS "${_header_dir}")
    file(
      GLOB_RECURSE _headers
      CONFIGURE_DEPENDS
      "${_header_dir}/*.h"
      "${_header_dir}/*.hh"
      "${_header_dir}/*.hpp"
      "${_header_dir}/*.hxx"
    )
    list(APPEND PUBLIC_HEADERS ${_headers})
  endif()
endforeach()
list(REMOVE_DUPLICATES PUBLIC_HEADERS)

foreach(_header IN LISTS PUBLIC_HEADERS)
  set_source_files_properties("${_header}" PROPERTIES HEADER_FILE_ONLY TRUE)
endforeach()

# ----------------------------------------------------------------------------
# Target Clothoids
# ----------------------------------------------------------------------------
add_library(Clothoids ${SOURCES} ${PUBLIC_HEADERS})
add_library(Clothoids::Clothoids ALIAS Clothoids)

clothoids_set_target_output_dirs(Clothoids lib/lib)

target_compile_features(Clothoids PUBLIC cxx_std_20)
target_compile_definitions(Clothoids PUBLIC FMT_HEADER_ONLY)

# Include directories: pubblico per gli utenti, privato per i sorgenti
target_include_directories(
  Clothoids
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE
    "${SOURCE_DIR}"
)

# Collegamento alle dipendenze
target_link_libraries(
  Clothoids
  PUBLIC
    ${ROOTS}
    ${GC}
    ${UTILS}
)

# Flag di compilazione
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
  target_compile_options(Clothoids PRIVATE -Wall -Wextra)
  if(CLOTHOIDS_STRICT_WARNINGS)
    target_compile_options(
      Clothoids PRIVATE
      -Wconversion -Wsign-conversion -Wshadow -Wdouble-promotion
    )
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(Clothoids PUBLIC -stdlib=libc++)
  target_link_options(Clothoids PUBLIC -stdlib=libc++)
endif()

if(MSVC)
  target_compile_options(Clothoids PUBLIC /utf-8 /bigobj)
endif()

# Imposta versioni
set_target_properties(
  Clothoids PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Clothoids
)

# Imposta directory di output per le dipendenze (vanno in lib3rd/lib)
clothoids_set_dependency_output_dirs(lib3rd/lib
  GenericContainer
  ${UTILS}
  ${ROOTS}
)


# ----------------------------------------------------------------------------
# Copia locale degli header dopo la compilazione
# ----------------------------------------------------------------------------
# Le librerie sono prodotte direttamente in lib/lib e lib3rd/lib durante il
# build; quindi gli header devono essere copiati con la stessa semantica, senza
# dipendere solo dal target install.
set(CLOTHOIDS_HEADER_COPY_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/copy_clothoids_headers.cmake")
file(WRITE "${CLOTHOIDS_HEADER_COPY_SCRIPT}" [=[
function(_clothoids_copy_header_tree source_dir destination_dir label)
  if(NOT source_dir OR NOT EXISTS "${source_dir}")
    message(STATUS "Skipping ${label} headers: ${source_dir} not found")
    return()
  endif()

  file(MAKE_DIRECTORY "${destination_dir}")
  file(COPY "${source_dir}/" DESTINATION "${destination_dir}"
    FILES_MATCHING
      PATTERN "*.h"
      PATTERN "*.hh"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
  )
endfunction()

file(MAKE_DIRECTORY
  "@CLOTHOIDS_LOCAL_INCLUDE_DIR@"
  "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@"
)

_clothoids_copy_header_tree("@CMAKE_CURRENT_SOURCE_DIR@/include" "@CLOTHOIDS_LOCAL_INCLUDE_DIR@" "Clothoids include")
_clothoids_copy_header_tree("@SOURCE_DIR@" "@CLOTHOIDS_LOCAL_INCLUDE_DIR@" "Clothoids src")
_clothoids_copy_header_tree("@GenericContainer_SOURCE_DIR@/include" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "GenericContainer include")
_clothoids_copy_header_tree("@GenericContainer_SOURCE_DIR@/src" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "GenericContainer src")
_clothoids_copy_header_tree("@UtilsLite_SOURCE_DIR@/include" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "UtilsLite include")
_clothoids_copy_header_tree("@UtilsLite_SOURCE_DIR@/src" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "UtilsLite src")
_clothoids_copy_header_tree("@quarticRootsFlocke_SOURCE_DIR@/include" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "quarticRootsFlocke include")
_clothoids_copy_header_tree("@quarticRootsFlocke_SOURCE_DIR@/src" "@CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR@" "quarticRootsFlocke src")
]=])
configure_file(
  "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  @ONLY
)

add_custom_target(Clothoids_copy_headers ALL
  COMMAND "${CMAKE_COMMAND}" -P "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  COMMENT "Copying Clothoids headers to lib/include and dependency headers to lib3rd/include"
  VERBATIM
)
add_dependencies(Clothoids Clothoids_copy_headers)

# ----------------------------------------------------------------------------
# Alias con suffisso di piattaforma (come in Splines e GenericContainer)
# ----------------------------------------------------------------------------
if(APPLE)
  execute_process(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
  set(CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE})
  set(CLOTHOIDS_PLATFORM_SUFFIX osx)
elseif(WIN32)
  if(MINGW)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(CLOTHOIDS_PLATFORM_SUFFIX mingw64)
    else()
      set(CLOTHOIDS_PLATFORM_SUFFIX mingw32)
    endif()
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CLOTHOIDS_PLATFORM_SUFFIX win64)
  else()
    set(CLOTHOIDS_PLATFORM_SUFFIX win32)
  endif()
elseif(UNIX)
  set(CLOTHOIDS_PLATFORM_SUFFIX linux)
else()
  string(TOLOWER "${CMAKE_SYSTEM_NAME}" CLOTHOIDS_PLATFORM_SUFFIX)
endif()

function(clothoids_set_platform_output_name TARGET_NAME BASE_NAME)
  if(BUILD_SHARED_LIBS)
    set(_output_name "${BASE_NAME}_${CLOTHOIDS_PLATFORM_SUFFIX}")
  else()
    set(_output_name "${BASE_NAME}_${CLOTHOIDS_PLATFORM_SUFFIX}_static")
  endif()

  set_target_properties(
    ${TARGET_NAME} PROPERTIES
    OUTPUT_NAME         "${_output_name}"
    ARCHIVE_OUTPUT_NAME "${_output_name}"
    LIBRARY_OUTPUT_NAME "${_output_name}"
    RUNTIME_OUTPUT_NAME "${_output_name}"
  )
endfunction()

function(clothoids_add_generic_library_symlink TARGET_NAME BASE_NAME)
  if(BUILD_SHARED_LIBS)
    set(REAL_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}_${CLOTHOIDS_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(LINK_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  else()
    set(REAL_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}_${CLOTHOIDS_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(LINK_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}${BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  endif()

  if(WIN32)
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${TARGET_NAME}>"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LINK_FILENAME}"
      VERBATIM
    )

    install(
      FILES "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LINK_FILENAME}"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  else()
    add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LINK_FILENAME}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${REAL_FILENAME}"
              "$<TARGET_FILE_DIR:${TARGET_NAME}>/${LINK_FILENAME}"
      VERBATIM
    )

    # The install-time symlink is handled in the CLOTHOIDS_INSTALL block
    # below, using explicit project-local paths.  Do not use
    # ${CMAKE_INSTALL_PREFIX} here: in install(CODE) it can generate
    # malformed relative paths such as lib/Users/... when the install cache
    # has stale or relative values.
  endif()
endfunction()

if(BUILD_SHARED_LIBS)
  set(CLOTHOIDS_REAL_LIBRARY_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}Clothoids_${CLOTHOIDS_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(CLOTHOIDS_GENERIC_LIBRARY_FILENAME "${CMAKE_SHARED_LIBRARY_PREFIX}Clothoids${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  set(CLOTHOIDS_REAL_LIBRARY_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}Clothoids_${CLOTHOIDS_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(CLOTHOIDS_GENERIC_LIBRARY_FILENAME "${CMAKE_STATIC_LIBRARY_PREFIX}Clothoids${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()
set(CLOTHOIDS_INSTALL_LIB_ABS "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib")

clothoids_set_platform_output_name(Clothoids Clothoids)
clothoids_add_generic_library_symlink(Clothoids Clothoids)

# ----------------------------------------------------------------------------
# Test
# ----------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  file(GLOB CLOTHOIDS_TEST_SOURCES CONFIGURE_DEPENDS
       "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/test*.cc")
  list(SORT CLOTHOIDS_TEST_SOURCES)
  add_custom_target("${PROJECT_NAME}_all_tests" ALL)

  foreach(TEST_SOURCE ${CLOTHOIDS_TEST_SOURCES})
    get_filename_component(S "${TEST_SOURCE}" NAME_WE)
    add_executable(${S} "${TEST_SOURCE}")
    target_link_libraries(${S} PRIVATE Clothoids::Clothoids ${CMAKE_DL_LIBS})
    set_target_properties(
      ${S} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    add_test(NAME "${S}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${S}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    add_dependencies("${PROJECT_NAME}_all_tests" ${S})
  endforeach()
endif()

# ----------------------------------------------------------------------------
# Benchmarks (opzionali)
# ----------------------------------------------------------------------------
if(CLOTHOIDS_BUILD_BENCHMARKS AND PROJECT_IS_TOP_LEVEL)
  # Aggiungi qui i target per i benchmark se presenti, ad esempio:
  # add_executable(bench_clothoids "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/bench.cc")
  # target_link_libraries(bench_clothoids PRIVATE Clothoids::Clothoids)
  # set_target_properties(...)
endif()

# ----------------------------------------------------------------------------
# Installazione
# ----------------------------------------------------------------------------
if(CLOTHOIDS_INSTALL)
  install(
    TARGETS Clothoids
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(CODE
    "set(_clothoids_libdir \"${CLOTHOIDS_INSTALL_LIB_ABS}\")\n
    file(MAKE_DIRECTORY \"\${_clothoids_libdir}\")\n
    execute_process(
      COMMAND \"${CMAKE_COMMAND}\" -E rm -f
              \"\${_clothoids_libdir}/${CLOTHOIDS_GENERIC_LIBRARY_FILENAME}\"
    )\n
    if(NOT WIN32)
      execute_process(
        COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
                \"${CLOTHOIDS_REAL_LIBRARY_FILENAME}\"
                \"\${_clothoids_libdir}/${CLOTHOIDS_GENERIC_LIBRARY_FILENAME}\"
      )
    endif()\n
    execute_process(
      COMMAND \"${CMAKE_COMMAND}\" -E rm -f
              \"${CMAKE_CURRENT_SOURCE_DIR}/lib/${CLOTHOIDS_REAL_LIBRARY_FILENAME}\"
              \"${CMAKE_CURRENT_SOURCE_DIR}/lib/${CLOTHOIDS_GENERIC_LIBRARY_FILENAME}\"
    )\n
    file(REMOVE_RECURSE \"${CMAKE_CURRENT_SOURCE_DIR}/lib/Users\")"
  )
  # Copia header del progetto in lib/include anche durante cmake --install.
  foreach(_header_dir IN LISTS CLOTHOIDS_PROJECT_HEADER_DIRS)
    if(EXISTS "${_header_dir}")
      install(
        DIRECTORY "${_header_dir}/"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
          PATTERN "*.h"
          PATTERN "*.hh"
          PATTERN "*.hpp"
          PATTERN "*.hxx"
      )
    endif()
  endforeach()

  # Ripeti la copia esplicita anche durante cmake --install.
  # Uso install(SCRIPT) per evitare problemi di quoting in install(CODE).
  install(SCRIPT "${CLOTHOIDS_HEADER_COPY_SCRIPT}")

  # CPack (opzionale)
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
  include(CPack)
endif()

# ----------------------------------------------------------------------------
# Riepilogo
# ----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Clothoids dependency summary:")
message(STATUS "  GenericContainer  : ${GenericContainer_SOURCE}")
message(STATUS "  UtilsLite         : ${UtilsLite_SOURCE}")
message(STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}")
message(STATUS "")