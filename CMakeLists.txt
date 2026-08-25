############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#                                                                          #
#  Clothoids -- self-contained modern CMake build.                         #
#                                                                          #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Department of Industrial Engineering                                #
#      University of Trento                                                #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.25)

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# On Apple platforms a single library can contain slices for several CPU
# architectures, for example arm64 and x86_64.  This is controlled by
# CMAKE_OSX_ARCHITECTURES.
#
# The variable must be set before project(), otherwise CMake may initialize
# compilers and targets before the requested architectures are known.
#
# Default behavior:
#   - ON  on Apple hosts or Apple cross-compilation targets.
#   - OFF on Linux/Windows and other platforms, where there is no equivalent
#     single-file universal library format.
#
# User overrides:
#   -DUTILS_BUILD_FAT=OFF
#   -DUTILS_FAT_ARCHITECTURES=arm64
#   -DCMAKE_OSX_ARCHITECTURES=arm64\;x86_64
set(UTILS_BUILD_FAT_DEFAULT OFF)
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  set(UTILS_BUILD_FAT_DEFAULT ON)
endif()
option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")

if(UTILS_BUILD_FAT AND NOT CMAKE_OSX_ARCHITECTURES)
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
endif()

# ----------------------------------------------------------------------------
# Project version discovery
# ----------------------------------------------------------------------------
# The project version is derived from the most recent Git tag when the source
# tree is a Git checkout and the Git executable is available.  This keeps CMake
# package metadata, install metadata, and generated archives aligned with the
# repository tag used for a release.
#
# If Git is unavailable, or if the source tree is a release archive without tag
# information, the script falls back to a neutral 0.0.0 version so configuration
# remains deterministic and does not fail for users building from snapshots.
#
# Tags may be written either as "1.2.3" or "v1.2.3".  The optional leading
# "v" is stripped before passing the version to project().
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
# Main user-facing build options
# ----------------------------------------------------------------------------
# Single-configuration generators such as Unix Makefiles and Ninja use
# CMAKE_BUILD_TYPE.  If the user did not choose a build type explicitly, default
# to Release so the produced local library is optimized.  Multi-configuration
# generators, such as Visual Studio and Ninja Multi-Config, ignore this variable
# and select the configuration at build time instead.
# ----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# BUILD_SHARED_LIBS is the standard CMake switch used by add_library() when no
# explicit STATIC or SHARED type is given.  Clothoids defaults to a static
# library, matching the existing local-library workflow.
option(BUILD_SHARED_LIBS "Build Clothoids as a shared library" OFF)

# Tests, install rules, and benchmarks default to ON only when this source tree
# is configured directly.  When Clothoids is pulled into another project with
# FetchContent or add_subdirectory(), these optional targets stay disabled by
# default to avoid surprising the parent build.
option(BUILD_TESTING "Build the Clothoids test executables" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_BUILD_BENCHMARKS "Build the micro-benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(CLOTHOIDS_POPULATE_TOOLBOX "Refresh toolbox/src during a normal build" OFF)
option(CLOTHOIDS_ALLOW_NETWORK_FETCH "Allow FetchContent network fallbacks for missing sibling projects" OFF)

# Strict warnings are useful for maintenance but are opt-in because they may be
# too aggressive for users building with newer compilers or vendor toolchains.
option(CLOTHOIDS_STRICT_WARNINGS "Compile with narrowing/conversion warnings" OFF)
option(
  CLOTHOIDS_ENABLE_GENERIC_CONTAINER_JSON
  "Enable the optional GenericContainer JSON bridge"
  OFF
)

# CheckIPOSupported is needed by the IPO/LTO probing below.
include(CheckIPOSupported)

# ----------------------------------------------------------------------------
# Compilation and optimization options
# ----------------------------------------------------------------------------
# Native CPU tuning is enabled for ordinary single-architecture builds, but
# disabled by default for Apple universal binaries: -march=native cannot
# describe both arm64 and x86_64 slices in the same artifact.
set(CLOTHOIDS_APPLE_MULTI_ARCH FALSE)
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
  list(LENGTH CMAKE_OSX_ARCHITECTURES CLOTHOIDS_OSX_ARCHITECTURE_COUNT)
  if(CLOTHOIDS_OSX_ARCHITECTURE_COUNT GREATER 1)
    set(CLOTHOIDS_APPLE_MULTI_ARCH TRUE)
  endif()
endif()

set(CLOTHOIDS_USE_NATIVE_ARCH_DEFAULT ON)
if(CLOTHOIDS_APPLE_MULTI_ARCH)
  set(CLOTHOIDS_USE_NATIVE_ARCH_DEFAULT OFF)
endif()

option(
  CLOTHOIDS_USE_NATIVE_ARCH
  "Optimize Release builds for the current CPU architecture"
  ${CLOTHOIDS_USE_NATIVE_ARCH_DEFAULT}
)
option(CLOTHOIDS_ENABLE_LTO "Enable interprocedural optimization / LTO" ON)
option(CLOTHOIDS_UNROLL_LOOPS "Enable explicit loop unrolling in Release builds" OFF)
option(
  CLOTHOIDS_ENABLE_FAST_MATH
  "Enable aggressive floating-point optimizations; may change IEEE-754 semantics"
  OFF
)

# Probe IPO/LTO support once and reuse the result for every Clothoids target.
set(CLOTHOIDS_IPO_SUPPORTED FALSE)
set(CLOTHOIDS_IPO_STATUS "disabled")
if(CLOTHOIDS_ENABLE_LTO)
  check_ipo_supported(
    RESULT CLOTHOIDS_IPO_SUPPORTED
    OUTPUT CLOTHOIDS_IPO_ERROR
    LANGUAGES CXX
  )
  if(CLOTHOIDS_IPO_SUPPORTED)
    set(CLOTHOIDS_IPO_STATUS "enabled")
  else()
    set(CLOTHOIDS_IPO_STATUS "not supported: ${CLOTHOIDS_IPO_ERROR}")
  endif()
endif()

# ----------------------------------------------------------------------------
# Install prefix and directory layout
# ----------------------------------------------------------------------------
# The historical project layout installs into a project-local lib/ tree rather
# than into a system prefix such as /usr/local.  When the user has not supplied
# a custom install prefix, use <source>/lib and define GNUInstallDirs-compatible
# subdirectories below it.
# ----------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Install include directory" FORCE)
set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Install library directory" FORCE)
set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "Install binary directory" FORCE)

# These absolute local paths are used by custom build/install helper commands.
# They make the Rake-style local development layout explicit:
#   - project headers are copied to lib/include;
#   - project libraries are produced in lib/lib;
#   - runtime executables are produced in lib/bin;
#   - third-party headers are copied to lib3rd/include.
#
# Keeping libraries only in lib/lib avoids duplicate artifacts directly under
# lib/, which was a common source of stale files in older build layouts.
set(CLOTHOIDS_LOCAL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/include")
set(CLOTHOIDS_LOCAL_LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib")
set(CLOTHOIDS_LOCAL_RUNTIME_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/bin")
set(CLOTHOIDS_LOCAL_3RD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib3rd/include")

# GNUInstallDirs defines standard installation directory variables such as
# CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_LIBDIR, and CMAKE_INSTALL_BINDIR.
# FetchContent provides the dependency population mechanism used below.
include(GNUInstallDirs)
include(FetchContent)

# ----------------------------------------------------------------------------
# Dependency helper: use a local checkout or fall back to FetchContent
# ----------------------------------------------------------------------------
# clothoids_use_local_or_fetch(NAME LOCAL_DIR EXPECTED_REF) centralizes the dependency
# loading policy used for GenericContainer, UtilsLite, and quarticRootsFlocke.
#
# Resolution order:
#   1. If LOCAL_DIR contains a CMakeLists.txt and either belongs to the same
#      Git worktree as Clothoids or resolves to EXPECTED_REF, use it locally.
#   2. Otherwise, use the FetchContent declaration already registered by the
#      caller and let CMake download/populate the dependency.
#
# The helper also suppresses dependency install rules while a dependency is
# included as part of this super-build.  This prevents nested projects from
# installing their own files when only the Clothoids install target should
# control the final layout.
# ----------------------------------------------------------------------------
function(clothoids_use_local_or_fetch NAME LOCAL_DIR EXPECTED_REF)
  # Normalize the user-supplied or default local path so existence checks and
  # diagnostic messages do not depend on the current working directory.
  get_filename_component(
    _local_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )

  # FetchContent override variables are conventionally upper-case, while the
  # generated _deps build directories use lower-case names.
  string(TOUPPER "${NAME}" _name_upper)
  string(TOLOWER "${NAME}" _name_lower)

  if(FETCHCONTENT_BASE_DIR)
    set(_fetchcontent_base_dir "${FETCHCONTENT_BASE_DIR}")
  else()
    set(_fetchcontent_base_dir "${CMAKE_BINARY_DIR}/_deps")
  endif()

  # Dependencies included through FetchContent may create their own install
  # targets.  Temporarily skip those rules while the dependency is made
  # available; Clothoids owns the final install layout.
  set(_dependency_install_stub "${_fetchcontent_base_dir}/${_name_lower}-build/cmake_install.cmake")
  set(_saved_skip_install_rules "${CMAKE_SKIP_INSTALL_RULES}")
  set(CMAKE_SKIP_INSTALL_RULES ON)

  set(_local_matches FALSE)
  set(_local_branch "not available")
  set(_local_state "not a compatible Git checkout")
  if(EXISTS "${_local_abs}/CMakeLists.txt" AND GIT_EXECUTABLE)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" -C "${_local_abs}" symbolic-ref --quiet --short HEAD
      RESULT_VARIABLE _branch_result
      OUTPUT_VARIABLE _local_branch
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(NOT _branch_result EQUAL 0)
      set(_local_branch "detached HEAD")
    endif()

    execute_process(
      COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}" rev-parse --show-toplevel
      RESULT_VARIABLE _project_root_result
      OUTPUT_VARIABLE _project_git_root
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" -C "${_local_abs}" rev-parse --show-toplevel
      RESULT_VARIABLE _local_root_result
      OUTPUT_VARIABLE _local_git_root
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )

    if(
      _project_root_result EQUAL 0 AND
      _local_root_result EQUAL 0 AND
      _project_git_root STREQUAL _local_git_root
    )
      set(_local_matches TRUE)
      set(_local_state "shared Git worktree, branch '${_local_branch}'")
    elseif(_branch_result EQUAL 0 AND _local_branch STREQUAL EXPECTED_REF)
      set(_local_matches TRUE)
      set(_local_state "branch '${_local_branch}'")
    else()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${_local_abs}" rev-parse HEAD
        RESULT_VARIABLE _head_result
        OUTPUT_VARIABLE _head_commit
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${_local_abs}" rev-parse "${EXPECTED_REF}^{commit}"
        RESULT_VARIABLE _ref_result
        OUTPUT_VARIABLE _ref_commit
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      if(_head_result EQUAL 0 AND _ref_result EQUAL 0 AND _head_commit STREQUAL _ref_commit)
        set(_local_matches TRUE)
        set(_local_state "HEAD at ref '${EXPECTED_REF}'")
      elseif(_head_result EQUAL 0 AND _ref_result EQUAL 0)
        execute_process(
          COMMAND "${GIT_EXECUTABLE}" -C "${_local_abs}" merge-base --is-ancestor "${_ref_commit}" "${_head_commit}"
          RESULT_VARIABLE _ref_is_ancestor
          ERROR_QUIET
        )
        if(_ref_is_ancestor EQUAL 0)
          set(_local_matches TRUE)
          set(_local_state "branch '${_local_branch}' contains ref '${EXPECTED_REF}'")
        else()
          set(_local_state "branch '${_local_branch}'")
        endif()
      else()
        set(_local_state "branch '${_local_branch}'")
      endif()
    endif()
  endif()

  if(_local_matches)
    message(STATUS "Clothoids: ${NAME} -> local ${_local_abs} (${_local_state}, requested '${EXPECTED_REF}')")
    set(FETCHCONTENT_SOURCE_DIR_${_name_upper} "${_local_abs}" CACHE PATH "Local ${NAME} checkout" FORCE)
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "local: ${_local_abs} [${_local_state}]" PARENT_SCOPE)
    set("${NAME}_SOURCE_DIR" "${_local_abs}" PARENT_SCOPE)
  else()
    unset(FETCHCONTENT_SOURCE_DIR_${_name_upper} CACHE)
    if(EXISTS "${_local_abs}/CMakeLists.txt")
      message(STATUS "Clothoids: ${NAME} -> ignoring local ${_local_abs} (${_local_state}, requested '${EXPECTED_REF}')")
      set(_local_problem "${_local_state} does not match requested ref '${EXPECTED_REF}'")
    else()
      set(_local_problem "checkout not found at ${_local_abs}")
    endif()
    if(NOT CLOTHOIDS_ALLOW_NETWORK_FETCH)
      message(FATAL_ERROR
        "Cannot use local ${NAME}: ${_local_problem}. "
        "Enable CLOTHOIDS_ALLOW_NETWORK_FETCH to fetch '${EXPECTED_REF}'.")
    endif()
    message(STATUS "Clothoids: ${NAME} -> FetchContent (requested '${EXPECTED_REF}')")
    FetchContent_MakeAvailable(${NAME})
    set("${NAME}_SOURCE" "FetchContent: ${EXPECTED_REF}" PARENT_SCOPE)
    if(DEFINED "${_name_lower}_SOURCE_DIR")
      set("${NAME}_SOURCE_DIR" "${${_name_lower}_SOURCE_DIR}" PARENT_SCOPE)
    endif()
  endif()

  # Some dependency projects or generated install steps expect a
  # cmake_install.cmake file in their build directory.  Write a harmless stub
  # when install rules have been suppressed, keeping the super-build robust.
  file(MAKE_DIRECTORY "${_fetchcontent_base_dir}/${_name_lower}-build")
  if(NOT EXISTS "${_dependency_install_stub}")
    file(WRITE "${_dependency_install_stub}" "# Install disabled for dependency ${NAME} inside Clothoids.\n")
  endif()

  set(CMAKE_SKIP_INSTALL_RULES "${_saved_skip_install_rules}")
endfunction()

# Dependencies loaded by a nested project cannot reliably identify their own
# branch when all sources have been copied into one aggregate Git worktree.
# Preselect a local source when it belongs to Clothoids' worktree, or when its
# standalone checkout resolves to the requested ref.  FetchContent then passes
# the source directly to the nested project, bypassing its branch-name check.
function(clothoids_override_nested_source NAME LOCAL_DIR EXPECTED_REF)
  get_filename_component(
    _nested_abs
    "${LOCAL_DIR}"
    ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  if(NOT EXISTS "${_nested_abs}/CMakeLists.txt" OR NOT GIT_EXECUTABLE)
    return()
  endif()

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}" rev-parse --show-toplevel
    RESULT_VARIABLE _project_root_result
    OUTPUT_VARIABLE _project_git_root
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${_nested_abs}" rev-parse --show-toplevel
    RESULT_VARIABLE _nested_root_result
    OUTPUT_VARIABLE _nested_git_root
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${_nested_abs}" rev-parse HEAD
    RESULT_VARIABLE _nested_head_result
    OUTPUT_VARIABLE _nested_head
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${_nested_abs}" rev-parse "${EXPECTED_REF}^{commit}"
    RESULT_VARIABLE _nested_ref_result
    OUTPUT_VARIABLE _nested_ref
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  set(_nested_matches FALSE)
  set(_nested_state "does not match '${EXPECTED_REF}'")
  if(
    _project_root_result EQUAL 0 AND
    _nested_root_result EQUAL 0 AND
    _project_git_root STREQUAL _nested_git_root
  )
    set(_nested_matches TRUE)
    set(_nested_state "shared Git worktree")
  elseif(
    _nested_head_result EQUAL 0 AND
    _nested_ref_result EQUAL 0 AND
    _nested_head STREQUAL _nested_ref
  )
    set(_nested_matches TRUE)
    set(_nested_state "HEAD at ref '${EXPECTED_REF}'")
  endif()

  if(_nested_matches)
    string(TOUPPER "${NAME}" _nested_name_upper)
    set(
      FETCHCONTENT_SOURCE_DIR_${_nested_name_upper}
      "${_nested_abs}"
      CACHE PATH "Local source override for nested dependency ${NAME}"
      FORCE
    )
    message(STATUS "Clothoids: nested ${NAME} -> local ${_nested_abs} (${_nested_state})")
  endif()
endfunction()

# ----------------------------------------------------------------------------
# Local dependency checkout paths
# ----------------------------------------------------------------------------
# These cache variables allow a developer to keep sibling repositories next to
# Clothoids and reuse them without network access.  Each value can be overridden
# from the command line with -D<VAR>=<path>.  If the path is absent or invalid,
# the corresponding FetchContent declaration below is used instead.
# ----------------------------------------------------------------------------
set(CLOTHOIDS_GC_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../GC" CACHE PATH "Local GenericContainer checkout")
set(CLOTHOIDS_UTILSLITE_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../UtilsLite" CACHE PATH "Local UtilsLite checkout")
set(CLOTHOIDS_QUARTIC_ROOTS_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../quarticRootsFlocke" CACHE PATH "Local quarticRootsFlocke checkout")
set(CLOTHOIDS_FKYAML_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../fkYAML" CACHE PATH "Local fkYAML checkout used by GenericContainer")
set(CLOTHOIDS_GC_REF "GC2.0" CACHE STRING "GenericContainer branch/tag to fetch")
set(CLOTHOIDS_UTILSLITE_REF "Utils2.0" CACHE STRING "UtilsLite branch/tag to fetch")
set(CLOTHOIDS_QUARTIC_ROOTS_REF "1.2.0" CACHE STRING "quarticRootsFlocke branch/tag to fetch")

# ----------------------------------------------------------------------------
# Dependency 1: GenericContainer
# ----------------------------------------------------------------------------
# GenericContainer provides heterogeneous container types and YAML/JSON/TOML/Lua
# bridge targets used by higher-level examples and serialization paths.  Its
# install rules and examples are disabled here because Clothoids is the outer
# project controlling the build tree and install tree.
# ----------------------------------------------------------------------------
set(GENERIC_CONTAINER_INSTALL OFF CACHE BOOL "Disable GenericContainer install rules inside Clothoids" FORCE)
set(GENERIC_CONTAINER_BUILD_EXAMPLES OFF CACHE BOOL "Disable GenericContainer examples inside Clothoids" FORCE)
# Clothoids links only the YAML bridge. Disable the unused JSON, TOML and Lua
# bridges before adding GenericContainer so their dependencies are neither
# searched nor downloaded. Force YAML on because its target is required below.
set(GENERIC_CONTAINER_ENABLE_YAML ON CACHE BOOL "Enable GenericContainer YAML bridge inside Clothoids" FORCE)
set(
  GENERIC_CONTAINER_ENABLE_JSON
  ${CLOTHOIDS_ENABLE_GENERIC_CONTAINER_JSON}
  CACHE BOOL "Enable GenericContainer JSON bridge inside Clothoids"
  FORCE
)
set(GENERIC_CONTAINER_ENABLE_TOML OFF CACHE BOOL "Disable GenericContainer TOML bridge inside Clothoids" FORCE)
set(GENERIC_CONTAINER_ENABLE_LUA OFF CACHE BOOL "Disable GenericContainer Lua bridge inside Clothoids" FORCE)

clothoids_override_nested_source(fkYAML "${CLOTHOIDS_FKYAML_LOCAL_DIR}" "v0.4.3")
if(FETCHCONTENT_SOURCE_DIR_FKYAML)
  FetchContent_Declare(fkYAML SOURCE_DIR "${FETCHCONTENT_SOURCE_DIR_FKYAML}")
  FetchContent_MakeAvailable(fkYAML)
endif()

FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/ebertolazzi/GenericContainer.git
  GIT_TAG        ${CLOTHOIDS_GC_REF}
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(GenericContainer "${CLOTHOIDS_GC_LOCAL_DIR}" "${CLOTHOIDS_GC_REF}")

# Clang compatibility adjustment.
# Some dependency targets may not propagate the standard C++ library choice in
# the same way as the top-level target.  When using Clang, explicitly propagate
# libc++ to the GenericContainer targets that Clothoids may link against.  The
# code distinguishes INTERFACE targets from compiled targets because CMake uses
# different visibility keywords for those cases.
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
# Dependency 2: UtilsLite
# ----------------------------------------------------------------------------
# UtilsLite contains common utility headers and implementation code shared by
# this project family.  Tests and local install rules are disabled because they
# are not needed when UtilsLite is embedded as a dependency of Clothoids.
# ----------------------------------------------------------------------------
set(UTILS_INSTALL_LOCAL OFF CACHE BOOL "Disable UtilsLite install rules inside Clothoids" FORCE)
set(UTILS_ENABLE_TESTS OFF CACHE BOOL "Disable UtilsLite test suite" FORCE)

FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/ebertolazzi/UtilsLite.git
  GIT_TAG        ${CLOTHOIDS_UTILSLITE_REF}
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(UtilsLite "${CLOTHOIDS_UTILSLITE_LOCAL_DIR}" "${CLOTHOIDS_UTILSLITE_REF}")

# Treat UtilsLite include directories as SYSTEM include directories.
# This prevents third-party or dependency warnings from being reported as if
# they were warnings in Clothoids itself, especially when strict warning flags
# are enabled on the top-level target.
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
# Dependency 3: quarticRootsFlocke
# ----------------------------------------------------------------------------
# quarticRootsFlocke supplies polynomial root solvers used internally by the
# geometric algorithms.  It is declared through FetchContent but may be resolved
# to a local sibling checkout by clothoids_use_local_or_fetch().
# ----------------------------------------------------------------------------
FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        ${CLOTHOIDS_QUARTIC_ROOTS_REF}
  GIT_SHALLOW    TRUE
)
clothoids_use_local_or_fetch(quarticRootsFlocke "${CLOTHOIDS_QUARTIC_ROOTS_LOCAL_DIR}" "${CLOTHOIDS_QUARTIC_ROOTS_REF}")

# Canonical dependency target names used by the Clothoids target.
# These variables keep the final target_link_libraries() call compact and make
# the expected imported or embedded target names explicit in one place.
set(ROOTS utils_quarticRootsFlocke_Static)
set(UTILS utils::UtilsLite_Static)
set(GC GenericContainer::Yaml)

# ----------------------------------------------------------------------------
# Output directory helper for project targets
# ----------------------------------------------------------------------------
# This helper forces a target's archives, shared libraries, and runtime
# executables into the local project layout.  It is used for the Clothoids target
# so generated artifacts end up in lib/lib and bin rather than in generator-
# specific directories inside the build tree.
# ----------------------------------------------------------------------------
function(clothoids_set_target_output_dirs TARGET_NAME LIB_DIR)
  set_target_properties(${TARGET_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${LIB_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

# ----------------------------------------------------------------------------
# Output directory helper for dependency targets
# ----------------------------------------------------------------------------
# Some dependencies are represented through ALIAS targets.  Since CMake
# properties cannot be set directly on aliases, this function resolves each
# alias to its real target before applying output-directory properties.  It
# skips INTERFACE libraries because they do not produce binary artifacts.
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
# Source discovery
# ----------------------------------------------------------------------------
# Collect all .cc files below src/.  The non-recursive glob preserves the
# original top-level behavior, while the recursive glob also captures optional
# implementation files stored in subdirectories.  Duplicates are removed after
# combining both lists.
#
# CONFIGURE_DEPENDS asks CMake to re-check the glob when the build system is
# regenerated.  This is convenient for this repository layout, although explicit
# source lists remain the most reproducible approach for packaged libraries.
# ----------------------------------------------------------------------------
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB SOURCES "${SOURCE_DIR}/*.cc")
# Also include implementation files in nested source directories, if present.
file(GLOB_RECURSE SOURCES_RECURSE "${SOURCE_DIR}/*.cc")
list(APPEND SOURCES ${SOURCES_RECURSE})
list(REMOVE_DUPLICATES SOURCES)

# Public header discovery.
# Clothoids historically exposes headers from both include/ and src/.  Both
# trees are therefore scanned and later copied to lib/include, preserving the
# relative directory layout expected by downstream projects and existing local
# build scripts.
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
# Main Clothoids library target
# ----------------------------------------------------------------------------
# The Clothoids target is the primary build artifact.  It is created from the
# discovered implementation files plus public headers, and exported through the
# namespaced alias Clothoids::Clothoids so downstream CMake code can link to a
# stable target name independent of the physical library filename.
# ----------------------------------------------------------------------------
add_library(Clothoids ${SOURCES} ${PUBLIC_HEADERS})
add_library(Clothoids::Clothoids ALIAS Clothoids)

clothoids_set_target_output_dirs(Clothoids lib/lib)

target_compile_features(Clothoids PUBLIC cxx_std_20)
target_compile_definitions(Clothoids PUBLIC FMT_HEADER_ONLY)

# Common target configuration: diagnostics, definitions and optimization policy.
# Flags are target-local (PRIVATE) so they do not leak into consumers that link
# this project.
function(clothoids_configure_target TARGET_NAME)
  target_compile_options(
    ${TARGET_NAME}
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wextra>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wpedantic>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wshadow>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wnon-virtual-dtor>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Woverloaded-virtual>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wcast-align>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wunused>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wformat=2>
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wimplicit-fallthrough>
      $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-cond>
      $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-branches>
      $<$<CXX_COMPILER_ID:GNU>:-Wlogical-op>
      $<$<CXX_COMPILER_ID:GNU>:-Wnull-dereference>
      $<$<CXX_COMPILER_ID:GNU>:-Wmisleading-indentation>
      $<$<CXX_COMPILER_ID:MSVC>:/W4>
      $<$<CXX_COMPILER_ID:MSVC>:/permissive->
      $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
  )

  target_compile_options(
    ${TARGET_NAME}
    PRIVATE
      $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3>
      $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/O2>
      $<$<AND:$<CONFIG:RelWithDebInfo>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O2>
      $<$<AND:$<CONFIG:RelWithDebInfo>,$<CXX_COMPILER_ID:MSVC>>:/O2>
  )

  target_compile_definitions(
    ${TARGET_NAME}
    PRIVATE
      $<$<CONFIG:Release>:NDEBUG>
      $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
  )

  if(CLOTHOIDS_USE_NATIVE_ARCH)
    if(CLOTHOIDS_APPLE_MULTI_ARCH)
      message(STATUS
        "Skipping native CPU flags for ${TARGET_NAME}: Apple universal build uses "
        "CMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}"
      )
    else()
      target_compile_options(
        ${TARGET_NAME}
        PRIVATE
          $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-march=native>
          $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-mtune=native>
          $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:IntelLLVM>>:-xHost>
          $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/arch:AVX2>
      )
    endif()
  endif()

  if(CLOTHOIDS_UNROLL_LOOPS)
    target_compile_options(
      ${TARGET_NAME}
      PRIVATE
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-funroll-loops>
    )
  endif()

  # Deliberately OFF by default because it permits reassociation and other
  # transformations that weaken strict IEEE-754 behavior.
  if(CLOTHOIDS_ENABLE_FAST_MATH)
    target_compile_options(
      ${TARGET_NAME}
      PRIVATE
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-ffast-math>
        $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/fp:fast>
    )
  endif()

  if(CLOTHOIDS_IPO_SUPPORTED)
    set_property(TARGET ${TARGET_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set_property(TARGET ${TARGET_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
  endif()

  if(APPLE AND UTILS_BUILD_FAT)
    set_target_properties(
      ${TARGET_NAME}
      PROPERTIES
        OSX_ARCHITECTURES "${CMAKE_OSX_ARCHITECTURES}"
    )
  endif()
endfunction()

clothoids_configure_target(Clothoids)

# Include directories.
# BUILD_INTERFACE paths are used while building from the source tree.  The
# INSTALL_INTERFACE path is used by consumers after installation.  The src/
# directory is also listed privately so internal implementation files can keep
# their historical include style.
target_include_directories(
  Clothoids
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE
    "${SOURCE_DIR}"
)

# Link against the resolved dependency targets.  These are PUBLIC because the
# Clothoids headers and linked library may expose types, compile definitions, or
# transitive usage requirements from the dependency stack.
target_link_libraries(
  Clothoids
  PUBLIC
    ${ROOTS}
    ${GC}
    ${UTILS}
)

# Compiler warnings and platform-specific compile options.
# Extra narrowing and conversion diagnostics are guarded by
# CLOTHOIDS_STRICT_WARNINGS because they can be too noisy for regular
# development builds.  The baseline warning set comes from
# clothoids_configure_target() above.
if(CLOTHOIDS_STRICT_WARNINGS)
  target_compile_options(
    Clothoids PRIVATE
    -Wconversion -Wsign-conversion -Wshadow -Wdouble-promotion
  )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(Clothoids PUBLIC -stdlib=libc++)
  target_link_options(Clothoids PUBLIC -stdlib=libc++)
endif()

if(MSVC)
  target_compile_options(Clothoids PUBLIC /utf-8 /bigobj)
endif()

# Target metadata used by shared-library versioning and CMake export metadata.
# VERSION/SOVERSION affect shared objects on platforms that support sonames;
# EXPORT_NAME controls how the target appears when exported to consumers.
set_target_properties(
  Clothoids PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Clothoids
)

# Place compiled dependency artifacts under lib3rd/lib.  This keeps project
# artifacts and third-party artifacts separate while preserving a predictable
# local tree for the historical build workflow.
clothoids_set_dependency_output_dirs(lib3rd/lib
  GenericContainer
  ${UTILS}
  ${ROOTS}
)


# ----------------------------------------------------------------------------
# Local header and toolbox source copy target
# ----------------------------------------------------------------------------
# Libraries are produced directly in lib/lib and lib3rd/lib during the normal
# build, so headers must be copied with the same build-time semantics.  Relying
# only on cmake --install would leave the local development layout incomplete
# after a regular build.
#
# The script generated below copies project headers to lib/include and selected
# dependency headers to lib3rd/include.  It also collects the C/C++ headers and
# implementation files used by Clothoids into toolbox/src.  Dependency source
# directories come from the helper above, so this works both for local checkouts
# and for sources populated by FetchContent.  A CMake script is used instead of
# shell commands to remain portable across Unix, macOS, and Windows generators.
# ----------------------------------------------------------------------------
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

function(_clothoids_copy_toolbox_tree source_dir destination_dir label)
  if(NOT source_dir OR NOT EXISTS "${source_dir}")
    message(STATUS "Skipping ${label} toolbox sources: ${source_dir} not found")
    return()
  endif()

  file(MAKE_DIRECTORY "${destination_dir}")
  file(COPY "${source_dir}/" DESTINATION "${destination_dir}"
    FILES_MATCHING
      PATTERN "*.c"
      PATTERN "*.cc"
      PATTERN "*.cpp"
      PATTERN "*.cxx"
      PATTERN "*.h"
      PATTERN "*.hh"
      PATTERN "*.hpp"
      PATTERN "*.hxx"
      PATTERN "*.inc"
      PATTERN "*.inl"
      PATTERN "*.ipp"
      PATTERN "*.tcc"
      PATTERN "*.tpp"
  )
endfunction()

# Eigen exposes a number of public module headers without a file extension,
# such as Eigen/Core, Eigen/Dense and unsupported/Eigen/AutoDiff.  These files
# cannot be selected reliably with the suffix-based filter above, therefore the
# resolved Eigen trees vendored by UtilsLite must be copied in their entirety.
function(_clothoids_copy_complete_tree source_dir destination_dir label)
  if(NOT source_dir OR NOT EXISTS "${source_dir}")
    message(STATUS "Skipping ${label}: ${source_dir} not found")
    return()
  endif()

  file(MAKE_DIRECTORY "${destination_dir}")
  file(COPY "${source_dir}/" DESTINATION "${destination_dir}")
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

# Merge the public and implementation trees into the flat toolbox source root
# only when explicitly requested. A normal build must not rewrite source files.
# file(COPY) preserves nested directories (for example Clothoids/, Utils/ and
# GenericContainer/) while files located at the tree root remain directly in
# toolbox/src, matching the historical populate_toolbox.rb layout.
if(@CLOTHOIDS_POPULATE_TOOLBOX@)
file(MAKE_DIRECTORY "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src")
_clothoids_copy_toolbox_tree("@CMAKE_CURRENT_SOURCE_DIR@/include" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "Clothoids include")
_clothoids_copy_toolbox_tree("@SOURCE_DIR@" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "Clothoids src")
_clothoids_copy_toolbox_tree("@GenericContainer_SOURCE_DIR@/include" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "GenericContainer include")
_clothoids_copy_toolbox_tree("@GenericContainer_SOURCE_DIR@/src" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "GenericContainer src")
_clothoids_copy_toolbox_tree("@UtilsLite_SOURCE_DIR@/include" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "UtilsLite include")
_clothoids_copy_toolbox_tree("@UtilsLite_SOURCE_DIR@/src" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "UtilsLite src")
_clothoids_copy_complete_tree("@UtilsLite_SOURCE_DIR@/src/Utils/3rd/Eigen" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src/Utils/3rd/Eigen" "Eigen headers")
_clothoids_copy_complete_tree("@UtilsLite_SOURCE_DIR@/src/Utils/3rd/unsupported" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src/Utils/3rd/unsupported" "Eigen unsupported headers")
_clothoids_copy_toolbox_tree("@quarticRootsFlocke_SOURCE_DIR@/include" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "quarticRootsFlocke include")
_clothoids_copy_toolbox_tree("@quarticRootsFlocke_SOURCE_DIR@/src" "@CMAKE_CURRENT_SOURCE_DIR@/toolbox/src" "quarticRootsFlocke src")
endif()
]=])
configure_file(
  "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  @ONLY
)

add_custom_target(Clothoids_copy_headers ALL
  COMMAND "${CMAKE_COMMAND}" -P "${CLOTHOIDS_HEADER_COPY_SCRIPT}"
  COMMENT "Copying Clothoids/dependency headers (toolbox refresh is opt-in)"
  VERBATIM
)
add_dependencies(Clothoids Clothoids_copy_headers)

# ----------------------------------------------------------------------------
# Platform-specific library naming and compatibility aliases
# ----------------------------------------------------------------------------
# The real library filename carries the platform suffix, for example:
#   libClothoids_osx_static.a
#   libClothoids_linux_static.a
#
# A generic compatibility filename is also created:
#   libClothoids.a
#
# On Unix-like platforms the generic file is a symlink to the real artifact.  On
# Windows it is copied because symlink behavior depends on user privileges and
# filesystem settings.
#
# CMAKE_OSX_ARCHITECTURES is intentionally not set in this block.  Apple
# universal/fat architectures must be selected before project(), in the FAT
# setup block at the top of this file.
# ----------------------------------------------------------------------------
if(APPLE)
  # The platform suffix names the produced artifact.  Do not set
  # CMAKE_OSX_ARCHITECTURES here: Apple universal/fat architectures must be
  # selected before project(), otherwise CMake may have already initialized the
  # compiler and target defaults.
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
  # Compute the real artifact name.  The suffix records both platform and, for
  # static builds, the static-library kind.  This matches the naming convention
  # used by the related Splines, GenericContainer, and UtilsLite projects.
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
  # Create the compatibility filename expected by legacy consumers.  The real
  # artifact keeps the explicit platform suffix; the generic filename points to
  # or copies that artifact after each build.
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

# Precompute the real and generic filenames used by both build-time aliasing and
# install-time cleanup/symlink generation.  Keeping the names in variables avoids
# repeating platform-specific filename logic in install(CODE).
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
# Tests
# ----------------------------------------------------------------------------
# When BUILD_TESTING is enabled and this is the top-level project, every
# src_tests/test*.cc file becomes an executable linked against Clothoids.  The
# executables are emitted into bin/ and registered with CTest so they can be run
# through ctest or through the generated all-tests convenience target.
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
# Optional benchmarks
# ----------------------------------------------------------------------------
# This block is intentionally empty unless benchmark sources are added later.
# Keeping the option and block in place gives the project a stable extension
# point without affecting regular library or test builds.
# ----------------------------------------------------------------------------
if(CLOTHOIDS_BUILD_BENCHMARKS AND PROJECT_IS_TOP_LEVEL)
  # Add benchmark targets here if benchmark sources are introduced later, for example:
  # add_executable(bench_clothoids "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/bench.cc")
  # target_link_libraries(bench_clothoids PRIVATE Clothoids::Clothoids)
  # set_target_properties(...)
endif()

# ----------------------------------------------------------------------------
# Install rules
# ----------------------------------------------------------------------------
# Install the library target, project headers, and the generated header-copy
# script output.  The install-time symlink logic mirrors the build-time
# compatibility alias so both regular builds and cmake --install leave the
# project-local lib/lib directory in the same usable state.
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
  # Copy project headers during cmake --install as well.  This duplicates the
  # build-time header-copy behavior for users who rely exclusively on install.
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

  # Reuse the generated copy script during installation.  install(SCRIPT) avoids
  # fragile quoting in large install(CODE) strings and keeps the copy logic in a
  # single generated file.
  install(SCRIPT "${CLOTHOIDS_HEADER_COPY_SCRIPT}")

  # Optional CPack configuration for source and binary archives generated from
  # the same project metadata used by project() and install().
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
  set(CPACK_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/")
  include(CPack)
endif()

# ----------------------------------------------------------------------------
# Configuration summary
# ----------------------------------------------------------------------------
# Print the dependency resolution summary at the end of configuration.  These
# messages make it immediately visible whether each dependency came from a local
# checkout or from FetchContent, which is useful when debugging reproducibility
# and offline builds.
# ----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Clothoids dependency summary:")
message(STATUS "  GenericContainer  : ${GenericContainer_SOURCE}")
message(STATUS "  UtilsLite         : ${UtilsLite_SOURCE}")
message(STATUS "  quarticRootsFlocke: ${quarticRootsFlocke_SOURCE}")
message(STATUS "  Compiler          : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "  Build shared      : ${BUILD_SHARED_LIBS}")
message(STATUS "  Native arch       : ${CLOTHOIDS_USE_NATIVE_ARCH}")
message(STATUS "  Apple multi-arch  : ${CLOTHOIDS_APPLE_MULTI_ARCH}")
message(STATUS "  LTO requested     : ${CLOTHOIDS_ENABLE_LTO}")
message(STATUS "  LTO status        : ${CLOTHOIDS_IPO_STATUS}")
message(STATUS "  Loop unrolling    : ${CLOTHOIDS_UNROLL_LOOPS}")
message(STATUS "  Fast math         : ${CLOTHOIDS_ENABLE_FAST_MATH}")

set(CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:-O3" "RelWithDebInfo:-O2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:/O2" "RelWithDebInfo:/O2")
endif()
if(CLOTHOIDS_USE_NATIVE_ARCH)
  if(CLOTHOIDS_APPLE_MULTI_ARCH)
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "native-arch:disabled-for-Apple-universal-build")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:-march=native" "Release:-mtune=native")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:-xHost")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:/arch:AVX2")
  endif()
endif()
if(CLOTHOIDS_UNROLL_LOOPS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:-funroll-loops")
endif()
if(CLOTHOIDS_ENABLE_FAST_MATH)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:-ffast-math")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:/fp:fast")
  endif()
endif()
if(CLOTHOIDS_IPO_SUPPORTED)
  list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release/RelWithDebInfo:IPO/LTO")
endif()
list(APPEND CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS "Release:NDEBUG" "Release:EIGEN_NO_DEBUG")
list(JOIN CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS " " CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS_TEXT)
message(STATUS "  Optimization flags: ${CLOTHOIDS_EFFECTIVE_OPTIMIZATION_FLAGS_TEXT}")
message(STATUS "")
