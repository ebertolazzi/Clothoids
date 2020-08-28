##
# PYTHON Wrapper for Clothoids
# 
# License MIT - See LICENSE file
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
#      Enrico Bertolazzi, Marco Frego
#

cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

project(G2lib VERSION 1.0)

find_package(pybind11 REQUIRED)
find_package(PkgConfig)
find_package(Eigen3 3.3 NO_MODULE)

if (PkgConfig_FOUND)
  pkg_check_modules(ipopt ipopt)
endif()

set(G2LIB_SOURCES
	"src/python-Triangle2D.cc"
	"src/python-LineSegment.cc"
	"src/python-G2lib.cc"
	"src/python-Clothoid.cc"
	"src/python-CircleArc.cc"
	"src/python-Biarc.cc"
	"src/python-BaseCurve.cc"
	"src/python-AABBtree.cc"
  "src/python-ClothoidSpline-Solver.cc"
  "src/python-ClothoidSpline-LMSolver.cc"
  "src/python-ClothoidSpline-IpoptSolver.cc"
  "src/python-ClothoidSpline-Interpolation.cc")

include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists-QuarticRootsFlocke.txt")
include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists-Clothoids.txt")
pybind11_add_module(${PROJECT_NAME} MODULE ${G2LIB_SOURCES} NO_EXTRAS)

target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
target_link_libraries(${PROJECT_NAME} PRIVATE quartic_roots_flocke)
target_link_libraries(${PROJECT_NAME} PRIVATE clothoids)
if (ipopt_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE IPOPT_CLOTHOID_SPLINE)
  target_include_directories(${PROJECT_NAME} PRIVATE ${ipopt_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${ipopt_LIBRARIES})
  target_compile_definitions(${PROJECT_NAME} PRIVATE HAVE_CSTDDEF)
endif()
if (TARGET Eigen3::Eigen)
  target_compile_definitions(${PROJECT_NAME} PRIVATE LMSOLVE_CLOTHOID_SPLINE)
  target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

  target_compile_options(${PROJECT_NAME} PRIVATE 
    -std=c++14 -Wall -Wno-float-equal -Wno-padded
    $<$<CONFIG:Release>: -funroll-loops -O2 -g0 >
    $<$<CONFIG:Debug>: -O0 -g3 >)

  target_compile_definitions(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Debug>: MECHATRONIX_DEBUG DEBUG >)

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  target_compile_options(${PROJECT_NAME} PRIVATE
    -std=c++14 -stdlib=libc++  -Weverything -Wno-weak-vtables 
    -Wno-implicit-fallthrough -Wno-documentation-unknown-command 
    -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic 
    -Wno-padded -Wno-reserved-id-macro -Wno-documentation
    $<$<CONFIG:Release>: -funroll-loops -O2 -g0 >
    $<$<CONFIG:Debug>: -O0 -g3 >)

  target_compile_definitions(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Debug>: " MECHATRONIX_DEBUG DEBUG " >)

elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

  target_compile_options(${PROJECT_NAME} PRIVATE
    /nologo /GS /W3 /WX- /Gm- /fp:precise /EHsc 
    $<$<CONFIG:Release>: /O2 /MD >
    $<$<CONFIG:Debug>: /Od /Ob0 /MDd /Z7 /RTC1 >)

  target_compile_definitions(${PROJECT_NAME} PRIVATE
    _WINDOWS _SCL_SECURE_NO_WARNINGS HAVE_STRING_H NO_GETTIMEOFDAY
    $<$<CONFIG:Debug>: _DEBUG MECHATRONIX_DEBUG >)

else()
  
	message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")

endif()


set(INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib)
install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_PREFIX})
