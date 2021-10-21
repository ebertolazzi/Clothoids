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

find_package(pybind11 REQUIRED
  HINTS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/pybind11/share/cmake)
find_package(PkgConfig)
find_package(Eigen3 3.3 NO_MODULE)

if (PkgConfig_FOUND)
  pkg_check_modules(ipopt ipopt)
endif()

include(./CMakeLists-Flags.txt)
add_subdirectory(.. Clothoids)

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

pybind11_add_module(${PROJECT_NAME} MODULE ${G2LIB_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
target_link_libraries(${PROJECT_NAME} PRIVATE ClothoidsStatic UtilsStatic QuarticStatic)
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
target_assign_flags(${PROJECT_NAME})

set(INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib)
install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_PREFIX})
