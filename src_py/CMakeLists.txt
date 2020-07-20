##
# PYTHON Wrapper for Clothoids
# 
# License MIT - See LICENSE file
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
#      Enrico Bertolazzi, Marco Frego
#

cmake_minimum_required(VERSION 3.10.2)

project(G2lib VERSION 1.0)

find_package(pybind11 REQUIRED)

find_package(pybind11)

	set(CLOTHOIDS "${CMAKE_CURRENT_SOURCE_DIR}/../lib")

	set(SOURCES
	  "src/python-Triangle2D.cc"
		"src/python-LineSegment.cc"
		"src/python-G2lib.cc"
		"src/python-Clothoid.cc"
		"src/python-CircleArc.cc"
		"src/python-Biarc.cc"
		"src/python-BaseCurve.cc"
		"src/python-AABBtree.cc")

pybind11_add_module(${PROJECT_NAME} MODULE ${SOURCES} NO_EXTRAS)

	target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

	target_include_directories(${PROJECT_NAME}
			PUBLIC "${CLOTHOIDS}/include"
			PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

	set_property(TARGET ${PROJECT_NAME} PROPERTY
			POSITION_INDEPENDENT_CODE 1)

if(APPLE)
  target_link_libraries(${PROJECT_NAME} 
    PUBLIC "${CLOTHOIDS}/lib/libClothoids_osx_static.a")
ELSEIF(LINUX OR UNIX)
  target_link_libraries(${PROJECT_NAME} 
    PUBLIC "${CLOTHOIDS}/lib/libClothoids_linux_static.a")
else()
  target_link_libraries(${PROJECT_NAME} 
    PUBLIC "${CLOTHOIDS}/lib/Clothoids_win__static.lib")
endif()

set(INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib)
install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_PREFIX})
