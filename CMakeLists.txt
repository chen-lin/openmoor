cmake_minimum_required(VERSION 3.9)

if(APPLE)
	set(CMAKE_C_COMPILER "gcc")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_CXX_STANDARD 11)
endif()

# set(CMAKE_CXX_FLAGS "$CMAKE_CXX_FLAGS -std=c++11 -stdlib=libstdc++")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(OpenMOOR)

# The version number
set(OpenMOOR_VERSION_MAJOR 0)
set(OpenMOOR_VERSION_MINOR 1)

# Compile dynamic linked library or executables: comment the followed line if 
# compile executables
set(OpenMOOR_API 1)

# set OpenMP for parallel computing
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Configure a header file to pass some of the CMake settings to the source code
configure_file(
	"${PROJECT_SOURCE_DIR}/moorconfig.h.in"
	"${PROJECT_BINARY_DIR}/moorconfig.h"
	)
	
# Add the binary tree to the search path for include file
# So that we will find moorconfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")

# Add_subdirectory (include)
# set (EXTRA_LIBS ${EXTRA_LIBS} include)

# Add the excutable from the src folder
if (OpenMOOR_API)
	add_library(MoorApi SHARED
		src/cable.cc
		src/catenary.cc
		src/input.cc
		src/main.cc
		src/meshoptimizer.cc
		src/moorerror.cc
		src/node.cc
		src/platform.cc
		src/reader.cc
		src/setting.cc
		src/simulation.cc
		src/solver.cc
		src/utility.cc
		src/writer.cc
		src/moorapi.cc
		include/Eigen/Dense
		include/Eigen/Polynomials
		include/Eigen/Geometry
		include/Eigen/QR
		include/boost/numeric/odeint.hpp
		include/rapidxml-1.13/rapidxml_print.hpp
		include/rapidxml-1.13/rapidxml_utils.hpp
		include/rapidxml-1.13/rapidxml.hpp
		)
# Add the excutable from the src folder		
else()
	add_executable(OpenMOOR 
		src/cable.cc
		src/catenary.cc
		src/input.cc
		src/main.cc
		src/meshoptimizer.cc
		src/moorerror.cc
		src/node.cc
		src/platform.cc
		src/reader.cc
		src/setting.cc
		src/simulation.cc
		src/solver.cc
		src/utility.cc
		src/writer.cc
		include/Eigen/Dense
		include/Eigen/Polynomials
		include/Eigen/Geometry
		include/Eigen/QR
		include/boost/numeric/odeint.hpp
		include/rapidxml-1.13/rapidxml_print.hpp
		include/rapidxml-1.13/rapidxml_utils.hpp
		include/rapidxml-1.13/rapidxml.hpp
		)
endif()