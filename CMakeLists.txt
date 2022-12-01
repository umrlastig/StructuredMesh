# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( compute-LOD2 )

set (CMAKE_CXX_STANDARD 17)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# GDAL and its components
find_package( GDAL CONFIG REQUIRED )

if ( NOT GDAL_FOUND )

  message(STATUS "This project requires the GDAL library, and will not be compiled.")

  return()  

endif()


# Eigen and its components
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if ( NOT Eigen3_FOUND )

  message(STATUS "This project requires the Eigen3 library, and will not be compiled.")

  return()  

endif()


# include for local directory

# include for local package


# Creating entries for target: compute-LOD2
# ############################

add_executable( compute-LOD2  main.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS compute-LOD2 )

# Link the executable to CGAL and third-party libraries
target_link_libraries(compute-LOD2 PRIVATE CGAL::CGAL GDAL::GDAL Eigen3::Eigen)
