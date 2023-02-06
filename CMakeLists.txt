# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( compute-LOD2 )

set (CMAKE_CXX_STANDARD 17)

# CGAL and its components
find_package( CGAL 5.4.1 QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# GDAL and its components
find_package( GDAL 3.5.0 CONFIG REQUIRED )

if ( NOT GDAL_FOUND )

  message(STATUS "This project requires the GDAL library, and will not be compiled.")

  return()  

endif()


# Ceres and its components
find_package(Ceres 2.1.0 REQUIRED)

if ( NOT Ceres_FOUND )

  message(STATUS "This project requires the Ceres library, and will not be compiled.")

  return()  

endif()


# include for local directory

# include for local package


# Creating entries for target: compute-LOD2
# ############################

add_executable( compute-LOD2  main.cpp raster.cpp code.cpp edge_collapse.cpp path.cpp bridge.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS compute-LOD2 )

# Link the executable to CGAL and third-party libraries
target_link_libraries(compute-LOD2 PRIVATE CGAL::CGAL GDAL::GDAL Ceres::ceres)
