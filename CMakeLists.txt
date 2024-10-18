# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( compute-LOD2 )

set (CMAKE_CXX_STANDARD 17)

#add_compile_options(-O0)

if(NOT CMAKE_BUILD_TYPE MATCHES "Release")
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

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


# Eigen3 and its components
find_package(Eigen3 3.4.0 REQUIRED)

if ( NOT Eigen3_FOUND )

  message(STATUS "This project requires the Eigen3 library, and will not be compiled.")

  return()  

endif()


# include for local directory

# include for local package

add_library (EdgeCollapse edge_collapse.cpp)
target_link_libraries(EdgeCollapse PRIVATE CGAL::CGAL Eigen3::Eigen)

# Creating entries for target: compute-LOD2
# ############################

add_executable( compute-LOD2  main.cpp raster.cpp code.cpp path.cpp bridge.cpp)

target_compile_options(compute-LOD2 PRIVATE -Wall -Wextra -Wpedantic)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS compute-LOD2 )

# Link the executable to CGAL and third-party libraries
target_link_libraries(compute-LOD2 PRIVATE CGAL::CGAL GDAL::GDAL Ceres::ceres Eigen3::Eigen EdgeCollapse)

# Creating entries for target: do-edge-collapse
# ############################

add_executable( do-edge-collapse  main_edge_collapse.cpp)

target_compile_options(do-edge-collapse PRIVATE -Wall -Wextra -Wpedantic)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS do-edge-collapse )

# Link the executable to CGAL and third-party libraries
target_link_libraries(do-edge-collapse PRIVATE CGAL::CGAL GDAL::GDAL Eigen3::Eigen EdgeCollapse)
