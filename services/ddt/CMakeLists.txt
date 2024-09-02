# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
# Then modified by hand to add Eigen3.


project( wasure )


# ---- CMAKE stuffs ----
cmake_minimum_required(VERSION 3.5)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
cmake_policy(VERSION 2.8.4)
else()
cmake_policy(VERSION 2.6)
endif()
endif()

SET(GLOBAL_BUILD_DIR "$ENV{GLOBAL_BUILD_DIR}")
SET(EXTERN_PROJECT_SRC_DIR "$ENV{EXTERN_PROJECT_SRC_DIR}")
SET(DDT_MAIN_DIR "$ENV{DDT_MAIN_DIR}")
SET(CMAKE_FIND_DIR "$ENV{CMAKE_FIND_DIR}")

#add_definitions("-DDDT_CGAL_TRAITS_D=3")
#add_definitions("-DDDT_CGAL_TRAITS_D=2")

if(DDT_TRAITS STREQUAL "2")
  add_definitions("-DDDT_CGAL_TRAITS_2=1")
  message(STATUS "Using 2D cgal traits")
else()
  if(DDT_TRAITS STREQUAL "3")
    add_definitions("-DDDT_CGAL_TRAITS_3=1")
    message(STATUS "Using 3D cgal traits")
  else()
    string(SUBSTRING ${DDT_TRAITS} 1 -1 DDT_TRAITS_D)
    add_definitions("-DDDT_CGAL_TRAITS_D=${DDT_TRAITS_D}")
    message(STATUS "Using dD cgal traits (d=${DDT_TRAITS_D})")
  endif()
endif()


# ---- Add /Defines depts ----
SET(INCLUDE_DIR ${wasure_SOURCE_DIR}/include/)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/bin)


# Includes contrib
#INCLUDE(${CMAKE_FIND_DIR}/FindCRITERION.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindANN.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindDDTLIB.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindEIGEN.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindGRAPHCUT.cmake)
#INCLUDE(${CMAKE_FIND_DIR}/FindPUGIXML.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindTINYPLY.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindCIMG.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindHDFS.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindDouble-conversion.cmake)
#INCLUDE(${CMAKE_FIND_DIR}/FindCUT-PURSUIT.cmake)


add_definitions("-DDDT_USE_HDFS=1")



find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

find_package(CGAL QUIET)
if ( CGAL_FOUND )
  set(CGAL_OPTIMISATION_NO_PRECONDITIONS TRUE)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Don't override flags")
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
else()
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
endif()

find_package(Eigen3)
if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )	
else()
    message(STATUS "NOTICE: This program requires the Eigen3 library, and will not be compiled.")
endif()


set(WASURE_EXE_LIB_DEP
  ${DDTLIB_LIBRARIES}
  ${ANN_LIBRARIES}
  ${EIGEN_LIBRARIES}
  ${GRAPHCUT_LIBRARIES}
  ${TINYPLY_LIBRARIES}
  ${DOUBLE_CONV_LIBRARIES}
  ${CGAL_LIBRARIES} 
  ${CGAL_3RD_PARTY_LIBRARIES}
  ${CIMG_LIB}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)




# ---- Defin global flags and versions ----
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -Wuninitialized -Wall -O0 -g -fPIC -pthread -std=gnu++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -O3 -DNDEBUG -pthread  -fPIC -std=gnu++11")


set(CIMG_LIB -L/usr/X11R6/lib -lm -lpthread -lX11)


# ---- Add subprojects ------

add_subdirectory (${wasure_SOURCE_DIR}/src/exe/)




