include(cmake_params.txt)
project( ddt )


 # ---- CMAKE configuration -------
cmake_minimum_required(VERSION 3.5)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set (CMAKE_CXX_STANDARD 20)

set(CPP_CORE_LOCAL_DIR "src/core/")
set(SERVICE_LOCAL_DIR "services/")
SET(EXTERN_PROJECT_SRC_DIR "${ddt_SOURCE_DIR}/services/extern/")
SET(CMAKE_FIND_DIR "${ddt_SOURCE_DIR}/services/cmake/")


# ---- Add /Defines depts -----
# Local includ dir
SET(INCLUDE_DIR
  ${ddt_SOURCE_DIR}/src/core/
  ${ddt_SOURCE_DIR}/build/conf_header/
  )

# ---- Global system dependencies : Booost Egein CGAL ----
# BOOST
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
# EIGEN
find_package(Eigen3)
if(EIGEN3_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCGAL_EIGEN3_ENABLED")
  if(EIGEN3_USE_FILE)
    include(${EIGEN3_USE_FILE})
  endif()
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message(STATUS "NOTICE: Eigen library is not found.")
endif()
# CGAL
find_package(CGAL QUIET)
if ( CGAL_FOUND )
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Don't override flags")
  set(CGAL_OPTIMISATION_NO_PRECONDITIONS TRUE)
  include( ${CGAL_USE_FILE} )
  message(STATUS "CGAL version : ${CGAL_VERSION}")
else()
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
endif()
find_package( TBB QUIET )
if( TBB_FOUND )
  if(NOT COMMAND CGAL_target_use_TBB)
    include(${TBB_USE_FILE})
    add_definitions("-DTBB_IMPLEMENT_CPP0X=0")
    macro(CGAL_target_use_TBB target)
    endmacro()
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()
endif()

# ----- Thirs part libs hardcoded -------
INCLUDE(${CMAKE_FIND_DIR}/FindHDFS.cmake)
INCLUDE(${CMAKE_FIND_DIR}/FindDouble-conversion.cmake)


# ---- Defin global flags and versions ----
#set(DDT_TRAITS "3" CACHE STRING "DDT Traits (2,3,D2,D3...)")
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

add_definitions("-Wunused -O3")
if(DDT_DEBUG)
  add_definitions("-DDDT_DEBUG=1")
endif()
if(DDT_USE_THREADS)
  add_definitions("-DDDT_USE_THREADS=1")
endif()
# if(DDT_USE_HDFS)
#   add_definitions("-DDDT_USE_HDFS=1")
# endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0  -std=gnu++11")



# --------- Start Lib dependances
include_directories(${INCLUDE_DIR})
#include_directories(${HDFS_INCLUDE_DIR})
include_directories(${DOUBLE_CONV_INCLUDE_DIR})
include(${CGAL_USE_FILE})

# ---------- Make libs ----------
set(libddt_SOURCES
  ${CPP_CORE_LOCAL_DIR}/algorithm.hpp
  ${CPP_CORE_LOCAL_DIR}/DDT.hpp
  ${CPP_CORE_LOCAL_DIR}/bbox.hpp
  ${CPP_CORE_LOCAL_DIR}/tile.hpp

  ${CPP_CORE_LOCAL_DIR}/io/logging.cpp
  ${CPP_CORE_LOCAL_DIR}/io/logging_stream.cpp
  ${CPP_CORE_LOCAL_DIR}/io/logging.hpp
  ${CPP_CORE_LOCAL_DIR}/io/logging_stream.hpp
  ${CPP_CORE_LOCAL_DIR}/io/read.hpp
  ${CPP_CORE_LOCAL_DIR}/io/read_stream.hpp
  ${CPP_CORE_LOCAL_DIR}/io/stream_api.hpp
  ${CPP_CORE_LOCAL_DIR}/io/write.hpp
  ${CPP_CORE_LOCAL_DIR}/io/write_vrt.hpp
  ${CPP_CORE_LOCAL_DIR}/io/write_vrt.cpp
  ${CPP_CORE_LOCAL_DIR}/io/write_ply.cpp
  ${CPP_CORE_LOCAL_DIR}/io/write_ply.hpp
  ${CPP_CORE_LOCAL_DIR}/io/write_stream.hpp
  ${CPP_CORE_LOCAL_DIR}/io/number_parser.hpp
  ${CPP_CORE_LOCAL_DIR}/io/number_parser.cpp
  ${CPP_CORE_LOCAL_DIR}/io/base64.hpp
  ${CPP_CORE_LOCAL_DIR}/io/base64.cpp

  ${CPP_CORE_LOCAL_DIR}/iterator/Cell_const_iterator.hpp
  ${CPP_CORE_LOCAL_DIR}/iterator/Facet_const_iterator.hpp
  ${CPP_CORE_LOCAL_DIR}/iterator/Map_iterators.hpp
  ${CPP_CORE_LOCAL_DIR}/iterator/Vertex_const_iterator.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/Facet_const_iterator_d.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/Facet_const_iterator_2.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/Facet_const_iterator_3.hpp

  ${CPP_CORE_LOCAL_DIR}/scheduler/thread_pool.hpp
  ${CPP_CORE_LOCAL_DIR}/scheduler/safe.hpp
  ${CPP_CORE_LOCAL_DIR}/scheduler/multithread_scheduler.hpp
  ${CPP_CORE_LOCAL_DIR}/scheduler/sequential_scheduler.hpp
  ${CPP_CORE_LOCAL_DIR}/scheduler/scheduler.hpp

  ${CPP_CORE_LOCAL_DIR}/partitioner/random_partitioner.hpp
  ${CPP_CORE_LOCAL_DIR}/partitioner/grid_partitioner.hpp
  ${CPP_CORE_LOCAL_DIR}/partitioner/const_partitioner.hpp

  ${CPP_CORE_LOCAL_DIR}/traits/cgal_traits_2.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/cgal_traits_3.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/cgal_traits_base.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/cgal_traits_d.hpp
  ${CPP_CORE_LOCAL_DIR}/traits/traits.hpp
  )

add_library( libddt SHARED ${libddt_SOURCES} )
set_property( TARGET libddt PROPERTY CXX_STANDARD 20)  
set_target_properties( libddt PROPERTIES PREFIX "" )


configure_file(${CPP_CORE_LOCAL_DIR}/conf/conf.hpp.in conf_header/conf.hpp)


# Add a custom target
add_custom_target("format" COMMAND
  astyle
  --style=allman --keep-one-line-blocks --recursive  --suffix=none
  ${CMAKE_CURRENT_SOURCE_DIR}/${CPP_CORE_LOCAL_DIR}/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${CPP_CORE_LOCAL_DIR}/*.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${SERVICES_LOCAL_DIR}/ddt/src/exe/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${SERVICES_LOCAL_DIR}/ddt/include/*.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${SERVICES_LOCAL_DIR}/wasure/src/exe/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${SERVICES_LOCAL_DIR}/wasure/src/lib/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/${SERVICES_LOCAL_DIR}/wasure/include/*.hpp
  VERBATIM
  )


