
cmake_minimum_required (VERSION 3.5)

# --------- Includes -----------
include(${CGAL_USE_FILE})
include_directories(${INCLUDE_DIR})
include_directories(SYSTEM ${TINYPLY_INCLUDE_DIR})
#include_directories(SYSTEM ${PUGIXML_INCLUDE_DIR})
include_directories(SYSTEM ${DDTLIB_INCLUDE_DIR})
include_directories(SYSTEM ${ANN_INCLUDE_DIR})
include_directories(SYSTEM ${GRAPHCUT_INCLUDE_DIR})
include_directories(SYSTEM ${CIMG_INCLUDE_DIR})
include_directories(SYSTEM ${GLOBAL_BUILD_DIR})
include_directories(SYSTEM ${HDFS_INCLUDE_DIR})
include_directories(SYSTEM ${DOUBLE_CONV_INCLUDE_DIR})
find_package(ZLIB)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(ddt_stream_SOURCES
  main_ddt_stream.cpp
  )
add_executable(ddt-stream-exe ${ddt_stream_SOURCES})
target_compile_options(ddt-stream-exe PRIVATE -Wall -Wextra )
set_property(TARGET ddt-stream-exe PROPERTY CXX_STANDARD 11)
target_link_libraries(ddt-stream-exe
  # intern dependencies
  ${WASURE_EXE_LIB_DEP}
  ${HDFS_LIBRARIES}
  ${DOUBLE_CONV_LIBRARIES}
  )



set(ddt_stream_SOURCES
  identity.cpp
  )
add_executable(identity-exe ${ddt_stream_SOURCES})
set_property(TARGET identity-exe PROPERTY CXX_STANDARD 11)
target_link_libraries(identity-exe
  # intern dependencies
  ${WASURE_EXE_LIB_DEP}
  ${HDFS_LIBRARIES}
  )
 
if(CMAKE_VERBOSITY GREATER 0)
message(STATUS "Added target wasure (exe)")
endif(CMAKE_VERBOSITY GREATER 0)
