project (tinyply)

cmake_minimum_required(VERSION 3.5)
#cmake_minimum_required(VERSION 3.12)
# src CMakeLists.txt

SET(GLOBAL_BUILD_DIR "$ENV{GLOBAL_BUILD_DIR}")
SET(EXTERN_PROJECT_SRC_DIR "$ENV{EXTERN_PROJECT_SRC_DIR}")

SET(TINYPLY_SRC_DIR ${tinyply_SOURCE_DIR}/source/)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib)
SET(INCLUDE_DIR ${tinyply_SOURCE_DIR}/include)

message(STATUS "include dir => ${INCLUDE_DIR}")
message(STATUS "libs dir => ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -lboost_thread -lboost_system -lboost_filesystem -fPIC -std=c++0x" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

include_directories (${INCLUDE_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


ADD_LIBRARY(tinyply STATIC
         ${TINYPLY_SRC_DIR}/tinyply.cpp
	   )

INSTALL(CODE "MESSAGE(\"Installing headers...\")")
