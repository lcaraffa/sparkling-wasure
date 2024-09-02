project (gc1)

cmake_minimum_required(VERSION 2.8)
# src CMakeLists.txt
set (OPTGC_SRC_DIR ./src)

SET(GLOBAL_BUILD_DIR "$ENV{GLOBAL_BUILD_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib/)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib/)

SET(INCLUDE_DIR ./include)
message(STATUS "include dir => ${YOLO}")

set(CMAKE_CXX_FLAGS "-O3 -lboost_thread -lboost_system -fPIC -lboost_filesystem")


include_directories (${INCLUDE_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

ADD_LIBRARY(gc1 STATIC
  
         ${OPTGC_SRC_DIR}/graph.cpp
	 ${OPTGC_SRC_DIR}/maxflow.cpp
       )

#TARGET_LINK_LIBRARIES(opt ${LIBRARIES_SOFTWARE} )
INSTALL(CODE "MESSAGE(\"Installing headers...\")")
#INSTALL(FILES ${OPTGC_SRC_DIR}/QPBO.h   DESTINATION ${INCLUDE_DIR})
