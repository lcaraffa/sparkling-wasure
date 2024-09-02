cmake_minimum_required(VERSION 2.8)

project (qpbo)
# src CMakeLists.txt
set (QPBO_SRC_DIR src)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib/)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_BUILD_DIR}/lib/)
SET(INCLUDE_DIR ${qpbo_SOURCE_DIR}/include)

SET(GLOBAL_BUILD_DIR "$ENV{GLOBAL_BUILD_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -std=c++0x")

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

INCLUDE_DIRECTORIES(
    ${INCLUDE_DIR}
    )

ADD_LIBRARY(qpbo STATIC
         ${QPBO_SRC_DIR}/QPBO.cpp
	 ${QPBO_SRC_DIR}/QPBO_extra.cpp
	 ${QPBO_SRC_DIR}/QPBO_maxflow.cpp
	 ${QPBO_SRC_DIR}/QPBO_postprocessing.cpp
       )

#TARGET_LINK_LIBRARIES(opt ${LIBRARIES_SOFTWARE} )

#INSTALL(CODE "MESSAGE(\"Installing headers...\")")
#INSTALL(FILES ${QPBO_SRC_DIR}/QPBO.h ${QPBO_SRC_DIR}/block.h DESTINATION ${INCLUDE_DIR})



