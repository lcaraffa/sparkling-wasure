cmake_minimum_required (VERSION 3.5)

set (CMAKE_CXX_STANDARD 20)


# --------- Includes -----------
include(${CGAL_USE_FILE})
include_directories(${INCLUDE_DIR})
include_directories(SYSTEM ${ANN_INCLUDE_DIR})
include_directories(SYSTEM ${GRAPHCUT_INCLUDE_DIR})
include_directories(SYSTEM ${CIMG_INCLUDE_DIR})
include_directories(SYSTEM ${TINYPLY_INCLUDE_DIR})
include_directories(SYSTEM ${PUGIXML_INCLUDE_DIR})
include_directories(SYSTEM ${DOUBLE_CONV_INCLUDE_DIR})
include_directories(${DDTLIB_INCLUDE_DIR})
include_directories(${CUT-PURSUIT_INCLUDE_DIR})

include_directories(../../../ddt/include)

include_directories(${DDT_MAIN_DIR}/include/)
include_directories(../../../../build-spark-Release/)


# ---------- Make lib ----------
set(wasureLib_SOURCES
wasure_algo.cpp
simplex_data_wasure.cpp
simplex_data.cpp
maths.cpp
)
add_library(tbmrf SHARED ${wasureLib_SOURCES} )

target_link_libraries(tbmrf
  ${WASURE_EXE_LIB_DEP}
  ${CIMG_LIB}
  )


