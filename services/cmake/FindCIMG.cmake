FIND_PATH(CIMG_INCLUDE_DIR CImg.h 
  ${EXTERN_PROJECT_SRC_DIR}/cimg/
)


IF (CIMG_INCLUDE_DIR)
  SET(CIMG_INCLUDE_FOUND "YES")
  MESSAGE(STATUS "Found CIMG include dir: ${CIMG_INCLUDE_DIR}")
ELSE(CIMG_INCLUDE_DIR)
 MESSAGE(STATUS "/!\ /!\ Not Found : CIMG include dir : ${CIMG_INCLUDE_DIR}")
ENDIF (CIMG_INCLUDE_DIR)

IF(CIMG_INCLUDE_FOUND)
  SET(CIMG_FOUND "YES")
ENDIF(CIMG_INCLUDE_FOUND)