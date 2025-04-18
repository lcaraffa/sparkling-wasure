FIND_PATH(HDFS_INCLUDE_DIR hdfs.h 
  ${EXTERN_PROJECT_SRC_DIR}/libhdfs/include/
)


FIND_LIBRARY(HDFS_LIBRARIES libhdfs.so 
  ${EXTERN_PROJECT_SRC_DIR}/libhdfs/lib/ 
)

  MESSAGE(STATUS "====> :   ${EXTERN_PROJECT_SRC_DIR}/libhdfs/source/ ")

FIND_LIBRARY(JVM_LIBRARIES
  NAMES libjvm.so
  HINTS  /usr/lib/jvm/java-8-openjdk-amd64/jre/lib/amd64/server/
)



IF (HDFS_INCLUDE_DIR)
  SET(HDFS_INCLUDE_FOUND "YES")
  MESSAGE(STATUS "Found HDFS include dir: ${HDFS_INCLUDE_DIR}")
ELSE(HDFS_INCLUDE_DIR)
 MESSAGE(STATUS "/!\ /!\ Not Found : HDFS include dir ")
ENDIF (HDFS_INCLUDE_DIR)

IF (HDFS_LIBRARIES)
  SET(HDFS_LIB_FOUND "YES")
  MESSAGE(STATUS "Found HDFS lib: ${HDFS_LIBRARIES}")
ELSE(HDFS_LIBRARIES) 
 MESSAGE(STATUS "/!\ /!\ Not Found : HDFS libs")
ENDIF (HDFS_LIBRARIES)

IF (JVM_LIBRARIES)
  SET(JVM_LIB_FOUND "YES")
  MESSAGE(STATUS "Found JVM lib: ${JVM_LIBRARIES}")
ELSE(JVM_LIBRARIES) 
 MESSAGE(STATUS "/!\ /!\ Not Found : JVM libs")
ENDIF (JVM_LIBRARIES)


IF(HDFS_LIB_FOUND AND HDFS_INCLUDE_FOUND)
  SET(HDFS_FOUND "YES")
ENDIF(HDFS_LIB_FOUND AND HDFS_INCLUDE_FOUND)
