file(GLOB MY_PUBLIC_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.qml
  )
file(COPY ${MY_PUBLIC_HEADERS}  DESTINATION $ENV{HOME}/.qgis2/processing/scripts/)
  




