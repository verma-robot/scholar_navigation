cmake_minimum_required(VERSION 2.8.3)
project(scholar_openslam_gmapping)

find_package(catkin REQUIRED)

include(GenerateExportHeader)
set(EXPORT_HEADER_DIR "${CATKIN_DEVEL_PREFIX}/include")
file(MAKE_DIRECTORY "${EXPORT_HEADER_DIR}")

catkin_package(  INCLUDE_DIRS    include    ${EXPORT_HEADER_DIR}  LIBRARIES utils sensor_base sensor_odometry sensor_range   scanmatcher log gridfastslam)

include_directories(  include  ${EXPORT_HEADER_DIR})

add_library(utils
  utils/stat.cpp
  )

add_library(sensor_base
  sensor/sensor_base/sensor.cpp
  sensor/sensor_base/sensorreading.cpp)

add_library(sensor_odometry
  sensor/sensor_odometry/odometrysensor.cpp
  sensor/sensor_odometry/odometryreading.cpp)

target_link_libraries(sensor_odometry  sensor_base)

add_library(sensor_range
  sensor/sensor_range/rangesensor.cpp
  sensor/sensor_range/rangereading.cpp)

target_link_libraries(sensor_range  sensor_base)

add_library(log
  log/sensorlog.cpp
)
 
target_link_libraries(log  sensor_range sensor_odometry sensor_base )

add_library(scanmatcher
  scanmatcher/smmap.cpp
  scanmatcher/scholar_scanmatcher.cpp
  )


target_link_libraries(scanmatcher  log sensor_range sensor_odometry sensor_base utils )

add_library(gridfastslam
  gridfastslam/scholar_gridslamprocessor_tree.cpp
  gridfastslam/scholar_motionmodel.cpp
  gridfastslam/scholar_gridslamprocessor.cpp
  )


target_link_libraries(gridfastslam  scanmatcher log sensor_range sensor_odometry sensor_base utils )

install(TARGETS utils sensor_base sensor_odometry sensor_range  log scanmatcher gridfastslam 
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)

install(DIRECTORY include/gmapping
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h*"
  PATTERN ".svn" EXCLUDE
)

generate_export_header(utils
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/utils/utils_export.h)
generate_export_header(sensor_base
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/sensor/sensor_base/sensor_base_export.h)
generate_export_header(sensor_odometry
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/sensor/sensor_odometry/sensor_odometry_export.h)
generate_export_header(sensor_range
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/sensor/sensor_range/sensor_range_export.h)
generate_export_header(log
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/log/log_export.h)

generate_export_header(scanmatcher
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/scanmatcher/scanmatcher_export.h)
generate_export_header(gridfastslam
  EXPORT_FILE_NAME ${EXPORT_HEADER_DIR}/gmapping/gridfastslam/gridfastslam_export.h)
