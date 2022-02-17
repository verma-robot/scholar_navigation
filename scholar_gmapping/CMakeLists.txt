cmake_minimum_required(VERSION 3.0.2)
project(scholar_gmapping)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  nav_msgs
  nodelet
  scholar_openslam_gmapping
  roscpp
  rospy
  sensor_msgs
  std_msgs
  tf
)
find_package(Boost REQUIRED)

catkin_package()
include_directories(src ${catkin_INCLUDE_DIRS})


add_executable(scholar_gmapping src/scholar_gmapping.cpp src/main.cpp)
target_link_libraries(scholar_gmapping ${Boost_LIBRARIES} ${catkin_LIBRARIES})
add_dependencies(scholar_gmapping ${catkin_EXPORTED_TARGETS})
