cmake_minimum_required(VERSION 2.8.3)
project(integrated_nav)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  cv_bridge
  geometry_msgs
  image_transport
  message_generation
  nav_msgs
  roscpp
  rospy
  sensor_msgs
  std_msgs
  tf
)

find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Sophus REQUIRED)
FIND_PACKAGE(fast REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system)

include_directories(
	include
	${catkin_INCLUDE_DIRS} 
	${CERES_INCLUDE_DIRS}
	${OpenCV_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
	${Sophus_INCLUDE_DIRS}
	${fast_INCLUDE_DIRS}
)

list(APPEND SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/global.cpp
    ${PROJECT_SOURCE_DIR}/src/fusion_handler.cpp
    ${PROJECT_SOURCE_DIR}/src/ins_gnss_simulator/ins_gnss_simulator_node.cpp
    ${PROJECT_SOURCE_DIR}/src/ins_gnss_simulator/simulator.cpp
)

# Set link libraries
list(APPEND LINK_LIBS 
  ${OpenCV_LIBS}
  ${Sophus_LIBRARIES}
  ${Boost_LIBRARIES}
  ${catkin_LIBRARIES}
  ${CERES_LIBRARIES}
  ${fast_LIBRARIES}
)  

## Generate messages in the 'msg' folder
add_message_files(
  FILES
  GNSS_ECEF.msg
)


## Generate added messages and services with any dependencies listed here
generate_messages(
  DEPENDENCIES
  std_msgs
)

catkin_package(
  CATKIN_DEPENDS
  message_runtime
  std_msgs
  geometry_msgs
)

catkin_package()

add_library(sensor_fusion SHARED ${SOURCE_FILES})
add_dependencies(sensor_fusion ${${PROJECT_NAME}_EXPORTED_TARGETS}})
target_link_libraries(sensor_fusion ${LINK_LIBS})

add_executable(main_node src/main_node.cpp)
target_link_libraries(main_node sensor_fusion)

add_executable(ins_gnss_simulator_node src/ins_gnss_simulator/ins_gnss_simulator_node.cpp)
target_link_libraries(ins_gnss_simulator_node sensor_fusion)

