cmake_minimum_required(VERSION 3.1)

project(pgo4imu)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

find_package(Eigen3 REQUIRED)

set (CMAKE_CXX_STANDARD 17)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)

link_directories(
)

add_executable(attopt
    src/attopt.cpp
)

target_link_libraries(attopt
    Eigen3::Eigen
)

add_executable(posopt
    src/posopt.cpp
)

target_link_libraries(posopt
    Eigen3::Eigen
)

add_executable(se3test
    src/se3test.cpp
)

target_link_libraries(se3test
    Eigen3::Eigen
)

add_executable(monocam_loop_closure_test
    src/monocam_loop_closure_test.cpp
    src/imutrj_optimizer.cpp
    include/imutrj_optimizer.hpp
)

target_link_libraries(monocam_loop_closure_test
    Eigen3::Eigen
)

add_executable(pointcloud_registration_test
        src/pointcloud_registration_test.cpp
)

target_link_libraries(pointcloud_registration_test
        Eigen3::Eigen
)

add_executable(foot_mounted_test
        src/imutrj_optimizer.cpp
        src/foot_mounted_test.cpp
)

