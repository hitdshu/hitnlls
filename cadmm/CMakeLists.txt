cmake_minimum_required(VERSION 3.1.3)
project(hitcadmm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 REQUIRED)

include_directories(
${PROJECT_SOURCE_DIR}
${EIGEN3_INCLUDE_DIR}
)

add_executable(Test test.cc)
