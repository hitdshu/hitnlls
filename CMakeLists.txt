cmake_minimum_required(VERSION 3.1.3)
project(hitnlls C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

include_directories(
${PROJECT_SOURCE_DIR}
)

file(GLOB COMMON_CC "${PROJECT_SOURCE_DIR}/common/*.cc")
file(GLOB MATRIX_CC "${PROJECT_SOURCE_DIR}/matrix/*.cc")
file(GLOB NODE_CC "${PROJECT_SOURCE_DIR}/node/*.cc")
file(GLOB FACTOR_CC "${PROJECT_SOURCE_DIR}/factor/*.cc")
file(GLOB GRAPH_CC "${PROJECT_SOURCE_DIR}/graph/*.cc")
file(GLOB SOLVER_CC "${PROJECT_SOURCE_DIR}/solver/*.cc")

add_library(${PROJECT_NAME} SHARED
${COMMON_CC}
${MATRIX_CC}
${NODE_CC}
${FACTOR_CC}
${GRAPH_CC}
${SOLVER_CC}
)

target_link_libraries(${PROJECT_NAME}
)

file(GLOB TEST_CC "${PROJECT_SOURCE_DIR}/test/*.cc")

add_executable(RunAllTest
${TEST_CC})
target_link_libraries(RunAllTest
${PROJECT_NAME})