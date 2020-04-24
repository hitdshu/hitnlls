cmake_minimum_required(VERSION 3.1.3)
project(hitnlls C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

include_directories(
${PROJECT_SOURCE_DIR}
)

file(GLOB TYPE_CC "${PROJECT_SOURCE_DIR}/type/*.cc")
file(GLOB NLLS_CC "${PROJECT_SOURCE_DIR}/nlls/*.cc")

add_library(${PROJECT_NAME} SHARED
${TYPE_CC}
${NLLS_CC}
)
target_link_libraries(${PROJECT_NAME})

