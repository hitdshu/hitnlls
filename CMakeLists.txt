cmake_minimum_required(VERSION 3.1.3)
project(hitnlls C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

include_directories(
${PROJECT_SOURCE_DIR}
)

file(GLOB ADMM_CC "${PROJECT_SOURCE_DIR}/admm/*.cc")
file(GLOB CAMERA_CC "${PROJECT_SOURCE_DIR}/camera/*.cc")
file(GLOB COMMON_CC "${PROJECT_SOURCE_DIR}/common/*.cc")
file(GLOB ILS_CC "${PROJECT_SOURCE_DIR}/ils/*.cc")
file(GLOB TYPE_CC "${PROJECT_SOURCE_DIR}/type/*.cc")
add_library(${PROJECT_NAME} SHARED
${ADMM_CC}
${COMMON_CC}
${CAMERA_CC}
${ILS_CC}
${TYPE_CC}
)
target_link_libraries(${PROJECT_NAME}
)

