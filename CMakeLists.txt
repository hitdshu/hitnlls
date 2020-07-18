cmake_minimum_required(VERSION 3.1.3)
project(hitnlls LANGUAGES CXX VERSION 2.1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -ffast-math")
set(CMAKE_BUILD_TYPE Release)
set(HITNLLS_HEADER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/hitnlls)

file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/*/*.cc")
add_library(${PROJECT_NAME} SHARED ${SOURCES})

install(
    TARGETS ${PROJECT_NAME}
    EXPORT hitnllsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

add_subdirectory(admm)
add_subdirectory(camera)
add_subdirectory(common)
add_subdirectory(geometry)
add_subdirectory(ils)
add_subdirectory(matrix)
add_subdirectory(type)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    hitnllsConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    EXPORT hitnllsTargets
    FILE hitnllsTargets.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
)

configure_file(${PROJECT_SOURCE_DIR}/cmake/hitnllsConfig.cmake.in hitnllsConfig.cmake @ONLY)
install(
    FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/hitnllsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/hitnllsConfigVersion.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME}
)
