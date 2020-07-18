file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable(HitNllsTest ${SOURCES})
target_link_libraries(
    HitNllsTest
    ${PROJECT_NAME}
)