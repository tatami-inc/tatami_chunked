cmake_minimum_required(VERSION 3.24)

project(tatami_chunked
    VERSION 2.0.0
    DESCRIPTION "tatami bindings for chunked matrices"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami_chunked INTERFACE)
add_library(tatami::tatami_chunked ALIAS tatami_chunked)

option(TATAMI_CHUNKED_FETCH_EXTERN "Automatically fetch tatami_chunked's external dependencies." ON)
if(TATAMI_CHUNKED_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 3.0.0 CONFIG REQUIRED)
endif()

target_link_libraries(tatami_chunked INTERFACE tatami::tatami)

# Switch between include directories depending on whether the downstream is
# using the build directly or is using the installed package.
include(GNUInstallDirs)
target_include_directories(tatami_chunked 
    INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_chunked>"
)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(TATAMI_CHUNKED_TESTS "Build tatami_chunked's test suite." ON)
else()
    option(TATAMI_CHUNKED_TESTS "Build tatami_chunked's test suite." OFF)
endif()

if(TATAMI_CHUNKED_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Installing for find_package.
include(CMakePackageConfigHelpers)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_chunked)

install(TARGETS tatami_chunked 
    EXPORT tatami_chunkedTargets)

install(EXPORT tatami_chunkedTargets
    FILE tatami_tatami_chunkedTargets.cmake
    NAMESPACE tatami::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_chunked)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_chunkedConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_chunked)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_chunkedConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_chunkedConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_chunkedConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_chunked)


