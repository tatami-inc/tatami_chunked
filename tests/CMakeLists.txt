include(FetchContent)
FetchContent_Declare(
    tatami_test
    GIT_REPOSITORY https://github.com/tatami-inc/tatami_test
    GIT_TAG master
)
FetchContent_MakeAvailable(tatami_test)

add_executable(
    libtest 
    src/LruSlabCache.cpp
    src/OracularSlabCache.cpp
    src/OracularVariableSlabCache.cpp
    src/OracularSubsettedSlabCache.cpp
    src/ChunkDimensionStats.cpp
    src/SlabCacheStats.cpp
    src/CustomDenseChunkedMatrix.cpp
    src/CustomSparseChunkedMatrix.cpp
)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

target_compile_options(libtest PRIVATE -Wall -Wextra -Wpedantic -Werror)

target_link_libraries(libtest tatami_chunked tatami_test)

# Making the tests discoverable.
include(GoogleTest)
gtest_discover_tests(libtest)
