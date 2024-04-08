include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/afd902e992b720d1b3e106bc5e425a5768872265.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Avoid installing GoogleTest when installing this project.
option(INSTALL_GTEST "Enable installation of googletest." OFF)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    libtest 
    src/LruSlabCache.cpp
    src/OracleSlabCache.cpp
    src/SubsettedOracleSlabCache.cpp
    src/mock_dense_chunk.cpp
    src/mock_sparse_chunk.cpp
#    src/CustomDenseChunkedMatrix.cpp
#    src/CustomSparseChunkedMatrix.cpp
)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

target_compile_options(libtest PRIVATE -Wall -Wextra -Wpedantic -Werror)

target_link_libraries(libtest tatami_chunked gtest_main)

# Making the tests discoverable.
include(GoogleTest)
gtest_discover_tests(libtest)
