# tatami bindings for chunked matrices

![Unit tests](https://github.com/tatami-inc/tatami_chunked/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/tatami-inc/tatami_chunked/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/tatami-inc/tatami_chunked/branch/master/graph/badge.svg?token=Z189ORCLLR)](https://codecov.io/gh/tatami-inc/tatami_chunked)

## Overview

Consider the situation where a matrix is represented as a grid of equally-sized independent chunks,
where extraction of a particular matrix row/column requires retrieval of data from all overlapping chunks.
This is a common paradigm in random-access file formats like HDF5 where entire chunks are read from file at once.
It can also be used for in-memory data structures where the data for each chunk can be efficiently compressed to reduce memory usage,
e.g., the `ChunkedRleArraySeed` for run-length-encoding-compressed chunks from the [**DelayedArray**](https://bioconductor.org/packages/DelayedArray) R package.
The **tatami_chunked** library implements some common functionality for **tatami** extension developers to create their own chunked matrix representations.

## Features

Given a rectangular grid of chunks that make up a chunked matrix,
we define a "slab" as the set of chunks that overlap a single row or column (or some subset/contiguous block thereof).
We typically want to load and cache an entire slab at once, ensuring that future requests to adjacent row/columns can just use the cached values rather than re-reading or decompressing the same chunks.
The **tatami_chunked** library provides the `LruSlabCache` and `OracularSlabCache` classes to facilitate caching of the slabs in `tatami::Matrix` extractors.
The `TypicalSlabCacheWorkspace` class allows developers to easily switch between caching strategies, depending on whether an oracle is provided to predict the future access pattern.

The `CustomDenseChunkedMatrix` and `CustomSparseChunkedMatrix` classes implement the `tatami::Matrix` interface on top of a matrix of custom chunks.
These classes automatically perform slab caching given a set of options including the maximum cache size (see the `CustomDenseChunkedMatrixOptions` and `CustomSparseChunkedMatrixOptions` classes).
Developers can use this to quickly create matrix representations with arbitrary chunk compression schemes that can reduce the memory footprint, e.g., DEFLATE, run length encodings.
Obviously, this comes at the cost of speed whereby the chunks must be unpacked to extract the relevant data -
developers are expected to define an appropriate extraction method for dense/sparse chunks.

In simple cases, chunk extraction is "atomic", i.e., the entire chunk must be unpacked to extract a subset of data.
Developers can then use the `SimpleDenseChunkWrapper` and `SimpleSparseChunkWrapper` to wrap these simple chunks for use in the `Custom*ChunkedMatrix` classes.
These wrappers only need a method to inflate the entire chunk; they will automatically handle the extraction of the desired block/subset from each chunk.
(More advanced developers may prefer to write their own extraction methods that avoid inflating the entire chunk, in which case these wrappers are not necessary.)

Still confused?
Read the [documentation](https://tatami-inc.github.io/tatami_chunked).

## Building with CMake

### CMake using `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  tatami_chunked
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_chunked
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(tatami_chunked)
```

Then you can link to **tatami_chunked** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe tatami_chunked)

# For libaries
target_link_libraries(mylib INTERFACE tatami_chunked)
```

### CMake using `find_package()`

You can install the library by cloning a suitable version of this repository and running the following commands:

```sh
mkdir build && cd build
cmake .. -DTATAMI_CHUNKED_TESTS=OFF
cmake --build . --target install
```

Then you can use `find_package()` as usual:

```cmake
find_package(tatami_tatami_chunked CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE tatami::tatami_chunked)
```

### Manual

If you're not using CMake, the simple approach is to just copy the files the `include/` subdirectory -
either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This will also require the dependencies listed in the [`extern/CMakeLists.txt`](extern/CMakeLists.txt) file, namely the core [**tatami**](https://github.com/tatami-inc/tatami) library.
