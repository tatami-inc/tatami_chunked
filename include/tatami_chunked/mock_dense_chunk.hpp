#ifndef TATAMI_CHUNKED_MOCK_DENSE_CHUNK_HPP
#define TATAMI_CHUNKED_MOCK_DENSE_CHUNK_HPP

#include <vector>

/**
 * @file mock_dense_chunk.hpp
 * @brief Dense chunk interface to use in a `CustomDenseChunkedMatrix`.
 */

namespace tatami_chunked {

/**
 * @cond
 */
namespace MockDenseChunk_internal {

template<typename InflatedValue_>
using Workspace =  std::vector<InflatedValue_>; 

template<class Blob_>
struct Core {
    Core() = default;

    Core(Blob_ c) : chunk(std::move(c)) {}

private:
    Blob_ chunk;

    typedef typename Blob_::value_type value_type;

public:
    auto get_target_chunkdim(bool row) const {
        if (row) {
            return chunk.nrow();
        } else {
            return chunk.ncol();
        }
    }

    auto get_non_target_chunkdim(bool row) const {
        if (row) {
            return chunk.ncol();
        } else {
            return chunk.nrow();
        }
    }

public:
    template<typename Index_>
    void extract(bool row, Index_ target_start, Index_ target_length, Index_ non_target_start, Index_ non_target_length, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);
        output += static_cast<size_t>(target_start) * stride; // cast to size_t to avoid overflow.

        if (Blob_::row_major == row) {
            size_t non_target_chunkdim = get_non_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(target_start) * non_target_chunkdim + static_cast<size_t>(non_target_start);
            for (Index_ p = 0; p < target_length; ++p) {
                std::copy_n(srcptr, non_target_length, output);
                srcptr += non_target_chunkdim;
                output += stride;
            }

        } else {
            size_t target_chunkdim = get_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(non_target_start) * target_chunkdim + static_cast<size_t>(target_start);
            for (Index_ p = 0; p < target_length; ++p) {
                auto copy_srcptr = srcptr;
                auto copy_output = output;
                for (Index_ s = 0; s < non_target_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += target_chunkdim;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

    template<typename Index_>
    void extract(bool row, Index_ target_start, Index_ target_length, const std::vector<Index_>& non_target_indices, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);
        output += static_cast<size_t>(target_start) * stride; // cast to size_t to avoid overflow.

        if (Blob_::row_major == row) {
            size_t non_target_chunkdim = get_non_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(target_start) * non_target_chunkdim;
            for (Index_ p = 0; p < target_length; ++p) {
                auto copy_output = output;
                for (auto x : non_target_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
                srcptr += non_target_chunkdim;
                output += stride;
            }

        } else {
            size_t target_chunkdim = get_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(target_start);
            for (Index_ p = 0; p < target_length; ++p) {
                auto copy_output = output;
                for (size_t x : non_target_indices) {
                    *copy_output = srcptr[x * target_chunkdim];
                    ++copy_output;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

public:
    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& target_indices, Index_ non_target_start, Index_ non_target_length, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);

        if (Blob_::row_major == row) {
            size_t non_target_chunkdim = get_non_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            size_t offset = non_target_start;
            for (size_t p : target_indices) {
                auto srcptr = work.data() + p * non_target_chunkdim + offset;
                std::copy_n(srcptr, non_target_length, output + p * stride);
            }

        } else {
            size_t target_chunkdim = get_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(non_target_start) * target_chunkdim;
            for (size_t p : target_indices) {
                auto copy_srcptr = srcptr + p;
                auto copy_output = output + p * stride;
                for (Index_ s = 0; s < non_target_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += target_chunkdim;
                }
            }
        }
    }

    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& target_indices, const std::vector<Index_>& non_target_indices, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);

        if (Blob_::row_major == row) {
            size_t non_target_chunkdim = get_non_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            for (size_t p : target_indices) {
                auto srcptr = work.data() + p * non_target_chunkdim;
                auto copy_output = output + p * stride;
                for (auto x : non_target_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
            }

        } else {
            size_t target_chunkdim = get_target_chunkdim(row); // use size_t to avoid integer overflow with Index_.
            for (size_t p : target_indices) {
                auto srcptr = work.data() + p;
                auto copy_output = output + p * stride;
                for (size_t x : non_target_indices) {
                    *copy_output = srcptr[x * target_chunkdim];
                    ++copy_output;
                }
            }
        }
    }
};

struct MockBlob {
    MockBlob(std::vector<double> c, size_t nr, size_t nc) : data(std::move(c)), NR(nr), NC(nc) {}
    MockBlob() = default;

    typedef double value_type;

    static constexpr bool row_major = true;

private:
    std::vector<double> data;
    size_t NR;
    size_t NC;

public:
    size_t nrow() const {
        return NR;
    }

    size_t ncol() const {
        return NC;
    }

    void inflate(std::vector<double>& buffer) const {
        buffer.resize(data.size());
        std::copy(data.begin(), data.end(), buffer.begin());
    }
};

}
/**
 * @endcond
 */

/**
 * @brief Mock a simple dense chunk for a `CustomDenseChunkedMatrix`.
 *
 * Mock a simple dense chunk for use inside a `CustomDenseChunkedMatrix`.
 * Each chunk should represent a 2-dimensional array of numeric values.
 * The interface is "simple" as any extraction of data from the chunk retrieves the full extent of the target dimension, 
 * with no attempt at optimization if only a subset of dimension elements are of interest.
 */
struct MockSimpleDenseChunk {
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Temporary workspace for extracting data from the chunk.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations may use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        MockDenseChunk_internal::Workspace<value_type> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the target dimension.
     * This should be set to `false`, otherwise a `MockSubsettedDenseChunk` is expected.
     */
    static constexpr bool use_subset = false;

public:
    /**
     * @cond
     */
    // You can construct this however you like, I don't care.
    MockSimpleDenseChunk() = default;
    MockSimpleDenseChunk(std::vector<double> c, size_t nr, size_t nc) : core(MockDenseChunk_internal::MockBlob(std::move(c), nr, nc)) {}
    /**
     * @endcond
     */

private:
    MockDenseChunk_internal::Core<MockDenseChunk_internal::MockBlob> core;

public:
    /**
     * Extract all elements of the target dimension into an output buffer.
     * For each element, this method will extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param non_target_start Index of the start of the continguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * P`,
     * where `P` is the number of rows (if `row = true`) or columns (otherwise) in this chunk.
     * @param stride Distance between corresponding values from adjacent elements of the target dimension when they are being stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_length`.
     *
     * If `row = true`, we would extract data for all rows and a block of columns `[non_target_start, non_target_start + non_target_length)`;
     * conversely, if `row = false`, we would extract data for all columns and a block of rows.
     * For a target dimension index `p` and non-target dimension index `non_target_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * The `stride` option allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomDenseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, Index_ non_target_start, Index_ non_target_length, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<Index_>(row, 0, core.get_target_chunkdim(row), non_target_start, non_target_length, work.work, output, stride);
    }

    /**
     * Extract all elements of the target dimension into an output buffer.
     * For each element, this method will extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param non_target_indices Indexed subset of the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * P`,
     * where `P` is the number of rows (if `row = true`) or columns (otherwise) in this chunk.
     * @param stride Distance between corresponding values from consecutive target dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_indices.size()`.
     *
     * If `row = true`, we would extract all rows and a subset of columns in `non_target_indices`;
     * conversely, if `row = false`, we would extract data for all columns and a subset of rows.
     * For a target dimension index `p` and non-target dimension index `non_target_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * The `stride` option allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomDenseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& non_target_indices, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<Index_>(row, 0, core.get_target_chunkdim(row), non_target_indices, work.work, output, stride);
    }
};

/**
 * @brief Create a dense chunk for a `CustomDenseChunkedMatrix`.
 *
 * Wraps a dense blob in a simple chunk interface for use inside a `CustomDenseChunkedMatrix`.
 * Each blob should hold a (possibly compressed) 2-dimensional array of numeric values.
 * The wrapper satisfies the `MockSimpleDenseChunk` interface, but is even simpler;
 * extraction of any data involves realization of the entire blob, with no optimization for subsets of interest along either dimension.
 *
 * The `Blob_` class should provide the following:
 *
 * - `static constexpr bool row_major`, specifying whether the array is row-major.
 * - `typedef value_type`, specifying the type of the value in the array.
 * - A `nrow() const` method, defining the number of rows in the array.
 * - A `ncol() const` method, defining the number of columns in the array.
 * - A `void inflate(std::vector<value_type>& buffer) const` method that fills `buffer` with the contents of the array.
 *   This should be filled in row-major format if `row_major = true` and in column-major format otherwise.
 *
 * @tparam Blob_ Class to represent a simple chunk.
 */
template<class Blob_>
struct SimpleDenseChunkWrapper {
    /**
     * @cond
     */
    typedef typename Blob_::value_type value_type;

    typedef MockDenseChunk_internal::Workspace<value_type> Workspace;

    static constexpr bool use_subset = false;

    SimpleDenseChunkWrapper() = default;

    SimpleDenseChunkWrapper(Blob_ c) : core(std::move(c)) {}
    /**
     * @endcond
     */

private:
    MockDenseChunk_internal::Core<Blob_> core;

public:
    /**
     * @cond
     */
    template<typename Index_>
    void extract(bool row, Index_ non_target_start, Index_ non_target_length, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, 0, core.get_target_chunkdim(row), non_target_start, non_target_length, work, output, stride);
    }

    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& non_target_indices, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, 0, core.get_target_chunkdim(row), non_target_indices, work, output, stride);
    }
    /**
     * @endcond
     */
};

/**
 * @brief Mock a subsettable dense chunk for a `CustomDenseChunkedMatrix`.
 *
 * Mock a subsettable dense chunk for use inside a `CustomDenseChunkedMatrix`.
 * Each chunk should represent a (possible compressed) 2-dimensional array of numeric values.
 * The interface is smarter as it only extracts elements of interest along the target dimension.
 * The elements of interest may be either a contiguous block or a indexed subset,
 * as predicted for each chunk from the `OracularSubsettedSlabCache`.
 * This provides some opportunities for optimization if the chunk supports partial reads.
 */
struct MockSubsettedDenseChunk {
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Workspace for chunk extraction.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations may use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        MockDenseChunk_internal::Workspace<value_type> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the target dimension.
     * This should be set to `false`, otherwise a `MockSimpleDenseChunk` is expected.
     */
    static constexpr bool use_subset = true;

public:
    /**
     * @cond
     */
    MockSubsettedDenseChunk() = default;
    MockSubsettedDenseChunk(std::vector<double> c, size_t nr, size_t nc) : core(MockDenseChunk_internal::MockBlob(std::move(c), nr, nc)) {}
    /**
     * @endcond
     */

private:
    MockDenseChunk_internal::Core<MockDenseChunk_internal::MockBlob> core;

public:
    /**
     * Extract a contiguous block of the target dimension into an output buffer.
     * For each element, this method will only extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_start Index of the start of the contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the first row, otherwise it is the first column.
     * @param target_length Length of the contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param non_target_start Index of the start of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * target_length + non_target_length`.
     * @param stride Distance between corresponding values from consecutive target dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_length`.
     *
     * If `row = true`, we would extract a block of rows `[target_start, target_start + length)` and a block of columns `[non_target_start, non_target_start + non_target_length)`;
     * conversely, if `row = false`, we would extract a block of columns as the target and the block of rows as the secondary.
     * For a target dimension index `target_start + p` and non-target dimension index `non_target_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, Index_ target_start, Index_ target_length, Index_ non_target_start, Index_ non_target_length, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, target_start, target_length, non_target_start, non_target_length, work.work, output, stride);
    }

    /**
     * Extract a contiguous block of the target dimension into an output buffer.
     * For each element, this method will only extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_start Index of the first element on the target dimension to be extracted.
     * If `row = true`, this is the first row, otherwise it is the first column.
     * @param target_length Number of elements on the target dimension to be extracted.
     * If `row = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param non_target_indices Indices of the elements on the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (target_start + target_length)`.
     * @param stride Distance between corresponding values from consecutive target dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_indices.size()`.
     *
     * If `row = true`, we would extract a block of rows `[target_start, target_start + length)` and a subset of columns in `non_target_indices`;
     * conversely, if `row = false`, we would extract a block of columns as the target and a subset of rows instead.
     * For a target dimension index `p` and non-target dimension index `non_target_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, Index_ target_start, Index_ target_length, const std::vector<Index_>& non_target_indices, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, target_start, target_length, non_target_indices, work.work, output, stride);
    }

public:
    /**
     * Extract an indexed subset of the target dimension into an output buffer.
     * For each element, this method will only extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_indices Indices of the elements of the target dimension to be extracted.
     * If `row = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param non_target_start Index of the start of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (target_indices.back() + 1)`.
     * @param stride Distance between corresponding values from consecutive target dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_length`.
     *
     * If `row = true`, we would extract a subset of rows in `target_indices` and a block of columns `[non_target_start, non_target_start + non_target_length)`,
     * conversely, if `row = false`, we would extract a subset of columns as the target and the block of rows as the secondary.
     * For a target dimension index `p` and non-target dimension index `non_target_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& target_indices, Index_ non_target_start, Index_ non_target_length, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, target_indices, non_target_start, non_target_length, work.work, output, stride);
    }

    /**
     * Extract an indexed subset of the target dimension into an output buffer.
     * For each element, this method will only extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_indices Indices of the elements on the target dimension to be extracted.
     * If `row = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param non_target_indices Indices of the elements on the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (target_indices.back() + 1)`.
     * @param stride Distance between corresponding values from consecutive target dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `non_target_indices.size()`.
     *
     * If `row = true`, we would extract a subset of rows in `target_indices` and a subset columns in `non_target_indices`.
     * conversely, if `row = false`, we would extract a subset of columns as the target and the subset of rows as the secondary.
     * For a target dimension index `p` and non-target dimension index `non_target_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `row` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<typename Index_>
    void extract(bool row, const std::vector<Index_>& target_indices, const std::vector<Index_>& non_target_indices, Workspace& work, value_type* output, size_t stride) const {
        core.extract(row, target_indices, non_target_indices, work.work, output, stride);
    }
};

}

#endif
