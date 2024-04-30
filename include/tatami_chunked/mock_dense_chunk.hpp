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
    template<bool accrow_>
    auto get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.nrow();
        } else {
            return chunk.ncol();
        }
    }

    template<bool accrow_>
    auto get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.ncol();
        } else {
            return chunk.nrow();
        }
    }

public:
    template<bool accrow_, typename Index_>
    void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);
        output += static_cast<size_t>(primary_start) * stride; // cast to size_t to avoid overflow.

        if constexpr(Blob_::row_major == accrow_) {
            size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(primary_start) * secondary_chunkdim + static_cast<size_t>(secondary_start);
            for (Index_ p = 0; p < primary_length; ++p) {
                std::copy_n(srcptr, secondary_length, output);
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(secondary_start) * primary_chunkdim + static_cast<size_t>(primary_start);
            for (Index_ p = 0; p < primary_length; ++p) {
                auto copy_srcptr = srcptr;
                auto copy_output = output;
                for (Index_ s = 0; s < secondary_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += primary_chunkdim;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

    template<bool accrow_, typename Index_>
    void extract(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);
        output += static_cast<size_t>(primary_start) * stride; // cast to size_t to avoid overflow.

        if constexpr(Blob_::row_major == accrow_) {
            size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(primary_start) * secondary_chunkdim;
            for (Index_ p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(primary_start);
            for (Index_ p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (size_t x : secondary_indices) {
                    *copy_output = srcptr[x * primary_chunkdim];
                    ++copy_output;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

public:
    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);

        if constexpr(Blob_::row_major == accrow_) {
            size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            size_t offset = secondary_start;
            for (size_t p : primary_indices) {
                auto srcptr = work.data() + p * secondary_chunkdim + offset;
                std::copy_n(srcptr, secondary_length, output + p * stride);
            }

        } else {
            size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            auto srcptr = work.data() + static_cast<size_t>(secondary_start) * primary_chunkdim;
            for (size_t p : primary_indices) {
                auto copy_srcptr = srcptr + p;
                auto copy_output = output + p * stride;
                for (Index_ s = 0; s < secondary_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += primary_chunkdim;
                }
            }
        }
    }

    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& primary_indices, const std::vector<Index_>& secondary_indices, Workspace<value_type>& work, value_type* output, size_t stride) const {
        chunk.inflate(work);

        if constexpr(Blob_::row_major == accrow_) {
            size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            for (size_t p : primary_indices) {
                auto srcptr = work.data() + p * secondary_chunkdim;
                auto copy_output = output + p * stride;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
            }

        } else {
            size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
            for (size_t p : primary_indices) {
                auto srcptr = work.data() + p;
                auto copy_output = output + p * stride;
                for (size_t x : secondary_indices) {
                    *copy_output = srcptr[x * primary_chunkdim];
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
 * The interface is "simple" as extraction of any data involves realization of the entire blob along the primary dimension
 * (i.e., the dimension used to create instances of the various `tatami::DenseExtractor` classes),
 * with no optimization for subsets of interest along that dimension.
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
     * Whether to extract a subset of elements on the primary dimension.
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
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * P`,
     * where `P` is the number of rows (if `accrow_ = true`) or columns (otherwise) in this chunk.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_length`.
     *
     * If `accrow_ = true`, we would extract data for all rows and a block of columns `[secondary_start, secondary_start + secondary_length)`;
     * conversely, if `accrow_ = false`, we would extract data for all columns and a block of rows.
     * For a primary dimension index `p` and secondary dimension index `secondary_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * The `stride` option allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomDenseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(Index_ secondary_start, Index_ secondary_length, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_, Index_>(0, core.template get_primary_chunkdim<accrow_>(), secondary_start, secondary_length, work.work, output, stride);
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * P`,
     * where `P` is the number of rows (if `accrow_ = true`) or columns (otherwise) in this chunk.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_indices.size()`.
     *
     * If `accrow_ = true`, we would extract all rows and a subset of columns in `secondary_indices`;
     * conversely, if `accrow_ = false`, we would extract data for all columns and a subset of rows.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * The `stride` option allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomDenseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& secondary_indices, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_, Index_>(0, core.template get_primary_chunkdim<accrow_>(), secondary_indices, work.work, output, stride);
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
    template<bool accrow_, typename Index_>
    void extract(Index_ secondary_start, Index_ secondary_length, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(0, core.template get_primary_chunkdim<accrow_>(), secondary_start, secondary_length, work, output, stride);
    }

    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& secondary_indices, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(0, core.template get_primary_chunkdim<accrow_>(), secondary_indices, work, output, stride);
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
 * The interface is smarter as it only extracts elements of interest along the primary dimension
 * (i.e., the dimension used to create instances of the various `tatami::DenseExtractor` classes).
 * The elements of interest may be either as a contiguous block or a indexed subset,
 * as predicted for each chunk from the `SubsettedtedOracleSlabCache`.
 * This provides some opportunities for optimization if the chunk can be partially read.
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
     * Whether to extract a subset of elements on the primary dimension.
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
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * primary_length + secondary_length`.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_length`.
     *
     * If `accrow_ = true`, we would extract a block of rows `[primary_start, primary_start + length)` and a block of columns `[secondary_start, secondary_start + secondary_length)`;
     * conversely, if `accrow_ = false`, we would extract a block of columns as the primary and the block of rows as the secondary.
     * For a primary dimension index `primary_start + p` and secondary dimension index `secondary_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(primary_start, primary_length, secondary_start, secondary_length, work.work, output, stride);
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param primary_start Index of the first element on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the first row, otherwise it is the first column.
     * @param primary_length Number of elements on the primary dimension to be extracted.
     * If `accrow_ = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (primary_start + primary_length)`.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_indices.size()`.
     *
     * If `accrow_ = true`, we would extract a block of rows `[primary_start, primary_start + length)` and a subset of columns in `secondary_indices`;
     * conversely, if `accrow_ = false`, we would extract a block of columns as the primary and a subset of rows instead.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(primary_start, primary_length, secondary_indices, work.work, output, stride);
    }

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param primary_indices Indices of the elements on the primary dimension to be extracted.
     * If `accrow_ = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (primary_indices.back() + 1)`.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_length`.
     *
     * If `accrow_ = true`, we would extract a subset of rows in `primary_indices` and a block of columns `[secondary_start, secondary_start + secondary_length)`,
     * conversely, if `accrow_ = false`, we would extract a subset of columns as the primary and the block of rows as the secondary.
     * For a primary dimension index `p` and secondary dimension index `secondary_start + i`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(primary_indices, secondary_start, secondary_length, work.work, output, stride);
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomDenseChunkedMatrix`.
     *
     * @param primary_indices Indices of the elements on the primary dimension to be extracted.
     * If `accrow_ = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output Pointer to an output array of length no less than `stride * (primary_indices.back() + 1)`.
     * @param stride Distance between corresponding values from consecutive primary dimension elements when stored in `output`.
     * This is guaranteed to be greater than or equal to `secondary_indices.size()`.
     *
     * If `accrow_ = true`, we would extract a subset of rows in `primary_indices` and a subset columns in `secondary_indices`.
     * conversely, if `accrow_ = false`, we would extract a subset of columns as the primary and the subset of rows as the secondary.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be stored in `output[p * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomDenseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(const std::vector<Index_>& primary_indices, const std::vector<Index_>& secondary_indices, Workspace& work, value_type* output, size_t stride) const {
        core.template extract<accrow_>(primary_indices, secondary_indices, work.work, output, stride);
    }
};

}

#endif
