#ifndef TATAMI_CHUNKED_SIMPLE_CHUNK_WRAPPERS_HPP
#define TATAMI_CHUNKED_SIMPLE_CHUNK_WRAPPERS_HPP

#include <vector>

/**
 * @file simple_chunk_wrappers.hpp
 * @brief Wrappers for simple chunks to use in a chunked matrix.
 */

namespace tatami_chunked {

/**
 * @brief Wrap a simple dense chunk for a `CustomChunkedDenseMatrix`.
 *
 * Implements a simple wrapper around a dense chunk for use inside a `CustomChunkedDenseMatrix`.
 * Each chunk should hold a (typically compressed) 2-dimensional array of numeric values.
 * This chunk is considered to be "simple", as extraction of any data involves realization of the entire chunk.
 *
 * The `SimpleChunk_` class should provide the following:
 *
 * - A `static constexpr bool row_major`, specifying whether the array is row-major.
 * - A `typedef value_type`, specifying the type of the value in the array.
 * - A `nrow() const` method, defining the number of rows in the array.
 * - A `ncol() const` method, defining the number of columns in the array.
 * - A `void inflate(std::vector<value_type>& buffer) const` method that fills `buffer` with the contents of the array.
 *   This should be filled in row-major format if `row_major = true` and in column-major format otherwise.
 *
 * @tparam SimpleChunk_ Class to represent a simple chunk.
 */
template<class SimpleChunk_>
struct SimpleDenseChunkWrapper {
    /**
     * Type of the value stored in this chunk.
     */
    typedef typename SimpleChunk_::value_type value_type;

    /**
     * Workspace for chunk extraction.
     * This can be used in multiple `fetch()` calls, possibly across different chunks.
     */
    typedef std::vector<value_type> Workspace;

public:
    /**
     * Default constructor.
     */
    SimpleDenseChunkWrapper() = default;

    /**
     * @param c Chunk to be wrapped.
     */
    SimpleDenseChunkWrapper(SimpleChunk_ c) : chunk(std::move(c)) {}

private:
    SimpleChunk_ chunk;

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
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam Output_ Numeric type for the output.
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
     * @param[out] output Pointer to an output array of length no less than `primary_length * stride`.
     * @param stride Stride separating corresponding values from consecutive elements on the primary dimension.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns `[secondary_start, secondary_start + secondary_length)`.
     * For a primary dimension index `p` and secondary dimension index `s`, the value from the chunk should be stored in `output[(p - primary_start) * stride + (s - secondary_start)]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename Output_>
    void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output, size_t stride) const {
        chunk.inflate(work);
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>();

        if constexpr(SimpleChunk_::row_major == accrow_) {
            auto srcptr = work.data() + primary_start * secondary_chunkdim + secondary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                std::copy(srcptr, srcptr + secondary_length, output);
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            auto srcptr = work.data() + secondary_start * primary_chunkdim + primary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_srcptr = srcptr;
                auto copy_output = output;
                for (size_t s = 0; s < secondary_length; ++s) {
                    *copy_output = *copy_srcptr;
                    ++copy_output;
                    copy_srcptr += primary_chunkdim;
                }
                ++srcptr;
                output += stride;
            }
        }
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam Output_ Numeric type for the output.
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
     * @param[out] output Pointer to an output array of length no less than `primary_length * stride`.
     * @param stride Stride separating corresponding values from consecutive elements on the primary dimension.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns in `secondary_indices`.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be stored in `output[(p - primary_start) * stride + i]`.
     * This layout allows concatenation of multiple chunks into a single contiguous array for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename Output_>
    void extract(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const {
        chunk.inflate(work);
        size_t primary_chunkdim = get_primary_chunkdim<accrow_>(); // use size_t to avoid integer overflow with Index_.
        size_t secondary_chunkdim = get_secondary_chunkdim<accrow_>(); 

        if constexpr(SimpleChunk_::row_major == accrow_) {
            auto srcptr = work.data() + primary_start * secondary_chunkdim;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x];
                    ++copy_output;
                }
                srcptr += secondary_chunkdim;
                output += stride;
            }

        } else {
            auto srcptr = work.data() + primary_start;
            for (size_t p = 0; p < primary_length; ++p) {
                auto copy_output = output;
                for (auto x : secondary_indices) {
                    *copy_output = srcptr[x * primary_chunkdim];
                    ++copy_output;
                }
                ++srcptr;
                output += stride;
            }
        }
    }
};

/**
 * @brief Wrap a simple sparse chunk for a `CustomChunkedSparseMatrix`.
 *
 * Implements a simple wrapper around a simple sparse chunk for use inside a `CustomChunkedSparseMatrix`.
 * Each chunk should hold a (typically compressed) 2-dimensional sparse submatrix of numeric values.
 * This chunk is considered to be "simple", as extraction of any data involves realization of the entire chunk.
 *
 * The `SimpleChunk_` class should provide the following:
 *
 * - A `static constexpr bool row_major`, specifying whether the submatrix is in the compressed sparse row layout.
 * - A `typedef value_type`, specifying the type of the value in the submatrix.
 * - A `typedef index_type`, specifying the type of the index in the submatrix.
 * - A `nrow() const` method, defining the number of rows in the submatrix.
 * - A `ncol() const` method, defining the number of columns in the submatrix.
 * - A `void inflate(std::vector<value_type>& values, std::vector<index_type>& indices, std::vector<size_t>& pointers) const` method that fills `values`, `indices` and `pointers` with the standard compressed sparse data.
 *   This should be a CSR matrix if `row_major = true` and a CSC matrix otherwise.
 *
 * @tparam SimpleChunk_ Class to represent the chunk.
 */
template<class SimpleChunk_>
struct SimpleSparseChunkWrapper {
    /**
     * Type of the value stored in this chunk.
     */
    typedef typename SimpleChunk_::value_type value_type;

    /**
     * Type of the index stored in this chunk.
     */
    typedef typename SimpleChunk_::index_type index_type;

    /**
     * @brief Workspace for chunk extraction.
     *
     * This can be used in multiple `fetch()` calls, possibly across different chunks.
     */
    struct Workspace {
        /**
         * @cond
         */
        std::vector<value_type> values;
        std::vector<index_type> indices;
        std::vector<size_t> indptrs;
        /**
         * @endcond
         */
    };

public:
    /**
     * Default constructor.
     */
    SimpleSparseChunkWrapper() = default;

    /**
     * @param c Chunk to be wrapped.
     */
    SimpleSparseChunkWrapper(SimpleChunk_ c) : chunk(std::move(c)) {}

private:
    SimpleChunk_ chunk;

    template<bool accrow_>
    size_t get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.nrow();
        } else {
            return chunk.ncol();
        }
    }

    template<bool accrow_>
    size_t get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk.ncol();
        } else {
            return chunk.nrow();
        }
    }

    template<typename Index_>
    static void refine_start_and_end(size_t& start, size_t& end, Index_ desired_start, Index_ desired_end, Index_ max_end, const std::vector<typename SimpleChunk_::index_type>& indices) {
        if (desired_start) {
            auto it = indices.begin();
            start = std::lower_bound(it + start, it + end, static_cast<typename SimpleChunk_::index_type>(desired_start)) - it;
        }

        if (desired_end != max_end) {
            if (desired_end == desired_start + 1) {
                if (start != end && indices[start] == desired_start) {
                    end = start + 1;
                } else {
                    end = start;
                }
            } else {
                auto it = indices.begin();
                end = std::lower_bound(it + start, it + end, static_cast<typename SimpleChunk_::index_type>(desired_end)) - it;
            }
        }
    }

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam OutputValue_ Numeric type for the output values.
     * @tparam OutputIndex_ Integer type for the output indices.
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
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns `[secondary_start, secondary_start + secondary_length)`.
     * For a primary dimension index `p` and secondary dimension index `s`, the value from the chunk should be appended to `output_values[p - primary_start]`.
     * Similarly, the secondary index should be increased by `shift` and then appended to `output_indices[p - primary_start]`.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename OutputValue_, typename OutputIndex_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        Index_ secondary_start,
        Index_ secondary_length,
        Workspace& work,
        std::vector<std::vector<OutputValue_> >& output_values,
        std::vector<std::vector<OutputIndex_> >& output_indices,
        OutputIndex_ shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ primary_end = primary_start + primary_length;
        Index_ secondary_end = secondary_start + secondary_length;

        if constexpr(SimpleChunk_::row_major == accrow_) {
            for (Index_ p = primary_start; p < primary_end; ++p) {
                auto start = work.indptrs[p], end = work.indptrs[p + 1];

                if (start < end) {
                    refine_start_and_end(start, end, secondary_start, secondary_end, secondary_chunkdim, work.indices);

                    auto& current_values = output_values[p - primary_start];
                    current_values.insert(current_values.end(), work.values.begin() + start, work.values.begin() + end);

                    auto& current_indices = output_indices[p - primary_start];
                    for (size_t i = start; i < end; ++i) {
                        current_indices.push_back(work.indices[i] + shift);
                    }
                }
            }

        } else {
            for (size_t s = secondary_start; s < secondary_end; ++s) {
                auto start = work.indptrs[s], end = work.indptrs[s + 1];
                refine_start_and_end(start, end, primary_start, primary_end, primary_chunkdim, work.indices);

                for (size_t i = start; i < end; ++i) {
                    auto p = work.indices[i] - primary_start;
                    output_values[p].push_back(work.values[i]);
                    output_indices[p].push_back(s + shift);
                }
            }
        }
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the chunk.
     * @tparam OutputValue_ Numeric type for the output values.
     * @tparam OutputIndex_ Integer type for the output indices.
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
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_length`; each inner vector corresponds to an element of the primary dimension, starting at `primary_start`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * This method extracts the specified values from the chunk into `output`.
     * For example, if `accrow_ = true`, we would extract rows `[primary_start, primary_start + length)` and columns in `secondary_indices`.
     * For a primary dimension index `p` and secondary dimension index `secondary_indices[i]`, the value from the chunk should be appended to `output_values[p - primary_start]`.
     * Similarly, the secondary index should be increased by `shift` and then appended to `output_indices[p - primary_start]`.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomChunkedMatrix`.
     */
    template<bool accrow_, typename Index_, typename OutputValue_, typename OutputIndex_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        const std::vector<Index_>& secondary_indices,
        Workspace& work,
        std::vector<std::vector<OutputValue_> >& output_values,
        std::vector<std::vector<OutputIndex_> >& output_indices,
        OutputIndex_ shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ primary_end = primary_start + primary_length;

        if constexpr(SimpleChunk_::row_major == accrow_) {
            for (Index_ p = primary_start; p < primary_end; ++p) {
                auto start = work.indptrs[p], end = work.indptrs[p + 1];

                if (start < end) {
                    if (secondary_indices.front()) {
                        auto it = work.indices.begin();
                        start = std::lower_bound(it + start, it + end, static_cast<typename SimpleChunk_::index_type>(secondary_indices.front())) - it;
                    }

                    auto sIt = secondary_indices.begin();
                    auto& current_values = output_values[p - primary_start];
                    auto& current_indices = output_indices[p - primary_start];

                    for (size_t i = start; i < end; ++i) {
                        Index_ target = work.indices[i];
                        while (sIt != secondary_indices.end() && *sIt < target) {
                            ++sIt;
                        }
                        if (sIt == secondary_indices.end()) {
                            break;
                        }
                        if (*sIt == target) {
                            current_values.push_back(work.values[i]);
                            current_indices.push_back(work.indices[i] + shift);
                            ++sIt;
                        }
                    }
                }
            }

        } else {
            for (auto s : secondary_indices) {
                auto start = work.indptrs[s], end = work.indptrs[s + 1];
                refine_start_and_end(start, end, primary_start, primary_end, primary_chunkdim, work.indices);

                for (size_t i = start; i < end; ++i) {
                    auto p = work.indices[i] - primary_start;
                    output_values[p].push_back(work.values[i]);
                    output_indices[p].push_back(s + shift);
                }
            }
        }
    }
};

}

#endif
