#ifndef TATAMI_CHUNKED_MOCK_SPARSE_CHUNK_HPP
#define TATAMI_CHUNKED_MOCK_SPARSE_CHUNK_HPP

#include <vector>
#include <cstdint>

/**
 * @file mock_sparse_chunk.hpp
 * @brief Sparse chunk interface to use in a `CustomSparseChunkedMatrix`.
 */

namespace tatami_chunked {

/**
 * @cond
 */
namespace MockSparseChunk_internal {

template<typename InflatedValue_, typename InflatedIndex_>
struct Workspace {
    // Standard compressed sparse members:
    std::vector<InflatedValue_> values;
    std::vector<InflatedIndex_> indices;
    std::vector<size_t> indptrs;

    // Allocation to allow for O(1) mapping of requested indices to sparse indices.
    // This mimics what is done in the indexed sparse primary extractors in tatami proper.
    std::vector<uint8_t> remap;
};

template<class Blob_>
struct Core {
    Core() = default;

    Core(Blob_ c) : chunk(std::move(c)) {}

private:
    Blob_ chunk;

    typedef typename Blob_::value_type value_type;

    typedef typename Blob_::index_type index_type;

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

    template<typename Index_>
    static void refine_start_and_end(size_t& start, size_t& end, Index_ desired_start, Index_ desired_end, Index_ max_end, const std::vector<index_type>& indices) {
        if (desired_start) {
            auto it = indices.begin();
            start = std::lower_bound(it + start, it + end, static_cast<index_type>(desired_start)) - it;
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
                end = std::lower_bound(it + start, it + end, static_cast<index_type>(desired_end)) - it;
            }
        }
    }

private:
    // Building a present/absent mapping for the requested indices to allow O(1) look-up from the sparse matrix indices.
    template<typename Index_>
    static void configure_remap(Workspace<value_type, index_type>& work, const std::vector<Index_>& indices, size_t full) {
        work.remap.resize(full);
        for (auto i : indices) {
            work.remap[i] = 1;
        }
    }

    // Resetting just the affected indices so we can avoid a fill operation over the entire array.
    template<typename Index_>
    static void reset_remap(Workspace<value_type, index_type>& work, const std::vector<Index_>& indices) {
        for (auto i : indices) {
            work.remap[i] = 0;
        }
    }

private:
    template<bool is_block_, typename Index_>
    void fill_primary(
        Index_ p, 
        Index_ secondary_start, 
        Index_ secondary_end, 
        Index_ secondary_chunkdim,
        Workspace<value_type, index_type>& work, 
        std::vector<value_type>& current_values, 
        std::vector<index_type>& current_indices,
        index_type shift)
    const {
        auto start = work.indptrs[p], end = work.indptrs[p + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, secondary_start, secondary_end, secondary_chunkdim, work.indices);

        if constexpr(is_block_) {
            current_values.insert(current_values.end(), work.values.begin() + start, work.values.begin() + end);
            for (size_t i = start; i < end; ++i) {
                current_indices.push_back(work.indices[i] + shift);
            }
        } else {
            // Assumes that work.remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = work.indices[i];
                if (work.remap[target]) {
                    current_values.push_back(work.values[i]);
                    current_indices.push_back(target + shift);
                }
            }
        }
    }

    template<bool is_block_, typename Index_>
    void fill_secondary(
        Index_ s,
        Index_ primary_start, 
        Index_ primary_end, 
        Index_ primary_chunkdim,
        Workspace<value_type, index_type>& work, 
        std::vector<std::vector<value_type> >& output_values, 
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        auto start = work.indptrs[s], end = work.indptrs[s + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, primary_start, primary_end, primary_chunkdim, work.indices);

        if constexpr(is_block_) {
            for (size_t i = start; i < end; ++i) {
                auto p = work.indices[i];
                output_values[p].push_back(work.values[i]);
                output_indices[p].push_back(s + shift);
            }
        } else {
            // Assumes that work.remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = work.indices[i];
                if (work.remap[target]) {
                    output_values[target].push_back(work.values[i]);
                    output_indices[target].push_back(s + shift);
                }
            }
        }
    }

public:
    template<bool accrow_, typename Index_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        Index_ secondary_start,
        Index_ secondary_length,
        Workspace<value_type, index_type>& work,
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_end = primary_start + primary_length;
        Index_ secondary_end = secondary_start + secondary_length;

        if constexpr(Blob_::row_major == accrow_) {
            Index_ secondary_chunkdim = get_secondary_chunkdim<accrow_>();
            for (Index_ p = primary_start; p < primary_end; ++p) {
                fill_primary<true>(p, secondary_start, secondary_end, secondary_chunkdim, work, output_values[p], output_indices[p], shift);
            }
        } else {
            Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
            for (Index_ s = secondary_start; s < secondary_end; ++s) {
                fill_secondary<true>(s, primary_start, primary_end, primary_chunkdim, work, output_values, output_indices, shift);
            }
        }
    }

    template<bool accrow_, typename Index_>
    void extract(
        Index_ primary_start,
        Index_ primary_length,
        const std::vector<Index_>& secondary_indices,
        Workspace<value_type, index_type>& work,
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ primary_end = primary_start + primary_length;

        if constexpr(Blob_::row_major == accrow_) {
            // secondary_indices is guaranteed to be non-empty, see contracts below.
            auto secondary_start = secondary_indices.front();
            auto secondary_end = secondary_indices.back() + 1; // need 1 past end.
            auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();

            configure_remap(work, secondary_indices, secondary_chunkdim);
            for (Index_ p = primary_start; p < primary_end; ++p) {
                fill_primary<false>(p, secondary_start, secondary_end, secondary_chunkdim, work, output_values[p], output_indices[p], shift);
            }
            reset_remap(work, secondary_indices);

        } else {
            Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
            for (auto s : secondary_indices) {
                fill_secondary<true>(s, primary_start, primary_end, primary_chunkdim, work, output_values, output_indices, shift);
            }
        }
    }

public:
    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& primary_indices,
        Index_ secondary_start,
        Index_ secondary_length,
        Workspace<value_type, index_type>& work,
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ secondary_end = secondary_start + secondary_length;

        if constexpr(Blob_::row_major == accrow_) {
            Index_ secondary_chunkdim = get_secondary_chunkdim<accrow_>();
            for (auto p : primary_indices) {
                fill_primary<true>(p, secondary_start, secondary_end, secondary_chunkdim, work, output_values[p], output_indices[p], shift);
            }

        } else {
            // primary_indices is guaranteed to be non-empty, see contracts below.
            auto primary_start = primary_indices.front();
            auto primary_end = primary_indices.back() + 1; // need 1 past end.
            auto primary_chunkdim = get_primary_chunkdim<accrow_>();

            configure_remap(work, primary_indices, primary_chunkdim);
            for (Index_ s = secondary_start; s < secondary_end; ++s) {
                fill_secondary<false>(s, primary_start, primary_end, primary_chunkdim, work, output_values, output_indices, shift);
            }
            reset_remap(work, primary_indices);
        }
    }

    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& primary_indices,
        const std::vector<Index_>& secondary_indices,
        Workspace<value_type, index_type>& work,
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        chunk.inflate(work.values, work.indices, work.indptrs);

        if constexpr(Blob_::row_major == accrow_) {
            // secondary_indices is guaranteed to be non-empty, see contracts below.
            auto secondary_start = secondary_indices.front();
            auto secondary_end = secondary_indices.back() + 1; // need 1 past end.
            auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();

            configure_remap(work, secondary_indices, secondary_chunkdim);
            for (auto p : primary_indices) {
                fill_primary<false>(p, secondary_start, secondary_end, secondary_chunkdim, work, output_values[p], output_indices[p], shift);
            }
            reset_remap(work, secondary_indices);

        } else {
            // primary_indices is guaranteed to be non-empty, see contracts below.
            auto primary_start = primary_indices.front();
            auto primary_end = primary_indices.back() + 1; // need 1 past end.
            auto primary_chunkdim = get_primary_chunkdim<accrow_>();

            configure_remap(work, primary_indices, primary_chunkdim);
            for (auto s : secondary_indices) {
                fill_secondary<false>(s, primary_start, primary_end, primary_chunkdim, work, output_values, output_indices, shift);
            }
            reset_remap(work, primary_indices);
        }
    }
};

struct MockBlob {
    typedef int index_type;
    typedef double value_type;
    static constexpr bool row_major = true;

private:
    int nrows, ncols;
    std::vector<double> vcontents;
    std::vector<int> icontents;
    std::vector<size_t> pcontents;

public:
    MockBlob() = default;

    MockBlob(int nr, int nc, std::vector<double> v, std::vector<int> i, std::vector<size_t> p) : 
        nrows(nr), ncols(nc), vcontents(std::move(v)), icontents(std::move(i)), pcontents(std::move(p)) {}

    int nrow() const {
        return nrows;
    }

    int ncol() const {
        return ncols;
    }

    void inflate(std::vector<double>& vbuffer, std::vector<int>& ibuffer, std::vector<size_t>& pbuffer) const {
        vbuffer.resize(vcontents.size());
        std::copy(vcontents.begin(), vcontents.end(), vbuffer.begin());
        ibuffer.resize(icontents.size());
        std::copy(icontents.begin(), icontents.end(), ibuffer.begin());
        pbuffer.resize(pcontents.size());
        std::copy(pcontents.begin(), pcontents.end(), pbuffer.begin());
    }
};

}
/**
 * @endcond
 */

/**
 * @brief Mock a simple sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Mock a simple sparse chunk for use inside a `CustomSparseChunkedMatrix`.
 * Each chunk should represent a 2-dimensional array of numeric values.
 * The interface is "simple" as extraction of any data involves realization of the entire blob along the primary dimension
 * (i.e., the dimension used to create instances of the various `tatami::SparseExtractor` classes),
 * with no optimization for subsets of interest along that dimension.
 */
struct MockSimpleSparseChunk {
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Type of the index stored in this chunk.
     * Implementations can use any integer type. 
     */
    typedef int index_type;

    /**
     * Temporary workspace for extracting data from the chunk.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations may use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        // Hiding this here to avoid giving the impression that we NEED to implement this.
        MockSparseChunk_internal::Workspace<value_type, index_type> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the primary dimension.
     * This should be set to `false`, otherwise a `MockSubsettedSparseChunk` is expected.
     */
    static constexpr bool use_subset = false;

public:
    /**
     * @cond
     */
    // You can construct this however you like, I don't care.
    MockSimpleSparseChunk() = default;
    MockSimpleSparseChunk(int nr, int nc, std::vector<double> x, std::vector<int> i, std::vector<size_t> p) : 
        core(MockSparseChunk_internal::MockBlob(nr, nc, std::move(x), std::move(i), std::move(p))) {}
    /**
     * @endcond
     */

private:
    MockSparseChunk_internal::Core<MockSparseChunk_internal::MockBlob> core;

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param secondary_start Index of the first element on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the first column, otherwise it is the first row.
     * @param secondary_length Number of elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract data for all rows and a block of columns `[secondary_start, secondary_start + secondary_length)`;
     * conversely, if `accrow_ = false`, we would extract data for all columns and a block of rows.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        Index_ secondary_start, 
        Index_ secondary_length, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_, Index_>(
            0, 
            core.template get_primary_chunkdim<accrow_>(), 
            secondary_start, 
            secondary_length, 
            work.work, 
            output_values, 
            output_indices, 
            shift
        );
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract all rows and a subset of columns in `secondary_indices`;
     * conversely, if `accrow_ = false`, we would extract data for all columns and a subset of rows.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& secondary_indices,
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_, Index_>(
            0, 
            core.template get_primary_chunkdim<accrow_>(), 
            secondary_indices, 
            work.work, 
            output_values, 
            output_indices,
            shift
        );
    }
};

/**
 * @brief Create a sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Wraps a sparse blob in a simple chunk interface for use inside a `CustomSparseChunkedMatrix`.
 * Each blob should hold a (possibly compressed) 2-dimensional sparse array of numeric values.
 * The wrapper satisfies the `MockSimpleSparseChunk` interface, but is even simpler;
 * extraction of any data involves realization of the entire blob, with no optimization for subsets of interest along either dimension.
 *
 * The `Blob_` class should provide the following:
 *
 * - `static constexpr bool row_major`, specifying whether the array is row-major.
 * - `typedef value_type`, specifying the type of the value in the array.
 * - `typedef index_type`, specifying the type of the index in the submatrix.
 * - A `nrow() const` method, defining the number of rows in the array.
 * - A `ncol() const` method, defining the number of columns in the array.
 * - A `void inflate(std::vector<value_type>& values, std::vector<index_type>& indices, std::vector<size_t>& pointers) const` method that fills `values`, `indices` and `pointers`,
 *   each with the corresponding field for compressed sparse matrix format.
 *   This should constitute a CSR matrix if `row_major = true` and a CSC matrix otherwise.
 *
 * @tparam Blob_ Class to represent a simple chunk.
 */
template<class Blob_>
struct SimpleSparseChunkWrapper {
    /**
     * @cond
     */
    typedef typename Blob_::value_type value_type;

    typedef typename Blob_::index_type index_type;

    typedef MockSparseChunk_internal::Workspace<value_type, index_type> Workspace;

    static constexpr bool use_subset = false;

    SimpleSparseChunkWrapper() = default;

    SimpleSparseChunkWrapper(Blob_ c) : core(std::move(c)) {}
    /**
     * @endcond
     */

private:
    MockSparseChunk_internal::Core<Blob_> core;

public:
    /**
     * @cond
     */
    template<bool accrow_, typename Index_>
    void extract(
        Index_ secondary_start, 
        Index_ secondary_length, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_, Index_>(
            0, 
            core.template get_primary_chunkdim<accrow_>(), 
            secondary_start, 
            secondary_length, 
            work, 
            output_values, 
            output_indices, 
            shift
        );
    }

    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& secondary_indices,
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_, Index_>(
            0, 
            core.template get_primary_chunkdim<accrow_>(), 
            secondary_indices, 
            work, 
            output_values, 
            output_indices,
            shift
        );
    }
    /**
     * @endcond
     */
};

/**
 * @brief Mock a subsettable sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Mock a subsettable sparse chunk for use inside a `CustomSparseChunkedMatrix`.
 * Each chunk should represent a (possible compressed) 2-dimensional array of numeric values.
 * The interface is smarter as it only extracts elements of interest along the primary dimension
 * (i.e., the dimension used to create instances of the various `tatami::SparseExtractor` classes).
 * The elements of interest may be either as a contiguous block or a indexed subset,
 * as predicted for each chunk from the `SubsettedtedOracleSlabCache`.
 * This provides some opportunities for optimization if the chunk can be partially read.
 */
struct MockSubsettedSparseChunk {
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Type of the index stored in this chunk.
     * Implementations can use any integer type. 
     */
    typedef int index_type;

    /**
     * Workspace for chunk extraction.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations maye use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        // Hiding this here to avoid giving the impression that we NEED to implement this.
        MockSparseChunk_internal::Workspace<value_type, index_type> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the primary dimension.
     * This should be set to `true`, otherwise a `MockSimpleSparseChunk` is expected.
     */
    static constexpr bool use_subset = true;

public:
    /**
     * @cond
     */
    MockSubsettedSparseChunk() = default;
    MockSubsettedSparseChunk(int nr, int nc, std::vector<double> x, std::vector<int> i, std::vector<size_t> p) : 
        core(MockSparseChunk_internal::MockBlob(nr, nc, std::move(x), std::move(i), std::move(p))) {}
    /**
     * @endcond
     */

private:
    MockSparseChunk_internal::Core<MockSparseChunk_internal::MockBlob> core;

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
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
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract a block of rows `[primary_start, primary_start + length)` and a block of columns `[secondary_start, secondary_start + secondary_length)`;
     * conversely, if `accrow_ = false`, we would extract a block of columns as the primary and the block of rows as the secondary.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        Index_ primary_start, 
        Index_ primary_length, 
        Index_ secondary_start, 
        Index_ secondary_length, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_>(
            primary_start, 
            primary_length, 
            secondary_start, 
            secondary_length, 
            work.work,
            output_values,
            output_indices,
            shift
        );
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
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
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_start + primary_length`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract a block of rows `[primary_start, primary_start + length)` and a subset of columns in `secondary_indices`;
     * conversely, if `accrow_ = false`, we would extract a block of columns as the primary and a subset of rows instead.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        Index_ primary_start, 
        Index_ primary_length, 
        const std::vector<Index_>& secondary_indices, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_>(
            primary_start, 
            primary_length, 
            secondary_indices, 
            work.work, 
            output_values, 
            output_indices, 
            shift
        );
    }

public:
    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
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
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_indices.back() + 1`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_indices.back() + 1`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract a subset of rows in `primary_indices` and a block of columns `[secondary_start, secondary_start + secondary_length)`,
     * conversely, if `accrow_ = false`, we would extract a subset of columns as the primary and the block of rows as the secondary.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& primary_indices, 
        Index_ secondary_start, 
        Index_ secondary_length, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_>(
            primary_indices, 
            secondary_start, 
            secondary_length, 
            work.work, 
            output_values,
            output_indices,
            shift
        );
    }

    /**
     * @tparam accrow_ Whether the rows are the primary dimension.
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param primary_indices Indices of the elements on the primary dimension to be extracted.
     * If `accrow_ = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param secondary_indices Indices of the elements on the secondary dimension to be extracted.
     * If `accrow_ = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of vectors in which to store the output values.
     * The outer vector is of length no less than `primary_indices.back() + 1`; each inner vector corresponds to an element of the primary dimension.
     * @param[out] output_indices Vector of vectors in which to store the output indices.
     * The outer vector is of length no less than `primary_indices.back() + 1`; each inner vector corresponds to an element of the primary dimension.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `accrow_ = true`, we would extract a subset of rows in `primary_indices` and a subset columns in `secondary_indices`.
     * conversely, if `accrow_ = false`, we would extract a subset of columns as the primary and the subset of rows as the secondary.
     * For a non-zero entry in primary dimension index `p`, the value from the chunk should be appended to `output_values[p]`.
     * The secondary index for this non-zero entry should be increased by `shift` and then appended to `output_indices[p]`.
     * The method should maintain a strictly increasing order among the appended secondary indices.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here.
     * Only the `accrow_` template parameter is explicitly passed when this method is called by `CustomSparseChunkedMatrix`.
     */
    template<bool accrow_, typename Index_>
    void extract(
        const std::vector<Index_>& primary_indices, 
        const std::vector<Index_>& secondary_indices, 
        Workspace& work, 
        std::vector<std::vector<value_type> >& output_values,
        std::vector<std::vector<index_type> >& output_indices,
        index_type shift)
    const {
        core.template extract<accrow_>(
            primary_indices, 
            secondary_indices, 
            work.work, 
            output_values,
            output_indices,
            shift 
        );
    }
};

}

#endif
