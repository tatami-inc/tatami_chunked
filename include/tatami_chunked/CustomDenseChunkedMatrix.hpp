#ifndef TATAMI_CHUNKED_CUSTOM_DENSE_CHUNKED_MATRIX_HPP
#define TATAMI_CHUNKED_CUSTOM_DENSE_CHUNKED_MATRIX_HPP

#include "tatami/tatami.hpp"
#include "custom_chunk_coordinator.hpp"

#include <vector>

/**
 * @file CustomDenseChunkedMatrix.hpp
 * @brief Custom dense chunked matrix.
 */

namespace tatami_chunked {

/**
 * @brief Options for custom dense chunk extraction.
 */
struct CustomDenseChunkedOptions {
    /**
     * Size of the in-memory cache in bytes.
     * Larger caches improve access speed at the cost of memory usage.
     * Small values may be ignored if `require_minimum_cache` is `true`.
     */
    size_t maximum_cache_size = 100000000;

    /**
     * Whether to automatically enforce a minimum size for the cache, regardless of `maximum_cache_size`.
     * This minimum is chosen to ensure that a single slab can be retained in memory,
     * so that the same chunks are not repeatedly re-read when iterating over consecutive rows/columns of the matrix.
     */
    bool require_minimum_cache = true;
};

/**
 * @cond
 */
namespace CustomChunkedMatrix_internal {

const Value_* process_dense_slab(const std::pair<const Slab*, Index_>& fetched, Value_* buffer, size_t secondary_length) {
    auto ptr = fetched.first->data() + static_cast<size_t>(fetched.second) * secondary_length; // cast to size_t to avoid overflow.
    std::copy_n(ptr, secondary_length, buffer);
    return buffer;
}

template<bool accrow_, typename Value_, typename Index_, typename Chunk_>
struct DenseBaseSolo {
private:
    const ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>::DenseSlab Slab;

    // These two instances are not fully allocated Slabs; rather, tmp_solo just
    // holds the content for a single chunk, while final_solo holds the content
    // across chunks but only for the requested dimension element. Both cases
    // are likely to be much smaller than a full Slab, so we're already more
    // memory-efficient than trying to use num_slabs_in_cache = 1.
    Slab tmp_solo;
    Slab final_solo;

public:
    DenseBaseSolo(const ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>& coordinator, Index_ secondary_length) : 
        coordinator(coordinator),
        tmp_solo(static_cast<size_t>(coordinator.get_chunk_nrow()) * static_cast<size_t>(coordinator.get_chunk_ncol())),
        final_solo(secondary_length)
    {}

    ~DenseBaseSolo() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(Index_ i, Args_&& ... args) {
        return coordinator.template fetch_solo<accrow_>(i, std::forward<Args_>(args)..., chunk_workspace, tmp_solo, final_solo);
    }
};

template<bool accrow_, typename Value_, typename Index_, typename Chunk_>
struct DenseBaseMyopic {
protected:
    const ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>::DenseSlab Slab;
    LruSlabCache<Index_, Slab_> cache;

public:
    DenseBaseMyopic(const ChunkCoordinator<Index_, Chunk_::use_subset, Chunk_, int>& coordinator, size_t max_slabs_in_cache) :
        coordinator(coordinator), cache(max_slabs_in_cache) {}

    ~DenseBaseMyopic() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(Index_ i, Args_&& ... args) {
        return this->coordinator.template fetch_myopic<accrow_>(i, std::forward<Args_>(args)..., chunk_workspace, cache);
    }
};

template<bool accrow_, typename Value_, typename Index_, typename Chunk_>
struct DenseBaseOracular {
protected:
    const ChunkCoordinator<Index_, false, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, false, Chunk_, int>::DenseSlab Slab;
    typename std::conditional<subset_, OracularSubsettedSlabCache<Index_, Index_, Slab>, OracularSlabCache<Index_, Index_, Slab> >::type cache;

public:
    DenseBaseOracular(const ChunkCoordinator<Index_, false, Chunk_, int>& coordinator, size_t max_slabs_in_cache, tatami::MaybeOracle<oracle_, Index_> oracle) : 
        coordinator(coordinator), cache(std::move(oracle), max_slabs_in_cache) {}

    ~DenseBaseOracular() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(Index_ i, Args_&& ... args) {
        return this->coordinator.template fetch_oracle<accrow_>(i, std::forward<Args_>(args)..., chunk_workspace, cache);
    }
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
using DenseBase = typename std::conditional<solo_, 
      DenseBaseSolo<accrow_, Value_, Index_, Chunk_>,
      typename std::conditional<oracle_,
          DenseBaseOracular<accrow_, Value_, Index_, Chunk_>,
          DenseBaseMyopic<accrow_, Value_, Index_, Chunk_>
      >::type
>::type;

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DenseFull : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    DenseFull(
        const ChunkCoordinator<Index_, false, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle) :
        DenseBase<accrow_, oracle_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            coordinator.template get_secondary_dim<accrow_>()
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto secondary_len = this->coordinator.template get_secondary_dim<accrow_>();
        auto fetched = this->coordinator.template fetch<accrow_, oracle_>(
            i, 
            0, 
            secondary_len,
            this->cache_workspace,
            this->chunk_workspace, 
            this->tmp_solo,
            this->final_solo
        );
        return this->process_dense_slab(fetched, buffer, secondary_len);
    }
};

template<bool accrow_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DenseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, oracle_, Value_, Index_, Chunk_> {
    DenseBlock(
        const ChunkCoordinator<Index_, false, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length) : 
        DenseBase<accrow_, oracle_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            block_length
        ),
        block_start(block_start),
        block_length(block_length)
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = this->coordinator.template fetch<accrow_, oracle_>(
            i, 
            block_start, 
            block_length,
            this->cache_workspace, 
            this->chunk_workspace, 
            this->tmp_solo,
            this->final_solo
        );
        return this->process_dense_slab(fetched, buffer, block_length);
    }

private:
    Index_ block_start, block_length;
};

template<bool accrow_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DenseIndex : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, oracle_, Value_, Index_, Chunk_> {
    DenseIndex(
        const ChunkCoordinator<Index_, false, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> idx_ptr) :
        DenseBase<accrow_, oracle_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            idx_ptr->size()
        ),
        indices_ptr(std::move(idx_ptr))
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = this->coordinator.template fetch<accrow_, oracle_>(
            i, 
            *indices_ptr,
            tmp_indices,
            this->cache_workspace, 
            this->chunk_workspace, 
            this->tmp_solo,
            this->final_solo
        );
        return this->process_dense_slab(fetched, buffer, indices_ptr->size());
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
    std::vector<Index_> tmp_indices;
};

}
/**
 * @endcond
 */

/**
 * @brief Matrix of custom dense chunks.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Chunk_ Class of the chunk, implementing either the `MockSimpleDenseChunk` or `MockSubsetDenseChunk` interfaces.
 *
 * Implements a `Matrix` subclass where data is contained in dense rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage; on access, each chunk is decompressed and the desired values are extracted.
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero - the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<typename Value_, typename Index_, typename Chunk_>
class CustomDenseChunkedMatrix : public tatami::Matrix<Value_, Index_> {
public:
    /**
     * @param mat_nrow Number of rows in the matrix.
     * @param mat_ncol Number of columns in the matrix.
     * @param chunk_nrow Number of rows in each chunk.
     * @param chunk_ncol Number of columns in each chunk.
     * @param chunks Vector containing a 2D array of chunks that cover the entire matrix.
     * This should have length equal to the product of the number of chunks along the rows and columns of the matrix, i.e., `ceil(mat_nrow / chunk_nrow) * ceil(mat_ncol / chunk_ncol)`.
     * @param row_major Whether `chunks` is in row-major format.
     * @param opt Further options for chunked extraction.
     */
    CustomDenseChunkedMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomDenseChunkedOptions& opt) : 
        coordinator(ChunkDimensionStats<Index_>(mat_nrow, chunk_nrow), ChunkDimensionStats<Index_>(mat_ncol, chunk_ncol), std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / sizeof(typename Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {}

private:
    CustomChunkedMatrix_internal::ChunkCoordinator<Index_, false, Chunk_, int> coordinator;
    size_t cache_size_in_elements;
    bool require_minimum_cache;

public:
    Index_ nrow() const { 
        return coordinator.get_nrow(); 
    }

    Index_ ncol() const { 
        return coordinator.get_ncol(); 
    }

    bool prefer_rows() const { 
        return coordinator.prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(coordinator.prefer_rows_internal());
    }

    bool sparse() const {
        return false;
    }

    double sparse_proportion() const {
        return 0;
    }

    using tatami::Matrix<Value_, Index_>::dense;

    using tatami::Matrix<Value_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options&) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::DenseFull<true, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseFull<false, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle)
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options&) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::DenseBlock<true, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                block_start, 
                block_length
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseBlock<false, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle),
                block_start, 
                block_length
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options&) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::DenseIndex<true, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                std::move(indices_ptr)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseIndex<false, oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                std::move(indices_ptr)
            );
        }
    }

public:
    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
private:
    Index_ get_secondary_dim(bool row) const {
        return row ? coordinator.get_ncol() : coordinator.get_nrow();
    }

public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const {
        return std::make_unique<tatami::FullSparsifiedWrapper<false, Value_, Index_> >(dense(row, opt), get_secondary_dim(row), opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::make_unique<tatami::BlockSparsifiedWrapper<false, Value_, Index_> >(dense(row, block_start, block_length, opt), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        auto d = dense(row, indices_ptr, opt);
        return std::make_unique<tatami::IndexSparsifiedWrapper<false, Value_, Index_> >(std::move(d), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        const tatami::Options& opt)
    const {
        return std::make_unique<tatami::FullSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(oracle), opt), get_secondary_dim(row), opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return std::make_unique<tatami::BlockSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(oracle), block_start, block_length, opt), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt)
    const {
        auto d = dense(row, std::move(oracle), indices_ptr, opt);
        return std::make_unique<tatami::IndexSparsifiedWrapper<true, Value_, Index_> >(std::move(d), std::move(indices_ptr), opt);
    }
};

}

#endif
