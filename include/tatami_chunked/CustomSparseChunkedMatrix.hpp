#ifndef TATAMI_CHUNKED_CUSTOM_SPARSE_CHUNKED_MATRIX_HPP
#define TATAMI_CHUNKED_CUSTOM_SPARSE_CHUNKED_MATRIX_HPP

#include "tatami/tatami.hpp"
#include "custom_internals.hpp"
#include "LruSlabCache.hpp"
#include "OracularSlabCache.hpp"
#include "OracularSubsettedSlabCache.hpp"

#include <vector>

/**
 * @file CustomSparseChunkedMatrix.hpp
 * @brief Custom sparse chunked matrix.
 */

namespace tatami_chunked {

/**
 * @brief Options for custom sparse chunk extraction.
 */
struct CustomSparseChunkedOptions {
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

/*********************
 **** Base classes ***
 *********************/

template<bool accrow_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct SparseBaseSolo {
protected:
    const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, true, Chunk_, int>::SparseSlab Slab;

    tatami::MaybeOracle<oracle_, Index_> oracle;
    typename std::conditional<oracle_, size_t, bool>::type counter = 0;

    // These two instances are not fully allocated Slabs; rather, tmp_solo just
    // holds the content for a single chunk, while final_solo holds the content
    // across chunks but only for the requested dimension element. Both cases
    // are likely to be much smaller than a full Slab, so we're already more
    // memory-efficient than trying to use num_slabs_in_cache = 1.
    Slab tmp_solo;
    Slab final_solo;

public:
    SparseBaseSolo(
        const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator, 
        [[maybe_unused]] size_t max_slabs_in_cache, // for consistency with the other base classes.
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ secondary_length) :
        coordinator(coordinator),
        oracle(std::move(ora)),
        tmp_solo(static_cast<size_t>(coordinator.get_chunk_nrow()) * static_cast<size_t>(coordinator.get_chunk_ncol())),
        final_solo(secondary_length)
    {}

    ~SparseBaseSolo() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(Index_ i, Args_&& ... args) {
        if constexpr(oracle_) {
            i = oracle->get(counter++);
        }
        return coordinator.template fetch_single<accrow_>(i, std::forward<Args_>(args)..., chunk_workspace, tmp_solo, final_solo);
    }
};

template<bool accrow_, typename Value_, typename Index_, typename Chunk_>
struct SparseBaseMyopic {
protected:
    const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, true, Chunk_, int>::SparseSlab Slab;
    LruSlabCache<Index_, Slab> cache;

public:
    SparseBaseMyopic(
        const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator,
        size_t max_slabs_in_cache, 
        [[maybe_unused]] tatami::MaybeOracle<false, Index_> ora, // for consistency with the other base classes
        [[maybe_unused]] Index_ secondary_length) : // ditto
        coordinator(coordinator),
        cache(max_slabs_in_cache) 
    {}

    ~SparseBaseMyopic() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(Index_ i, Args_&& ... args) {
        return this->coordinator.template fetch_myopic<accrow_>(i, std::forward<Args_>(args)..., chunk_workspace, cache);
    }
};

template<bool accrow_, typename Value_, typename Index_, typename Chunk_>
struct SparseBaseOracular {
protected:
    const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, true, Chunk_, int>::SparseSlab Slab;
    typename std::conditional<Chunk_::use_subset, OracularSubsettedSlabCache<Index_, Index_, Slab>, OracularSlabCache<Index_, Index_, Slab> >::type cache;

public:
    SparseBaseOracular(
        const ChunkCoordinator<Index_, true, Chunk_, int>& coordinator,
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<true, Index_> ora,
        [[maybe_unused]] Index_ secondary_length) : // for consistency with the other base classes
        coordinator(coordinator), 
        cache(std::move(ora), max_slabs_in_cache) 
    {}

    ~SparseBaseOracular() = default;

protected:
    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw([[maybe_unused]] Index_ i, Args_&& ... args) {
        return this->coordinator.template fetch_oracular<accrow_>(std::forward<Args_>(args)..., chunk_workspace, cache);
    }
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
using SparseBase = typename std::conditional<solo_, 
      SparseBaseSolo<accrow_, oracle_, Value_, Index_, Chunk_>,
      typename std::conditional<oracle_,
          SparseBaseOracular<accrow_, Value_, Index_, Chunk_>,
          SparseBaseMyopic<accrow_, Value_, Index_, Chunk_>
      >::type
>::type;

/***********************
 **** Sparse classes ***
 ***********************/

template<class Slab_, typename Index_, typename Value_>
tatami::SparseRange<Value_, Index_> process_sparse_slab(const std::pair<const Slab_*, Index_>& fetched, Value_* vbuffer, Index_* ibuffer, bool needs_value, bool needs_index) {
    const auto& values = fetched.first->values[fetched.second];
    const auto& indices = fetched.first->indices[fetched.second];

    if (needs_value) {
        std::copy(values.begin(), values.end(), vbuffer);
    } else {
        vbuffer = NULL;
    }

    if (needs_index) {
        std::copy(indices.begin(), indices.end(), ibuffer);
    } else {
        ibuffer = NULL;
    }

    return tatami::SparseRange<Value_, Index_>(values.size(), vbuffer, ibuffer);
}

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct SparseFull : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    SparseFull(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        const tatami::Options& opt) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), coordinator.template get_secondary_dim<accrow_>()),
        needs_value(opt.sparse_extract_value),
        needs_index(opt.sparse_extract_index)
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        Index_ secondary_len = this->coordinator.template get_secondary_dim<accrow_>();
        auto fetched = this->fetch_raw(i, 0, secondary_len);
        return process_sparse_slab(fetched, vbuffer, ibuffer, needs_value, needs_index); 
    }

private:
    bool needs_value, needs_index;
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct SparseBlock : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    SparseBlock(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), block_length),
        block_start(block_start),
        block_length(block_length),
        needs_value(opt.sparse_extract_value),
        needs_index(opt.sparse_extract_index)
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto fetched = this->fetch_raw(i, block_start, block_length);
        return process_sparse_slab(fetched, vbuffer, ibuffer, needs_value, needs_index);
    }

private:
    Index_ block_start, block_length;
    bool needs_value, needs_index;
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct SparseIndex : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    SparseIndex(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        tatami::VectorPtr<Index_> idx_ptr, 
        const tatami::Options& opt) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), idx_ptr->size()),
        indices_ptr(std::move(idx_ptr)),
        needs_value(opt.sparse_extract_value),
        needs_index(opt.sparse_extract_index)
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto fetched = this->fetch_raw(i, *indices_ptr, tmp_indices);
        return process_sparse_slab(fetched, vbuffer, ibuffer, needs_value, needs_index);
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
    std::vector<Index_> tmp_indices;
    bool needs_value, needs_index;
};

/**************************
 **** Densified classes ***
 **************************/

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedFull : public tatami::DenseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    DensifiedFull(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        const tatami::Options&) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), coordinator.template get_secondary_dim<accrow_>())
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto secondary_length = this->coordinator.template get_secondary_dim<accrow_>();
        auto contents = this->fetch_raw(i, 0, secondary_length);
        const auto& values = contents.first->values[contents.second];
        const auto& indices = contents.first->indices[contents.second];
        std::fill_n(buffer, secondary_length, 0);
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            buffer[indices[i]] = values[i];
        }
        return buffer;
    }
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedBlock : public tatami::DenseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    DensifiedBlock(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora,
        Index_ block_start, 
        Index_ block_length,
        const tatami::Options&) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), block_length),
        block_start(block_start),
        block_length(block_length)
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto contents = this->fetch_raw(i, block_start, block_length);
        const auto& values = contents.first->values[contents.second];
        const auto& indices = contents.first->indices[contents.second];
        std::fill_n(buffer, block_length, 0);
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            buffer[indices[i] - block_start] = values[i];
        }
        return buffer;
    }

private:
    Index_ block_start, block_length;
};

template<bool accrow_, bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedIndex : public tatami::DenseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_> {
    DensifiedIndex(
        const ChunkCoordinator<Index_, true, Chunk_>& coordinator, 
        size_t max_slabs_in_cache,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        tatami::VectorPtr<Index_> idx_ptr,
        const tatami::Options&) :
        SparseBase<accrow_, solo_, oracle_, Value_, Index_, Chunk_>(coordinator, max_slabs_in_cache, std::move(ora), idx_ptr->size()),
        indices_ptr(std::move(idx_ptr))
    {
        const auto& indices = *indices_ptr;
        if (!indices.empty()) {
            remap_offset = indices.front();
            size_t alloc = indices.back() - remap_offset + 1;
            remap.resize(alloc);
            Index_ counter = 0;
            for (auto i : indices) {
                remap[i - remap_offset] = counter;
                ++counter;
            }
        }
    }

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto nidx = indices_ptr->size();
        auto contents = this->fetch_raw(i, *indices_ptr, tmp_indices);
        const auto& values = contents.first->values[contents.second];
        const auto& indices = contents.first->indices[contents.second];
        std::fill_n(buffer, nidx, 0);
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            buffer[remap[indices[i] - remap_offset]] = values[i];
        }
        return buffer;
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
    Index_ remap_offset = 0;
    std::vector<Index_> remap;
    std::vector<Index_> tmp_indices;
};

}
/**
 * @endcond
 */

/**
 * @brief Matrix of custom sparse chunks.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Chunk_ Class of the chunk, implementing either the `MockSimpleSparseChunk` or `MockSubsetSparseChunk` interfaces.
 *
 * Implements a `Matrix` subclass where data is contained in sparse rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage; on access, each chunk is decompressed and the desired values are extracted.
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero - the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<typename Value_, typename Index_, typename Chunk_>
class CustomSparseChunkedMatrix : public tatami::Matrix<Value_, Index_> {
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
    CustomSparseChunkedMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomSparseChunkedOptions& opt) : 
        coordinator(ChunkDimensionStats<Index_>(mat_nrow, chunk_nrow), ChunkDimensionStats<Index_>(mat_ncol, chunk_ncol), std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / (sizeof(typename Chunk_::value_type) + sizeof(typename Chunk_::index_type))),
        require_minimum_cache(opt.require_minimum_cache)
    {}

private:
    CustomChunkedMatrix_internal::ChunkCoordinator<Index_, true, Chunk_> coordinator;
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
        return true;
    }

    double sparse_proportion() const {
        return 1;
    }

    using tatami::Matrix<Value_, Index_>::dense;

    using tatami::Matrix<Value_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<
        template<bool, typename, typename> class Interface_, 
        bool oracle_, 
        template<bool, bool, bool, typename, typename, class> class Extractor_, 
        typename ... Args_
    >
    std::unique_ptr<Interface_<oracle_, Value_, Index_> > raw_internal(bool row, Index_ secondary_length, Args_&& ... args) const {
        if (row) {
            SlabCacheStats stats(coordinator.get_chunk_nrow(), secondary_length, cache_size_in_elements, require_minimum_cache);
            if (stats.num_slabs_in_cache > 0) {
                return std::make_unique<Extractor_<true, false, oracle_, Value_, Index_, Chunk_> >(coordinator, stats.num_slabs_in_cache, std::forward<Args_>(args)...);
            } else {
                return std::make_unique<Extractor_<true, true, oracle_, Value_, Index_, Chunk_> >(coordinator, stats.num_slabs_in_cache, std::forward<Args_>(args)...);
            }
        } else {
            SlabCacheStats stats(coordinator.get_chunk_ncol(), secondary_length, cache_size_in_elements, require_minimum_cache);
            if (stats.num_slabs_in_cache > 0) {
                return std::make_unique<Extractor_<false, false, oracle_, Value_, Index_, Chunk_> >(coordinator, stats.num_slabs_in_cache, std::forward<Args_>(args)...);
            } else {
                return std::make_unique<Extractor_<false, true, oracle_, Value_, Index_, Chunk_> >(coordinator, stats.num_slabs_in_cache, std::forward<Args_>(args)...);
            }
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, tatami::MaybeOracle<oracle_, Index_> oracle, const tatami::Options& opt) const {
        auto secondary = (row ? coordinator.get_ncol() : coordinator.get_nrow());
        return raw_internal<tatami::DenseExtractor, oracle_, CustomChunkedMatrix_internal::DensifiedFull>(row, secondary, std::move(oracle), opt);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return raw_internal<tatami::DenseExtractor, oracle_, CustomChunkedMatrix_internal::DensifiedBlock>(row, block_length, std::move(oracle), block_start, block_length, opt);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        auto num_indices = indices_ptr->size();
        return raw_internal<tatami::DenseExtractor, oracle_, CustomChunkedMatrix_internal::DensifiedIndex>(row, num_indices, std::move(oracle), std::move(indices_ptr), opt);
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

    /***********************
     *** Oracular dense ***
     ***********************/
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
    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(bool row, tatami::MaybeOracle<oracle_, Index_> oracle, const tatami::Options& opt) const {
        auto secondary = (row ? coordinator.get_ncol() : coordinator.get_nrow());
        return raw_internal<tatami::SparseExtractor, oracle_, CustomChunkedMatrix_internal::SparseFull>(row, secondary, std::move(oracle), opt);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return raw_internal<tatami::SparseExtractor, oracle_, CustomChunkedMatrix_internal::SparseBlock>(row, block_length, std::move(oracle), block_start, block_length, opt);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        auto num_indices = indices_ptr->size();
        return raw_internal<tatami::SparseExtractor, oracle_, CustomChunkedMatrix_internal::SparseIndex>(row, num_indices, std::move(oracle), std::move(indices_ptr), opt);
    }

public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const {
        return sparse_internal<false>(row, false, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return sparse_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        return sparse_internal<false>(row, false, std::move(indices_ptr), opt);
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
        return sparse_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return sparse_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        return sparse_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }
};

}

#endif
