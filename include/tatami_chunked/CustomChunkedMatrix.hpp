#ifndef TATAMI_CHUNKED_CUSTOM_CHUNKED_MATRIX_HPP
#define TATAMI_CHUNKED_CUSTOM_CHUNKED_MATRIX_HPP

#include "tatami/tatami.hpp"
#include "typical_slab_cache.hpp"

#include <vector>

/**
 * @file CustomChunkedMatrix.hpp
 * @brief Custom chunked matrix.
 */

namespace tatami_chunked {

/**
 * @brief Options for custom chunk extraction.
 */
struct CustomChunkedOptions : public TypicalSlabCacheOptions {};

/**
 * @cond
 */
namespace CustomChunkedMatrix_internal {

template<
    typename Index_, 
    class Chunk_, 
    typename ChunkIndex_ = typename Chunk_::index_type // need to put this here as the dense Chunk_ contract doesn't have an index_type.
>
class ChunkCoordinator {
public:
    ChunkCoordinator(Index_ mat_nr, Index_ mat_nc, Index_ chunk_nr, Index_ chunk_nc, std::vector<Chunk_> chunks, bool rm) :
        mat_nrow(mat_nr),
        mat_ncol(mat_nc),
        chunk_nrow(chunk_nr), 
        chunk_ncol(chunk_nc),
        num_chunks_per_row(integer_ceil(mat_ncol, chunk_ncol)),
        num_chunks_per_column(integer_ceil(mat_nrow, chunk_nrow)),
        last_row_chunk_rows(mat_nrow ? (mat_nrow - (num_chunks_per_column - 1) * chunk_nrow) : 0),
        last_column_chunk_cols(mat_ncol ? (mat_ncol - (num_chunks_per_row - 1) * chunk_ncol) : 0),
        chunk_array(std::move(chunks)),
        row_major(rm)
    {
        if (static_cast<size_t>(num_chunks_per_row) * static_cast<size_t>(num_chunks_per_column) != chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    Index_ mat_nrow, mat_ncol;
    Index_ chunk_nrow, chunk_ncol;
    Index_ num_chunks_per_row, num_chunks_per_column;
    Index_ last_row_chunk_rows, last_column_chunk_cols; // number of rows in the chunks overlapping the last row, or columns in the chunks overlapping the last column.

    std::vector<Chunk_> chunk_array;
    bool row_major;

public:
    Index_ get_nrow() const {
        return mat_nrow;
    }

    Index_ get_ncol() const {
        return mat_ncol;
    }

    bool prefer_rows_internal() const {
        // Prefer rows if we have to extract fewer chunks per row.
        return num_chunks_per_column > num_chunks_per_row; 
    }

public:
    template<bool accrow_>
    Index_ get_primary_dim() const {
        if constexpr(accrow_) {
            return mat_nrow;
        } else {
            return mat_ncol;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_dim() const {
        if constexpr(accrow_) {
            return mat_ncol;
        } else {
            return mat_nrow;
        }
    }

    template<bool accrow_>
    Index_ get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_nrow;
        } else {
            return chunk_ncol;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return chunk_ncol;
        } else {
            return chunk_nrow;
        }
    }

    template<bool accrow_>
    Index_ get_primary_num_chunks() const {
        if constexpr(accrow_) {
            return num_chunks_per_column;
        } else {
            return num_chunks_per_row;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_num_chunks() const {
        if constexpr(accrow_) {
            return num_chunks_per_row;
        } else {
            return num_chunks_per_column;
        }
    }

    // Overload that handles the truncated chunk at the bottom/right edges of each matrix.
    template<bool accrow_>
    Index_ get_primary_chunkdim(Index_ chunk_id) const {
        if (chunk_id + 1 == get_primary_num_chunks<accrow_>()) {
            if constexpr(accrow_) {
                return last_row_chunk_rows;
            } else {
                return last_column_chunk_cols;
            }
        } else {
            return get_primary_chunkdim<accrow_>();
        }
    }

private:
    static Index_ integer_ceil(Index_ left, Index_ right) {
        return left / right + (left % right > 0); // avoids overflow.
    }

private:
    typedef typename Chunk_::value_type ChunkValue;
    typedef typename Chunk_::Workspace ChunkWork;

public:
    typedef std::vector<ChunkValue> DenseSlab;

    struct SparseSlab {
        SparseSlab() = default;
        SparseSlab(size_t primary_dim) : values(primary_dim), indices(primary_dim) {}

        std::vector<std::vector<ChunkValue> > values;
        std::vector<std::vector<ChunkIndex_> > indices;

        void resize(size_t primary_dim) {
            values.resize(primary_dim);
            indices.resize(primary_dim);
        }
    };

private:
    template<bool sparse_>
    using Slab = typename std::conditional<sparse_, SparseSlab, DenseSlab>::type;

private:
    template<bool accrow_>
    std::pair<size_t, size_t> offset_and_increment(Index_ chunk_id) const {
        if (row_major) {
            if constexpr(accrow_) {
                return std::pair<size_t, size_t>(static_cast<size_t>(chunk_id) * num_chunks_per_row, 1);
            } else {
                return std::pair<size_t, size_t>(chunk_id, num_chunks_per_row);
            }
        } else {
            if constexpr(accrow_) {
                return std::pair<size_t, size_t>(chunk_id, num_chunks_per_column);
            } else {
                return std::pair<size_t, size_t>(static_cast<size_t>(chunk_id) * num_chunks_per_column, 1);
            }
        }
    }

    template<bool sparse_>
    static auto configure_slab(Slab<sparse_>& slab) {
        if constexpr(!sparse_) {
            return slab.data();
        } else {
            for (auto& x : slab.indices) {
                x.clear();
            }
            for (auto& x : slab.values) {
                x.clear();
            }
            return false;
        }
    }

    template<bool accrow_, bool sparse_, class ExtractFunction_>
    void extract_secondary_block(
        Index_ chunk_id, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab<sparse_>& slab, 
        ExtractFunction_ extract)
    const {
        auto slab_ptr = configure_slab<sparse_>();

        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ start_chunk_index = secondary_block_start / secondary_chunkdim;
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;
        Index_ secondary_block_end = secondary_block_start + secondary_block_length;
        Index_ end_chunk_index = integer_ceil(secondary_block_end, secondary_chunkdim);

        auto oi = offset_and_increment(chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        for (Index_ c = start_chunk_index; c < end_chunk_index; ++c) {
            const auto& chunk = chunk_array[offset];
            Index_ from = (c == start_chunk_index ? secondary_block_start - secondary_start_pos : 0);
            Index_ to = (c + 1 == end_chunk_index ? secondary_block_end - secondary_start_pos : secondary_chunkdim);
            Index_ len = to - from;

            // No need to protect against a zero length, as it should be impossible
            // here (otherwise, start_chunk_index == end_chunk_index and we'd never iterate).
            if constexpr(sparse_) {
                extract(from, len, slab.values, slab.indices, secondary_start_pos);
            } else {
                extract(from, len, slab_ptr);
            }

            secondary_start_pos += len;
            offset += increment;
            if constexpr(!sparse_) {
                slab_ptr += len;
            }
        }
    }

    template<bool accrow_, bool sparse_, class ExtractFunction_>
    void extract_secondary_index(
        Index_ chunk_id, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices, 
        Slab<sparse_>& slab, 
        ExtractFunction_ extract)
    const {
        auto slab_ptr = configure_slab<sparse_>();
        if (secondary_indices.empty()) {
            return;
        }

        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ start_chunk_index = secondary_indices.front() / secondary_chunkdim; // 'secondary_indices' is guaranteed to be non-empty at this point.
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;

        auto oi = offset_and_increment(chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        auto secondary_dim = get_secondary_dim<accrow_>();
        auto iIt = secondary_indices.begin();
        auto iEnd = secondary_indices.end();
        while (iIt != iEnd) {
            const auto& chunk = chunk_array[offset];

            Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // avoid overflow.
            chunk_indices.clear();
            while (iIt != iEnd && *iIt < secondary_end_pos) {
                chunk_indices.push_back(*iIt - secondary_start_pos);
                ++iIt;
            }

            if (!chunk_indices.empty()) {
                if constexpr(sparse_) {
                    extract(chunk, chunk_indices, slab.values, slab.indices, secondary_start_pos);
                } else {
                    extract(chunk, chunk_indices, slab_ptr);
                }
            }

            secondary_start_pos = secondary_end_pos;
            offset += increment;
            if constexpr(!sparse_) {
                slab_ptr += chunk_indices.size();
            }
        }
    }

private:
    // Extract a single element of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_single(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace)
    const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(chunk_offset, from, len, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(chunk_offset, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                }
            );
        }
    }

    // Extract a single element of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_single(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        const std::vector<Index_>& secondary_indices, 
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_indices, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_indices, chunk_workspace, slab_ptr);
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_block(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_block(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        const std::vector<Index_>& secondary_indices, 
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                }
            );
        }
    }

    // Extract an indexed subset of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_index(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                }
            );
        }
    }

    // Extract an indexed subset of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_, bool sparse_>
    void extract_index(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices,
        const std::vector<Index_>& secondary_indices,
        Slab<sparse_>& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex_>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                }
            );
        }
    }

public:
    // Obtain the slab containing the 'i'-th element of the primary dimension.
    // The 'args...' determines whether a contiguous block or indexed subset is
    // used on the secondary dimension.
    template<bool accrow_, bool sparse_, bool oracle_, bool subset_, typename ... Args_>
    std::pair<const Slab<sparse_>*, Index_> fetch(
        Index_ i, 
        typename Chunk_::Workspace& chunk_workspace,
        TypicalSlabCacheWorkspace<oracle_, subset_, Index_, Slab<sparse_> >& cache_workspace,
        Slab<sparse_>& solo,
        Index_ secondary_length,
        const Args_& ... args) 
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t alloc = sparse_ ? primary_chunkdim : primary_chunkdim * static_cast<size_t>(secondary_length); // use size_t to avoid integer overflow.

        if constexpr(!oracle_) {
            auto chunk_id = i / primary_chunkdim;
            auto chunk_offset = i % primary_chunkdim;

            if (cache_workspace.num_slabs_in_cache == 0) {
                extract_single<accrow_, sparse_>(chunk_id, chunk_offset, args..., solo, chunk_workspace);
                return std::make_pair(&solo, static_cast<Index_>(0));

            } else {
                auto& cache = cache_workspace.cache.find(
                    chunk_id,
                    /* create = */ [&]() -> Slab<sparse_> {
                        return Slab<sparse_>(alloc);
                    },
                    /* populate = */ [&](Index_ id, Slab<sparse_>& slab) -> void {
                        extract_block<accrow_, sparse_>(id, 0, get_primary_chunkdim<accrow_>(id), args..., slab, chunk_workspace);
                    }
                );
                return std::make_pair(&cache, chunk_offset);
            }

        } else if constexpr(subset_) {
            auto out = cache_workspace.cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab<sparse_> {
                    return Slab<sparse_>(alloc);
                },
                /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                    for (const auto& p : in_need) {
                        auto ptr = data[p.second];
                        switch (ptr->subset.selection) {
                            case SubsetSelection::FULL:
                                extract_block<accrow_, sparse_>(p.first, 0, get_primary_chunkdim<accrow_>(p.first), args..., ptr->contents, chunk_workspace);
                                break;
                            case SubsetSelection::BLOCK:
                                extract_block<accrow_, sparse_>(p.first, ptr->subset.block_start, ptr->subset.block_length, args..., ptr->contents, chunk_workspace);
                                break;
                            case SubsetSelection::INDEX:
                                extract_index<accrow_, sparse_>(p.first, ptr->subset.indices, args..., ptr->contents, chunk_workspace);
                                break;
                        }
                    }
                }
            );
            return std::make_pair(&(out.first->contents), out.second);

        } else {
            return cache_workspace.cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab<sparse_> {
                    return Slab<sparse_>(alloc);
                },
                /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                    for (const auto& p : in_need) {
                        extract_block<accrow_, sparse_>(p.first, 0, get_primary_chunkdim<accrow_>(p.first), args..., p.second->contents, chunk_workspace);
                    }
                }
            );
        }
    }
};

}
/**
 * @endcond
 */

/**
 * @cond
 */
namespace CustomChunkedMatrix_internal {

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DenseBase {
protected:
    const ChunkCoordinator<Index_, Chunk_, int>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, Chunk_, int>::DenseSlab Slab;
    TypicalSlabCacheWorkspace<oracle_, subset_, Index_, Slab> cache_workspace;
    Slab solo;

public:
    DenseBase(
        const ChunkCoordinator<Index_, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ secondary_length) : 
        coordinator(coordinator),
        cache_workspace(
            coordinator.template get_primary_chunkdim<accrow_>(), 
            secondary_length, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle)
        )
    {
        if (cache_workspace.num_slabs_in_cache == 0) {
            solo.resize(secondary_length);
        }
    }

protected:
    const Value_* process_dense_slab(const std::pair<const Slab*, Index_>& fetched, Value_* buffer, size_t secondary_length) {
        auto ptr = fetched.first->data() + fetched.second * secondary_length; // already cast to size_t to avoid overflow.
        std::copy_n(ptr, secondary_length, buffer);
        return buffer;
    }
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DenseFull : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DenseFull(
        const ChunkCoordinator<Index_, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle) :
        DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            coordinator.template get_secondary_dim<accrow_>()
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto secondary_len = this->coordinator.template get_secondary_dim<accrow_>();
        auto fetched = this->coordinator.template fetch<accrow_, false, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            secondary_len, 
            0, 
            secondary_len
        );
        return this->process_dense_slab(fetched, buffer, secondary_len);
    }
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DenseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DenseBlock(
        const ChunkCoordinator<Index_, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length) : 
        DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
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
        auto fetched = this->coordinator.template fetch<accrow_, false, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            block_length, 
            block_start, 
            block_length
        );
        return this->process_dense_slab(fetched, buffer, block_length);
    }

private:
    Index_ block_start, block_length;
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DenseIndex : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DenseIndex(
        const ChunkCoordinator<Index_, Chunk_, int>& coordinator, 
        size_t cache_size_in_elements, 
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> idx_ptr) :
        DenseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            idx_ptr->size()
        ),
        indices_ptr(std::move(idx_ptr))
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto nidx = indices_ptr->size();
        auto fetched = this->coordinator.template fetch<accrow_, false, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            nidx, 
            *indices_ptr
        );
        return this->process_dense_slab(fetched, buffer, nidx);
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
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
 * @tparam Chunk_ Class of the chunk.
 * @tparam use_subsetted_oracle_ Whether to extract a subset of the primary dimension from each chunk during oracle predictions.
 *
 * Implements a `Matrix` subclass where data is contained in dense rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage;
 * on access, each chunk is decompressed and the desired values are extracted.
 * The `Chunk_` class should provide the following:
 *
 * - A `typedef value_type` specifying the type of the decompressed chunk data.
 * - A nested `Workspace` class that allocates memory for decompression.
 *   This should be default-constructible and may be re-used for decompressing multiple chunks.
 * - A `void extract<accrow_>(Index_ primary, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output) const` method,
 *   which extracts a single primary element from the chunk and stores it in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(Index_ primary, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts a single primary element from the chunk and stores it in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts a contiguous block of primary elements from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(Index_ primary_start, Index_ primary_length, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts a contiguous block of primary elements from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 * 
 * If `use_subsetted_oracle_ = true`, this class will use a `SubsettedOracleCache` to extract subsets of the primary dimension for each chunk when an `Oracle` is supplied.
 * This may improve performance if the chunk is capable of providing optimized access to subsets along the primary dimension.
 * In such cases, we expect the following additional methods:
 *
 * - A `void extract<accrow_>(const std::vector<Index_>& primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts an indexed subset of primary elements from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(const std::vector<Index_>& primary_indices, const std::vector<Index_>& secondary_indices, Workspace& work, Output_* output, size_t stride) const` method,
 *   which extracts an indexed subset of primary elements from the chunk and stores them in `output`.
 *   See `SimpleDenseChunkWrapper::extract()` for more details.
 *
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero -
 * the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<bool use_subsetted_oracle_, typename Value_, typename Index_, typename Chunk_>
class CustomChunkedDenseMatrix : public tatami::Matrix<Value_, Index_> {
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
    CustomChunkedDenseMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        coordinator(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / sizeof(typename Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {}

private:
    CustomChunkedMatrix_internal::ChunkCoordinator<Index_, Chunk_, int> coordinator;
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
public:
    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options&) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::DenseFull<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseFull<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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
            return std::make_unique<CustomChunkedMatrix_internal::DenseBlock<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                block_start, 
                block_length
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseBlock<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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
            return std::make_unique<CustomChunkedMatrix_internal::DenseIndex<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                std::move(indices_ptr)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DenseIndex<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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

/**
 * @cond
 */
namespace CustomChunkedMatrix_internal {

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct SparseBase {
protected:
    const ChunkCoordinator<Index_, Chunk_>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, Chunk_>::SparseSlab Slab;
    TypicalSlabCacheWorkspace<oracle_, subset_, Index_, Slab> cache_workspace;
    bool needs_value, needs_index;
    Slab solo;

public:
    SparseBase(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache,
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ secondary_length, 
        const tatami::Options& opt) : 
        coordinator(coordinator),
        cache_workspace(
            coordinator.template get_primary_chunkdim<accrow_>(), 
            secondary_length, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle)
        ),
        needs_value(opt.sparse_extract_value),
        needs_index(opt.sparse_extract_index)
    {
        if (cache_workspace.num_slabs_in_cache == 0) {
            solo.resize(1);
        }
    }

protected:
    tatami::SparseRange<Value_, Index_> process_sparse_slab(const std::pair<const Slab*, Index_>& fetched, Value_* vbuffer, Index_* ibuffer) {
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
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct SparseFull : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    SparseFull(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options& opt) : 
        SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            coordinator.template get_secondary_dim<accrow_>(),
            opt
        )
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto secondary_len = this->coordinator.template get_secondary_dim<accrow_>();
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            secondary_len, 
            0, 
            secondary_len
        );
        return this->process_sparse_slab(contents, vbuffer, ibuffer); 
    }
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct SparseBlock : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    SparseBlock(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) :
        SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            block_length, 
            opt
        ),
        block_start(block_start),
        block_length(block_length)
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            block_length, 
            block_start, 
            block_length
        );
        return this->process_sparse_slab(contents, vbuffer, ibuffer);
    }

private:
    Index_ block_start, block_length;
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct SparseIndex : public tatami::SparseExtractor<oracle_, Value_, Index_>, public SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    SparseIndex(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> idx_ptr, 
        const tatami::Options& opt) :
        SparseBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache,
            std::move(oracle), 
            idx_ptr->size(),
            opt
        ),
        indices_ptr(std::move(idx_ptr))
    {}

    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            indices_ptr->size(), 
            *indices_ptr
        );
        return this->process_sparse_slab(contents, vbuffer, ibuffer);
    }

private:
    tatami::VectorPtr<Index_> indices_ptr;
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedBase {
protected:
    const ChunkCoordinator<Index_, Chunk_>& coordinator;
    typename Chunk_::Workspace chunk_workspace;
    typedef typename ChunkCoordinator<Index_, Chunk_>::SparseSlab Slab;
    TypicalSlabCacheWorkspace<oracle_, subset_, Index_, Slab> cache_workspace;
    Slab solo;

public:
    DensifiedBase(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ secondary_length) : 
        coordinator(coordinator),
        cache_workspace(
            coordinator.template get_primary_chunkdim<accrow_>(), 
            secondary_length, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle)
        )
    {
        if (cache_workspace.num_slabs_in_cache == 0) {
            solo.resize(1);
        }
    }
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedFull : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DensifiedFull(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle) : 
        DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle),
            coordinator.template get_secondary_dim<accrow_>()
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto secondary_length = this->coordinator.template get_secondary_dim<accrow_>();
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            secondary_length, 
            0, 
            secondary_length
        );

        const auto& values = contents.first->values[contents.second];
        const auto& indices = contents.first->indices[contents.second];
        std::fill_n(buffer, secondary_length, 0);
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            buffer[indices[i]] = values[i];
        }
        return buffer;
    }

private:
    Index_ secondary_length;
};

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedBlock : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DensifiedBlock(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start, 
        Index_ block_length) : 
        DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
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
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            block_length, 
            block_start, 
            block_length
        );

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

template<bool accrow_, bool oracle_, bool subset_, typename Value_, typename Index_, typename Chunk_>
struct DensifiedIndex : public tatami::DenseExtractor<oracle_, Value_, Index_>, public DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_> {
    DensifiedIndex(
        const ChunkCoordinator<Index_, Chunk_>& coordinator, 
        size_t cache_size_in_elements,
        bool require_minimum_cache, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> idx_ptr) :
        DensifiedBase<accrow_, oracle_, subset_, Value_, Index_, Chunk_>(
            coordinator, 
            cache_size_in_elements, 
            require_minimum_cache, 
            std::move(oracle), 
            idx_ptr->size()
        ),
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
        auto contents = this->coordinator.template fetch<accrow_, true, oracle_, subset_>(
            i, 
            this->chunk_workspace, 
            this->cache_workspace, 
            this->solo, 
            nidx, 
            *indices_ptr
        );

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
 * @tparam Chunk_ Class of the chunk.
 * @tparam use_subsetted_oracle_ Whether to report the subset of each chunk during oracle predictions, see `SubsettedOracleCache` for details.
 *
 * Implements a `Matrix` subclass where data is contained in sparse rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage;
 * on access, each chunk is decompressed and the desired values are extracted.
 * The `Chunk_` class should provide the following:
 *
 * - A `typedef value_type` specifying the type of the decompressed chunk data.
 * - A `typedef index_type` specifying the type of the indices in the decompressed chunk.
 * - A nested `Workspace` class that allocates memory for decompression.
 *   This should be default-constructible and may be re-used for decompressing multiple chunks.
 * - A `void extract<accrow_>(Index_ primary, Index_ secondary_start, Index_ secondary_length, Workspace& work, std::vector<value_type>& output_values, std::vector<index_type>& output_indices) const` method,
 *   which extracts a single primary element from the chunk and stores it in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(Index_ primary, const std::vector<Index_>& secondary_indices, Workspace& work, std::vector<value_type>& output_values, std::vector<index_type>& output_indices) const` method,
 *   which extracts a single primary element from the chunk and stores it in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 * - A `void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, 
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts contiguous ranges from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 * - A `void extract(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace& work, 
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts an indexed subset from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 *
 * If `use_subsetted_oracle_ = true`, this class will use a `SubsettedOracleCache` to extract subsets of the primary dimension for each chunk when an `Oracle` is supplied.
 * This may improve performance if the chunk is capable of providing optimized access to subsets along the primary dimension.
 * In such cases, we expect the following additional methods:
 *
 * - A `void extract<accrow_>(const std::vector<Index_>& primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace& work,
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts an indexed subset of primary elements from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 * - A `void extract<accrow_>(const std::vector<Index_>& primary_indices, const std::vector<Index_>& secondary_indices, Workspace& work,
 *   std::vector<std::vector<value_type> >& output_values, std::vector<std::vector<index_type> >& output_indices) const` method,
 *   which extracts an indexed subset of primary elements from the chunk and stores them in `output_values` and `output_indices`.
 *   See `SimpleSparseChunkWrapper::extract()` for more details.
 *
 * All chunks should have the same dimensions, i.e., covering the same shape/area of the matrix.
 * The matrix should be partitioned at regular intervals starting from zero -
 * the first chunk should start at (0, 0), the next chunk should be immediately adjacent in one of the dimensions, and so on.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<bool use_subsetted_oracle_, typename Value_, typename Index_, typename Chunk_>
class CustomChunkedSparseMatrix : public tatami::Matrix<Value_, Index_> {
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
    CustomChunkedSparseMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomChunkedOptions& opt) : 
        coordinator(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / (sizeof(typename Chunk_::value_type) + sizeof(typename Chunk_::index_type))),
        require_minimum_cache(opt.require_minimum_cache)
    {}

private:
    CustomChunkedMatrix_internal::ChunkCoordinator<Index_, Chunk_> coordinator;
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
public:
    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options&) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedFull<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedFull<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedBlock<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                block_start, 
                block_length
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedBlock<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedIndex<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                std::move(indices_ptr)
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::DensifiedIndex<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
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
public:
    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options& opt) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::SparseFull<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                opt
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::SparseFull<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                opt
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::SparseBlock<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                block_start, 
                block_length, 
                opt
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::SparseBlock<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle), 
                block_start, 
                block_length, 
                opt
            );
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > sparse_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        if (row) {
            return std::make_unique<CustomChunkedMatrix_internal::SparseIndex<true, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle),
                std::move(indices_ptr), 
                opt
            );
        } else {
            return std::make_unique<CustomChunkedMatrix_internal::SparseIndex<false, oracle_, use_subsetted_oracle_, Value_, Index_, Chunk_> >(
                coordinator,
                cache_size_in_elements,
                require_minimum_cache,
                std::move(oracle),
                std::move(indices_ptr), 
                opt
            );
        }
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
