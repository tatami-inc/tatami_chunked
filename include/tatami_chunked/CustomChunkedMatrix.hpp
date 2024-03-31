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

template<typename Index_, bool sparse_, class Chunk_>
class ChunkCoordinator {
protected:
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
    {}

protected:
    Index_ mat_nrow, mat_ncol;
    Index_ chunk_nrow, chunk_ncol;
    Index_ num_chunks_per_row, num_chunks_per_column;
    Index_ leftover_chunks_per_row, leftover_chunks_per_column;

    std::vector<Chunk_> chunk_array;
    bool row_major;

    bool prefer_rows_internal() const {
        // Prefer rows if we have to extract fewer chunks per row.
        return num_chunks_per_column > num_chunks_per_row; 
    }

protected:
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

    template<bool accrow_>
    Index_ get_primary_chunkdim(Index_ chunk_id) {
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

protected:
    static Index_ integer_ceil(Index_ left, Index_ right) {
        return left / right + (left % right > 0); // avoids overflow.
    }

protected:
    typedef typename Chunk_::value_type ChunkValue;
    typedef typename Chunk_::index_type ChunkIndex;

    typedef std::vector<ChunkValue> DenseSlab;

    struct SparseSlab {
        SparseSlab() = default;
        SparseSlab(size_t primary_dim) : values(primary_dim), indices(primary_dim) {}

        std::vector<std::vector<ChunkValue> > values;
        std::vector<std::vector<ChunkIndex> > indices;

        void resize(size_t primary_dim) {
            values.resize(primary_dim);
            indices.resize(primary_dim);
        }
    };

    typedef typename std::conditional<sparse_, SparseSlab, DenseSlab>::type Slab;

protected:
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

    static auto configure_slab(Slab& slab) {
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

    template<bool accrow_, class ExtractFunction_>
    void extract_secondary_block(
        Index_ chunk_id, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        ExtractFunction_ extract)
    const {
        auto slab_ptr = configure_slab();

        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ start_chunk_index = secondary_block_start / secondary_chunkdim;
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;
        Index_ secondary_block_end = secondary_block_start + secondary_block_length;
        Index_ end_chunk_index = integer_ceil(block_end, secondary_chunkdim);

        auto oi = offset_and_increment(chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        for (Index_ c = start_chunk_index; c < end_chunk_index; ++c) {
            const auto& chunk = chunk_array[offset];
            Index_ from = (c == start_chunk_index ? secondary_block_start - secondary_start_pos : 0);
            Index_ to = (c + 1 == end_chunk_index ? secondary_block_end - secondary_start_pos : secondary_chunkdim);

            // No need to protect against a zero length, as it should be impossible
            // here (otherwise, start_chunk_index == end_chunk_index and we'd never iterate).
            if constexpr(sparse_) {
                extract(from, to - from, slab.values, slab.indices, secondary_start_pos);
            } else {
                extract(from, to - from, slab_ptr);
            }

            secondary_start_pos += secondary_chunkdim;
            offset += increment;
            if constexpr(!sparse_) {
                slab_ptr += (to - from);
            }
        }
    }

    template<bool accrow_, class ExtractFunction_>
    void extract_secondary_index(
        Index_ chunk_id, 
        const std::vector<Index_>& secondary_indices, 
        typename Chunk_::Workspace& workspace, 
        std::vector<Index_>& chunk_indices, 
        Slab& slab, 
        ExtractFunction_ extract)
    const {
        auto slab_ptr = configure_slab();
        if (secondary_indices.empty()) {
            return;
        }

        auto secondary_num_chunks = get_secondary_num_chunks<accrow_>();
        Index_ start_chunk_index = ext->indices.front() / secondary_chunkdim; // 'indices' is guaranteed to be non-empty at this point.
        offset += static_cast<size_t>(start_chunk_index) * increment; // use size_t to avoid integer overflow.
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;

        auto oi = offset_and_increment(chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        auto iIt = secondary_indices.begin();
        auto iEnd = secondary_indices.end();
        for (Index_ c = start_chunk_index; c < secondary_num_chunks; ++c) {
            const auto& chunk = chunk_array[offset];

            Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // avoid overflow.
            chunk_indices.clear();
            while (iIt != iEnd && *iIt < secondary_end_pos) {
                chunk_indices.push_back(*iIt - secondary_start_pos);
                ++iIt;
            }

            if (!chunk_indices.empty()) {
                if constexpr(sparse_) {
                    extract(chunk, from, to - from, slab.values, slab.indices, secondary_start_pos);
                } else {
                    extract(chunk, from, to - from, slab_ptr);
                }
            }
        }
    }

protected:
    // Extract a single element of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_>
    void extract_single_and_block(Index_ chunk_id, Index_ chunk_offset, Index_ secondary_block_start, Index_ secondary_block_length, Slab& slab, typename Chunk_::Workspace& chunk_workspace) const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
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
    template<bool accrow_>
    void extract_single_and_index(Index_ chunk_id, Index_ chunk_offset, const std::vector<Index_>& secondary_indices, Slab& slab, typename Chunk_::Workspace& chunk_workspace) const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_indices, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(chunk_offset, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_>
    void extract_block_and_block(
        Index_ chunk_id,
        Index_ primary_block_start, 
        Index_ primary_block_length, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        typename Chunk_::Workspace& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(primary_block_start, primary_block_length, from, len, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(primary_block_start, primary_block_length, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_>
    void extract_block_and_index(
        Index_ chunk_id,
        Index_ primary_block_start, 
        Index_ primary_block_length, 
        const std::vector<Index_>& secondary_indices,
        Slab& slab, 
        typename Chunk_::Workspace& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
                    chunk.template extract<accrow_>(primary_block_start, primary_block_length, chunk_indices, chunk_workspace, slab_values, slab_indices, secondary_start_pos);
                }
            );
        } else {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, ChunkValue* slab_ptr) {
                    chunk.template extract<accrow_>(primary_block_start, primary_block_length, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                }
            );
        }
    }

    // Extract an indexed subset of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_>
    void extract_index_and_block(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        typename Chunk_::Workspace& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, slab,
                [&](const Chunk_& chunk, Index_ from, Index_ len, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
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
    template<bool accrow_>
    void extract_index_and_index(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices,
        const std::vector<Index_>& secondary_indices,
        Slab& slab, 
        typename Chunk_::Workspace& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, slab,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, std::vector<ChunkValue>& slab_values, std::vector<ChunkIndex>& slab_indices, Index_ secondary_start_pos) {
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
    // Obtain the slab containing the 'i'-th element of the primary dimension and a contiguous block on the secondary dimension.
    template<bool accrow_, bool use_subsetted_oracle_>
    std::pair<const Slab*, Index_> fetch_block(Index_ i, Index_ block_start, Index_ block_length, typename Chunk_::Workspace& chunk_workspace) const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t alloc = sparse_ ? primary_chunkdim : primary_chunkdim * static_cast<size_t>(block_length); // use size_t to avoid integer overflow.

        if (cache_workspace.oracle_cache) {
            if constexpr(use_subsetted_oracle_) {
                auto out = cache_workspace.oracle_cache->next(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                    },
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                        for (const auto& p : in_need) {
                            auto ptr = data[p.second];
                            switch (ptr->subset.selection) {
                                case SubsetSelection::FULL:
                                    extract_block_and_block(p.first, 0, get_primary_chunkdim<accrow_>(p.first), block_start, block_length, ptr->contents, chunk_workspace);
                                    break;
                                case SubsetSelection::BLOCK:
                                    extract_block_and_block(p.first, ptr->subset.block_start, ptr->subset.block_length, block_start, block_length, ptr->contents, chunk_workspace);
                                    break;
                                case SubsetSelection::INDEX:
                                    extract_index_and_block(p.first, ptr->subset.indices, block_start, block_length, ptr->contents, chunk_workspace);
                                    break;
                            }
                        }
                    }
                );
                return std::make_pair(&(out.first->contents), out.second);

            } else {
                return cache_workspace.oracle_cache->next(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                    },
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                        for (const auto& p : in_need) {
                            extract_block_and_block(p.first, 0, get_primary_chunkdim<accrow_>(p.first), block_start, block_length, p.second->contents, chunk_workspace);
                        }
                    }
                );
            }

        } else {
            auto chunk_id = i / primary_chunkdim;
            auto chunk_offset = i % primary_chunkdim;

            if (cache_workspace.num_slabs_in_cache == 0) {
                extract_single_and_block(chunk_id, chunk_offset, block-start, block_length, ptr->contents, chunk_workspace);
                return std::make_pair(&(ext->solo), static_cast<Index_>(0));

            } else {
                auto& cache = cache_workspace.lru_cache->find(
                    chunk_id,
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate = */ [&](Index_ id, Slab& slab) -> void {
                        extract_block_and_block(id, 0, get_primary_chunkdim<accrow_>(id), block_start, block_length, slab, chunk_workspace);
                    }
                );
                return std::make_pair(&cache, chunk_offset);
            }
        }
    }

    // Obtain the slab containing the 'i'-th element of the primary dimension and an indexed subset on the secondary dimension.
    template<bool accrow_, bool use_subsetted_oracle_>
    std::pair<const Slab*, Index_> fetch_block(Index_ i, const std::vector<Index_>& indices, typename Chunk_::Workspace& chunk_workspace) const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t alloc = sparse_ ? primary_chunkdim : primary_chunkdim * static_cast<size_t>(indices.size()); // use size_t to avoid integer overflow.

        if (cache_workspace.oracle_cache) {
            if constexpr(use_subsetted_oracle_) {
                auto out = cache_workspace.oracle_cache->next(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                    },
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                        for (const auto& p : in_need) {
                            auto ptr = data[p.second];
                            switch (ptr->subset.selection) {
                                case SubsetSelection::FULL:
                                    extract_block_and_index(p.first, 0, get_primary_chunkdim<accrow_>(p.first), indices, ptr->contents, chunk_workspace);
                                    break;
                                case SubsetSelection::BLOCK:
                                    extract_block_and_index(p.first, ptr->subset.block_start, ptr->subset.block_length, indices, ptr->contents, chunk_workspace);
                                    break;
                                case SubsetSelection::INDEX:
                                    extract_index_and_index(p.first, ptr->subset.indices, indices, ptr->contents, chunk_workspace);
                                    break;
                            }
                        }
                    }
                );
                return std::make_pair(&(out.first->contents), out.second);

            } else {
                return cache_workspace.oracle_cache->next(
                    /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                        return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                    },
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                        for (const auto& p : in_need) {
                            extract_block_and_index(p.first, 0, get_primary_chunkdim<accrow_>(p.first), indices, p.second->contents, chunk_workspace);
                        }
                    }
                );
            }

        } else {
            auto chunk_id = i / primary_chunkdim;
            auto chunk_offset = i % primary_chunkdim;

            if (cache_workspace.num_slabs_in_cache == 0) {
                extract_single_and_index(chunk_id, chunk_offset, indices, ptr->contents, chunk_workspace);
                return std::make_pair(&(ext->solo), static_cast<Index_>(0));

            } else {
                auto& cache = cache_workspace.lru_cache->find(
                    chunk_id,
                    /* create = */ [&]() -> Slab {
                        return Slab(alloc);
                    },
                    /* populate = */ [&](Index_ id, Slab& slab) -> void {
                        extract_block_and_index(id, 0, get_primary_chunkdim<accrow_>(id), indices, slab, chunk_workspace);
                    }
                );
                return std::make_pair(&cache, chunk_offset);
            }
        }
    }
};

template<bool oracle_, typename Value_, typename Index_>
struct Dense : public tatami::DenseExtractor<oracle_, Value_, Index_> {
    CustomExtractor(const CustomChunkedDenseMatrix* p) : parent(p) {
        auto len = tatami::extracted_length<selection_, Index_>(*this);

        cache_workspace = TypicalSlabCacheWorkspace<Index_, Slab, use_subsetted_oracle_>(
            accrow_ ? parent->chunk_nrow : parent->chunk_ncol,
            len,
            parent->cache_size_in_elements,
            parent->require_minimum_cache
        );

        if (cache_workspace.num_slabs_in_cache == 0) {
            solo.resize(len);
        }
    }

    CustomExtractor(const CustomChunkedDenseMatrix* p, Index_ bs, Index_ bl) : parent(p) {
        if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
            this->block_start = bs;
            this->block_length = bl;
        }
        initialize_cache();
    }

    CustomExtractor(const CustomChunkedDenseMatrix* p, std::vector<Index_> idx) : parent(p) {
        if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
            indices = std::move(idx);
            this->index_length = indices.size();
        }
        initialize_cache();
    }

private:
    const CustomChunkedDenseMatrix* parent;

    typedef typename CustomChunkedMatrixMethods<Index_, false, Chunk_>::Slab Slab;
    typename Chunk_::Workspace chunk_workspace;
    TypicalSlabCacheWorkspace<Index_, Slab, use_subsetted_oracle_> cache_workspace;
    Slab solo;

    void initialize_cache() {
        }
    }

public:
    friend class CustomChunkedMatrixMethods<Index_, false, Chunk_>;

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = parent->template fetch_cache<accrow_, selection_, use_subsetted_oracle_>(i, this);
        size_t len = tatami::extracted_length<selection_, Index_>(*this); // size_t to avoid overflow.
        auto ptr = fetched.first->data() + fetched.second * len;
        std::copy(ptr, ptr + len, buffer);
        return buffer;
    }
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
template<bool use_subsetted_oracle_ = false, typename Value_, typename Index_, typename Chunk_>
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
        CustomChunkedMatrixMethods<Index_, false, Chunk_>(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / sizeof(typename Chunk_::value_type)),
        require_minimum_cache(opt.require_minimum_cache)
    {
        if (static_cast<size_t>(this->num_chunks_per_column) * static_cast<size_t>(this->num_chunks_per_row) != this->chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    size_t cache_size_in_elements;
    bool require_minimum_cache;

public:
    Index_ nrow() const { return this->mat_nrow; }

    Index_ ncol() const { return this->mat_ncol; }

    bool prefer_rows() const { 
        return this->prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(this->prefer_rows_internal());
    }

    using tatami::Matrix<Value_, Index_>::dense_row;

    using tatami::Matrix<Value_, Index_>::dense_column;

    using tatami::Matrix<Value_, Index_>::sparse_row;

    using tatami::Matrix<Value_, Index_>::sparse_column;

private:
public:
    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_row(const tatami::Options& ) const {
        auto ptr = new CustomExtractor<true, tatami::DimensionSelectionType::FULL>(this);
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const tatami::Options& ) const {
        auto ptr = new CustomExtractor<true, tatami::DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const tatami::Options& ) const {
        auto ptr = new CustomExtractor<true, tatami::DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_column(const tatami::Options& ) const {
        auto ptr = new CustomExtractor<false, tatami::DimensionSelectionType::FULL>(this);
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const tatami::Options& ) const {
        auto ptr = new CustomExtractor<false, tatami::DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const tatami::Options& ) const {
        auto ptr = new CustomExtractor<false, tatami::DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(ptr);
    }
};

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
template<typename Value_, typename Index_, typename Chunk_, bool use_subsetted_oracle_ = false>
class CustomChunkedSparseMatrix : public tatami::Matrix<Value_, Index_>, public CustomChunkedMatrixMethods<Index_, true, Chunk_> {
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
        CustomChunkedMatrixMethods<Index_, true, Chunk_>(mat_nrow, mat_ncol, chunk_nrow, chunk_ncol, std::move(chunks), row_major),
        cache_size_in_elements(opt.maximum_cache_size / (sizeof(typename Chunk_::value_type) + sizeof(typename Chunk_::index_type))),
        require_minimum_cache(opt.require_minimum_cache)
    {
        if (static_cast<size_t>(this->num_chunks_per_column) * static_cast<size_t>(this->num_chunks_per_row) != this->chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    size_t cache_size_in_elements;
    bool require_minimum_cache;

public:
    Index_ nrow() const { return this->mat_nrow; }

    Index_ ncol() const { return this->mat_ncol; }

    bool prefer_rows() const { 
        return this->prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(this->prefer_rows_internal());
    }

    bool sparse() const {
        return true;
    }

    double sparse_proportion() const {
        return 1;
    }

    using tatami::Matrix<Value_, Index_>::dense_row;

    using tatami::Matrix<Value_, Index_>::dense_column;

    using tatami::Matrix<Value_, Index_>::sparse_row;

    using tatami::Matrix<Value_, Index_>::sparse_column;

private:
    template<bool accrow_, bool sparse_, tatami::DimensionSelectionType selection_>
    struct CustomExtractorBase : public tatami::Extractor<selection_, sparse_, Value_, Index_> {
        CustomExtractorBase(const CustomChunkedSparseMatrix* p) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
                this->full_length = parent->template get_secondary_dim<accrow_>();
            }
            initialize_cache();
        }

        CustomExtractorBase(const CustomChunkedSparseMatrix* p, Index_ bs, Index_ bl) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                this->block_start = bs;
                this->block_length = bl;
            }
            initialize_cache();
        }

        CustomExtractorBase(const CustomChunkedSparseMatrix* p, std::vector<Index_> idx) : parent(p) {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
            }
            initialize_cache();
        }

    protected:
        const CustomChunkedSparseMatrix* parent;
        typename std::conditional<selection_ == tatami::DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type indices;
        typename std::conditional<selection_ == tatami::DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type chunk_indices;

        typedef typename CustomChunkedMatrixMethods<Index_, true, Chunk_>::Slab Slab;
        typename Chunk_::Workspace chunk_workspace;
        TypicalSlabCacheWorkspace<Index_, Slab, use_subsetted_oracle_> cache_workspace;
        Slab solo;

        void initialize_cache() {
            cache_workspace = TypicalSlabCacheWorkspace<Index_, Slab, use_subsetted_oracle_>(
                accrow_ ? parent->chunk_nrow : parent->chunk_ncol,
                tatami::extracted_length<selection_, Index_>(*this),
                parent->cache_size_in_elements,
                parent->require_minimum_cache
            );

            if (cache_workspace.num_slabs_in_cache == 0) {
                solo.resize(1);
            }
        }

    public:
        const Index_* index_start() const {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<tatami::Oracle<Index_> > o) {
            cache_workspace.set_oracle(std::move(o));
            return;
        }

        friend class CustomChunkedMatrixMethods<Index_, true, Chunk_>;
    };

private:
    template<bool accrow_, tatami::DimensionSelectionType selection_>
    struct CustomDenseExtractor : public CustomExtractorBase<accrow_, false, selection_> {
    public:
        template<typename ... Args_>
        CustomDenseExtractor(const CustomChunkedSparseMatrix* p, Args_&& ... args) : CustomExtractorBase<accrow_, false, selection_>(p, std::forward<Args_>(args)...) {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                const auto& indices = this->indices;
                if (!indices.empty()) {
                    remap.resize(indices.back() + 1);
                    for (Index_ i = 0, end = indices.size(); i < end; ++i) {
                        remap[indices[i]] = i;
                    }
                }
            }
        }

    public:
        const Value_* fetch(Index_ i, Value_* buffer) {
            auto contents = this->parent->template fetch_cache<accrow_, selection_, use_subsetted_oracle_>(i, this);
            const auto& values = contents.first->values[contents.second];
            const auto& indices = contents.first->indices[contents.second];

            auto len = tatami::extracted_length<selection_, Index_>(*this);
            std::fill(buffer, buffer + len, 0);

            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    buffer[remap[indices[i]]] = values[i];
                }
            } else {
                Index_ start = 0;
                if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                    start = this->block_start;
                }

                for (size_t i = 0, end = indices.size(); i < end; ++i) {
                    buffer[indices[i] - start] = values[i];
                }
            }

            return buffer;
        }

    private:
        typename std::conditional<selection_ == tatami::DimensionSelectionType::INDEX, std::vector<Index_>, bool>::type remap;
    };

public:
    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_row(const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<true, tatami::DimensionSelectionType::FULL>(this);
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_row(Index_ block_start, Index_ block_length, const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<true, tatami::DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_row(std::vector<Index_> indices, const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<true, tatami::DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> > dense_column(const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<false, tatami::DimensionSelectionType::FULL>(this);
        return std::unique_ptr<tatami::FullDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> > dense_column(Index_ block_start, Index_ block_length, const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<false, tatami::DimensionSelectionType::BLOCK>(this, block_start, block_length);
        return std::unique_ptr<tatami::BlockDenseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> > dense_column(std::vector<Index_> indices, const tatami::Options&) const {
        auto ptr = new CustomDenseExtractor<false, tatami::DimensionSelectionType::INDEX>(this, std::move(indices));
        return std::unique_ptr<tatami::IndexDenseExtractor<Value_, Index_> >(ptr);
    }

private:
    template<bool accrow_, tatami::DimensionSelectionType selection_>
    struct CustomSparseExtractor : public CustomExtractorBase<accrow_, true, selection_> {
    public:
        template<typename ... Args_>
        CustomSparseExtractor(const CustomChunkedSparseMatrix* p, bool ev, bool ei, Args_&& ... args) : 
            CustomExtractorBase<accrow_, true, selection_>(p, std::forward<Args_>(args)...), 
            extract_value(ev), 
            extract_index(ei)
        {}

    public:
        tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* vbuffer, Index_* ibuffer) {
            auto contents = this->parent->template fetch_cache<accrow_, selection_, use_subsetted_oracle_>(i, this);
            const auto& values = contents.first->values[contents.second];
            const auto& indices = contents.first->indices[contents.second];

            if (extract_value) {
                std::copy(values.begin(), values.end(), vbuffer);
            } else {
                vbuffer = NULL;
            }

            if (extract_index) {
                std::copy(indices.begin(), indices.end(), ibuffer);
            } else {
                ibuffer = NULL;
            }

            return tatami::SparseRange<Value_, Index_>(values.size(), vbuffer, ibuffer);
        }

    private:
        bool extract_value, extract_index;
    };

public:
    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_row(const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, tatami::DimensionSelectionType::FULL>(this, opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_row(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, tatami::DimensionSelectionType::BLOCK>(this, opt.sparse_extract_value, opt.sparse_extract_index, block_start, block_length);
        return std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_row(std::vector<Index_> indices, const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<true, tatami::DimensionSelectionType::INDEX>(this, opt.sparse_extract_value, opt.sparse_extract_index, std::move(indices));
        return std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> > sparse_column(const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, tatami::DimensionSelectionType::FULL>(this, opt.sparse_extract_value, opt.sparse_extract_index);
        return std::unique_ptr<tatami::FullSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> > sparse_column(Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, tatami::DimensionSelectionType::BLOCK>(this, opt.sparse_extract_value, opt.sparse_extract_index, block_start, block_length);
        return std::unique_ptr<tatami::BlockSparseExtractor<Value_, Index_> >(ptr);
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> > sparse_column(std::vector<Index_> indices, const tatami::Options& opt) const {
        auto ptr = new CustomSparseExtractor<false, tatami::DimensionSelectionType::INDEX>(this, opt.sparse_extract_value, opt.sparse_extract_index, std::move(indices));
        return std::unique_ptr<tatami::IndexSparseExtractor<Value_, Index_> >(ptr);
    }
};

}

#endif
