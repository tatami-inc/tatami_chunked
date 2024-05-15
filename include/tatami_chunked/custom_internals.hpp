#ifndef TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP
#define TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP

#include "tatami/tatami.hpp"
#include "OracularSlabCache.hpp"
#include "OracularSubsettedSlabCache.hpp"
#include "ChunkDimensionStats.hpp"

#include <vector>

namespace tatami_chunked {

namespace CustomChunkedMatrix_internal {

/*******************
 *** Dense slabs ***
 *******************/

template<typename CachedValue_>
struct DenseSlabFactory {
    DenseSlabFactory(size_t primary_chunkdim, size_t secondary_length, size_t num_slabs) : 
        slab_size(primary_chunkdim * secondary_length), 
        pool(num_slabs * slab_size) 
    {}

    // Delete the copy constructors as we're passing out pointers.
    DenseSlabFactory(const DenseSlabFactory&) = delete;
    DenseSlabFactory& operator=(const DenseSlabFactory&) = delete;

    // Move constructors are okay though.
    DenseSlabFactory(DenseSlabFactory&&) = default;
    DenseSlabFactory& operator=(DenseSlabFactory&&) = default;

private:
    size_t offset = 0, slab_size;
    std::vector<CachedValue_> pool;

public:
    struct Slab {
        CachedValue_* data;
    };

    Slab create() {
        Slab output;
        output.data = pool.data() + offset;
        offset += slab_size;
        return output;
    }
};

template<typename CachedValue_>
using DenseSingleWorkspace = std::vector<CachedValue_>;

/********************
 *** Sparse slabs ***
 ********************/

template<typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseSlabFactory {
    SparseSlabFactory(size_t primary_chunkdim, size_t secondary_length, size_t num_slabs, bool needs_value, bool needs_index) : 
        primary_chunkdim(primary_chunkdim),
        secondary_length(secondary_length),
        slab_size(primary_chunkdim * secondary_length), 
        needs_value(needs_value),
        needs_index(needs_index),
        number_pool(num_slabs * primary_chunkdim)
    {
        size_t total_size = num_slabs * slab_size;
        if (needs_value) {
            value_pool.resize(total_size);
        }
        if (needs_index) {
            index_pool.resize(total_size);
        }
    }

    // Delete the copy constructors as we're passing out pointers.
    SparseSlabFactory(const SparseSlabFactory&) = delete;
    SparseSlabFactory& operator=(const SparseSlabFactory&) = delete;

    // Move constructors are okay though.
    SparseSlabFactory(SparseSlabFactory&&) = default;
    SparseSlabFactory& operator=(SparseSlabFactory&&) = default;

private:
    size_t offset_slab = 0, offset_number = 0;
    size_t primary_chunkdim, secondary_length, slab_size;
    bool needs_value, needs_index;
    std::vector<CachedValue_> value_pool;
    std::vector<CachedIndex_> index_pool;
    std::vector<Index_> number_pool;

public:
    struct Slab {
        Index_* number;
        std::vector<CachedValue_*> values;
        std::vector<CachedIndex_*> indices;
    };

    Slab create() {
        Slab output;
        output.number = number_pool.data() + offset_number;
        offset_number += primary_chunkdim;

        if (needs_value) {
            output.values.reserve(primary_chunkdim);
            auto vptr = value_pool.data() + offset_slab;
            for (size_t p = 0; p < primary_chunkdim; ++p, vptr += secondary_length) {
                output.values.push_back(vptr);
            }
        }

        if (needs_index) {
            output.indices.reserve(primary_chunkdim);
            auto iptr = index_pool.data() + offset_slab;
            for (size_t p = 0; p < primary_chunkdim; ++p, iptr += secondary_length) {
                output.indices.push_back(iptr);
            }
        }

        offset_slab += slab_size;
        return output;
    }
};

template<typename Index_, typename CachedValue_, typename CachedIndex_>
struct SparseSingleWorkspace {
    SparseSingleWorkspace(size_t primary_chunkdim, size_t secondary_chunkdim, bool needs_value, bool needs_index) : number(primary_chunkdim) {
        size_t total_size = primary_chunkdim * secondary_chunkdim;
        if (needs_value) {
            value_pool.resize(total_size);
            values.reserve(primary_chunkdim);
            auto vptr = value_pool.data();
            for (size_t p = 0; p < primary_chunkdim; ++p, vptr += secondary_chunkdim) {
                values.push_back(vptr);
            }
        }
        if (needs_index) {
            index_pool.resize(total_size);
            indices.reserve(primary_chunkdim);
            auto iptr = index_pool.data();
            for (size_t p = 0; p < primary_chunkdim; ++p, iptr += secondary_chunkdim) {
                indices.push_back(iptr);
            }
        }
    }

    // Delete the copy constructors as we're passing out pointers.
    SparseSingleWorkspace(const SparseSingleWorkspace&) = delete;
    SparseSingleWorkspace& operator=(const SparseSingleWorkspace&) = delete;

    // Move constructors are okay though.
    SparseSingleWorkspace(SparseSingleWorkspace&&) = default;
    SparseSingleWorkspace& operator=(SparseSingleWorkspace&&) = default;

private:
    std::vector<CachedValue_> value_pool;
    std::vector<CachedIndex_> index_pool;

public:
    std::vector<CachedValue_*> values;
    std::vector<CachedIndex_*> indices;
    std::vector<Index_> number;
};

/*******************
 *** Coordinator ***
 *******************/

template<
    typename Index_, 
    bool sparse_,   
    class Chunk_, 
    typename ChunkIndex_ = typename Chunk_::index_type // need to put this here as the dense Chunk_ contract doesn't have an index_type.
>
class ChunkCoordinator {
public:
    ChunkCoordinator(ChunkDimensionStats<Index_> rows, ChunkDimensionStats<Index_> cols, std::vector<Chunk_> chunks, bool rm) :
        row_stats(std::move(rows)), col_stats(std::move(cols)), chunk_array(std::move(chunks)), row_major(rm)
    {
        if (static_cast<size_t>(row_stats.num_chunks) * static_cast<size_t>(col_stats.num_chunks) != chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    ChunkDimensionStats<Index_> row_stats;
    ChunkDimensionStats<Index_> col_stats;
    std::vector<Chunk_> chunk_array;
    bool row_major;

public:
    // Number of chunks along the rows is equal to the number of chunks for
    // each column, and vice versa; hence the flipped definitions.
    Index_ get_num_chunks_per_row() const {
        return col_stats.num_chunks;
    }

    Index_ get_num_chunks_per_column() const {
        return row_stats.num_chunks;
    }

    Index_ get_nrow() const {
        return row_stats.dimension_extent;
    }

    Index_ get_ncol() const {
        return col_stats.dimension_extent;
    }

    bool prefer_rows_internal() const {
        // Prefer rows if we have to extract fewer chunks per row.
        return get_num_chunks_per_column() > get_num_chunks_per_row(); 
    }

    Index_ get_chunk_nrow() const {
        return row_stats.chunk_length;
    }

    Index_ get_chunk_ncol() const {
        return col_stats.chunk_length;
    }

public:
    template<bool accrow_>
    Index_ get_primary_dim() const {
        if constexpr(accrow_) {
            return row_stats.dimension_extent;
        } else {
            return col_stats.dimension_extent;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_dim() const {
        if constexpr(accrow_) {
            return col_stats.dimension_extent;
        } else {
            return row_stats.dimension_extent;
        }
    }

    template<bool accrow_>
    Index_ get_primary_chunkdim() const {
        if constexpr(accrow_) {
            return row_stats.chunk_length;
        } else {
            return col_stats.chunk_length;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_chunkdim() const {
        if constexpr(accrow_) {
            return col_stats.chunk_length;
        } else {
            return row_stats.chunk_length;
        }
    }

    template<bool accrow_>
    Index_ get_primary_num_chunks() const {
        if constexpr(accrow_) {
            return row_stats.num_chunks;
        } else {
            return col_stats.num_chunks;
        }
    }

    template<bool accrow_>
    Index_ get_secondary_num_chunks() const {
        if constexpr(accrow_) {
            return col_stats.num_chunks;
        } else {
            return row_stats.num_chunks;
        }
    }

    // Overload that handles the truncated chunk at the bottom/right edges of each matrix.
    template<bool accrow_>
    Index_ get_primary_chunkdim(Index_ chunk_id) const {
        if constexpr(accrow_) {
            return row_stats.get_chunk_length(chunk_id);
        } else {
            return col_stats.get_chunk_length(chunk_id);
        }
    }

private:
    template<bool accrow_>
    std::pair<size_t, size_t> offset_and_increment(Index_ chunk_id) const {
        if (row_major) {
            if constexpr(accrow_) {
                return std::pair<size_t, size_t>(static_cast<size_t>(chunk_id) * static_cast<size_t>(get_num_chunks_per_row()), 1); // use size_t to avoid overflow.
            } else {
                return std::pair<size_t, size_t>(chunk_id, get_num_chunks_per_row());
            }
        } else {
            if constexpr(accrow_) {
                return std::pair<size_t, size_t>(chunk_id, get_num_chunks_per_column());
            } else {
                return std::pair<size_t, size_t>(static_cast<size_t>(chunk_id) * static_cast<size_t>(get_num_chunks_per_column()), 1); // ditto.
            }
        }
    }

    template<bool accrow_, class ExtractFunction_>
    void extract_secondary_block(
        Index_ chunk_id, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        ExtractFunction_ extract)
    const {
        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ start_chunk_index = secondary_block_start / secondary_chunkdim;
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;
        Index_ secondary_block_end = secondary_block_start + secondary_block_length;
        Index_ end_chunk_index = secondary_block_end / secondary_chunkdim + (secondary_block_end % secondary_chunkdim > 0); // i.e., integer ceiling.

        auto oi = offset_and_increment<accrow_>(chunk_id);
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
                extract(chunk, from, len, secondary_start_pos);
            } else {
                extract(chunk, from, len);
            }

            secondary_start_pos += to; // yes, this is deliberate; '+ to' means that either we add 'secondary_chunkdim' or set it to 'secondary_block_end', the latter of which avoids overflow.
            offset += increment;
        }
    }

    template<bool accrow_, class ExtractFunction_>
    void extract_secondary_index(
        Index_ chunk_id, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer, 
        ExtractFunction_ extract)
    const {
        if (secondary_indices.empty()) {
            return;
        }

        auto secondary_chunkdim = get_secondary_chunkdim<accrow_>();
        Index_ start_chunk_index = secondary_indices.front() / secondary_chunkdim; // 'secondary_indices' is guaranteed to be non-empty at this point.
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;

        auto oi = offset_and_increment<accrow_>(chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        auto secondary_dim = get_secondary_dim<accrow_>();
        auto iIt = secondary_indices.begin();
        auto iEnd = secondary_indices.end();
        while (iIt != iEnd) {
            const auto& chunk = chunk_array[offset];

            Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // this convoluted method avoids overflow.
            chunk_indices_buffer.clear();
            while (iIt != iEnd && *iIt < secondary_end_pos) {
                chunk_indices_buffer.push_back(*iIt - secondary_start_pos);
                ++iIt;
            }

            if (!chunk_indices_buffer.empty()) {
                if constexpr(sparse_) {
                    extract(chunk, chunk_indices_buffer, secondary_start_pos);
                } else {
                    extract(chunk, chunk_indices_buffer);
                }
            }

            secondary_start_pos = secondary_end_pos;
            offset += increment;
        }
    }

    typedef typename Chunk_::Workspace ChunkWork;
    typedef typename Chunk_::value_type ChunkValue;
    typedef ChunkIndex_ ChunkIndex;
    typedef typename std::conditional<sparse_, typename SparseSlabFactory<Index_, ChunkValue, ChunkIndex>::Slab, typename DenseSlabFactory<ChunkValue>::Slab>::type Slab;
    typedef typename std::conditional<sparse_, SparseSingleWorkspace<Index_, ChunkValue, ChunkIndex>, DenseSingleWorkspace<ChunkValue> >::type SingleWorkspace;

public:
    // Extract a single element of the primary dimension, using a contiguous
    // block on the secondary dimension. 
    //
    // Unfortunately, we can't just re-use the fetch_block() functions with a
    // length of 1, because some chunks do not support partial extraction; this
    // requires special handling to extract the full chunk into 'tmp_work', and
    // then pull out what we need into 'final_slab'.
    //
    // Even the use_subset=true chunks that do support partial extraction
    // require a workspace involving the full chunk size, so we end up needing
    // the 'tmp_work' workspace anyway.
    template<bool accrow_>
    std::pair<const Slab*, Index_> fetch_single(
        Index_ i,
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        ChunkWork& chunk_workspace,
        SingleWorkspace& tmp_work,
        Slab& final_slab)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;

        if constexpr(sparse_) {
            auto& final_num = *final_slab.number;
            final_num = 0;
            bool needs_value = !final_slab.values.empty();
            bool needs_index = !final_slab.indices.empty();

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    std::fill_n(tmp_work.number.begin(), primary_chunkdim, 0);

                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, from, len, chunk_workspace, tmp_work.values, tmp_work.indices, tmp_work.number.data(), secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, tmp_work.values, tmp_work.indices, tmp_work.number.data(), secondary_start_pos);
                    }

                    auto count = tmp_work.number[chunk_offset];
                    if (needs_value) {
                        std::copy_n(tmp_work.values[chunk_offset], count, final_slab.values[0] + final_num);
                    }
                    if (needs_index) {
                        std::copy_n(tmp_work.indices[chunk_offset], count, final_slab.indices[0] + final_num);
                    }
                    final_num += count;
                }
            );

        } else {
            auto final_slab_ptr = final_slab.data;
            auto tmp_buffer_ptr = tmp_work.data();

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, from, len, chunk_workspace, tmp_buffer_ptr, len);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, tmp_buffer_ptr, len);
                    }

                    size_t tmp_offset = static_cast<size_t>(len) * static_cast<size_t>(chunk_offset);
                    std::copy_n(tmp_buffer_ptr + tmp_offset, len, final_slab_ptr);
                    final_slab_ptr += len;
                }
            );
        }

        return std::make_pair(&final_slab, static_cast<Index_>(0));
    }

    // Extract a single element of the primary dimension, using an indexed
    // subset on the secondary dimension.
    template<bool accrow_>
    std::pair<const Slab*, Index_> fetch_single(
        Index_ i,
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer,
        ChunkWork& chunk_workspace,
        SingleWorkspace& tmp_work,
        Slab& final_slab)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;

        if constexpr(sparse_) {
            auto& final_num = *final_slab.number;
            final_num = 0;
            bool needs_value = !final_slab.values.empty();
            bool needs_index = !final_slab.indices.empty();

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    std::fill_n(tmp_work.number.begin(), primary_chunkdim, 0);

                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, chunk_indices, chunk_workspace, tmp_work.values, tmp_work.indices, tmp_work.number.data(), secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, tmp_work.values, tmp_work.indices, tmp_work.number.data(), secondary_start_pos);
                    }

                    auto count = tmp_work.number[chunk_offset];
                    if (needs_value) {
                        std::copy_n(tmp_work.values[chunk_offset], count, final_slab.values[0] + final_num);
                    }
                    if (needs_index) {
                        std::copy_n(tmp_work.indices[chunk_offset], count, final_slab.indices[0] + final_num);
                    }
                    final_num += count;
                }
            );

        } else {
            auto final_slab_ptr = final_slab.data;
            auto tmp_buffer_ptr = tmp_work.data();

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    size_t nidx = chunk_indices.size();
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, chunk_indices, chunk_workspace, tmp_buffer_ptr, nidx);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, tmp_buffer_ptr, nidx);
                    }

                    size_t tmp_offset = nidx * static_cast<size_t>(chunk_offset);
                    std::copy_n(tmp_buffer_ptr + tmp_offset, nidx, final_slab_ptr);
                    final_slab_ptr += nidx;
                }
            );
        }

        return std::make_pair(&final_slab, static_cast<Index_>(0));
    }

private:
    // Extract a contiguous block of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_>
    void fetch_block(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim<accrow_>(), 0);

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_block_length;

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += len;
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_>
    void fetch_block(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer,
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim<accrow_>(), 0);

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_indices.size();

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += chunk_indices.size();
                }
            );
        }
    }

private:
    // Extract an indexed subset of the primary dimension, using a contiguous block on the secondary dimension.
    template<bool accrow_>
    void fetch_index(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim<accrow_>(), 0);

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_block_length;

            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += len;
                }
            );
        }
    }

    // Extract an indexed subset of the primary dimension, using an indexed subset on the secondary dimension.
    template<bool accrow_>
    void fetch_index(
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices,
        const std::vector<Index_>& secondary_indices,
        std::vector<Index_>& chunk_indices_buffer,
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim<accrow_>(), 0);

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_indices.size();

            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer, 
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += chunk_indices.size();
                }
            );
        }
    }

public:
    // Obtain the slab containing the 'i'-th element of the primary dimension.
    template<bool accrow_, class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_myopic(
        Index_ i, 
        Index_ block_start,
        Index_ block_length,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;
        auto& out = cache.find(
            chunk_id,
            /* create = */ [&]() -> Slab {
                return factory.create();
            },
            /* populate = */ [&](Index_ id, Slab& slab) -> void {
                fetch_block<accrow_>(id, 0, get_primary_chunkdim<accrow_>(id), block_start, block_length, slab, chunk_workspace);
            }
        );
        return std::make_pair(&out, chunk_offset);
    }

    template<bool accrow_, class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_myopic(
        Index_ i, 
        const std::vector<Index_>& indices,
        std::vector<Index_>& tmp_indices,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;
        auto& out = cache.find(
            chunk_id,
            /* create = */ [&]() -> Slab {
                return factory.create();
            },
            /* populate = */ [&](Index_ id, Slab& slab) -> void {
                fetch_block<accrow_>(id, 0, get_primary_chunkdim<accrow_>(id), indices, tmp_indices, slab, chunk_workspace);
            }
        );
        return std::make_pair(&out, chunk_offset);
    }

private:
    template<bool accrow_, class Cache_, class Factory_, typename FetchBlock_, typename FetchIndex_>
    std::pair<const Slab*, Index_> fetch_oracular_core(
        Cache_& cache,
        Factory_& factory,
        FetchBlock_ fetch_block,
        FetchIndex_ fetch_index)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        if constexpr(Chunk_::use_subset) {
            return cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::pair<Index_, Index_>(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab {
                    return factory.create();
                },
                /* populate =*/ [&](std::vector<std::tuple<Index_, Slab*, OracularSubsettedSlabCacheSelectionDetails<Index_>*> >& in_need) -> void {
                    for (const auto& p : in_need) {
                        auto id = std::get<0>(p);
                        auto ptr = std::get<1>(p);
                        auto sub = std::get<2>(p);
                        switch (sub->selection) {
                            case OracularSubsettedSlabCacheSelectionType::FULL:
                                fetch_block(id, 0, get_primary_chunkdim<accrow_>(id), *ptr);
                                break;
                            case OracularSubsettedSlabCacheSelectionType::BLOCK:
                                fetch_block(id, sub->block_start, sub->block_length, *ptr);
                                break;
                            case OracularSubsettedSlabCacheSelectionType::INDEX:
                                fetch_index(id, sub->indices, *ptr);
                                break;
                        }
                    }
                }
            );

        } else {
            return cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::pair<Index_, Index_>(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab {
                    return factory.create();
                },
                /* populate =*/ [&](std::vector<std::pair<Index_, Slab*> >& to_populate) -> void {
                    for (auto& p : to_populate) {
                        fetch_block(p.first, 0, get_primary_chunkdim<accrow_>(p.first), *(p.second));
                    }
                }
            );
        }
    }

public:
    template<bool accrow_, class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_oracular(
        Index_ block_start,
        Index_ block_length,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        return fetch_oracular_core<accrow_>(
            cache,
            factory,
            [&](Index_ pid, Index_ pstart, Index_ plen, Slab& slab) {
                fetch_block<accrow_>(pid, pstart, plen, block_start, block_length, slab, chunk_workspace);
            },
            [&](Index_ pid, const std::vector<Index_>& pindices, Slab& slab) {
                fetch_index<accrow_>(pid, pindices, block_start, block_length, slab, chunk_workspace);
            }
        );
    }

    template<bool accrow_, class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_oracular(
        const std::vector<Index_>& indices,
        std::vector<Index_>& chunk_indices_buffer,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        return fetch_oracular_core<accrow_>(
            cache,
            factory,
            [&](Index_ pid, Index_ pstart, Index_ plen, Slab& slab) {
                fetch_block<accrow_>(pid, pstart, plen, indices, chunk_indices_buffer, slab, chunk_workspace);
            },
            [&](Index_ pid, const std::vector<Index_>& pindices, Slab& slab) {
                fetch_index<accrow_>(pid, pindices, indices, chunk_indices_buffer, slab, chunk_workspace);
            }
        );
    }
};

}

}

#endif
