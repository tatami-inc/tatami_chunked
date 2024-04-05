#ifndef TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP
#define TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP

#include "tatami/tatami.hpp"
#include "typical_slab_cache.hpp"

#include <vector>

namespace tatami_chunked {

namespace CustomChunkedMatrix_internal {

template<
    typename Index_, 
    bool sparse_,   
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

    Index_ get_chunk_nrow() const {
        return chunk_nrow;
    }

    Index_ get_chunk_ncol() const {
        return chunk_ncol;
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
    typedef typename std::conditional<sparse_, SparseSlab, DenseSlab>::type Slab;

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

    static void flush_slab(Slab& slab) {
        if constexpr(sparse_) {
            for (auto& x : slab.indices) {
                x.clear();
            }
            for (auto& x : slab.values) {
                x.clear();
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
        Index_ end_chunk_index = integer_ceil(secondary_block_end, secondary_chunkdim);

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

            Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // avoid overflow.
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

private:
    // Extract a single element of the primary dimension, using a contiguous
    // block on the secondary dimension. 
    //
    // Unfortunately, we can't just re-use the fetch_block() functions with a
    // length of 1, because some chunks do not support partial extraction; this
    // requires special handling to extract the full chunk into 'tmp_slab', and
    // then pull out what we need into 'final_slab'.
    //
    // Even the use_subset=true chunks that do support partial extraction
    // require a workspace involving the full chunk size, so we end up needing
    // the 'tmp_slab' anyway.
    template<bool accrow_>
    void fetch_single(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& tmp_slab, 
        Slab& final_slab,
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            flush_slab(final_slab);
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    flush_slab(tmp_slab);
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, from, len, chunk_workspace, tmp_slab.values, tmp_slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, tmp_slab.values, tmp_slab.indices, secondary_start_pos);
                    }
                    final_slab.values.insert(final_slab.values.end(), tmp_slab.values[chunk_offset].begin(), tmp_slab.values[chunk_offset].end());
                    final_slab.indices.insert(final_slab.indices.end(), tmp_slab.indices[chunk_offset].begin(), tmp_slab.indices[chunk_offset].end());
                }
            );
        } else {
            auto final_slab_ptr = final_slab.data();
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    auto tmp_slab_ptr = tmp_slab.data();
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, from, len, chunk_workspace, tmp_slab_ptr, len);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, tmp_slab_ptr, len);
                    }
                    auto tmp_start = tmp_slab_ptr + len * chunk_offset;
                    std::copy(tmp_start, tmp_start + len, final_slab_ptr);
                    final_slab_ptr += len;
                }
            );
        }
    }

    // Extract a single element of the primary dimension, using an indexed
    // subset on the secondary dimension.
    template<bool accrow_>
    void fetch_single(
        Index_ chunk_id, 
        Index_ chunk_offset, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer,
        Slab& tmp_slab, 
        Slab& final_slab,
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            flush_slab(final_slab);
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    flush_slab(tmp_slab);
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, chunk_indices, chunk_workspace, tmp_slab.values, tmp_slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, tmp_slab.values, tmp_slab.indices, secondary_start_pos);
                    }
                    final_slab.values.insert(final_slab.values.end(), tmp_slab.values[chunk_offset].begin(), tmp_slab.values[chunk_offset].end());
                    final_slab.indices.insert(final_slab.indices.end(), tmp_slab.indices[chunk_offset].begin(), tmp_slab.indices[chunk_offset].end());
                }
            );
        } else {
            auto final_slab_ptr = final_slab.data();
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    auto tmp_slab_ptr = tmp_slab.data();
                    size_t nidx = chunk_indices.size();
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, 1, chunk_indices, chunk_workspace, tmp_slab_ptr, nidx);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, tmp_slab_ptr, nidx);
                    }
                    auto tmp_start = tmp_slab_ptr + nidx * chunk_offset;
                    std::copy(tmp_start, tmp_start + nidx, final_slab_ptr);
                    final_slab_ptr += nidx;
                }
            );
        }
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
            flush_slab(slab);
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    }
                }
            );
        } else {
            auto slab_ptr = slab.data();
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab_ptr, secondary_block_length);
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
            flush_slab(slab);
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    }
                }
            );
        } else {
            auto slab_ptr = slab.data();
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
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
            flush_slab(slab);
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    }
                }
            );
        } else {
            auto slab_ptr = slab.data();
            extract_secondary_block<accrow_>(
                chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, from, len, chunk_workspace, slab_ptr, secondary_block_length);
                    } else {
                        chunk.template extract<accrow_>(from, len, chunk_workspace, slab_ptr, secondary_block_length);
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
            flush_slab(slab);
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab.values, slab.indices, secondary_start_pos);
                    }
                }
            );
        } else {
            auto slab_ptr = slab.data();
            extract_secondary_index<accrow_>(
                chunk_id, secondary_indices, chunk_indices_buffer, 
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.template extract<accrow_>(primary_indices, chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                    } else {
                        chunk.template extract<accrow_>(chunk_indices, chunk_workspace, slab_ptr, secondary_indices.size());
                    }
                    slab_ptr += chunk_indices.size();
                }
            );
        }
    }

public:
    // Obtain the slab containing the 'i'-th element of the primary dimension.
    template<bool accrow_, bool oracle_, typename FetchSingle_, typename FetchBlock_, typename FetchIndex_>
    std::pair<const Slab*, Index_> fetch_core(
        Index_ i, 
        TypicalSlabCacheWorkspace<oracle_, Chunk_::use_subset, Index_, Slab>& cache_workspace,
        Slab& tmp_solo,
        Slab& final_solo,
        Index_ secondary_length,
        FetchSingle_ fetch_single,
        FetchBlock_ fetch_block,
        FetchIndex_ fetch_index)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim<accrow_>();
        size_t alloc = sparse_ ? primary_chunkdim : primary_chunkdim * static_cast<size_t>(secondary_length); // use size_t to avoid integer overflow.

        if (cache_workspace.num_slabs_in_cache == 0) {
            if constexpr(oracle_) {
                i = cache_workspace.cache.next();
            }
            fetch_single(i / primary_chunkdim, i % primary_chunkdim, tmp_solo, final_solo);
            return std::make_pair(&final_solo, static_cast<Index_>(0));

        } else if constexpr(!oracle_) {
            auto chunk_id = i / primary_chunkdim;
            auto chunk_offset = i % primary_chunkdim;
            auto& cache = cache_workspace.cache.find(
                chunk_id,
                /* create = */ [&]() -> Slab {
                    return Slab(alloc);
                },
                /* populate = */ [&](Index_ id, Slab& slab) -> void {
                    fetch_block(id, 0, get_primary_chunkdim<accrow_>(id), slab);
                }
            );
            return std::make_pair(&cache, chunk_offset);

        } else if constexpr(Chunk_::use_subset) {
            auto out = cache_workspace.cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::make_pair(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab {
                    return Slab(alloc);
                },
                /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                    for (const auto& p : in_need) {
                        auto& ptr = data[p.second];
                        switch (ptr->subset.selection) {
                            case SubsetSelection::FULL:
                                fetch_block(p.first, 0, get_primary_chunkdim<accrow_>(p.first), ptr->contents);
                                break;
                            case SubsetSelection::BLOCK:
                                fetch_block(p.first, ptr->subset.block_start, ptr->subset.block_length, ptr->contents);
                                break;
                            case SubsetSelection::INDEX:
                                fetch_index(p.first, ptr->subset.indices, ptr->contents);
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
                /* create = */ [&]() -> Slab {
                    return Slab(alloc);
                },
                /* populate =*/ [&](const std::vector<std::pair<Index_, Index_> >& in_need, auto& data) -> void {
                    for (const auto& p : in_need) {
                        fetch_block(p.first, 0, get_primary_chunkdim<accrow_>(p.first), *(data[p.second]));
                    }
                }
            );
        }
    }

    template<bool accrow_, bool oracle_>
    std::pair<const Slab*, Index_> fetch(
        Index_ i, 
        Index_ block_start,
        Index_ block_length,
        TypicalSlabCacheWorkspace<oracle_, Chunk_::use_subset, Index_, Slab>& cache_workspace,
        ChunkWork& chunk_workspace,
        Slab& tmp_solo,
        Slab& final_solo)
    const {
        return fetch_core<accrow_, oracle_>(
            i, 
            cache_workspace, 
            tmp_solo, 
            final_solo, 
            block_length,
            [&](Index_ pid, Index_ pstart, Slab& tmp_solo, Slab& final_solo) {
                fetch_single<accrow_>(pid, pstart, block_start, block_length, tmp_solo, final_solo, chunk_workspace);
            },
            [&](Index_ pid, Index_ pstart, Index_ plen, Slab& slab) {
                fetch_block<accrow_>(pid, pstart, plen, block_start, block_length, slab, chunk_workspace);
            },
            [&](Index_ pid, const std::vector<Index_>& pindices, Slab& slab) {
                fetch_index<accrow_>(pid, pindices, block_start, block_length, slab, chunk_workspace);
            }
        );
    }

    template<bool accrow_, bool oracle_>
    std::pair<const Slab*, Index_> fetch(
        Index_ i, 
        const std::vector<Index_>& indices,
        std::vector<Index_>& chunk_indices_buffer,
        TypicalSlabCacheWorkspace<oracle_, Chunk_::use_subset, Index_, Slab>& cache_workspace,
        ChunkWork& chunk_workspace,
        Slab& tmp_solo,
        Slab& final_solo)
    const {
        return fetch_core<accrow_, oracle_>(
            i, 
            cache_workspace, 
            tmp_solo, 
            final_solo, 
            indices.size(),
            [&](Index_ pid, Index_ pstart, Slab& tmp_solo, Slab& final_solo) {
                fetch_single<accrow_>(pid, pstart, indices, chunk_indices_buffer, tmp_solo, final_solo, chunk_workspace);
            },
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
