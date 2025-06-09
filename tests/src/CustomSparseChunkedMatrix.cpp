#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"

#include "tatami_chunked/CustomSparseChunkedMatrix.hpp"

typedef double ChunkValue_;
typedef int Index_;

struct MockSparseChunk {
    // Standard compressed sparse members:
    std::vector<ChunkValue_> values;
    std::vector<Index_> indices;
    std::vector<std::size_t> indptrs;
};

struct MockSparseChunkData { 
    tatami_chunked::ChunkDimensionStats<Index_> row_stats, col_stats;
    std::vector<MockSparseChunk> chunks;
};

class MockSparseChunkWorkspace final : public tatami_chunked::CustomSparseChunkedMatrixWorkspace<ChunkValue_, Index_> {
public:
    MockSparseChunkWorkspace(const MockSparseChunkData& data) : my_data(data) {}

private:
    const MockSparseChunkData& my_data;

    // Allocation to allow for O(1) mapping of requested indices to sparse indices.
    // This mimics what is done in the indexed sparse extractors in tatami proper.
    std::vector<unsigned char> remap;

private:
    auto get_target_chunkdim(bool row) const {
        if (row) {
            return my_data.row_stats.chunk_length;
        } else {
            return my_data.col_stats.chunk_length;
        }
    }

    auto get_non_target_chunkdim(bool row) const {
        if (row) {
            return my_data.col_stats.chunk_length;
        } else {
            return my_data.row_stats.chunk_length;
        }
    }

    static void refine_start_and_end(std::size_t& start, std::size_t& end, Index_ desired_start, Index_ desired_end, Index_ max_end, const std::vector<Index_>& indices) {
        if (desired_start) {
            auto it = indices.begin();
            // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
            start = std::lower_bound(it + start, it + end, desired_start, [](Index_ a, Index_ b) -> bool { return a < b; }) - it;
        }

        if (desired_end != max_end) {
            if (desired_end == desired_start + 1) {
                if (start != end && static_cast<Index_>(indices[start]) == desired_start) {
                    end = start + 1;
                } else {
                    end = start;
                }
            } else {
                auto it = indices.begin();
                end = std::lower_bound(it + start, it + end, desired_end, [](Index_ a, Index_ b) -> bool { return a < b; }) - it;
            }
        }
    }

    // Building a present/absent mapping for the requested indices to allow O(1) look-up from the sparse matrix indices.
    void configure_remap(const std::vector<Index_>& indices, std::size_t full) {
        remap.resize(full);
        for (auto i : indices) {
            remap[i] = 1;
        }
    }

    // Resetting just the affected indices so we can avoid a fill operation over the entire array.
    void reset_remap(const std::vector<Index_>& indices) {
        for (auto i : indices) {
            remap[i] = 0;
        }
    }

    template<bool is_block_>
    void fill_target(
        Index_ p, 
        const MockSparseChunk& chunk, 
        Index_ non_target_start, 
        Index_ non_target_end, 
        Index_ non_target_chunkdim,
        const std::vector<ChunkValue_*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        size_t start = chunk.indptrs[p], end = chunk.indptrs[p + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, non_target_start, non_target_end, non_target_chunkdim, chunk.indices);

        auto& current_number = output_number[p];
        const bool needs_value = !output_values.empty();
        auto vptr = needs_value ? output_values[p] + current_number : NULL;
        const bool needs_index = !output_indices.empty();
        auto iptr = needs_index ? output_indices[p] + current_number : NULL;

        if constexpr(is_block_) {
            if (needs_value) {
                std::copy(chunk.values.begin() + start, chunk.values.begin() + end, vptr);
            }
            if (needs_index) {
                for (size_t i = start; i < end; ++i, ++iptr) {
                    *iptr = chunk.indices[i] + shift;
                }
            }
            current_number += end - start;

        } else {
            // Assumes that chunk.remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = chunk.indices[i];
                if (remap[target]) {
                    if (needs_value) {
                        *vptr = chunk.values[i];
                        ++vptr;
                    }
                    if (needs_index) {
                        *iptr = target + shift;
                        ++iptr;
                    }
                    ++current_number;
                }
            }
        }
    }

    template<bool is_block_>
    void fill_secondary(
        Index_ s,
        const MockSparseChunk& chunk, 
        Index_ target_start, 
        Index_ target_end, 
        Index_ target_chunkdim,
        const std::vector<ChunkValue_*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        auto start = chunk.indptrs[s], end = chunk.indptrs[s + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, target_start, target_end, target_chunkdim, chunk.indices);

        bool needs_value = !output_values.empty();
        bool needs_index = !output_indices.empty();

        if constexpr(is_block_) {
            for (size_t i = start; i < end; ++i) {
                auto p = chunk.indices[i];
                auto& num = output_number[p];
                if (needs_value) {
                    output_values[p][num] = chunk.values[i];
                }
                if (needs_index) {
                    output_indices[p][num] = s + shift;
                }
                ++num;
            }

        } else {
            // Assumes that remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = chunk.indices[i];
                if (remap[target]) {
                    auto& num = output_number[target];
                    if (needs_value) {
                        output_values[target][num] = chunk.values[i];
                    }
                    if (needs_index) {
                        output_indices[target][num] = s + shift;
                    }
                    ++num;
                }
            }
        }
    }

public:
    void extract(
        Index_ chunk_row_id,
        Index_ chunk_column_id,
        bool row,
        Index_ target_start, 
        Index_ target_length, 
        Index_ non_target_start, 
        Index_ non_target_length, 
        const std::vector<ChunkValue_*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        const auto& current_chunk = my_data.chunks[chunk_row_id * my_data.col_stats.num_chunks + chunk_column_id];
        Index_ target_end = target_start + target_length;
        Index_ non_target_end = non_target_start + non_target_length;

        if (row) {
            Index_ non_target_chunkdim = get_non_target_chunkdim(row);
            for (Index_ p = target_start; p < target_end; ++p) {
                fill_target<true>(p, current_chunk, non_target_start, non_target_end, non_target_chunkdim, output_values, output_indices, output_number, shift);
            }

        } else {
            Index_ target_chunkdim = get_target_chunkdim(row);
            for (Index_ s = non_target_start; s < non_target_end; ++s) {
                fill_secondary<true>(s, current_chunk, target_start, target_end, target_chunkdim, output_values, output_indices, output_number, shift);
            }
        }
    }

    void extract(
        Index_ chunk_row_id,
        Index_ chunk_column_id,
        bool row,
        Index_ target_start, 
        Index_ target_length, 
        const std::vector<Index_>& non_target_indices, 
        const std::vector<ChunkValue_*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        const auto& current_chunk = my_data.chunks[chunk_row_id * my_data.col_stats.num_chunks + chunk_column_id];
        Index_ target_end = target_start + target_length;

        if (row) {
            // non_target_indices is guaranteed to be non-empty, see contracts below.
            auto non_target_start = non_target_indices.front();
            auto non_target_end = non_target_indices.back() + 1; // need 1 past end.
            auto non_target_chunkdim = get_non_target_chunkdim(row);

            configure_remap(non_target_indices, non_target_chunkdim);
            for (Index_ p = target_start; p < target_end; ++p) {
                fill_target<false>(p, current_chunk, non_target_start, non_target_end, non_target_chunkdim, output_values, output_indices, output_number, shift);
            }
            reset_remap(non_target_indices);

        } else {
            Index_ target_chunkdim = get_target_chunkdim(row);
            for (auto s : non_target_indices) {
                fill_secondary<true>(s, current_chunk, target_start, target_end, target_chunkdim, output_values, output_indices, output_number, shift);
            }
        }
    }

    void extract(
        Index_ chunk_row_id,
        Index_ chunk_column_id,
        bool row,
        const std::vector<Index_>& target_indices, 
        Index_ non_target_start, 
        Index_ non_target_length, 
        const std::vector<ChunkValue_*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        const auto& current_chunk = my_data.chunks[chunk_row_id * my_data.col_stats.num_chunks + chunk_column_id];
        Index_ non_target_end = non_target_start + non_target_length;

        if (row) {
            Index_ non_target_chunkdim = get_non_target_chunkdim(row);
            for (auto p : target_indices) {
                fill_target<true>(p, current_chunk, non_target_start, non_target_end, non_target_chunkdim, output_values, output_indices, output_number, shift);
            }

        } else {
            // target_indices is guaranteed to be non-empty, see contracts below.
            auto target_start = target_indices.front();
            auto target_end = target_indices.back() + 1; // need 1 past end.
            auto target_chunkdim = get_target_chunkdim(row);

            configure_remap(target_indices, target_chunkdim);
            for (Index_ s = non_target_start; s < non_target_end; ++s) {
                fill_secondary<false>(s, current_chunk, target_start, target_end, target_chunkdim, output_values, output_indices, output_number, shift);
            }
            reset_remap(target_indices);
        }
    }

    void extract(
        Index_ chunk_row_id,
        Index_ chunk_column_id,
        bool row,
        const std::vector<Index_>& target_indices, 
        const std::vector<Index_>& non_target_indices, 
        const std::vector<ChunkValue_*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    {
        const auto& current_chunk = my_data.chunks[chunk_row_id * my_data.col_stats.num_chunks + chunk_column_id];

        if (row) {
            // non_target_indices is guaranteed to be non-empty, see contracts below.
            auto non_target_start = non_target_indices.front();
            auto non_target_end = non_target_indices.back() + 1; // need 1 past end.
            auto non_target_chunkdim = get_non_target_chunkdim(row);

            configure_remap(non_target_indices, non_target_chunkdim);
            for (auto p : target_indices) {
                fill_target<false>(p, current_chunk, non_target_start, non_target_end, non_target_chunkdim, output_values, output_indices, output_number, shift);
            }
            reset_remap(non_target_indices);

        } else {
            // target_indices is guaranteed to be non-empty, see contracts below.
            auto target_start = target_indices.front();
            auto target_end = target_indices.back() + 1; // need 1 past end.
            auto target_chunkdim = get_target_chunkdim(row);

            configure_remap(target_indices, target_chunkdim);
            for (auto s : non_target_indices) {
                fill_secondary<false>(s, current_chunk, target_start, target_end, target_chunkdim, output_values, output_indices, output_number, shift);
            }
            reset_remap(target_indices);
        }
    }
};

class MockSparseChunkManager final : public tatami_chunked::CustomSparseChunkedMatrixManager<ChunkValue_, Index_> {
public:
    MockSparseChunkManager(MockSparseChunkData data) : my_data(std::move(data)) {}

    std::unique_ptr<tatami_chunked::CustomSparseChunkedMatrixWorkspace<ChunkValue_, Index_> > new_workspace() const {
        return std::make_unique<MockSparseChunkWorkspace>(my_data);
    }

    bool prefer_rows() const {
        return true;
    }

    const tatami_chunked::ChunkDimensionStats<Index_>& row_stats() const {
        return my_data.row_stats;
    }

    const tatami_chunked::ChunkDimensionStats<Index_>& column_stats() const {
        return my_data.col_stats;
    }

private:
    MockSparseChunkData my_data; 
};

class CustomSparseChunkedMatrixCore {
public:
    typedef std::tuple<
        std::pair<int, int>, // matrix dimensions
        std::pair<int, int>, // chunk dimensions
        double // cache fraction
    > SimulationParameters;

protected:
    inline static std::unique_ptr<tatami::Matrix<double, int> > ref, simple_mat, subset_mat;
    inline static SimulationParameters last_params;

    static void assemble(const SimulationParameters& params) {
        if (ref && params == last_params) {
            return;
        }
        last_params = params;

        auto matdim = std::get<0>(params);
        auto chunkdim = std::get<1>(params);
        double cache_fraction = std::get<2>(params);

        auto full = tatami_test::simulate_compressed_sparse<double, int>(matdim.second, matdim.first, [&]{
            tatami_test::SimulateCompressedSparseOptions opt;
            opt.density = 0.1;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = matdim.first * matdim.second + chunkdim.first * chunkdim.second + 100 * cache_fraction;
            return opt;
        }());

        ref.reset(new tatami::CompressedSparseColumnMatrix<double, int>(
            matdim.first,
            matdim.second,
            std::move(full.data),
            std::move(full.index),
            std::move(full.indptr)
        ));

        MockSparseChunkData data;
        data.row_stats = tatami_chunked::ChunkDimensionStats<Index_>(matdim.first, chunkdim.first);
        data.col_stats = tatami_chunked::ChunkDimensionStats<Index_>(matdim.second, chunkdim.second);
        data.chunks.resize(data.row_stats.num_chunks * data.col_stats.num_chunks);

        for (int r = 0; r < data.row_stats.num_chunks; ++r) {
            for (int c = 0; c < data.col_stats.num_chunks; ++c) {
                auto cstart = c * chunkdim.second;
                auto cend = std::min(cstart + chunkdim.second, matdim.second);
                auto clen = cend - cstart;

                auto rstart = r * chunkdim.first;
                auto rend = std::min(rstart + chunkdim.first, matdim.first);
                auto rlen = rend - rstart;

                MockSparseChunk chunk;
                chunk.indptrs.resize(1);
                auto ext = ref->sparse_row(cstart, clen);
                std::vector<double> vbuffer(clen);
                std::vector<int> ibuffer(clen);

                for (int r2 = 0; r2 < rlen; ++r2) {
                    auto range = ext->fetch(r2 + rstart, vbuffer.data(), ibuffer.data());
                    chunk.values.insert(chunk.values.end(), range.value, range.value + range.number);
                    for (int i = 0; i < range.number; ++i) {
                        chunk.indices.push_back(range.index[i] - cstart);
                    }
                    chunk.indptrs.push_back(chunk.indptrs.back() + range.number);
                }

                auto offset = r * data.col_stats.num_chunks + c;
                data.chunks[offset] = std::move(chunk);
            }
        }

        tatami_chunked::CustomSparseChunkedMatrixOptions opt;
        std::size_t cache_size = static_cast<double>(matdim.first) * static_cast<double>(matdim.second) * cache_fraction * static_cast<double>(sizeof(double) + sizeof(int));
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = (cache_size > 0);

        auto manager = std::make_shared<MockSparseChunkManager>(std::move(data));
        simple_mat.reset(new tatami_chunked::CustomSparseChunkedMatrix<double, int, double>(manager, opt));

        opt.cache_subset = true;
        subset_mat.reset(new tatami_chunked::CustomSparseChunkedMatrix<double, int, double>(manager, opt));
    }
};

/*******************************************************/

class CustomSparseChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<typename CustomSparseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions> >, 
    public CustomSparseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomSparseChunkedMatrixFullTest, Basic) {
    auto opt = tatami_test::convert_test_access_options(std::get<1>(GetParam()));
    tatami_test::test_full_access(*simple_mat, *ref, opt);
    tatami_test::test_full_access(*subset_mat, *ref, opt);
}

INSTANTIATE_TEST_SUITE_P(
    CustomSparseChunkedMatrix,
    CustomSparseChunkedMatrixFullTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values( // matrix dimensions
                std::make_pair(200, 50),
                std::make_pair(100, 300)
            ),

            ::testing::Values( // chunk dimensions
                std::make_pair(1, 20),
                std::make_pair(20, 1),
                std::make_pair(11, 13) // odd numbers
            ),

            ::testing::Values(0, 0.01, 0.1) // cache fraction
        ),

        tatami_test::standard_test_access_options_combinations()
    )
);

/*******************************************************/

class CustomSparseChunkedMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<typename CustomSparseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CustomSparseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomSparseChunkedMatrixBlockTest, Basic) {
    auto tparam = GetParam();
    auto opt = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto block = std::get<2>(tparam);
    tatami_test::test_block_access(*simple_mat, *ref, block.first, block.second, opt);
    tatami_test::test_block_access(*subset_mat, *ref, block.first, block.second, opt);
}

INSTANTIATE_TEST_SUITE_P(
    CustomSparseChunkedMatrix,
    CustomSparseChunkedMatrixBlockTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values( // matrix dimensions
                std::make_pair(201, 67),
                std::make_pair(123, 372)
            ),

            ::testing::Values( // chunk dimensions
                std::make_pair(1, 20),
                std::make_pair(20, 1),
                std::make_pair(10, 10)
            ),

            ::testing::Values(0, 0.01, 0.1) // cache fraction
        ),

        tatami_test::standard_test_access_options_combinations(),

        ::testing::Values( // block boundaries
            std::make_pair(0.0, 0.35),
            std::make_pair(0.15, 0.71),
            std::make_pair(0.38, 0.62)
        )
    )
);

/*******************************************************/

class CustomSparseChunkedMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<typename CustomSparseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CustomSparseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomSparseChunkedMatrixIndexTest, Basic) {
    auto tparam = GetParam();
    auto opt = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto index = std::get<2>(tparam);
    tatami_test::test_indexed_access(*simple_mat, *ref, index.first, index.second, opt);
    tatami_test::test_indexed_access(*subset_mat, *ref, index.first, index.second, opt);
}

INSTANTIATE_TEST_SUITE_P(
    CustomSparseChunkedMatrix,
    CustomSparseChunkedMatrixIndexTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values( // matrix dimensions
                std::make_pair(198, 67),
                std::make_pair(187, 300)
            ),

            ::testing::Values( // chunk dimensions
                std::make_pair(1, 20),
                std::make_pair(20, 1),
                std::make_pair(7, 13)
            ),

            ::testing::Values(0, 0.01, 0.1) // cache fraction
        ),

        tatami_test::standard_test_access_options_combinations(),

        ::testing::Values( // index information.
            std::make_pair(0.0, 0.15),
            std::make_pair(0.2, 0.24),
            std::make_pair(0.7, 0.3)
        )
    )
);
