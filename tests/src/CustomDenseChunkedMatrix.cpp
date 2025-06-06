#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"

#include "tatami_chunked/CustomDenseChunkedMatrix.hpp"

typedef double ChunkValue_;
typedef int Index_;

struct MockDenseChunkData { 
    tatami_chunked::ChunkDimensionStats<Index_> row_stats, col_stats;
    std::vector<std::vector<ChunkValue_> > chunks;
};

class MockDenseChunkWorkspace final : public tatami_chunked::CustomDenseChunkedMatrixWorkspace<ChunkValue_, Index_> {
public:
    MockDenseChunkWorkspace(const MockDenseChunkData& data) : my_data(data) {}

    void extract(
        Index_ chunk_row,
        Index_ chunk_column,
        bool row,
        Index_ target_start,
        Index_ target_length,
        Index_ non_target_start,
        Index_ non_target_length,
        ChunkValue_* output,
        Index_ stride
    ) {
        const auto& curchunk = my_data.chunks[chunk_row * my_data.col_stats.num_chunks + chunk_column];
        if (row) {
            for (Index_ tidx = target_start, tend = target_start + target_length; tidx < tend; ++tidx) {
                std::size_t in_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(non_target_start);
                std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride);
                std::copy_n(curchunk.data() + in_offset, non_target_length, output + out_offset);
            }
        } else {
            // Could be a bit more efficient by using tatami::transpose(), but whatever.
            for (Index_ tidx = target_start, tend = target_start + target_length; tidx < tend; ++tidx) {
                for (Index_ nidx = 0; nidx < non_target_length; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(nidx + non_target_start) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(tidx);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        }
    }

    void extract(
        Index_ chunk_row,
        Index_ chunk_column,
        bool row,
        Index_ target_start,
        Index_ target_length,
        const std::vector<Index_>& non_target_indices,
        ChunkValue_* output,
        Index_ stride
    ) {
        const auto& curchunk = my_data.chunks[chunk_row * my_data.col_stats.num_chunks + chunk_column];
        auto ntsize = non_target_indices.size();
        if (row) {
            for (Index_ tidx = target_start, tend = target_start + target_length; tidx < tend; ++tidx) {
                for (decltype(ntsize) nidx = 0; nidx < ntsize; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(non_target_indices[nidx]);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        } else {
            for (Index_ tidx = target_start, tend = target_start + target_length; tidx < tend; ++tidx) {
                for (decltype(ntsize) nidx = 0; nidx < ntsize; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(non_target_indices[nidx]) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(tidx);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        }
    }

    void extract(
        Index_ chunk_row,
        Index_ chunk_column,
        bool row,
        const std::vector<Index_>& target_indices,
        Index_ non_target_start,
        Index_ non_target_length,
        ChunkValue_* output,
        Index_ stride
    ) {
        const auto& curchunk = my_data.chunks[chunk_row * my_data.col_stats.num_chunks + chunk_column];
        if (row) {
            for (auto tidx : target_indices) {
                for (Index_ nidx = 0; nidx < non_target_length; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(nidx + non_target_start);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        } else {
            for (auto tidx : target_indices) {
                for (Index_ nidx = 0; nidx < non_target_length; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(nidx + non_target_start) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(tidx);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        }
    }

    void extract(
        Index_ chunk_row,
        Index_ chunk_column,
        bool row,
        const std::vector<Index_>& target_indices,
        const std::vector<Index_>& non_target_indices,
        ChunkValue_* output,
        Index_ stride
    ) {
        const auto& curchunk = my_data.chunks[chunk_row * my_data.col_stats.num_chunks + chunk_column];
        auto ntsize = non_target_indices.size();
        if (row) {
            for (auto tidx : target_indices) {
                for (decltype(ntsize) nidx = 0; nidx < ntsize; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(non_target_indices[nidx]);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        } else {
            for (auto tidx : target_indices) {
                for (decltype(ntsize) nidx = 0; nidx < ntsize; ++nidx) {
                    std::size_t in_offset = static_cast<std::size_t>(non_target_indices[nidx]) * static_cast<std::size_t>(my_data.col_stats.chunk_length) + static_cast<std::size_t>(tidx);
                    std::size_t out_offset = static_cast<std::size_t>(tidx) * static_cast<std::size_t>(stride) + static_cast<std::size_t>(nidx);
                    output[out_offset] = curchunk[in_offset];
                }
            }
        }
    }

private:
    const MockDenseChunkData& my_data;
};

class MockDenseChunkManager final : public tatami_chunked::CustomDenseChunkedMatrixManager<ChunkValue_, Index_> {
public:
    MockDenseChunkManager(MockDenseChunkData data) : my_data(std::move(data)) {}

    std::unique_ptr<tatami_chunked::CustomDenseChunkedMatrixWorkspace<ChunkValue_, Index_> > new_workspace() const {
        return std::make_unique<MockDenseChunkWorkspace>(my_data);
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
    MockDenseChunkData my_data; 
};

struct CustomDenseChunkedMatrixCore {
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

        auto full = tatami_test::simulate_vector<double>(matdim.first * matdim.second, [&]{
            tatami_test::SimulateVectorOptions opt;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = matdim.first * matdim.second + chunkdim.first * chunkdim.second + 100 * cache_fraction;
            return opt;
        }());
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        MockDenseChunkData data;
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

                auto& contents = data.chunks[r * data.col_stats.num_chunks + c];
                contents.resize(chunkdim.first * chunkdim.second);
                auto ccptr = contents.data();
                auto ext = ref->dense_row(cstart, clen);
                for (int r2 = 0; r2 < rlen; ++r2) {
                    auto ptr = ext->fetch(r2 + rstart, ccptr);
                    tatami::copy_n(ptr, clen, ccptr);
                    ccptr += chunkdim.second;
                }
            }
        }

        tatami_chunked::CustomDenseChunkedMatrixOptions opt;
        std::size_t cache_size = static_cast<double>(matdim.first) * static_cast<double>(matdim.second) * cache_fraction * static_cast<double>(sizeof(double));
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = (cache_size > 0);

        auto manager = std::make_shared<MockDenseChunkManager>(std::move(data));
        simple_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, double>(manager, opt));

        opt.cache_subset = true;
        subset_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, double>(manager, opt));
    }
};

/*******************************************************/

class CustomDenseChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<typename CustomDenseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions> >, 
    public CustomDenseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomDenseChunkedMatrixFullTest, Basic) {
    auto opts = tatami_test::convert_test_access_options(std::get<1>(GetParam()));
    tatami_test::test_full_access(*simple_mat, *ref, opts);
    tatami_test::test_full_access(*subset_mat, *ref, opts);
}

INSTANTIATE_TEST_SUITE_P(
    CustomDenseChunkedMatrix,
    CustomDenseChunkedMatrixFullTest,
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

class CustomDenseChunkedMatrixBlockTest :
    public ::testing::TestWithParam<std::tuple<typename CustomDenseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CustomDenseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomDenseChunkedMatrixBlockTest, Basic) {
    auto tparam = GetParam();
    auto opt = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto block = std::get<2>(tparam);
    tatami_test::test_block_access(*simple_mat, *ref, block.first, block.second, opt);
    tatami_test::test_block_access(*subset_mat, *ref, block.first, block.second, opt);
}

INSTANTIATE_TEST_SUITE_P(
    CustomDenseChunkedMatrix,
    CustomDenseChunkedMatrixBlockTest,
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
            std::make_pair(0.15, 0.72),
            std::make_pair(0.38, 0.62)
        )
    )
);

/*******************************************************/

class CustomDenseChunkedMatrixIndexTest :
    public ::testing::TestWithParam<std::tuple<typename CustomDenseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessOptions, std::pair<double, double> > >, 
    public CustomDenseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomDenseChunkedMatrixIndexTest, Basic) {
    auto tparam = GetParam();
    auto opt = tatami_test::convert_test_access_options(std::get<1>(tparam));
    auto index = std::get<2>(tparam);
    tatami_test::test_indexed_access(*simple_mat, *ref, index.first, index.second, opt);
    tatami_test::test_indexed_access(*subset_mat, *ref, index.first, index.second, opt);
}

INSTANTIATE_TEST_SUITE_P(
    CustomDenseChunkedMatrix,
    CustomDenseChunkedMatrixIndexTest,
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
            std::make_pair(0.0, 0.1),
            std::make_pair(0.2, 0.2),
            std::make_pair(0.7, 0.3)
        )
    )
);
