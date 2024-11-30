#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"

#include "tatami_chunked/CustomSparseChunkedMatrix.hpp"
#include "tatami_chunked/mock_sparse_chunk.hpp"

#include "mock_blob.h"

class CustomSparseChunkedMatrixCore {
public:
    typedef std::tuple<
        std::pair<int, int>, // matrix dimensions
        std::pair<int, int>, // chunk dimensions
        bool, // row major chunks
        double // cache fraction
    > SimulationParameters;

protected:
    inline static std::unique_ptr<tatami::Matrix<double, int> > ref, mock_mat, simple_mat, subset_mat;

    typedef tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> > SChunk;
    typedef tatami_chunked::MockSimpleSparseChunk MockSimple;
    typedef tatami_chunked::MockSubsettedSparseChunk MockSubsetted;

    inline static SimulationParameters last_params;

    static void assemble(const SimulationParameters& params) {
        if (ref && params == last_params) {
            return;
        }
        last_params = params;

        auto matdim = std::get<0>(params);
        auto chunkdim = std::get<1>(params);
        bool rowmajor = std::get<2>(params);
        double cache_fraction = std::get<3>(params);

        auto full = tatami_test::simulate_compressed_sparse<double, int>(matdim.second, matdim.first, [&]{
            tatami_test::SimulateCompressedSparseOptions opt;
            opt.density = 0.1;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor + 100 * cache_fraction;
            return opt;
        }());

        ref.reset(new tatami::CompressedSparseColumnMatrix<double, int>(
            matdim.first,
            matdim.second,
            std::move(full.data),
            std::move(full.index),
            std::move(full.indptr)
        ));

        auto num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        auto num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<SChunk> mock_chunks(num_chunks_per_row * num_chunks_per_column);
        std::vector<MockSimple> simple_chunks(mock_chunks.size());
        std::vector<MockSubsetted> subset_chunks(mock_chunks.size());

        for (int r = 0; r < num_chunks_per_column; ++r) {
            for (int c = 0; c < num_chunks_per_row; ++c) {
                auto cstart = c * chunkdim.second;
                auto cend = std::min(cstart + chunkdim.second, matdim.second);
                auto clen = cend - cstart;

                auto rstart = r * chunkdim.first;
                auto rend = std::min(rstart + chunkdim.first, matdim.first);
                auto rlen = rend - rstart;

                auto offset = rowmajor ? (r * num_chunks_per_row + c) : (c * num_chunks_per_column + r);

                {
                    std::vector<double> vcontents;
                    std::vector<int> icontents;
                    std::vector<size_t> pcontents(1);

                    auto ext = ref->sparse_column(rstart, rlen);
                    std::vector<double> vbuffer(rlen);
                    std::vector<int> ibuffer(rlen);

                    for (int c2 = 0; c2 < clen; ++c2) {
                        auto range = ext->fetch(c2 + cstart, vbuffer.data(), ibuffer.data());
                        vcontents.insert(vcontents.end(), range.value, range.value + range.number);
                        for (int i = 0; i < range.number; ++i) {
                            icontents.push_back(range.index[i] - rstart);
                        }
                        pcontents.push_back(pcontents.back() + range.number);
                    }

                    mock_chunks[offset] = SChunk(MockSparseBlob<false>(rlen, clen, std::move(vcontents), std::move(icontents), std::move(pcontents)));
                }

                {
                    std::vector<double> vcontents;
                    std::vector<int> icontents;
                    std::vector<size_t> pcontents(1);

                    auto ext = ref->sparse_row(cstart, clen);
                    std::vector<double> vbuffer(clen);
                    std::vector<int> ibuffer(clen);

                    for (int r2 = 0; r2 < rlen; ++r2) {
                        auto range = ext->fetch(r2 + rstart, vbuffer.data(), ibuffer.data());
                        vcontents.insert(vcontents.end(), range.value, range.value + range.number);
                        for (int i = 0; i < range.number; ++i) {
                            icontents.push_back(range.index[i] - cstart);
                        }
                        pcontents.push_back(pcontents.back() + range.number);
                    }

                    simple_chunks[offset] = MockSimple(rlen, clen, vcontents, icontents, pcontents);
                    subset_chunks[offset] = MockSubsetted(rlen, clen, std::move(vcontents), std::move(icontents), std::move(pcontents));
                }
            }
        }

        tatami_chunked::CustomSparseChunkedMatrixOptions opt;
        size_t cache_size = static_cast<double>(matdim.first) * static_cast<double>(matdim.second) * cache_fraction * static_cast<double>(sizeof(double) + sizeof(int));
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = (cache_size > 0);

        mock_mat.reset(new tatami_chunked::CustomSparseChunkedMatrix<double, int, SChunk>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(mock_chunks), rowmajor, opt
        ));

        simple_mat.reset(new tatami_chunked::CustomSparseChunkedMatrix<double, int, MockSimple>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(simple_chunks), rowmajor, opt
        ));

        subset_mat.reset(new tatami_chunked::CustomSparseChunkedMatrix<double, int, MockSubsetted>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(subset_chunks), rowmajor, opt
        ));
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
    tatami_test::test_full_access(*mock_mat, *ref, opt);
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

            ::testing::Values(true, false), // row major
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
    tatami_test::test_block_access(*mock_mat, *ref, block.first, block.second, opt);
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

            ::testing::Values(true, false), // row major
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
    tatami_test::test_indexed_access(*mock_mat, *ref, index.first, index.second, opt);
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

            ::testing::Values(true, false), // row major 
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
