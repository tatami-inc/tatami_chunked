#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"

#include "tatami_chunked/CustomDenseChunkedMatrix.hpp"
#include "tatami_chunked/mock_dense_chunk.hpp"

#include "mock_blob.h"

struct CustomDenseChunkedMatrixCore {
    typedef std::tuple<
        std::pair<int, int>, // matrix dimensions
        std::pair<int, int>, // chunk dimensions
        bool, // row major chunks
        double // cache fraction
    > SimulationParameters;

protected:
    inline static std::unique_ptr<tatami::Matrix<double, int> > ref, mock_mat, simple_mat, subset_mat;

    typedef tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> > DChunk;
    typedef tatami_chunked::MockSimpleDenseChunk MockSimple;
    typedef tatami_chunked::MockSubsettedDenseChunk MockSubsetted;

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

        auto full = tatami_test::simulate_vector<double>(matdim.first * matdim.second, [&]{
            tatami_test::SimulateVectorOptions opt;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor + 100 * cache_fraction;
            return opt;
        }());
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        auto num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        auto num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<DChunk> mock_chunks(num_chunks_per_row * num_chunks_per_column);
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

                std::vector<double> contents(chunkdim.second * chunkdim.first);
                auto ccptr = contents.data();
                auto ext = ref->dense_row(cstart, clen);
                for (int r2 = 0; r2 < rlen; ++r2) {
                    auto ptr = ext->fetch(r2 + rstart, ccptr);
                    tatami::copy_n(ptr, clen, ccptr);
                    ccptr += chunkdim.second;
                }

                auto offset = rowmajor ? (r * num_chunks_per_row + c) : (c * num_chunks_per_column + r);
                mock_chunks[offset] = DChunk(MockDenseBlob<true>(chunkdim.first, chunkdim.second, contents));
                simple_chunks[offset] = MockSimple(contents, chunkdim.first, chunkdim.second);
                subset_chunks[offset] = MockSubsetted(std::move(contents), chunkdim.first, chunkdim.second);
            }
        }

        tatami_chunked::CustomDenseChunkedMatrixOptions opt;
        size_t cache_size = static_cast<double>(matdim.first) * static_cast<double>(matdim.second) * cache_fraction * static_cast<double>(sizeof(double));
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = (cache_size > 0);

        mock_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, DChunk>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(mock_chunks), rowmajor, opt
        ));

        simple_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, MockSimple>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(simple_chunks), rowmajor, opt
        ));

        subset_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, MockSubsetted>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(subset_chunks), rowmajor, opt
        ));
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
    tatami_test::test_full_access(*mock_mat, *ref, opts);
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

            ::testing::Values(true, false), // row major
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
    tatami_test::test_block_access(*mock_mat, *ref, block.first, block.second, opt);
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

            ::testing::Values(true, false), // row major
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
    tatami_test::test_indexed_access(*mock_mat, *ref, index.first, index.second, opt);
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

            ::testing::Values(true, false), // row major 
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
