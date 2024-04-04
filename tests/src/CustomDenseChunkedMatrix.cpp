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
        int   // cache size
    > SimulationParameters;

protected:
    inline static std::unique_ptr<tatami::Matrix<double, int> > ref, mock_mat, simple_mat, subset_mat;

    typedef tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> > DChunk;
    typedef tatami_chunked::MockSimpleDenseChunk MockSimple;
    typedef tatami_chunked::MockSubsetDenseChunk MockSubset;

    inline static SimulationParameters last_params;

    static void assemble(const SimulationParameters& params) {
        if (ref && params == last_params) {
            return;
        }
        last_params = params;

        auto matdim = std::get<0>(params);
        auto chunkdim = std::get<1>(params);
        bool rowmajor = std::get<2>(params);
        size_t cache_size = std::get<3>(params);

        auto full = tatami_test::simulate_dense_vector<double>(matdim.first * matdim.second, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor + cache_size);
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        auto num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        auto num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<DChunk> mock_chunks(num_chunks_per_row * num_chunks_per_column);
        std::vector<MockSimple> simple_chunks(mock_chunks.size());
        std::vector<MockSubset> subset_chunks(mock_chunks.size());

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
                subset_chunks[offset] = MockSubset(std::move(contents), chunkdim.first, chunkdim.second);
            }
        }

        tatami_chunked::CustomDenseChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mock_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, DChunk>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(mock_chunks), rowmajor, opt
        ));

        simple_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, MockSimple>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(simple_chunks), rowmajor, opt
        ));

        subset_mat.reset(new tatami_chunked::CustomDenseChunkedMatrix<double, int, MockSubset>(
            matdim.first, matdim.second, chunkdim.first, chunkdim.second, std::move(subset_chunks), rowmajor, opt
        ));
    }
};

/*******************************************************/

class CustomDenseChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<typename CustomDenseChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessParameters> >, 
    public CustomDenseChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomDenseChunkedMatrixFullTest, Basic) {
    auto params = tatami_test::convert_access_parameters(std::get<1>(GetParam()));
    tatami_test::test_full_access(params, mock_mat.get(), ref.get());
    tatami_test::test_full_access(params, simple_mat.get(), ref.get());
    tatami_test::test_full_access(params, subset_mat.get(), ref.get());
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
            ::testing::Values(0, 1000, 10000) // cache size
        ),

        tatami_test::standard_test_access_parameter_combinations()
    )
);

///*******************************************************/
//
//class CustomDenseChunkedMatrixBlockTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
//    public CustomDenseChunkedMatrixMethods {};
//
//TEST_P(CustomDenseChunkedMatrixBlockTest, Row) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    auto bounds = std::get<5>(param);
//
//    assemble(sparse, matdim, chunkdim, rowmajor, cache_size);
//
//    bool FORWARD = true;
//    size_t JUMP = 1;
//    int FIRST = bounds.first * ref->ncol();
//    int LAST = bounds.second * ref->ncol();
//
//    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
//    tatami_test::test_sliced_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
//}
//
//TEST_P(CustomDenseChunkedMatrixBlockTest, Column) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    auto bounds = std::get<5>(param);
//
//    assemble(sparse, matdim, chunkdim, rowmajor, cache_size);
//
//    bool FORWARD = true;
//    size_t JUMP = 1;
//    int FIRST = bounds.first * ref->nrow();
//    int LAST = bounds.second * ref->nrow();
//
//    tatami_test::test_sliced_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
//    tatami_test::test_sliced_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, LAST);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomDenseChunkedMatrix,
//    CustomDenseChunkedMatrixBlockTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(201, 67),
//            std::make_pair(123, 372)
//        ),
//
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(10, 10)
//        ),
//
//        ::testing::Values(true, false), // row major
//        ::testing::Values(false, true), // sparse chunks
//        ::testing::Values(0, 1000, 10000), // cache size
//
//        ::testing::Values( // block boundaries
//            std::make_pair(0.0, 0.35),
//            std::make_pair(0.15, 0.87),
//            std::make_pair(0.38, 1.0)
//        )
//    )
//);
//
///*******************************************************/
//
//class CustomDenseChunkedMatrixIndexTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
//    public CustomDenseChunkedMatrixMethods {};
//
//TEST_P(CustomDenseChunkedMatrixIndexTest, Row) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    auto bounds = std::get<5>(param);
//
//    assemble(sparse, matdim, chunkdim, rowmajor, cache_size);
//
//    bool FORWARD = true;
//    size_t JUMP = 1;
//    int FIRST = bounds.first * ref->ncol(), STEP = bounds.second;
//
//    tatami_test::test_indexed_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
//    tatami_test::test_indexed_row_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
//}
//
//TEST_P(CustomDenseChunkedMatrixIndexTest, Column) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    auto bounds = std::get<5>(param);
//
//    assemble(sparse, matdim, chunkdim, rowmajor, cache_size);
//
//    bool FORWARD = true;
//    size_t JUMP = 1;
//    int FIRST = bounds.first * ref->nrow(), STEP = bounds.second;
//
//    tatami_test::test_indexed_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
//    tatami_test::test_indexed_column_access(mat.get(), ref.get(), FORWARD, JUMP, FIRST, STEP);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomDenseChunkedMatrix,
//    CustomDenseChunkedMatrixIndexTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(198, 67),
//            std::make_pair(187, 300)
//        ),
//
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(7, 13)
//        ),
//
//        ::testing::Values(true, false), // row major 
//        ::testing::Values(false, true), // sparse chunks
//        ::testing::Values(0, 1000, 10000), // cache size
//
//        ::testing::Values( // index information.
//            std::make_pair(0.0, 10),
//            std::make_pair(0.2, 5),
//            std::make_pair(0.7, 3)
//        )
//    )
//);
//
///*******************************************************/
//
//class CustomDenseChunkedMatrixOracleFullTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool> >, 
//    public CustomDenseChunkedMatrixMethods {};
//
//TEST_P(CustomDenseChunkedMatrixOracleFullTest, Row) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true);
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false);
//}
//
//TEST_P(CustomDenseChunkedMatrixOracleFullTest, Column) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true);
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomDenseChunkedMatrix,
//    CustomDenseChunkedMatrixOracleFullTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(200, 50),
//            std::make_pair(100, 300)
//        ),
//
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(11, 13) // odd numbers
//        ),
//
//        ::testing::Values(true, false), // row major
//        ::testing::Values(false, true), // sparse chunks
//        ::testing::Values(1000, 10000), // cache size
//        ::testing::Values(false, true)  // subsetted oracle
//    )
//);
//
///*******************************************************/
//
//class CustomDenseChunkedMatrixOracleBlockTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool, std::pair<double, double> > >, 
//    public CustomDenseChunkedMatrixMethods {};
//
//TEST_P(CustomDenseChunkedMatrixOracleBlockTest, Row) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//    auto bounds = std::get<6>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//
//    int FIRST = bounds.first * ref->ncol();
//    int LAST = bounds.second * ref->ncol();
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true, FIRST, LAST - FIRST);
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false, FIRST, LAST - FIRST);
//}
//
//TEST_P(CustomDenseChunkedMatrixOracleBlockTest, Column) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//    auto bounds = std::get<6>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//
//    int FIRST = bounds.first * ref->nrow();
//    int LAST = bounds.second * ref->nrow();
//
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true, FIRST, LAST - FIRST);
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false, FIRST, LAST - FIRST);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomDenseChunkedMatrix,
//    CustomDenseChunkedMatrixOracleBlockTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(201, 67),
//            std::make_pair(123, 372)
//        ),
//
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(10, 10)
//        ),
//
//        ::testing::Values(true, false), // row major
//        ::testing::Values(false, true), // sparse chunks
//        ::testing::Values(1000, 10000), // cache size
//        ::testing::Values(false, true), // subsetted oracle
//
//        ::testing::Values( // block boundaries
//            std::make_pair(0.0, 0.35),
//            std::make_pair(0.15, 0.87),
//            std::make_pair(0.38, 1.0)
//        )
//    )
//);
//
///*******************************************************/
//
//class CustomDenseChunkedMatrixOracleIndexTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool, std::pair<double, double> > >, 
//    public CustomDenseChunkedMatrixMethods {};
//
//TEST_P(CustomDenseChunkedMatrixOracleIndexTest, Row) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//    auto bounds = std::get<6>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//
//    int NC = ref->ncol();
//    int FIRST = bounds.first * NC, STEP = bounds.second;
//    std::vector<int> indices;
//    {
//        int counter = FIRST;
//        while (counter < NC) {
//            indices.push_back(counter);
//            counter += STEP;
//        }
//    }
//
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ true, indices);
//    tatami_test::test_oracle_row_access(mat.get(), ref.get(), /* random = */ false, indices);
//}
//
//TEST_P(CustomDenseChunkedMatrixOracleIndexTest, Column) {
//    auto param = GetParam();
//    auto matdim = std::get<0>(param);
//    auto chunkdim = std::get<1>(param);
//    bool rowmajor = std::get<2>(param);
//    bool sparse = std::get<3>(param);
//    auto cache_size = std::get<4>(param);
//    bool subset_oracle = std::get<5>(param);
//    auto bounds = std::get<6>(param);
//
//    assemble2(sparse, subset_oracle, matdim, chunkdim, rowmajor, cache_size);
//
//    int NR = ref->nrow();
//    int FIRST = bounds.first * NR, STEP = bounds.second;
//    std::vector<int> indices;
//    {
//        int counter = FIRST;
//        while (counter < NR) {
//            indices.push_back(counter);
//            counter += STEP;
//        }
//    }
//
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ true, indices);
//    tatami_test::test_oracle_column_access(mat.get(), ref.get(), /* random = */ false, indices);
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    CustomDenseChunkedMatrix,
//    CustomDenseChunkedMatrixOracleIndexTest,
//    ::testing::Combine(
//        ::testing::Values( // matrix dimensions
//            std::make_pair(198, 67),
//            std::make_pair(187, 300)
//        ),
//
//        ::testing::Values( // chunk dimensions
//            std::make_pair(1, 20),
//            std::make_pair(20, 1),
//            std::make_pair(7, 13)
//        ),
//
//        ::testing::Values(true, false), // row major 
//        ::testing::Values(false, true), // sparse chunks
//        ::testing::Values(1000, 10000), // cache size
//        ::testing::Values(false, true), // subsetted oracle
//
//        ::testing::Values( // index information.
//            std::make_pair(0.0, 10),
//            std::make_pair(0.2, 5),
//            std::make_pair(0.7, 3)
//        )
//    )
//);
