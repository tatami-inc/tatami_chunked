#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"
#include "tatami_chunked/CustomChunkedMatrix.hpp"
#include "tatami_chunked/simple_chunk_wrappers.hpp"

#include "mock_chunk.h"

class CustomChunkedMatrixCore {
protected:
    static std::unique_ptr<tatami::Matrix<double, int> > ref, mat, subset_mat;

    typedef std::tuple<
        bool, // sparse
        std::pair<int, int>, // matrix dimensions
        std::pair<int, int>, // chunk dimensions
        bool, // row major chunks
        int   // cache size
    > SimulationParameters;

    typedef tatami_chunked::SimpleDenseChunkWrapper<MockDenseChunk<true> > DChunk;
    typedef tatami_chunked::SimpleSparseChunkWrapper<MockSparseChunk<false> > SChunk;

    static void dense_assemble(const SimulationParameters& params) {
        auto matdim = std::get<1>(params);
        auto chunkdim = std::get<2>(params);
        bool rowmajor = std::get<3>(params);
        size_t cache_size = std::get<4>(params);

        auto full = tatami_test::simulate_dense_vector<double>(matdim.first * matdim.second, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor + cache_size);
        ref.reset(new tatami::DenseRowMatrix<double, int>(matdim.first, matdim.second, std::move(full)));

        auto num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        auto num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<DChunk> chunks(num_chunks_per_row * num_chunks_per_column);
        std::vector<DChunk> sub_chunks = chunks;

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
                chunks[offset] = DChunk(MockDenseChunk<true>(chunkdim.first, chunkdim.second, contents));
                sub_chunks[offset] = DChunk(MockDenseChunk<true>(chunkdim.first, chunkdim.second, std::move(contents)));
            }
        }

        tatami_chunked::CustomChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mat.reset(new tatami_chunked::CustomChunkedDenseMatrix<false, double, int, DChunk>(
            matdim.first,
            matdim.second,
            chunkdim.first,
            chunkdim.second,
            std::move(chunks),
            rowmajor,
            opt
        ));

        subset_mat.reset(new tatami_chunked::CustomChunkedDenseMatrix<true, double, int, DChunk>(
            matdim.first,
            matdim.second,
            chunkdim.first,
            chunkdim.second,
            std::move(sub_chunks),
            rowmajor,
            opt
        ));
    }

    template<bool use_subsetted_oracle_ = false>
    static void sparse_assemble(const SimulationParameters& params) {
        auto matdim = std::get<1>(params);
        auto chunkdim = std::get<2>(params);
        bool rowmajor = std::get<3>(params);
        size_t cache_size = std::get<4>(params);

        auto full = tatami_test::simulate_sparse_compressed<double>(matdim.first, matdim.second, 0.1, -10, 10, 
            /* seed = */ matdim.first * matdim.second + chunkdim.first * chunkdim.second + rowmajor + cache_size);
        ref.reset(new tatami::CompressedSparseColumnMatrix<double, int>(matdim.first, matdim.second, std::move(full.value), std::move(full.index), std::move(full.ptr)));

        auto num_chunks_per_row = (matdim.second + chunkdim.second - 1) / chunkdim.second;
        auto num_chunks_per_column = (matdim.first + chunkdim.first - 1) / chunkdim.first;
        std::vector<SChunk> chunks(num_chunks_per_row * num_chunks_per_column);
        std::vector<SChunk> sub_chunks = chunks;

        for (int r = 0; r < num_chunks_per_column; ++r) {
            for (int c = 0; c < num_chunks_per_row; ++c) {
                auto cstart = c * chunkdim.second;
                auto cend = std::min(cstart + chunkdim.second, matdim.second);
                auto clen = cend - cstart;

                auto rstart = r * chunkdim.first;
                auto rend = std::min(rstart + chunkdim.first, matdim.first);
                auto rlen = rend - rstart;

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

                auto offset = rowmajor ? (r * num_chunks_per_row + c) : (c * num_chunks_per_column + r);
                chunks[offset] = SChunk(MockSparseChunk<false>(chunkdim.first, chunkdim.second, vcontents, icontents, pcontents));
                sub_chunks[offset] = SChunk(MockSparseChunk<false>(chunkdim.first, chunkdim.second, std::move(vcontents), std::move(icontents), std::move(pcontents)));
            }
        }

        tatami_chunked::CustomChunkedOptions opt;
        opt.maximum_cache_size = cache_size;
        opt.require_minimum_cache = false;

        mat.reset(new tatami_chunked::CustomChunkedSparseMatrix<false, double, int, SChunk>(
            matdim.first,
            matdim.second,
            chunkdim.first,
            chunkdim.second,
            std::move(chunks),
            rowmajor,
            opt
        ));

        subset_mat.reset(new tatami_chunked::CustomChunkedSparseMatrix<true, double, int, SChunk>(
            matdim.first,
            matdim.second,
            chunkdim.first,
            chunkdim.second,
            std::move(sub_chunks),
            rowmajor,
            opt
        ));
    }

    static SimulationParameters last_params;

    static void assemble(const SimulationParameters& params) {
        if (ref && params == last_params) {
            return;
        }
        last_params = params;
        if (std::get<0>(params)) {
            sparse_assemble(params);
        } else {
            dense_assemble(params);
        }
    }
};

/*******************************************************/

class CustomChunkedMatrixFullTest :
    public ::testing::TestWithParam<std::tuple<typename CustomChunkedMatrixCore::SimulationParameters, tatami_test::StandardTestAccessParameters> >, 
    public CustomChunkedMatrixCore {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(CustomChunkedMatrixFullTest, Column) {
    auto params = tatami_test::convert_access_parameters(std::get<1>(GetParam()));
    tatami_test::test_full_access(params, mat.get(), ref.get());
    if (params.use_oracle) {
        tatami_test::test_full_access(params, subset_mat.get(), ref.get());
    }
}

INSTANTIATE_TEST_SUITE_P(
    CustomChunkedMatrix,
    CustomChunkedMatrixFullTest,
    ::testing::Combine(
        ::testing::Combine(
            ::testing::Values(false, true), // sparse chunks

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
//class CustomChunkedMatrixBlockTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
//    public CustomChunkedMatrixMethods {};
//
//TEST_P(CustomChunkedMatrixBlockTest, Row) {
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
//TEST_P(CustomChunkedMatrixBlockTest, Column) {
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
//    CustomChunkedMatrix,
//    CustomChunkedMatrixBlockTest,
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
//class CustomChunkedMatrixIndexTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, std::pair<double, double> > >, 
//    public CustomChunkedMatrixMethods {};
//
//TEST_P(CustomChunkedMatrixIndexTest, Row) {
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
//TEST_P(CustomChunkedMatrixIndexTest, Column) {
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
//    CustomChunkedMatrix,
//    CustomChunkedMatrixIndexTest,
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
//class CustomChunkedMatrixOracleFullTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool> >, 
//    public CustomChunkedMatrixMethods {};
//
//TEST_P(CustomChunkedMatrixOracleFullTest, Row) {
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
//TEST_P(CustomChunkedMatrixOracleFullTest, Column) {
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
//    CustomChunkedMatrix,
//    CustomChunkedMatrixOracleFullTest,
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
//class CustomChunkedMatrixOracleBlockTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool, std::pair<double, double> > >, 
//    public CustomChunkedMatrixMethods {};
//
//TEST_P(CustomChunkedMatrixOracleBlockTest, Row) {
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
//TEST_P(CustomChunkedMatrixOracleBlockTest, Column) {
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
//    CustomChunkedMatrix,
//    CustomChunkedMatrixOracleBlockTest,
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
//class CustomChunkedMatrixOracleIndexTest :
//    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, bool, bool, int, bool, std::pair<double, double> > >, 
//    public CustomChunkedMatrixMethods {};
//
//TEST_P(CustomChunkedMatrixOracleIndexTest, Row) {
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
//TEST_P(CustomChunkedMatrixOracleIndexTest, Column) {
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
//    CustomChunkedMatrix,
//    CustomChunkedMatrixOracleIndexTest,
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
