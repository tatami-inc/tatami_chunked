#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"
#include "tatami_chunked/mock_sparse_chunk.hpp"
#include "mock_blob.h"

static std::vector<int>& remove_shift(std::vector<int>& buffer, int shift) {
    for (auto& i : buffer) {
        i -= shift;
    }
    return buffer;
}

static tatami::VectorPtr<int> create_indices(int dim, const std::pair<double, int>& config) {
    auto output = std::make_shared<std::vector<int> >();
    auto& indices = *output;
    int start = config.first * dim, step = config.second;
    while (start < dim) {
        indices.push_back(start);
        start += step;
    } 
    return output;
}

/**********************************************************
 **********************************************************/

class MockSimpleSparseChunkUtils {
protected:
    static inline std::unique_ptr<tatami::Matrix<double, int> > ref;

    static inline tatami_chunked::MockSimpleSparseChunk mock;
    static inline tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> > drm_chunk;
    static inline tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> > dcm_chunk;

    static inline std::pair<int, int> last_params;

    static void assemble(const std::pair<int, int>& params) {
        if (ref && last_params == params) {
            return;
        }
        last_params = params;

        const auto& dim = params;
        auto full = tatami_test::simulate_sparse_vector<double>(dim.first * dim.second, 0.2, -10, 10, /* seed = */ dim.first * dim.second);
        tatami::DenseRowMatrix<double, int> tmp(dim.first, dim.second, std::move(full));

        auto compressed = tatami::retrieve_compressed_sparse_contents<true, double, int>(&tmp, true);
        mock = tatami_chunked::MockSimpleSparseChunk(
            dim.first, 
            dim.second, 
            compressed.value, 
            compressed.index, 
            compressed.pointers
        );

        drm_chunk = tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >(MockSparseBlob<true>(
            dim.first, 
            dim.second, 
            compressed.value,
            compressed.index,
            compressed.pointers
        ));

        ref.reset(new tatami::CompressedSparseRowMatrix<double, int>(dim.first, dim.second, std::move(compressed.value), std::move(compressed.index), std::move(compressed.pointers)));

        auto compressed2 = tatami::retrieve_compressed_sparse_contents<false, double, int>(ref.get(), true);
        dcm_chunk = tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >(MockSparseBlob<false>(
            dim.first, 
            dim.second, 
            std::move(compressed2.value),
            std::move(compressed2.index),
            std::move(compressed2.pointers)
        ));
    }
};

/**********************************************************
 **********************************************************/

class MockSimpleSparseChunkBlockTest : public MockSimpleSparseChunkUtils, public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double> > > {
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSimpleSparseChunkBlockTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto block = std::get<1>(params);

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    // First extracting row-wise blocks.
    int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;
    {
        std::vector<std::vector<double> > mock_values(dim.first), drm_values(dim.first), dcm_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first), drm_indices(dim.first), dcm_indices(dim.first);
        mock.template extract<true>(c_first, c_len, mock_work, mock_values, mock_indices, 7);
        drm_chunk.template extract<true>(c_first, c_len, drm_work, drm_values, drm_indices, 3);
        dcm_chunk.template extract<true>(c_first, c_len, dcm_work, dcm_values, dcm_indices, 5);

        auto ext = this->ref->sparse_row(c_first, c_len);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.value, drm_values[r]);
            EXPECT_EQ(refrow.value, dcm_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 7));
            EXPECT_EQ(refrow.index, remove_shift(drm_indices[r], 3));
            EXPECT_EQ(refrow.index, remove_shift(dcm_indices[r], 5));
        }
    }

    // Then extracting column-wise blocks.
    int r_first = block.first * dim.first, r_last = block.second * dim.first, r_len = r_last - r_first;
    {
        std::vector<std::vector<double> > mock_values(dim.second), drm_values(dim.second), dcm_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second), drm_indices(dim.second), dcm_indices(dim.second);
        mock.template extract<false>(r_first, r_len, mock_work, mock_values, mock_indices, 1);
        drm_chunk.template extract<false>(r_first, r_len, drm_work, drm_values, drm_indices, 6);
        dcm_chunk.template extract<false>(r_first, r_len, dcm_work, dcm_values, dcm_indices, 10);

        auto ext = this->ref->sparse_column(r_first, r_len);
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.value, drm_values[c]);
            EXPECT_EQ(refcol.value, dcm_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 1));
            EXPECT_EQ(refcol.index, remove_shift(drm_indices[c], 6));
            EXPECT_EQ(refcol.index, remove_shift(dcm_indices[c], 10));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    MockSimpleSparseChunk,
    MockSimpleSparseChunkBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23),
            std::make_pair(190, 150)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 1.0),
            std::make_pair(0.25, 0.75),
            std::make_pair(0.33, 0.66)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSimpleSparseChunkIndexTest : public MockSimpleSparseChunkUtils, public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int> > > {
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSimpleSparseChunkIndexTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto index = std::get<1>(params);

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    // First extracting row-wise blocks.
    {
        auto c_indices = create_indices(dim.second, index);

        std::vector<std::vector<double> > mock_values(dim.first), drm_values(dim.first), dcm_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first), drm_indices(dim.first), dcm_indices(dim.first);
        mock.template extract<true>(*c_indices, mock_work, mock_values, mock_indices, 17);
        drm_chunk.template extract<true>(*c_indices, drm_work, drm_values, drm_indices, 72);
        dcm_chunk.template extract<true>(*c_indices, dcm_work, dcm_values, dcm_indices, 5);

        auto ext = this->ref->sparse_row(c_indices);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.value, drm_values[r]);
            EXPECT_EQ(refrow.value, dcm_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 17));
            EXPECT_EQ(refrow.index, remove_shift(drm_indices[r], 72));
            EXPECT_EQ(refrow.index, remove_shift(dcm_indices[r], 5));
        }
    }

    // Then extracting column-wise blocks.
    {
        auto r_indices = create_indices(dim.first, index);

        std::vector<std::vector<double> > mock_values(dim.second), drm_values(dim.second), dcm_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second), drm_indices(dim.second), dcm_indices(dim.second);
        mock.template extract<false>(*r_indices, mock_work, mock_values, mock_indices, 17);
        drm_chunk.template extract<false>(*r_indices, drm_work, drm_values, drm_indices, 72);
        dcm_chunk.template extract<false>(*r_indices, dcm_work, dcm_values, dcm_indices, 5);

        auto ext = this->ref->sparse_column(r_indices);
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.value, drm_values[c]);
            EXPECT_EQ(refcol.value, dcm_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 17));
            EXPECT_EQ(refcol.index, remove_shift(drm_indices[c], 72));
            EXPECT_EQ(refcol.index, remove_shift(dcm_indices[c], 5));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    MockSimpleSparseChunk,
    MockSimpleSparseChunkIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23),
            std::make_pair(190, 150)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 10),
            std::make_pair(0.25, 5),
            std::make_pair(0.33, 3)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsetSparseChunkUtils {
protected:
    static inline std::unique_ptr<tatami::Matrix<double, int> > ref;

    static inline tatami_chunked::MockSubsetSparseChunk mock;

    static inline std::pair<int, int> last_params;

    static void assemble(const std::pair<int, int>& params) {
        if (ref && last_params == params) {
            return;
        }
        last_params = params;

        const auto& dim = params;
        auto full = tatami_test::simulate_sparse_vector<double>(dim.first * dim.second, 0.2, -10, 10, /* seed = */ dim.first * dim.second + 1);
        tatami::DenseRowMatrix<double, int> tmp(dim.first, dim.second, std::move(full));

        auto compressed = tatami::retrieve_compressed_sparse_contents<true, double, int>(&tmp, true);
        mock = tatami_chunked::MockSubsetSparseChunk(
            dim.first, 
            dim.second, 
            compressed.value, 
            compressed.index, 
            compressed.pointers
        );

        ref.reset(new tatami::CompressedSparseRowMatrix<double, int>(dim.first, dim.second, std::move(compressed.value), std::move(compressed.index), std::move(compressed.pointers)));
    }
};

/**********************************************************
 **********************************************************/

class MockSubsetSparseChunkBlockBlockTest : 
    public MockSubsetSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetSparseChunkBlockBlockTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsetSparseChunk::Workspace work;

    // Row extraction.
    {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;

        std::vector<std::vector<double> > mock_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first);
        mock.template extract<true>(r_first, r_len, c_first, c_len, work, mock_values, mock_indices, 10);

        auto ext = ref->sparse_row(c_first, c_len);
        for (int r = r_first; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 10));
        }
    }

    // Column extraction.
    {
        int r_first = block.first * dim.first,  r_last = block.second * dim.first,  r_len = r_last - r_first;
        int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;

        std::vector<std::vector<double> > mock_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second);
        mock.template extract<false>(c_first, c_len, r_first, r_len, work, mock_values, mock_indices, 22);

        auto ext = ref->sparse_column(r_first, r_len);
        for (int c = c_first; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 22));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetSparseChunk,
    MockSubsetSparseChunkBlockBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 1.0),
            std::make_pair(0.2, 0.8),
            std::make_pair(0.4, 0.6)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 1.0),
            std::make_pair(0.25, 0.75),
            std::make_pair(0.33, 0.66)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsetSparseChunkBlockIndexTest : 
    public MockSubsetSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, int> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetSparseChunkBlockIndexTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto iparam = std::get<2>(params);

    typename tatami_chunked::MockSubsetSparseChunk::Workspace work;

    {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        auto c_indices = create_indices(dim.second, iparam);

        std::vector<std::vector<double> > mock_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first);
        mock.template extract<true>(r_first, r_len, *c_indices, work, mock_values, mock_indices, 10);

        auto ext = ref->sparse_row(c_indices);
        for (int r = r_first; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 10));
        }
    }

    {
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;
        auto r_indices = create_indices(dim.first, iparam);

        std::vector<std::vector<double> > mock_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second);
        mock.template extract<false>(c_first, c_len, *r_indices, work, mock_values, mock_indices, 5);

        auto ext = ref->sparse_column(r_indices);
        for (int c = c_first; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 5));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetSparseChunk,
    MockSubsetSparseChunkBlockIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 1.0),
            std::make_pair(0.2, 0.8),
            std::make_pair(0.4, 0.6)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 10),
            std::make_pair(0.25, 5),
            std::make_pair(0.33, 3)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsetSparseChunkIndexBlockTest : 
    public MockSubsetSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetSparseChunkIndexBlockTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsetSparseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, bounds);
        int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;

        std::vector<std::vector<double> > mock_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first);
        mock.template extract<true>(*r_indices, c_first, c_len, work, mock_values, mock_indices, 32);

        auto ext = ref->sparse_row(c_first, c_len);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 32));
        }
    }

    {
        auto c_indices = create_indices(dim.second, bounds);
        int r_first = block.first * dim.first,  r_last = block.second * dim.first,  r_len = r_last - r_first;

        std::vector<std::vector<double> > mock_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second);
        mock.template extract<false>(*c_indices, r_first, r_len, work, mock_values, mock_indices, 99);

        auto ext = ref->sparse_column(r_first, r_len);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 99));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetSparseChunk,
    MockSubsetSparseChunkIndexBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary indices 
            std::make_pair(0.0, 10),
            std::make_pair(0.2, 5),
            std::make_pair(0.4, 3)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 1.0),
            std::make_pair(0.25, 0.75),
            std::make_pair(0.33, 0.66)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsetSparseChunkIndexIndexTest : 
    public MockSubsetSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int>, std::pair<double, int> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetSparseChunkIndexIndexTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto iparam1 = std::get<1>(params);
    auto iparam2 = std::get<2>(params);

    typename tatami_chunked::MockSubsetSparseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, iparam1);
        auto c_indices = create_indices(dim.second, iparam2);

        std::vector<std::vector<double> > mock_values(dim.first);
        std::vector<std::vector<int> > mock_indices(dim.first);
        mock.template extract<true>(*r_indices, *c_indices, work, mock_values, mock_indices, 9);

        auto ext = ref->sparse_row(*c_indices);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_values[r]);
            EXPECT_EQ(refrow.index, remove_shift(mock_indices[r], 9));
        }
    }

    {
        auto c_indices = create_indices(dim.second, iparam1);
        auto r_indices = create_indices(dim.first, iparam2);

        std::vector<std::vector<double> > mock_values(dim.second);
        std::vector<std::vector<int> > mock_indices(dim.second);
        mock.template extract<false>(*c_indices, *r_indices, work, mock_values, mock_indices, 19);

        auto ext = ref->sparse_column(*r_indices);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_values[c]);
            EXPECT_EQ(refcol.index, remove_shift(mock_indices[c], 19));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetSparseChunk,
    MockSubsetSparseChunkIndexIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary indices
            std::make_pair(0.0, 10),
            std::make_pair(0.3, 5),
            std::make_pair(0.5, 3)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 10),
            std::make_pair(0.25, 5),
            std::make_pair(0.33, 3)
        )
    )
);
