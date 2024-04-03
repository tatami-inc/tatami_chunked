#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"
#include "tatami_chunked/mock_dense_chunk.hpp"
#include "mock_blob.h"

static std::vector<double> slice_and_fill(std::vector<double>& buffer, size_t primary, size_t stride, size_t len) {
    auto it = buffer.begin() + primary * stride;
    std::vector<double> output(it, it + len);
    std::fill(it, it + len, 0);
    return output;
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

static bool is_empty(const std::vector<double>& buffer) {
    for (auto& x : buffer) {
        if (x) {
            return false;
        }
    }
    return true;
}

/**********************************************************
 **********************************************************/

class MockSimpleDenseChunkUtils {
protected:
    static inline std::unique_ptr<tatami::Matrix<double, int> > ref;

    static inline tatami_chunked::MockSimpleDenseChunk mock;
    static inline tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> > drm_chunk;
    static inline tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<false> > dcm_chunk;

    static inline std::pair<int, int> last_params;

    void assemble(const std::pair<int, int>& params) {
        if (ref && last_params == params) {
            return;
        }
        last_params = params;

        const auto& dim = params;
        auto full = tatami_test::simulate_dense_vector<double>(dim.first * dim.second, -10, 10, /* seed = */ dim.first * dim.second);
        mock = tatami_chunked::MockSimpleDenseChunk(full, dim.first, dim.second);
        drm_chunk = tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> >(MockDenseBlob<true>(dim.first, dim.second, full));
        ref.reset(new tatami::DenseRowMatrix<double, int>(dim.first, dim.second, std::move(full)));

        std::vector<double> replacement(dim.first * dim.second);
        auto ext = ref->dense_column();
        for (int c = 0; c < dim.second; ++c) {
            auto dest = replacement.data() + c * dim.first;
            auto ptr = ext->fetch(c, dest);
            tatami::copy_n(ptr, dim.first, dest);
        }
        dcm_chunk = tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<false> >(MockDenseBlob<false>(dim.first, dim.second, std::move(replacement)));
    }
};

/**********************************************************
 **********************************************************/

class MockSimpleDenseChunkBlockTest : public MockSimpleDenseChunkUtils, public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double> > > {
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSimpleDenseChunkBlockTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto block = std::get<1>(params);

    typename tatami_chunked::MockSimpleDenseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<false> >::Workspace dcm_work;

    // First extracting row-wise blocks.
    int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;
    {
        int stride = std::max(100, dim.second);
        size_t buffer_size = dim.first * stride;
        std::vector<double> mock_buffer(buffer_size), drm_buffer(buffer_size), dcm_buffer(buffer_size);
        mock.template extract<true>(c_first, c_len, mock_work, mock_buffer.data(), stride);
        drm_chunk.template extract<true>(c_first, c_len, drm_work, drm_buffer.data(), stride);
        dcm_chunk.template extract<true>(c_first, c_len, dcm_work, dcm_buffer.data(), stride);

        auto ext = ref->dense_row(c_first, c_len);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow, slice_and_fill(mock_buffer, r, stride, c_len));
            EXPECT_EQ(refrow, slice_and_fill(drm_buffer, r, stride, c_len));
            EXPECT_EQ(refrow, slice_and_fill(dcm_buffer, r, stride, c_len));
        }

        EXPECT_TRUE(is_empty(mock_buffer));
        EXPECT_TRUE(is_empty(drm_buffer));
        EXPECT_TRUE(is_empty(dcm_buffer));
    }

    // Then extracting column-wise blocks.
    int r_first = block.first * dim.first, r_last = block.second * dim.first, r_len = r_last - r_first;
    {
        int stride = std::max(100, dim.first);
        size_t buffer_size = dim.second * stride;
        std::vector<double> mock_buffer(buffer_size), drm_buffer(buffer_size), dcm_buffer(buffer_size);
        mock.template extract<false>(r_first, r_len, mock_work, mock_buffer.data(), stride);
        drm_chunk.template extract<false>(r_first, r_len, drm_work, drm_buffer.data(), stride);
        dcm_chunk.template extract<false>(r_first, r_len, dcm_work, dcm_buffer.data(), stride);

        auto ext = ref->dense_column(r_first, r_len);
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol, slice_and_fill(mock_buffer, c, stride, r_len));
            EXPECT_EQ(refcol, slice_and_fill(drm_buffer, c, stride, r_len));
            EXPECT_EQ(refcol, slice_and_fill(dcm_buffer, c, stride, r_len));
        }

        EXPECT_TRUE(is_empty(mock_buffer));
        EXPECT_TRUE(is_empty(drm_buffer));
        EXPECT_TRUE(is_empty(dcm_buffer));
    }
}

INSTANTIATE_TEST_SUITE_P(
    MockSimpleDenseChunk,
    MockSimpleDenseChunkBlockTest,
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

class MockSimpleDenseChunkIndexTest : public MockSimpleDenseChunkUtils, public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int> > > {
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSimpleDenseChunkIndexTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto index = std::get<1>(params);

    typename tatami_chunked::MockSimpleDenseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleDenseChunkWrapper<MockDenseBlob<false> >::Workspace dcm_work;

    // First extracting row-wise blocks.
    {
        auto c_indices = create_indices(dim.second, index);

        int stride = std::max(100, dim.second);
        size_t buffer_size = dim.first * stride;
        std::vector<double> mock_buffer(buffer_size), drm_buffer(buffer_size), dcm_buffer(buffer_size);
        mock.template extract<true>(*c_indices, mock_work, mock_buffer.data(), stride);
        drm_chunk.template extract<true>(*c_indices, drm_work, drm_buffer.data(), stride);
        dcm_chunk.template extract<true>(*c_indices, dcm_work, dcm_buffer.data(), stride);

        auto ext = ref->dense_row(c_indices);
        size_t c_len = c_indices->size();
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow, slice_and_fill(mock_buffer, r, stride, c_len));
            EXPECT_EQ(refrow, slice_and_fill(drm_buffer, r, stride, c_len));
            EXPECT_EQ(refrow, slice_and_fill(dcm_buffer, r, stride, c_len));
        }

        EXPECT_TRUE(is_empty(mock_buffer));
        EXPECT_TRUE(is_empty(drm_buffer));
        EXPECT_TRUE(is_empty(dcm_buffer));
    }

    // Then extracting column-wise blocks.
    {
        auto r_indices = create_indices(dim.first, index);

        int stride = std::max(100, dim.first);
        size_t buffer_size = dim.second * stride;
        std::vector<double> mock_buffer(buffer_size), drm_buffer(buffer_size), dcm_buffer(buffer_size);
        mock.template extract<false>(*r_indices, mock_work, mock_buffer.data(), stride);
        drm_chunk.template extract<false>(*r_indices, drm_work, drm_buffer.data(), stride);
        dcm_chunk.template extract<false>(*r_indices, dcm_work, dcm_buffer.data(), stride);

        auto ext = ref->dense_column(r_indices);
        size_t r_len = r_indices->size();
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol, slice_and_fill(mock_buffer, c, stride, r_len));
            EXPECT_EQ(refcol, slice_and_fill(drm_buffer, c, stride, r_len));
            EXPECT_EQ(refcol, slice_and_fill(dcm_buffer, c, stride, r_len));
        }

        EXPECT_TRUE(is_empty(mock_buffer));
        EXPECT_TRUE(is_empty(drm_buffer));
        EXPECT_TRUE(is_empty(dcm_buffer));
    }
}

INSTANTIATE_TEST_SUITE_P(
    MockSimpleDenseChunk,
    MockSimpleDenseChunkIndexTest,
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

class MockSubsetDenseChunkUtils {
protected:
    static inline std::unique_ptr<tatami::Matrix<double, int> > ref;

    static inline tatami_chunked::MockSubsetDenseChunk mock;

    static inline std::pair<int, int> last_params;

    void assemble(const std::pair<int, int>& params) {
        if (ref && last_params == params) {
            return;
        }
        last_params = params;

        const auto& dim = params;
        auto full = tatami_test::simulate_dense_vector<double>(dim.first * dim.second, -10, 10, /* seed = */ dim.first * dim.second + 1);
        mock = tatami_chunked::MockSubsetDenseChunk(full, dim.first, dim.second);
        ref.reset(new tatami::DenseRowMatrix<double, int>(dim.first, dim.second, std::move(full)));
    }
};

/**********************************************************
 **********************************************************/

class MockSubsetDenseChunkBlockBlockTest : 
    public MockSubsetDenseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetDenseChunkBlockBlockTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsetDenseChunk::Workspace work;

    // Row extraction.
    {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;

        int stride = std::min(100, dim.second);
        std::vector<double> buffer(dim.first * stride);
        mock.template extract<true>(r_first, r_len, c_first, c_len, work, buffer.data(), stride);

        auto ext = ref->dense_row(c_first, c_len);
        for (int r = r_first; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow, slice_and_fill(buffer, r, stride, c_len));
        }

        EXPECT_TRUE(is_empty(buffer));
    }

    // Column extraction.
    {
        int r_first = block.first * dim.first,  r_last = block.second * dim.first,  r_len = r_last - r_first;
        int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;

        int stride = std::min(98, dim.first);
        std::vector<double> buffer(dim.second * stride);
        mock.template extract<false>(c_first, c_len, r_first, r_len, work, buffer.data(), stride);

        auto ext = ref->dense_column(r_first, r_len);
        for (int c = c_first; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol, slice_and_fill(buffer, c, stride, r_len));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetDenseChunk,
    MockSubsetDenseChunkBlockBlockTest,
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

class MockSubsetDenseChunkBlockIndexTest : 
    public MockSubsetDenseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, int> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetDenseChunkBlockIndexTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto iparam = std::get<2>(params);

    typename tatami_chunked::MockSubsetDenseChunk::Workspace work;

    {
        int r_first = bounds.first * dim.first,  r_last = bounds.second * dim.first,  r_len = r_last - r_first;
        auto c_indices = create_indices(dim.second, iparam);

        int stride = std::max(dim.second, 100);
        std::vector<double> buffer(dim.first * stride);
        mock.template extract<true>(r_first, r_len, *c_indices, work, buffer.data(), stride);

        auto ext = ref->dense_row(c_indices);
        for (int r = r_first; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_indices->size());
            EXPECT_EQ(refrow, slice_and_fill(buffer, r, stride, c_indices->size()));
        }

        EXPECT_TRUE(is_empty(buffer));
    }

    {
        int c_first = bounds.first * dim.second, c_last = bounds.second * dim.second, c_len = c_last - c_first;
        auto r_indices = create_indices(dim.first, iparam);

        int stride = std::max(dim.second, 87);
        std::vector<double> buffer(dim.second * stride);
        mock.template extract<false>(c_first, c_len, *r_indices, work, buffer.data(), stride);

        auto ext = ref->dense_column(r_indices);
        for (int c = c_first; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_indices->size());
            EXPECT_EQ(refcol, slice_and_fill(buffer, c, stride, r_indices->size()));
        }

        EXPECT_TRUE(is_empty(buffer));
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetDenseChunk,
    MockSubsetDenseChunkBlockIndexTest,
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

class MockSubsetDenseChunkIndexBlockTest : 
    public MockSubsetDenseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetDenseChunkIndexBlockTest, Dense) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsetDenseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, bounds);
        int c_first = block.first * dim.second, c_last = block.second * dim.second, c_len = c_last - c_first;

        int stride = std::min(dim.second, 100);
        std::vector<double> buffer(dim.first * stride);
        mock.template extract<true>(*r_indices, c_first, c_len, work, buffer.data(), stride);

        auto ext = ref->dense_row(c_first, c_len);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_len);
            EXPECT_EQ(refrow, slice_and_fill(buffer, r, stride, c_len));
        }

        EXPECT_TRUE(is_empty(buffer));
    }

    {
        auto c_indices = create_indices(dim.second, bounds);
        int r_first = block.first * dim.first,  r_last = block.second * dim.first,  r_len = r_last - r_first;

        int stride = std::min(dim.first, 98);
        std::vector<double> buffer(dim.second * stride);
        mock.template extract<false>(*c_indices, r_first, r_len, work, buffer.data(), stride);

        auto ext = ref->dense_column(r_first, r_len);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_len);
            EXPECT_EQ(refcol, slice_and_fill(buffer, c, stride, r_len));
        }

        EXPECT_TRUE(is_empty(buffer));
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetDenseChunk,
    MockSubsetDenseChunkIndexBlockTest,
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

class MockSubsetDenseChunkIndexIndexTest : 
    public MockSubsetDenseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, int>, std::pair<double, int> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsetDenseChunkIndexIndexTest, Dense) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto iparam1 = std::get<1>(params);
    auto iparam2 = std::get<2>(params);

    typename tatami_chunked::MockSubsetDenseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, iparam1);
        auto c_indices = create_indices(dim.second, iparam2);

        int stride = std::max(dim.second, 100);
        std::vector<double> buffer(dim.first * stride);
        mock.template extract<true>(*r_indices, *c_indices, work, buffer.data(), stride);

        auto ext = ref->dense_row(*c_indices);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(ext.get(), r, c_indices->size());
            EXPECT_EQ(refrow, slice_and_fill(buffer, r, stride, c_indices->size()));
        }

        EXPECT_TRUE(is_empty(buffer));
    }

    {
        auto c_indices = create_indices(dim.second, iparam1);
        auto r_indices = create_indices(dim.first, iparam2);

        int stride = std::max(dim.first, 98);
        std::vector<double> buffer(dim.second * stride);
        mock.template extract<false>(*c_indices, *r_indices, work, buffer.data(), stride);

        auto ext = ref->dense_column(*r_indices);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(ext.get(), c, r_indices->size());
            EXPECT_EQ(refcol, slice_and_fill(buffer, c, stride, r_indices->size()));
        }

        EXPECT_TRUE(is_empty(buffer));
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsetDenseChunk,
    MockSubsetDenseChunkIndexIndexTest,
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
