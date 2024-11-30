#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_test/tatami_test.hpp"
#include "tatami_chunked/mock_sparse_chunk.hpp"
#include "mock_blob.h"

static tatami::VectorPtr<int> create_indices(int extent, const std::pair<double, double>& config) {
    return tatami_test::create_indexed_subset(
        extent,
        config.first,
        config.second,
        extent + 997 * config.first + 2003 * config.second
    );
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
        auto full = tatami_test::simulate_vector<double>(dim.first * dim.second, [&]{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = dim.first * dim.second + 17;
            return opt;
        }());

        tatami::DenseRowMatrix<double, int> tmp(dim.first, dim.second, std::move(full));

        auto compressed = tatami::retrieve_compressed_sparse_contents<double, int>(&tmp, true, true);
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

        auto compressed2 = tatami::retrieve_compressed_sparse_contents<double, int>(ref.get(), false, true);
        dcm_chunk = tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >(MockSparseBlob<false>(
            dim.first, 
            dim.second, 
            std::move(compressed2.value),
            std::move(compressed2.index),
            std::move(compressed2.pointers)
        ));
    }
};

struct ExtractSpace {
    ExtractSpace(size_t primary, size_t secondary) : value_pool(primary * secondary), index_pool(primary * secondary) {
        number.resize(primary);
        values.reserve(primary);
        indices.reserve(primary);
        auto vptr = value_pool.data();
        auto iptr = index_pool.data();
        for (size_t p = 0; p < primary; ++p, vptr += secondary, iptr += secondary) {
            values.push_back(vptr);
            indices.push_back(iptr);
        }
    }

private:
    std::vector<double> value_pool;
    std::vector<int> index_pool;

public:
    std::vector<int> number;
    std::vector<double*> values;
    std::vector<int*> indices;

public:
    std::vector<double> values_as_vector(size_t p, size_t skip = 0) const {
        return std::vector<double>(values[p] + skip, values[p] + number[p]);
    }

    std::vector<int> indices_as_vector(size_t p, int shift = 0, size_t skip = 0) const {
        std::vector<int> output(indices[p] + skip, indices[p] + number[p]);
        for (auto& o : output) {
            o -= shift;
        }
        return output;
    }
};

/**********************************************************
 **********************************************************/

class MockSparseChunkGeneralTest : 
    public MockSimpleSparseChunkUtils,
    public ::testing::Test {
protected:
    void SetUp() {
        assemble({ 20, 20 });
    }
};

TEST_F(MockSparseChunkGeneralTest, AppendBlock) {
    auto dim = last_params;

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    // Injecting some initial 'number', to make sure the chunks append properly.
    ExtractSpace mock_data(dim.first, dim.second + 10);
    std::fill(mock_data.number.begin(), mock_data.number.end(), 10);
    mock.extract(true, 0, dim.second, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

    ExtractSpace drm_data(dim.first, dim.second + 20);
    std::fill(drm_data.number.begin(), drm_data.number.end(), 20);
    drm_chunk.extract(true, 0, dim.second, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

    ExtractSpace dcm_data(dim.first, dim.second + 5);
    std::fill(dcm_data.number.begin(), dcm_data.number.end(), 5);
    dcm_chunk.extract(true, 0, dim.second, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

    auto ext = this->ref->sparse_row();
    for (int r = 0; r < dim.first; ++r) {
        auto refrow = tatami_test::fetch(*ext, r, dim.second);
        EXPECT_EQ(refrow.value, mock_data.values_as_vector(r, 10));
        EXPECT_EQ(refrow.value, drm_data.values_as_vector(r, 20));
        EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r, 5));
        EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 0, 10));
        EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r, 0, 20));
        EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r, 0, 5));
    }
}

TEST_F(MockSparseChunkGeneralTest, AppendIndex) {
    auto dim = last_params;

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    auto c_indices = create_indices(dim.second, { 0.1, 3 });

    // Injecting some initial 'number', to make sure the chunks append properly.
    ExtractSpace mock_data(dim.first, dim.second + 10);
    std::fill(mock_data.number.begin(), mock_data.number.end(), 10);
    mock.extract(true, *c_indices, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

    ExtractSpace drm_data(dim.first, dim.second + 20);
    std::fill(drm_data.number.begin(), drm_data.number.end(), 20);
    drm_chunk.extract(true, *c_indices, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

    ExtractSpace dcm_data(dim.first, dim.second + 5);
    std::fill(dcm_data.number.begin(), dcm_data.number.end(), 5);
    dcm_chunk.extract(true, *c_indices, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

    auto ext = this->ref->sparse_row(std::move(c_indices));
    for (int r = 0; r < dim.first; ++r) {
        auto refrow = tatami_test::fetch(*ext, r, dim.second);
        EXPECT_EQ(refrow.value, mock_data.values_as_vector(r, 10));
        EXPECT_EQ(refrow.value, drm_data.values_as_vector(r, 20));
        EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r, 5));
        EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 0, 10));
        EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r, 0, 20));
        EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r, 0, 5));
    }
}

TEST_F(MockSparseChunkGeneralTest, SkippingBlock) {
    auto dim = last_params;

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    // Skipping values.
    {
        ExtractSpace mock_data(dim.first, dim.second);
        mock_data.values.clear();
        mock.extract(true, 0, dim.second, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_data.values.clear();
        drm_chunk.extract(true, 0, dim.second, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_data.values.clear();
        dcm_chunk.extract(true, 0, dim.second, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

        auto ext = this->ref->sparse_row();
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, dim.second);
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r));
            EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r));
            EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r));
        }
    }

    // Skipping indices.
    {
        ExtractSpace mock_data(dim.first, dim.second);
        mock_data.indices.clear();
        mock.extract(true, 0, dim.second, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_data.indices.clear();
        drm_chunk.extract(true, 0, dim.second, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_data.indices.clear();
        dcm_chunk.extract(true, 0, dim.second, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

        auto ext = this->ref->sparse_row();
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, dim.second);
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, drm_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r));
        }
    }
}

TEST_F(MockSparseChunkGeneralTest, SkippingIndices) {
    auto dim = last_params;

    typename tatami_chunked::MockSimpleSparseChunk::Workspace mock_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<true> >::Workspace drm_work;
    typename tatami_chunked::SimpleSparseChunkWrapper<MockSparseBlob<false> >::Workspace dcm_work;

    auto c_indices = create_indices(dim.second, { 0.1, 3 });

    // Skipping values.
    {
        ExtractSpace mock_data(dim.first, dim.second);
        mock_data.values.clear();
        mock.extract(true, *c_indices, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_data.values.clear();
        drm_chunk.extract(true, *c_indices, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_data.values.clear();
        dcm_chunk.extract(true, *c_indices, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

        auto ext = this->ref->sparse_row(c_indices);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, dim.second);
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r));
            EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r));
            EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r));
        }
    }

    // Skipping indices.
    {
        ExtractSpace mock_data(dim.first, dim.second);
        mock_data.indices.clear();
        mock.extract(true, *c_indices, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 0);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_data.indices.clear();
        drm_chunk.extract(true, *c_indices, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 0);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_data.indices.clear();
        dcm_chunk.extract(true, *c_indices, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 0);

        auto ext = this->ref->sparse_row(c_indices);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, dim.second);
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, drm_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r));
        }
    }
}

/**********************************************************
 **********************************************************/

class MockSimpleSparseChunkBlockTest : 
    public MockSimpleSparseChunkUtils,
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double> > > {
protected:
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
    int c_first = block.first * dim.second, c_len = block.second * dim.second;
    {
        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, c_first, c_len, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 7);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_chunk.extract(true, c_first, c_len, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 3);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_chunk.extract(true, c_first, c_len, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 5);

        auto ext = this->ref->sparse_row(c_first, c_len);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, c_len);
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, drm_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 7));
            EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r, 3));
            EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r, 5));
        }
    }

    // Then extracting column-wise blocks.
    int r_first = block.first * dim.first, r_len = block.second * dim.first;
    {
        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, r_first, r_len, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 1);

        ExtractSpace drm_data(dim.second, dim.first);
        drm_chunk.extract(false, r_first, r_len, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 6);

        ExtractSpace dcm_data(dim.second, dim.first);
        dcm_chunk.extract(false, r_first, r_len, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 10);

        auto ext = this->ref->sparse_column(r_first, r_len);
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(*ext, c, r_len);
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.value, drm_data.values_as_vector(c));
            EXPECT_EQ(refcol.value, dcm_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 1));
            EXPECT_EQ(refcol.index, drm_data.indices_as_vector(c, 6));
            EXPECT_EQ(refcol.index, dcm_data.indices_as_vector(c, 10));
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
            std::make_pair(0.0, 0.6),
            std::make_pair(0.25, 0.55),
            std::make_pair(0.33, 0.44)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSimpleSparseChunkIndexTest : 
    public MockSimpleSparseChunkUtils,
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double> > > {
protected:
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

        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, *c_indices, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 17);

        ExtractSpace drm_data(dim.first, dim.second);
        drm_chunk.extract(true, *c_indices, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 72);

        ExtractSpace dcm_data(dim.first, dim.second);
        dcm_chunk.extract(true, *c_indices, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 5);

        auto ext = this->ref->sparse_row(c_indices);
        for (int r = 0; r < dim.first; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, drm_data.values_as_vector(r));
            EXPECT_EQ(refrow.value, dcm_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 17));
            EXPECT_EQ(refrow.index, drm_data.indices_as_vector(r, 72));
            EXPECT_EQ(refrow.index, dcm_data.indices_as_vector(r, 5));
        }
    }

    // Then extracting column-wise blocks.
    {
        auto r_indices = create_indices(dim.first, index);

        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, *r_indices, mock_work, mock_data.values, mock_data.indices, mock_data.number.data(), 17);

        ExtractSpace drm_data(dim.second, dim.first);
        drm_chunk.extract(false, *r_indices, drm_work, drm_data.values, drm_data.indices, drm_data.number.data(), 72);

        ExtractSpace dcm_data(dim.second, dim.first);
        dcm_chunk.extract(false, *r_indices, dcm_work, dcm_data.values, dcm_data.indices, dcm_data.number.data(), 5);

        auto ext = this->ref->sparse_column(r_indices);
        for (int c = 0; c < dim.second; ++c) {
            auto refcol = tatami_test::fetch(*ext, c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.value, drm_data.values_as_vector(c));
            EXPECT_EQ(refcol.value, dcm_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 17));
            EXPECT_EQ(refcol.index, drm_data.indices_as_vector(c, 72));
            EXPECT_EQ(refcol.index, dcm_data.indices_as_vector(c, 5));
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

class MockSubsettedSparseChunkUtils {
protected:
    static inline std::unique_ptr<tatami::Matrix<double, int> > ref;

    static inline tatami_chunked::MockSubsettedSparseChunk mock;

    static inline std::pair<int, int> last_params;

    static void assemble(const std::pair<int, int>& params) {
        if (ref && last_params == params) {
            return;
        }
        last_params = params;

        const auto& dim = params;
        auto full = tatami_test::simulate_vector<double>(dim.first * dim.second, [&]{
            tatami_test::SimulateVectorOptions opt;
            opt.density = 0.2;
            opt.lower = -10;
            opt.upper = 10;
            opt.seed = dim.first * dim.second + 1333;
            return opt;
        }());

        tatami::DenseRowMatrix<double, int> tmp(dim.first, dim.second, std::move(full));

        auto compressed = tatami::retrieve_compressed_sparse_contents<double, int>(&tmp, true, true);
        mock = tatami_chunked::MockSubsettedSparseChunk(
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

class MockSubsettedSparseChunkBlockBlockTest : 
    public MockSubsettedSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsettedSparseChunkBlockBlockTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsettedSparseChunk::Workspace work;

    // Row extraction.
    {
        int r_first = bounds.first * dim.first, r_len = bounds.second * dim.first;
        int c_first = bounds.first * dim.second, c_len = bounds.second * dim.second;

        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, r_first, r_len, c_first, c_len, work, mock_data.values, mock_data.indices, mock_data.number.data(), 10);

        auto ext = ref->sparse_row(c_first, c_len);
        for (int r = r_first, r_last = r_first + r_len; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, c_len);
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 10));
        }
    }

    // Column extraction.
    {
        int r_first = block.first * dim.first, r_len = block.second * dim.first;
        int c_first = block.first * dim.second, c_len = block.second * dim.second;

        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, c_first, c_len, r_first, r_len, work, mock_data.values, mock_data.indices, mock_data.number.data(), 22);

        auto ext = ref->sparse_column(r_first, r_len);
        for (int c = c_first, c_last = c_first + c_len; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(*ext, c, r_len);
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 22));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsettedSparseChunk,
    MockSubsettedSparseChunkBlockBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 0.7),
            std::make_pair(0.2, 0.6),
            std::make_pair(0.4, 0.5)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 0.51),
            std::make_pair(0.25, 0.52),
            std::make_pair(0.33, 0.53)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsettedSparseChunkBlockIndexTest : 
    public MockSubsettedSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsettedSparseChunkBlockIndexTest, Basic) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto iparam = std::get<2>(params);

    typename tatami_chunked::MockSubsettedSparseChunk::Workspace work;

    {
        int r_first = bounds.first * dim.first, r_len = bounds.second * dim.first;
        auto c_indices = create_indices(dim.second, iparam);

        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, r_first, r_len, *c_indices, work, mock_data.values, mock_data.indices, mock_data.number.data(), 10);

        auto ext = ref->sparse_row(c_indices);
        for (int r = r_first, r_last = r_first + r_len; r < r_last; ++r) {
            auto refrow = tatami_test::fetch(*ext, r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 10));
        }
    }

    {
        int c_first = bounds.first * dim.second, c_len = bounds.second * dim.second;
        auto r_indices = create_indices(dim.first, iparam);

        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, c_first, c_len, *r_indices, work, mock_data.values, mock_data.indices, mock_data.number.data(), 5);

        auto ext = ref->sparse_column(r_indices);
        for (int c = c_first, c_last = c_first + c_len; c < c_last; ++c) {
            auto refcol = tatami_test::fetch(*ext, c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 5));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsettedSparseChunk,
    MockSubsettedSparseChunkBlockIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary range
            std::make_pair(0.0, 0.6),
            std::make_pair(0.2, 0.6),
            std::make_pair(0.4, 0.6)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 0.15),
            std::make_pair(0.25, 0.2),
            std::make_pair(0.33, 0.25)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsettedSparseChunkIndexBlockTest : 
    public MockSubsettedSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsettedSparseChunkIndexBlockTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto bounds = std::get<1>(params);
    auto block = std::get<2>(params);

    typename tatami_chunked::MockSubsettedSparseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, bounds);
        int c_first = block.first * dim.second, c_len = block.second * dim.second;

        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, *r_indices, c_first, c_len, work, mock_data.values, mock_data.indices, mock_data.number.data(), 32);

        auto ext = ref->sparse_row(c_first, c_len);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(*ext, r, c_len);
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 32));
        }
    }

    {
        auto c_indices = create_indices(dim.second, bounds);
        int r_first = block.first * dim.first, r_len = block.second * dim.first;

        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, *c_indices, r_first, r_len, work, mock_data.values, mock_data.indices, mock_data.number.data(), 99);

        auto ext = ref->sparse_column(r_first, r_len);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(*ext, c, r_len);
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 99));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsettedSparseChunk,
    MockSubsettedSparseChunkIndexBlockTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary indices 
            std::make_pair(0.0, 0.2),
            std::make_pair(0.2, 0.2),
            std::make_pair(0.4, 0.2)
        ),
        ::testing::Values( // secondary block
            std::make_pair(0.0, 0.5),
            std::make_pair(0.25, 0.45),
            std::make_pair(0.33, 0.65)
        )
    )
);

/**********************************************************
 **********************************************************/

class MockSubsettedSparseChunkIndexIndexTest : 
    public MockSubsettedSparseChunkUtils, 
    public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<double, double>, std::pair<double, double> > > {
protected:
    void SetUp() {
        assemble(std::get<0>(GetParam()));
    }
};

TEST_P(MockSubsettedSparseChunkIndexIndexTest, Sparse) {
    auto params = GetParam();
    auto dim = std::get<0>(params);
    auto iparam1 = std::get<1>(params);
    auto iparam2 = std::get<2>(params);

    typename tatami_chunked::MockSubsettedSparseChunk::Workspace work;

    {
        auto r_indices = create_indices(dim.first, iparam1);
        auto c_indices = create_indices(dim.second, iparam2);

        ExtractSpace mock_data(dim.first, dim.second);
        mock.extract(true, *r_indices, *c_indices, work, mock_data.values, mock_data.indices, mock_data.number.data(), 9);

        auto ext = ref->sparse_row(*c_indices);
        for (auto r : *r_indices) {
            auto refrow = tatami_test::fetch(*ext, r, c_indices->size());
            EXPECT_EQ(refrow.value, mock_data.values_as_vector(r));
            EXPECT_EQ(refrow.index, mock_data.indices_as_vector(r, 9));
        }
    }

    {
        auto c_indices = create_indices(dim.second, iparam1);
        auto r_indices = create_indices(dim.first, iparam2);

        ExtractSpace mock_data(dim.second, dim.first);
        mock.extract(false, *c_indices, *r_indices, work, mock_data.values, mock_data.indices, mock_data.number.data(), 19);

        auto ext = ref->sparse_column(*r_indices);
        for (auto c : *c_indices) {
            auto refcol = tatami_test::fetch(*ext, c, r_indices->size());
            EXPECT_EQ(refcol.value, mock_data.values_as_vector(c));
            EXPECT_EQ(refcol.index, mock_data.indices_as_vector(c, 19));
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    MockSubsettedSparseChunk,
    MockSubsettedSparseChunkIndexIndexTest,
    ::testing::Combine(
        ::testing::Values( // chunk dimensions
            std::make_pair(1, 50),
            std::make_pair(50, 1),
            std::make_pair(13, 42),
            std::make_pair(27, 23) 
        ),
        ::testing::Values( // primary indices
            std::make_pair(0.0, 0.2),
            std::make_pair(0.3, 0.15),
            std::make_pair(0.5, 0.35)
        ),
        ::testing::Values( // secondary indices
            std::make_pair(0.0, 0.1),
            std::make_pair(0.25, 0.3),
            std::make_pair(0.33, 0.2)
        )
    )
);
