#include <gtest/gtest.h>
#include "tatami_chunked/SubsettedOracleSlabCache.hpp"

#include <random>

class SubsettedOracleSlabCacheTestMethods {
protected:
    struct TestSlab {
        unsigned char chunk_id;
        int populate_number;
    };

    template<class Cache_>
    auto next(Cache_& cache, int& counter, int& nalloc) {
        return cache.next(
            [](int i) -> std::pair<unsigned char, int> {
                return std::make_pair<unsigned char, int>(i / 10, i % 10);
            },
            [&]() -> TestSlab {
                ++nalloc;
                return TestSlab();
            },
            [&](std::vector<std::pair<unsigned char, int> >& in_need, auto& data) -> void {
                for (auto& x : in_need) {
                    auto& current = data[x.second];
                    current->contents.chunk_id = x.first;
                    current->contents.populate_number = counter++;
                }
            }
        );
    }

    template<class SubsetDetails_>
    void confirm_mapping_integrity(const SubsetDetails_& subset) {
        EXPECT_EQ(subset.mapping.size(), subset.indices.size());
        for (size_t i = 0; i < subset.indices.size(); ++i) {
            auto it = subset.mapping.find(subset.indices[i]);
            ASSERT_TRUE(it != subset.mapping.end());
            EXPECT_EQ(it->second, i);
        }
    }
};

class SubsettedOracleSlabCacheTest : public ::testing::Test, public SubsettedOracleSlabCacheTestMethods {};

TEST_F(SubsettedOracleSlabCacheTest, Consecutive) {
    std::vector<int> predictions{
        11, // Cycle 1
        12,

        22,
        21,
        20,

        33,
        34,
        32,
        35,
        
        44, // Cycle 2
        46,
        45,

        55,
        58,
        50,
        52,

        66,
        66,

        77, // Cycle 3
        88,
        99
    };

    tatami_chunked::SubsettedOracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 100, 3);
    int counter = 0;
    int nalloc = 0;

    for (int i = 0; i < 2; ++i) {
        auto out = next(cache, counter, nalloc); // extracting 11, 12.
        EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(1));
        EXPECT_EQ(out.first->contents.populate_number, 0);
        EXPECT_EQ(out.second, i + 1);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::BLOCK);
        EXPECT_EQ(out.first->subset.block_start, 1);
        EXPECT_EQ(out.first->subset.block_length, 2);
    }

    for (int i = 0; i < 3; ++i) { // extracting 22, 21, 20.
        auto out = next(cache, counter, nalloc); 
        EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(2));
        EXPECT_EQ(out.first->contents.populate_number, 1);
        EXPECT_EQ(out.second, 2 - i);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::BLOCK);
        EXPECT_EQ(out.first->subset.block_start, 0);
        EXPECT_EQ(out.first->subset.block_length, 3);
    }

    {
        std::vector<int> offsets { 3, 4, 2, 5 };
        for (int i = 0; i < 4; ++i) { // extracting 33, 34, 32, 35
            auto out = next(cache, counter, nalloc); 
            EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(3));
            EXPECT_EQ(out.first->contents.populate_number, 2);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 2);
            EXPECT_EQ(out.first->subset.block_length, 4);
        }
    }

    {
        std::vector<int> offsets { 4, 6, 5 };
        std::vector<int> sorted { 4, 5, 6 };

        for (int i = 0; i < 3; ++i) { // extracting 44, 46, 45
            auto out = next(cache, counter, nalloc); 
            EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(4));
            EXPECT_EQ(out.first->contents.populate_number, 3);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::INDEX);
            EXPECT_EQ(out.first->subset.indices, sorted);
            confirm_mapping_integrity(out.first->subset);
        }
    }

    {
        std::vector<int> offsets { 5, 8, 0, 2 };
        std::vector<int> sorted { 0, 2, 5, 8 };

        for (int i = 0; i < 4; ++i) { // extracting 55, 58, 50, 52
            auto out = next(cache, counter, nalloc); 
            EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(5));
            EXPECT_EQ(out.first->contents.populate_number, 4);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::INDEX);
            EXPECT_EQ(out.first->subset.indices, sorted);
            confirm_mapping_integrity(out.first->subset);
        }
    }

    {
        for (int i = 0; i < 2; ++i) { // extracting 66, 66
            auto out = next(cache, counter, nalloc); 
            EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(6));
            EXPECT_EQ(out.first->contents.populate_number, 5);
            EXPECT_EQ(out.second, 6);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 6);
            EXPECT_EQ(out.first->subset.block_length, 1);
        }
    }

    {
        for (int i = 0; i < 3; ++i) { // extracting 77, 88, 99
            auto out = next(cache, counter, nalloc); 
            EXPECT_EQ(out.first->contents.chunk_id, static_cast<unsigned char>(7 + i));
            EXPECT_EQ(out.first->contents.populate_number, 6 + i);
            EXPECT_EQ(out.second, 7 + i);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::SubsetSelection::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 7 + i);
            EXPECT_EQ(out.first->subset.block_length, 1);
        }
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

//class SubsettedOracleSlabCacheStressTest : public ::testing::TestWithParam<std::tuple<int, int> >, public SubsettedOracleSlabCacheTestMethods {};
//
//TEST_P(SubsettedOracleSlabCacheStressTest, Stressed) {
//    auto param = GetParam();
//    auto max_pred = std::get<0>(param);
//    auto cache_size = std::get<1>(param);
//
//    std::mt19937_64 rng(max_pred * cache_size + cache_size + 1);
//    std::vector<int> predictions(10000);
//    for (size_t i = 0; i < predictions.size(); ++i) {
//        predictions[i] = rng() % 50 + 10;
//    }
//
//    // Using limited predictions to force more cache interations.
//    tatami_chunked::SubsettedOracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedSubsettedOracle<int> >(predictions.data(), predictions.size()), max_pred, cache_size);
//    int counter = 0;
//    int nalloc = 0;
//
//    for (size_t i = 0; i < predictions.size(); ++i) {
//        auto out = next(cache, counter, nalloc);
//        EXPECT_EQ(out.first->chunk_id, predictions[i] / 10);
//        EXPECT_EQ(out.second, predictions[i] % 10);
//    }
//
//    EXPECT_EQ(nalloc, std::min({ 5, max_pred, cache_size }));
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    SubsettedOracleSlabCache,
//    SubsettedOracleSlabCacheStressTest,
//    ::testing::Combine(
//        ::testing::Values(3, 5, 10), // max predictions
//        ::testing::Values(3, 5, 10)  // max cache size
//    )
//);
