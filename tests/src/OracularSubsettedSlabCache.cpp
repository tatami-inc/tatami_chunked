#include <gtest/gtest.h>
#include "tatami_chunked/OracularSubsettedSlabCache.hpp"

#include <random>
#include <algorithm>
#include <unordered_set>
#include <vector>

class OracularSubsettedSlabCacheTestMethods {
protected:
    struct TestSlab {
        tatami_chunked::OracularSubsettedSlabCacheSelectionDetails<int> subset;
        unsigned char chunk_id;
        int populate_number;
        int cycle;
    };

    template<class Cache_>
    auto next(Cache_& cache, int& counter, int& nalloc, int& cycle) {
        return cache.next(
            [](int i) -> std::pair<unsigned char, int> {
                return std::make_pair<unsigned char, int>(i / 10, i % 10);
            },
            [&]() -> TestSlab {
                ++nalloc;
                return TestSlab();
            },
            [&](std::vector<std::tuple<unsigned char, TestSlab*, const tatami_chunked::OracularSubsettedSlabCacheSelectionDetails<int>*> >& in_need) -> void {
                EXPECT_FALSE(in_need.empty());
                std::sort(in_need.begin(), in_need.end());
                for (auto& x : in_need) {
                    auto current = std::get<1>(x);
                    current->chunk_id = std::get<0>(x);
                    current->populate_number = counter++;
                    current->subset = *std::get<2>(x);
                    current->cycle = cycle;
                }
                ++cycle;
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

class OracularSubsettedSlabCacheTest : public ::testing::Test, public OracularSubsettedSlabCacheTestMethods {};

TEST_F(OracularSubsettedSlabCacheTest, Consecutive) {
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

    tatami_chunked::OracularSubsettedSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 3);
    EXPECT_EQ(cache.get_max_slabs(), 3);
    EXPECT_EQ(cache.get_num_slabs(), 0);
    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (int i = 0; i < 2; ++i) {
        auto out = next(cache, counter, nalloc, cycle); // extracting 11, 12.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
        EXPECT_EQ(out.first->populate_number, 0);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, i + 1);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
        EXPECT_EQ(out.first->subset.block_start, 1);
        EXPECT_EQ(out.first->subset.block_length, 2);
    }

    for (int i = 0; i < 3; ++i) { // extracting 22, 21, 20.
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
        EXPECT_EQ(out.first->populate_number, 1);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 2 - i);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
        EXPECT_EQ(out.first->subset.block_start, 0);
        EXPECT_EQ(out.first->subset.block_length, 3);
    }

    {
        std::vector<int> offsets { 3, 4, 2, 5 };
        for (int i = 0; i < 4; ++i) { // extracting 33, 34, 32, 35
            auto out = next(cache, counter, nalloc, cycle); 
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->cycle, 1);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 2);
            EXPECT_EQ(out.first->subset.block_length, 4);
        }
    }

    {
        std::vector<int> offsets { 4, 6, 5 };
        std::vector<int> sorted { 4, 5, 6 };

        for (int i = 0; i < 3; ++i) { // extracting 44, 46, 45
            auto out = next(cache, counter, nalloc, cycle); 
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->cycle, 2);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::INDEX);
            EXPECT_EQ(out.first->subset.indices, sorted);
            confirm_mapping_integrity(out.first->subset);
        }
    }

    {
        std::vector<int> offsets { 5, 8, 0, 2 };
        std::vector<int> sorted { 0, 2, 5, 8 };

        for (int i = 0; i < 4; ++i) { // extracting 55, 58, 50, 52
            auto out = next(cache, counter, nalloc, cycle); 
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(5));
            EXPECT_EQ(out.first->populate_number, 4);
            EXPECT_EQ(out.first->cycle, 2);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::INDEX);
            EXPECT_EQ(out.first->subset.indices, sorted);
            confirm_mapping_integrity(out.first->subset);
        }
    }

    {
        for (int i = 0; i < 2; ++i) { // extracting 66, 66
            auto out = next(cache, counter, nalloc, cycle); 
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(6));
            EXPECT_EQ(out.first->populate_number, 5);
            EXPECT_EQ(out.first->cycle, 2);
            EXPECT_EQ(out.second, 6);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 6);
            EXPECT_EQ(out.first->subset.block_length, 1);
        }
    }

    {
        for (int i = 0; i < 3; ++i) { // extracting 77, 88, 99
            auto out = next(cache, counter, nalloc, cycle); 
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(7 + i));
            EXPECT_EQ(out.first->populate_number, 6 + i);
            EXPECT_EQ(out.first->cycle, 3);
            EXPECT_EQ(out.second, 7 + i);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 7 + i);
            EXPECT_EQ(out.first->subset.block_length, 1);
        }
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
    EXPECT_EQ(cache.get_num_slabs(), 3);
}

TEST_F(OracularSubsettedSlabCacheTest, FullFallback) {
    std::vector<int> predictions{
        11, // Cycle 1
        22,
        33,
        12,
        24,
        36,

        55, // Cycle 2
        13, 
        45,

        28, // Cycle 3
        29,
        27,
        40, 
        56,
        45
    };

    tatami_chunked::OracularSubsettedSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 3);
    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 11.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
        EXPECT_EQ(out.first->populate_number, 0);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 1);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // because it's used in the next cycle.
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 22.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
        EXPECT_EQ(out.first->populate_number, 1);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 2);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::INDEX);
        std::vector<int> expected { 2, 4 };
        EXPECT_EQ(out.first->subset.indices, expected);
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 33.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
        EXPECT_EQ(out.first->populate_number, 2);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 3);

        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::INDEX);
        std::vector<int> expected { 3, 6 };
        EXPECT_EQ(out.first->subset.indices, expected);
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 12.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
        EXPECT_EQ(out.first->populate_number, 0);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 2);
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 24.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
        EXPECT_EQ(out.first->populate_number, 1);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 4);
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 36.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
        EXPECT_EQ(out.first->populate_number, 2);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, 6);
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 55.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(5));
        EXPECT_EQ(out.first->populate_number, 4); // even though 55 occurs before 45, it sorts afterward in our populate() function, so it is populated later.
        EXPECT_EQ(out.first->cycle, 2);
        EXPECT_EQ(out.second, 5);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // used in the next cycle.
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 13.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
        EXPECT_EQ(out.first->populate_number, 0);
        EXPECT_EQ(out.first->cycle, 1); // reused from the previous cycle.
        EXPECT_EQ(out.second, 3);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // not used in the next cycle, but it's already there, so we don't reallocate.
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 45.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
        EXPECT_EQ(out.first->populate_number, 3); // see above for 55's explanation.
        EXPECT_EQ(out.first->cycle, 2);
        EXPECT_EQ(out.second, 5);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // used in the next cycle.
    }

    {
        std::vector<int> offsets { 8, 9, 7 };
        for (int i = 0; i < 3; ++i) { // extracting 28, 29, 27
            auto out = next(cache, counter, nalloc, cycle);
            EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
            EXPECT_EQ(out.first->populate_number, 5);
            EXPECT_EQ(out.first->cycle, 3);
            EXPECT_EQ(out.second, offsets[i]);

            EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK);
            EXPECT_EQ(out.first->subset.block_start, 7);
            EXPECT_EQ(out.first->subset.block_length, 3);
        }
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 40.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
        EXPECT_EQ(out.first->populate_number, 3);
        EXPECT_EQ(out.first->cycle, 2); // reused from the previous cycle.
        EXPECT_EQ(out.second, 0);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // not used in the next cycle, but it's already extracted.
    }

    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 56.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(5));
        EXPECT_EQ(out.first->populate_number, 4);
        EXPECT_EQ(out.first->cycle, 2); // reused from the previous cycle
        EXPECT_EQ(out.second, 6);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // not used in the next cycle, but it's already extracted.
    }


    {
        auto out = next(cache, counter, nalloc, cycle); // extracting 40.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
        EXPECT_EQ(out.first->populate_number, 3);
        EXPECT_EQ(out.first->cycle, 2); // reused from the previous cycle
        EXPECT_EQ(out.second, 5);
        EXPECT_EQ(out.first->subset.selection, tatami_chunked::OracularSubsettedSlabCacheSelectionType::FULL); // not used in the next cycle, but it's already extracted.
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

class OracularSubsettedSlabCacheStressTest : public ::testing::TestWithParam<int>, public OracularSubsettedSlabCacheTestMethods {};

TEST_P(OracularSubsettedSlabCacheStressTest, Stressed) {
    auto cache_size = GetParam();
    int expected_max_cache = std::min(5, cache_size);

    std::mt19937_64 rng(cache_size + 1);
    std::vector<int> predictions(10000);
    for (size_t i = 0; i < predictions.size(); ++i) {
        predictions[i] = rng() % 50 + 10;
    }

    // Using limited predictions to force more cache interations.
    tatami_chunked::OracularSubsettedSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), cache_size);
    int counter = 0;
    int nalloc = 0;
    int cycle = 1;
    std::unordered_set<int> previous_chunks;

    for (size_t i = 0; i < predictions.size(); ++i) {
        int last_cycle = cycle;

        auto out = next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, predictions[i] / 10);
        EXPECT_EQ(out.second, predictions[i] % 10);

        // Checking the subset.
        const auto& sub = out.first->subset;
        if (sub.selection == tatami_chunked::OracularSubsettedSlabCacheSelectionType::BLOCK) {
            EXPECT_TRUE(out.second >= sub.block_start);
            EXPECT_TRUE(out.second < sub.block_start + sub.block_length);
        } else if (sub.selection == tatami_chunked::OracularSubsettedSlabCacheSelectionType::INDEX) {
            confirm_mapping_integrity(sub);
            EXPECT_TRUE(sub.mapping.find(out.second) != sub.mapping.end());
        } 

        auto cid = out.first->chunk_id;
        if (cycle != last_cycle && last_cycle != 1) {
            EXPECT_FALSE(previous_chunks.find(cid) != previous_chunks.end());
            EXPECT_EQ(expected_max_cache, previous_chunks.size());
            previous_chunks.clear();
        }
        previous_chunks.insert(cid);
        EXPECT_LE(previous_chunks.size(), cache_size);
    }

    EXPECT_EQ(nalloc, expected_max_cache);
}

INSTANTIATE_TEST_SUITE_P(
    OracularSubsettedSlabCache,
    OracularSubsettedSlabCacheStressTest,
    ::testing::Values(3, 5, 10)  // max cache size
);
