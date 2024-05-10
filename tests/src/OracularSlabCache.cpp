#include <gtest/gtest.h>
#include "tatami_chunked/OracularSlabCache.hpp"

#include <random>

class OracularSlabCacheTestMethods {
protected:
    struct TestSlab {
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
            [&](std::vector<std::pair<unsigned char, TestSlab*> >& in_need) -> void {
                for (auto& x : in_need) {
                    auto& current = *(x.second);
                    current.chunk_id = x.first;
                    current.populate_number = counter++;
                    current.cycle = cycle;
                }
                ++cycle;
            }
        );
    }
};

class OracularSlabCacheTest : public ::testing::Test, public OracularSlabCacheTestMethods {};

TEST_F(OracularSlabCacheTest, Consecutive) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        33,
        44, // Cycle 2
        55,
        66,
        77, // Cycle 3
        88,
        99
    };

    tatami_chunked::OracularSlabCache<unsigned char, int, TestSlab> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 3);
    EXPECT_EQ(cache.get_max_slabs(), 3);
    EXPECT_EQ(cache.get_num_slabs(), 0);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 3; ++i) { // 11 to 33
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    for (size_t i = 3; i < 6; ++i) { // 44 to 66
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->cycle, 2);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    for (size_t i = 6; i < predictions.size(); ++i) { // 77 to 99
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->cycle, 3);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
    EXPECT_EQ(cache.get_num_slabs(), 3);
}

TEST_F(OracularSlabCacheTest, Reuse) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        12, 
        31,
        23, 
        14, 
        28,
        45, // Cycle 2
        11, 
        36,
        32,
        42,
        24, // Cycle 3
        15
    };

    tatami_chunked::OracularSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 3);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 7; ++i) { // 11 to 28.
        auto out = next(cache, counter, nalloc, cycle); // fetching 11.
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, predictions[i] % 10);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
        } else if (out.first->chunk_id == 2) {
            EXPECT_EQ(out.first->populate_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
        }
    }

    for (size_t i = 7; i < 12; ++i) { // 45 to 42.
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.second, predictions[i] % 10);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->cycle, 1);
        } else if (out.first->chunk_id == 3) {
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->cycle, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->cycle, 2);
        }
    }

    for (size_t i = 12; i < predictions.size(); ++i) { // 24 to 15.
        auto out = next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.second, predictions[i] % 10);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->cycle, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 4); // slab 2 needs to be reloaded.
            EXPECT_EQ(out.first->cycle, 3);
        }
    }
}

TEST_F(OracularSlabCacheTest, ConsecutiveReuse) {
    // This just checks that the short-circuiting works correctly for
    // consecutive requests to the same slab.

    std::vector<int> predictions{
        11, // Cycle 1.
        12, 
        32, 
        33, 
        21,
        23,
        10, 
        12, 
        44, // Cycle 2.
        46, 
        23, 
        24
    };

    tatami_chunked::OracularSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 3);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 8; ++i) { // fetching 11 to 12
        auto out = next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.second, predictions[i] % 10);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
        } else if (out.first->chunk_id == 3) { // yes, as slab 3 is created before slab 2.
            EXPECT_EQ(out.first->populate_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
        } 
    }

    for (size_t i = 8; i < predictions.size(); ++i) { // fetching 44 to 24.
        auto out = next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.second, predictions[i] % 10);

        if (out.first->chunk_id == 4) {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->cycle, 2);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->cycle, 1);
        } 
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

class OracularSlabCacheStressTest : public ::testing::TestWithParam<int>, public OracularSlabCacheTestMethods {};

TEST_P(OracularSlabCacheStressTest, Stressed) {
    auto cache_size = GetParam();

    std::mt19937_64 rng(cache_size + 1);
    std::vector<int> predictions(10000);
    for (size_t i = 0; i < predictions.size(); ++i) {
        predictions[i] = rng() % 50 + 10;
    }

    tatami_chunked::OracularSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), cache_size);
    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < predictions.size(); ++i) {
        auto out = next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, predictions[i] / 10);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    EXPECT_EQ(nalloc, std::min({ 5, cache_size }));
}

INSTANTIATE_TEST_SUITE_P(
    OracularSlabCache,
    OracularSlabCacheStressTest,
    ::testing::Values(3, 5, 10)  // max cache size
);
