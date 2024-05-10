#include <gtest/gtest.h>
#include "tatami_chunked/OracularVariableSlabCache.hpp"

#include <random>

class OracularVariableSlabCacheTestMethods {
protected:
    struct TestSlab {
        unsigned char chunk_id;
        int cycle;
        int populate_number;
        int reuse_number;
    };

    template<class Cache_>
    auto simple_next(Cache_& cache, int& counter, int& nalloc, int& cycle) {
        return cache.next(
            [](int i) -> std::pair<unsigned char, int> {
                return std::make_pair<unsigned char, int>(i / 10, i % 10);
            },
            [](int i) -> int {
                return i * 10;
            },
            [](int i, const TestSlab&) -> int {
                return i * 10;
            },
            [&]() -> TestSlab {
                ++nalloc;
                return TestSlab();
            },
            [&](std::vector<std::pair<unsigned char, size_t> >& in_need, std::vector<std::pair<unsigned char, size_t> >& to_reuse, std::vector<TestSlab>& all_slabs) {
                for (auto& x : in_need) {
                    auto& current = all_slabs[x.second];
                    current.chunk_id = x.first;
                    current.cycle = cycle;
                    current.populate_number = counter++;
                    current.reuse_number = 0;
                }

                for (auto& y : to_reuse) {
                    auto& current = all_slabs[y.second];
                    current.cycle = cycle;
                    ++(current.reuse_number);
                }

                ++cycle;
            }
        );
    }

    template<class Cache_>
    auto complex_next(Cache_& cache, int& counter, int& nalloc, int& cycle) {
        return cache.next(
            [](int i) -> std::pair<unsigned char, int> {
                return std::make_pair<unsigned char, int>(i / 10, i % 10);
            },
            [](int i) -> int {
                return i * 10;
            },
            [](int, const TestSlab& slab) -> int {
                return slab.chunk_id * 5; // difference in the sizes here!
            },
            [&]() -> TestSlab {
                ++nalloc;
                return TestSlab();
            },
            [&](std::vector<std::pair<unsigned char, size_t> >& in_need, std::vector<std::pair<unsigned char, size_t> >& to_reuse, std::vector<TestSlab>& all_slabs) {
                for (auto& x : in_need) {
                    auto& current = all_slabs[x.second];
                    current.chunk_id = x.first;
                    current.cycle = cycle;
                    current.populate_number = counter++;
                    current.reuse_number = 0;
                }

                for (auto& y : to_reuse) {
                    auto& current = all_slabs[y.second];
                    current.cycle = cycle;
                    ++(current.reuse_number);
                }

                ++cycle;
            }
        );
    }
};

class OracularVariableSlabCacheTest : public ::testing::Test, public OracularVariableSlabCacheTestMethods {};

TEST_F(OracularVariableSlabCacheTest, Consecutive) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        33,
        44,
        55, // Cycle 2
        66,
        77, // Cycle 3
        88, // Cycle 4 
        99  // Cycle 5 
    };

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 120);
    EXPECT_EQ(cache.get_max_size(), 120);
    EXPECT_EQ(cache.get_used_size(), 0);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 4; ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 100);
    }

    for (size_t i = 4; i < 6; ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 2);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 110);
    }

    for (size_t i = 6; i < predictions.size(); ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 3 + i - 6);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), out.first->chunk_id * 10);
    }

    // Max of 4 slabs in circulation at any given time.
    EXPECT_EQ(nalloc, 4);
}

TEST_F(OracularVariableSlabCacheTest, Reverse) {
    std::vector<int> predictions{
        66, // Cycle 1
        55, // Cycle 2
        44,
        33, // Cycle 3
        22, 
        11 
    };

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 100);
    EXPECT_EQ(cache.get_max_size(), 100);
    EXPECT_EQ(cache.get_used_size(), 0);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 1; ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 60);
    }

    for (size_t i = 1; i < 3; ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 2);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 90);
    }

    for (size_t i = 3; i < predictions.size(); ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 3);
        EXPECT_EQ(out.first->populate_number, i);
        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 60);
    }

    // Max of 3 slabs in circulation at any given time.
    EXPECT_EQ(nalloc, 3);
}

TEST_F(OracularVariableSlabCacheTest, Reused) {
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

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 80);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 7; ++i) { // fetching 11 to 28.
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
        } else if (out.first->chunk_id == 2) {
            EXPECT_EQ(out.first->populate_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
        } 

        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 60);
    }

    for (size_t i = 7; i < 12; ++i) { // fetching 45 to 42.
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 2);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->reuse_number, 1);
        } else if (out.first->chunk_id == 3) {
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->reuse_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->reuse_number, 0);
        } 

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 80);
    }

    for (size_t i = 12, end = predictions.size(); i < end; ++i) { // fetching 24 to 15.
        auto out = simple_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 3);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->reuse_number, 2);
        } else {
            EXPECT_EQ(out.first->populate_number, 4); // slab 2 needs to be reloaded.
            EXPECT_EQ(out.first->reuse_number, 0);
        }

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 30);
    }

    // Max of 3 slabs in circulation at any given time.
    EXPECT_EQ(nalloc, 3);
}

TEST_F(OracularVariableSlabCacheTest, ReusedContracted) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        12, 
        31, // Cycle 2 (reuses slabs 1 and 2 at half price)
        23, 
        14, 
        28,
        45, // Cycle 3 (reuses slab 1 at half price)
        11, 
        36, // Cycle 4 (reuses slab 4 at half price)
        32,
        42,
        24, // Cycle 5 (no slab reuse).
        15
    };

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 50);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 3; ++i) { // fetching 11 to 12.
        auto out = complex_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
        } else {
            EXPECT_EQ(out.first->populate_number, 1);
        } 

        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 30);
    }

    for (size_t i = 3; i < 7; ++i) { // fetching 31 to 28.
        auto out = complex_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 2);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->reuse_number, 1);
        } else if (out.first->chunk_id == 2) {
            EXPECT_EQ(out.first->populate_number, 1);
            EXPECT_EQ(out.first->reuse_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->reuse_number, 0);
        } 

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 45);
    }

    for (size_t i = 7; i < 9; ++i) { // fetching 45 to 11.
        auto out = complex_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 3);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
            EXPECT_EQ(out.first->reuse_number, 2);
        } else {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->reuse_number, 0);
        } 

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 45);
    }

    for (size_t i = 9; i < 12; ++i) { // fetching 36 to 42.
        auto out = complex_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 4);

        if (out.first->chunk_id == 4) {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->reuse_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 4);
            EXPECT_EQ(out.first->reuse_number, 0);
        } 

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 50);
    }

    for (size_t i = 12, end = predictions.size(); i < end; ++i) { // fetching 24 to 15.
        auto out = complex_next(cache, counter, nalloc, cycle); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 5);

        if (out.first->chunk_id == 2) {
            EXPECT_EQ(out.first->populate_number, 5); // slab 1 needs to be reloaded.
            EXPECT_EQ(out.first->reuse_number, 0);
        } else {
            EXPECT_EQ(out.first->populate_number, 6); // slab 2 needs to be reloaded.
            EXPECT_EQ(out.first->reuse_number, 0);
        }

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 30);
    }

    // Max of 3 slabs in circulation at any given time.
    EXPECT_EQ(nalloc, 3);
}


TEST_F(OracularVariableSlabCacheTest, ConsecutiveReuse) {
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

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_shared<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), 60);

    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < 8; ++i) { // fetching 11 to 12
        auto out = simple_next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 1);

        if (out.first->chunk_id == 1) {
            EXPECT_EQ(out.first->populate_number, 0);
        } else if (out.first->chunk_id == 3) { // yes, as slab 3 is created before slab 2.
            EXPECT_EQ(out.first->populate_number, 1);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
        } 

        EXPECT_EQ(out.first->reuse_number, 0);
        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 60);
    }

    for (size_t i = 8; i < predictions.size(); ++i) { // fetching 44 to 24.
        auto out = simple_next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(predictions[i] / 10));
        EXPECT_EQ(out.first->cycle, 2);

        if (out.first->chunk_id == 4) {
            EXPECT_EQ(out.first->populate_number, 3);
            EXPECT_EQ(out.first->reuse_number, 0);
        } else {
            EXPECT_EQ(out.first->populate_number, 2);
            EXPECT_EQ(out.first->reuse_number, 1);
        } 

        EXPECT_EQ(out.second, predictions[i] % 10);
        EXPECT_EQ(cache.get_used_size(), 60);
    }
}

class OracularVariableSlabCacheStressTest : public ::testing::TestWithParam<int>, public OracularVariableSlabCacheTestMethods {};

TEST_P(OracularVariableSlabCacheStressTest, Stressed) {
    auto cache_size = GetParam();

    std::mt19937_64 rng(cache_size + 1);
    std::vector<int> predictions(10000);
    for (size_t i = 0; i < predictions.size(); ++i) {
        predictions[i] = rng() % 50 + 10;
    }

    tatami_chunked::OracularVariableSlabCache<unsigned char, int, TestSlab, int> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), cache_size);
    int counter = 0;
    int nalloc = 0;
    int cycle = 1;

    for (size_t i = 0; i < predictions.size(); ++i) {
        auto out = simple_next(cache, counter, nalloc, cycle);
        EXPECT_EQ(out.first->chunk_id, predictions[i] / 10);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    EXPECT_TRUE(nalloc <= 5);
}

INSTANTIATE_TEST_SUITE_P(
    OracularVariableSlabCache,
    OracularVariableSlabCacheStressTest,
    ::testing::Values(30, 50, 100)  // max cache size
);
