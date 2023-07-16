#include <gtest/gtest.h>
#include "tatami_chunked/OracleSlabCache.hpp"

#include <random>

class OracleSlabCacheTest : public ::testing::Test {
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
            [&](std::vector<std::pair<unsigned char, TestSlab*> >& chunks_in_need) -> void {
                for (auto& x : chunks_in_need) {
                    x.second->chunk_id = x.first;
                    x.second->populate_number = counter++;
                }
            }
        );
    }
};

TEST_F(OracleSlabCacheTest, Consecutive) {
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

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 100, 3);
    int counter = 0;
    int nalloc = 0;

    for (int i = 1; i < 9; ++i) {
        auto out = next(cache, counter, nalloc); 
        EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(i));
        EXPECT_EQ(out.first->populate_number, i - 1);
        EXPECT_EQ(out.second, i);
    }

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

TEST_F(OracleSlabCacheTest, AllPredictions) {
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
        42, // Cycle 3
        24,
        15
    };

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 100, 3);
    int counter = 0;
    int nalloc = 0;

    auto out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 22.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 31.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 14.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 28.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 8);

    out = next(cache, counter, nalloc); // fetching 45.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->populate_number, 3);
    EXPECT_EQ(out.second, 5);

    out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 36.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 6);

    out = next(cache, counter, nalloc); // fetching 32.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 42.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->populate_number, 3);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 24.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 4);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 15.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 5);

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

TEST_F(OracleSlabCacheTest, LimitedPredictions) {
    std::vector<int> predictions{
        11, // Cycle 1
        22, 
        12, 
        31, // Cycle 2 (forced by limited predictions)
        23, 
        14, 
        45, // Cycle 3
        51, 
        32, 
        34, // Cycle 4 (forced by limited predictions)
        15,
        12, 
        23, // Cycle 5 (forced by limited predictions)
        25
    };

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedOracle<int> >(predictions.data(), predictions.size()), 3, 3);
    int counter = 0;
    int nalloc = 0;

    auto out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 22.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 31.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 14.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 45.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->populate_number, 3);
    EXPECT_EQ(out.second, 5);

    out = next(cache, counter, nalloc); // fetching 51.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(5));
    EXPECT_EQ(out.first->populate_number, 4);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 32.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 34.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 15.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 5);
    EXPECT_EQ(out.second, 5);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 5);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 6);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 25.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 6);
    EXPECT_EQ(out.second, 5);

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}
