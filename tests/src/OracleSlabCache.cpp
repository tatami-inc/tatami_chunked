#include <gtest/gtest.h>
#include "tatami_chunked/OracleSlabCache.hpp"

#include <random>

class OracleSlabCacheTestMethods {
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
            [&](std::vector<std::pair<unsigned char, TestSlab*> >& in_need) -> void {
                for (auto& x : in_need) {
                    auto& current = *(x.second);
                    current.chunk_id = x.first;
                    current.populate_number = counter++;
                }
            }
        );
    }
};

class OracleSlabCacheTest : public ::testing::Test, public OracleSlabCacheTestMethods {};

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

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_shared<tatami::FixedVectorOracle<int> >(std::move(predictions)), 3);
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

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedVectorOracle<int> >(std::move(predictions)), 3);
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

TEST_F(OracleSlabCacheTest, ShortCircuit) {
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

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedVectorOracle<int> >(std::move(predictions)), 3);
    int counter = 0;
    int nalloc = 0;

    auto out = next(cache, counter, nalloc); // fetching 11.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 32.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 33.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(3));
    EXPECT_EQ(out.first->populate_number, 1);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 21.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 1);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 10.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 0);

    out = next(cache, counter, nalloc); // fetching 12.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(1));
    EXPECT_EQ(out.first->populate_number, 0);
    EXPECT_EQ(out.second, 2);

    out = next(cache, counter, nalloc); // fetching 44.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->populate_number, 3);
    EXPECT_EQ(out.second, 4);

    out = next(cache, counter, nalloc); // fetching 46.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(4));
    EXPECT_EQ(out.first->populate_number, 3);
    EXPECT_EQ(out.second, 6);

    out = next(cache, counter, nalloc); // fetching 23.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 3);

    out = next(cache, counter, nalloc); // fetching 24.
    EXPECT_EQ(out.first->chunk_id, static_cast<unsigned char>(2));
    EXPECT_EQ(out.first->populate_number, 2);
    EXPECT_EQ(out.second, 4);

    EXPECT_EQ(nalloc, 3); // respects the max cache size.
}

class OracleSlabCacheStressTest : public ::testing::TestWithParam<int>, public OracleSlabCacheTestMethods {};

TEST_P(OracleSlabCacheStressTest, Stressed) {
    auto cache_size = GetParam();

    std::mt19937_64 rng(cache_size + 1);
    std::vector<int> predictions(10000);
    for (size_t i = 0; i < predictions.size(); ++i) {
        predictions[i] = rng() % 50 + 10;
    }

    tatami_chunked::OracleSlabCache<unsigned char, int, TestSlab> cache(std::make_unique<tatami::FixedViewOracle<int> >(predictions.data(), predictions.size()), cache_size);
    int counter = 0;
    int nalloc = 0;

    for (size_t i = 0; i < predictions.size(); ++i) {
        auto out = next(cache, counter, nalloc);
        EXPECT_EQ(out.first->chunk_id, predictions[i] / 10);
        EXPECT_EQ(out.second, predictions[i] % 10);
    }

    EXPECT_EQ(nalloc, std::min({ 5, cache_size }));
}

INSTANTIATE_TEST_SUITE_P(
    OracleSlabCache,
    OracleSlabCacheStressTest,
    ::testing::Values(3, 5, 10)  // max cache size
);
