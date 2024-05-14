#include <gtest/gtest.h>
#include "tatami/tatami.hpp"
#include "tatami_chunked/SlabCacheStats.hpp"

TEST(SlabCacheStats, Basic) {
    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, 1000, false);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 5);
    }

    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, 100, false);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 0);
    }

    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, 100, true);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 1);
    }

    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, -1, false);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 50);
    }
}

TEST(SlabCacheStats, ElementSize) {
    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, 20000, 10, false);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 10);
    }

    {
        tatami_chunked::SlabCacheStats stats(10, 20, 50, 20000, 0, false);
        EXPECT_EQ(stats.slab_size_in_elements, 200);
        EXPECT_EQ(stats.num_slabs_in_cache, 50);
    }
}
