#include <gtest/gtest.h>
#include "tatami_chunked/ChunkDimensionStats.hpp"

TEST(ChunkDimensionStats, Basic) {
    {
        tatami_chunked::ChunkDimensionStats<int> stuff(100, 20);
        EXPECT_EQ(stuff.dimension_extent, 100);
        EXPECT_EQ(stuff.chunk_length, 20);
        EXPECT_EQ(stuff.num_chunks, 5);
        EXPECT_EQ(stuff.last_chunk_length, 20);
        EXPECT_EQ(tatami_chunked::get_chunk_length(stuff, 0), 20);
        EXPECT_EQ(tatami_chunked::get_chunk_length(stuff, 4), 20);
    }

    {
        tatami_chunked::ChunkDimensionStats<int> stuff(55, 20);
        EXPECT_EQ(stuff.dimension_extent, 55);
        EXPECT_EQ(stuff.chunk_length, 20);
        EXPECT_EQ(stuff.num_chunks, 3);
        EXPECT_EQ(stuff.last_chunk_length, 15);
        EXPECT_EQ(tatami_chunked::get_chunk_length(stuff, 0), 20);
        EXPECT_EQ(tatami_chunked::get_chunk_length(stuff, 2), 15);
    }

    {
        tatami_chunked::ChunkDimensionStats<int> stuff(19, 40);
        EXPECT_EQ(stuff.dimension_extent, 19);
        EXPECT_EQ(stuff.chunk_length, 40);
        EXPECT_EQ(stuff.num_chunks, 1);
        EXPECT_EQ(stuff.last_chunk_length, 19);
        EXPECT_EQ(tatami_chunked::get_chunk_length(stuff, 0), 19);
    }

    {
        tatami_chunked::ChunkDimensionStats<int> stuff;
        EXPECT_EQ(stuff.dimension_extent, 0);
        EXPECT_EQ(stuff.chunk_length, 0);
        EXPECT_EQ(stuff.num_chunks, 0);
        EXPECT_EQ(stuff.last_chunk_length, 0);
    }
}
