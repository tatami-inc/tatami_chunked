#ifndef TATAMI_CHUNKED_CHUNK_DIMENSION_STATS_HPP
#define TATAMI_CHUNKED_CHUNK_DIMENSION_STATS_HPP

namespace tatami_chunked {

/**
 * @tparam Index_ Integer type.
 * 
 * Obtain the integer ceiling of `left/right`.
 *
 * @param left A non-negative number. 
 * @param right A non-negative number. 
 * 
 * @return The integer ceiling if `right > 0`, otherwise 0.
 */
template<typename Index_>
Index_ integer_ceil(Index_ left, Index_ right) {
    if (right) {
        return left / right + (left % right > 0); // avoids overflow.
    } else {
        return 0;
    }
}

/**
 * @brief Statistics for regular chunks along a dimension.
 *
 * This class holds useful statistics for contiguous equilength chunks along a dimension.
 * The first chunk starts at position zero. 
 * All chunks have the same length, except for (potentially) the last chunk, which may contain fewer elements than the chunk length.
 *
 * @tparam Index_ Integer type for the various dimensions.
 */
template<typename Index_>
struct ChunkDimensionStats {
    /**
     * @param dimension_extent Full extent of the dimension.
     * @param chunk_length Length of each chunk.
     */
    ChunkDimensionStats(Index_ dimension_extent, Index_ chunk_length) :
        dimension_extent(dimension_extent),
        chunk_length(chunk_length),
        num_chunks(integer_ceil(dimension_extent, chunk_length)),
        last_chunk_length(num_chunks ? (dimension_extent - (num_chunks - 1) * chunk_length) : 0)
    {}

    /**
     * Default constructor.
     */
    ChunkDimensionStats() : ChunkDimensionStats(0, 0) {}

    /**
     * Extent of the dimension.
     */
    Index_ dimension_extent;

    /**
     * Length of the chunks, except for (possibly) the last chunk.
     */
    Index_ chunk_length;

    /**
     * Number of chunks along this dimension. 
     */
    Index_ num_chunks;

    /**
     * Length of the last chunk.
     */
    Index_ last_chunk_length;

    /**
     * @param i Zero-based index of the chunk of interest.
     * @return Length of chunk `i`.
     */
    Index_ get_chunk_length(Index_ i) const {
        return (i + 1 == num_chunks ? last_chunk_length : chunk_length);
    }
};

}

#endif
