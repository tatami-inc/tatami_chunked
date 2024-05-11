#ifndef TATAMI_CHUNKED_SLAB_CACHE_STATS_HPP
#define TATAMI_CHUNKED_SLAB_CACHE_STATS_HPP

#include <algorithm>

/**
 * @file SlabCacheStats.hpp
 * @brief Slab cache statistics.
 */

namespace tatami_chunked {

/**
 * @brief Typical statistics for slab cache construction.
 */
struct SlabCacheStats {
    /**
     * Size of each slab, in terms of the number of elements.
     */
    size_t slab_size_in_elements;

    /**
     * Number of slabs that can fit in the cache.
     * This is used as `num_slabs` in `OracularSlabCache` and `OracularSubsettedSlabCache`, and as `m` in `LruSlabCache`.
     */
    size_t num_slabs_in_cache;

    /**
     * @param primary_length Length of the primary dimension of each slab.
     * The primary dimension contains the elements to be extracted from each cached slab.
     * For example, if we were iterating through rows of a matrix, `primary_length` would be the number of rows spanned by each slab.
     * @param secondary_length Length of the secondary dimension of each slab.
     * This is the dimension that is not the primary, e.g., if we were iterating through matrix rows, the `secondary_length` would be the number of columns spanned by each slab.
     * @param cache_size_in_elements Total size of the cache, in terms of the number of elements.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction of consecutive dimension elements, even if it exceeds `cache_size_in_cache`.
     */
    SlabCacheStats(size_t primary_length, size_t secondary_length, size_t cache_size_in_elements, bool require_minimum_cache) :
        slab_size_in_elements(primary_length * secondary_length),
        num_slabs_in_cache(compute_num_slabs_in_cache(slab_size_in_elements, cache_size_in_elements, require_minimum_cache))
    {}

private:
    static size_t compute_num_slabs_in_cache(size_t slab_size_in_elements, size_t cache_size_in_elements, bool require_minimum_cache) {
        if (slab_size_in_elements == 0) {
            return 0;
        }

        auto tmp = cache_size_in_elements / slab_size_in_elements;
        if (tmp == 0 && require_minimum_cache) {
            return 1;
        } 

        return tmp;
    }
};

}

#endif
