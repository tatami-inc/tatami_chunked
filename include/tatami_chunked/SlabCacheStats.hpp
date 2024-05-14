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
     * This is used as `max_slabs` in `LruSlabCache`, `OracularSlabCache` and friends.
     */
    size_t num_slabs_in_cache;

    /**
     * @param target_length Length of the target dimension of each slab.
     * For example, if we were iterating through rows of a matrix, `target_length` would be the number of rows spanned by each slab.
     * @param non_target_length Length of the non-target dimension of each slab.
     * For example, if we were iterating through matrix rows, the `non_target_length` would be the number of columns spanned by each slab.
     * @param target_num_slabs Number of slabs required to span the full extent of the target dimension of the matrix.
     * This is used as an upper bound on the value of `num_slabs_in_cache`.
     * @param cache_size_in_elements Total size of the cache, in terms of the number of elements.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction of consecutive dimension elements, even if it exceeds `cache_size_in_cache`.
     */
    SlabCacheStats(size_t target_length, size_t non_target_length, size_t target_num_slabs, size_t cache_size_in_elements, bool require_minimum_cache) :
        slab_size_in_elements(target_length * non_target_length),
        num_slabs_in_cache(std::min(compute_num_slabs_in_cache(slab_size_in_elements, cache_size_in_elements, require_minimum_cache), target_num_slabs))
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
