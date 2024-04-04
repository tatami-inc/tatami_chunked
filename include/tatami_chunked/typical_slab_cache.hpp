#ifndef TATAMI_CHUNKED_TYPICAL_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_TYPICAL_SLAB_CACHE_HPP

#include <memory>
#include <type_traits>

#include "OracleSlabCache.hpp"
#include "SubsettedOracleSlabCache.hpp"
#include "LruSlabCache.hpp"

/**
 * @file typical_slab_cache.hpp
 * @brief Classes for typical caching of slabs. 
 */

namespace tatami_chunked {

/**
 * @brief Workspace for typical slab extraction.
 *
 * Implements a workspace to initialize the slab caches (i.e., the `LruSlabCache` and `OracleSlabCache`) for extraction from a chunked matrix representation.
 * This is intended to be a member of an `Extractor` class, allowing extraction to switch between different caches, e.g., when `ExtractorBase::set_oracle()` is called.
 * It also handles the calculation of various cache size statistics.
 *
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Slab_ Class that contains a single slab.
 * @tparam subset_ Whether to report the requested subset from each slab.
 * Only relevant when an oracle is available.
 */
template<bool oracle_, bool subset_, typename Index_, class Slab_>
struct TypicalSlabCacheWorkspace {
    /**
     * Default constructor.
     * Instances constructed in this manner should only be used for copy/move-assignment from instances created with the other constructor.
     */
    TypicalSlabCacheWorkspace() = default;

    /**
     * @param primary_length Length of the primary dimension of each slab.
     * The primary dimension contains the elements to be extracted from each cached slab.
     * For example, if we were iterating through rows of a matrix, `primary_length` would be the number of rows spanned by each slab.
     * @param secondary_length Length of the secondary dimension of each slab.
     * This is the dimension that is not the primary, e.g., if we were iterating through matrix rows, the `secondary_length` would be the number of columns spanned by each slab.
     * @param cache_size_in_elements Total size of the cache, in terms of the number of elements.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction of consecutive dimension elements, even if it exceeds `cache_size_in_cache`.
     * @param oracle Oracle containing the predicted accesses.
     * Only relevant if `oracle_ = true`, otherwise a bool should be passed in (which is ignored).
     */
    TypicalSlabCacheWorkspace(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache, tatami::MaybeOracle<oracle_, Index_> oracle) {
        slab_size_in_elements = static_cast<size_t>(primary_length) * static_cast<size_t>(secondary_length);

        if (!slab_size_in_elements) {
            num_slabs_in_cache = 0;
        } else {
            num_slabs_in_cache = (slab_size_in_elements ? cache_size_in_elements / slab_size_in_elements : 1);
            if (num_slabs_in_cache == 0 && require_minimum_cache) {
                num_slabs_in_cache = 1;
            }
        }

        if constexpr(!oracle_) {
            cache = LruSlabCache<Index_, Slab_>(num_slabs_in_cache);
        } else if constexpr(!subset_) {
            cache = OracleSlabCache<Index_, Index_, Slab_>(std::move(oracle), 10000, num_slabs_in_cache);
        } else {
            cache = SubsettedOracleSlabCache<Index_, Index_, Slab_>(std::move(oracle), 10000, num_slabs_in_cache);
        }
    }

public:
    /**
     * Size of each slab, in terms of the number of elements.
     */
    size_t slab_size_in_elements;

    /**
     * Number of slabs that can fit in the cache.
     */
    size_t num_slabs_in_cache;

    /**
     * Type of the slab cache, depending on `oracle_` and `subset_`.
     */
    typedef typename std::conditional<
        oracle_,
        typename std::conditional<
            subset_, 
            SubsettedOracleSlabCache<Index_, Index_, Slab_>, 
            OracleSlabCache<Index_, Index_, Slab_> 
        >::type,
        LruSlabCache<Index_, Slab_>
    >::type SlabCache;

    /**
     * Cache of to-be-used slabs, based on an `Oracle`'s predictions.
     * This may be NULL, see `set_oracle()` for more details.
     */
    SlabCache cache;
};

}

#endif
