#ifndef TATAMI_CHUNKED_ORACLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include "tatami/tatami.hpp"

/**
 * @file OracleSlabCache.hpp
 * @brief Create a oracle-aware cache for slabs.
 */

namespace tatami_chunked {

/**
 * @brief Oracle-aware cache for slabs.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slabs.
 * Each slab is defined as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during iteration through a `tatami::Matrix`.
 * This cache can be used for `Matrix` representations where the data is costly to load (e.g., from file) and a `tatami::Oracle` is provided to predict future accesses.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 */
template<typename Id_, typename Index_, class Slab_> 
class OracleSlabCache {
    tatami::OracleStream<Index_> prediction_stream;
    std::vector<std::pair<Id_, Index_> > predictions_made;
    size_t predictions_fulfilled = 0;
    size_t max_predictions;

    size_t max_slabs;
    std::unordered_map<Id_, Index_> cache_exists, next_cache_exists;
    std::vector<Slab_> cache_data, next_cache_data;
    std::vector<std::pair<Id_, Index_> > slabs_in_need; 

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_slabs Maximum number of slabs to store.
     */
    OracleSlabCache(std::unique_ptr<tatami::Oracle<Index_> > oracle, size_t per_iteration, size_t num_slabs) :
        prediction_stream(std::move(oracle)), 
        max_predictions(per_iteration),
        max_slabs(num_slabs), 
        cache_data(max_slabs),
        next_cache_data(max_slabs)
    {
        slabs_in_need.reserve(max_slabs);
    } 

    /**
     * @cond
     */
    // For testing only.
    OracleSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Sfunction_ Function to swap two slabs' contents.
     * @tparam Rfunction_ Function to determine whether a slab object is a non-mock instance.
     * @tparam Afunction_ Function to allocate memory to a mock slab, to turn it into a non-mock instance.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the slab containing `i`, and (2) the index of row/column `i` inside that slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2 and an offset of 1.
     * @param swap Function that accepts two `Slab_&` and swaps their contents.
     * @param ready Function that accepts a `const Slab_&` and returns a boolean indicating whether it has already been allocated.
     * This should return `true` for objects that have been used in `allocate()`, and `false` otherwise.
     * @param allocate Function that accepts a single default-initialized `Slab_` object,
     * and allocates sufficient memory to it in order to hold a slab's contents when used in `populate()`.
     * @param populate Function that accepts two arguments, `slabs_in_need` and `slab_data`.
     * (1) `slabs_in_need` is a `const std::vector<std::pair<Id_, Index_> >&` specifying the slabs to be populated.
     * The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     * The second `Index_` element specifies the index of the `Slab_` in `slab_data` in which to store the contents of each slab.
     * (2) `slab_data` is a `std::vector<Slab_>&` containing the cached slab contents to be populated.
     * This function should iterate over the `slabs_in_need` and populate the corresponding entries in `slab_data`.
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Sfunction_, class Rfunction_, class Afunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Sfunction_ swap, Rfunction_ ready, Afunction_ allocate, Pfunction_ populate) {
        if (predictions_made.size() > predictions_fulfilled) {
            const auto& chosen = predictions_made[predictions_fulfilled++];
            return std::make_pair(cache_data.data() + chosen.first, chosen.second);
        }

        next_cache_exists.clear();
        slabs_in_need.clear();
        size_t used = 0;

        predictions_made.clear();
        predictions_made.reserve(max_predictions);
        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!prediction_stream.next(current)) {
                break;
            }

            auto slab_id = identify(current);
            auto curslab = slab_id.first;
            auto curindex = slab_id.second;

            auto it = next_cache_exists.find(curslab);
            if (it == next_cache_exists.end()) {
                if (used == max_slabs) {
                    prediction_stream.back();
                    break;
                }

                next_cache_exists[curslab] = used;
                predictions_made.emplace_back(used, curindex);

                auto it2 = cache_exists.find(curslab);
                if (it2 != cache_exists.end()) {
                    swap(next_cache_data[used], cache_data[it2->second]);
                } else {
                    slabs_in_need.emplace_back(curslab, used);
                }

                ++used;
            } else {
                predictions_made.emplace_back(it->second, curindex);
            }
        }

        /**
         * Doing a linear scan across slabs to find the allocated but unused cache
         * elements. This is the simplest and safest approach; trying to keep
         * track of the unused caches would require a scan to fill a map/set
         * anyway, and deleting an iterator from cache_exists only works if
         * cache_exists actually contains all cache elements, which it might
         * not be if we reached max_predictions without filling up the cache.
         *
         * In any case, a linear scan should be pretty good for consecutive
         * access; the first cache elements would be the oldest, so there
         * wouldn't be any wasted iterations to find available cache elements.
         */
        size_t search = 0;
        for (const auto& c : slabs_in_need) {
            if (!ready(next_cache_data[c.second])) {
                while (search < cache_data.size() && !ready(cache_data[search])) {
                    ++search;
                }
                if (search < cache_data.size()) {
                    swap(next_cache_data[c.second], cache_data[search]);
                    ++search;
                } else {
                    /*
                     * This should be called no more than 'max_slabs' times
                     * across the lifetime of this Cache object. At any given
                     * point in time, allocated slabs will be interspersed
                     * between 'cache_data' and 'next_cache_data', so either
                     * 'next_cache_data[c.second]' is already ready, or we'll
                     * find a ready slab from the linear scan through
                     * 'cache_data'.
                     */
                    allocate(next_cache_data[c.second]);
                }
            }
        }

        populate(slabs_in_need, next_cache_data);

        cache_data.swap(next_cache_data);
        cache_exists.swap(next_cache_exists);
        predictions_fulfilled = 1; // well, because we just used one.
        const auto& chosen = predictions_made.front(); // assuming at least one prediction was made, otherwise, why was this function even called?
        return std::make_pair(cache_data.data() + chosen.first, chosen.second);
    }
};

}

#endif
