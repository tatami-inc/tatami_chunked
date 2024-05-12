#ifndef TATAMI_CHUNKED_ORACULAR_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACULAR_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file OracularSlabCache.hpp
 * @brief Create a oracle-aware cache for slabs.
 */

namespace tatami_chunked {

/**
 * @brief Oracular-aware cache for slabs.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slabs.
 * Each slab is defined as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during iteration through a `tatami::Matrix`.
 * This cache can be used for `Matrix` representations where the data is costly to load (e.g., from file) and a `tatami::Oracle` is provided to predict future accesses.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 *
 * It is assumed that each slab has the same size such that `Slab_` instances can be effectively reused between slabs without requiring reallocations.
 * For variable-sized slabs, consider using `OracularVariableSlabCache` instead.
 */
template<typename Id_, typename Index_, class Slab_> 
class OracularSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > oracle;
    size_t total;
    size_t counter = 0;

    Index_ last_slab_id = 0;
    Slab_* last_slab = NULL;

    size_t max_slabs;
    std::vector<Slab_> all_slabs;
    std::unordered_map<Id_, Slab_*> current_cache, future_cache;
    std::vector<std::pair<Id_, Slab_*> > to_populate;
    std::vector<Id_> in_need;
    size_t refresh_point = 0;

public:
    /**
     * @param ora Pointer to an `tatami::Oracle` to be used for predictions.
     * @param max_slabs Maximum number of slabs to store in the cache.
     */
    OracularSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > ora, size_t max_slabs) : 
        oracle(std::move(ora)), 
        total(oracle->total()),
        max_slabs(max_slabs) 

    {
        all_slabs.reserve(max_slabs);
        current_cache.reserve(max_slabs);
        future_cache.reserve(max_slabs);
    } 

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularSlabCache(const OracularSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularSlabCache& operator=(const OracularSlabCache&) = delete;

    /**
     * @cond
     */
    // Move operators are still okay as pointers still point to the moved vectors.
    // see https://stackoverflow.com/questions/43988553/stdvector-stdmove-and-pointer-invalidation.
    OracularSlabCache& operator=(OracularSlabCache&&) = default;
    OracularSlabCache(OracularSlabCache&&) = default;

    // Might as well define this.
    ~OracularSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method is intended to be called when `num_slabs = 0`, to provide callers with the oracle predictions for non-cached extraction of data.
     * Calls to this method should not be intermingled with calls to its overload below; the latter should only be called when `num_slabs > 0`.
     *
     * @return The next prediction from the oracle.
     */
    Index_ next() {
        return oracle->get(counter++);
    }

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     * This method should only be called if `num_slabs > 0` in the constructor; otherwise, no slabs are actually available and cannot be returned.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Cfunction_ Function to create a new slab.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the slab containing `i`, and (2) the index of row/column `i` inside that slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2 and an offset of 1.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts a `std::vector<std::pair<Id_, Slab_*> >&` specifying the slabs to be populated.
     * The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     * The second `Slab_*` element contains a pointer to a `Slab_` returned by `create()`.
     * This function should iterate over the vector and populate each slab.
     * Note that the vector is not guaranteed to be sorted. 
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == last_slab_id && last_slab) {
            return std::make_pair(last_slab, slab_info.second);
        }
        last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (counter - 1 == refresh_point) {
            // Note that, for any given populate cycle, the first prediction's
            // slab cannot already be in the cache, otherwise it would have
            // incorporated into the previous cycle. So we can skip some code.
            future_cache[slab_info.first] = NULL;
            in_need.push_back(slab_info.first);
            size_t used_slabs = 1;
            auto last_future_slab_id = slab_info.first;

            while (++refresh_point < total) {
                auto future_index = oracle->get(refresh_point);
                auto future_slab_info = identify(future_index);
                if (last_future_slab_id == future_slab_info.first) {
                    continue;
                }

                last_future_slab_id = future_slab_info.first;
                if (future_cache.find(future_slab_info.first) != future_cache.end()) {
                    continue;
                }

                if (used_slabs == max_slabs) {
                    break;
                } 
                ++used_slabs;

                auto ccIt = current_cache.find(future_slab_info.first);
                if (ccIt != current_cache.end()) {
                    auto slab_ptr = ccIt->second;
                    future_cache[future_slab_info.first] = slab_ptr;
                    current_cache.erase(ccIt);
                } else {
                    future_cache[future_slab_info.first] = NULL;
                    in_need.push_back(future_slab_info.first);
                }
            }

            auto cIt = current_cache.begin();
            for (auto a : in_need) {
                if (cIt != current_cache.end()) {
                    to_populate.emplace_back(a, cIt->second);
                    future_cache[a] = cIt->second;
                    ++cIt;
                } else {
                    // We reserved all_slabs so further push_backs() should not 
                    // trigger any reallocation or invalidation of the pointers.
                    all_slabs.push_back(create());
                    auto slab_ptr = &(all_slabs.back());
                    to_populate.emplace_back(a, slab_ptr);
                    future_cache[a] = slab_ptr;
                }
            }
            in_need.clear();

            populate(to_populate);
            to_populate.clear();

            // We always fill future_cache to the brim so every entry of
            // all_slabs should be referenced by a pointer in future_cache.
            // There shouldn't be any free cache entries remaining in
            // current_cache i.e., at this point, cIt should equal
            // current_cache.end(), as we transferred everything to
            // future_cache. Thus it is safe to clear current_cache without
            // worrying about leaking memory. The only exception is if we're at
            // the end of the predictions, in which case it doesn't matter.
            current_cache.clear();
            current_cache.swap(future_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = current_cache.find(slab_info.first);
        last_slab = ccIt->second;
        return std::make_pair(last_slab, slab_info.second);
    }

public:
    /**
     * @return Maximum number of slabs in the cache.
     */
    size_t get_max_slabs() const {
        return max_slabs;
    }

    /**
     * @return Number of slabs currently in the cache.
     */
    size_t get_num_slabs() const {
        return current_cache.size();
    }
};

}

#endif