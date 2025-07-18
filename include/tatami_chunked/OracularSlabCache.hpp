#ifndef TATAMI_CHUNKED_ORACULAR_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACULAR_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include <type_traits>
#include <memory>
#include <cstddef>

#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"

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
 * @tparam track_reuse_ Whether to track slabs in the cache that are re-used.
 *
 * Implement an oracle-aware cache for slabs.
 * Each slab is defined as the set of chunks required to read an element of the target dimension (or a contiguous block/indexed subset thereof) from a `tatami::Matrix`.
 * This cache can be used for `Matrix` representations where the data is costly to load (e.g., from file) and a `tatami::Oracle` is provided to predict future accesses on the target dimension.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 *
 * It is assumed that each slab has the same size such that `Slab_` instances can be effectively reused between slabs without requiring any reallocation of memory.
 * For variable-sized slabs, consider using `OracularVariableSlabCache` instead.
 */
template<typename Id_, typename Index_, class Slab_, bool track_reuse_ = false> 
class OracularSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > my_oracle;
    tatami::PredictionIndex my_total;
    tatami::PredictionIndex my_counter = 0;

    Index_ my_last_slab_id = 0;
    Slab_* my_last_slab = NULL;

    typedef std::vector<Slab_> SlabPool;
    typename SlabPool::size_type my_max_slabs;
    SlabPool my_all_slabs;

    std::unordered_map<Id_, Slab_*> my_current_cache, my_future_cache;
    std::vector<std::pair<Id_, Slab_*> > my_to_populate;
    std::vector<Id_> my_in_need;
    tatami::PredictionIndex my_refresh_point = 0;

    typename std::conditional<track_reuse_, std::vector<std::pair<Id_, Slab_*> >, bool>::type my_to_reuse;

public:
    /**
     * @tparam MaxSlabs_ Integer type of the maximum number of slabs.
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param max_slabs Maximum number of slabs to store in the cache.
     */
    template<typename MaxSlabs_>
    OracularSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > oracle, MaxSlabs_ max_slabs) : 
        my_oracle(std::move(oracle)), 
        my_total(my_oracle->total()),
        my_max_slabs(sanisizer::cast<decltype(my_max_slabs)>(max_slabs)) 
    {
        my_all_slabs.reserve(max_slabs);
        my_current_cache.reserve(max_slabs);
        my_future_cache.reserve(max_slabs);
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
        return my_oracle->get(my_counter++);
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
     * @param identify Function that accepts `i`, an `Index_` containing the predicted index of a single element on the target dimension.
     * This should return a pair containing:
     * 1. An `Id_`, the identifier of the slab containing `i`.
     *    This is typically defined as the index of the slab on the target dimension.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2.
     * 2. An `Index_`, the index of row/column `i` inside that slab.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would yield an offset of 1.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function where the arguments depend on `track_reuse_`.
     * The return value is ignored.
     * - If `track_reuse = false`, this function accepts a single `std::vector<std::pair<Id_, Slab_*> >&` specifying the slabs to be populated.
     *   The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     *   The second `Slab_*` element contains a pointer to a `Slab_` returned by `create()`.
     *   This function should iterate over the vector and populate each slab.
     *   The vector is guaranteed to be non-empty but is not guaranteed to be sorted. 
     * - If `track_reuse_ = true`, this function accepts two `std::vector<std::pair<Id_, Slab_*> >&` arguments.
     *   The first vector (`to_populate`) specifies the slabs to be populated, identical to the sole expected argument when `track_reuse_ = false`.
     *   The second vector (`to_reuse`) specifies the existing slabs in the cache to be reused.
     *   This function should iterate over `to_populate` and populate each slab. 
     *   The function may also iterate over `to_reuse` to perform some housekeeping on the existing slabs (e.g., defragmentation).
     *   `to_populate` is guaranteed to be non-empty but `to_reuse` is not.
     *   Neither vector is guaranteed to be sorted. 
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == my_last_slab_id && my_last_slab) {
            return std::make_pair(my_last_slab, slab_info.second);
        }
        my_last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (my_counter - 1 == my_refresh_point) {
            // Note that, for any given populate cycle, the first prediction's
            // slab cannot already be in the cache, otherwise it would have
            // incorporated into the previous cycle. So we can skip some code.
            my_future_cache[slab_info.first] = NULL;
            my_in_need.push_back(slab_info.first);
            decltype(my_max_slabs) used_slabs = 1;
            auto last_future_slab_id = slab_info.first;

            while (++my_refresh_point < my_total) {
                auto future_index = my_oracle->get(my_refresh_point);
                auto future_slab_info = identify(future_index);
                if (last_future_slab_id == future_slab_info.first) {
                    continue;
                }

                last_future_slab_id = future_slab_info.first;
                if (my_future_cache.find(future_slab_info.first) != my_future_cache.end()) {
                    continue;
                }

                if (used_slabs == my_max_slabs) {
                    break;
                } 
                ++used_slabs;

                auto ccIt = my_current_cache.find(future_slab_info.first);
                if (ccIt == my_current_cache.end()) {
                    my_future_cache[future_slab_info.first] = NULL;
                    my_in_need.push_back(future_slab_info.first);

                } else {
                    auto slab_ptr = ccIt->second;
                    my_future_cache[future_slab_info.first] = slab_ptr;
                    my_current_cache.erase(ccIt);
                    if constexpr(track_reuse_) {
                        my_to_reuse.emplace_back(future_slab_info.first, slab_ptr);
                    }
                }
            }

            auto cIt = my_current_cache.begin();
            for (auto a : my_in_need) {
                if (cIt != my_current_cache.end()) {
                    my_to_populate.emplace_back(a, cIt->second);
                    my_future_cache[a] = cIt->second;
                    ++cIt;
                } else {
                    // We reserved my_all_slabs so further push_backs() should not 
                    // trigger any reallocation or invalidation of the pointers.
                    my_all_slabs.push_back(create());
                    auto slab_ptr = &(my_all_slabs.back());
                    my_to_populate.emplace_back(a, slab_ptr);
                    my_future_cache[a] = slab_ptr;
                }
            }
            my_in_need.clear();

            if constexpr(track_reuse_) {
                populate(my_to_populate, my_to_reuse);
            } else {
                populate(my_to_populate);
            }

            my_to_populate.clear();
            if constexpr(track_reuse_) {
                my_to_reuse.clear();
            }

            // We always fill my_future_cache to the brim so every entry of
            // my_all_slabs should be referenced by a pointer in
            // my_future_cache.  There shouldn't be any free cache entries
            // remaining in my_current_cache i.e., at this point, cIt should
            // equal my_current_cache.end(), as we transferred everything to
            // my_future_cache. Thus it is safe to clear my_current_cache
            // without worrying about leaking memory. The only exception is if
            // we run out of predictions, in which case it doesn't matter.
            my_current_cache.clear();
            my_current_cache.swap(my_future_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = my_current_cache.find(slab_info.first);
        my_last_slab = ccIt->second;
        return std::make_pair(my_last_slab, slab_info.second);
    }

public:
    /**
     * @return Maximum number of slabs in the cache.
     * The type is an unsigned integer defined inin `std::vector::size_type`.
     */
    auto get_max_slabs() const {
        return my_max_slabs;
    }

    /**
     * @return Number of slabs currently in the cache.
     * The type is an unsigned integer defined in `std::vector::size_type`.
     */
    auto get_num_slabs() const {
        return my_current_cache.size();
    }
};

}

#endif
