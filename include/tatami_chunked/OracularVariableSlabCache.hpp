#ifndef TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file OracularVariableSlabCache.hpp
 * @brief Oracle-aware cache for variable-size slabs.
 */

namespace tatami_chunked {

/**
 * @brief Oracle-aware cache for variable-size slabs.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Size_ Numeric type for the slab size.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for variable-size slabs.
 * This is similar to `OracleSlabCache` but enables improved cache utilization when the slabs vary in size.
 * For example, the number of non-zero entries in a sparse matrix might vary between slabs,
 * so the cache could be optimized to fit more slabs into memory when they have fewer non-zeros.
 *
 * The size of each slab is defined by `Size_`, which can be any non-negative measure of slab size.
 * This could be the number of non-zero elements, or the number of dimension elements, or the size of the slab in bytes, etc.,
 * as long as its interpretation is consistent between slabs and with the `max_size` used in the constructor.
 * Users can also differentiate between the estimated and actual size of the slab, if the latter is not known until after the slab has been loaded into memory.
 *
 * When implementing `Slab_`,  we generally suggest using a common memory pool that is referenced by each `Slab_` instance.
 * This guarantees that the actual cache size does not exceed the limit associated with `max_size` when `Slab_` instances are re-used for different slabs.
 * (Otherwise, if each `Slab_` allocates its own memory, re-use of an instance may cause its allocation to inflate to the size of the largest encountered slab.)
 * A side effect of this implementation is that callers may need to occasionally defragment the pool to ensure that enough memory is available for loading new slabs.
 */
template<typename Id_, typename Index_, class Slab_, typename Size_> 
class OracularVariableSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > oracle;
    size_t total;
    size_t counter = 0;

    Index_ last_slab_id = 0;
    Slab_* last_slab = NULL;

    size_t max_size;
    std::vector<Slab_> all_slabs;
    std::unordered_map<Id_, Slab_*> current_cache, future_cache;
    std::vector<std::pair<Id_, Slab_*> > to_populate;
    std::vector<Id_> to_reuse;
    std::vector<Id_> in_need;
    std::vector<Slab_*> free_pool;
    size_t refresh_point = 0;

public:
    /**
     * @param ora Pointer to an `tatami::Oracle` to be used for predictions.
     * @param max_size Total size of all slabs to store in the cache.
     * This may be zero, in which case no caching should be performed.
     */
    OracularVariableSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > ora, size_t max_size) : 
        oracle(std::move(ora)), 
        total(oracle->total()),
        max_size(max_size) 
    {} 

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularVariableSlabCache(const OracularVariableSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularVariableSlabCache& operator=(const OracularVariableSlabCache&) = delete;

    /**
     * @cond
     */
    // Move operators are still okay as pointers still point to the moved vectors.
    // see https://stackoverflow.com/questions/43988553/stdvector-stdmove-and-pointer-invalidation.
    OracularVariableSlabCache& operator=(OracularVariableSlabCache&&) = default;
    OracularVariableSlabCache(OracularVariableSlabCache&&) = default;

    // Might as well define this.
    ~OracularVariableSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method is intended to be called when `max_size = 0`, to provide callers with the oracle predictions for non-cached extraction of data.
     * Calls to this method should not be intermingled with calls to its overload below; the latter should only be called when `max_size > 0`.
     *
     * @return The next prediction from the oracle.
     */
    Index_ next() {
        return oracle->get(counter++);
    }

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     * This method should only be called if `max_size > 0` in the constructor; otherwise, no slabs are actually available and cannot be returned.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Efunction_ Function to compute the estimated size of a slab.
     * @tparam Afunction_ Function to compute the actual size of a slab.
     * @tparam Cfunction_ Function to create a new slab.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the slab containing `i`, and (2) the index of row/column `i` inside that slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2 and an offset of 1.
     * @param estimated_size Function that accepts `j`, an `Id_` containing the slab identifier.
     * It should return the size of the slab as a non-negative `Size_`.
     * @param actual_size Function that accepts `j`, an `Id_` containing the slab identifier; and `slab`, a populated `const Slab_&` instance corresponding to `j`.
     * It should return the actual size of the slab as a non-negative `Size_` that is no greater than `estimated_size(j)`.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts two arguments - `to_populate` and `to_reuse`.
     * - The `to_populate` argument is a `std::vector<std::pair<Id_, Slab_*> >&` specifying the slabs to be populated.
     *   The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     *   The second `Slab_*` element contains a pointer to a `Slab_` returned by `create()`.
     * - The `to_reuse` argument is a `std::vector<std::pair<Id_, Slab_*> >&` specifying the cached slabs that were re-used in the upcoming set of predictions.
     *   The elements of each pair are interpreted in the same manner as `to_populate`. 
     * .
     * The `populate` function should iterate over `to_populate` and fill each `Slab_` with the contents of the corresponding slab.
     * It may optionally iterate over `to_reuse` to defragment the cache in order to free up enough space for the new contents in `to_populate` -
     * this is typically required for `Slab_` implementations that point into a single shared buffer, as repeated variable-size allocations can fragment that buffer.
     * Note that neither `to_populate` or `to_reuse` vector is guaranteed to be sorted. 
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Efunction_, class Afunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Efunction_ estimated_size, Afunction_ actual_size, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == last_slab_id && last_slab) {
            return std::make_pair(last_slab, slab_info.second);
        }
        last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (counter - 1 == refresh_point) {
            Size_ used_size;

            auto ccIt = current_cache.find(slab_info.first);
            if (ccIt != current_cache.end()) {
                auto slab_ptr = ccIt->second;
                used_size = actual_size(slab_info.first, *slab_ptr);
                future_cache[slab_info.first] = slab_ptr;
                to_reuse.emplace_back(slab_info.first, slab_ptr);
                current_cache.erase(ccIt);
            } else {
                used_size = estimated_size(slab_info.first);
                requisition_new_slab(slab_info.first);
            }

            auto last_future_slab_id = slab_info.first;
            to_reuse.clear();

            while (++refresh_point < total) {
                auto future_index = oracle->get(refresh_point);
                auto future_slab_info = identify(future_index);
                if (last_future_slab_id != future_slab_info.first) {
                    if (future_cache.find(future_slab_info.first) == future_cache.end()) {
                        auto ccIt = current_cache.find(future_slab_info.first);
                        if (ccIt != current_cache.end()) {
                            auto slab_ptr = ccIt->second;
                            used_size += actual_size(slab_info.first, *slab_ptr);
                            if (used_size > max_size) {
                                break;
                            } 
                            future_cache[future_slab_info.first] = slab_ptr;
                            to_reuse.emplace_back(future_slab_info.first, slab_ptr);
                            current_cache.erase(ccIt);
                        } else {
                            used_size += estimated_size(future_slab_info.first);
                            if (used_size > max_size) {
                                break;
                            } 
                            requisition_new_slab(future_slab_info.first);
                        }
                    }
                }
            }

            auto cIt = current_cache.begin();
            for (auto a : in_need) {
                if (cIt != current_cache.end()) {
                    to_populate.emplace_back(a, cIt->second);
                    future_cache[a] = cIt->second;
                    ++cIt;
                } else {
                    all_slabs.push_back(create());
                    to_populate.emplace_back(a, &(all_slabs.back()));
                }
            }
            in_need.clear();

            for (; cIt != current_cache.end(); ++cIt) {
                free_pool.emplace_back(cIt->second);
            }

            populate(to_populate, to_reuse);
            to_populate.clear();
            to_reuse.clear();

            current_cache.clear();
            current_cache.swap(future_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = current_cache.find(slab_info.first);
        last_slab = ccIt->second;
        return std::make_pair(last_slab, slab_info.second);
    }

private:
    void requisition_new_slab(Id_ slab_id) {
        if (!free_pool.empty()) {
            auto slab_ptr = free_pool.back();
            future_cache[slab_id] = slab_ptr;
            free_pool.pop_back();
            to_populate.emplace_back(slab_id, slab_ptr);
        } else {
            future_cache[slab_id] = NULL;
            in_need.push_back(slab_id);
        }
    }

public:
    /**
     * @return Maximum total size of all slabs in the cache.
     */
    size_t get_max_size() const {
        return max_size;
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
