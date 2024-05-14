#ifndef TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include <type_traits>
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
 * Each slab is defined as the set of chunks required to read an element of the target dimension (or a contiguous block/indexed subset thereof) from a `tatami::Matrix`.
 * This cache is similar to `OracularSlabCache` but enables improved cache utilization when the slabs vary in size.
 * For example, the number of non-zero entries in a sparse matrix might vary between slabs,
 * so the cache could be optimized to fit more slabs into memory when they have fewer non-zeros.
 *
 * The size of each slab is defined by `Size_`, which can be any non-negative measure of slab size.
 * This could be the number of non-zero elements, or the number of dimension elements, or the size of the slab in bytes, etc.,
 * as long as its interpretation is consistent between slabs and with the `max_size` used in the constructor.
 * Users can also differentiate between the estimated and actual size of the slab, if the latter is not known until after the slab has been loaded into memory,
 * e.g., the number of non-zero entries in a file-backed sparse matrix.
 *
 * When implementing `Slab_`,  we generally suggest using a common memory pool that is referenced by each `Slab_` instance.
 * This guarantees that the actual cache size does not exceed the limit associated with `max_size` when `Slab_` instances are re-used for different slabs.
 * (Otherwise, if each `Slab_` allocates its own memory, re-use of an instance may cause its allocation to increase to the size of the largest encountered slab.)
 * Callers may need to occasionally defragment the pool to ensure that enough memory is available for loading new slabs.
 */
template<typename Id_, typename Index_, class Slab_, typename Size_> 
class OracularVariableSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > oracle;
    size_t total;
    size_t counter = 0;

    Index_ last_slab_id = 0;
    size_t last_slab_num = -1;

    Size_ max_size, used_size = 0;
    std::vector<Slab_> all_slabs;

    // We need to hold an offset into 'all_slabs' rather than a pointer, as
    // 'all_slabs' might be reallocated upon addition of new slabs, given that
    // we don't know the maximum number of slabs ahead of time.
    std::unordered_map<Id_, size_t> current_cache, future_cache;
    std::vector<std::pair<Id_, size_t> > to_populate, to_reuse;
    std::vector<Id_> in_need;
    std::vector<size_t> free_pool;
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
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted index of a single element on the target dimension.
     * This should return a pair containing:
     * 1. An `Id_`, the identifier of the slab containing `i`.
     *    This is typically defined as the index of the slab on the target dimension.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2.
     * 2. An `Index_`, the index of row/column `i` inside that slab.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would yield an offset of 1.
     * @param estimated_size Function that accepts `j`, an `Id_` containing the slab identifier.
     * It should return the size of the slab as a non-negative `Size_`.
     * @param actual_size Function that accepts `j`, an `Id_` containing the slab identifier; and `slab`, a populated `const Slab_&` instance corresponding to `j`.
     * It should return the actual size of the slab as a non-negative `Size_` that is no greater than `estimated_size(j)`.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts three arguments - `to_populate`, `to_reuse` and `all_slabs`.
     * - The `to_populate` argument is a `std::vector<std::pair<Id_, size_t> >&` specifying the slabs to be populated.
     *   The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     *   The second `size_t` element specifies the entry of `all_slabs` containing the corresponding `Slab_` instance, as returned by `create()`.
     *   This argument can be modified in any manner.
     *   It is not guaranteed to be sorted.
     * - The `to_reuse` argument is a `std::vector<std::pair<Id_, size_t> >&` specifying the cached slabs that were re-used in the upcoming set of predictions.
     *   The elements of each pair are interpreted in the same manner as `to_populate`. 
     *   This argument can be modified in any manner.
     *   It is not guaranteed to be sorted.
     * - The `all_slabs` argument is a `std::vector<Slab_>&` containing all slabs in the cache.
     *   This may include instances that are not referenced by `to_populate` or `to_reuse`.
     *   Each element of this argument can be modified but the length should not change.
     * .
     * The `populate` function should iterate over `to_populate` and fill each `Slab_` with the contents of the corresponding slab.
     * Optionally, callers may use `to_reuse` to defragment the already-in-use parts of the cache, in order to free up enough space for new data from `to_populate`.
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Efunction_, class Afunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Efunction_ estimated_size, Afunction_ actual_size, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == last_slab_id && last_slab_num != static_cast<size_t>(-1)) {
            return std::make_pair(all_slabs.data() + last_slab_num, slab_info.second);
        }
        last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (counter - 1 == refresh_point) {
            // Note that, for any given populate cycle, the first prediction's
            // slab cannot already be in the cache, otherwise it would have
            // incorporated into the previous cycle. So we can skip some code.
            used_size = estimated_size(slab_info.first);
            requisition_new_slab(slab_info.first);

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

                auto ccIt = current_cache.find(future_slab_info.first);
                if (ccIt != current_cache.end()) {
                    size_t slab_num = ccIt->second;
                    auto candidate = used_size + actual_size(future_slab_info.first, all_slabs[slab_num]);
                    if (candidate > max_size) {
                        break;
                    } 
                    used_size = candidate;
                    future_cache[future_slab_info.first] = slab_num;
                    to_reuse.emplace_back(future_slab_info.first, slab_num);
                    current_cache.erase(ccIt);
                } else {
                    auto candidate = used_size + estimated_size(future_slab_info.first);
                    if (candidate > max_size) {
                        break;
                    } 
                    used_size = candidate;
                    requisition_new_slab(future_slab_info.first);
                }
            }

            auto cIt = current_cache.begin();
            for (auto a : in_need) {
                if (cIt != current_cache.end()) {
                    size_t slab_num = cIt->second;
                    to_populate.emplace_back(a, slab_num);
                    future_cache[a] = slab_num;
                    ++cIt;
                } else {
                    size_t slab_num = all_slabs.size();
                    all_slabs.push_back(create());
                    to_populate.emplace_back(a, slab_num);
                    future_cache[a] = slab_num;
                }
            }
            in_need.clear();

            for (; cIt != current_cache.end(); ++cIt) {
                free_pool.emplace_back(cIt->second);
            }

            populate(to_populate, to_reuse, all_slabs);
            to_populate.clear();
            to_reuse.clear();

            current_cache.clear();
            current_cache.swap(future_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = current_cache.find(slab_info.first);
        last_slab_num = ccIt->second;
        return std::make_pair(all_slabs.data() + last_slab_num, slab_info.second);
    }

private:
    void requisition_new_slab(Id_ slab_id) {
        if (!free_pool.empty()) {
            auto slab_num = free_pool.back();
            future_cache[slab_id] = slab_num;
            free_pool.pop_back();
            to_populate.emplace_back(slab_id, slab_num);
        } else {
            future_cache[slab_id] = 0;
            in_need.push_back(slab_id);
        }
    }

public:
    /**
     * @return Maximum total size of the cache.
     * This is the same as the `max_size` used in the constructor.
     */
    size_t get_max_size() const {
        return max_size;
    }

    /**
     * @return Current usage across all slabs in the cache.
     * This should be interpreted as an upper bound on usage if there is a difference between estimated and actual slab sizes.
     */
    size_t get_used_size() const {
        return used_size;
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
