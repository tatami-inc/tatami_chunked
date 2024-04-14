#ifndef TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file SubsettedOracleSlabCache.hpp
 * @brief Create a oracle-aware cache with subsets.
 */

namespace tatami_chunked {

/**
 * Type of subset selection.
 * Used to determine the subsets to extract in `SubsettedOracleSlabCache::SubsettedOracleSlabSubset`.
 * 
 * - `FULL`: all rows/columns. 
 * - `BLOCK`: a contiguous block of rows/columns. 
 * - `INDEX`: an indexed subset of rows/columns.
 */
enum class SubsettedOracleSelection : char { FULL, BLOCK, INDEX };

/**
 * @brief Details on the subset to extract in `SubsettedOracleSlabCache`.
 * @tparam Index_ Type of row/column index produced by the oracle.
 */
template<typename Index_>
struct SubsettedOracleSlabSubset {
    /**
     * Type of subset selection to extract from the slab.
     */
    SubsettedOracleSelection selection;

    /**
     * Row/column index representing the start of the block within the slab to be extracted.
     * Only used if `selection` is set to `SubsettedOracleSelection::BLOCK`.
     *
     * Note that the index is relative to the start of the slab, not to the matrix containing the slab,
     * i.e., if the slab consists of rows 10-20 and we want to extract row 11, this will be reported here as an index of 1.
     */
    Index_ block_start;

    /**
     * Length of the block within the slab to be extracted.
     * Only used if `selection` is set to `SubsettedOracleSelection::BLOCK`.
     */
    Index_ block_length;

    /**
     * Indices within the slab to be extracted.
     * Guaranteed to be sorted and unique.
     * Only used if `selection` is set to `SubsettedOracleSelection::INDEX`.
     *
     * Note that all indices is relative to the start of the slab, not to the matrix containing the slab,
     * i.e., if the slab consists of rows 10-20 and we want to extract row 11, this will be reported here as an index of 1.
     */
    std::vector<Index_> indices;

    /**
     * Mapping of indices to be extracted to their positions inside `indices`.
     * All values of `indices` are present as keys here where `mapping[indices[i]] = i`.
     * Only used if `selection` is set to `SubsettedOracleSelection::INDEX`.
     */
    std::unordered_map<Index_, Index_> mapping;

private:
    Index_ block_end = 0;

    void fill_mapping() {
        for (size_t i = 0, end = indices.size(); i < end; ++i) {
            mapping[indices[i]] = i;
        }
    }

public:
    /**
     * @cond
     */
    void set(Index_ i) {
        selection = SubsettedOracleSelection::BLOCK;
        block_start = i;
        block_end = i + 1;
        indices.clear();
        mapping.clear();
    }

    void add(Index_ i) {
        if (selection == SubsettedOracleSelection::FULL) {
            return;
        }

        if (selection == SubsettedOracleSelection::BLOCK) {
            if (i == block_end) {
                block_end = i + 1;
                return;

            } else if (i + 1 == block_start) {
                block_start = i;
                return;

            } else if (i >= block_start && i < block_end) {
                return;
            }

            selection = SubsettedOracleSelection::INDEX;
            indices.resize(block_end - block_start);
            std::iota(indices.begin(), indices.end(), block_start);
            fill_mapping();
        }

        if (mapping.find(i) == mapping.end()) {
            mapping[i] = indices.size();
            indices.push_back(i);
        }
    }

    void finalize() {
        if (selection == SubsettedOracleSelection::BLOCK) {
            block_length = block_end - block_start;
        } else if (selection == SubsettedOracleSelection::INDEX) {
            if (!std::is_sorted(indices.begin(), indices.end())) {
                std::sort(indices.begin(), indices.end());
                fill_mapping();
            }
        }
    }
    /**
     * @endcond
     */
};

/**
 * @brief Oracle-aware cache for slabs, plus subsets.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slab subsets.
 * This is similar to the `OracleSlabCache` except that it remembers the subset of rows/columns that were requested for each slab.
 * Slab extractors can use this information to optimize the slab population process by ignoring rows/columns that are not needed.
 */
template<typename Id_, typename Index_, class Slab_> 
class SubsettedOracleSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > oracle;
    size_t total;
    size_t counter = 0;

    Index_ last_slab_id = 0;
    Slab_* last_slab = NULL;

    size_t max_slabs;
    std::vector<Slab_> all_slabs;
    std::unordered_map<Id_, Slab_*> current_cache, future_cache;

    std::vector<SubsettedOracleSlabSubset<Index_> > all_subset_details;
    std::vector<SubsettedOracleSlabSubset<Index_>*> free_subset_details;
    std::unordered_map<Id_, SubsettedOracleSlabSubset<Index_>*> close_future_subset_cache, far_future_subset_cache;
    size_t close_refresh_point = 0;
    size_t far_refresh_point = 0;
    Id_ far_slab_id;
    Index_ far_slab_offset;

    std::vector<std::pair<Id_, SubsettedOracleSlabSubset<Index_>*> > to_reassign;
    std::vector<std::tuple<Id_, Slab_*, SubsettedOracleSlabSubset<Index_>*> > to_populate;

public:
    /**
     * @param ora Pointer to an `tatami::Oracle` to be used for predictions.
     * @param num_slabs Maximum number of slabs to store.
     */
    SubsettedOracleSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > ora, size_t num_slabs) :
        oracle(std::move(ora)), 
        total(oracle->total()),
        max_slabs(num_slabs)
    {
        all_slabs.reserve(max_slabs);
        current_cache.reserve(max_slabs);
        future_cache.reserve(max_slabs);
        close_future_subset_cache.reserve(max_slabs);
        far_future_subset_cache.reserve(max_slabs);

        all_subset_details.resize(max_slabs * 2);
        for (auto& as : all_subset_details) {
            free_subset_details.push_back(&as);
        }
    }

    /**
     * Deleted as the cache holds persistent pointers.
     */
    SubsettedOracleSlabCache(const SubsettedOracleSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent pointers.
     */
    SubsettedOracleSlabCache& operator=(const SubsettedOracleSlabCache&) = delete;

    /**
     * @cond
     */
    // Move operators are still okay as pointers still point to the moved vectors,
    // see https://stackoverflow.com/questions/43988553/stdvector-stdmove-and-pointer-invalidation.
    SubsettedOracleSlabCache(SubsettedOracleSlabCache&&) = delete;
    SubsettedOracleSlabCache& operator=(SubsettedOracleSlabCache&&) = delete;

    // Might as well define this.
    ~SubsettedOracleSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method is intended to be called when `num_slabs = 0`, to provide callers with the oracle predictions for non-cached extraction of data.
     * Calls to this method should not be intermingled with calls to its overload below; the latter should only be called when `max_slabs > 0`.
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
     * @param populate Function that accepts a `std::vector<std::pair<Id_, Slab_*, SubsettedOracleSlabSubset<Index_>*> >&` specifying the slabs to be populated.
     * The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     * The second `Slab_*` element specifies the object which to store the contents of each slab.
     * The thid `SubsettedOracleSlabSubset<Index_>*` element contains information about the subset of each slab that is required.
     * This function should iterate over the vector and populate the desired subset of each slab.
     * Note that the vector is not guaranteed to be sorted. 
     *
     * @return Pair containing (1) a pointer to a cached slab and (2) the index of the next predicted row/column inside the retrieved slab.
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
        if (counter - 1 == close_refresh_point) {
            if (all_slabs.empty()) {
                // This section only runs once, at the start, to populate the close_future_subset_cache.
                requisition_subset_close(slab_info.first, slab_info.second);
                size_t used_slabs = 1;

                while (++close_refresh_point < total) {
                    auto future_index = oracle->get(close_refresh_point);
                    auto future_slab_info = identify(future_index);
                    auto cfcIt = close_future_subset_cache.find(future_slab_info.first);
                    if (cfcIt != close_future_subset_cache.end()) {
                        cfcIt->second->add(future_slab_info.second);
                    } else if (used_slabs < max_slabs) {
                        requisition_subset_close(future_slab_info.first, future_slab_info.second);
                        ++used_slabs;
                    } else {
                        far_slab_id = future_slab_info.first;
                        far_slab_offset = future_slab_info.second;
                        break;
                    }
                }

                far_refresh_point = close_refresh_point;
            } else {
                close_refresh_point = far_refresh_point;
            }

            // Populating the far future cache. 
            if (far_refresh_point < total) {
                requisition_subset_far(far_slab_id, far_slab_offset);
                size_t used_slabs = 1;

                while (++far_refresh_point < total) {
                    auto future_index = oracle->get(far_refresh_point);
                    auto future_slab_info = identify(future_index);
                    auto ffcIt = far_future_subset_cache.find(future_slab_info.first);
                    if (ffcIt != far_future_subset_cache.end()) {
                        ffcIt->second->add(future_slab_info.second);
                    } else if (used_slabs < max_slabs) {
                        requisition_subset_far(future_slab_info.first, future_slab_info.second);
                        ++used_slabs;
                    } else {
                        far_slab_id = future_slab_info.first;
                        far_slab_offset = future_slab_info.second;
                        break;
                    }
                }
            }

            // Reusing slabs from current_cache; these should all have FULL selections already.
            for (auto& cf : close_future_subset_cache) {
                auto cIt = current_cache.find(cf.first);
                if (cIt == current_cache.end()) {
                    to_reassign.emplace_back(cf.first, cf.second);
                } else {
                    future_cache[cf.first] = cIt->second;
                    current_cache.erase(cIt);
                }
            }

            // Creating new slabs for everything that's left.
            auto cIt = current_cache.begin();
            for (auto a : to_reassign) {
                Slab_* slab_ptr;
                if (cIt == current_cache.end()) {
                    all_slabs.emplace_back(create());
                    slab_ptr = &(all_slabs.back());
                } else {
                    slab_ptr = cIt->second;
                    ++cIt;
                }
                future_cache[a.first] = slab_ptr;
                to_populate.emplace_back(a.first, slab_ptr, a.second);
            }
            to_reassign.clear();

            for (auto p : to_populate) {
                std::get<2>(p)->finalize();
            }
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

            // Putting the no-longer-used subset pointers back in the free pool
            // before we swap the close and far futures.
            for (auto& cfc : close_future_subset_cache) {
                free_subset_details.push_back(cfc.second);
            }
            close_future_subset_cache.clear();
            close_future_subset_cache.swap(far_future_subset_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = current_cache.find(slab_info.first);
        last_slab = ccIt->second;
        return std::make_pair(last_slab, slab_info.second);
    }

private:
    void requisition_subset_close(Id_ slab_id, Index_ slab_offset) {
        auto selected = free_subset_details.back();
        selected->set(slab_offset);
        close_future_subset_cache[slab_id] = selected;
        free_subset_details.pop_back();
    }

    void requisition_subset_far(Id_ slab_id, Index_ slab_offset) {
        auto selected = free_subset_details.back();
        selected->set(slab_offset);
        far_future_subset_cache[slab_id] = selected;
        free_subset_details.pop_back();

        // If a slab is still being used in the far future, it might continue
        // to be used in an even further future, in which case we need to do a
        // FULL extraction just to be safe.
        auto cfcIt = close_future_subset_cache.find(slab_id);
        if (cfcIt != close_future_subset_cache.end()) {
            selected->selection = SubsettedOracleSelection::FULL;
            cfcIt->second->selection = SubsettedOracleSelection::FULL;
        }
    }
};

}

#endif
