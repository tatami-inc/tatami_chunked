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
    size_t max_predictions;
    size_t max_slabs;

public:
    /**
     * @brief A cached slab.
     */
    struct CachedSlab {
        /**
         * Identity of the slab.
         */
        Id_ self;

        /**
         * Contents of the slab.
         */
        Slab_ contents;

        /**
         * Indices inside the slab to be retrieved.
         * If empty, all indices should be retrieved.
         * If non-empty, this is guaranteed to be sorted and unique.
         */
        std::vector<Index_> indices;

    public:
        /**
         * @cond
         */
        CachedSlab(Id_ s, Slab_ c) : self(s), contents(std::move(c)) {}
        /**
         * @endcond
         */

    private:
        std::unordered_set<Index_> present;
    };

private:
    std::list<CachedSlab> slab_cache;
    std::list<CachedSlab> free_cache;
    typedef typename std::list<CachedSlab>::iterator cache_iterator;
    std::unordered_map<Id_, std::pair<cache_iterator, bool> > slab_exists;

    std::list<CachedSlab> mock_slab_cache;
    std::list<CachedSlab> mock_free_cache;

    std::vector<std::pair<const Slab_*, Index_> > predictions_made;
    size_t predictions_fulfilled = 0;

    struct PredictionRound {
        std::vector<std::pair<const Slab_*, Index_> > predictions_made;
        std::vector<Id_> slabs_used;
        void clear() {
            slabs_used.clear();
            predictions_made.clear();
        }
    };
    std::list<PredictionRound> next_rounds;
    std::list<CachedSlab> free_rounds;

    std::vector<CachedSlab*> slabs_to_fill; 

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_slabs Maximum number of slabs to store.
     */
    OracleSlabCache(std::unique_ptr<tatami::Oracle<Index_> > oracle, size_t per_iteration, size_t num_slabs) :
        prediction_stream(std::move(oracle)), 
        max_predictions(per_iteration),
        max_slabs(num_slabs)
    {
        slabs_to_fill.reserve(max_slabs);
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
     * @tparam Afunction_ Function to allocate memory to a mock slab, to turn it into a non-mock instance.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the slab containing `i`, and (2) the index of row/column `i` inside that slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2 and an offset of 1.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * @param populate Function that accepts a vector of `CachedSlab` pointers and fills them with the contents of the relevant slabs.
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Sfunction_, class Rfunction_, class Afunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Afunction_ allocate, Pfunction_ populate) {
        if (predictions_made.size() > predictions_fulfilled) {
            return predictions_made[predictions_fulfilled++];
        }

        if (next_rounds.empty()) {
            // If we're just starting, then we fill it up the first prediction round.
            predictions_made.reserve(max_predictions);
            next_rounds.resize(1);
            auto& baseline = next_rounds.front().slabs_used;
            size_t used = 0;

            for (size_t p = 0; p < max_predictions; ++p) {
                Index_ current;
                if (!prediction_stream.next(current)) {
                    break;
                }

                auto slab_id = identify(current);
                auto curslab = slab_id.first;
                auto curindex = slab_id.second;

                auto it = slab_exists.find(curslab);
                if (it == slab_exists.end()) {
                    if (used == max_slabs) {
                        prediction_stream.back();
                        break;
                    }

                    slab_cache.emplace_back(create());
                    it = slab_cache.end();
                    --it;
                    slab_exists[curslab] = std::make_pair(scIt, true);
                    baseline.push_back(curslab);
                }

                next_predictions_made.emplace_back(&(((it->second).first)->contents), curindex);
            }

        } else {
            // Otherwise we cycle out the previous round's information and we move the next round's information in.
            auto& previous_round = next_rounds.front();
            for (auto x : previous_round.slabs_used) {
                slab_exists[x].second = false;
            }
            previous_round.clear();
            free_rounds.splice(free_rounds.end(), next_rounds, next_rounds.begin());

            auto& current_round = next_rounds.front();
            for (auto x : current_round.slabs_used) {
                slab_exists[x].second = true;
            }
            predictions_made.swap(current_round.predictions_made);
        }

        constexpr int max_rounds = 10;
        for (int round = 0; round < max_rounds; ++round) {
            if (free_rounds.empty()) {
                free_rounds.resize(1);
            }
            next_rounds.splice(next_rounds.end(), free_rounds, free_rounds.begin());
            auto& next_baseline = next_rounds.front().slabs_used;
            auto& next_predictions_made = next_rounds.front().predictions_made;
            size_t used = 0;
            bool overlaps = false;
            bool finished = false;

            for (size_t p = 0; p < max_predictions; ++p) {
                Index_ current;
                if (!prediction_stream.next(current)) {
                    finished = true;
                    break;
                }

                auto slab_id = identify(current);
                auto curslab = slab_id.first;
                auto curindex = slab_id.second;

                auto it = slab_exists.find(curslab);
                if (it == slab_exists.end()) {
                    if (used == max_slabs) {
                        prediction_stream.back();
                        break;
                    }

                    slab_cache.emplace_back(create()); // OOPS. What do we do here?
                    it = slab_cache.end();
                    --it;
                    slab_exists[curslab] = std::make_pair(scIt, false);
                    next_baseline.push_back(curslab);
                } else if (!overlaps) {
                    overlaps = (it->second).second;
                }

                next_predictions_made.emplace_back(&(((it->second).first)->contents), curindex);
            }

            if (!overlaps || finished) {
                break;
            }
        }

    }
};

}

#endif
