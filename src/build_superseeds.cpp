/*
 * superblocks.cpp
 *
 *  Created on: Nov 7, 2018
 *      Author: sarah
 */

#include "build_superseeds.hpp"
#include "seed.hpp"
#include "presence_checker.hpp"

#include <queue>
#include <iostream>
#include <vector>

class MyComparator {
public:
	int operator()(const std::tuple<size_t, size_t, double>& t1, const std::tuple<size_t, size_t, double>& t2) {
		return std::get<2>(t1) > std::get<2>(t2);
	}
};

Superseed mergeSuperseeds(Superseed& s1, Superseed& s2, size_t nTax, const Options& Options) {
	const std::vector<Seed>& s1Seeds = s1.getMySeeds();
	const std::vector<Seed>& s2Seeds = s2.getMySeeds();
	Superseed merged(nTax, s1Seeds[0]);
	for (size_t i = 1; i < s1Seeds.size(); ++i) {
		merged.addSeed(s1Seeds[i]);
	}
	for (size_t i = 0; i < s2Seeds.size(); ++i) {
		merged.addSeed(s2Seeds[i]);
	}
	s1.clear();
	s2.clear();
	return merged;
}

std::vector<Superseed> buildSuperseeds(const std::vector<Seed>& seeds, const std::string& T, const PresenceChecker& presenceChecker,
		size_t nTax, const Options& options) {
	size_t revCompStartIdx = T.size();
	if (options.reverseComplement) {
		revCompStartIdx /= 2;
	}
	// first, put each superseed into the priority queue.
	std::cout << "seeds.size() before building superseeds: " << seeds.size() << "\n";
	std::vector<Superseed> superseeds;
	superseeds.reserve(seeds.size());
	for (size_t i = 0; i < seeds.size(); ++i) {
		Superseed ss(nTax, seeds[i]);
		superseeds.push_back(ss);
	}
	std::vector<bool> active(superseeds.size(), true);

	std::priority_queue<std::tuple<size_t, size_t, double>, std::vector<std::tuple<size_t, size_t, double>>, MyComparator> pq;
	// initially fill the priority queue
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < superseeds.size(); ++i) {
		if (!active[i]) {
			continue;
		}
		for (size_t j = i + 1; j < superseeds.size(); ++j) {
			if (!active[j]) {
				continue;
			}
			double dist = superseeds[i].score(superseeds[j], revCompStartIdx, options.maxAllowedSuperseedDistance, options.minSharedSuperseedTax);
			if (dist < std::numeric_limits<double>::infinity()) { // can be merged
				std::tuple<size_t, size_t, double> entry = std::make_tuple(i, j, dist);
#pragma omp critical
				pq.push(entry);
			}
		}
	}

	while (!pq.empty()) {
		std::tuple<size_t, size_t, double> entry = pq.top();
		if (active[std::get<0>(entry)] && active[std::get<1>(entry)]) {
			// merge the two blocks
			size_t firstIdx = std::get<0>(entry);
			size_t secondIdx = std::get<1>(entry);
			std::cout << "Merging blocks with indices " << firstIdx << " and " << secondIdx << "...\n";
			Superseed mergedSeed = mergeSuperseeds(superseeds[firstIdx], superseeds[secondIdx], nTax, options);
			active[firstIdx] = false;
			active[secondIdx] = false;
			superseeds.push_back(mergedSeed);
			active.push_back(true);

			// update pq entries
			size_t i = superseeds.size() - 1;
			for (size_t j = 0; j < superseeds.size() - 1; ++j) {
				if (!active[j]) {
					continue;
				}
				double dist = superseeds[i].score(superseeds[j], revCompStartIdx, options.maxAllowedSuperseedDistance, options.minSharedSuperseedTax);
				if (dist < std::numeric_limits<double>::infinity()) { // can be merged
					std::tuple<size_t, size_t, double> entry = std::make_tuple(i, j, dist);
					pq.push(entry);
				}
			}

		}
		pq.pop();
	}

	std::vector<Superseed> res;
	for (size_t i = 0; i < superseeds.size(); ++i) {
		if (active[i]) {
			res.push_back(superseeds[i]);
		}
	}
	std::cout << "Number of superseeds: " << res.size() << "\n";
	std::cout << "Finished building superseeds.\n";
	return res;
}
