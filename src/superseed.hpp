/*
 * superseed.hpp
 *
 *  Created on: Nov 6, 2018
 *      Author: sarah
 */

#pragma once

#include <unordered_set>
#include "seed.hpp"
#include "msa_wrapper.hpp"

/* stores a set of merged seed regions. */

class Superseed {
public:
	Superseed(size_t nTax, Seed& mySeed) {
		taxonCoords.resize(nTax);
		for (size_t i = 0; i < nTax; ++i) {
			taxonCoords[i] = mySeed.getTaxonCoords(i);
			if (mySeed.hasTaxon(i)) {
				taxIDs.insert(i);
			}
		}
	}
	void addSeed(Seed& seedToAdd) {
		for (size_t i = 0; i < taxonCoords.size(); ++i) {
			if (seedToAdd.hasTaxon(i)) {
				if (taxonCoords[i].first == std::string::npos) {
					taxonCoords[i] = seedToAdd.getTaxonCoords(i);
					taxIDs.insert(i);
				} else {
					taxonCoords[i].first = std::min(taxonCoords[i].first, seedToAdd.getTaxonCoords(i).first);
					taxonCoords[i].second = std::max(taxonCoords[i].second, seedToAdd.getTaxonCoords(i).second);
				}
			}
		}
	}
	std::pair<size_t, size_t> getTaxonCoords(size_t taxID) const {
		return taxonCoords[taxID];
	}
	size_t getNTaxInBlock() const {
		return taxIDs.size();
	}
	bool hasTaxon(size_t taxID) const {
		return (taxIDs.find(taxID) != taxIDs.end());
	}
private:
	std::vector<Seed> mySeeds;
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	std::unordered_set<size_t> taxIDs;
};
