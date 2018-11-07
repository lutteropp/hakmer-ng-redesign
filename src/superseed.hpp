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
	const std::vector<Seed>& getMySeeds() const {
		return mySeeds;
	}

	const bool orderCompatible(const Superseed& otherSeed, size_t revCompStartIdx) const {
		// This has to also take into account reverse complement matches.
		if (nSharedTax(otherSeed) == 0) {
			return true;
		}
		bool orderSet = false;
		bool thisComesFirst;
		for (size_t i = 0; i < taxonCoords.size(); ++i) {
			if (hasTaxon(i) && otherSeed.hasTaxon(i)) {
				bool thisIsFirstNow;
				if (taxonCoords[i].first < revCompStartIdx && otherSeed.taxonCoords[i].first < revCompStartIdx) { // both are on forward strand
					if (taxonCoords[i].first < otherSeed.taxonCoords[i].first) {
						thisIsFirstNow = true;
					} else {
						thisIsFirstNow = false;
					}
				} else if (taxonCoords[i].first >= revCompStartIdx && otherSeed.taxonCoords[i].first >= revCompStartIdx) { // both are on reverse-complement strand
					if (taxonCoords[i].first > otherSeed.taxonCoords[i].first) {
						thisIsFirstNow = true;
					} else {
						thisIsFirstNow = false;
					}
				} else if (taxonCoords[i].first < revCompStartIdx && otherSeed.taxonCoords[i].first >= revCompStartIdx) { // this on forward strand, the other one on reverse-complement strand
					size_t rcBackToForward = revCompStartIdx - (otherSeed.taxonCoords[i].first - revCompStartIdx);
					if (taxonCoords[i].first < rcBackToForward) {
						thisIsFirstNow = true;
					} else {
						thisIsFirstNow = false;
					}
				} else { // this is on reverse-complement strand, the other one is on forward strand
					size_t rcBackToForward = revCompStartIdx - (taxonCoords[i].first - revCompStartIdx);
					if (rcBackToForward < otherSeed.taxonCoords[i].first) {
						thisIsFirstNow = true;
					} else {
						thisIsFirstNow = false;
					}
				}

				if (!orderSet) {
					thisComesFirst = thisIsFirstNow;
					orderSet = true;
				} else {
					if (thisIsFirstNow != thisComesFirst) {
						return false;
					}
				}
			}
		}
		return true;
	}

	size_t distance(const Superseed& otherSeed, size_t revCompStartIdx) const {
		if (nSharedTax(otherSeed) == 0) {
			return std::numeric_limits<size_t>::infinity();
		} else {
			if (!orderCompatible(otherSeed, revCompStartIdx)) {
				return std::numeric_limits<size_t>::infinity();
			}
			size_t maxDist = 0;
			for (size_t i = 0; i < taxonCoords.size(); ++i) {
				if (this->hasTaxon(i) && otherSeed.hasTaxon(i)) {
					size_t actDist;
					size_t thisFirstForward = taxonCoords[i].first;
					size_t thisSecondForward = taxonCoords[i].second;
					if (thisFirstForward >= revCompStartIdx) { // this is on reverse-complement strand
						thisFirstForward = revCompStartIdx - (thisFirstForward - revCompStartIdx);
						thisSecondForward = revCompStartIdx - (thisSecondForward - revCompStartIdx);
					}
					size_t otherFirstForward = otherSeed.taxonCoords[i].first;
					size_t otherSecondForward = otherSeed.taxonCoords[i].second;
					if (otherFirstForward >= revCompStartIdx) { // other is on reverse-complement strand
						otherFirstForward = revCompStartIdx - (otherFirstForward - revCompStartIdx);
						otherSecondForward = revCompStartIdx - (otherSecondForward - revCompStartIdx);
					}
					if (thisSecondForward < otherFirstForward) {
						actDist = otherFirstForward - thisSecondForward;
					} else if (otherSecondForward < thisFirstForward) {
						actDist = thisFirstForward - otherSecondForward;
					} else {
						throw std::runtime_error("This should not happen - do the seeds overlap?");
					}
					maxDist = std::max(maxDist, actDist);
				}
			}
			return maxDist;
		}
	}
	size_t nSharedTax(const Superseed& otherSeed) const {
		size_t count = 0;
		for (size_t i = 0; i < taxonCoords.size(); ++i) {
			if (this->hasTaxon(i) && otherSeed.hasTaxon(i)) {
				count++;
			}
		}
		return count;
	}
	double pSharedTax(const Superseed& otherSeed) const {
		return (double) 2 * nSharedTax(otherSeed) / (this->getNTaxInBlock() + otherSeed.getNTaxInBlock());
	}
private:
	std::vector<Seed> mySeeds;
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	std::unordered_set<size_t> taxIDs;
};
