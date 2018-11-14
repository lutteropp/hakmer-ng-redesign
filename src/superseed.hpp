/*
 * superseed.hpp
 *
 *  Created on: Nov 6, 2018
 *      Author: sarah
 */

#pragma once

#include <unordered_set>
#include <iostream>
#include "seed.hpp"
#include "alignment/msa_wrapper.hpp"

/* stores a set of merged seed regions. */

class Superseed {
public:
	Superseed(size_t nTax, const Seed& mySeed) {
		taxonCoords.resize(nTax);
		for (size_t i = 0; i < nTax; ++i) {
			taxonCoords[i] = mySeed.getTaxonCoords(i);
			if (mySeed.hasTaxon(i)) {
				taxIDs.insert(i);
			}
		}
		mySeeds.push_back(mySeed);
	}
	void addSeed(const Seed& seedToAdd) { // TODO: The interaction between reverse-complement and the merged taxon coordinates seems to be still wrong.
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
		mySeeds.push_back(seedToAdd);
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

	bool orderCompatible(const Superseed& otherSeed, size_t revCompStartIdx) const {
		if (nSharedTax(otherSeed) == 0) {
			return true;
		}
		for (size_t i = 0; i < this->mySeeds.size(); ++i) {
			for (size_t j = 0; j < otherSeed.mySeeds.size(); ++j) {
				if (!mySeeds[i].orderCompatible(otherSeed.mySeeds[j], revCompStartIdx)) {
					return false;
				}
			}
		}
		return true;
	}

	bool overlap(const Superseed& otherSeed, size_t revCompStartIdx) const {
		if (nSharedTax(otherSeed) == 0) {
			return false;
		}
		for (size_t i = 0; i < this->mySeeds.size(); ++i) {
			for (size_t j = 0; j < otherSeed.mySeeds.size(); ++j) {
				if (mySeeds[i].overlap(otherSeed.mySeeds[j], revCompStartIdx)) {
					return true;
				}
			}
		}
		return false;
	}

	double score(const Superseed& otherSeed, size_t revCompStartIdx, size_t maxAllowedDistance, size_t minSharedTax) const {
		size_t sharedTaxa = nSharedTax(otherSeed);
		if (sharedTaxa < minSharedTax) {
			return std::numeric_limits<double>::infinity();
		}
		size_t dist = distance(otherSeed, revCompStartIdx, maxAllowedDistance, sharedTaxa);
		if (dist == std::numeric_limits<size_t>::infinity()) {
			return std::numeric_limits<double>::infinity();
		} else {
			size_t nTaxA = this->getNTaxInBlock();
			size_t nTaxB = otherSeed.getNTaxInBlock();
			return (1.0 - (2.0 * sharedTaxa / (nTaxA + nTaxB))) * dist;
		}
	}

	std::pair<size_t, size_t> computeForwardStrandCoordinates(const std::pair<size_t, size_t>& coords, size_t revCompStartPos) const {
		size_t thisForwardFirst = coords.first;
		size_t thisForwardSecond = coords.second;
		if (thisForwardFirst >= revCompStartPos) {
			thisForwardFirst = revCompStartPos - (thisForwardFirst - revCompStartPos);
			thisForwardSecond = revCompStartPos - (thisForwardSecond - revCompStartPos);
		}
		return std::make_pair(std::min(thisForwardFirst, thisForwardSecond), std::max(thisForwardFirst, thisForwardSecond));
	}

	size_t distance(const Superseed& otherSeed, size_t revCompStartIdx, size_t maximumDistance, size_t sharedTax) const {
		if (sharedTax == 0) {
			return std::numeric_limits<size_t>::infinity();
		} else {
			if (!orderCompatible(otherSeed, revCompStartIdx) || overlap(otherSeed, revCompStartIdx)) {
				return std::numeric_limits<size_t>::infinity();
			}

			size_t maxDist = 0;
			for (size_t i = 0; i < taxonCoords.size(); ++i) {
				if (this->hasTaxon(i) && otherSeed.hasTaxon(i)) {
					size_t actDist;
					std::pair<size_t, size_t> thisForward = computeForwardStrandCoordinates(taxonCoords[i], revCompStartIdx);
					std::pair<size_t, size_t> otherForward = computeForwardStrandCoordinates(otherSeed.taxonCoords[i], revCompStartIdx);

					if (thisForward.second < otherForward.first) {
						actDist = otherForward.first - thisForward.second;
					} else if (otherForward.second < thisForward.first) {
						actDist = thisForward.first - otherForward.second;
					} else {
						actDist = 0;
						//throw std::runtime_error("This should not happen - do the seeds overlap?");
						std::cout << "my seeds for this taxon:" << "\n";
						for (size_t j = 0; j < mySeeds.size(); ++j) {
							if (mySeeds[j].hasTaxon(i)) {
								std::pair<size_t, size_t> coords = computeForwardStrandCoordinates(mySeeds[j].getTaxonCoords(i),
										revCompStartIdx);
								std::cout << coords.first << " - " << coords.second << "\n";
							}
						}
						std::cout << "other seeds for this taxon: " << "\n";
						for (size_t j = 0; j < otherSeed.mySeeds.size(); ++j) {
							if (mySeeds[j].hasTaxon(i)) {
								std::pair<size_t, size_t> coords = computeForwardStrandCoordinates(otherSeed.mySeeds[j].getTaxonCoords(i),
										revCompStartIdx);
								std::cout << coords.first << " - " << coords.second << "\n";
							}
						}
					}
					maxDist = std::max(maxDist, actDist);
					if (maxDist >= maximumDistance) {
						return std::numeric_limits<size_t>::infinity();
					}
				}
			}
			return maxDist;
		}
	}
	size_t nSharedTax(const Superseed& otherSeed) const {
		size_t cnt = 0;
		for (size_t tID: taxIDs) {
			if (otherSeed.taxIDs.find(tID) != otherSeed.taxIDs.end()) {
				cnt++;
			}
		}
		return cnt;
	}
	double pSharedTax(const Superseed& otherSeed) const {
		return (double) 2 * nSharedTax(otherSeed) / (this->getNTaxInBlock() + otherSeed.getNTaxInBlock());
	}
	void clear() {
		mySeeds.clear();
		mySeeds.shrink_to_fit();
		taxonCoords.clear();
		taxonCoords.shrink_to_fit();
		taxIDs.clear();
	}
private:
	std::vector<Seed> mySeeds;
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	std::unordered_set<size_t> taxIDs;
};
