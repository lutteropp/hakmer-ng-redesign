/*
 * seeded_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <vector>

#include "alignment/simple_coords.hpp"

/*
 * Store seed coordinates for each included taxon
 */

class Seed {
public:
	Seed(size_t nTax);
	void addTaxon(size_t saPos, size_t taxID, size_t firstCoord, size_t lastCoord);
	void blockTaxon(size_t taxID);
	void removeTaxon(size_t taxID);
	size_t getNTaxInBlock() const;
	SimpleCoords getSeedCoords(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	double getAverageSeedSize() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	void increaseAllTaxonCoordsRight();
	void decreaseAllTaxonCoordsLeft();
	void increaseTaxonCoordsRight(size_t taxID);
	void decreaseTaxonCoordsLeft(size_t taxID);
	void decreaseTaxonCoordRight(size_t taxID);
	void increaseTaxonCoordLeft(size_t taxID);
	bool operator <(const Seed& str) const {
		/*if (k >= 42 && str.k < 42) {
			// this is always smaller
			return 0;
		} else if (str.k >= 42 && k < 42) {
			// this is always larger
			return 1;
		}*/

		if (taxIDs.size() == str.taxIDs.size()) {
			if (subRate == str.subRate) {
				return k < str.k;
			} else {
				return subRate > str.subRate;
			}
		} else {
			return taxIDs.size() < str.taxIDs.size();
		}
	}
	bool operator >(const Seed& str) const {
		/*if (k >= 42 && str.k < 42) {
			// this is always larger
			return 1;
		} else if (str.k >= 42 && k < 42) {
			// this is always smaller
			return 0;
		}*/

		if (taxIDs.size() == str.taxIDs.size()) {
			if (subRate == str.subRate) {
				return k > str.k;
			} else {
				return subRate < str.subRate;
			}
		} else {
			return taxIDs.size() > str.taxIDs.size();
		}
	}
	bool orderCompatible(const Seed& other, size_t revCompStartPos) const;
	bool overlap(const Seed& other, size_t revCompStartPos) const;
	size_t distance(const Seed& other, size_t revCompStartPos) const;
	bool allSeedsSameSize() const;
	size_t getSeedSize(size_t taxonID) const;

	void addGapLeft(size_t taxonID);
	void addGapRight(size_t taxonID);
	void removeGapLeft(size_t taxonID);
	void removeGapRight(size_t taxonID);
	std::vector<SimpleCoords> getSeedCoords() const;
	size_t getOriginalK() const;
	size_t getFirstSAPos() const;
	void setFirstSAPos(size_t pos);
	std::vector<size_t> getSAPositions() const;

	void setSubRate(double rate);
	double getSubRate() const;
	const std::vector<size_t> getBlockedTaxa() const;
private:
	void cleanTaxIDs();
	std::vector<SimpleCoords> seedCoords;
	std::vector<size_t> taxIDs;
	std::vector<size_t> saPositions;
	std::vector<size_t> blockedTaxa;
	size_t k;
	size_t originalK;
	size_t firstSAPos;
	double subRate;
};
