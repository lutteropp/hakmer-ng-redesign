/*
 * seeded_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <string>

#include "alignment/simple_msa.hpp"

/*
 * Store seed coordinates for each included taxon
 */

class Seed {
public:
	Seed(size_t nTax);
	void addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord);
	void removeTaxon(size_t taxID);
	size_t getNTaxInBlock() const;
	SimpleCoords getSeedCoords(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	double getAverageSeedSize() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	void increaseAllTaxonCoordsRight();
	void decreaseAllTaxonCoordsLeft();
	void decreaseAllTaxonCoordsRight();
	void increaseAllTaxonCoordsLeft();
	void increaseTaxonCoordsRight(size_t taxID);
	void decreaseTaxonCoordsLeft(size_t taxID);
	void decreaseTaxonCoordRight(size_t taxID);
	void increaseTaxonCoordLeft(size_t taxID);
	bool operator <(const Seed& str) const {
		if (n == str.n) {
			return k < str.k;
		} else {
			return n < str.n;
		}
	}
	bool operator >(const Seed& str) const {
		if (n == str.n) {
			return k > str.k;
		} else {
			return n > str.n;
		}
	}
	bool orderCompatible(const Seed& other, size_t revCompStartPos) const;
	bool overlap(const Seed& other, size_t revCompStartPos) const;
	size_t distance(const Seed& other, size_t revCompStartPos) const;
	bool allSeedsSameSize() const;
	size_t getSeedSize(size_t taxonID) const;

	void addGapLeft(size_t taxonID);
	void addGapRight(size_t taxonID);
	std::vector<SimpleCoords> getSeedCoords() const;
private:
	std::vector<SimpleCoords> seedCoords;
	std::vector<size_t> taxIDs;
	size_t n;
	size_t k;
};
