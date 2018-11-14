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

/*
 * Store seed coordinates for each included taxon
 */

class Seed {
public:
	Seed(size_t nTax);
	void addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord);
	size_t getNTaxInBlock() const;
	std::pair<size_t, size_t> getTaxonCoords(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	size_t getSeedSize() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	void increaseTaxonCoordsRight();
	void decreaseTaxonCoordsLeft();
	void decreaseTaxonCoordsRight();
	void increaseTaxonCoordsLeft();
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
private:
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	std::vector<size_t> taxIDs;
	size_t n;
	size_t k;
};
