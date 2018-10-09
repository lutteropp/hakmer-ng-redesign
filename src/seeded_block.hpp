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

class SeededBlock {
public:
	SeededBlock(size_t nTax);
	void addTaxon(size_t taxID, size_t firstCoord, size_t lastCoord);
	size_t getNTaxInBlock() const;
	std::pair<size_t, size_t> getTaxonCoords(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	size_t getSeedSize() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	void increaseTaxonCoordsRight();
	void decreaseTaxonCoordsLeft();
private:
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	size_t n;
	size_t k;
};
