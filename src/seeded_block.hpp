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
	size_t getN() const;
	std::vector<std::pair<size_t, size_t> > getTaxonCoords() const;
	bool hasTaxon(size_t taxID) const;
private:
	std::vector<std::pair<size_t, size_t> > taxonCoords;
	size_t n;
};
