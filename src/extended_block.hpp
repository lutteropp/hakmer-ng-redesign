/*
 * extended_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <utility>
#include <vector>

#include "alignment/simple_coords.hpp"
#include "seed.hpp"

class ExtendedBlock {
public:
	ExtendedBlock(const Seed& seededBlock, size_t nTax);
	double getAverageLeftFlankSize() const;
	double getAverageRightFlankSize() const;
	void setLeftFlankSize(size_t val);
	void setRightFlankSize(size_t val);
	void growLeftFlank(size_t tID);
	void growRightFlank(size_t tID);
	void addLeftFlankGap(size_t tID);
	void addRightFlankGap(size_t tID);
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID) const;
	std::pair<size_t, size_t> getTaxonCoordsWithoutFlanks(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
	double getAverageSeedSize() const;
	size_t getNTaxInBlock() const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	size_t getTotalBasesUsed() const;
	void removeTaxon(size_t tID);

	std::vector<SimpleCoords> getLeftFlankCoords() const;
	std::vector<SimpleCoords> getSeedCoords() const;
	std::vector<SimpleCoords> getRightFlankCoords() const;
private:
	Seed mySeededBlock;
	std::vector<SimpleCoords> leftFlankCoords;
	std::vector<SimpleCoords> rightFlankCoords;
};
