/*
 * extended_block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include "seeded_block.hpp"

class ExtendedBlock {
public:
	ExtendedBlock(const SeededBlock& seededBlock, size_t nTax);
	size_t getLeftFlankSize(size_t taxID) const;
	size_t getRightFlankSize(size_t taxID) const;
	void setLeftFlankSize(size_t taxID, size_t val);
	void setRightFlankSize(size_t taxID, size_t val);
	void incrementAllLeftFlanks();
	void incrementAllRightFlanks();
	void decrementAllLeftFlanks();
	void decrementAllRightFlanks();
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID) const;
	bool hasTaxon(size_t taxID) const;
private:
	SeededBlock mySeededBlock;
	std::vector<size_t> leftFlankSizes;
	std::vector<size_t> rightFlankSizes;
};
