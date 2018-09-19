/*
 * aligned_block.hpp
 *
 *  Created on: Sep 4, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <string>

#include "extended_block.hpp"

class AlignedBlock {
public:
	AlignedBlock(const ExtendedBlock& extendedBlock, size_t nTax);
	bool isAligned() const;
	bool hasTaxon(size_t taxID) const;
	void align(const std::string& T);
	std::vector<std::string> getAlignment() const;
	size_t getAlignmentWidth() const;
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID) const;
	std::vector<size_t> getTaxonIDsInBlock() const;
private:
	bool aligned;
	ExtendedBlock myBlock;
	size_t nTax;
	std::vector<std::string> alignment;
};
