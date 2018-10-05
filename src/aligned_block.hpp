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
#include "options.hpp"

class AlignedBlock {
public:
	AlignedBlock(const ExtendedBlock& extendedBlock, size_t nTax);
	bool isAligned() const;
	bool hasTaxon(size_t taxID) const;
	void alignMAFFT(const std::string& T, const Options& options);
	std::vector<std::string> getAlignment() const;
	void setAlignment(const std::vector<std::string>& alignmentOfPresentTaxa);
	size_t getAlignmentWidth() const;
	std::pair<size_t, size_t> getTaxonCoordsWithFlanks(size_t taxID) const;
	std::vector<size_t> getTaxonIDsInBlock() const;
	size_t getSeedSize() const;
	std::vector<double> getPairwiseNormalizedDistances(const Options& options);
private:
	bool aligned;
	ExtendedBlock myBlock;
	size_t nTax;
	std::vector<std::string> alignment;
};
