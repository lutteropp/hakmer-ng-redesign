/*
 * IndexedConcat.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

#include "options.hpp"
#include "indexing/suffix_array_classic.hpp"

class IndexedTaxonCoords {
	std::string label;
	size_t firstCoord;
	size_t lastCoord;
	std::vector<std::pair<size_t, size_t> > contigCoords;
public:
	std::string getLabel() const;
	size_t getFirstCoord() const;
	size_t getLastCoord() const;
	std::pair<size_t, size_t> getContigCoords(size_t contigIdx) const;
	size_t getContigStart(size_t contigIdx) const;
	size_t getContigEnd(size_t contigIdx) const;
	size_t nContigs() const;
	IndexedTaxonCoords(const std::string& label, const std::vector<std::string>& contigs, size_t coordOffset, const Options& options);
};

// TODO: Re-integrate external indexing via FM index

class IndexedConcatenatedSequence {
	SuffixArrayClassic sa;
	std::string concatenatedSeq;
	std::vector<IndexedTaxonCoords> taxonCoords;
public:
	IndexedConcatenatedSequence(const std::string& seq, const std::vector<IndexedTaxonCoords>& coords, const Options& options); // internal-indexing variant
	const std::vector<size_t>& getSuffixArray() const;
	const std::vector<size_t>& getLcpArray() const;
	size_t getConcatSize() const;
	size_t nTax() const;
	IndexedTaxonCoords getTaxonCoords(size_t taxonIdx) const;
	const std::string& getConcatenatedSeq() const;
};
