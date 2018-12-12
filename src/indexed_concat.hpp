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
	bool contains(size_t pos) const;
	IndexedTaxonCoords(const std::string& label, const std::vector<std::string>& contigs, size_t coordOffset);
	size_t getTotalLength() const;
};

// TODO: Re-integrate external indexing via FM index

class IndexedConcatenatedSequence {
	std::string concatenatedSeq;
	std::vector<IndexedTaxonCoords> taxonCoords;
	std::vector<std::string> taxonLabels;
	size_t sequenceDataSize;
	std::vector<size_t> suffixArray;
	std::vector<size_t> lcpArray;
	bool proteinData;
	std::vector<uint16_t> posToTaxonArray;
	size_t posToTaxonInternal(size_t pos, bool revComp) const;
public:
	IndexedConcatenatedSequence(const std::string& seq, const std::vector<IndexedTaxonCoords>& coords, bool protein, const Options& options); // internal-indexing variant
	const std::vector<size_t>& getSuffixArray() const;
	const std::vector<size_t>& getLcpArray() const;
	size_t getConcatSize() const;
	size_t nTax() const;
	size_t getSequenceDataSize() const;
	IndexedTaxonCoords getTaxonCoords(size_t taxonIdx) const;
	const std::string& getConcatenatedSeq() const;
	const std::vector<std::string>& getTaxonLabels() const;
	const std::vector<IndexedTaxonCoords>& getTaxonCoords() const;
	bool isProtein() const;
	size_t posToTaxon(size_t pos) const;
};
