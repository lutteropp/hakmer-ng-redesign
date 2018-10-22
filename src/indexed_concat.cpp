#include "indexed_concat.hpp"

IndexedTaxonCoords::IndexedTaxonCoords(const std::string& label, const std::vector<std::string>& contigs, size_t coordOffset, bool revComp) {
	this->label = label;

	size_t nC = contigs.size();
	if (revComp) {
		nC *= 2;
	}
	contigCoords.resize(nC);

	contigCoords[0] = std::make_pair(coordOffset, contigs[0].size() - 1 + coordOffset);
	for (size_t i = 1; i < contigs.size(); ++i) {
		contigCoords[i] = std::make_pair(contigCoords[i - 1].second + 2, contigCoords[i - 1].second + 2 + contigs[i].size() - 1);
	}
	if (revComp) {
		for (size_t i = 0; i < contigs.size(); ++i) {
			contigCoords[i + contigs.size()] = std::make_pair(contigCoords[i + contigs.size() - 1].second + 2, contigCoords[i + contigs.size() - 1].second + 2 + contigs[i].size() - 1);
		}
	}

	firstCoord = contigCoords[0].first;
	lastCoord = contigCoords[0].second;
	for (size_t i = 1; i < contigCoords.size(); ++i) {
		firstCoord = std::min(firstCoord, contigCoords[i].first);
		lastCoord = std::max(lastCoord, contigCoords[i].second);
	}
}
std::string IndexedTaxonCoords::getLabel() const {
	return label;
}
size_t IndexedTaxonCoords::getFirstCoord() const {
	return firstCoord;
}
size_t IndexedTaxonCoords::getLastCoord() const {
	return lastCoord;
}

size_t IndexedTaxonCoords::getContigStart(size_t contigIdx) const {
	return contigCoords[contigIdx].first;
}
size_t IndexedTaxonCoords::getContigEnd(size_t contigIdx) const {
	return contigCoords[contigIdx].second;
}

std::pair<size_t, size_t> IndexedTaxonCoords::getContigCoords(size_t contigIdx) const {
	return contigCoords[contigIdx];
}

bool IndexedTaxonCoords::contains(size_t pos) const {
	for (size_t i = 0; i < contigCoords.size(); ++i) {
		if (pos >= contigCoords[i].first && pos <= contigCoords[i].second) {
			return true;
		}
	}
	return false;
}

IndexedConcatenatedSequence::IndexedConcatenatedSequence(const std::string& seq, const std::vector<IndexedTaxonCoords>& coords,
		const Options& options) {
	concatenatedSeq = seq;
	taxonCoords = coords;
	sa.buildSuffixArray(seq, seq.size(), options);
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		taxonLabels.push_back(taxonCoords[i].getLabel());
	}
}

const std::vector<size_t>& IndexedConcatenatedSequence::getSuffixArray() const {
	return sa.getSA();
}
const std::vector<size_t>& IndexedConcatenatedSequence::getLcpArray() const {
	return sa.getLCP();
}

size_t IndexedConcatenatedSequence::nTax() const {
	return taxonCoords.size();
}
IndexedTaxonCoords IndexedConcatenatedSequence::getTaxonCoords(size_t taxonIdx) const {
	return taxonCoords[taxonIdx];
}
const std::string& IndexedConcatenatedSequence::getConcatenatedSeq() const {
	return concatenatedSeq;
}

size_t IndexedConcatenatedSequence::getConcatSize() const {
	return concatenatedSeq.size();
}

const std::vector<std::string>& IndexedConcatenatedSequence::getTaxonLabels() const {
	return taxonLabels;
}

const std::vector<IndexedTaxonCoords>& IndexedConcatenatedSequence::getTaxonCoords() const {
	return taxonCoords;
}
