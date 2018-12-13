#include "indexed_concat.hpp"

IndexedTaxonCoords::IndexedTaxonCoords(const std::string& label, const std::vector<std::string>& contigs, size_t coordOffset) {
	this->label = label;

	size_t nC = contigs.size();
	contigCoords.resize(nC);

	contigCoords[0] = std::make_pair(coordOffset, contigs[0].size() - 1 + coordOffset);
	for (size_t i = 1; i < contigs.size(); ++i) {
		contigCoords[i] = std::make_pair(contigCoords[i - 1].second + 2, contigCoords[i - 1].second + 2 + contigs[i].size() - 1);
	}

	firstCoord = contigCoords[0].first;
	lastCoord = contigCoords[0].second;
	for (size_t i = 1; i < contigCoords.size(); ++i) {
		firstCoord = std::min(firstCoord, contigCoords[i].first);
		lastCoord = std::max(lastCoord, contigCoords[i].second);
	}
}
size_t IndexedTaxonCoords::getTotalLength() const {
	return 1 + lastCoord - firstCoord;
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

size_t longestCommonPrefix(const std::string& seq, size_t start1, size_t start2) {
	size_t res = 0;
	for (size_t i = 0; i < seq.size(); ++i) {
		if (start1 + i >= seq.size() || start2 + i >= seq.size()) {
			break;
		}
		if (ambiguousMatch(seq[start1 + i], seq[start2 + i])) {
			res++;
		} else {
			break;
		}
	}
	return res;
}

std::pair<std::vector<size_t>, std::vector<size_t> > IndexedConcatenatedSequence::shrinkArrays(size_t wantedTaxon) {
	std::vector<size_t> resSA;
	std::vector<size_t> resLCP;

	size_t lastITaken = 0;

	bool recomputeNeeded = false;
	for (size_t i = 0; i < suffixArray.size(); ++i) {
		if (taxonCoords[wantedTaxon].contains(suffixArray[i])) {
			resSA.push_back(suffixArray[i]);
			size_t lcpVal = lcpArray[i];
			if (recomputeNeeded) {
				// recompute lcpVal
				lcpVal = longestCommonPrefix(concatenatedSeq, resSA[resSA.size() - 2], resSA[resSA.size() - 1]);
				recomputeNeeded = false;
			}
			resLCP.push_back(lcpVal);

			lastITaken = i;
		} else {
			recomputeNeeded = true;
		}
	}
	return std::make_pair(resSA, resLCP);
}

size_t IndexedConcatenatedSequence::posToTaxonInternal(size_t pos, bool revComp) const {
	if (revComp && pos >= concatenatedSeq.size() / 2) {
		pos = concatenatedSeq.size() - pos - 1;
	}
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (taxonCoords[i].contains(pos)) {
			return i;
		}
	}
	return std::numeric_limits<uint16_t>::infinity();
}

IndexedConcatenatedSequence::IndexedConcatenatedSequence(const std::string& seq, const std::vector<IndexedTaxonCoords>& coords,
		bool protein, const Options& options) {
	concatenatedSeq = seq;
	taxonCoords = coords;
	{
		SuffixArrayClassic sa;
		sa.buildSuffixArray(seq, seq.size(), options);
		suffixArray = sa.getSA();
		lcpArray = sa.getLCP();
	}

	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		taxonLabels.push_back(taxonCoords[i].getLabel());
	}
	sequenceDataSize = 0;
	for (size_t i = 0; i < seq.size(); ++i) {
		if (seq[i] != '$') {
			sequenceDataSize++;
		}
	}
	proteinData = protein;
	// precompute posToTaxonArray
	std::cout << "Precomputing posToTaxon array...\n";
	posToTaxonArray.resize(suffixArray.size());
#pragma omp parallel
	{
#pragma omp for
		for (size_t i = 0; i < suffixArray.size(); ++i) {
			size_t tID = posToTaxonInternal(suffixArray[i], options.reverseComplement);
			posToTaxonArray[suffixArray[i]] = tID;
		}
	}

	std::cout << "Computing per-taxon SA and LCP arrays...\n";
	perTaxonArrays.resize(taxonCoords.size());
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		perTaxonArrays[i] = shrinkArrays(i);
	}
}

size_t IndexedConcatenatedSequence::getSequenceDataSize() const {
	return sequenceDataSize;
}

const std::vector<size_t>& IndexedConcatenatedSequence::getSuffixArray() const {
	return suffixArray;
}
const std::vector<size_t>& IndexedConcatenatedSequence::getLcpArray() const {
	return lcpArray;
}

const std::vector<size_t>& IndexedConcatenatedSequence::getSuffixArray(size_t taxonID) const {
	return perTaxonArrays[taxonID].first;
}

const std::vector<size_t>& IndexedConcatenatedSequence::getLcpArray(size_t taxonID) const {
	return perTaxonArrays[taxonID].second;
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

bool IndexedConcatenatedSequence::isProtein() const {
	return proteinData;
}

size_t IndexedConcatenatedSequence::posToTaxon(size_t pos) const {
	return posToTaxonArray[pos];
}

