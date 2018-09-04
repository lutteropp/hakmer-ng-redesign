/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include <unordered_set>

size_t longestCommonPrefix(const std::string& seq, size_t start1, size_t start2, unsigned int lTop) {
	size_t res = 0;
	for (size_t i = 0; i < seq.size(); ++i) {
		if (start1 + i >= seq.size() || start2 + i >= seq.size()) {
			break;
		}
		if (seq[start1 + i] == seq[start2 + i]) {
			res++;
			if (res >= lTop) {
				return res; // bail because lTop is largest k-mer size we're gonna search
			}
		} else {
			break;
		}
	}
	return res;
}

size_t posToTaxon(size_t pos, const std::vector<std::pair<size_t, size_t> >& taxonCoords) {
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (pos >= taxonCoords[i].first && pos <= taxonCoords[i].second) {
			return i;
		}
	}
	return std::string::npos;
}

bool canStay(size_t pos, const std::vector<std::pair<size_t, size_t> >& taxonCoords, const std::vector<size_t>& wantedTaxa) {
	for (size_t i = 0; i < wantedTaxa.size(); ++i) {
		if (pos >= taxonCoords[i].first && pos <= taxonCoords[i].second) {
			return true;
		}
	}
	return false;
}

std::pair<std::vector<size_t>, std::vector<size_t> > shrinkArrays(const std::string& T, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, const std::vector<std::pair<size_t, size_t> >& taxonCoords, const std::vector<size_t>& wantedTaxa,
		const Options& options) {
	std::vector<size_t> resSA;
	std::vector<size_t> resLCP;

	bool recomputeNeeded = false;

	for (size_t i = 0; i < SA.size(); ++i) {
		if (canStay(i, taxonCoords, wantedTaxa)) {
			resSA.push_back(SA[i]);
			size_t lcpVal = lcp[i];
			if (recomputeNeeded) {
				// recompute lcpVal
				size_t lTop = std::min(T.size(), options.maxK);
				lcpVal = longestCommonPrefix(T, SA[i - 1], SA[i], lTop);
				recomputeNeeded = false;
			}
			resLCP.push_back(lcpVal);
		} else {
			recomputeNeeded = true;
		}
	}
	return std::make_pair(resSA, resLCP);
}

bool acceptSeed(size_t actSAPos, size_t matchCount, size_t k, size_t nTax, const std::vector<size_t>& SA, PresenceChecker& presenceChecker,
		const std::vector<std::pair<size_t, size_t> >& taxonCoords, const Options& options) {
	if (matchCount > nTax) { // easy test for paralogy
		return false;
	}
	if (matchCount < options.minTaxaPerBlock) { // easy test for not enough taxa
		return false;
	}
	// more complicated check for paralogy
	std::unordered_set<size_t> takenTaxa;
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		size_t taxID = posToTaxon(SA[i], taxonCoords);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}
	// check for presence/absence
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		if (!presenceChecker.isFree(SA[i], SA[i] + k - 1)) {
			return false;
		}
	}
	return true;
}

size_t countMatches(size_t actSAPos, const std::vector<size_t>& lcp, size_t k) {
	size_t matchCount = 1;
	for (size_t i = actSAPos + 1; i < lcp.size(); ++i) {
		if (lcp[i] >= k) {
			matchCount++;
		} else {
			break;
		}
	}
	return matchCount;
}

// TODO: Re-add mismatches and indels in seeds
SeededBlock nextSeededBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	SeededBlock block(nTax);
	bool foundBlock = false;
	if (actSAPos >= SA.size()) {
		return block;
	}

	size_t startPos = SA[actSAPos];
	size_t k = options.minK;
	if (startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1)) { // early stop
		actSAPos++;
		block = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	}

	k = options.minK;
	size_t matchCount = countMatches(actSAPos, lcp, k);

	while (matchCount >= options.minTaxaPerBlock) {
		if (acceptSeed(actSAPos, matchCount, k, nTax, SA, presenceChecker, taxonCoords, options)) {
			foundBlock = true;
			for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
				block.addTaxon(posToTaxon(SA[i], taxonCoords), SA[i], SA[i] + k - 1);
			}
			//presenceChecker.reserveSeededBlock(block);
			break;
		} else {
			if (k == options.maxK) { // no further extension of seed length
				break;
			}
			if (startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // newly added character would be already taken anyway
				break;
			}
			k++;
			matchCount = countMatches(actSAPos, lcp, k);
		}
	}

	if (!foundBlock) {
		actSAPos++;
		block = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	} else {
		actSAPos++;
	}

	return block;
}

std::vector<SeededBlock> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	std::vector<SeededBlock> res;
	size_t actSAPos = 0;
	while (actSAPos < SA.size()) {
		res.push_back(nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options));
	}
	return res;
}

// TODO: Re-add dynamic extension of flanking regions
ExtendedBlock extendBlock(const SeededBlock& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options) {
	ExtendedBlock block(seededBlock, nTax);
	if (!options.dynamicFlanks) {
		size_t flankWidth = options.flankWidth;
		for (size_t i = 0; i < nTax; ++i) {
			if (block.hasTaxon(i)) {
				block.setLeftFlankSize(i, flankWidth);
				block.setRightFlankSize(i, flankWidth);
			}
		}
	} else { // TODO: Implement me.

	}
	return block;
}

std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	std::vector<ExtendedBlock> res;
	size_t actSAPos = 0;
	while (actSAPos < SA.size()) {
		SeededBlock seededBlock = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
		ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
		// check if the extended block can still be accepted.
		if (presenceChecker.isFine(extendedBlock)) {
			presenceChecker.reserveExtendedBlock(extendedBlock);
			res.push_back(extendedBlock);
		}
	}
	return res;
}
