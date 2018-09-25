/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"

#include <unordered_set>
#include <cmath>

#include "distance_estimator.hpp"
#include "indexing/suffix_array_classic.hpp"

size_t posToTaxon(size_t pos, const std::vector<std::pair<size_t, size_t> >& taxonCoords, size_t concatSize, bool revComp) {
	if (pos >= concatSize && revComp) {
		pos = concatSize - 1 - pos;
	}
	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (pos >= taxonCoords[i].first && pos <= taxonCoords[i].second) {
			return i;
		}
	}
	return std::string::npos;
}

bool canStay(size_t pos, const std::vector<std::pair<size_t, size_t> >& taxonCoords, const std::vector<size_t>& wantedTaxa) {
	for (size_t i = 0; i < wantedTaxa.size(); ++i) {
		if (pos >= taxonCoords[wantedTaxa[i]].first && pos <= taxonCoords[wantedTaxa[i]].second) {
			return true;
		}
	}
	return false;
}

// Shrink arrays... needed for quartet-based stuff.
std::pair<std::vector<size_t>, std::vector<size_t> > shrinkArrays(const IndexedConcatenatedSequence& concat,
		const std::vector<std::pair<size_t, size_t> >& taxonCoords, const std::vector<size_t>& wantedTaxa, const Options& options) {
	std::vector<size_t> resSA;
	std::vector<size_t> resLCP;

	bool recomputeNeeded = false;

	for (size_t i = 0; i < concat.getSuffixArray().size(); ++i) {
		if (canStay(concat.getSuffixArray()[i], taxonCoords, wantedTaxa)) {
			resSA.push_back(concat.getSuffixArray()[i]);
			size_t lcpVal = concat.getLcpArray()[i];
			if (recomputeNeeded) {
				// recompute lcpVal
				size_t lTop = std::min(concat.getConcatSize(), options.maxK);
				lcpVal = longestCommonPrefix(concat.getConcatenatedSeq(), resSA[resSA.size() - 2], resSA[resSA.size() - 1], lTop);
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
		const std::vector<std::pair<size_t, size_t> >& taxonCoords, size_t concatSize, const Options& options) {
	if (matchCount > nTax /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
		return false;
	}
	if (matchCount < options.minTaxaPerBlock) { // easy test for not enough taxa
		return false;
	}
	// more complicated check for paralogy
	std::unordered_set<size_t> takenTaxa;
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		size_t taxID = posToTaxon(SA[i], taxonCoords, concatSize, options.reverseComplement);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}

	size_t flankOffset = 0;
	if (!options.dynamicFlanks) {
		flankOffset = options.flankWidth;
	}

	// check for presence/absence
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		if (flankOffset > SA[i] || SA[i] + k - 1 + flankOffset >= concatSize) {
			return false;
		}

		if (!presenceChecker.isFree(SA[i] - flankOffset, SA[i] + k - 1 + flankOffset)) {
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
// TODO: Something here is going wrong. Maybe get rid of all these performance-killing recursive calls?
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
	// find the next suitable SA pos
	while ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1)) && actSAPos < SA.size()) {
		actSAPos++;
		startPos = SA[actSAPos];
	}
	size_t matchCount = countMatches(actSAPos, lcp, k);
	while (matchCount >= options.minTaxaPerBlock) {
		if (acceptSeed(actSAPos, matchCount, k, nTax, SA, presenceChecker, taxonCoords, T.size(), options)) {
			foundBlock = true;

			for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
				block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
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

	actSAPos++;
	if (!foundBlock) {
		block = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	}
	return block;
}

std::vector<SeededBlock> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	std::vector<SeededBlock> res;
	size_t actSAPos = 0;
	while (actSAPos < SA.size()) {
		SeededBlock bl = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
		if (bl.getSeedSize() > 0) {
			res.push_back(nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options));
		} else {
			break;
		}
	}
	return res;
}

void swap(std::vector<double>& vec, size_t i, size_t j) {
	double tmp = vec[i];
	vec[i] = vec[j];
	vec[j] = tmp;
}

double deltaScore(const std::array<double, 6>& pairwiseDist, const Options& options) {
	// This happens if two sequences are too divergent, this is, relative uncorrected distance > 0.75
	for (size_t i = 0; i < pairwiseDist.size(); ++i) {
		if (std::isnan(pairwiseDist[i]) || std::isinf(pairwiseDist[i])) {
			return 1;
		}
	}

	// Check for 3 or more identical sequences, as this leads to star topology
	if (pairwiseDist[0] == 0 && pairwiseDist[1] == 0) { // d(a,b) and d(a,c) -> a,b,c are equal
		return 1;
	} else if (pairwiseDist[0] == 0 && pairwiseDist[2] == 0) { // d(a,b) and d(a,d) -> a,b,d are equal
		return 1;
	} else if (pairwiseDist[1] == 0 && pairwiseDist[2] == 0) { // d(a,c) and d(a,d) -> a,c,d are equal
		return 1;
	} else if (pairwiseDist[3] == 0 && pairwiseDist[4] == 0) { // d(b,c) and d(b,d) -> b,c,d are equal
		return 1;
	}

	double d12_34 = pairwiseDist[0] + pairwiseDist[5];
	double d13_24 = pairwiseDist[1] + pairwiseDist[4];
	double d14_23 = pairwiseDist[2] + pairwiseDist[3];

	std::vector<double> distances = { d12_34, d13_24, d14_23 };
	if (distances[0] > distances[1])
		swap(distances, 0, 1);
	if (distances[0] > distances[2])
		swap(distances, 0, 2);
	if (distances[1] > distances[2])
		swap(distances, 1, 2);

	double score;
	if (distances[0] == distances[2]) { // star topology
		score = 1.0;
	} else {
		score = (distances[2] - distances[1]) / (distances[2] - distances[0]);
	}
	if (std::isnan(score)) {
		throw std::runtime_error("The delta score is nan!!!\n");
	}
	return score;
}

double averageDeltaScore(std::vector<std::vector<HammingDistanceEstimator> >& pairwiseDistanceEstimators, size_t nTaxBlock,
		const Options& options) {
	double sum = 0.0;
	size_t numQuartets = nTaxBlock * (nTaxBlock - 1) * (nTaxBlock - 2) * (nTaxBlock - 3) / 24;
	for (size_t i = 0; i < nTaxBlock - 3; ++i) {
		for (size_t j = i + 1; j < nTaxBlock - 2; ++j) {
			for (size_t k = j + 1; k < nTaxBlock - 1; ++k) {
				for (size_t l = k + 1; l < nTaxBlock; ++l) {
					std::array<double, 6> pairwiseDist = { pairwiseDistanceEstimators[i][j].distance(),
							pairwiseDistanceEstimators[i][k].distance(), pairwiseDistanceEstimators[i][l].distance(),
							pairwiseDistanceEstimators[j][k].distance(), pairwiseDistanceEstimators[j][l].distance(),
							pairwiseDistanceEstimators[k][l].distance() };
					sum += deltaScore(pairwiseDist, options);
				}
			}
		}
	}
	return sum / numQuartets;
}

bool canGoLeftAll(const ExtendedBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoordsWithoutFlanks(i).first - offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool canGoRightAll(const ExtendedBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoordsWithoutFlanks(i).second + offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

std::pair<size_t, double> findPerfectFlankSize(const ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	double bestScore = 1;

	// for the pairwise distances:
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<std::vector<HammingDistanceEstimator> > est;
	est.resize(nTaxBlock - 1);
	for (size_t i = 0; i < nTaxBlock - 1; ++i) {
		est[i].resize(nTaxBlock - i - 1);
	}

	for (size_t i = 1; i <= options.maximumExtensionWidth; ++i) {
		if (bestScore == 0) {
			break;
		}
		if (directionRight) {
			if (!canGoRightAll(block, presenceChecker, nTax, i)) {
				break;
			}
		} else {
			if (!canGoLeftAll(block, presenceChecker, nTax, i)) {
				break;
			}
		}

		// Update sequences in est
		for (size_t j = 0; j < est.size(); ++j) {
			for (size_t k = 0; k < est[i].size(); ++k) {
				size_t coordA;
				if (directionRight) {
					coordA = block.getTaxonCoordsWithoutFlanks(j).second + i;
				} else {
					coordA = block.getTaxonCoordsWithoutFlanks(j).first - i;
				}
				size_t idB = j + 1 + k;
				size_t coordB;
				if (directionRight) {
					coordB = block.getTaxonCoordsWithoutFlanks(idB).second + i;
				} else {
					coordB = block.getTaxonCoordsWithoutFlanks(idB).first - i;
				}
				est[j][k].addChars(T[coordA], T[coordB]);
			}
		}

		double score = averageDeltaScore(est, nTaxBlock, options);
		if (score < bestScore) {
			bestScore = score;
			bestSize = i;
		}
	}

	return std::make_pair(bestSize, bestScore);
}

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
	} else {
		std::pair<size_t, double> bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false);
		std::pair<size_t, double> bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true);
		for (size_t i = 0; i < nTax; ++i) {
			if (block.hasTaxon(i)) {
				block.setLeftFlankSize(i, bestLeft.first);
				block.setRightFlankSize(i, bestRight.first);
			}
		}
	}
	return block;
}

ExtendedBlock nextExtendedBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	SeededBlock seededBlock = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
	// check if the extended block can still be accepted.
	if (presenceChecker.isFine(extendedBlock)) {
		return extendedBlock;
	} else {
		throw std::runtime_error("The extended block could not be accepted");
	}
}

std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	std::vector<ExtendedBlock> res;
	size_t actSAPos = 0;

	while (actSAPos < SA.size()) {
		if (actSAPos == 474898) {
			std::cout << "'well hello, nice little world' said the rabbit.\n";
		}

		SeededBlock seededBlock = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
		if (seededBlock.getSeedSize() == 0) { // no more seeded blocks found
			return res;
		}
		ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
		// check if the extended block can still be accepted.
		if (presenceChecker.isFine(extendedBlock)) {
			presenceChecker.reserveExtendedBlock(extendedBlock);
			res.push_back(extendedBlock);
		}
	}
	return res;
}

AlignedBlock nextAlignedBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	ExtendedBlock extBlock = nextExtendedBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	AlignedBlock bl(extBlock, nTax);
	if (bl.getSeedSize() != 0) {
		bl.align(T, options);
	}
	return bl;
}

std::vector<AlignedBlock> extractAlignedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	std::vector<AlignedBlock> res;
	std::vector<ExtendedBlock> extBlocks = extractExtendedBlocks(T, nTax, SA, lcp, presenceChecker, taxonCoords, options);

	for (size_t i = 0; i < extBlocks.size(); ++i) {
		if (extBlocks[i].getSeedSize() == 0) {
			break;
		}
		AlignedBlock bl(extBlocks[i], nTax);
		res.push_back(bl);
	}
#pragma omp parallel for
	for (size_t i = 0; i < res.size(); ++i) {
		res[i].align(T, options);
	}
	return res;
}
