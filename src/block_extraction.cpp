/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include "block_helper_functions.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"

#include "hmm_classification.hpp"

#include <unordered_set>
#include <cmath>
#include <queue>
#include <algorithm>

#include "indexing/suffix_array_classic.hpp"

size_t posToTaxon(size_t pos, const std::vector<IndexedTaxonCoords>& taxonCoords, size_t concatSize, bool revComp) {
	if (revComp && pos >= concatSize / 2) {
		pos = concatSize - pos - 1;
	}

	for (size_t i = 0; i < taxonCoords.size(); ++i) {
		if (taxonCoords[i].contains(pos)) {
			return i;
		}
	}
	return std::string::npos;
}

bool canStay(size_t pos, const std::vector<IndexedTaxonCoords>& taxonCoords, const std::vector<size_t>& wantedTaxa) {
	for (size_t i = 0; i < wantedTaxa.size(); ++i) {
		if (taxonCoords[wantedTaxa[i]].contains(pos)) {
			return true;
		}
	}
	return false;
}

// Shrink arrays... needed for quartet-based stuff.
std::pair<std::vector<size_t>, std::vector<size_t> > shrinkArrays(const IndexedConcatenatedSequence& concat,
		const std::vector<IndexedTaxonCoords>& taxonCoords, const std::vector<size_t>& wantedTaxa, const Options& options) {
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

/**
 * Compute the entropy of a given k-mer.
 * @param S The entire sequence data.
 * @param startPos Index of the first base in the k-mer.
 * @param k Size of the k-mer.
 * @param NCount Variable to store the number of 'N' bases in the k-mer.
 */
float info_content(const std::string& seed, unsigned int *NCount) {
	char histo[64] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	float x = 0.0;
	*NCount = 0;
	if (seed.size() == 3) {
		return 0.0; // degenerate case
	}

	size_t p = 0, q = 0;

	for (unsigned int i = 0; i <= seed.size() - 3; i++) {
		unsigned int I = 0;
		q = p;
		for (int j = 0; j < 3; j++) {
			I = (I << 2);
			switch (seed[q]) {
			case 'A':
				break;
			case 'C':
				I |= 1;
				break;
			case 'G':
				I |= 2;
				break;
			case 'T':
				I |= 3;
				break;
			case 'U': // Added by Sarah
				I |= 3;
				break;
			case 'N':
				++*NCount;
				break; // implicitly treating an N as an A...
			default:
				break; // implicitly treating all other chars as A also... :)
			}
			++q;
		}
		++p;
		++histo[I];
	}
	for (unsigned int i = 0; i < 64; i++) {
		x += histo[i] * (histo[i] - 1) / 2.0;
	}
	x /= (seed.size() - 3);
	return x;
}

/**
 * Checks whether a k-mer is of high enough complexity. For example, k-mers like AAAAAAAAAAA are bad and will fail the test.
 * @param S The input sequence data.
 * @param Position of the first base of the k-mer in the sequence data.
 * @param kmerSize Size of the k-mer.
 * @param Ncutoff Maximum number of 'N' bases allowed in a k-mer.
 * @param lowComplexityCutoff The k-mer will only be accepted if its complexity is strictly greater than this value.
 * @param gLowComplexityKmers Variable to keep track of the number of rejected k-mers due to low complexity.
 */
bool highComplexity(const std::string& seed, unsigned int Ncutoff, double lowComplexityCutoff) {
	unsigned int NCount;
	float lowComplexity = info_content(seed, &NCount);

	if (NCount > Ncutoff || lowComplexity > lowComplexityCutoff) {
		return false;
	} else {
		return true;
	}
}

bool acceptSeed(size_t actSAPos, size_t matchCount, const std::vector<std::pair<size_t, size_t> >& extraOccs, size_t k, size_t nTax,
		const std::vector<size_t>& SA, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const std::string& T, const Options& options) {
	if (matchCount + extraOccs.size() > nTax /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
		return false;
	}
	if (matchCount + extraOccs.size() < options.minTaxaPerBlock) { // easy test for not enough taxa
		return false;
	}

	if (!options.lowComplexity) {
		std::string seed = T.substr(SA[actSAPos], k);
		if (!highComplexity(seed, k / 8, 1.0)) {
			//std::cout << "rejecting seed complexity: " << seed << "\n";
			return false;
		}
		//std::cout << "accepting seed complexity: " << seed << "\n";
	}

	size_t concatSize = T.size();

	// more complicated check for paralogy
	std::unordered_set<size_t> takenTaxa;
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		size_t taxID = posToTaxon(SA[i], taxonCoords, concatSize, options.reverseComplement);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}
	for (size_t i = 0; i < extraOccs.size(); ++i) {
		size_t taxID = posToTaxon(extraOccs[i].first, taxonCoords, concatSize, options.reverseComplement);
		if (takenTaxa.find(taxID) != takenTaxa.end()) { // multiple match in taxon
			return false;
		}
		takenTaxa.insert(taxID);
	}

	size_t flankOffset = 0;
	if (!options.dynamicFlanks && options.fixedFlanks) {
		flankOffset = options.flankWidth;
	}

	// check for presence/absence of the whole region, including flanks
	for (size_t i = actSAPos; i < actSAPos + matchCount; ++i) {
		if (flankOffset > SA[i] || SA[i] + k - 1 + flankOffset >= concatSize) {
			return false;
		}
		if (!presenceChecker.isFree(SA[i] - flankOffset, SA[i] + k - 1 + flankOffset)) {
			return false;
		}
	}

	for (size_t i = 0; i < extraOccs.size(); ++i) {
		if (flankOffset > extraOccs[i].first || extraOccs[i].second + flankOffset >= concatSize) {
			return false;
		}
		if (!presenceChecker.isFree(extraOccs[i].first - flankOffset, extraOccs[i].second + flankOffset)) {
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

bool canGoLeftAll(const SeededBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoords(i).first - offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool canGoRightAll(const SeededBlock& block, const PresenceChecker& presenceChecker, size_t nTax, size_t offset) {
	bool canGo = true;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			if (!presenceChecker.isFree(block.getTaxonCoords(i).second + offset)) {
				canGo = false;
				break;
			}
		}
	}
	return canGo;
}

bool allLeftSame(const SeededBlock& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char leftChar = T[seededBlock.getTaxonCoords(taxIDs[0]).first - offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).first - offset];
		if (!ambiguousMatch(actChar, leftChar))
			return false;
	}
	return true;
}

bool allRightSame(const SeededBlock& seededBlock, const std::string& T, size_t nTax, size_t offset = 1) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char rightChar = T[seededBlock.getTaxonCoords(taxIDs[0]).second + offset];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).second + offset];
		if (!ambiguousMatch(actChar, rightChar))
			return false;
	}
	return true;
}

std::pair<size_t, size_t> computeBestCaseMaxSizes(SeededBlock& seededBlock, const std::string& T, PresenceChecker& presenceChecker,
		size_t nTax) {
	size_t maxSizeLeft = 0;
	size_t maxSizeRight = 0;
	bool canLeft = true;
	bool canRight = true;
	size_t leftOffset = 1;
	size_t rightOffset = 1;
	while (canLeft || canRight) {
		if (!canGoLeftAll(seededBlock, presenceChecker, nTax, leftOffset)) {
			canLeft = false;
		}
		if (!canGoRightAll(seededBlock, presenceChecker, nTax, rightOffset)) {
			canRight = false;
		}
		if (canLeft && !allLeftSame(seededBlock, T, nTax, leftOffset)) {
			canLeft = false;
		}
		if (canRight && !allRightSame(seededBlock, T, nTax, rightOffset)) {
			canRight = false;
		}
		if (canLeft) {
			maxSizeLeft++;
			leftOffset++;
		}
		if (canRight) {
			maxSizeRight++;
			rightOffset++;
		}
	}
	return std::make_pair(maxSizeLeft, maxSizeRight);
}

void trivialExtension(SeededBlock& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	// first, perform trivial extension of the seeded block
	std::pair<size_t, size_t> bestCaseMaxSize = computeBestCaseMaxSizes(seededBlock, T, presenceChecker, nTax);
	for (size_t i = 1; i <= bestCaseMaxSize.first; ++i) {
		if (!canGoLeftAll(seededBlock, presenceChecker, nTax, 1)) {
			break;
		}
		seededBlock.decreaseTaxonCoordsLeft();
	}

	for (size_t i = 1; i <= bestCaseMaxSize.second; ++i) {
		if (!canGoRightAll(seededBlock, presenceChecker, nTax, 1)) {
			break;
		}
		seededBlock.increaseTaxonCoordsRight();
	}
}

// TODO: Re-add mismatches and indels in seeds
std::vector<SeededBlock> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const Options& options) {
	std::vector<SeededBlock> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < SA.size(); ++sIdx) {
		size_t startPos = SA[sIdx];
		size_t k = options.minK;
		if ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, lcp, k);
		std::vector<std::pair<size_t, size_t> > extraOccs;

		while (matchCount + extraOccs.size() >= options.minTaxaPerBlock) {
			//we can stop if at least one of the exact matches we've found occurs before the current match in the suffix array... if we don't do that, we don't find all matches!
			bool stopEarly = false;
			if (lcp[sIdx] >= k) {
				stopEarly = true;
			}

			if (!stopEarly && acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options)) {
				if (options.largeSeeds) {
					// try to find largest seed size that still gets accepted
					size_t bestK = k;
					while (acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options)) {
						bestK = k;
						if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
							break;
						}
						k++;
						matchCount = countMatches(sIdx, lcp, k);
					}
					k = bestK;
					matchCount = countMatches(sIdx, lcp, k);
				}

				if (options.maxMismatches > 0) {
					std::string pattern = T.substr(startPos, k);
					extraOccs = approxMatcher.findOccurrences(T, SA, presenceChecker, pattern, options.maxMismatches, 1, false);

					// TODO: Maybe only add those approximate matches that don't collide with the exact matches we already have?

					if (!acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options)) {
						if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
							break;
						}
						k++;
						matchCount = countMatches(sIdx, lcp, k);
						continue;
					}
				}

				SeededBlock block(nTax);
				for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
					block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
				}
				for (size_t i = 0; i < extraOccs.size(); ++i) {
					block.addTaxon(posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement), extraOccs[i].first,
							extraOccs[i].second);
				}
				computeBestCaseMaxSizes(block, T, presenceChecker, nTax);
#pragma omp critical
				{
					res.push_back(block);
					if (options.verboseDebug) {
						std::cout << "Pushing back a seeded block with " << block.getNTaxInBlock() << " taxa and seed size "
								<< block.getSeedSize() << "\n";
					}
				}
				if (!options.quartetFlavor) {
					double progress = (double) 100 * sIdx / SA.size(); // TODO: Fix this, this looks kinda wrong in parallel mode
					if (progress > lastP + 1) {
#pragma omp critical
						{
							if (progress > lastP + 1) {
								std::cout << progress << " %\n";
								lastP = progress;
							}
						}
					}
				}
				break;
			} else {
				if (k == options.maxK || startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // no further extension of seed length, or newly added character would be already taken anyway
					break;
				}
				k++;
				matchCount = countMatches(sIdx, lcp, k);
			}
		}
	}
	return res;
}

void swap(std::vector<double>& vec, size_t i, size_t j) {
	double tmp = vec[i];
	vec[i] = vec[j];
	vec[j] = tmp;
}

double deltaScore(const std::array<double, 6>& pairwiseDist) {
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

double averageDeltaScore(ExtendedBlock& block, size_t nTaxBlock, const Options& options) {
	double sum = 0.0;
	size_t numQuartets;

	if (options.quickDelta) {
		numQuartets = nTaxBlock - 3;
		for (size_t i = 0; i < nTaxBlock - 3; ++i) {
			size_t j = i + 1;
			size_t k = i + 2;
			size_t l = i + 3;
			std::array<double, 6> pairwiseDist = { block.getPairwiseNormalizedDistance(i, j, options), block.getPairwiseNormalizedDistance(
					i, k, options), block.getPairwiseNormalizedDistance(i, l, options), block.getPairwiseNormalizedDistance(j, k, options),
					block.getPairwiseNormalizedDistance(j, l, options), block.getPairwiseNormalizedDistance(k, l, options) };
			sum += deltaScore(pairwiseDist);
		}
	} else {
		numQuartets = nTaxBlock * (nTaxBlock - 1) * (nTaxBlock - 2) * (nTaxBlock - 3) / 24;
		for (size_t i = 0; i < nTaxBlock - 3; ++i) {
			for (size_t j = i + 1; j < nTaxBlock - 2; ++j) {
				for (size_t k = j + 1; k < nTaxBlock - 1; ++k) {
					for (size_t l = k + 1; l < nTaxBlock; ++l) {
						std::array<double, 6> pairwiseDist = { block.getPairwiseNormalizedDistance(i, j, options),
								block.getPairwiseNormalizedDistance(i, k, options), block.getPairwiseNormalizedDistance(i, l, options),
								block.getPairwiseNormalizedDistance(j, k, options), block.getPairwiseNormalizedDistance(j, l, options),
								block.getPairwiseNormalizedDistance(k, l, options) };
						sum += deltaScore(pairwiseDist);
					}
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

std::pair<size_t, double> findPerfectFlankSizeStatic(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	double bestScore = 1.0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	size_t finalFlankSize = 0;
	for (size_t i = 1; i <= options.flankWidth; ++i) {
		if (directionRight) {
			if (!canGoRightAll(block, presenceChecker, nTax, i)) {
				break;
			}
		} else {
			if (!canGoLeftAll(block, presenceChecker, nTax, i)) {
				break;
			}
		}
		std::vector<char> charsToAdd(taxIDs.size());
		for (size_t j = 0; j < taxIDs.size(); ++j) {
			size_t coord;
			if (directionRight) {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
			} else {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
			}
			if (coord >= T.size()) {
				throw std::runtime_error("This should not happen! Coord is too large.");
			}
			charsToAdd[j] = T[coord];
		}

		if (directionRight) {
			block.msaWrapper.addCharsRight(charsToAdd);
		} else {
			block.msaWrapper.addCharsLeft(charsToAdd);
		}
		finalFlankSize = i;
	}
	double score = averageDeltaScore(block, nTaxBlock, options);
	return std::make_pair(finalFlankSize, score);
}

std::pair<size_t, double> findPerfectFlankSizeDelta(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	double bestScore = 1.0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();

	for (size_t i = 1; i <= options.maximumExtensionWidth; ++i) {
		if (bestScore <= options.maxDelta) {
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

		std::vector<char> charsToAdd(taxIDs.size());
		for (size_t j = 0; j < taxIDs.size(); ++j) {
			size_t coord;
			if (directionRight) {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
			} else {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
			}
			if (coord >= T.size()) {
				throw std::runtime_error("This should not happen! Coord is too large.");
			}
			charsToAdd[j] = T[coord];
		}

		if (directionRight) {
			block.msaWrapper.addCharsRight(charsToAdd);
		} else {
			block.msaWrapper.addCharsLeft(charsToAdd);
		}

		double score = averageDeltaScore(block, nTaxBlock, options);
		if (score < bestScore) {
			bestScore = score;
			bestSize = i;
			//std::cout << "found a better score: " << score << "\n" << " with flank size: " << i << "\n";
		} else {
			if (score > bestScore && i - bestSize >= options.earlyStopCount) {
				break;
			}
		}
	}

	return std::make_pair(bestSize, bestScore);
}

std::pair<size_t, double> findPerfectFlankSizeHMM(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	size_t bestSize = 0;
	size_t nTaxBlock = block.getNTaxInBlock();
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();

	for (size_t i = 1; i <= options.maximumExtensionWidth; ++i) {
		if (directionRight) {
			if (!canGoRightAll(block, presenceChecker, nTax, i)) {
				break;
			}
		} else {
			if (!canGoLeftAll(block, presenceChecker, nTax, i)) {
				break;
			}
		}
		std::vector<char> charsToAdd(taxIDs.size());
		for (size_t j = 0; j < taxIDs.size(); ++j) {
			size_t coord;
			if (directionRight) {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
			} else {
				coord = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
			}
			if (coord >= T.size()) {
				throw std::runtime_error("This should not happen! Coord is too large.");
			}
			charsToAdd[j] = T[coord];
		}
		if (directionRight) {
			block.msaWrapper.addCharsRight(charsToAdd);
		} else {
			block.msaWrapper.addCharsLeft(charsToAdd);
		}
	}

	std::pair<size_t, size_t> seedCoords; // TODO: This currently only works with exactly matching seeds I think...
	std::string seed;
	for (size_t i = 0; i < nTax; ++i) {
		if (block.hasTaxon(i)) {
			seed = T.substr(block.getTaxonCoordsWithoutFlanks(i).first, block.getSeedSize());
			seedCoords.first = block.msaWrapper.assembleMSA()[0].find(seed);
			seedCoords.second = seedCoords.first + block.getSeedSize() - 1;
			break;
		}
	}

	Params hmmParams = prepareHmmParams(T, options);
	bestSize = findNumGoodSitesMSA(block.msaWrapper, directionRight, seedCoords, hmmParams);

	return std::make_pair(bestSize, 1);
}

std::pair<size_t, double> findPerfectFlankSizeBigGaps(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	throw std::runtime_error("Not implemented yet");
}

std::pair<size_t, double> findPerfectFlankSize(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
		const std::string& T, const Options& options, bool directionRight) {
	if (!options.dynamicFlanks) {
		return findPerfectFlankSizeStatic(block, nTax, presenceChecker, T, options, directionRight);
	} else if (!options.useHMM) {
		return findPerfectFlankSizeDelta(block, nTax, presenceChecker, T, options, directionRight);
	} else if (!options.useBigGapsCriterion) {
		return findPerfectFlankSizeHMM(block, nTax, presenceChecker, T, options, directionRight);
	} else {
		return findPerfectFlankSizeBigGaps(block, nTax, presenceChecker, T, options, directionRight);
	}
}

ExtendedBlock extendBlock(const SeededBlock& seededBlock, const std::string& T, size_t nTax, PresenceChecker& presenceChecker,
		const Options& options) {
	ExtendedBlock block(seededBlock, nTax, options.noIndels);

	std::string seedSequence = extractTaxonSequence(seededBlock, seededBlock.getTaxonIDsInBlock()[0], T);
	block.msaWrapper.init(block.getNTaxInBlock());
	block.msaWrapper.setSeeds(seedSequence);

	std::pair<size_t, double> bestLeft;
	std::pair<size_t, double> bestRight;

	if (!options.dynamicFlanks && options.fixedFlanks) {
		bestLeft.first = options.flankWidth;
		bestRight.first = options.flankWidth;
		// we still need to add the chars to the MSA
		std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
		for (size_t i = 1; i <= options.flankWidth; ++i) {
			std::vector<char> charsToAddLeft(taxIDs.size());
			std::vector<char> charsToAddRight(taxIDs.size());
			for (size_t j = 0; j < taxIDs.size(); ++j) {
				size_t coordRight = 0;
				size_t coordLeft = 0;
				coordRight = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).second + i;
				coordLeft = block.getTaxonCoordsWithoutFlanks(taxIDs[j]).first - i;
				if (coordRight >= T.size() || coordLeft >= T.size()) {
					throw std::runtime_error("This should not happen! Coord is too large.");
				}
				charsToAddLeft[j] = T[coordLeft];
				charsToAddRight[j] = T[coordRight];
			}
			block.msaWrapper.addCharsRight(charsToAddRight);
			block.msaWrapper.addCharsLeft(charsToAddLeft);
		}
	} else {
		bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false);
		bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true);
	}

	block.setLeftFlankSize(bestLeft.first);
	block.setRightFlankSize(bestRight.first);

	block.msaWrapper.shrinkDownToLeftFlank(bestLeft.first);
	block.msaWrapper.shrinkDownToRightFlank(bestRight.first);

	block.msaWrapper.clearMSADataStructures();

	return block;
}

std::vector<SeededBlock> filterSeededBlocks(const std::vector<SeededBlock>& seededBlocks, PresenceChecker& presenceChecker,
		const Options& options) {
	std::vector<SeededBlock> res;
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		if (presenceChecker.isFine(seededBlocks[i])) {
			res.push_back(seededBlocks[i]);
			presenceChecker.reserveSeededBlock(seededBlocks[i]);
		}
	}
	return res;
}

std::vector<ExtendedBlock> extractExtendedBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const Options& options) {
	std::vector<ExtendedBlock> res;
	if (!options.quartetFlavor)
		std::cout << "Extracting seeded blocks...\n";
	std::vector<SeededBlock> seededBlocks = extractSeededBlocks(T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	std::sort(seededBlocks.begin(), seededBlocks.end(), std::greater<SeededBlock>());

	if (options.preselectSeeds) {
		seededBlocks = filterSeededBlocks(seededBlocks, presenceChecker, options);
	}

	if (!options.quartetFlavor)
		std::cout << "Assembling extended blocks...\n";
	double lastP = 0;
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		SeededBlock seededBlock = seededBlocks[i];
		if (!options.preselectSeeds && !presenceChecker.isFine(seededBlock))
			continue;
		trivialExtension(seededBlock, T, presenceChecker, nTax);
		ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
		// check if the extended block can still be accepted.
		if ((!options.preselectSeeds && presenceChecker.isFine(extendedBlock))
				|| (options.preselectSeeds && presenceChecker.isFineWithoutSeed(extendedBlock))) {
			presenceChecker.reserveExtendedBlock(extendedBlock);
			std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
			if (options.verboseDebug) {
				std::cout << "Pushing back a block with alignment: \n";
				for (size_t i = 0; i < msa.size(); ++i) {
					std::cout << msa[i] << "\n";
				}
			}
			res.push_back(extendedBlock);
		}
		if (!options.quartetFlavor) {
			double progress = (double) 100 * i / seededBlocks.size();
			if (progress > lastP + 1) {
				std::cout << progress << " %\n";
				lastP = progress;
			}
		}
	}
	return res;
}
