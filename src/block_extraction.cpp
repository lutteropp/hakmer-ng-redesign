/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include "block_helper_functions.hpp"

#include <unordered_set>
#include <cmath>

#include "indexing/suffix_array_classic.hpp"

size_t posToTaxon(size_t pos, const std::vector<std::pair<size_t, size_t> >& taxonCoords, size_t concatSize, bool revComp) {
	if (pos >= concatSize/2 && revComp) {
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

bool acceptSeed(size_t actSAPos, size_t matchCount, size_t k, size_t nTax, const std::vector<size_t>& SA, PresenceChecker& presenceChecker,
		const std::vector<std::pair<size_t, size_t> >& taxonCoords, const std::string& T, const Options& options) {
	if (matchCount > nTax /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
		return false;
	}
	if (matchCount < options.minTaxaPerBlock) { // easy test for not enough taxa
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

bool allLeftSame(const SeededBlock& seededBlock, const std::string& T, size_t nTax) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char leftChar = T[seededBlock.getTaxonCoords(taxIDs[0]).first - 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).first - 1];
		if (actChar != leftChar)
			return false;
	}
	return true;
}

bool allRightSame(const SeededBlock& seededBlock, const std::string& T, size_t nTax) {
	std::vector<size_t> taxIDs = seededBlock.getTaxonIDsInBlock();
	char rightChar = T[seededBlock.getTaxonCoords(taxIDs[0]).second + 1];
	for (size_t i = 1; i < taxIDs.size(); ++i) {
		char actChar = T[seededBlock.getTaxonCoords(taxIDs[i]).second + 1];
		if (actChar != rightChar)
			return false;
	}
	return true;
}

void trivialExtension(SeededBlock& seededBlock, const std::string& T, PresenceChecker& presenceChecker, size_t nTax) {
	// first, perform trivial extension of the seeded block
	std::vector<size_t> taxonIDs = seededBlock.getTaxonIDsInBlock();
	bool canLeft = true;
	bool canRight = true;
	while (canLeft || canRight) {
		if (!canGoLeftAll(seededBlock, presenceChecker, nTax, 1)) {
			canLeft = false;
		}
		if (!canGoRightAll(seededBlock, presenceChecker, nTax, 1)) {
			canRight = false;
		}
		if (canLeft && !allLeftSame(seededBlock, T, nTax)) {
			canLeft = false;
		}
		if (canRight && !allRightSame(seededBlock, T, nTax)) {
			canRight = false;
		}
		if (canLeft) {
			seededBlock.decreaseTaxonCoordsLeft();
		}
		if (canRight) {
			seededBlock.increaseTaxonCoordsRight();
		}
	}
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

	size_t lastPos = SA.size() - 1;
	for (size_t sIdx = actSAPos; sIdx < SA.size(); ++sIdx) {
		size_t startPos = SA[sIdx];
		size_t k = options.minK;
		if ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, lcp, k);

		while (matchCount >= options.minTaxaPerBlock) {
			if (acceptSeed(sIdx, matchCount, k, nTax, SA, presenceChecker, taxonCoords, T, options)) {
				foundBlock = true;
				for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
					block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
				}
				trivialExtension(block, T, presenceChecker, nTax); // TODO: Improve this, do the trivial extension already before
				break;
			} else {
				if (k == options.maxK) { // no further extension of seed length
					break;
				}
				if (startPos + k + 1 >= T.size() || !presenceChecker.isFree(startPos + k)) { // newly added character would be already taken anyway
					break;
				}
				k++;
				matchCount = countMatches(sIdx, lcp, k);
			}
		}
		if (foundBlock) {
			lastPos = sIdx;
			break;
		}
	}
	actSAPos = lastPos + 1;
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
	size_t numQuartets = nTaxBlock * (nTaxBlock - 1) * (nTaxBlock - 2) * (nTaxBlock - 3) / 24;
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

std::pair<size_t, double> findPerfectFlankSize(ExtendedBlock& block, size_t nTax, const PresenceChecker& presenceChecker,
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
		if (options.noIndels) {
			if (directionRight) {
				block.noGapsMSA.addCharsRight(charsToAdd);
			} else {
				block.noGapsMSA.addCharsLeft(charsToAdd);
			}
		} else {
			if (directionRight) {
				block.starMSA.addCharsRight(charsToAdd);
			} else {
				block.starMSA.addCharsLeft(charsToAdd);
			}
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
		std::string seedSequence = extractTaxonSequence(seededBlock, seededBlock.getTaxonIDsInBlock()[0], T);
		if (options.noIndels) {
			block.noGapsMSA.init(block.getNTaxInBlock());
			block.noGapsMSA.setSeeds(seedSequence);
		} else {
			block.starMSA.init(block.getNTaxInBlock());
			block.starMSA.setSeeds(seedSequence);
		}

		std::pair<size_t, double> bestLeft = findPerfectFlankSize(block, nTax, presenceChecker, T, options, false);
		std::pair<size_t, double> bestRight = findPerfectFlankSize(block, nTax, presenceChecker, T, options, true);
		for (size_t i = 0; i < nTax; ++i) {
			if (block.hasTaxon(i)) {
				block.setLeftFlankSize(i, bestLeft.first);
				block.setRightFlankSize(i, bestRight.first);
			}
		}
		if (options.noIndels) {
			block.noGapsMSA.shrinkDownToLeftFlank(bestLeft.first);
			block.noGapsMSA.shrinkDownToRightFlank(bestRight.first);
		} else {
			block.starMSA.shrinkDownToLeftFlank(bestLeft.first);
			block.starMSA.shrinkDownToRightFlank(bestRight.first);
		}
	}
	return block;
}

ExtendedBlock nextExtendedBlock(size_t& actSAPos, const std::string& T, size_t nTax, const std::vector<size_t>& SA,
		const std::vector<size_t>& lcp, PresenceChecker& presenceChecker, const std::vector<std::pair<size_t, size_t> >& taxonCoords,
		const Options& options) {
	SeededBlock seededBlock = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
	if (seededBlock.getSeedSize() == 0) {
		ExtendedBlock bl(seededBlock, nTax);
		return bl;
	}
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

	double lastP = 0;

	while (actSAPos < SA.size()) {
		SeededBlock seededBlock = nextSeededBlock(actSAPos, T, nTax, SA, lcp, presenceChecker, taxonCoords, options);
		if (seededBlock.getSeedSize() == 0) { // no more seeded blocks found
			return res;
		}
		ExtendedBlock extendedBlock = extendBlock(seededBlock, T, nTax, presenceChecker, options);
		// check if the extended block can still be accepted.
		if (presenceChecker.isFine(extendedBlock)) {
			presenceChecker.reserveExtendedBlock(extendedBlock);

			if (options.verboseDebug) {
				std::cout << "Pushing back a block with alignment: \n";
				std::vector<std::string> msa;
				if (options.noIndels) {
					msa = extendedBlock.noGapsMSA.assembleMSA();
				} else {
					msa = extendedBlock.starMSA.assembleMSA();
				}
				for (size_t i = 0; i < msa.size(); ++i) {
					std::cout << msa[i] << "\n";
				}
			}

			res.push_back(extendedBlock);
		}

		double progress = (double) 100 * actSAPos / SA.size();
		if (progress > lastP + 1) {
			std::cout << progress << " %\n";
			lastP = progress;
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
		if (options.noIndels) {
			bl.setAlignment(extBlock.noGapsMSA.assembleMSA());
		} else {
			bl.setAlignment(extBlock.starMSA.assembleMSA());
		}
		// bl.alignMAFFT(T, options);
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
		if (options.noIndels) {
			res[i].setAlignment(extBlocks[i].noGapsMSA.assembleMSA());
		} else {
			res[i].setAlignment(extBlocks[i].starMSA.assembleMSA());
		}
		// res[i].alignMAFFT(T, options);
	}
	return res;
}
