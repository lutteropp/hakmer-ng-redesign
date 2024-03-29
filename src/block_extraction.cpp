/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "alignment/simple_coords.hpp"
#include "alignment/simple_msa.hpp"
#include "block_extension.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"

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

bool acceptSeedComplexity(size_t actSAPos, size_t k, const IndexedConcatenatedSequence& concat, const Options& options) {
	if (!options.lowComplexity) {
		std::string seed = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[actSAPos], k);
		if (!highComplexity(seed, k / 8, 1.0)) {
			return false;
		}
	}
	return true;
}

void trimSeededBlockExtra(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	// Further trimming: Prune seed nucleotide occurrences which are kept in too few taxa.
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	if (taxIDs.size() < options.minTaxaPerBlock) {
		return;
	}
	std::vector<unsigned int> nFine(
			block.getSeedSize(taxIDs[0]) + block.getSeedCoords(taxIDs[0]).leftGapSize + block.getSeedCoords(taxIDs[0]).rightGapSize, 0);
	for (size_t i = 0; i < nFine.size(); ++i) {
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize <= i && block.getSeedCoords(tID).leftGapSize + block.getSeedSize(tID) >= i) {
				nFine[i]++;
			}
		}
	}
	// left side...
	for (size_t i = 0; i < nFine.size(); ++i) {
		if (nFine[i] >= options.minSeedTaxInColumn) {
			break;
		}
		// now we must remove the leftmost characters
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize <= i) {
				block.increaseTaxonCoordLeft(tID);
				block.addGapLeft(tID);
			}
		}
	}
	// right side...
	for (size_t i = 0; i < nFine.size(); ++i) {
		if (nFine[nFine.size() - i - 1] >= options.minSeedTaxInColumn) {
			break;
		}
		// now we must remove the rightmost characters
		for (size_t tID : taxIDs) {
			if (block.getSeedCoords(tID).leftGapSize + block.getSeedSize(tID) <= nFine.size() - i - 1) {
				block.decreaseTaxonCoordRight(tID);
				block.addGapRight(tID);
			}
		}
	}
	for (size_t tID : taxIDs) {
		if (!block.hasTaxon(tID)) {
			block.removeTaxon(tID);
		}
	}

	// Remove taxa with not enough kept seed sequence data from the seed
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		if (block.hasTaxon(taxIDs[i]) && block.getSeedCoords(taxIDs[i]).size() < options.minSeedSitesKept) {
			block.removeTaxon(taxIDs[i]);
		}
	}

	// reduce gap sizes to reasonable amount
	// TODO: This, of course, changes the k value
	taxIDs = block.getTaxonIDsInBlock();

	if (taxIDs.empty()) {
		return;
	}

	size_t minLeftGapSize = std::numeric_limits<size_t>::max();
	size_t minRightGapSize = std::numeric_limits<size_t>::max();
	for (size_t tID : taxIDs) {
		minLeftGapSize = std::min(minLeftGapSize, block.getSeedCoords(tID).leftGapSize);
		minRightGapSize = std::min(minRightGapSize, block.getSeedCoords(tID).rightGapSize);
	}
	while (minLeftGapSize > 0) {
		for (size_t tID : taxIDs) {
			block.removeGapLeft(tID);
		}
		minLeftGapSize--;
	}
	while (minRightGapSize > 0) {
		for (size_t tID : taxIDs) {
			block.removeGapRight(tID);
		}
		minRightGapSize--;
	}
}

void trimSeededBlock(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
	// prune positions which are overlapping with other, already taken blocks
	for (size_t i = 0; i < taxIDs.size(); ++i) {
		size_t tID = taxIDs[i];
		while (block.hasTaxon(tID)) {
			size_t leftCoord = block.getSeedCoords(tID).first;
			if (!presenceChecker.isFree(leftCoord)) {
				block.increaseTaxonCoordLeft(tID);
				block.addGapLeft(tID);
			} else {
				break;
			}
		}
		while (block.hasTaxon(tID)) {
			size_t rightCoord = block.getSeedCoords(tID).second;
			if (!presenceChecker.isFree(rightCoord)) {
				block.decreaseTaxonCoordRight(tID);
				block.addGapRight(tID);
			} else {
				break;
			}
		}
		if (!block.hasTaxon(tID)) {
			block.removeTaxon(tID);
		}
	}
	assert(block.getTaxonIDsInBlock().size() == block.getNTaxInBlock());
}

ExtendedBlock trimmedExtendedBlock(const Seed& oldSeed, const std::vector<std::pair<size_t, size_t>>& extraOccs,
		const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, const Options& options) {
	// ensure that the old seed does not contain a $ sign
	for (size_t i = 0; i < oldSeed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = oldSeed.getSeedCoords(oldSeed.getTaxonIDsInBlock()[i]).first;
				j <= oldSeed.getSeedCoords(oldSeed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the old seed before trimmings!");
			}
		}
	}

	// make a new seed out of the extended block.
	std::vector<size_t> taxonUsage(concat.nTax(), 0);
	Seed seed(concat.nTax());

	size_t k = oldSeed.getOriginalK();
	for (size_t tID : oldSeed.getTaxonIDsInBlock()) {
		size_t saPos = oldSeed.getSAPositions()[tID];
		size_t firstCoord = concat.getSuffixArray()[saPos];
		size_t lastCoord = firstCoord + k - 1;
		seed.addTaxon(saPos, tID, firstCoord, lastCoord);
		taxonUsage[tID]++;
	}
	for (size_t i = 0; i < extraOccs.size(); ++i) {
		size_t taxID = concat.posToTaxon(extraOccs[i].first);
		if (!seed.hasTaxon(taxID)) {
			seed.addTaxon(std::numeric_limits<size_t>::max(), taxID, extraOccs[i].first, extraOccs[i].second);
		}
		taxonUsage[taxID]++;
	}
	if (options.discardParalogMismatches) {
		for (size_t i = 0; i < concat.nTax(); ++i) {
			if (taxonUsage[i] > 1) {
				seed.removeTaxon(i);
			}
		}
	}

	// ensure that the new seed does not contain a $ sign
	for (size_t i = 0; i < seed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).first;
				j <= seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the new seed before trimmings!");
			}
		}
	}

	trimSeededBlock(seed, presenceChecker, options);

	// ensure that the new seed does not contain a $ sign
	for (size_t i = 0; i < seed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).first;
				j <= seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the seed after trimming 1!");
			}
		}
	}

	trimSeededBlockExtra(seed, presenceChecker, options);

	// ensure that the new seed does not contain a $ sign
	for (size_t i = 0; i < seed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).first;
				j <= seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the seed after trimming 2!");
			}
		}
	}

	trivialExtensionPartial(seed, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);

	// ensure that the new seed does not contain a $ sign
	for (size_t i = 0; i < seed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).first;
				j <= seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the seed after trivial extension!");
			}
		}
	}

	ExtendedBlock res = extendBlock(seed, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, k, options);

	// ensure that the extended block does not contain a $ sign
	for (size_t i = 0; i < res.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = res.getTaxonCoordsWithFlanks(res.getTaxonIDsInBlock()[i]).first;
				j <= res.getTaxonCoordsWithFlanks(res.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign after computing extended block!");
			}
		}
	}

	return res;
}

double jukesCantorCorrection(double dist) {
	// TODO: Jukes Cantor Correction doesn't work if dist >= 0.75. In this case, it will return infinity.
	return -0.75 * log(1 - (4.0 / 3) * dist);
}

double averageSubstitutionRate(const std::vector<std::string>& msa) {
	std::vector<size_t> validSeqIndices;
	for (size_t i = 0; i < msa.size(); ++i) {
		for (size_t j = 0; j < msa[i].size(); ++j) {
			if (std::find(validSeqIndices.begin(), validSeqIndices.end(), i) != validSeqIndices.end()) {
				break;
			}
			if (msa[i][j] != '?' && msa[i][j] != '-') {
				validSeqIndices.push_back(i);
				break;
			}
		}
	}
	//size_t nPairs = msa[0].size() * (validSeqIndices.size() * (validSeqIndices.size() - 1) / 2);
	size_t nPairs = 0;
	size_t nDiff = 0;
	for (size_t i = 0; i < msa[0].size(); ++i) {
		for (size_t j = 0; j < validSeqIndices.size(); ++j) {
			for (size_t k = j + 1; k < validSeqIndices.size(); ++k) {
				/*if (msa[validSeqIndices[j]][i] != '-' && msa[validSeqIndices[k]][i] != '-' && msa[validSeqIndices[j]][i] != '?'
				 && msa[validSeqIndices[k]][i] != '?' && !ambiguousMatch(msa[validSeqIndices[j]][i], msa[validSeqIndices[k]][i])) {
				 nDiff++;
				 }*/
				if (msa[validSeqIndices[j]][i] != '-' && msa[validSeqIndices[k]][i] != '-' && msa[validSeqIndices[j]][i] != '?'
						&& msa[validSeqIndices[k]][i] != '?') {
					nPairs++;
					if (!ambiguousMatch(msa[validSeqIndices[j]][i], msa[validSeqIndices[k]][i])) {
						nDiff++;
					}
				}
			}
		}
	}
	double subRate = (double) nDiff / nPairs;
	return jukesCantorCorrection(subRate);
}

size_t processExtendedBlockBuffer(std::vector<ExtendedBlock>& extendedBlockBuffer, const Options& options, SummaryStatistics& stats,
		BlockWriter& writer, const IndexedConcatenatedSequence& concat, ApproximateMatcher& approxMatcher,
		PresenceChecker& presenceChecker) {
	size_t newlyAddedBlocks = 0;
#pragma omp parallel for reduction(+:newlyAddedBlocks)
	for (size_t i = 0; i < extendedBlockBuffer.size(); ++i) {
		ExtendedBlock block = extendedBlockBuffer[i];
		std::vector<std::pair<size_t, size_t> > extraOccs;

		if (block.getNTaxInBlock() < concat.nTax()) {
			//double subRate = std::min(block.getMySeededBlock().getSubRate(), options.maxErrorRate);

			double subRate = options.maxErrorRate;

//#pragma omp critical
			//		std::cout << "subRate: " << subRate << "\n";
			size_t k = block.getMySeededBlock().getOriginalK();
			size_t maxMismatches = std::ceil(k * subRate);
			maxMismatches = std::min(maxMismatches, options.maxMismatches);

			if (maxMismatches > 0) { // augment the block with approximate matches
				std::string pattern = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[block.getMySeededBlock().getFirstSAPos()],
						k);
				std::vector<size_t> taxPresence(concat.nTax(), 0);
				for (size_t tID : block.getTaxonIDsInBlock()) {
					taxPresence[tID] += 2; // in order not to search for approx matches in the taxa we already have with exact matches
				}
				for (size_t tID : block.getBlockedTaxa()) {
					taxPresence[tID] += 2;
				}
				// do it on a per-taxon basis

				for (size_t taxonID = 0; taxonID < concat.nTax(); taxonID++) {
					if (taxPresence[taxonID] >= 2) {
						continue;
					}
					std::vector<std::pair<size_t, size_t> > eOccs = approxMatcher.findFewOccurrences(concat.getConcatenatedSeq(),
							concat.getSuffixArray(taxonID), concat.getLcpArray(taxonID), presenceChecker, pattern, maxMismatches, 1, false,
							concat, taxPresence);
					for (size_t eIdx = 0; eIdx < eOccs.size(); ++eIdx) {
						extraOccs.push_back(eOccs[eIdx]);
					}
				}
			}
		}

		if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
#pragma omp critical
			{
				presenceChecker.freeExtendedBlock(block); // TODO: Here is the bug!!! Because we also free the $ positions!!!
				block = trimmedExtendedBlock(block.getMySeededBlock(), extraOccs, concat, presenceChecker, options);
				if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
					presenceChecker.reserveExtendedBlock(block);
				}
			}

			if (block.getNTaxInBlock() >= options.minTaxaPerBlock) {
				std::vector<std::string> msa = computeMSA(block, concat.getConcatenatedSeq(), concat.nTax(), options);
#pragma omp critical
				{
					newlyAddedBlocks++;
					stats.updateSummaryStatistics(block, concat.nTax());
					if (!options.outpath.empty()) {
						try {
							writer.writeTemporaryBlockMSA(msa, concat.nTax());
						} catch (std::runtime_error& e) {
							//cout << e.what() << "\n";
							throw e;
						}
					}
				}
			}
		}
	}
	return newlyAddedBlocks;
}

std::vector<ExtendedBlock> processSeedBuffer(std::vector<Seed>& seedInfoBuffer, const IndexedConcatenatedSequence& concat,
		PresenceChecker& presenceChecker, const Options& options, ApproximateMatcher& approxMatcher) {
	std::vector<ExtendedBlock> extendedBlocks;

#pragma omp parallel for
	for (size_t i = 0; i < seedInfoBuffer.size(); ++i) {
		Seed block = seedInfoBuffer[i];

		// first check if the current seed candidate is hopeless
		trimSeededBlock(block, presenceChecker, options);
		if (block.getNTaxInBlock() < options.minTaxaPerBlock) {
			continue;
		}

		size_t k = block.getSeedCoords(block.getTaxonIDsInBlock()[0]).size()
				+ block.getSeedCoords(block.getTaxonIDsInBlock()[0]).leftGapSize
				+ block.getSeedCoords(block.getTaxonIDsInBlock()[0]).rightGapSize;

		// check some little assertion...
		for (size_t tID : block.getTaxonIDsInBlock()) {
			assert(k == block.getSeedCoords(tID).size() + block.getSeedCoords(tID).leftGapSize + block.getSeedCoords(tID).rightGapSize);
		}

		// Final trimming and adding
		size_t oldSeedSize = block.getAverageSeedSize();
		trimSeededBlock(block, presenceChecker, options);
		trimSeededBlockExtra(block, presenceChecker, options);

		if (block.getNTaxInBlock() >= options.minTaxaPerBlock) { // we need this extra check because in parallel mode, our taxa in the block could get invalidated all the time
			// Partially extend the seed
			trivialExtensionPartial(block, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
			ExtendedBlock extendedBlock = extendBlock(block, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, oldSeedSize,
					options);

			bool discardMe = (options.discardUninformativeBlocks && extendedBlock.getAverageLeftFlankSize() == 0
					&& extendedBlock.getAverageRightFlankSize() == 0);

			if (!discardMe) {
#pragma omp critical
				{
					presenceChecker.reserveExtendedBlock(extendedBlock);
					extendedBlocks.push_back(extendedBlock);
				}
			}
		}
	}
	return extendedBlocks;
}

void printMSA(const std::vector<std::string>& msa) {
	std::cout << "block MSA:\n";
	for (size_t i = 0; i < msa.size(); ++i) {
		std::cout << msa[i] << "\n";
	}
}

double estimateSubRate(const Seed& unextendedSeed, const IndexedConcatenatedSequence& concat, const Options& options) {
	ExtendedBlock block(unextendedSeed, concat.nTax());
	block.setLeftFlankSize(unextendedSeed.getOriginalK());
	block.setRightFlankSize(unextendedSeed.getOriginalK());
	std::vector<std::string> msa = computeMSA(block, concat.getConcatenatedSeq(), concat.nTax(), options);
	printMSA(msa);
	double subRate = averageSubstitutionRate(msa);
	return subRate;
}

// Just assume ungapped MSA for subRate estimation here
double estimateSubRateQuick(const Seed& unextendedSeed, const IndexedConcatenatedSequence& concat, const Options& options) {
	size_t k = unextendedSeed.getOriginalK();
	size_t nPairs = 0;
	size_t nDiff = 0;
	for (size_t tID : unextendedSeed.getTaxonIDsInBlock()) {
		size_t firstA = unextendedSeed.getSeedCoords(tID).first;
		size_t lastA = unextendedSeed.getSeedCoords(tID).second;
		for (size_t tID2 : unextendedSeed.getTaxonIDsInBlock()) {
			if (tID == tID2) {
				continue;
			}
			size_t firstB = unextendedSeed.getSeedCoords(tID2).first;
			size_t lastB = unextendedSeed.getSeedCoords(tID2).second;
			for (size_t i = 1; i <= k; ++i) { // left flank
				if (concat.getConcatenatedSeq()[firstA - i] == '$' || concat.getConcatenatedSeq()[firstB - i] == '$') {
					break;
				}
				if (firstA < i || firstB < i) {
					break;
				}
				nPairs++;
				if (!ambiguousMatch(concat.getConcatenatedSeq()[firstA - i], concat.getConcatenatedSeq()[firstB - i])) {
					nDiff++;
				}
			}
			for (size_t i = 1; i <= k; ++i) { // right flank
				if (concat.getConcatenatedSeq()[lastA + i] == '$' || concat.getConcatenatedSeq()[lastB + i] == '$') {
					break;
				}
				if (lastA + i >= concat.getConcatSize() || lastB + i >= concat.getConcatSize()) {
					break;
				}
				nPairs++;
				if (!ambiguousMatch(concat.getConcatenatedSeq()[lastA + i], concat.getConcatenatedSeq()[lastB + i])) {
					nDiff++;
				}
			}
		}
	}

	for (size_t i = 0; i < unextendedSeed.getNTaxInBlock(); ++i) {
		for (size_t j = i + 1; j < unextendedSeed.getNTaxInBlock(); ++j) {
			nPairs += k;
		}
	}

	//nPairs += k * unextendedSeed.getNTaxInBlock() * (unextendedSeed.getNTaxInBlock() - 1) / 2;
	double subRate = (double) nDiff / nPairs;
	subRate = jukesCantorCorrection(subRate);

	/*if (subRate > options.maxAvgSubstitutionRate) {
	 #pragma omp critical
	 {
	 double subRateNew = estimateSubRate(unextendedSeed, concat, options);
	 std::cout << "old: " << subRate << "; new: " << subRateNew << "\n";
	 subRate = subRateNew;
	 }
	 }*/

	return subRate;
}

Seed findSeed(size_t saPos, const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, size_t minK,
		const Options& options) {
	Seed emptySeed(concat.nTax());
	size_t kStart = minK; //std::max(minK, concat.getLcpArray()[saPos] + 1);
	if ((concat.getSuffixArray()[saPos] + kStart >= concat.getConcatenatedSeq().size())
			|| (!presenceChecker.isFree(concat.getSuffixArray()[saPos], concat.getSuffixArray()[saPos] + kStart - 1))) {
		return emptySeed;
	}

	std::vector<size_t> taxCounts(concat.nTax(), 0);
	// at least nMin next entries have to be included
	for (size_t i = 1; i < options.minTaxaPerBlock; ++i) {
		size_t lcpVal = concat.getLcpArray()[saPos + i];
		if (lcpVal < kStart) { // not enough occurrences found
			return emptySeed;
		}
	}

	// how much more can we go to the right without getting smaller than our kStart?
	taxCounts[concat.posToTaxon(concat.getSuffixArray()[saPos])]++;
	size_t lastIdx = saPos;
	for (size_t sIdx = saPos + 1; sIdx < concat.getSuffixArray().size(); ++sIdx) {
		if (concat.getLcpArray()[sIdx] < kStart) {
			break;
		} else {
			taxCounts[concat.posToTaxon(concat.getSuffixArray()[sIdx])]++;
			lastIdx = sIdx;
		}
	}
	size_t nGood = 0;
	for (size_t i = 0; i < concat.nTax(); ++i) {
		if (taxCounts[i] == 1) {
			nGood++;
		}
	}
	if (nGood < options.minTaxaPerBlock) { // not enough non-paralogous occurrences found
		return emptySeed;
	}

	// Now, it gets interesting... what's the largest value for k we can use for the non-paralogous taxa we still have?
	size_t k = std::numeric_limits<size_t>::max();
	size_t lastCheckedPos = saPos;
	for (size_t i = saPos + 1; i <= lastIdx; ++i) {
		size_t taxID = concat.posToTaxon(concat.getSuffixArray()[i]);
		if (taxCounts[taxID] == 1) { // I think that here was the problem... we were skipping some suffix array positions!
			for (size_t j = lastCheckedPos + 1; j <= i; ++j) {
				k = std::min(k, concat.getLcpArray()[j]);
			}
			lastCheckedPos = i;
		}
	}

	if (k < minK) {
		return emptySeed;
	}

	// reject instances where the same seed size led to some taxa showing paralog occurrences
	for (size_t i = saPos + 1; i <= lastIdx; ++i) {
		if (taxCounts[concat.posToTaxon(concat.getSuffixArray()[i])] > 1) {
			if (concat.getLcpArray()[i] >= k) {
				return emptySeed;
			}
		}
	}

	Seed seed(concat.nTax());
	// Now, we add those taxa to the seed.
	for (size_t i = saPos; i <= lastIdx; ++i) {
		size_t taxID = concat.posToTaxon(concat.getSuffixArray()[i]);
		if (taxCounts[taxID] == 1) {
			size_t firstPos = concat.getSuffixArray()[i];
			size_t lastPos = firstPos + k - 1;
			seed.addTaxon(i, taxID, firstPos, lastPos);
		} else if (taxCounts[taxID] > 1) {
			seed.blockTaxon(taxID);
		}
	}
	seed.setFirstSAPos(saPos);
	//double subRate = estimateSubRateQuick(seed, concat, options);
	double subRate = options.maxAvgSubstitutionRate;
	seed.setSubRate(subRate);

	/*
	// ensure that the newly discovered seed does not contain a $ sign
	for (size_t i = 0; i < seed.getTaxonIDsInBlock().size(); ++i) {
		for (size_t j = seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).first;
				j <= seed.getSeedCoords(seed.getTaxonIDsInBlock()[i]).second; ++j) {
			if (concat.getConcatenatedSeq()[j] == '$') {
				throw std::runtime_error("Encountered a $ sign in the extracted seed!");
			}
		}
	}
	*/

//#pragma omp critical
	//std::cout << "k: " << k << "; n: " << seed.getNTaxInBlock() << "; Estimated sub rate: " << subRate << "\n";

	if ((options.discardUninformativeBlocks && subRate == 0) || subRate > options.maxAvgSubstitutionRate) {
		return emptySeed;
	}

	return seed;
}

std::vector<Seed> extractSeeds(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, size_t minK, size_t maxK,
		const Options& options) {
	std::vector<Seed> res;
	size_t actSAPos = 0;
	double lastP = 0;

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < concat.getSuffixArray().size() - options.minTaxaPerBlock; ++sIdx) {
		bool skip = false;

		Seed seed = findSeed(sIdx, concat, presenceChecker, minK, options);
		if (seed.getNTaxInBlock() == 0) {
			skip = true;
		}

		if (!skip && acceptSeedComplexity(sIdx, seed.getAverageSeedSize(), concat, options)) {
			if (!canGoLeftAll(seed, presenceChecker, concat.nTax())
					|| !allLeftSame(seed, concat.getConcatenatedSeq(), seed.getTaxonIDsInBlock())) {
#pragma omp critical
				res.push_back(seed);
				if (options.verboseDebug) {
#pragma omp critical
					std::cout << "Pushing back a seeded occurrence with " << seed.getNTaxInBlock() << " taxa and seed size "
							<< seed.getAverageSeedSize() << "\n";
				}
			}
		}
		double progress = (double) 100 * sIdx / concat.getSuffixArray().size();
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
	return res;
}

size_t elbowMethod(const std::vector<std::pair<size_t, size_t> >& seedSizeCounts, const Options& options) {
// see https://www.linkedin.com/pulse/finding-optimal-number-clusters-k-means-through-elbow-asanka-perera/
// see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
// we need to find the point with the largest distance to the line from the first to the last point; this point corresponds to our chosen minK value.
	static size_t minReasonableCount = 1000;

	// find the entry with the highest count, this one will be at minIdx
	size_t minIdx = 0;
	for (size_t i = 1; i < seedSizeCounts.size(); ++i) {
		if (seedSizeCounts[i].second > seedSizeCounts[minIdx].second) {
			minIdx = i;
		}
	}

	size_t lastIdx = seedSizeCounts.size() - 1;
	while (seedSizeCounts[lastIdx].second < minReasonableCount && lastIdx > minIdx) {
		lastIdx--;
	}

	int maxDist = 0;
	size_t maxDistIdx = minIdx;

	int x1 = seedSizeCounts[minIdx].first;
	int y1 = seedSizeCounts[minIdx].second;
	int x2 = seedSizeCounts[lastIdx].first;
	int y2 = seedSizeCounts[lastIdx].second;
	for (size_t i = minIdx + 1; i <= lastIdx; ++i) { // because the endpoints trivially have distance 0
		int x0 = seedSizeCounts[i].first;
		int y0 = seedSizeCounts[i].second;
		int d = std::abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1);
		if (d > maxDist) {
			maxDist = d;
			maxDistIdx = i;
		}
	}
	return seedSizeCounts[maxDistIdx].first;
}

void printSeedSizeHistogram(const std::vector<std::pair<size_t, size_t> >& seedSizes) {
	std::cout << "Seed size histogram:\n";
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		std::cout << seedSizes[i].first << ": " << seedSizes[i].second << "\n";
	}
}

std::vector<std::pair<size_t, size_t> > estimateSeedSizes(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker,
		size_t minK, size_t maxK, const Options& options) {
	std::vector<std::pair<size_t, size_t> > resSizes;
	std::vector<size_t> seedSizes(1000, 0);

	size_t actSAPos = 0;
	double lastP = 0;

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < concat.getSuffixArray().size() - options.minTaxaPerBlock; ++sIdx) {
		bool skip = false;

		Seed seed = findSeed(sIdx, concat, presenceChecker, minK, options);
		if (seed.getNTaxInBlock() == 0) {
			skip = true;
		}

		if (!skip && acceptSeedComplexity(sIdx, seed.getAverageSeedSize(), concat, options)) {
			if (!canGoLeftAll(seed, presenceChecker, concat.nTax())
					|| !allLeftSame(seed, concat.getConcatenatedSeq(), seed.getTaxonIDsInBlock())) {
				size_t k = seed.getAverageSeedSize();
#pragma omp critical
				{
					if (k >= seedSizes.size()) {
						seedSizes.resize(k + 1);
					}
					seedSizes[k]++;
				}
			}
		}
		double progress = (double) 100 * sIdx / concat.getSuffixArray().size();
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

	for (size_t i = 0; i < seedSizes.size(); ++i) {
		if (seedSizes[i] > 0) {
			std::pair<size_t, size_t> p = std::make_pair(i, seedSizes[i]);
			resSizes.push_back(p);
		}
	}

	return resSizes;
}

size_t selectAndProcessSeeds(const std::vector<Seed>& seeds, ApproximateMatcher& approxMatcher, const IndexedConcatenatedSequence& concat,
		PresenceChecker& presenceChecker, BlockWriter& writer, SummaryStatistics& stats, const Options& options, size_t minK, size_t maxK) {
	size_t newlyAddedBlocks = 0;
	std::vector<Seed> buffer;
	size_t lastN = seeds[0].getNTaxInBlock();
	buffer.push_back(seeds[0]);
	for (size_t i = 1; i < seeds.size(); ++i) {
		if (seeds[i].getNTaxInBlock() == lastN) {
			if (seeds[i].getAverageSeedSize() >= minK && seeds[i].getAverageSeedSize() <= maxK) {
				buffer.push_back(seeds[i]);
			}
		} else {
			std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
			std::vector<ExtendedBlock> extendedBlockBuffer = processSeedBuffer(buffer, concat, presenceChecker, options, approxMatcher);
			std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
			newlyAddedBlocks += processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat, approxMatcher,
					presenceChecker);
			buffer.clear();
			lastN = seeds[i].getNTaxInBlock();
			if (seeds[i].getAverageSeedSize() >= minK && seeds[i].getAverageSeedSize() <= maxK) {
				buffer.push_back(seeds[i]);
			}
		}
	}
// process the last buffer
	std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
	std::vector<ExtendedBlock> extendedBlockBuffer = processSeedBuffer(buffer, concat, presenceChecker, options, approxMatcher);
	std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
	newlyAddedBlocks += processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat, approxMatcher, presenceChecker);
	std::cout << "Newly added blocks: " << newlyAddedBlocks << "\n";
	return newlyAddedBlocks;
}

void printHypotheticalBestCaseTaxonCoverage(const IndexedConcatenatedSequence& concat, const std::vector<Seed>& seeds) {
// if we would ignore overlaps and stuff like this
	std::vector<size_t> taxUsage(concat.nTax(), 0);
	for (size_t i = 0; i < seeds.size(); ++i) {
		size_t k = seeds[i].getAverageSeedSize();
		for (size_t tID : seeds[i].getTaxonIDsInBlock()) {
			taxUsage[tID] += k;
		}
	}

	std::cout << "\nHypothetical best-case taxon usages if we wouldn't care about overlaps and stuff like this:\n";
	for (size_t i = 0; i < concat.nTax(); ++i) {
		double p = ((double) taxUsage[i] * 100) / concat.getTaxonCoords(i).getTotalLength();
		std::cout << concat.getTaxonLabels()[i] << ": " << p << " %\n";
	}
	std::cout << "\n";
}

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t maxK) {
	PresenceChecker seedingPresenceChecker(presenceChecker);
	size_t initialMinK = minK;

	std::cout << "Estimating seed sizes...\n";
	std::vector<std::pair<size_t, size_t> > seedSizes = estimateSeedSizes(concat, presenceChecker, minK, maxK, options);
	printSeedSizeHistogram(seedSizes);
	size_t newMinK = elbowMethod(seedSizes, options);
	seedSizes.clear();
	seedSizes.shrink_to_fit();
	std::cout << "New chosen minK by using the elbow method: " << newMinK << ". Ignoring all seeds with smaller k than this value.\n";
	minK = newMinK;

	std::cout << "Extracting seeded blocks...\n";
	std::vector<Seed> seeds = extractSeeds(concat, seedingPresenceChecker, minK, maxK, options);
	std::cout << "seeded block infos.size(): " << seeds.size() << "\n";
	std::cout << "Processing seeds...\n";
	std::sort(seeds.begin(), seeds.end(), std::greater<Seed>());
	ApproximateMatcher approxMatcher(options.mismatchesOnly);

	size_t newlyAddedBlocks = selectAndProcessSeeds(seeds, approxMatcher, concat, presenceChecker, writer, stats, options, minK, maxK);
	double seqDataUsed = stats.getAmountSeqDataUsed(concat.getSequenceDataSize());
	std::cout << "seqDataUsed: " << seqDataUsed << "\n";

// TODO: Remove me again, this is just out of curiosity
//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
//std::cout << "We have made " << superseeds.size() << " superseeds out of " << seededBlocks.size() << " seeded blocks.\n";
}
