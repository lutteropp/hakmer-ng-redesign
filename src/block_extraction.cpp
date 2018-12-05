/*
 * block_extraction.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include "block_extraction.hpp"
#include "dna_functions.hpp"
#include "indexing/approx_matching.hpp"

#include <unordered_set>
#include <cmath>
#include <queue>
#include <algorithm>

#include "indexing/suffix_array_classic.hpp"

#include "build_superseeds.hpp"
#include "block_extension.hpp"
#include "seed_info.hpp"

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

bool acceptSeedComplexity(size_t actSAPos, size_t matchCount, size_t k, size_t nTax, const std::vector<size_t>& SA,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const std::string& T,
		const Options& options) {
	if (!options.lowComplexity) {
		std::string seed = T.substr(SA[actSAPos], k);
		if (!highComplexity(seed, k / 8, 1.0)) {
			//std::cout << "rejecting seed complexity: " << seed << "\n";
			return false;
		}
		//std::cout << "accepting seed complexity: " << seed << "\n";
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

void processExtendedBlockBuffer(std::vector<ExtendedBlock>& extendedBlockBuffer, const Options& options, SummaryStatistics& stats,
		BlockWriter& writer, const IndexedConcatenatedSequence& concat) {
#pragma omp parallel for
	for (size_t i = 0; i < extendedBlockBuffer.size(); ++i) {
		ExtendedBlock extendedBlock = extendedBlockBuffer[i];
#pragma omp critical
		stats.updateSummaryStatistics(extendedBlock, concat.nTax());
		if (!options.outpath.empty()) {
			writer.writeTemporaryBlockMSA(extendedBlock, concat.getConcatenatedSeq(), concat.nTax(), options);
		}
	}
}

std::vector<ExtendedBlock> processSeedInfoBuffer(std::vector<SeedInfo>& seedInfoBuffer, const IndexedConcatenatedSequence& concat,
		PresenceChecker& presenceChecker, const Options& options, ApproximateMatcher& approxMatcher, size_t flankWidth,
		const std::vector<uint16_t>& posTaxonArray) {
	std::vector<ExtendedBlock> extendedBlocks;

#pragma omp parallel for
	for (size_t i = 0; i < seedInfoBuffer.size(); ++i) {
		SeedInfo seedInfo = seedInfoBuffer[i];
		// make a SeededBlock out of the seedInfo
		Seed block(concat.nTax());
		size_t sIdx = seedInfo.saPos;
		size_t k = seedInfo.k;
		size_t matchCount = seedInfo.n;
		for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
			//size_t tID = posToTaxon(concat.getSuffixArray()[i], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement);
			size_t tID = posTaxonArray[i];
			block.addTaxon(tID, concat.getSuffixArray()[i], concat.getSuffixArray()[i] + k - 1);
		}
		for (std::pair<size_t, size_t> extraOcc : seedInfo.extraOccs) {
			block.addTaxon(
					posToTaxon(extraOcc.first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement),
					extraOcc.first, extraOcc.second);
		}

		// first check if the current seed candidate is hopeless
		trimSeededBlock(block, presenceChecker, options);
		if (block.getNTaxInBlock() < options.minTaxaPerBlock) {
			continue;
		}

		k = block.getSeedCoords(block.getTaxonIDsInBlock()[0]).size() + block.getSeedCoords(block.getTaxonIDsInBlock()[0]).leftGapSize
				+ block.getSeedCoords(block.getTaxonIDsInBlock()[0]).rightGapSize;

		// check some little assertion...
		for (size_t tID : block.getTaxonIDsInBlock()) {
			assert(k == block.getSeedCoords(tID).size() + block.getSeedCoords(tID).leftGapSize + block.getSeedCoords(tID).rightGapSize);
		}

		// now we decided that we'd like to take this seeded block. Augment it with mismatches, then trim it again.
		// Augment the seed with mismatches
		if (options.maxMismatches > 0 && matchCount < concat.nTax()) {
			std::unordered_set<size_t> taxIDs;
			for (size_t t = 0; t < matchCount; ++t) {
				//size_t tID = posToTaxon(concat.getSuffixArray()[sIdx + t], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement);
				size_t tID = posTaxonArray[sIdx + t];
				taxIDs.insert(tID);
			}

			std::string pattern = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[sIdx], k);
			std::vector<std::pair<size_t, size_t> > extraOccs = approxMatcher.findOccurrences(concat.getConcatenatedSeq(),
					concat.getSuffixArray(), presenceChecker, pattern, options.maxMismatches, 1, false);
			for (size_t i = 0; i < extraOccs.size(); ++i) {
				size_t taxID = posToTaxon(extraOccs[i].first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
						options.reverseComplement);
				if (taxID < concat.nTax() && taxIDs.find(taxID) == taxIDs.end()) {
					// check if the extra occurrence is fine
					if (extraOccs[i].second >= concat.getConcatenatedSeq().size()) {
						continue;
					}
					block.addTaxon(taxID, extraOccs[i].first, extraOccs[i].second);
					taxIDs.insert(taxID);
				}
			}
			trimSeededBlock(block, presenceChecker, options);

			if (block.getNTaxInBlock() < options.minTaxaPerBlock) {
				continue; // this can happen in parallel mode, because some other block could have been taken in the meantime.
			}
		}

		// TODO: If we want to extend the block already here, then we need a trimExtendedBlock function.

		// Final trimming and adding, this time in critical mode
#pragma omp critical
		{
			trimSeededBlock(block, presenceChecker, options);
			trimSeededBlockExtra(block, presenceChecker, options);

			if (block.getNTaxInBlock() >= options.minTaxaPerBlock) { // we need this extra check because in parallel mode, our taxa in the block could get invalidated all the time
				// Partially extend the seed
				trivialExtensionPartial(block, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
				ExtendedBlock extendedBlock = extendBlock(block, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, flankWidth,
						options);

				// TODO: Correctly implement extended block trimming, this is needed because of the parallelization
				//trimExtendedBlock(extendedBlock, presenceChecker, options);
				presenceChecker.reserveExtendedBlock(extendedBlock);
				extendedBlocks.push_back(extendedBlock);
			}
		}
	}
	return extendedBlocks;
}

std::pair<size_t, size_t> findKAndNumExactMatches(size_t saPos, const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker,
		size_t minK, const Options& options, const std::vector<uint16_t>& posTaxonArray) {
	std::pair<size_t, size_t> emptyPair = std::make_pair(0, 0);
	size_t kStart = std::max(minK, concat.getLcpArray()[saPos] + 1);

	if ((concat.getSuffixArray()[saPos] + kStart >= concat.getConcatenatedSeq().size())
			|| (!presenceChecker.isFree(concat.getSuffixArray()[saPos], concat.getSuffixArray()[saPos] + kStart - 1))) {
		return emptyPair;
	}

	// at least nMin next entries have to be included
	std::unordered_set<size_t> taxIDs;
	size_t minNeighborK = std::numeric_limits<size_t>::max();
	for (size_t i = 0; i < options.minTaxaPerBlock; ++i) {
		if (i > 0 && concat.getLcpArray()[saPos + i] < kStart) {
			return emptyPair;
		}
		size_t tID = posTaxonArray[saPos + i];
		//size_t tID = posToTaxon(concat.getSuffixArray()[saPos + i], concat.getTaxonCoords(), concat.getConcatSize(), options.reverseComplement);
		taxIDs.insert(tID);
		if (i > 0) {
			minNeighborK = std::min(minNeighborK, concat.getLcpArray()[saPos + i]);
		}
	}
	if (taxIDs.size() < options.minTaxaPerBlock) {
		return emptyPair;
	}
	size_t actMatches = options.minTaxaPerBlock;
	// now, we are at a position where at least the minimum amount of taxa has been found. Let's check how much more we can go to the right without encountering any paralogs.
	size_t newMinK = kStart;
	for (size_t i = options.minTaxaPerBlock; i < concat.getSuffixArray().size() - saPos; ++i) {
		if (concat.getLcpArray()[saPos + i] < kStart) {
			break;
		}
		size_t tID = posTaxonArray[saPos + i];
		//size_t tID = posToTaxon(concat.getSuffixArray()[saPos + i], concat.getTaxonCoords(), concat.getConcatSize(), options.reverseComplement);
		taxIDs.insert(tID);
		if (taxIDs.size() != actMatches + 1) {
			newMinK = concat.getLcpArray()[saPos + i] + 1;
			break;
		}
		actMatches++;
	}

	if (minNeighborK < newMinK) {
		return emptyPair;
	}
	size_t k = newMinK; // now, we just need to check how many matches we still have and what's their min value

	size_t minValNewMatches = std::numeric_limits<size_t>::max();
	size_t nMatches = 1;
	for (size_t i = 1; i < concat.getSuffixArray().size() - saPos; ++i) {
		if (concat.getLcpArray()[saPos + i] < k) {
			break;
		}
		nMatches++;
		minValNewMatches = std::min(minValNewMatches, concat.getLcpArray()[saPos + i]);
	}
	k = minValNewMatches;

	return std::make_pair(k, nMatches);
}

std::vector<SeedInfo> extractSeedInfos(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, size_t minK,
		size_t maxK, const Options& options, const std::vector<uint16_t>& posTaxonArray) {
	std::vector<SeedInfo> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < concat.getSuffixArray().size() - options.minTaxaPerBlock; ++sIdx) {
		bool skip = false;
		std::pair<size_t, size_t> p = findKAndNumExactMatches(sIdx, concat, presenceChecker, minK, options, posTaxonArray);
		size_t k = p.first;
		size_t matchCount = p.second;
		if (k == 0) {
			skip = true;
		}

		if (!skip
				&& acceptSeedComplexity(sIdx, matchCount, k, concat.nTax(), concat.getSuffixArray(), presenceChecker,
						concat.getTaxonCoords(), concat.getConcatenatedSeq(), options)) {
			// perform simple extension of the seed
			Seed block(concat.nTax());
			for (size_t idx = 0; idx < matchCount; ++idx) {
				//size_t tID = posToTaxon(concat.getSuffixArray()[sIdx + idx], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement);
				size_t tID = posTaxonArray[sIdx + idx];
				block.addTaxon(tID, concat.getSuffixArray()[sIdx + idx], concat.getSuffixArray()[sIdx + idx] + k - 1);
			}
			trivialExtensionSimple(block, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
			k = block.getAverageSeedSize();

			SeedInfo info(sIdx, k, matchCount);
			/*
			 // Augment the seed with mismatches
			 if (options.maxMismatches > 0) {
			 std::unordered_set<size_t> taxIDs;
			 for (size_t t = 0; t < matchCount; ++t) {
			 //size_t tID = posToTaxon(concat.getSuffixArray()[sIdx + t], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(), options.reverseComplement);
			 size_t tID = posTaxonArray[sIdx + t];
			 taxIDs.insert(tID);
			 }

			 std::string pattern = concat.getConcatenatedSeq().substr(concat.getSuffixArray()[sIdx], k);
			 std::vector<std::pair<size_t, size_t> > extraOccs = approxMatcher.findOccurrences(concat.getConcatenatedSeq(),
			 concat.getSuffixArray(), presenceChecker, pattern, options.maxMismatches, 1, false);
			 for (size_t i = 0; i < extraOccs.size(); ++i) {
			 size_t taxID = posToTaxon(extraOccs[i].first, concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
			 options.reverseComplement);
			 if (taxID < concat.nTax() && taxIDs.find(taxID) == taxIDs.end()) {
			 // check if the extra occurrence is fine
			 if (extraOccs[i].second >= concat.getConcatenatedSeq().size()) {
			 continue;
			 }
			 info.addExtraOcc(extraOccs[i]);
			 taxIDs.insert(taxID);
			 }
			 }
			 assert(info.n + info.extraOccs.size() <= concat.nTax());
			 }*/

#pragma omp critical
			res.push_back(info);
			if (options.verboseDebug) {
#pragma omp critical
				std::cout << "Pushing back a seeded occurrence with " << info.n + info.extraOccs.size() << " taxa and seed size " << info.k
						<< "\n";
			}
		}
		double progress = (double) 100 * sIdx / concat.getSuffixArray().size(); // TODO: Fix this, this looks kinda wrong in parallel mode
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
	// see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
	// we need to find the point with the largest distance to the line from the first to the last point; this point corresponds to our chosen minK value.
	int maxDist = 0;
	size_t maxDistIdx = 0;
	int x1 = seedSizeCounts[0].first;
	int y1 = seedSizeCounts[0].second;
	int x2 = seedSizeCounts[seedSizeCounts.size() - 1].first;
	int y2 = seedSizeCounts[seedSizeCounts.size() - 1].second;
	for (size_t i = 1; i < seedSizeCounts.size() - 1; ++i) { // because the endpoints trivially have distance 0
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

std::vector<std::pair<size_t, size_t> > countSeedSizes(const std::vector<SeedInfo>& seedInfos, const Options& options) {
	std::vector<std::pair<size_t, size_t> > res;
	std::vector<size_t> seedSizes(1000, 0);
	for (size_t i = 0; i < seedInfos.size(); ++i) {
		seedSizes[seedInfos[i].k]++;
	}
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		if (seedSizes[i] > 0) {
			res.push_back(std::make_pair(i, seedSizes[i]));
		}
	}
	return res;
}

/*size_t rechooseMinK(const std::vector<SeedInfo>& seedInfos, const Options& options) {
	size_t maxCount = 0;
	std::vector<size_t> seedSizes(500, 0);
	for (size_t i = 0; i < seedInfos.size(); ++i) {
		seedSizes[seedInfos[i].k]++;
	}
	std::cout << "Seed size histogram:\n";
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		if (seedSizes[i] > 0) {
			maxCount = std::max(maxCount, seedSizes[i]);
			std::cout << i << ": " << seedSizes[i] << "\n";
		}
	}
	std::cout << "We'd suggest to discard all seeds with count >= " << maxCount / 2 << ".\n";
	for (size_t i = 0; i < seedSizes.size(); ++i) {
		if (seedSizes[i] > 0 && seedSizes[i] < maxCount / 2) {
			std::cout << "New suggested minK: " << i << ". Ignoring seeds with smaller k now.\n";
			return i;
		}
	}
	return options.minK;
}*/

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t maxK, size_t flankWidth) {
	PresenceChecker seedingPresenceChecker(presenceChecker);
	std::cout << "Precomputing posToTaxon array...\n";
	std::vector<uint16_t> posTaxonArray(concat.getSuffixArray().size(), 0);
#pragma omp parallel for
	for (size_t i = 0; i < concat.getSuffixArray().size(); ++i) {
		size_t tID = posToTaxon(concat.getSuffixArray()[i], concat.getTaxonCoords(), concat.getConcatenatedSeq().size(),
				options.reverseComplement);
		posTaxonArray[i] = tID;
	}

	std::cout << "Extracting seeded blocks...\n";
	std::vector<SeedInfo> seededBlockInfos = extractSeedInfos(concat, seedingPresenceChecker, minK, maxK, options, posTaxonArray);
	std::cout << "seeded block infos.size(): " << seededBlockInfos.size() << "\n";
	std::cout << "Processing seeds...\n";
	std::sort(seededBlockInfos.begin(), seededBlockInfos.end(), std::greater<SeedInfo>());

	//printSeedSizeHistogram(seededBlockInfos, options);
	if (options.overriddenK) {
		std::vector<std::pair<size_t, size_t> > seedSizes = countSeedSizes(seededBlockInfos, options);
		printSeedSizeHistogram(seedSizes);
		size_t newMinK = elbowMethod(seedSizes, options);
		std::cout << "New chosen minK by using the elbow method: " << newMinK << ". Ignoring all seeds with smaller k that this value.\n";

		//size_t newMinK = rechooseMinK(seededBlockInfos, options);
		minK = newMinK;
		flankWidth = newMinK; // TODO: maybe remove me again?
	}

	ApproximateMatcher approxMatcher(options.mismatchesOnly);
	std::vector<SeedInfo> buffer;
	size_t lastN = seededBlockInfos[0].n;
	buffer.push_back(seededBlockInfos[0]);
	for (size_t i = 1; i < seededBlockInfos.size(); ++i) {
		if (seededBlockInfos[i].n == lastN) {
			if (seededBlockInfos[i].k >= minK) {
				buffer.push_back(seededBlockInfos[i]);
			}
		} else {
			std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
			std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher,
					flankWidth, posTaxonArray);
			std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
			processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat);
			buffer.clear();
			lastN = seededBlockInfos[i].n;
			if (seededBlockInfos[i].k >= minK) {
				buffer.push_back(seededBlockInfos[i]);
			}
		}
	}
	// process the last buffer
	std::cout << "Number of seeded blocks with " << lastN << " exact matches: " << buffer.size() << "\n";
	std::vector<ExtendedBlock> extendedBlockBuffer = processSeedInfoBuffer(buffer, concat, presenceChecker, options, approxMatcher,
			flankWidth, posTaxonArray);
	std::cout << "Assembling " << extendedBlockBuffer.size() << " extended blocks with " << lastN << " exact taxa...\n";
	processExtendedBlockBuffer(extendedBlockBuffer, options, stats, writer, concat);

	double seqDataUsed = stats.getAmountSeqDataUsed(concat.getSequenceDataSize());
	std::cout << "seqDataUsed: " << seqDataUsed << "\n";
	if (seqDataUsed < options.minSeqDataUsage) {
		std::cout << "Current percentage of sequence data used: " << seqDataUsed * 100
				<< " %. This is too low. Trying to find more blocks with lower minimum kmer size.\n";
		size_t newMinK = std::max((size_t) 8, (size_t) minK - 1);
		if (newMinK != minK) {
			std::cout << "Using new value for minK: " << newMinK << "\n";
			extractExtendedBlocks(concat, presenceChecker, writer, stats, options, newMinK, maxK, newMinK);
		} else {
			std::cout << "Unfortunately, a smaller minK value makes not much sense. :-(\n";
		}
	}

	// TODO: Remove me again, this is just out of curiosity
	//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
	//std::cout << "We have made " << superseeds.size() << " superseeds out of " << seededBlocks.size() << " seeded blocks.\n";
}
