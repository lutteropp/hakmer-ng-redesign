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

#include <unordered_set>
#include <cmath>
#include <queue>
#include <algorithm>

#include "indexing/suffix_array_classic.hpp"

#include "build_superseeds.hpp"
#include "block_extension.hpp"

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

bool noOverlap(size_t actSAPos, size_t matchCount, const std::vector<std::pair<size_t, size_t> >& extraOccs, size_t k, size_t nTax,
		const std::vector<size_t>& SA, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const std::string& T, const Options& options, size_t maxMatches) {
	size_t concatSize = T.size();

	size_t flankOffset = 0;

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

bool acceptSeed(size_t actSAPos, size_t matchCount, const std::vector<std::pair<size_t, size_t> >& extraOccs, size_t k, size_t nTax,
		const std::vector<size_t>& SA, PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords,
		const std::string& T, const Options& options, size_t maxMatches) {
	if (matchCount + extraOccs.size() > nTax || matchCount + extraOccs.size() > maxMatches /*|| matchCount > options.maxTaxaPerBlock*/) { // easy test for paralogy
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

	if (!options.trimSeeds) {
		return noOverlap(actSAPos, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches);
	} else {
		return true;
	}
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

void trimSeededBlockSimple(Seed& block, PresenceChecker& presenceChecker) {
	// first, trim the left side
	while (block.getAverageSeedSize() > 0) {
		bool ok = true;
		for (size_t i = 0; i < block.getTaxonIDsInBlock().size(); ++i) {
			size_t coord = block.getSeedCoords(block.getTaxonIDsInBlock()[i]).first;
			if (!presenceChecker.isFree(coord)) {
				ok = false;
				break;
			}
		}
		if (ok) {
			break;
		} else {
			// trim back the left side by 1
			block.increaseAllTaxonCoordsLeft();
		}
	}
	// then, trim the right side (if seed.size() is not already zero)
	while (block.getAverageSeedSize() > 0) {
		bool ok = true;
		for (size_t i = 0; i < block.getTaxonIDsInBlock().size(); ++i) {
			size_t coord = block.getSeedCoords(block.getTaxonIDsInBlock()[i]).second;
			if (!presenceChecker.isFree(coord)) {
				ok = false;
				break;
			}
		}
		if (ok) {
			break;
		} else {
			// trim back the rightz side by 1
			block.decreaseAllTaxonCoordsRight();
		}
	}
}

void trimSeededBlock(Seed& block, PresenceChecker& presenceChecker, const Options& options) {
	trimSeededBlockSimple(block, presenceChecker);
	if (!options.simpleTrimming) {
		std::vector<size_t> taxIDs = block.getTaxonIDsInBlock();
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
}

std::vector<Seed> extractSeededBlocks(const std::string& T, size_t nTax, const std::vector<size_t>& SA, const std::vector<size_t>& lcp,
		PresenceChecker& presenceChecker, const std::vector<IndexedTaxonCoords>& taxonCoords, const Options& options, size_t minK,
		size_t minMatches, size_t maxMatches) {
	std::vector<Seed> res;
	size_t actSAPos = 0;
	double lastP = 0;

	ApproximateMatcher approxMatcher(true);

#pragma omp parallel for schedule(dynamic)
	for (size_t sIdx = 0; sIdx < SA.size(); ++sIdx) {
		size_t startPos = SA[sIdx];
		size_t k = minK;
		if ((startPos + k >= T.size() || !presenceChecker.isFree(startPos, startPos + k - 1))) {
			continue;
		}
		size_t matchCount = countMatches(sIdx, lcp, k);
		std::vector<std::pair<size_t, size_t> > extraOccs;

		while (matchCount >= std::max(options.minTaxaPerBlock, minMatches)) {
			//we can stop if at least one of the exact matches we've found occurs before the current match in the suffix array... if we don't do that, we don't find all matches!
			bool stopEarly = false;
			if (lcp[sIdx] >= k) {
				stopEarly = true;
			}

			if (!stopEarly && acceptSeed(sIdx, matchCount, extraOccs, k, nTax, SA, presenceChecker, taxonCoords, T, options, maxMatches)) {
				Seed block(nTax);
				for (size_t i = sIdx; i < sIdx + matchCount; ++i) {
					block.addTaxon(posToTaxon(SA[i], taxonCoords, T.size(), options.reverseComplement), SA[i], SA[i] + k - 1);
				}

				if (options.maxMismatches > 0) {
					std::string pattern = T.substr(startPos, k);
					extraOccs = approxMatcher.findOccurrences(T, SA, presenceChecker, pattern, options.maxMismatches, 1, false);

					// TODO: Maybe only add those approximate matches that don't collide with the exact matches we already have?
					for (size_t i = 0; i < extraOccs.size(); ++i) {
						size_t taxID = posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement);
						if (!block.hasTaxon(taxID)) {
							size_t flankOffset = 0;
							// check if the extra occurrence is fine
							if (flankOffset > extraOccs[i].first || extraOccs[i].second + flankOffset >= T.size()) {
								continue;
							}
							if (!presenceChecker.isFree(extraOccs[i].first - flankOffset, extraOccs[i].second + flankOffset)) {
								continue;
							}
							block.addTaxon(posToTaxon(extraOccs[i].first, taxonCoords, T.size(), options.reverseComplement),
									extraOccs[i].first, extraOccs[i].second);
						}
					}
				}

				if (options.trimSeeds) {
					trimSeededBlock(block, presenceChecker, options);
				}

				if (block.getAverageSeedSize() > 0 && block.getNTaxInBlock() >= options.minTaxaPerBlock) {
#pragma omp critical
					{
						res.push_back(block);
						if (options.verboseDebug) {
							std::cout << "Pushing back a seeded block with " << block.getNTaxInBlock() << " taxa and seed size "
									<< block.getAverageSeedSize() << "\n";
						}
					}

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

					break;
				} else {
					break;
				}
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

std::vector<Seed> filterSeededBlocks(std::vector<Seed>& seededBlocks, const std::string& T, size_t nTax,
		PresenceChecker& seedingPresenceChecker, const Options& options) {
	std::vector<Seed> res;
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		if (seedingPresenceChecker.isFine(seededBlocks[i])) {
			trivialExtension(seededBlocks[i], T, seedingPresenceChecker, nTax, options);
			res.push_back(seededBlocks[i]);
			seedingPresenceChecker.reserveSeededBlock(seededBlocks[i]);
		} else if (options.trimSeeds) {
			trimSeededBlock(seededBlocks[i], seedingPresenceChecker, options);
			if (seededBlocks[i].getNTaxInBlock() >= options.minTaxaPerBlock && seededBlocks[i].getAverageSeedSize() > 0
					&& seedingPresenceChecker.isFine(seededBlocks[i])) {
				trivialExtension(seededBlocks[i], T, seedingPresenceChecker, nTax, options);
				res.push_back(seededBlocks[i]);
				seedingPresenceChecker.reserveSeededBlock(seededBlocks[i]);
			}
		}
	}
	return res;
}

void processBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, PresenceChecker& seedingPresenceChecker,
		BlockWriter& writer, SummaryStatistics& stats, const Options& options, const std::vector<Seed>& seededBlocks,
		bool useSeedngPresenceChecker) {
	double lastP = 0;
	for (size_t i = 0; i < seededBlocks.size(); ++i) {
		Seed seededBlock = seededBlocks[i];
		if (!presenceChecker.isFine(seededBlock))
			continue;
		trivialExtension(seededBlock, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
		ExtendedBlock extendedBlock = extendBlock(seededBlock, concat.getConcatenatedSeq(), concat.nTax(), presenceChecker, options);
		// check if the extended block can still be accepted.
		if (presenceChecker.isFine(extendedBlock)) {
			presenceChecker.reserveExtendedBlock(extendedBlock);
			if (!useSeedngPresenceChecker) {
				seedingPresenceChecker.reserveExtendedBlock(extendedBlock);
			}
			std::vector<std::string> msa = extendedBlock.msaWrapper.assembleMSA();
			if (options.verboseDebug) {
				std::cout << "Pushing back a block with alignment: \n";
				for (size_t i = 0; i < msa.size(); ++i) {
					std::cout << msa[i] << "\n";
				}
			}
			stats.updateSummaryStatistics(extendedBlock, concat.nTax());
			writer.writeTemporaryBlockMSA(extendedBlock, concat.nTax());
		}
		double progress = (double) 100 * i / seededBlocks.size();
		if (progress > lastP + 1) {
			std::cout << progress << " %\n";
			lastP = progress;
		}
	}
}

void extractExtendedBlocks(const IndexedConcatenatedSequence& concat, PresenceChecker& presenceChecker, BlockWriter& writer,
		SummaryStatistics& stats, const Options& options, size_t minK, size_t nMin, size_t nMax) {
	std::cout << "Extracting seeded blocks...\n";
	std::vector<Seed> seededBlocks;
	PresenceChecker seedingPresenceChecker(presenceChecker);
	for (size_t i = nMax; i >= nMin; i--) {
		std::vector<Seed> actSeededBlocks = extractSeededBlocks(concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(),
				concat.getLcpArray(), seedingPresenceChecker, concat.getTaxonCoords(), options, minK, i, i);
		actSeededBlocks = filterSeededBlocks(actSeededBlocks, concat.getConcatenatedSeq(), concat.nTax(), seedingPresenceChecker, options);
		std::sort(actSeededBlocks.begin(), actSeededBlocks.end(), std::greater<Seed>());
		std::cout << "Found " << actSeededBlocks.size() << " new seeded blocks with " << i << " matches.\n";

		if (options.iterativeExtension) {
			processBlocks(concat, presenceChecker, seedingPresenceChecker, writer, stats, options, actSeededBlocks, false);
			std::cout << "Finished extension of the newly discovered seeded blocks with at least " << i << " matches.\n";
			std::cout << "Current extendedBlocks.size(): " << stats.getCurrentNBlocks() << "\n";
		} else {
			for (size_t j = 0; j < actSeededBlocks.size(); ++j) {
				seededBlocks.push_back(actSeededBlocks[j]);
			}
		}
	}

	if (!options.iterativeExtension) {
		std::cout << "seededBlocks.size(): " << seededBlocks.size() << "\n";
		std::cout << "Assembling extended blocks...\n";
		// TODO: Remove me again, this is just out of curiosity
		//std::vector<Superseed> superseeds = buildSuperseeds(seededBlocks, concat.getConcatenatedSeq(), presenceChecker, concat.nTax(), options);
		//std::cout << "We have made " << superseeds.size() << " superseeds out of " << seededBlocks.size() << " seeded blocks.\n";
		processBlocks(concat, presenceChecker, seedingPresenceChecker, writer, stats, options, seededBlocks, true);
	}

}
