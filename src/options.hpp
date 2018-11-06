/*
 * options.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <limits>

class Options {
public:
	std::string filepath = "";
	std::string outpath = "";
	std::string infopath = "";
	bool contigs = false;
	bool reportAndExit = false;
	bool noIndels = false;
	bool verbose = false;
	bool reverseComplement = false;
	bool proteinData = false;
	double maxDelta = 0.1;
	bool largeSeeds = false;

	bool noQuartets = true;
	bool redo = false;

	bool dynamicFlanks = false;
	size_t flankWidth = 20;
	size_t minTaxaPerBlock = 4;
	bool quickDelta = false;
	bool useHMM = false;
	bool hmm_gcCorrection = false;
	bool useBigGapsCriterion = false; // TODO: This still ahs to be implemented... stop alignment extension as soon as there are too big gaps.
	size_t maxMismatches = 0;

	// quartet-flavor-only parameters
	std::string speciesTreePath = "";
	std::string geneTreesPath = "";
	std::string multiSPAMPath = "";
	size_t minBlocksPerQuartet = 100;
	size_t maxBlocksPerQuartet = std::numeric_limits<size_t>::max();
	bool sampleQuartetBlocks = false;
	bool concatenatedDistance = false;
	bool concatenatedMSA = true;
	bool noFinalMSA = false;
	bool noPartitions = true;
	bool quartetFlavor = false;

	// internal parameters
	bool verboseDebug = false;
	size_t minK = 24;
	size_t maxK = std::numeric_limits<size_t>::max(); // if maxK == minK, then we don't do dynamic extension
	bool largeMem = false;
	bool jukesCantor = true;
	bool qicScoring = false;
	bool storeQICCounts = false;
	double minTimesLarger = 3.0;
	bool externalIndexing = false;
	bool discardNs = true;
	bool sameBlockCountAsMultiSPAM = false;
	size_t maxTaxaPerBlock = std::numeric_limits<size_t>::max();
	size_t maximumExtensionWidth = 1000;
	bool lowComplexity = false; // keep low complexity k-mers?
	size_t earlyStopCount = 50; // stop flank extension if earlyStopCount extra bases didn't improve the delta score
	bool fixedFlanks = false; // in non-dynamic extension mode: Enforce that the entire flank width is present?
	bool preselectSeeds = false;
};

