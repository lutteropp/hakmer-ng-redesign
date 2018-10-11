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
	bool contigs = false;
	bool reportAndExit = false;
	bool noIndels = false;
	bool verbose = false;
	bool reverseComplement = false;
	bool proteinData = false;
	double maxDelta = 0.1;

	bool noQuartets = true;
	bool redo = false;

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

	// internal parameters
	bool verboseDebug = false;
	size_t minK = 24;
	size_t maxK = std::numeric_limits<size_t>::max(); // if maxK == minK, then we don't do dynamic extension
	bool largeMem = false;
	bool jukesCantor = true;
	bool qicScoring = false;
	bool storeQICCounts = false;
	double minTimesLarger = 3.0;
	bool externalIndexing = true;
	bool discardNs = true;
	bool sameBlockCountAsMultiSPAM = false;
	size_t minTaxaPerBlock = 4;
	size_t maxTaxaPerBlock = std::numeric_limits<size_t>::max();
	size_t maximumExtensionWidth = 500;
	bool lowComplexity = false; // keep low complexity k-mers?
	size_t earlyStopCount = 50; // stop flank extension if earlyStopCount extra bases didn't improve the delta score

	bool dynamicFlanks = false;
	size_t flankWidth = 20;
};

