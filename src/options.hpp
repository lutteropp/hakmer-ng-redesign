/*
 * options.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: Sarah Lutteropp
 */

#ifndef SRC_OPTIONS_HPP_
#define SRC_OPTIONS_HPP_

#include "../DynamicEditDistance/src/dist_config.hpp"

class Options {
public:
	std::string filepath;
	std::string outpath;
	bool contigs;
	bool reportAndExit;
	bool noIndels;
	bool verbose;
	bool reverseComplement;
	bool proteinData;
	double maxDelta;

	bool noQuartets;
	bool redo;

	// quartet-flavor-only parameters
	std::string speciesTreePath;
	std::string geneTreesPath;
	std::string multiSPAMPath;
	size_t minBlocksPerQuartet;
	size_t maxBlocksPerQuartet;
	bool sampleQuartetBlocks;
	bool concatenatedDistance;
	bool concatenatedMSA;
	bool noFinalMSA;
	bool noPartitions;
	bool quickRaxml;

	// internal parameters
	bool verboseDebug = false;
	size_t minK = 8;
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
	size_t maximumExtensionWidth = 300;

	bool dynamicFlanks = true;
	size_t flankWidth = 25;

	DistConfig config;
};

#endif /* SRC_OPTIONS_HPP_ */
