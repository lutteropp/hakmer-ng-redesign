/*
 * options.hpp
 *
 *  Created on: Jun 6, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <stddef.h>
#include <limits>
#include <string>

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

	bool redo = false;

	size_t flankWidth = 50;
	size_t minTaxaPerBlock = 4;
	size_t maxMismatches = 1;

	// internal parameters
	bool verboseDebug = false;
	size_t minK = 0;
	size_t maxK = std::numeric_limits<size_t>::max(); // if maxK == minK, then we don't do dynamic extension
	bool externalIndexing = false;
	bool lowComplexity = false; // keep low complexity k-mers?
	size_t maxAllowedSuperseedDistance = 50;
	size_t minSharedSuperseedTax = 4;
	bool mismatchAugmentationOnly = true;
	bool discardNs = true;
	bool simpleExtension = false; // stop extending as soon as one taxon can't be further extended?
	bool mismatchesOnly = true; // allow only mismatches in the seeds - this means we don't need to run a MSA on them.
	size_t minSeedTaxInColumn = 4; // minimum number of taxa in column, needed for seed trimming
	size_t minSeedSitesKept = 8; // needed for seed trimming

	double minSeqDataUsage = 0.25; // minimum amount of overall sequence data to use, default value: 25 %
};

