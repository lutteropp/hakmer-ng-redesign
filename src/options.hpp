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

	bool redo = false;
	size_t minTaxaPerBlock = 4;

	// internal parameters
	bool proteinData = false;
	bool verboseDebug = false;
	bool externalIndexing = false;
	bool lowComplexity = false; // keep low complexity k-mers?
	size_t maxAllowedSuperseedDistance = 50;
	size_t minSharedSuperseedTax = 4;
	bool discardNs = true;
	bool mismatchesOnly = true; // allow only mismatches in the seeds - this means we don't need to run a MSA on them.
	size_t minSeedTaxInColumn = 4; // minimum number of taxa in column, needed for seed trimming
	size_t minSeedSitesKept = 1; // needed for seed trimming
	bool discardUninformativeBlocks = true;

	bool discardParalogMismatches = true;

	double maxAvgSubstitutionRate = 0.3;
	double maxErrorRate = 0.1;
	size_t maxMismatches = std::numeric_limits<size_t>::max();
};

