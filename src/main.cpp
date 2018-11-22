/*
 * main.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <limits>
#include <cassert>
#include <fstream>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "external/CLI11.hpp"

#include "options.hpp"
#include "block_extraction.hpp"
#include "io.hpp"
#include "indexed_concat.hpp"
#include "block_writer.hpp"

#include "summary_stats.hpp"

#include "indexing/suffix_array_fm.hpp"

size_t estimateMinK(const IndexedConcatenatedSequence& concat) {
	// as in Mauve, but we double the resulting number just to be sure
	double sum = 0;
	for (size_t i = 0; i < concat.nTax(); ++i) {
		sum += (double) concat.getTaxonCoords(i).getTotalLength() / concat.nTax();
	}
	return log(sum) * 2;
}

void matrixCallback(Options& options) {
	IndexedConcatenatedSequence concat = readConcat(options);
	size_t estK = estimateMinK(concat);
	std::cout << "Estimated minK: " << estK << "\n";
	if (options.minK == 0) {
		options.minK = estK;
	}
	if (options.minK < estK) {
		std::cout << "WARNING: The provided minK is smaller than the estimated minK value.\n";
	}

	PresenceChecker presenceChecker(concat, options.reverseComplement);
	SummaryStatistics stats;
	BlockWriter writer(concat.nTax(), options);
	extractExtendedBlocks(concat, presenceChecker, writer, stats, options, options.minK, options.minTaxaPerBlock, concat.nTax());
	stats.printSummaryStatistics(concat.nTax(), concat.getSequenceDataSize(), options);
	writer.assembleFinalSupermatrix(concat.getTaxonLabels(), options.outpath, options);
}

void reportCallback(const Options& options) {
	printInputStatistics(options.filepath, options.contigs);
}

int main(int argc, char* argv[]) {
	Options options;
	CLI::App app { "Program hakmer (\"homology aware k-mers\"): Reads a set of sequences with some homology between them\n"
			"and identifies \"k-mer blocks\" comprising sets of exact or nearly exact k-mer matches appearing only once\n"
			"per sequence, plus (nearby) upstream and downstream flanking sequences." };
	app.add_option("-f,--file", options.filepath,
			"Input a single file containing sequences in Fasta format (each sequence is treated as separate taxon)")->required()->check(
			CLI::ExistingFile);
	app.add_flag("-c,--contigs", options.contigs,
			"Expect that the file specified by -f will be a 2-column tab delimited file that points to a set of contig files.\n"
					"    The first column in each row contains a file name, the second contains a taxon label to be used for all contigs in that file.\n"
					"    Thus, all contigs within a file are treated as belonging to the same taxon.");
	bool reportMode = false;
	app.add_flag("--report", reportMode, "Report on input sequences and exit without analysis");

	app.add_flag("--gapfree,--noindels,--nogaps", options.noIndels, "Build gap-free alignments.");
	app.add_flag("-v,--verbose", options.verbose, "Print progress updates.");
	app.add_flag("--debug,--verboseDebug", options.verboseDebug, "Print debug output.");
	auto revCompOption = app.add_flag("--revcomp,-r", options.reverseComplement, "Also consider reverse-complement matches of DNA data.");
	app.add_flag("--protein", options.proteinData, "The sequences are protein data instead of DNA data.")->excludes(revCompOption);
	app.add_option("--kmin", options.minK, "Minimum kmer seed size.");
	app.add_option("--kmax", options.maxK, "Maximum kmer seed size.");
	app.add_option("-q,--mismatches", options.maxMismatches, "Maximum number of mismatches in a k-mer block seed.", true);

	size_t nThreads = 0;
	app.add_option("-t,--threads", nThreads, "Maximum number of threads to use.");

	app.add_flag("--redo", options.redo, "Redo the run, overwrite old result files.");
	app.add_option("-o,--outpath", options.outpath, "Path to the output file to be written.")->required();
	app.add_option("-i,--info,--infopath", options.infopath, "Path to the optional run-information file to be written.");

	app.add_option("--flankwidth", options.flankWidth,
			"Maximum size of flanking sequence kept on each side of k-mer. The side of a resulting k-mer block is at most 2*flankwidth+k.", true);

	app.add_option("--minTaxa", options.minTaxaPerBlock, "Minimum number of taxa per block.", true);

	CLI11_PARSE(app, argc, argv);

#ifdef WITH_OPENMP
	if (nThreads > 0) {
		omp_set_num_threads(nThreads);
	}
#endif
	if (!reportMode) {
		matrixCallback(options);
	} else {
		reportCallback(options);
	}
}
