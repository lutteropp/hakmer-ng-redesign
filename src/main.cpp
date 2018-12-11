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
#include <chrono>

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
	// as in Mauve
	double sum = 0;
	for (size_t i = 0; i < concat.nTax(); ++i) {
		sum += (double) concat.getTaxonCoords(i).getTotalLength() / concat.nTax();
	}
	return log(sum);
}

void matrixCallback(Options& options) {
	IndexedConcatenatedSequence concat = readConcat(options);
	size_t estK = estimateMinK(concat);
	std::cout << "Estimated minK: " << estK << "\n";

	PresenceChecker presenceChecker(concat, options.reverseComplement);
	SummaryStatistics stats(concat.nTax());
	BlockWriter writer(concat.nTax(), options);
	extractExtendedBlocks(concat, presenceChecker, writer, stats, options, estK, std::numeric_limits<size_t>::max());
	stats.printSummaryStatistics(concat, options);

	if (!options.outpath.empty()) {
		writer.assembleFinalSupermatrix(concat.getTaxonLabels(), options.outpath, options);
		std::cout << "Supermatrix written to: " << options.outpath << "\n";
	}
}

void reportCallback(const Options& options) {
	printInputStatistics(options.filepath, options.contigs);
}

int main(int argc, char* argv[]) {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	std::string commandStr = "";
	for (size_t i = 0; i < argc; ++i) {
		commandStr += argv[i];
		if (i < argc - 1) {
			commandStr += " ";
		}
	}
	std::cout << "Program called with the following arguments: " << commandStr << "\n";

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

	size_t nThreads = 0;
	app.add_option("-t,--threads", nThreads, "Maximum number of threads to use.");

	app.add_flag("--redo", options.redo, "Redo the run, overwrite old result files.");
	app.add_option("-o,--outpath", options.outpath, "Path to the output file to be written."); //->required();
	app.add_option("-i,--info,--infopath", options.infopath, "Path to the optional run-information file to be written.");

	app.add_option("--minTaxa", options.minTaxaPerBlock, "Minimum number of taxa per block.", true);
	app.add_option("-s,--minSeqData,--minSeqUsage", options.minSeqDataUsage,
			"Minimum percentage of sequence data to be used. Has to be a value between 0 and 1. If less sequence data has been used, minimum k-mer seed size gets reduced.",
			true);

	CLI11_PARSE(app, argc, argv);
	if (!options.outpath.empty() && options.infopath.empty()) {
		options.infopath = options.outpath + ".info";
	}

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
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "Program runtime: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " s" << std::endl;

}
