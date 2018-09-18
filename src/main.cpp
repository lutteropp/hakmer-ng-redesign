/*
 * main.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <limits>
#include <omp.h>

#include "external/CLI11.hpp"
#include "options.hpp"

void reportCallback(const Options& options) {

}

void evalQuartetsCallback(const Options& options) {

}

void quartetsCallback(const Options& options) {

}

void matrixCallback(const Options& options) {

}

int main(int argc, char* argv[]) {
	Options options;
	CLI::App app { "Program hakmer (\"homology aware k-mers\"): Reads a set of sequences with some homology between them\n"
			"and identifies \"k-mer blocks\" comprising sets of exact or nearly exact k-mer matches appearing only once\n"
			"per sequence, plus (nearby) upstream and downstream flanking sequences." };
	app.add_option("-f,--file", options.filepath,
			"Input a single file containing sequences in Fasta format (each sequence is treated as separate taxon)")->required()->check(
			CLI::ExistingFile);
	app.add_flag("-c,--contigs", options.contigs, "Expect that the file specified by -f will be a 2-column tab delimited file that points to a set of contig files.\n"
			"    The first column in each row contains a file name, the second contains a taxon label to be used for all contigs in that file.\n"
			"    Thus, all contigs within a file are treated as belonging to the same taxon.");
	auto reportMode = app.add_subcommand("report", "Report on input sequences and exit without analysis");

	app.add_flag("--gapfree,--noindels,--nogaps", options.noIndels, "Build gap-free alignments.");
	app.add_flag("-v,--verbose", options.verbose, "Print progress updates.");
	app.add_flag("--debug", options.verboseDebug, "Print debug output.");
	auto revCompOption = app.add_flag("--revcomp,-r", options.reverseComplement, "Also consider reverse-complement matches of DNA data.");
	app.add_flag("--protein", options.proteinData, "The sequences are protein data instead of DNA data.")->excludes(revCompOption);

	size_t nThreads = 0;
	app.add_option("-t,--threads", nThreads, "Maximum number of threads to use.");

	app.add_option("--redo", options.redo, "Redo the run, overwrite old result files.");
	app.add_option("-o,--outpath", options.outpath, "Path to the output file to be written.")->required();

	auto dynamicFlanksOption = app.add_flag("-d,--dynamic", options.dynamicFlanks, "Optional description");

	app.add_option("--flankwidth", options.flankWidth, "Optional description", true)->excludes(dynamicFlanksOption)->check(
			CLI::Range(0, std::numeric_limits<int>::max()));
	app.add_option("--maxdelta", options.maxDelta, "Maximum delta-score to be still considered tree-like.", true)->needs(dynamicFlanksOption)->check(CLI::Range(0, 1));

	auto evalMode = app.add_subcommand("evaluateQuartets", "Quartet evaluation mode");
	evalMode->add_option("--speciesTree", options.speciesTreePath, "Path to the trusted species tree topology.")->check(CLI::ExistingFile);
	evalMode->add_option("--geneTrees", options.geneTreesPath, "Path to the gene trees.")->check(CLI::ExistingFile);
	evalMode->add_option("--multispam", options.multiSPAMPath, "Path to the quartet topologies file from multi-SpaM.")->check(CLI::ExistingFile);

	auto quartetsMode = app.add_subcommand("quartets", "Quartets mode");
	quartetsMode->add_option("--minblocks", options.minBlocksPerQuartet, "Minimum number of blocks to be sampled for each quartet.", true);
	quartetsMode->add_option("--maxblocks", options.maxBlocksPerQuartet, "Maximum number of blocks to be sampled for each quartet.");
	quartetsMode->add_flag("-s,--sample", options.sampleQuartetBlocks, "Sample a small number of blocks per quartet.");
	quartetsMode->add_flag("--concatDist", options.concatenatedDistance, "Concatenate quartet block distances instead of majority voting the topology.");
	quartetsMode->add_flag("--concatMSA", options.concatenatedMSA, "Concatenate quartet blocks and then infer the quartet topology via RAxML.");
	quartetsMode->add_flag("--noPartitions", options.noPartitions, "Do not build a partitioned MSA, use a single partition instead.");

	auto supermatrixMode = app.add_subcommand("matrix", "Supermatrix mode");

	app.require_subcommand(1);
	CLI11_PARSE(app, argc, argv);

	if (nThreads > 0) {
		omp_set_num_threads(nThreads);
	}
	if (app.got_subcommand(evalMode)) {
		evalQuartetsCallback(options);
	} else if (app.got_subcommand(quartetsMode)) {
		quartetsCallback(options);
	} else if (app.got_subcommand(supermatrixMode)) {
		matrixCallback(options);
	} else if (app.got_subcommand(reportMode)) {
		reportCallback(options);
	}
}
