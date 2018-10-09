/*
 * main.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <limits>
#include <cassert>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "external/CLI11.hpp"
#include "external/ProgressBar.hpp"

#include "options.hpp"
#include "block_extraction.hpp"
#include "io.hpp"
#include "indexed_concat.hpp"

#include "quartet_topology.hpp"
#include "block_helper_functions.hpp"
#include "quartet_lookup_table.hpp"
#include "mafft_raxml_wrapper.hpp"
#include "quartet_topology_checker.hpp"

#include "quartet_indexer.hpp"

QuartetTopology inferQuartet(size_t a, size_t b, size_t c, size_t d, const IndexedConcatenatedSequence& concat, const Options& options) {
	std::vector<size_t> wantedTaxa = { a, b, c, d };
	std::vector<std::pair<size_t, size_t> > taxCoords;
	for (size_t i = 0; i < concat.nTax(); ++i) {
		taxCoords.push_back(std::make_pair(concat.getTaxonCoords(i).getFirstCoord(), concat.getTaxonCoords(i).getLastCoord()));
	}
	std::pair<std::vector<size_t>, std::vector<size_t> > shrunk = shrinkArrays(concat, taxCoords, wantedTaxa, options);
	PresenceChecker presenceChecker(concat, options.reverseComplement);

	std::vector<AlignedBlock> alignedBlocks = extractAlignedBlocks(concat.getConcatenatedSeq(), concat.nTax(), shrunk.first, shrunk.second,
			presenceChecker, taxCoords, options);
	if (options.concatenatedMSA) { // this is the one that uses raxml-ng for quartet inference
		std::array<std::string, 4> concatMSA = { "", "", "", "" };
		for (size_t i = 0; i < alignedBlocks.size(); ++i) {
			for (size_t j = 0; j < 4; ++j) {
				concatMSA[j] += extractTaxonSequence(alignedBlocks[i], j);
			}
		}
		std::string prefix = concat.getTaxonLabels()[a] + "_" + concat.getTaxonLabels()[b] + "_" + concat.getTaxonLabels()[c] + "_"
				+ concat.getTaxonLabels()[d];
		return inferTopology(prefix, concatMSA, options);
	} else if (options.concatenatedDistance) { // concatenated distances
		std::array<double, 6> concatDist = { 0, 0, 0, 0, 0, 0 };
		for (size_t i = 0; i < alignedBlocks.size(); ++i) {
			std::vector<double> pwdist = alignedBlocks[i].getPairwiseNormalizedDistances(options);
			for (size_t j = 0; j < 6; ++j) {
				concatDist[j] += pwdist[j];
			}
		}
		std::array<double, 6> pairwiseDistances =
				{ concatDist[0], concatDist[1], concatDist[2], concatDist[3], concatDist[4], concatDist[5] };
		return topologyFromDistances(pairwiseDistances);
	} else { // inferring distance-based topologies for each block separately, then doing a majority vote
		std::array<size_t, 3> counts = { 0, 0, 0 }; // ab|cd, ac|bd, ad|bc
		for (size_t i = 0; i < alignedBlocks.size(); ++i) {
			std::vector<double> pwdist = alignedBlocks[i].getPairwiseNormalizedDistances(options);
			std::array<double, 6> pairwiseDistances = { pwdist[0], pwdist[1], pwdist[2], pwdist[3], pwdist[4], pwdist[5] };
			QuartetTopology topo = topologyFromDistances(pairwiseDistances);
			switch (topo) {
			case QuartetTopology::AB_CD:
				counts[0]++;
				break;
			case QuartetTopology::AC_BD:
				counts[1]++;
				break;
			case QuartetTopology::AD_BC:
				counts[2]++;
				break;
			default:
				break;
			}
		}
		// perform a majority-voting for now...
		if (counts[0] > counts[1] && counts[0] > counts[2]) {
			return QuartetTopology::AB_CD;
		} else if (counts[1] > counts[0] && counts[1] > counts[2]) {
			return QuartetTopology::AC_BD;
		} else if (counts[2] > counts[0] && counts[2] > counts[1]) {
			return QuartetTopology::AD_BC;
		} else {
			return QuartetTopology::STAR;
		}
	}

	return QuartetTopology::STAR;
}

void quartetsCallback(const Options& options) {
#ifndef WITH_GENESIS
	if (!options.speciesTreePath.empty()) {
		throw std::runtime_error("Quartet evaluation mode requires genesis integration.");
	}
#endif
	IndexedConcatenatedSequence concat = readConcat(options);
	TopologyChecker topoChecker;

	QuartetIndexer quartetIndexer;
#ifdef WITH_OPENMP
	quartetIndexer.init(concat.nTax());
#endif

	if (!options.speciesTreePath.empty()) {
		topoChecker.init(concat, options.speciesTreePath, options.geneTreesPath, options.multiSPAMPath);
	}

	size_t nCorrect = 0;
	size_t nWrong = 0;
	size_t nStar = 0;
	size_t n = concat.nTax();
	QuartetLookupTable<size_t> quartetTable(n);
	size_t nQuartets = n * (n - 1) * (n - 2) * (n - 3) / 24;
	ProgressBar progressBar(nQuartets);

#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+:nCorrect,nWrong,nStar)
	for (size_t idx = 0; idx < nQuartets; ++idx) {
		std::array<size_t, 4> vals = quartetIndexer.numberToQuartet(idx);
		size_t a = vals[0];
		size_t b = vals[1];
		size_t c = vals[2];
		size_t d = vals[3];
#else
	for (size_t a = 0; a < n - 3; ++a) {
		for (size_t b = a + 1; b < n - 2; ++b) {
			for (size_t c = b + 1; c < n - 1; ++c) {
				for (size_t d = c + 1; d < n; ++d) {
#endif
					QuartetTopology topo = inferQuartet(a, b, c, d, concat, options);
					switch (topo) {
					case QuartetTopology::AB_CD:
						quartetTable.get_tuple(a, b, c, d)[quartetTable.tuple_index(a, b, c, d)]++;
						break;
					case QuartetTopology::AC_BD:
						quartetTable.get_tuple(a, c, b, d)[quartetTable.tuple_index(a, c, b, d)]++;
						break;
					case QuartetTopology::AD_BC:
						quartetTable.get_tuple(a, d, b, c)[quartetTable.tuple_index(a, d, b, c)]++;
						break;
					default:
						break;
					}

					if (!options.speciesTreePath.empty()) { // evaluation stuff
						bool correct = topoChecker.sameTopologyAsReference(a, b, c, d, topo);
						if (correct) {
							nCorrect++;
						} else if (topo == QuartetTopology::STAR) {
							nStar++;
						} else {
							nWrong++;
						}
					}
#ifdef WITH_OPENMP
#pragma omp critical
#endif
					progressBar.Update();
#ifndef WITH_OPENMP
				}
			}
		}
#endif
	}
	std::cout << "\n";
	if (!options.speciesTreePath.empty()) {
		std::cout << "nCorrect: " << nCorrect << "\n";
		std::cout << "nWrong: " << nWrong << "\n";
		std::cout << "nStar: " << nStar << "\n";
	}
	writeQuartets(quartetTable, concat.getTaxonLabels(), options.outpath);
}

void matrixCallback(const Options& options) {
	IndexedConcatenatedSequence concat = readConcat(options);
	std::vector<std::pair<size_t, size_t> > taxCoords;
	for (size_t i = 0; i < concat.nTax(); ++i) {
		taxCoords.push_back(std::make_pair(concat.getTaxonCoords(i).getFirstCoord(), concat.getTaxonCoords(i).getLastCoord()));
	}
	PresenceChecker presenceChecker(concat, options.reverseComplement);

	std::vector<AlignedBlock> alignedBlocks = extractAlignedBlocks(concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(),
			concat.getLcpArray(), presenceChecker, taxCoords, options);
	writeFASTASupermatrix(alignedBlocks, concat.getTaxonLabels(), options.outpath);

	/*
	 std::vector<std::string> msa(concat.nTax());
	 size_t actSAPos = 0;
	 double lastP = 0;
	 while (actSAPos < concat.getSuffixArray().size()) {
	 AlignedBlock bl = nextAlignedBlock(actSAPos, concat.getConcatenatedSeq(), concat.nTax(), concat.getSuffixArray(),
	 concat.getLcpArray(), presenceChecker, taxCoords, options);
	 if (bl.getSeedSize() == 0) {
	 break; // no blocks left
	 }
	 for (size_t i = 0; i < concat.nTax(); ++i) {
	 msa[i] += extractTaxonSequence(bl, i);
	 }
	 double progress = (double) 100 * actSAPos / concat.getSuffixArray().size();
	 if (progress > lastP + 1) {
	 std::cout << progress << "\n";
	 lastP = progress;
	 }
	 }
	 writeFASTASupermatrix(msa, concat.getTaxonLabels(), options.outpath);*/
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
	auto reportMode = app.add_subcommand("report", "Report on input sequences and exit without analysis");

	app.add_flag("--gapfree,--noindels,--nogaps", options.noIndels, "Build gap-free alignments.");
	app.add_flag("-v,--verbose", options.verbose, "Print progress updates.");
	app.add_flag("--debug", options.verboseDebug, "Print debug output.");
	auto revCompOption = app.add_flag("--revcomp,-r", options.reverseComplement, "Also consider reverse-complement matches of DNA data.");
	app.add_flag("--protein", options.proteinData, "The sequences are protein data instead of DNA data.")->excludes(revCompOption);
	app.add_option("--kmin", options.minK, "Minimum kmer seed size.");
	app.add_option("--kmax", options.maxK, "Maximum kmer seed size.");

	size_t nThreads = 0;
	app.add_option("-t,--threads", nThreads, "Maximum number of threads to use.");

	app.add_flag("--redo", options.redo, "Redo the run, overwrite old result files.");
	app.add_option("-o,--outpath", options.outpath, "Path to the output file to be written.")->required();

	auto dynamicFlanksOption = app.add_flag("-d,--dynamic", options.dynamicFlanks,
			"Dynamically extend the sequence regions flanking a kmer seed.");

	app.add_option("--flankwidth", options.flankWidth,
			"Length of flanking sequence kept on each side of k-mer. The side of a resulting k-mer block is 2*flankwidth+k.", true)->excludes(
			dynamicFlanksOption);
	app.add_option("--maxdelta", options.maxDelta, "Maximum delta-score to be still considered tree-like.", true)->needs(
			dynamicFlanksOption)->check(CLI::Range(0.0, 1.0));

	auto quartetsMode = app.add_subcommand("quartets", "Quartets mode");
	quartetsMode->add_option("--minblocks", options.minBlocksPerQuartet, "Minimum number of blocks to be sampled for each quartet.", true);
	quartetsMode->add_option("--maxblocks", options.maxBlocksPerQuartet, "Maximum number of blocks to be sampled for each quartet.");
	quartetsMode->add_flag("-s,--sample", options.sampleQuartetBlocks, "Sample a small number of blocks per quartet.");
	quartetsMode->add_flag("--concatDist", options.concatenatedDistance,
			"Concatenate quartet block distances instead of majority voting the topology.");
	quartetsMode->add_flag("--concatMSA", options.concatenatedMSA,
			"Concatenate quartet blocks and then infer the quartet topology via RAxML.");
	quartetsMode->add_flag("--noPartitions", options.noPartitions, "Do not build a partitioned MSA, use a single partition instead.");
	auto evalOption = quartetsMode->add_option("--speciesTree", options.speciesTreePath,
			"Path to the trusted species tree topology. Activates quartet evaluation mode.")->check(CLI::ExistingFile);
	quartetsMode->add_option("--geneTrees", options.geneTreesPath, "Path to the gene trees. Only used in quartet evaluation mode.")->check(
			CLI::ExistingFile)->needs(evalOption);
	quartetsMode->add_option("--multispam", options.multiSPAMPath,
			"Path to the quartet topologies file from multi-SpaM. Only used in quartet evaluation mode.")->check(CLI::ExistingFile)->needs(
			evalOption);

	auto supermatrixMode = app.add_subcommand("matrix", "Supermatrix mode");
	supermatrixMode->add_option("--minTaxa", options.minTaxaPerBlock, "Minimum number of taxa per block.", true);

	app.require_subcommand(1);
	CLI11_PARSE(app, argc, argv);

#ifdef WITH_OPENMP
	if (nThreads > 0) {
		omp_set_num_threads(nThreads);
	}
#endif
	if (app.got_subcommand(quartetsMode)) {
		options.maxTaxaPerBlock = 4;
		quartetsCallback(options);
	} else if (app.got_subcommand(supermatrixMode)) {
		matrixCallback(options);
	} else if (app.got_subcommand(reportMode)) {
		reportCallback(options);
	}
}
