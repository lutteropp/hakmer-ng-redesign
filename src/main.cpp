/*
 * main.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>
#include <limits>

#include "external/CLI11.hpp"
#include "options.hpp"

void evalQuartetsCallback(const Options& options) {

}

void quartetsCallback(const Options& options) {

}

void matrixCallback(const Options& options) {

}

int main(int argc, char* argv[]) {
	Options options;
	CLI::App app { "App description" };
	app.add_option("-f,--file", options.filepath,
			"Program hakmer (\"homology aware k-mers\"): Reads a set of sequences with some homology between them\n"
					"and identifies \"k-mer blocks\" comprising sets of exact or nearly exact k-mer matches appearing only once\n"
					"per sequence, plus (nearby) upstream and downstream flanking sequences.")->required()->check(CLI::ExistingFile);

	app.add_option("--redo", options.redo, "Redo the run");
	app.add_option("-o,--outpath", options.outpath, "Outpath")->required();

	auto dynamicFlanksOption = app.add_flag("-d,--dynamic", options.dynamicFlanks, "Optional description");

	app.add_option("--flankwidth", options.flankWidth, "Optional description", true)->excludes(dynamicFlanksOption)->check(
			CLI::Range(0, std::numeric_limits<int>::max()));
	app.add_option("--maxdelta", options.maxDelta, "Optional description", true)->needs(dynamicFlanksOption)->check(CLI::Range(0, 1));

	auto evalMode = app.add_subcommand("evaluateQuartets", "Quartet evaluation mode");
	evalMode->add_option("--speciesTree", options.speciesTreePath, "Path to the species tree")->check(CLI::ExistingFile);

	auto quartetsMode = app.add_subcommand("quartets", "Quartets mode");
	quartetsMode->add_option("--minblocks", options.minBlocksPerQuartet, "Optional description", true);

	auto supermatrixMode = app.add_subcommand("matrix", "Supermatrix mode");

	app.require_subcommand(1);
	CLI11_PARSE(app, argc, argv);

	if (app.got_subcommand(evalMode)) {
		evalQuartetsCallback(options);
	} else if (app.got_subcommand(quartetsMode)) {
		quartetsCallback(options);
	} else if (app.got_subcommand(supermatrixMode)) {
		matrixCallback(options);
	}
}
