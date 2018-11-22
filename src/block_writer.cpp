/*
 * block_writer.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: sarah
 */

#include "block_writer.hpp"

#include <stddef.h>
#include <cassert>
#include <cstdio>
#include <sstream>

#include "alignment/simple_msa.hpp"

std::string buildTempFilename(size_t taxID, const Options& options) {
	return options.filepath + "_temporary_taxon_alignment_file_" + std::to_string(taxID);
}

BlockWriter::BlockWriter(size_t nTax, const Options& options) :tempMSAFiles(nTax) {
	for (size_t i = 0; i < nTax; ++i) {
		std::string tempFilename = buildTempFilename(i, options);
		tempMSAFiles[i].open(tempFilename);
	}
}

void BlockWriter::writeTemporaryBlockMSA(ExtendedBlock& block, const std::string& T, size_t nTax, const Options& options) {
	std::vector<std::string> msa = computeMSA(block, T, nTax, options);
	for (size_t i = 0; i < nTax; ++i) {
		//tempMSAFiles[i] << extractTaxonSequence(block, i);
#pragma omp critical
		tempMSAFiles[i] << msa[i];
	}
}

std::string slurp(std::ifstream& in) {
    std::stringstream sstr;
    sstr << in.rdbuf();
    return sstr.str();
}

void BlockWriter::assembleFinalSupermatrix(const std::vector<std::string>& taxonLabels, const std::string& outpath, const Options& options) {
	assert(tempMSAFiles.size() == taxonLabels.size());
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		tempMSAFiles[i].close();
	}
	std::ofstream outfile(outpath);
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		outfile << ">" << taxonLabels[i] << "\n";
		std::string filename = buildTempFilename(i, options);
		std::ifstream tempAct(filename);
		outfile << slurp(tempAct) << "\n";
	}
	outfile.close();
	deleteTemporaryFiles(options);
}
void BlockWriter::deleteTemporaryFiles(const Options& options) {
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		std::string filename = buildTempFilename(i, options);
		std::remove(filename.c_str());
	}
	tempMSAFiles.clear();
	tempMSAFiles.shrink_to_fit();
}
