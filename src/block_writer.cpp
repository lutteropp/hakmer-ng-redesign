/*
 * block_writer.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: sarah
 */

#include "block_writer.hpp"
#include <cassert>
#include <sstream>
#include <cstdio>

#include "block_helper_functions.hpp"

std::string buildTempFilename(size_t taxID) {
	return "temporary_taxon_alignment_file_" + std::to_string(taxID);
}

BlockWriter::BlockWriter(size_t nTax) :tempMSAFiles(nTax) {
	for (size_t i = 0; i < nTax; ++i) {
		std::string tempFilename = buildTempFilename(i);
		tempMSAFiles[i].open(tempFilename);
	}
}

void BlockWriter::writeTemporaryBlockMSA(ExtendedBlock& block, size_t nTax) {
	for (size_t i = 0; i < nTax; ++i) {
		tempMSAFiles[i] << extractTaxonSequence(block, i);
	}
}

std::string slurp(std::ifstream& in) {
    std::stringstream sstr;
    sstr << in.rdbuf();
    return sstr.str();
}

void BlockWriter::assembleFinalSupermatrix(const std::vector<std::string>& taxonLabels, const std::string& outpath) {
	assert(tempMSAFiles.size() == taxonLabels.size());
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		tempMSAFiles[i].close();
	}
	std::ofstream outfile(outpath);
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		outfile << ">" << taxonLabels[i] << "\n";
		std::string filename = buildTempFilename(i);
		std::ifstream tempAct(filename);
		outfile << slurp(tempAct) << "\n";
	}
	outfile.close();
	deleteTemporaryFiles();
}
void BlockWriter::deleteTemporaryFiles() {
	for (size_t i = 0; i < tempMSAFiles.size(); ++i) {
		std::string filename = buildTempFilename(i);
		std::remove(filename.c_str());
	}
	tempMSAFiles.clear();
	tempMSAFiles.shrink_to_fit();
}
