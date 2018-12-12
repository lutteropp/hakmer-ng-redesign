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
#include <iostream>
#include <sstream>

#include "alignment/simple_msa.hpp"

std::string buildTempFilename(size_t taxID, const Options& options) {
	return options.filepath + "_temporary_taxon_alignment_file_" + std::to_string(taxID);
}

BlockWriter::BlockWriter(size_t nTax, const Options& options) :
		tempMSAFiles(nTax) {
	for (size_t i = 0; i < nTax; ++i) {
		std::string tempFilename = buildTempFilename(i, options);
		tempMSAFiles[i].open(tempFilename);
	}
}

void printBlockMSA(ExtendedBlock& block, const std::vector<std::string>& msa) {
	std::cout << "block MSA:\n";
	for (size_t i = 0; i < msa.size(); ++i) {
		std::cout << msa[i] << "\n";
	}
}

void BlockWriter::writeTemporaryBlockMSA(ExtendedBlock& block, const std::string& T, size_t nTax, const Options& options) {
	std::vector<std::string> msa = computeMSA(block, T, nTax, options);
/*#pragma omp critical
	printBlockMSA(block, msa);*/
#pragma omp critical
	for (size_t i = 0; i < nTax; ++i) {
		tempMSAFiles[i] << msa[i];
	}
}

void BlockWriter::writeTemporaryBlockMSA(const std::vector<std::string>& msa, size_t nTax) {
	for (size_t i = 0; i < nTax; ++i) {
		tempMSAFiles[i] << msa[i];
	}
}


std::string slurp(std::ifstream& in) {
	std::stringstream sstr;
	sstr << in.rdbuf();
	return sstr.str();
}

void BlockWriter::assembleFinalSupermatrix(const std::vector<std::string>& taxonLabels, const std::string& outpath,
		const Options& options) {
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
