/*
 * simple_msa.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#include <stdexcept>

#include "simple_msa.hpp"
#include "msa_wrapper.hpp"

std::vector<std::string> computeMSA(const std::vector<std::string>& seqs) {
	MSAWrapper wrapper(false);
	wrapper.init(seqs.size());
	wrapper.setSeeds(seqs);
	return wrapper.assembleMSA();
}

std::string createMissingDataString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}

std::vector<std::string> prepareSeqs(const std::vector<SimpleCoords>& seqCoords, const std::string& T) {
	std::vector<std::string> seqs;
	size_t concatSize = 0;
	for (size_t i = 0; i < seqCoords.size(); ++i) {
		if (seqCoords[i].size() == 0) {
			seqs.push_back("");
		} else {
			std::string leftTrim = createMissingDataString(seqCoords[i].leftGapSize);
			std::string middle = T.substr(seqCoords[i].first, seqCoords[i].second + 1 - seqCoords[i].first);
			std::string rightTrim = createMissingDataString(seqCoords[i].rightGapSize);
			std::string concat = leftTrim + middle + rightTrim;
			seqs.push_back(concat);
			if (concatSize > 0) {
				if (concat.size() != concatSize) {
					throw std::runtime_error("This should not happen");
				}
			} else {
				concatSize = concat.size();
			}
		}
	}
	for (size_t i = 0; i < seqs.size(); ++i) {
		if (seqs[i].empty()) {
			seqs[i] = createMissingDataString(concatSize);
		}
	}
	return seqs;
}

std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T) {
	std::vector<std::string> seqs;
	seqs = prepareSeqs(seqCoords, T);
	return computeMSA(seqs);
}

std::string createMissingString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}

std::vector<std::string> computeMSA(const ExtendedBlock& block, const std::string& T, size_t nTax, const Options& options) {
	std::vector<std::string> res;
	std::vector<std::string> leftFlankMSA = computeMSA(block.getLeftFlankCoords(), T);
	std::vector<std::string> seedMSA;
	if (options.mismatchesOnly) {
		seedMSA = prepareSeqs(block.getSeedCoords(), T);
	} else {
		seedMSA = computeMSA(block.getSeedCoords(), T);
	}
	std::vector<std::string> rightFlankMSA = computeMSA(block.getRightFlankCoords(), T);

	for (size_t i = 0; i < nTax; ++i) {
		res.push_back(leftFlankMSA[i] + seedMSA[i] + rightFlankMSA[i]);
	}
	return res;
}
