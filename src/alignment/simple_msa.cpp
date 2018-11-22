/*
 * simple_msa.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#include "simple_msa.hpp"

inline std::vector<std::string> computeMSA(const std::vector<std::string>& seqs) {
	MSAWrapper wrapper(false);
	wrapper.init(seqs.size());
	wrapper.setSeeds(seqs);
	return wrapper.assembleMSA();
}

std::string createGapString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '-';
	}
	return res;
}

inline std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T) {
	std::vector<std::string> seqs;
	for (size_t i = 0; i < seqCoords.size(); ++i) {
		if (seqCoords[i].size() == 0) {
			seqs.push_back("");
		} else {
			std::string leftTrim = createGapString(seqCoords[i].leftGapSize);
			std::string middle = T.substr(seqCoords[i].first, seqCoords[i].second + 1 - seqCoords[i].first);
			std::string rightTrim = createGapString(seqCoords[i].rightGapSize);
			seqs.push_back(leftTrim + middle + rightTrim);
		}
	}
	return computeMSA(seqs);
}

std::string createMissingString(size_t len) {
	std::string res;
	for (size_t i = 0; i < len; ++i) {
		res += '?';
	}
	return res;
}

inline std::vector<std::string> computeMSA(const ExtendedBlock& block, const std::string& T, size_t nTax) {
	std::vector<std::string> res;
	std::vector<std::string> leftFlankMSA = computeMSA(block.getLeftFlankCoords(), T);
	std::vector<std::string> seedMSA = computeMSA(block.getSeedCoords(), T);
	std::vector<std::string> rightFlankMSA = computeMSA(block.getRightFlankCoords(), T);

	for (size_t i = 0; i < nTax; ++i) {
		res.push_back(leftFlankMSA[i] + seedMSA[i] + rightFlankMSA[i]);
	}
	return res;
}
