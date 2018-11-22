/*
 * simple_msa.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#include "simple_msa.hpp"

std::vector<std::string> computeMSA(const std::vector<std::string>& seqs) {
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

std::vector<std::string> computeMSA(const std::vector<SimpleCoords>& seqCoords, const std::string& T) {
	std::vector<std::string> res;
	MSAWrapper wrapper(false);
	wrapper.init(seqCoords.size());
	std::vector<std::string> seqs;
	for (size_t i = 0; i < seqCoords.size(); ++i) {
		std::string leftTrim = createGapString(seqCoords[i].leftGapSize);
		std::string middle = T.substr(seqCoords[i].first, seqCoords[i].second + 1 - seqCoords[i].first);
		std::string rightTrim = createGapString(seqCoords[i].rightGapSize);
		seqs.push_back(leftTrim + middle + rightTrim);
	}
	wrapper.setSeeds(seqs);
	return wrapper.assembleMSA();
}
