/*
 * simple_msa.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: sarah
 */

#include "simple_msa.hpp"

std::vector<std::string> alignSequences(const std::vector<std::string>& seqs) {
	std::vector<std::string> res;
	MSAWrapper wrapper(false);
	wrapper.init(seqs.size());
	wrapper.setSeeds(seqs);
	return wrapper.assembleMSA();
}
