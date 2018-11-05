/*
 * hmm_classification.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: sarah
 */

#pragma once

#include "external/HomologyHMM/parameters.h"
#include "external/HomologyHMM/homology.h"
#include "options.hpp"
#include "alignment/msa_wrapper.hpp"

static char charmap[128];
inline char* getCharmap() {
	static bool initialized = false;
	if (initialized)
		return charmap;
	memset(charmap, 0, 128);
	charmap['a'] = 0;
	charmap['c'] = 1;
	charmap['g'] = 2;
	charmap['t'] = 3;
	charmap['-'] = 4;
	charmap['A'] = 0;
	charmap['C'] = 1;
	charmap['G'] = 2;
	charmap['T'] = 3;
	charmap['-'] = 4;
	initialized = true;
	return charmap;
}
// a mapping from pairwise alignment columns to HomologyHMM emission codes
// row/column indices are as given by the charmap above (ACGT- == 01234).
static char colmap[5][5] = {
//    A   C   G   T   -
		{ '1', '3', '4', '5', '7' },	// A
		{ '3', '2', '6', '4', '7' },  // C
		{ '4', '6', '2', '3', '7' },  // G
		{ '5', '4', '3', '1', '7' },  // T
		{ '7', '7', '7', '7', '\0' },  // -
		};

inline std::string predict(std::string& column_states, const Params& hmm_params) {
	std::string prediction;
	run(column_states, prediction, hmm_params);
	return prediction;
}

inline std::string pairwiseAlignmentToColumnStates(const std::string& s1Aligned, const std::string& s2Aligned) {
	std::string colStates;
	assert(s1Aligned.size() == s2Aligned.size());
	for (size_t i = 0; i < s1Aligned.size(); ++i) { // TODO: This only works with non-ambiguous DNA characters for now...
		colStates += colmap[charmap[s1Aligned[i]]][charmap[s2Aligned[i]]];
	}
	return colStates;
}

inline size_t findNumGoodSites(const std::string& s1Aligned, const std::string& s2Aligned, const Params& params) {
	std::string colStates = pairwiseAlignmentToColumnStates(s1Aligned, s2Aligned);
	std::string prediction = predict(colStates, params);
	for (size_t i = 0; i < prediction.size(); ++i) {
		if (prediction[i] == 'N') {
			return i;
		}
	}
	return prediction.size();
}

inline size_t findNumGoodSitesMSA(MSAWrapper& msaWrapper, const Params& params) {
	size_t goodSites = std::numeric_limits<size_t>::max();
	std::vector<std::string> msa = msaWrapper.assembleMSA();
	for (size_t i = 0; i < msa.size(); ++i) {
		for (size_t j = i+1; j < msa.size(); ++j) {
			goodSites = std::min(goodSites, findNumGoodSites(msa[i], msa[j], params));
		}
	}
	return goodSites;
}

inline Params prepareHmmParams(const std::string& concat, const Options& options) {
	Params params;
	if (options.hmm_gcCorrection) { // TODO: Also support ambiguity characters here and stuff like that... for now, only A,C,G,T characters are considered for computing the GC content.
		size_t nGC = 0;
		size_t nucCount = 0;
		for (size_t i = 0; i < concat.size(); ++i) {
			if (concat[i] == 'G' || concat[i] == 'C') {
				nGC++;
				nucCount++;
			} else if (concat[i] == 'A' || concat[i] == 'T') {
				nucCount++;
			}
		}
		if (nucCount == 0)
			nucCount = 1;
		double gcContent = (double) nGC / nucCount;
		params = getAdaptedHoxdMatrixParameters(gcContent);
	} else {
		params = getHoxdParams();
	}
	float pGoHomo = 0.0005f;
	float pGoUnrelated = 0.000001f;
	params.iGoHomologous = pGoHomo;
	params.iGoUnrelated = pGoUnrelated;
	params.iStartHomologous = 0.5;
	return params;
}
