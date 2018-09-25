#pragma once

#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "../indexing/suffix_array_construction/radix.hpp"
#include "../options.hpp"

inline size_t longestCommonPrefix(const std::string& seq, size_t start1, size_t start2, unsigned int lTop) {
	size_t res = 0;
	for (size_t i = 0; i < seq.size(); ++i) {
		if (start1 + i >= seq.size() || start2 + i >= seq.size()) {
			break;
		}
		if (seq[start1 + i] == seq[start2 + i]) {
			res++;
			if (res >= lTop) {
				return res; // bail because lTop is largest k-mer size we're gonna search
			}
		} else {
			break;
		}
	}
	return res;
}

class SuffixArrayClassic {
public:
	SuffixArrayClassic() {
		_nTotalSites = 0;
	}
	/**
	 * Build the suffix array according to Marius's code.
	 * @param seq The concatenated input sequences.
	 * @param nTotalSites
	 * @param ltop Maximum kmer size to consider.
	 */
	void buildSuffixArray(const std::string& seq, size_t nTotalSites, const Options& options) {
		// First heavy lifting: build the suffix array according to Marius's code
		std::cout << "Building suffix array...\n";

		unsigned int lTop = std::min(seq.size(), options.maxK);

		size_t* SA_radix = Radix<size_t>(seq, nTotalSites, lTop).build();
		for (size_t i = 0; i < nTotalSites; ++i) {
			SA.push_back(SA_radix[i]);
		}
		delete[] SA_radix;

		// Important to use ltop, the largest kmer among the range possible
		std::cout << "Suffix array built...Working...\n";
		_nTotalSites = nTotalSites;

		lcp.resize(nTotalSites);
		lcp[0] = 0;
		for (size_t i = 1; i < nTotalSites; ++i) {
			lcp[i] = longestCommonPrefix(seq, SA[i - 1], SA[i], lTop);
		}
	}

	void buildSuffixArrayFromFile(const std::string& filepath, const Options& options) {
		std::ifstream infile(filepath);
		std::string text { istreambuf_iterator<char>(infile), istreambuf_iterator<char>() };
		buildSuffixArray(text, text.size(), options);
	}

	inline size_t operator[](size_t pos) const {
		return SA[pos];
	}

	const std::vector<size_t>& getSA() const;
	const std::vector<size_t>& getLCP() const;

	size_t exactMatches(const std::string& pattern, const std::string& text, size_t firstIdx, std::vector<size_t> &matches);
	size_t exactMatches(size_t patternStartPos, unsigned int m, const std::string& text, size_t firstIdx, std::vector<size_t> &matches);
	size_t exactMatches(const std::string& pattern, const std::string& text, std::vector<size_t> &matches);
	size_t exactMatches(size_t patternStartPos, unsigned int m, const std::string& text, std::vector<size_t> &matches);
	size_t countExactMatches(size_t firstSAIndex, unsigned int m);
	size_t countExactMatches(const std::string& pattern, const std::string& text);
private:
	bool binarySearch3Prime(const std::string& pattern, std::pair<size_t, size_t>& res, const std::string& text);
	bool binarySearch3Prime(size_t patternSeqStartPos, unsigned int m, std::pair<size_t, size_t>& res, const std::string& text);
	std::vector<size_t> SA;

	size_t _nTotalSites;
	std::vector<size_t> lcp;
};
