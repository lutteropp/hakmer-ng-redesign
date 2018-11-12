#pragma once

#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "../options.hpp"

using namespace sdsl;

class SuffixArrayFM {
public:
	SuffixArrayFM() {
	}
	/**
	 * Build a compressed suffix array (FM index).
	 * @param seq The concatenated input sequences.
	 * @param nTotalSites
	 * @param ltop Maximum kmer size to consider.
	 */
	void buildSuffixArray(const std::string& seq, const Options& options) {
		std::cout << "Building FM index...\n";
		/* write seq to a file.
		 std::ofstream tempout("tempout.txt");
		 tempout << seq << "\n";
		 tempout.close();
		 // build the FM index.
		 construct(fm_index, "tempout.txt", 1);
		 */
		cache_config cc(false); // do not delete temp files after csa construction
		construct_im(fm_index, seq, 1);
		cc.delete_files = true; // delete temp files after lcp construction
		construct_im(lcp, seq, 1);
		std::cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << "\n";
	}
	inline size_t operator[](size_t pos) const {
		return fm_index[pos];
	}
	std::string extractText(size_t firstPos, size_t lastPos) {
		return sdsl::extract(fm_index, firstPos, lastPos);
	}
	size_t exactMatches(const std::string& pattern, std::vector<size_t>& matches) {
		auto locations = sdsl::locate(fm_index, pattern.begin(), pattern.end());
		std::sort(locations.begin(), locations.end());
		matches.reserve(locations.size());
		for (size_t i = 0; i < locations.size(); ++i) {
			matches.push_back(locations[i]);
		}
		return matches.size();
	}
	size_t exactMatches(size_t firstSAIndex, const std::string& pattern, std::vector<size_t>& matches) {
		size_t count = sdsl::count(fm_index, pattern.begin(), pattern.end());

		matches.reserve(count);
		for (size_t i = 0; i < count; ++i) {
			matches.push_back(fm_index[firstSAIndex + i]);
		}
		return matches.size();
	}
	size_t exactMatches(size_t patternStartPos, unsigned int m, std::vector<size_t>& matches) {
		std::string pattern = extractText(patternStartPos, patternStartPos + m - 1);
		return exactMatches(pattern, matches);
	}
	size_t countExactMatches(size_t firstSAIndex, unsigned int m) {
		size_t res = 1;
		for (size_t i = firstSAIndex+1; i < fm_index.size(); ++i) {
			if (lcp[i] >= m) {
				res++;
			} else {
				break;
			}
		}
		return res;
	}

	size_t getLCPEntry(size_t idx) {
		return lcp[idx];
	}
private:
	csa_wt<> fm_index;
	lcp_wt<> lcp;
};
