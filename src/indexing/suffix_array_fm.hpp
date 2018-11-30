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
private:
	csa_wt<> fm_index;
	lcp_wt<> lcp;
};
