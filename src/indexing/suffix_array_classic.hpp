#pragma once

#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "../indexing/suffix_array_construction/radix.hpp"
#include "../options.hpp"
#include "../dna_functions.hpp"

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
	void buildSuffixArray(const std::string& seq, size_t nTotalSites, unsigned int lTop, const std::string& sequencesPath,
			const Options& options) {
		std::string saPath = sequencesPath + ".sa";
		std::string lcpPath = sequencesPath + ".lcp";
		if (options.reverseComplement) {
			saPath += ".with_rc";
			lcpPath += ".with_rc";
		} else {
			saPath += ".without_rc";
			lcpPath += ".without_rc";
		}
		std::ifstream saFile(saPath);
		std::ifstream lcpFile(lcpPath);
		if (saFile.good() && lcpFile.good()) {
			std::cout << "Found SA and LCP array. No need to recompute, loading them instead...\n";
			readFromFiles(saFile, lcpFile);
			std::cout << "Finished reading SA and LCP from file.\n";
		} else {
			// First heavy lifting: build the suffix array according to Marius's code
			std::cout << "Building suffix array...\n";

			size_t* SA_radix = Radix<size_t>(seq, nTotalSites, lTop).build();
			for (size_t i = 0; i < nTotalSites; ++i) {
				SA.push_back(SA_radix[i]);
			}
			delete[] SA_radix;

			// Important to use ltop, the largest kmer among the range possible
			std::cout << "Suffix array built...Working...\n";
			_nTotalSites = nTotalSites;

			std::cout << "Computing longest common prefixes...\n";
			/*lcp.resize(nTotalSites);
			 lcp[0] = 0;
			 for (size_t i = 1; i < nTotalSites; ++i) {
			 lcp[i] = longestCommonPrefix(seq, SA[i - 1], SA[i], lTop);
			 }*/
			computeLCP(seq);
			std::cout << "Finished computation of longest common prefixes.\n";
			writeToFiles(saPath, lcpPath);
		}
	}

	void buildSuffixArray(const std::string& seq, size_t nTotalSites, const Options& options) {
		unsigned int lTop = seq.size();
		buildSuffixArray(seq, nTotalSites, lTop, options.filepath, options);
	}

	void buildSuffixArrayFromFile(const std::string& filepath, const Options& options) {
		std::ifstream infile(filepath);
		std::string text { istreambuf_iterator<char>(infile), istreambuf_iterator<char>() };
		buildSuffixArray(text, text.size(), options);
	}

	void buildSuffixArray(const std::string& seq) {
		// First heavy lifting: build the suffix array according to Marius's code
		std::cout << "Building suffix array...\n";

		size_t* SA_radix = Radix<size_t>(seq, seq.size(), seq.size()).build();
		for (size_t i = 0; i < seq.size(); ++i) {
			SA.push_back(SA_radix[i]);
		}
		delete[] SA_radix;

		// Important to use ltop, the largest kmer among the range possible
		std::cout << "Suffix array built...Working...\n";
		_nTotalSites = seq.size();

		std::cout << "Computing longest common prefixes...\n";
		computeLCP(seq);
		std::cout << "Finished computation of longest common prefixes.\n";
	}

	inline size_t operator[](size_t pos) const {
		return SA[pos];
	}

	const std::vector<size_t>& getSA() const {
		return SA;
	}

	const std::vector<size_t>& getLCP() const {
		return lcp;
	}

	bool checkLCP(const std::string& seq) const {
		// check whether the lcp array is fine
		for (size_t i = 0; i < SA.size() - 1; ++i) {
			for (size_t j = 1; j <= lcp[i]; ++j) {
				if (seq[SA[i] + j - 1] == '$' || seq[SA[i - 1] + j - 1] == '$') {
					return false;
				}
			}
		}
		std::cout << "Checked LCP and it was fine.\n";
		return true;
	}

private:
	void readFromFiles(std::ifstream& saFile, std::ifstream& lcpFile) {
		saFile >> _nTotalSites;
		lcpFile >> _nTotalSites;
		SA.resize(_nTotalSites);
		lcp.resize(_nTotalSites);
		for (size_t i = 0; i < _nTotalSites; ++i) {
			saFile >> SA[i];
			lcpFile >> lcp[i];
		}
		saFile.close();
		lcpFile.close();
	}

	void writeToFiles(const std::string& saPath, const std::string& lcpPath) {
		std::ofstream saOut(saPath);
		std::ofstream lcpOut(lcpPath);
		saOut << SA.size() << "\n";
		lcpOut << lcp.size() << "\n";
		for (size_t i = 0; i < SA.size(); ++i) {
			saOut << SA[i] << "\n";
			lcpOut << lcp[i] << "\n";
		}
		saOut.close();
		lcpOut.close();
	}

	// algorithm from http://web.cs.iastate.edu/~cs548/references/linear_lcp.pdf
	void computeLCP(const std::string& seq) {
		std::vector<size_t> rank(SA.size());
		for (size_t i = 0; i < SA.size(); ++i) {
			rank[SA[i]] = i; // stores for each position in the text where in the suffix array it occurs
		}
		lcp.resize(SA.size());
		size_t h = 0; // size of the current longest common prefix
		for (size_t i = 0; i < SA.size(); ++i) {
			size_t k = rank[i];
			size_t j;
			if (k == 0) {
				lcp[k] = 0;
			} else {
				j = SA[k - 1];
				while (i + h < SA.size() && j + h < SA.size() && ambiguousMatch(seq[i + h], seq[j + h])) {
					h++;
				}
				lcp[k] = h;
			}
			if (h > 0) {
				h--;
			}
		}
		/*
		// fix the lcp array...
#pragma omp parallel for schedule(dynamic)
		for (size_t k = 1; k < SA.size(); ++k) {
			// check lcp[k]
			for (size_t c1 = 0; c1 < lcp[k]; ++c1) {
				if (seq[SA[k] + c1] == '$' || seq[SA[k - 1] + c1] == '$') {
					lcp[k] = c1;
					break;
				}
			}
		}*/

		/*
		 // compute LCP, the normal way...
		 lcp.resize(SA.size());
		 lcp[0] = 0;
		 #pragma omp parallel for schedule(dynamic)
		 for (size_t i = 1; i < SA.size(); ++i) {
		 size_t h = 0;
		 while (ambiguousMatch(seq[SA[i] + h], seq[SA[i-1] + h])) {
		 h++;
		 }
		 lcp[i] = h;
		 }
		 */
	}
	std::vector<size_t> SA;

	size_t _nTotalSites;
	std::vector<size_t> lcp;
};
