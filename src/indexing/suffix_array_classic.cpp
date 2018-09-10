#include "../indexing/suffix_array_classic.hpp"

// Better version, as of Nov 10, 2015, handles boundary cases
/**
 * Perform a binary search in the suffix array to find matches.
 * @param patternSeqStartPos Index of the first base of the pattern in question in the input sequence data.
 * @param m Length of the pattern in question in the input sequence data.
 * @param kmerExtended3Prime Result giving the first and last index of the pattern occurrence in the suffix array as well as the total count.
 * @param S The input sequence data.
 */
bool SuffixArrayClassic::binarySearch3Prime(size_t patternSeqStartPos, unsigned int m, std::pair<size_t, size_t>& res,
		const std::string& text) {
	size_t l, h, mid, pos;
	l = 0;
	h = text.size() - 1;
	while (l < h) { // find lower bound
		mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (text.compare(patternSeqStartPos, m, text, pos, m) <= 0) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (text.compare(pos, m, text, patternSeqStartPos, m) != 0) {
		return false;
	} else {
		res.first = l;
	}

	l = 0;
	h = text.size() - 1;
	while (l < h) { // find upper bound
		if (h - l == 1) {
			mid = h; // treats integer rounding in mid calc for boundary case when interval is just [a,a+1]
		} else {
			mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		}
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (text.compare(patternSeqStartPos, m, text, pos, m) >= 0) {
			l = mid;
		} else {
			h = mid - 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (text.compare(pos, m, text, patternSeqStartPos, m) != 0) {
		return false;
	} else {
		res.second = l - 1;
	}
	return true;
}

// Better version, as of Nov 10, 2015, handles boundary cases
/**
 * Perform a binary search in the suffix array to find matches.
 * @param patternSeqStartPos Index of the first base of the pattern in question in the input sequence data.
 * @param m Length of the pattern in question in the input sequence data.
 * @param kmerExtended3Prime Result giving the first and last index of the pattern occurrence in the suffix array as well as the total count.
 * @param S The input sequence data.
 */
bool SuffixArrayClassic::binarySearch3Prime(const std::string& pattern, std::pair<size_t, size_t>& res,
		const std::string& text) {
	size_t l, h, mid, pos;
	size_t m = pattern.size();
	l = 0;
	h = text.size() - 1;
	while (l < h) { // find lower bound
		mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (pattern.compare(0, m, text, pos, m) <= 0) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (text.compare(pos, m, pattern, 0, m) != 0) {
		return false;
	} else {
		res.first = l;
	}

	l = 0;
	h = text.size() - 1;
	while (l < h) { // find upper bound
		if (h - l == 1) {
			mid = h; // treats integer rounding in mid calc for boundary case when interval is just [a,a+1]
		} else {
			mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		}
		pos = SA[mid];
		if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (pattern.compare(0, m, text, pos, m) >= 0) {
			l = mid;
		} else {
			h = mid - 1;
		}
	}
	pos = SA[l];
	if (pos > text.size() - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (text.compare(pos, m, pattern, 0, m) != 0) {
		return false;
	} else {
		res.second = l - 1;
	}
	return true;
}

size_t SuffixArrayClassic::exactMatches(const std::string& pattern, const std::string& text, std::vector<size_t> &matches) {
	std::pair<size_t, size_t> bl;
	binarySearch3Prime(pattern, bl, text);
	for (size_t i = bl.first; i <= bl.second; ++i) {
		matches.push_back(SA[i]);
	}
	return matches.size();
}
size_t SuffixArrayClassic::exactMatches(size_t patternStartPos, unsigned int m, const std::string& text,
		std::vector<size_t> &matches) {
	std::pair<size_t, size_t> bl;
	binarySearch3Prime(patternStartPos, m, bl, text);
	for (size_t i = bl.first; i <= bl.second; ++i) {
		matches.push_back(SA[i]);
	}
	return matches.size();
}

size_t SuffixArrayClassic::exactMatches(const std::string& pattern, const std::string& text, size_t firstIdx,
		std::vector<size_t> &matches) {
	throw std::runtime_error("Not implemented yet");
}

size_t SuffixArrayClassic::exactMatches(size_t patternStartPos, unsigned int m, const std::string& text, size_t firstIdx,
		std::vector<size_t> &matches) {
	size_t j = firstIdx + 1;
	size_t blockSize = 1;
	while (j <= _nTotalSites) {
		if (text.compare(SA[firstIdx], m, text, SA[j], m) == 0) {
			++blockSize;
			++j;
		} else {
			break;
		}
	}
	//fill the matches array
	matches.clear();
	matches.reserve(blockSize);
	for (size_t i = 0; i < blockSize; ++i) {
		matches.push_back(SA[firstIdx + i]);
	}
	return matches.size();
}

size_t SuffixArrayClassic::countExactMatches(size_t firstSAIndex, unsigned int m) {
	size_t res = 1;
	for (size_t i = firstSAIndex + 1; i < _nTotalSites; ++i) {
		if (lcp[i] >= m) {
			res++;
		} else {
			break;
		}
	}
	return res;
}

const std::vector<size_t>& SuffixArrayClassic::getSA() const {
	return SA;
}

const std::vector<size_t>& SuffixArrayClassic::getLCP() const {
	return lcp;
}
