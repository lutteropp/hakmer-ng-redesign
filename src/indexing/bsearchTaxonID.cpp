#include "../indexing/bsearchTaxonID.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#include "../seq_data.hpp"
#include "../utils.hpp"

// Better version, as of Nov 10, 2015, handles boundary cases
/**
 * Perform a binary search in the suffix array to find matches.
 * @param firstPos Index of the first base of the pattern in question in the input sequence data.
 * @param lastPos Index of the last base of the pattern in question in the input sequence data.
 * @param res Result giving the first and last index of the pattern occurrence in the suffix array.
 * @param S The input sequence data.
 * @param SA The suffix array.
 */
bool binarySearch3PrimeFirstLast(cint firstPos, cint lastPos, std::pair<cint, cint>& res, const SeqData& S,
		SuffixArrayWrapper& SA) {
	cint l, h, mid, pos;
	l = 0;
	h = S.nTotalSites - 1;
	size_t m = lastPos - firstPos + 1;
	while (l < h) { // find lower bound
		mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		pos = SA[mid];
		if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (S.seq.compare(firstPos, m, S.seq, pos, m) <= 0) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	pos = SA[l];
	if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (S.seq.compare(pos, m, S.seq, firstPos, m) != 0) {
		return false;
	} else {
		res.first = l;
	}

	l = 0;
	h = S.nTotalSites - 1;
	while (l < h) { // find upper bound
		if (h - l == 1) {
			mid = h; // treats integer rounding in mid calc for boundary case when interval is just [a,a+1]
		} else {
			mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		}
		pos = SA[mid];
		if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (S.seq.compare(firstPos, m, S.seq, pos, m) >= 0) {
			l = mid;
		} else {
			h = mid - 1;
		}
	}
	pos = SA[l];
	if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (S.seq.compare(pos, m, S.seq, firstPos, m) != 0) {
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
 * @param res Result giving the first and last index of the pattern occurrence in the suffix array.
 * @param S The input sequence data.
 * @param SA The suffix array.
 */
bool binarySearch3Prime(cint patternSeqStartPos, unsigned int m, std::pair<cint, cint>& res, const SeqData& S,
		SuffixArrayWrapper& SA) {
	cint l, h, mid, pos;
	l = 0;
	h = S.nTotalSites - 1;
	while (l < h) { // find lower bound
		mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		pos = SA[mid];
		if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (S.seq.compare(patternSeqStartPos, m, S.seq, pos, m) <= 0) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	pos = SA[l];
	if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (S.seq.compare(pos, m, S.seq, patternSeqStartPos, m) != 0) {
		return false;
	} else {
		res.first = l;
	}

	l = 0;
	h = S.nTotalSites - 1;
	while (l < h) { // find upper bound
		if (h - l == 1) {
			mid = h; // treats integer rounding in mid calc for boundary case when interval is just [a,a+1]
		} else {
			mid = l + (h - l) / 2; // calculate it this way to avoid overflows instead of (l+h)/2
		}
		pos = SA[mid];
		if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
			return false;
		}
		if (S.seq.compare(patternSeqStartPos, m, S.seq, pos, m) >= 0) {
			l = mid;
		} else {
			h = mid - 1;
		}
	}
	pos = SA[l];
	if (pos > S.nTotalSites - m - 1) { /*printf ("Seq array out of bounds in bsearch; bailing\n");*/
		return false;
	}
	if (S.seq.compare(pos, m, S.seq, patternSeqStartPos, m) != 0) {
		return false;
	} else {
		res.second = l - 1;
	}
	return true;
}

