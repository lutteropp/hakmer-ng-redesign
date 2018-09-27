/*
 * quartet_indexing.hpp
 *
 *  Created on: Sep 27, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <array>
#include <vector>
#include <cassert>
#include <algorithm>

class QuartetIndexer {
public:
	QuartetIndexer() {
	}

	void init(size_t n) {
		// fill nCr2, nCr3, and nCr4
		nCr2.resize(n);
		nCr3.resize(n);
		nCr4.resize(n);
		for (size_t i = 0; i < n; ++i) {
			nCr4[i] = binomialCoeff(i, 4);
			nCr3[i] = binomialCoeff(i, 3);
			nCr2[i] = binomialCoeff(i, 2);
		}
	}

// Returns value of Binomial Coefficient C(n, k), code from https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
	size_t binomialCoeff(size_t n, size_t k) {
		size_t res = 1;
		// Since C(n, k) = C(n, n-k)
		if (k > n - k) {
			k = n - k;
		}
		// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
		for (size_t i = 0; i < k; ++i) {
			res *= (n - i);
			res /= (i + 1);
		}
		return res;
	}

// using https://en.wikipedia.org/wiki/Combinatorial_number_system
	size_t quartetToNumber(size_t a, size_t b, size_t c, size_t d) {
		assert((a < b) && (b < c) && (c < d));
		size_t res = 0;
		res += binomialCoeff(d, 4);
		res += binomialCoeff(c, 3);
		res += binomialCoeff(b, 2);
		res += a;
		return res;
	}

// using https://en.wikipedia.org/wiki/Combinatorial_number_system
	std::array<size_t, 4> numberToQuartet(size_t number) {
		// find d maximal with nCr(d,4) <= N
		size_t d = (std::upper_bound(nCr4.begin(), nCr4.end(), number) - nCr4.begin()) - 1;
		// find c maximal with nCr(c,3) <= N - nCr(d,4)
		size_t c = (std::upper_bound(nCr3.begin(), nCr3.end(), number - nCr4[d]) - nCr3.begin()) - 1;
		// find b maximal with nCr(b,2) <= N - nCr(d,4) - nCr(c,3)
		size_t b = (std::upper_bound(nCr2.begin(), nCr2.end(), number - nCr4[d] - nCr3[c]) - nCr2.begin()) - 1;
		// find a maximal with nCr(a,1) <= N - nCr(d,4) - nCr(c,3) - nCr(b,2)
		size_t a = number - nCr4[d] - nCr3[c] - nCr2[b];
		std::array<size_t, 4> res = { a, b, c, d };
		return res;
	}
private:
	std::vector<size_t> nCr4;
	std::vector<size_t> nCr3;
	std::vector<size_t> nCr2;
};
