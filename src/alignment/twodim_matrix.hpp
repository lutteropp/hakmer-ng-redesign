/*
 * twodim_matrix.hpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>

template<typename T>
class TwoDimMatrix {
public:
	TwoDimMatrix() :
			n(0), m(0) {
	}
	void init(size_t n, size_t m) {
		this->n = n;
		this->m = m;
		entries.resize(n);
		for (size_t i = 0; i < m; ++i) {
			entries[i].resize(m);
		}
	}
	void addRow() {
		n++;
		entries.resize(n);
		entries[n - 1].resize(m);
	}
	void addColumn() {
		m++;
		for (size_t i = 0; i < entries.size(); ++i) {
			entries[i].resize(m);
		}
	}
	T& entryAt(size_t i, size_t j) {
		return entries[i][j];
	}
	size_t getN() const {
		return n;
	}
	size_t getM() const {
		return m;
	}
	void shrinkDownTo(size_t newN, size_t newM) {
		entries.resize(newN);
		for (size_t i = 0; i < newN; ++i) {
			entries[i].resize(newM);
			entries[i].shrink_to_fit();
		}
		entries.shrink_to_fit();
		n = newN;
		m = newM;
	}
	size_t size() const {
		return entries.size();
	}
	size_t size(size_t idx) const {
		return entries[idx].size();
	}
private:
	size_t n;
	size_t m;
	std::vector<std::vector<T> > entries;
};
