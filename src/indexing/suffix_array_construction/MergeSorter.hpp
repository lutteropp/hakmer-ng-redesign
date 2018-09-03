/*
 * MergeSorter.h
 *
 *  Created on: Jul 19, 2012
 *      Author: marius
 */

#pragma once

#include <vector>

#include "../../indexing/suffix_array_construction/utilsSA.hpp"
using namespace std;

template<class DataT, class KeyT, class S>
class MergeSorter {
private:
	static const S insertSortThreshold = 32;
	static const S tileSize = 0x4000;
	vector<KeyT> destK;
	vector<DataT> destD;

public:
	MergeSorter() {

	}

	~MergeSorter() {

	}

	void sort(S length, DataT *data, KeyT *key) {
		if (length <= insertSortThreshold) {
			insertSort(data, length, key);
			return;
		}

		mergeSort(length, data, key);
		return;
	}

	void mergeSort(S length, DataT *data, KeyT *key) {
		//	index tiles = (length / insertSortThreshold) + ((length
		//			% insertSortThreshold) ? 1 : 0);

		//	index insSortThr = (logNextPowOfTwo(tiles) & 1) // odd
		//	? (insertSortThreshold >> 1)
		//			: insertSortThreshold;
		S insSortThr = insertSortThreshold;

		S s = 0;
		for (; s + insSortThr < length; s += insSortThr)
		insertSort(data + s, insSortThr, key + s);
		if (length - s > 1)
		insertSort(data + s, length - s, key + s);

		iterativeMergeSort(length, data, key, insSortThr);
	}

	void iterativeMergeSort(S length, DataT *data, KeyT *key, S initialTileSize) {
		// iterative merge sort
		destD.reserve(length);
		destK.reserve(length);
		DataT *tD = &destD[0];
		KeyT *tK = &destK[0];
		for (S bSize = initialTileSize; bSize < length; bSize <<= 1) {
			int n = 0;
			for (S s = 0; s < length; s += (bSize << 1)) {
				S i = s;
				S m = s + bSize;
				if (m > length)
				m = length;
				S j = m;
				S e = m + bSize;
				if (e > length)
				e = length;
				for (; i < m && j < e; ++n) {
					if (key[i] <= key[j]) {
						tK[n] = key[i];
						tD[n] = data[i];
						++i;
					} else {
						tK[n] = key[j];
						tD[n] = data[j];
						++j;
					}
				}
				for (; i < m; ++n, ++i) {
					tK[n] = key[i];
					tD[n] = data[i];
				}
				for (; j < e; ++n, ++j) {
					tK[n] = key[j];
					tD[n] = data[j];
				}
			}

			DataT *tmp = data;
			data = tD;
			tD = tmp;

			KeyT *tmpk = key;
			key = tK;
			tK = tmpk;
		}
		if (data == &destD[0]) {
			memcpy(tK, key, length * sizeof(KeyT));
			memcpy(tD, data, length * sizeof(DataT));
		}
	}

};
