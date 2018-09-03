/*
 * RadixLSDCache.h
 *
 *  Created on: Jul 19, 2012
 *      Author: marius
 */

#pragma once

#include <vector>
#include <cstring>

#include "../../indexing/suffix_array_construction/MergeSorter.hpp"

template<class DataT, class KeyT, class S>
class RadixLSDCache {
	static const int charBuckets = 1 << (8 * sizeof(uchar));
	static const int shortBuckets = 1 << (8 * sizeof(ushort));
	static const S mergeSortThreshold = 512;
	static const S doubleCharSortThreshold = shortBuckets / 2;

	S bucketSizeShort[sizeof(KeyT) / sizeof(ushort)][shortBuckets];
	S bucketSizeChar[sizeof(KeyT) / sizeof(uchar)][charBuckets];
	vector<KeyT> tKey;
	vector<DataT> tData;
	MergeSorter<DataT, KeyT, S> ss;

public:

	void sort(S length, DataT *data, KeyT *key) {
		if (length < mergeSortThreshold)
		ss.sort(length, data, key);
		else if (length < doubleCharSortThreshold)
		sortP<uchar, charBuckets>(length, data, key, bucketSizeChar);
		else
		sortP<ushort, shortBuckets>(length, data, key, bucketSizeShort);
	}

	template<class T, int N>
	void sortP(S length, DataT *data, KeyT *key, DataT bucketSize[][N]) {
		tData.reserve(length);
		tKey.reserve(length);
		DataT *tD = &tData[0];
		KeyT *tK = &tKey[0];

		const int passes = sizeof(KeyT) / sizeof(T);

		for (int r = 0; r < passes; ++r)
		memset(bucketSize[r], 0, N * sizeof(DataT));

		S i = 0;
		for (T *k = (T*) key; i < length; ++i, k += passes)
		for (int r = passes; r; --r) {
			int o = getLSDOffset(passes, r);
			++bucketSize[r - 1][k[o]];
		}

		for (int r = passes; r; --r) {
			// convert bucketSizes into bucketStarts
			DataT *bStart = bucketSize[r - 1];
			inplaceShiftedBy1PrefixSum(bStart, N);

			int o = getLSDOffset(passes, r);
			T *k = ((T*) key) + o;
			for (i = 0; i < length; ++i, k += passes) {
				DataT b = bStart[*k]++;
				tD[b] = data[i];
				tK[b] = key[i];
			}

			DataT *t = tD;
			tD = data;
			data = t;

			KeyT *t2 = tK;
			tK = key;
			key = t2;
		}
	}
};
