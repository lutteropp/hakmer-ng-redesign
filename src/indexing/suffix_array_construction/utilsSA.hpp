/*
 * utils.h
 *
 *  Created on: Oct 17, 2011
 *      Author: marius
 */

#pragma once

#include <sys/time.h>
#include <sys/types.h>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>

using namespace std;

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned char uchar;

#ifdef LITTLE_ENDIAN_FLAG
#define getBufferStart(buffer, len) ((buffer) + (len) - 1)
#define getMSDOffset(increment, remainingKey) ((remainingKey) - 1)
#define getLSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#else
#define getBufferStart(buffer, len) (buffer)
#define getMSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#define getLSDOffset(increment, remainingKey) ((remainingKey) - 1)
#endif

inline int max(int a, int b) {
	return a > b ? a : b;
}
inline int min(int a, int b) {
	return a < b ? a : b;
}

/**
 * Returns the number of bits required to represent number n
 *
 * NOTE: if you want to know the number of bits required to
 * represent an alphabet of size n
 * you have to call bitsFor(n-1) since
 * n-1 is the largest character in the alphabet!
 */
template<class unum>
inline int bitsFor(unum n) {
	int b = 1;
	if (n >= 0x100000000L) {
		b += 32;
		n = ((unsigned long) n) >> 32;
	}
	if (n >= 0x10000) {
		b += 16;
		n >>= 16;
	}
	if (n >= 0x100) {
		b += 8;
		n >>= 8;
	}
	if (n >= 0x10) {
		b += 4;
		n >>= 4;
	}
	if (n >= 4) {
		b += 2;
		n >>= 2;
	}
	if (n >= 2)
		b += 1;
	return b;
}

template<class T>
inline void swapItems(T *a, T *b) {
	T tmp = *a;
	*a = *b;
	*b = tmp;
}

template<class T>
inline void prefixSum(T *s, T *d, int length) {
	*d++ = *s++;
	for (; --length; ++d, ++s)
		*d = *s + *(d - 1);
}

// Associates unique characters of the input with consecutive numerical indexes starting at 0.
// Returns the size of the input alphabet.
// 'charIndex' will map original characters to their numerical indexes.
// 'bitsPerChar' will be the number of bits necessary to store alphabet indexes.
template<class unum>
inline int indexAlphabet(const std::string& in, unum length, uint *charIndex, int& bitsPerChar) {
	const int maxAlpha = 256;
	memset(charIndex, 0, maxAlpha * sizeof(*charIndex));

	for (unum i = 0; i < length; ++i) {
		charIndex[(unsigned char) in[i]] = 1;
	}

	charIndex[0] -= 1;
	prefixSum(charIndex, charIndex, maxAlpha);

	int alphaSize = charIndex[maxAlpha - 1] + 1;
	bitsPerChar = bitsFor(alphaSize - 1);

	return alphaSize;
}

template<class unum>
inline unum max(unum *a, unum n) {
	unum m = *a;
	for (; --n;)
		if (*++a > m)
			m = *a;
	return m;
}

template<class DataT, class KeyT, class S>
inline void insertSort(DataT* data, S length, KeyT* key) {
	if (length == 2) {
		if (key[1] < key[0]) {
			swapItems(key, key + 1);
			swapItems(data, data + 1);
		}
		return;
	}

	for (S i = 1; i < length; ++i) {
		KeyT kp = key[i];
		DataT t = data[i];
		S j;
		for (j = i; j > 0 && kp < key[j - 1]; --j) {
			key[j] = key[j - 1];
			data[j] = data[j - 1];
		}
		key[j] = kp;
		data[j] = t;
	}
}

#define PRINT_LIMIT 500
inline void printCharData(const uchar *array, int nKey) {
	int e = (nKey > PRINT_LIMIT) ? PRINT_LIMIT : nKey;
	for (int i = 0; i < e; ++i) {
		uchar c = array[i];
		if (c == '\n')
			c = 'N';
		if (c == '\r')
			c = 'n';
		if (c == '\t')
			c = 'T';
		if (c == 0)
			printf("0");
		else
			printf("%c", c);
	}
	if (e < nKey)
		printf("... %d more", nKey - e);
}

template<typename T>
inline void FreeAll(T & t) {
	T tmp;
	t.swap(tmp);
}

inline double clock_diff_to_msec(long clock_diff) {
	return double(clock_diff) / CLOCKS_PER_SEC * 1000;
}

template<class Proc, class Arg>
inline double time_it(Proc proc, Arg a, int N) // returns time in microseconds
		{
	std::clock_t const start = std::clock();
	for (int i = 0; i < N; ++i)
		proc(a);
	std::clock_t const end = std::clock();
	if (clock_diff_to_msec(end - start) < 200)
		return time_it(proc, a, N * 5);
	return clock_diff_to_msec(end - start) / N;
}

template<class T>
inline void reverseArray(T *a, int n) {
	for (int i = 0, j = n - 1; i < j; ++i, --j) {
		T t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
}

/**
 * Takes as input an array a of length n and performs the following transformation:
 * a'[i] becomes (a[0] + a[1] + ... + a[i-1]) for all i = n-1 down to 1
 * a'[0] becomes 0
 *
 * The operation is done in place so a' overwrites a.
 */
template<class T, class S>
inline void inplaceShiftedBy1PrefixSum(T *a, S n) {
	T next = 0;
	for (S i = 0; i < n; ++i) {
		// we're doing in place prefix sum but the destination
		// is offset by 1 to the right;
		// so we save this value before we overwrite the cell
		T t = a[i];
		a[i] = next;
		next += t;
	}
}
