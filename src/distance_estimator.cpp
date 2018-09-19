/*
 * distance_estimator.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: Sarah Lutteropp
 */

#include <cmath>
#include "distance_estimator.hpp"

inline double jukesCantorCorrection(double dist) {
	// TODO: Jukes Cantor Correction doesn't work if dist >= 0.75. In this case, it will return infinity.
	return -0.75 * std::log(1 - (4.0 / 3) * dist);
}

void HammingDistanceEstimator::addChars(char c1, char c2) {
	a += c1;
	b += c2;
	if (c1 != c2) {
		dist++;
	}
}

double HammingDistanceEstimator::distance() {
	return jukesCantorCorrection(dist / a.size());
}
