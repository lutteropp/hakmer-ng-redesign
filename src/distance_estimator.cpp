/*
 * distance_estimator.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: Sarah Lutteropp
 */

#include "distance_estimator.hpp"

double DistanceEstimator::distance() {
	if (distValid) {
		return dist;
	} else {
		// TODO: Recompute distance
		distValid = true;
		return dist;
	}
}

void DistanceEstimator::addChars(char c1, char c2) {
	a += c1;
	b += c2;
}

void DistanceEstimator::addCharA(char c) {
	a += c;
}

void DistanceEstimator::addCharB(char c) {
	b += c;
}
