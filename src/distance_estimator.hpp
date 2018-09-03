/*
 * distance_estimator.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <string>

class DistanceConfig {
public:
	size_t insertionPenalty = 1;
	size_t deletionPenalty = 1;
	size_t substitutionPenalty = 1;
};

class DistanceEstimator { // currently uses simple edit distance
public:
	DistanceEstimator(const std::string& s1, const std::string& s2) : a(s1), b(s2), dist(0) {};
	double distance();
	void addChars(char c1, char c2);
private:
	std::string a;
	std::string b;
	double dist;
};
