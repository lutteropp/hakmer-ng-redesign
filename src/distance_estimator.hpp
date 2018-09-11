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

class DistanceEstimator {
public:
	DistanceEstimator() : a(""), b(""), dist(0), distValid(false) {};
	double distance(); // TODO: dist has to be between 0 (identical) and 1.
	void addChars(char c1, char c2);
	void addCharA(char c);
	void addCharB(char c);
private:
	std::string a;
	std::string b;
	double dist;
	bool distValid;
};
