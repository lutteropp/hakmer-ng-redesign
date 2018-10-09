/*
 * split_pairwise_alignment.cpp
 *
 *  Created on: Oct 9, 2018
 *      Author: Sarah Lutteropp
 */

#include <algorithm>
#include <iostream>

#include "split_pairwise_alignment.hpp"

SplitPairwiseAlignment::SplitPairwiseAlignment() {
	aliValid = false;
}
void SplitPairwiseAlignment::addCharsLeft(char a, char b) {
	leftFlank.addChars(a, b);
	aliValid = false;
}
void SplitPairwiseAlignment::addCharsRight(char a, char b) {
	rightFlank.addChars(a, b);
	aliValid = false;
}
void SplitPairwiseAlignment::setSeed(const std::string& seed1, const std::string& seed2) {
	if (seed1.size() != seed2.size()) {
		throw std::runtime_error("Different-sized seeds not supported yet");
	}
	for (size_t i = 0; i < seed1.size(); ++i) {
		seed.addChars(seed1[i], seed2[i]);
	}
	aliValid = false;
}
std::pair<std::string, std::string> SplitPairwiseAlignment::extractAlignment() {
	if (aliValid) {
		return ali;
	}
	std::string s1Aligned;
	std::string s2Aligned;
	std::pair<std::string, std::string> leftFlankAlignment = leftFlank.extractAlignment();

	//std::cout << "Left flank alignment: \n" << leftFlankAlignment.first << "\n" << leftFlankAlignment.second << "\n";

	std::reverse(leftFlankAlignment.first.begin(), leftFlankAlignment.first.end());
	std::reverse(leftFlankAlignment.second.begin(), leftFlankAlignment.second.end());
	s1Aligned += leftFlankAlignment.first;
	s2Aligned += leftFlankAlignment.second;
	std::pair<std::string, std::string> seedAlignment = seed.extractAlignment();
	s1Aligned += seedAlignment.first;
	s2Aligned += seedAlignment.second;
	std::pair<std::string, std::string> rightFlankAlignment = rightFlank.extractAlignment();
	s1Aligned += rightFlankAlignment.first;
	s2Aligned += rightFlankAlignment.second;
	ali = std::make_pair(s1Aligned, s2Aligned);
	aliValid = true;
	return ali;
}
void SplitPairwiseAlignment::printAlignment() {
	if (!aliValid) {
		extractAlignment();
	}
	std::cout << ali.first << "\n" << ali.second << "\n\n";
}
size_t SplitPairwiseAlignment::getAlignmentWidth() {
	return leftFlank.getAlignmentWidth() + seed.getAlignmentWidth() + rightFlank.getAlignmentWidth();
}
std::string SplitPairwiseAlignment::getS1() const {
	std::string s1Left = leftFlank.getS1();
	std::reverse(s1Left.begin(), s1Left.end());
	return s1Left + seed.getS1() + rightFlank.getS1();
}
std::string SplitPairwiseAlignment::getS2() const {
	std::string s2Left = leftFlank.getS2();
	std::reverse(s2Left.begin(), s2Left.end());
	return s2Left + seed.getS2() + rightFlank.getS2();
}
double SplitPairwiseAlignment::pairwiseDistance() {
	return leftFlank.pairwiseDistance() + seed.pairwiseDistance() + rightFlank.pairwiseDistance();
}
void SplitPairwiseAlignment::shrinkDownToLeftFlank(size_t newLeftFlankSize) {
	leftFlank.shrinkDownTo(newLeftFlankSize, newLeftFlankSize);
	aliValid = false;
}
void SplitPairwiseAlignment::shrinkDownToRightFlank(size_t newRightFlankSize) {
	rightFlank.shrinkDownTo(newRightFlankSize, newRightFlankSize);
	aliValid = false;
}
