/*
 * Block.hpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <string>

class BlockSequence {
public:
	BlockSequence(size_t taxonID, size_t firstCoord, size_t lastCoord, const std::string& seq);
	size_t getTaxonID();
	size_t getFirstCoord();
	size_t getLastCoord();
	std::string getSeq();
private:
	size_t taxonID;
	size_t firstCoord;
	size_t lastCoord;
	std::string seq;
};

BlockSequence::BlockSequence(size_t taxonID, size_t firstCoord, size_t lastCoord, const std::string& seq) :
		taxonID(taxonID), firstCoord(firstCoord), lastCoord(lastCoord), seq(seq) {
}

size_t BlockSequence::getTaxonID() {
	return taxonID;
}

size_t BlockSequence::getFirstCoord() {
	return firstCoord;
}

size_t BlockSequence::getLastCoord() {
	return lastCoord;
}

std::string BlockSequence::getSeq() {
	return seq;
}

class Block {
private:
	std::vector<BlockSequence> entries;
	size_t width;
};
