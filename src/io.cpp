/*
 * io.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: Sarah Lutteropp
 */

#include "io.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unordered_map>

class FASTARecord {
public:
	FASTARecord(const std::string& name, const std::string& seq) :
			name(name), seq(seq) {
	}
	std::string name;
	std::string seq;
};

class ContigRecord {
public:
	ContigRecord() = default;
	ContigRecord(const FASTARecord& fastaRecord);
	std::string name;
	std::vector<FASTARecord> contigs;
};

class InputReader {
public:
	InputReader(bool contigs);
	void openFile(const std::string& filepath);
	ContigRecord readNext();
	std::vector<ContigRecord> readAll();
	bool hasNext();
private:
	std::ifstream infile;
	bool contigs;
};

class ContigStats {
public:
	size_t length;
	size_t nContig;
	size_t Ns;
	size_t runs; // number of 'N' runs
	size_t maxRunLength;
	std::string name;
};

/**
 * Check if the file contains more records.
 */
bool hasNextRecord(std::ifstream& infile) {
	int len = infile.tellg(); // current position in the ifstream
	std::string line;
	getline(infile, line);
	while (line.empty() && infile.good()) { // skip newlines
		std::getline(infile, line);
	}
	bool res = infile.good();
	infile.seekg(len, std::ios_base::beg); // return back to position before reading the line
	return res;
}

bool fastaRecordEnded(std::ifstream& infile) {
	int len = infile.tellg(); // current position in the ifstream
	std::string line;
	getline(infile, line);
	infile.seekg(len, std::ios_base::beg); // return back to position before reading the line
	if (line.empty() || (line.size() > 0 && line[0] == '>')) {
		return true;
	} else {
		return false;
	}
}

FASTARecord readNextFASTA(std::ifstream& infile) {
	std::string line;
	std::getline(infile, line);
	while (line.empty() && infile.good()) { // skip newlines
		std::getline(infile, line);
	}
	if (line[0] != '>') {
		std::cout << line.size() << "\n";
		throw std::runtime_error("the sequence is not in FASTA format");
	}
	std::string name = line.substr(1, std::string::npos);
	std::getline(infile, line);
	std::string seq = line;
	while (!fastaRecordEnded(infile)) {
		std::getline(infile, line);
		seq += line;
	}
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	for (size_t i = 0; i < seq.size(); ++i) {
		/*if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
		 seq[i] = 'N';
		 }*/
	}
	return FASTARecord(name, seq);
}

ContigRecord::ContigRecord(const FASTARecord& fastaRecord) {
	name = fastaRecord.name;
	contigs.clear();
	contigs.push_back(fastaRecord);
}

ContigRecord readNextContig(std::ifstream& infile) {
	std::string line;
	std::getline(infile, line);
	size_t tabPos = line.find_first_of('\t');
	if (line.find_last_of('\t') != tabPos || tabPos == std::string::npos) {
		throw std::runtime_error("Something is wrong with your input file");
	}
	std::string contigFilePath = line.substr(0, tabPos);
	std::string taxonName = line.substr(tabPos + 1);
	ContigRecord contigRecord;
	contigRecord.name = taxonName;

	std::ifstream fastaFile;
	fastaFile.open(contigFilePath);
	if (!infile.good()) {
		throw std::runtime_error("Can not open file: " + contigFilePath);
	}

	while (hasNextRecord(fastaFile)) {
		FASTARecord fastaRecord = readNextFASTA(fastaFile);
		contigRecord.contigs.push_back(fastaRecord);
	}
	return contigRecord;
}

void InputReader::openFile(const std::string& filepath) {
	infile.open(filepath);
	if (!infile.good()) {
		throw std::runtime_error("Can not open file: " + filepath);
	}
}

bool InputReader::hasNext() {
	return hasNextRecord(infile);
}

ContigRecord InputReader::readNext() {
	if (!contigs) {
		FASTARecord fastaRecord = readNextFASTA(infile);
		return ContigRecord(fastaRecord);
	} else {
		return readNextContig(infile);
	}
}

std::vector<ContigRecord> InputReader::readAll() {
	infile.clear();
	infile.seekg(0, std::ios::beg);
	std::vector<ContigRecord> vec;
	while (hasNextRecord(infile)) {
		vec.push_back(readNext());
	}
	return vec;
}

InputReader::InputReader(bool contigs) :
		contigs(contigs) {
}

ContigStats computeContigStats(const ContigRecord& record) {
	ContigStats stats;
	stats.nContig = record.contigs.size();
	size_t length = 0;
	size_t numN = 0;
	size_t numNRun = 0;
	size_t maxNRun = 0;
	size_t actNRun = 0;
	for (const FASTARecord& contig : record.contigs) {
		length += contig.seq.size();
		for (size_t i = 0; i < contig.seq.size(); ++i) {
			if (contig.seq[i] == 'N') {
				numN++;
				actNRun++;
				if (i > 0 && contig.seq[i - 1] != 'N') {
					numNRun++;
				}
			} else {
				maxNRun = std::max(maxNRun, actNRun);
				actNRun = 0;
			}
		}

	}
	stats.Ns = numN;
	stats.maxRunLength = maxNRun;
	stats.runs = numNRun;
	stats.length = length;
	stats.name = "[" + record.name.substr(0, record.name.find(" ")) + "]";
	return stats;
}

void printInputStatistics(const std::string& filepath, bool contigs) {
	std::cout << "\nReport on input sequences:\n";
	const char separator = ' ';
	std::cout << std::left << std::setw(8) << std::setfill(separator) << "Seq" << "\t";
	std::cout << std::left << std::setw(14) << std::setfill(separator) << "Length" << "\t";
	std::cout << std::left << std::setw(8) << std::setfill(separator) << "nContig" << "\t";
	std::cout << std::left << std::setw(10) << std::setfill(separator) << "Ns" << "\t";
	std::cout << std::left << std::setw(10) << std::setfill(separator) << "Runs" << "\t";
	std::cout << std::left << std::setw(10) << std::setfill(separator) << "Max Run Length" << "\t";
	std::cout << std::left << std::setfill(separator) << "Name" << "\n";
	InputReader reader(contigs);
	reader.openFile(filepath);
	size_t i = 0;
	while (reader.hasNext()) {
		ContigRecord record = reader.readNext();
		ContigStats stats = computeContigStats(record);
		std::cout << std::left << std::setw(8) << std::setfill(separator) << i << "\t";
		std::cout << std::left << std::setw(14) << std::setprecision(14) << std::setfill(separator) << stats.length << "\t";
		std::cout << std::left << std::setw(8) << std::setprecision(8) << std::setfill(separator) << stats.nContig << "\t";
		std::cout << std::left << std::setw(10) << std::setprecision(10) << std::setfill(separator) << stats.Ns << "\t";
		std::cout << std::left << std::setw(10) << std::setprecision(10) << std::setfill(separator) << stats.runs << "\t";
		std::cout << std::left << std::setw(10) << std::setprecision(10) << std::setfill(separator) << stats.maxRunLength << "\t";
		std::cout << std::left << std::setfill(separator) << stats.name << "\n";
		++i;
	}
}

std::string revComp(const std::string& str) {
	static std::unordered_map<char, char> mapping; // TODO: Make this one static somehow
	mapping['A'] = 'T';
	mapping['T'] = 'A';
	mapping['C'] = 'G';
	mapping['G'] = 'C';
	mapping['N'] = 'N';
	// DNA ambiguity codes, see http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	mapping['Y'] = 'R';
	mapping['R'] = 'Y';
	mapping['W'] = 'W';
	mapping['S'] = 'S';
	mapping['K'] = 'M';
	mapping['M'] = 'K';
	mapping['D'] = 'H';
	mapping['V'] = 'B';
	mapping['H'] = 'D';
	mapping['B'] = 'V';
	mapping['X'] = 'X';

	std::string res;
	for (int i = str.size() - 1; i >= 0; --i) {
		if (mapping.find(str[i]) != mapping.end()) {
			res += mapping[str[i]];
		} else {
			std::string c;
			c += str[i];
			std::cout << "i: " << i << "str.size(): " << str.size() << "\n";
			throw std::runtime_error("Cannot reverse-complement the following base: " + c);
		}
	}
	return res;
}


IndexedConcatenatedSequence readConcat(const Options& options) {
	std::vector<IndexedTaxonCoords> coords;

	size_t coordOffset = 0;
	std::string concat = "";

	InputReader reader(options.contigs);
	reader.openFile(options.filepath);
	while (reader.hasNext()) {
		ContigRecord record = reader.readNext();
		std::vector<std::string> contigs(record.contigs.size());
		for (size_t i = 0; i < record.contigs.size(); ++i) {
			std::string seq = record.contigs[i].seq;
			if (options.discardNs) {
				std::string noNSeq = "";
				for (size_t j = 0; j < record.contigs[i].seq.size(); ++j) {
					if (record.contigs[i].seq[j] != 'N' && record.contigs[i].seq[j] != 'n' && record.contigs[i].seq[j] != '-') {
						noNSeq += record.contigs[i].seq[j];
					}
				}
				seq = noNSeq;
			}
			contigs[i] = seq;
		}

		IndexedTaxonCoords itc(record.name.substr(0, record.name.find(" ")), contigs, coordOffset, options);
		coords.push_back(itc);
		coordOffset = itc.getLastCoord() + 2;

		for (size_t i = 0; i < contigs.size(); ++i) {
			concat += contigs[i];
			concat += "$";
			if (options.reverseComplement) { // TODO: Is this where we want our reverse-complemented contigs to be? Might be better to add them at the end or so.
				concat += revComp(contigs[i]);
				concat += "$";
			}
		}
	}

	IndexedConcatenatedSequence res(concat, coords, options);
	return res;
}
