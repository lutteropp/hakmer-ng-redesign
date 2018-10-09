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
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "block_helper_functions.hpp"

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
	/*for (size_t i = 0; i < seq.size(); ++i) {
	 if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
	 seq[i] = 'N';
	 }
	 }*/
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
	std::locale loc("");
	infile.imbue(loc);
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

		IndexedTaxonCoords itc(record.name.substr(0, record.name.find(" ")), contigs, coordOffset);
		coords.push_back(itc);
		coordOffset = itc.getLastCoord() + 2;

		for (size_t i = 0; i < contigs.size(); ++i) {
			concat += contigs[i];
			concat += "$";
		}

		if (options.reverseComplement) {
			concat += "$";
			for (size_t i = contigs.size() - 1; i >= 1; --i) {
				concat += revComp(contigs[i]);
				concat += "$";
			}
			concat += revComp(contigs[0]);
		}
	}

	IndexedConcatenatedSequence res(concat, coords, options);
	return res;
}

void writeFASTASupermatrix(const std::vector<std::string>& msa, const std::vector<std::string>& taxonLabels, const std::string& filepath) {
	std::ofstream outfile(filepath);
	std::cout << "writing to file: " << filepath << "\n";
	for (size_t i = 0; i < taxonLabels.size(); ++i) {
		outfile << ">" + taxonLabels[i] << "\n";
		outfile << msa[i] << "\n";
	}
	outfile.close();
}

void writeFASTASupermatrix(const std::vector<AlignedBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath) {
	std::ofstream outfile(filepath);
	std::cout << "writing to file: " << filepath << "\n";
	for (size_t i = 0; i < taxonLabels.size(); ++i) {
		outfile << ">" + taxonLabels[i] << "\n";
		std::string concat = "";
		for (size_t j = 0; j < blocks.size(); ++j) {
			concat += extractTaxonSequence(blocks[j], i);
		}
		outfile << concat << "\n";
	}
	outfile.close();
}
void writeFASTASupermatrix(const std::vector<ExtendedBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath, const std::string& T) {
	std::ofstream outfile(filepath);
	for (size_t i = 0; i < taxonLabels.size(); ++i) {
		outfile << ">" + taxonLabels[i] << "\n";
		std::string concat = "";
		for (size_t j = 0; j < blocks.size(); ++j) {
			concat += extractTaxonSequence(blocks[j], i, T);
		}
		outfile << concat << "\n";
	}
	outfile.close();
}
void writeFASTASupermatrix(const std::vector<SeededBlock>& blocks, const std::vector<std::string>& taxonLabels, const std::string& filepath, const std::string& T) {
	std::ofstream outfile(filepath);
	for (size_t i = 0; i < taxonLabels.size(); ++i) {
		outfile << ">" + taxonLabels[i] << "\n";
		std::string concat = "";
		for (size_t j = 0; j < blocks.size(); ++j) {
			concat += extractTaxonSequence(blocks[j], i, T);
		}
		outfile << concat << "\n";
	}
	outfile.close();
}


std::string makeHeader(const std::vector<std::string>& taxonLabels) {
	std::string res;
	res += std::to_string(taxonLabels.size());
	res += "\n";
	for (size_t i = 0; i < taxonLabels.size(); ++i) {
		res += taxonLabels[i];
		res += "\n";
	}
	return res;
}
std::string makeTopologyString(size_t idx1, size_t idx2, size_t idx3, size_t idx4,
		const std::array<size_t, 3>& counts) {
	std::stringstream ss;
	ss << idx1 << "," << idx2 << "|" << idx3 << "," << idx4 << ":" << counts[0] << "\n";
	ss << idx1 << "," << idx3 << "|" << idx2 << "," << idx4 << ":" << counts[1] << "\n";
	ss << idx1 << "," << idx4 << "|" << idx2 << "," << idx3 << ":" << counts[2] << "\n";
	return ss.str();
}

void parseTopologyString(QuartetLookupTable<size_t>& table, const std::string& line) {
	size_t firstCommaPos = line.find(',', 1);
	size_t idx1 = std::stoul(line.substr(1, firstCommaPos - 1));
	size_t delimPos = line.find('|', firstCommaPos + 1);
	size_t idx2 = std::stoul(line.substr(firstCommaPos + 1, delimPos - 1 - (firstCommaPos + 1) + 1));
	size_t secondCommaPos = line.find(',', delimPos + 1);
	size_t idx3 = std::stoul(line.substr(delimPos + 1, secondCommaPos - 1 - delimPos));
	size_t colonBracketPos = line.find(':', secondCommaPos + 1);
	size_t idx4 = std::stoul(line.substr(secondCommaPos + 1, colonBracketPos - 1 - secondCommaPos));
	size_t count = std::stoul(line.substr(colonBracketPos + 1, line.size()));
	auto& tuple = table.get_tuple(idx1, idx2, idx3, idx4);
	size_t tupIndex = table.tuple_index(idx1, idx2, idx3, idx4);
	tuple[tupIndex] = count;
}

QuartetLookupTable<size_t> readQuartets(const std::string& filepath) {
	QuartetLookupTable<size_t> table;
	std::ifstream infile(filepath);
	if (!infile.good()) {
		throw std::runtime_error("Error reading file: " + filepath);
	}
	size_t n;
	infile >> n;
	table.init(n);
	std::string line;
	for (size_t i = 0; i < n; ++i) { // ignore the taxon labels
		infile >> line;
	}
	while (infile.good()) {
		infile >> line;
		parseTopologyString(table, line);
	}
	return table;
}

void writeQuartets(const QuartetLookupTable<size_t>& table, const std::vector<std::string>& taxonLabels,
		const std::string& filepath) {
	std::ofstream outfile(filepath);
	if (!outfile.good()) {
		throw std::runtime_error("Error writing file: " + filepath);
	}
	outfile << makeHeader(taxonLabels);
	for (size_t t1 = 0; t1 < taxonLabels.size() - 3; t1++) {
		for (size_t t2 = t1 + 1; t2 < taxonLabels.size() - 2; t2++) {
			for (size_t t3 = t2 + 1; t3 < taxonLabels.size() - 1; t3++) {
				for (size_t t4 = t3 + 1; t4 < taxonLabels.size(); t4++) {
					std::array < size_t, 3 > counts;
					auto& tuple = table.get_tuple(t1, t2, t3, t4);
					size_t t1t2t3t4 = table.tuple_index(t1, t2, t3, t4);
					size_t t1t3t2t4 = table.tuple_index(t1, t3, t2, t4);
					size_t t1t4t2t3 = table.tuple_index(t1, t4, t2, t3);
					counts[0] = tuple[t1t2t3t4];
					counts[1] = tuple[t1t3t2t4];
					counts[2] = tuple[t1t4t2t3];
					outfile << makeTopologyString(t1, t2, t3, t4, counts);
				}
			}
		}
	}
	outfile.close();
}
