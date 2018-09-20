#include "mafft_raxml_wrapper.hpp"
#include "io.hpp"

#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <unordered_set>

size_t countUniqueSequences(const std::array<std::string, 4>& sequences) {
	std::unordered_set<std::string> s;
	for (size_t i = 0; i < sequences.size(); ++i) {
		s.insert(sequences[i]);
	}
	return s.size();
}

std::vector<std::string> mafftAlign(const std::string& prefix, const std::vector<std::string>& sequences, const std::vector<std::string>& taxNames) {
	if (sequences[0].empty()) {
		return sequences;
	}

	std::vector<std::string> res;
	std::string filename = prefix + "_alignment.fasta";

	// write all sequences into a file
	std::ofstream outfile(filename);
	for (size_t i = 0; i < sequences.size(); ++i) {
		outfile << ">" << taxNames[i] << "\n";
		outfile << sequences[i] << "\n";
	}
	outfile.close();

	std::string alignedFilename = filename + ".aligned";

	// run mafft for aligning the sequences
	std::string mafftCall = "mafft-ginsi --quiet " + filename + " > " + alignedFilename;
	//std::string mafftCall = "mafft --retree 1 --quiet " + filename + " > " + alignedFilename;

	int status = std::system(mafftCall.c_str());
	if (!WIFEXITED(status)) {
		throw std::runtime_error("Something went wrong while running MAFFT!");
	}

	// read back the sequences
	std::ifstream infile(alignedFilename);
	for (size_t i = 0; i < sequences.size(); ++i) {
		FASTARecord rec = readNextFASTA(infile);
		res[i] = rec.seq;
		//std::transform(res[i].begin(), res[i].end(), res[i].begin(), ::toupper);
	}

	// remove the temporary files again
	std::remove(filename.c_str());
	std::remove(alignedFilename.c_str());

	return res;
}

std::array<std::string, 4> mafftAlign(const std::string& prefix, const std::array<std::string, 4>& sequences) {
	if (sequences[0].empty()) {
		return sequences;
	}

	std::array<std::string, 4> res;

	// TODO: randomly generate a file name
	std::string filename = prefix + "randomFilename.fasta";

	std::array<std::string, 4> taxNames = { "a", "b", "c", "d" };
	// write all sequences into a file
	std::ofstream outfile(filename);
	for (size_t i = 0; i < 4; ++i) {
		outfile << ">" << taxNames[i] << "\n";
		outfile << sequences[i] << "\n";
	}
	outfile.close();

	std::string alignedFilename = filename + ".aligned";

	// run mafft for aligning the sequences
	std::string mafftCall = "mafft-ginsi --quiet " + filename + " > " + alignedFilename;

	//std::string mafftCall = "mafft --retree 1 --quiet " + filename + " > " + alignedFilename;

	int status = std::system(mafftCall.c_str());
	if (!WIFEXITED(status)) {
		throw std::runtime_error("Something went wrong while running MAFFT!");
	}

	// read back the sequences
	std::ifstream infile(alignedFilename);
	for (size_t i = 0; i < sequences.size(); ++i) {
		FASTARecord rec = readNextFASTA(infile);
		res[i] = rec.seq;
		//std::transform(res[i].begin(), res[i].end(), res[i].begin(), ::toupper);
	}

	// remove the temporary files again
	std::remove(filename.c_str());
	std::remove(alignedFilename.c_str());

	return res;
}

std::array<double, 3> retrieveLogl(const std::string& logPath, const Options& options) {
	std::array<double, 3> res;
	std::ifstream infile(logPath);
	if (!infile.good()) {
		throw std::runtime_error("File does not exist: " + logPath);
	}

	std::string line;
	while (std::getline(infile, line)) {
		if (line.find("Evaluating 3 trees") != std::string::npos) {
			break;
		}
	}
	std::getline(infile, line);

	std::getline(infile, line);
	// parse logl of tree #1
	line = line.substr(line.find(": ") + 2, std::string::npos);
	res[0] = std::stod(line);
	std::getline(infile, line);
	std::getline(infile, line);
	// parse logl of tree #2
	line = line.substr(line.find(": ") + 2, std::string::npos);
	res[1] = std::stod(line);
	std::getline(infile, line);
	std::getline(infile, line);
	// parse logl of tree #3
	line = line.substr(line.find(": ") + 2, std::string::npos);
	res[2] = std::stod(line);

	return res;
}

QuartetTopology inferTopologyPlaceholder(const std::string& prefix, const std::vector<std::array<std::string, 4> >& alignedSequences,
		const Options& options) {
	return QuartetTopology::STAR;
}

QuartetTopology inferTopology(const std::string& prefix, const std::vector<std::array<std::string, 4> >& alignedSequences, const Options& options) {
	if (alignedSequences.size() == 0 || alignedSequences[0][0].size() == 0) {
#pragma omp critical
		std::cout << "WARNING: Could not find any blocks for quartet " << prefix << "\n";
		return QuartetTopology::STAR;
	}

	std::vector<size_t> partitionSizes;
	for (size_t i = 0; i < alignedSequences.size(); ++i) {
		partitionSizes.push_back(alignedSequences[i][0].size());
	}

	std::string msaFilepath = prefix + "_msa.txt";
	std::ofstream msaFile(msaFilepath);
	// write the MSA into a file
	std::array<std::string, 4> taxNames = { "a", "b", "c", "d" };
	std::array<std::string, 4> concatenatedSequences;
	for (size_t i = 0; i < 4; ++i) {
		msaFile << ">" << taxNames[i] << "\n";
		for (size_t j = 0; j < alignedSequences.size(); ++j) {
			msaFile << alignedSequences[j][i];
			concatenatedSequences[i] += alignedSequences[j][i];
		}
		msaFile << "\n";
	}
	msaFile.close();

	size_t numUnique = countUniqueSequences(concatenatedSequences);
	if (numUnique == 2) { // special case
		if (concatenatedSequences[0] == concatenatedSequences[1] && concatenatedSequences[2] == concatenatedSequences[3]) { // ab|cd
			return QuartetTopology::AB_CD;
		} else if (concatenatedSequences[0] == concatenatedSequences[2] && concatenatedSequences[1] == concatenatedSequences[3]) { // ac|bd
			return QuartetTopology::AC_BD;
		} else if (concatenatedSequences[0] == concatenatedSequences[3] && concatenatedSequences[1] == concatenatedSequences[2]) { // ad|bc
			return QuartetTopology::AD_BC;
		} else { // star topology
			return QuartetTopology::STAR;
		}
	}

	if (numUnique < 3) {
#pragma omp critical
		std::cout << "WARNING: Too many identical sequences for quartet " << prefix << "\n";
		return QuartetTopology::STAR;
	}

	std::string modelString;
	std::string partitionsFilepath = prefix + "_part.txt";
	if (!options.noPartitions) {
		std::ofstream partitionsFile(partitionsFilepath);
		size_t lastEndCoord = 0;
		for (size_t i = 0; i < alignedSequences.size(); ++i) {
			size_t startCoord = lastEndCoord + 1;
			size_t endCoord = startCoord + partitionSizes[i] - 1;
			partitionsFile << "GTR+G, p" << i << " = " << startCoord << "-" << endCoord << "\n";
			lastEndCoord = endCoord;
		}
		partitionsFile.close();

		modelString = "--model " + partitionsFilepath;
	} else {
		modelString = "--model GTR+G";
	}
	std::string outpath = prefix + "_out.txt";

	std::string raxmlCall = "./raxml-ng --evaluate --msa " + msaFilepath + " " + modelString
			+ " --tree quartet_trees --threads 1 --log INFO --nofiles > " + outpath;
	int status = std::system(raxmlCall.c_str());
	if (!WIFEXITED(status)) {
		throw std::runtime_error("Something went wrong while running RAxML!");
	}

	std::array<double, 3> logl = retrieveLogl(outpath, options);

	// read back the topology
	QuartetTopology topology = QuartetTopology::STAR;
	double maxLogl = std::max(logl[0], std::max(logl[1], logl[2]));
	double minLogl = std::min(logl[0], std::min(logl[1], logl[2]));
	if (minLogl != maxLogl) {
		if (logl[0] == maxLogl) {
			topology = QuartetTopology::AB_CD;
		} else if (logl[1] == maxLogl) {
			topology = QuartetTopology::AC_BD;
		} else {
			topology = QuartetTopology::AD_BC;
		}
		size_t countMaxLogl = 0;
		if (logl[0] == maxLogl)
			countMaxLogl++;
		if (logl[1] == maxLogl)
			countMaxLogl++;
		if (logl[2] == maxLogl)
			countMaxLogl++;
		if (countMaxLogl > 1) {
			topology = QuartetTopology::STAR;
		}
	}

	// remove the temporary files again
	std::remove(msaFilepath.c_str());
	if (!options.noPartitions) {
		std::remove(partitionsFilepath.c_str());
	}
	std::remove(outpath.c_str());

	return topology;
}

QuartetTopology inferTopology(const std::string& prefix, const std::array<std::string, 4>& alignedConcatenatedSequences, const Options& options) {
	if (alignedConcatenatedSequences[0].size() == 0) {
#pragma omp critical
		std::cout << "WARNING: Could not find any blocks for quartet " << prefix << "\n";
		return QuartetTopology::STAR;
	}

	std::string msaFilepath = prefix + "_msa.txt";
	std::ofstream msaFile(msaFilepath);
	// write the MSA into a file
	std::array<std::string, 4> taxNames = { "a", "b", "c", "d" };
	for (size_t i = 0; i < 4; ++i) {
		msaFile << ">" << taxNames[i] << "\n";
		msaFile << alignedConcatenatedSequences[i] << "\n";
	}
	msaFile.close();

	size_t numUnique = countUniqueSequences(alignedConcatenatedSequences);
	if (numUnique == 2) { // special case
		if (alignedConcatenatedSequences[0] == alignedConcatenatedSequences[1]
				&& alignedConcatenatedSequences[2] == alignedConcatenatedSequences[3]) { // ab|cd
			return QuartetTopology::AB_CD;
		} else if (alignedConcatenatedSequences[0] == alignedConcatenatedSequences[2]
				&& alignedConcatenatedSequences[1] == alignedConcatenatedSequences[3]) { // ac|bd
			return QuartetTopology::AC_BD;
		} else if (alignedConcatenatedSequences[0] == alignedConcatenatedSequences[3]
				&& alignedConcatenatedSequences[1] == alignedConcatenatedSequences[2]) { // ad|bc
			return QuartetTopology::AD_BC;
		} else { // star topology
			return QuartetTopology::STAR;
		}
	}

	if (numUnique < 3) {
#pragma omp critical
		std::cout << "WARNING: Too many identical sequences for quartet " << prefix << "\n";
		return QuartetTopology::STAR;
	}

	std::string modelString;
	modelString = "--model GTR+G";
	std::string outpath = prefix + "_out.txt";

	std::string raxmlCall = "./raxml-ng --evaluate --msa " + msaFilepath + " " + modelString
			+ " --tree quartet_trees --threads 1 --log INFO --nofiles > " + outpath;
	int status = std::system(raxmlCall.c_str());
	if (!WIFEXITED(status)) {
		throw std::runtime_error("Something went wrong while running RAxML!");
	}

	std::array<double, 3> logl = retrieveLogl(outpath, options);

	// read back the topology
	QuartetTopology topology = QuartetTopology::STAR;
	double maxLogl = std::max(logl[0], std::max(logl[1], logl[2]));
	double minLogl = std::min(logl[0], std::min(logl[1], logl[2]));
	if (minLogl != maxLogl) {
		if (logl[0] == maxLogl) {
			topology = QuartetTopology::AB_CD;
		} else if (logl[1] == maxLogl) {
			topology = QuartetTopology::AC_BD;
		} else {
			topology = QuartetTopology::AD_BC;
		}
		size_t countMaxLogl = 0;
		if (logl[0] == maxLogl)
			countMaxLogl++;
		if (logl[1] == maxLogl)
			countMaxLogl++;
		if (logl[2] == maxLogl)
			countMaxLogl++;
		if (countMaxLogl > 1) {
			topology = QuartetTopology::STAR;
		}
	}

	// remove the temporary files again
	std::remove(msaFilepath.c_str());
	std::remove(outpath.c_str());

	return topology;
}
