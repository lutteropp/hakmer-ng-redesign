#! /usr/env/python

import os
import sys
import math
from sets import Set
import os.path
import random
import shlex
from subprocess import check_output, check_call

from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
import sys

# function taken from the multiSPAM tool
def fix_tree(tree_file, id_map, out_file):
	import re
	with open(tree_file) as file, open(out_file, "w+") as out:
		for line in file:
			parts = re.split('(,|\(|\)|;)',line)
			for idx, p in enumerate(parts):
				if p.isdigit():
					value = id_map[int(p)]
					parts[idx] = value
			out.write(''.join(parts))

def readIDMap(lines):
	n = int(lines[0])
	id_map = []
	for i in range(n):
		id_map.append(lines[i + 1].split(" ")[0])
	return id_map

def make_tree(quartets_path, tree_path):
	print_all = False
	print_with_val = False
	with_weights = False
	print_counts = False

	lines = open(sys.argv[quartets_path]).readlines()
	id_map = readIDMap(lines)

	tmp_file = open("tempFile.txt", 'w')
	n = len(id_map)
	for i in range(n + 1, len(lines), 3):
		counts = []
		tops = []
		newline = ""
		splitted = lines[i].split(",")
		newline = splitted[0] + "," + splitted[1] + "," + splitted[2].split(":")[0]
		count = int(splitted[2].split(":")[1])
		counts.append(count)
		tops.append(newline)

		newline = ""
		splitted = lines[i+1].split(",")
		newline = splitted[0] + "," + splitted[1] + "," + splitted[2].split(":")[0]
		count = int(splitted[2].split(":")[1])
		counts.append(count)
		tops.append(newline)

		newline = ""
		splitted = lines[i+2].split(",")
		newline = splitted[0] + "," + splitted[1] + "," + splitted[2].split(":")[0]
		count = int(splitted[2].split(":")[1])
		counts.append(count)
		tops.append(newline)

		best = 0
		if counts[1] >= counts[0] and counts[1] >= counts[2]:
			best = 1
		if counts[2] >= counts[0] and counts[2] >= counts[1]:
			best = 2
		vals = []
		for j in range(3):
			if counts[j] == 0:
				vals.append(0)
			else:
				if print_counts:
					vals.append(counts[j])
				else:
					vals.append(float(counts[j]) / (counts[0] + counts[1] + counts[2]))
		
		if print_all:
			if (counts[0] > 0):
				if print_with_val:
					tmp_file.write(tops[0] + ":" + str(vals[0]) + "\n")
				else:
					tmp_file.write(tops[0] + "\n")
			if (counts[1] > 0):
				if print_with_val:
					tmp_file.write(tops[1] + ":" + str(vals[1]) + "\n")
				else:
					tmp_file.write(tops[1] + "\n")
			if (counts[2] > 0):
				if print_with_val:
					tmp_file.write(tops[2] + ":" + str(vals[2]) + "\n")
				else:
					tmp_file.write(tops[2] + "\n")
		elif counts[best] > 0:
			for j in range(3):
				if counts[j] == counts[best]:
					if print_with_val:
						tmp_file.write(tops[j] + ":" + str(vals[j]) + "\n")
					else:
						tmp_file.write(tops[j] + "\n")

	id_tree_file = "idTree.txt"
	if with_weights:
		check_call(["./max-cut-tree", "qrtt=tempFile.txt", "weights=on", "otre=" + id_tree_file])
	else:
		check_call(["./max-cut-tree", "qrtt=tempFile.txt", "weights=off", "otre=" + id_tree_file])		

	fix_tree(id_tree_file, id_map, tree_path)
	os.remove("tempFile.txt")
	os.remove(id_tree_file)

def compute_norm_rf_score(reference_tree_path, sampled_tree_path):
	out = check_output(["Rscript", "--vanilla", "rf_dist.r", reference_tree_path, sampled_tree_path])
	splitted = out.split(" ")
	return float(splitted[len(splitted) - 1])

def nexus_to_fasta(input_name, output_name):
	input_file = open(input_name, 'r')
	output_file = open(output_name, 'w')
	alignment = AlignIO.read(input_file, 'nexus')
	AlignIO.write(alignment, output_file, 'fasta')
	input_file.close()
	output_file.close()

def create_param_str(kmin, revcomp, largeseeds, flankwidth, flavor):
	param_str = "_kmin" + str(kmin) + "_revcomp" + str(revcomp) + "_largeseeds" + str(largeseeds) + "_flank" + str(flankwidth) + "_flavor" + str(flavor) + "_"
	return param_str

def do_experiment(kmin, revcomp, largeseeds, flankwidth, flavor, sequences_inpath, reference_tree_path):
	param_str = create_param_str(kmin, revcomp, largeseeds, flankwidth, flavor)
	hakmer_call = "./hakmer-ng-redesign -f " + sequences_inpath + " --kmin " + str(kmin)
	if flankwidth > 0:
		hakmer_call = hakmer_call + " --flankwidth " + str(flankwidth)
	else:
		hakmer_call = hakmer_call + " --dynamic"
	if revcomp:
		hakmer_call = hakmer_call + " --revcomp"
	if largeseeds:
		hakmer_call = hakmer_call + " --largeSeeds"

	supermatrix_path = sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + param_str
	hakmer_call = hakmer_call + " -o " + "hakmer_quartets/" + supermatrix_path

	hakmer_call = hakmer_call + " quartets"

	if flavor == 1:
		hakmer_call = hakmer_call + " --concatDist"
	elif flavor == 2:
		hakmer_call = hakmer_call + " --concatMSA"		

	print "calling hakmer-ng quartets..."
	print hakmer_call

	check_call(shlex.split(hakmer_call))
	score = -1
	try:
		# Run the tree inference via quartet max cut
		make_tree("hakmer_quartets/" + supermatrix_path, "hakmer_quartets/" + supermatrix_path + ".tree")
		sampled_tree_path = "hakmer_quartets/" + supermatrix_path + ".tree"
		score = compute_norm_rf_score(reference_tree_path, sampled_tree_path)
	except:
		pass
	return score

if __name__== "__main__":
	# Parse command line arguments
	reference_tree_path = ""
	sequences_inpath = ""
	if (len(sys.argv) == 3):
		reference_tree_path = sys.argv[1]
		sequences_inpath = sys.argv[2]
	else:
		print "Usage: python orig_hakmer_experiments.py PATH_TO_REFERENCE_TREE PATH_TO_SEQUENCES"
		sys.exit()
	if not os.path.exists("hakmer_quartets/"):
    		os.makedirs("hakmer_quartets/")

	kmin_values = [8, 16, 32, 64]
	revcomp_values = [False, True]
	flankwidth_values = [-1, 25, 50, 100, 200, 400]
	largeseeds_values = [False, True]
	flavor_values = [0,1,2]

	exp_outpath = sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + "_results.csv"

	if not os.path.exists(exp_outpath):
		exp_outfile = open(exp_outpath, 'w')
		exp_outfile.write("kmin;revComp;largeseeds;flanks;flavor;norm_rf\n")
		exp_outfile.close()

	for kmin in kmin_values:
		for revcomp in revcomp_values:
			for largeseeds in largeseeds_values:
				for flankwidth in flankwidth_values:
					for flavor in flavor_values:
						param_str = create_param_str(kmin, revcomp, largeseeds, flankwidth, flavor)
						supermatrix_path = "hakmer_quartets/" + sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + param_str
						if os.path.exists(supermatrix_path):
							print "Skipping " + supermatrix_path
							continue
						score = do_experiment(kmin, revcomp, largeseeds, flankwidth, flavor, sequences_inpath, reference_tree_path)
						exp_outfile = open(exp_outpath, 'a')
						exp_outfile.write(str(kmin) + ";" + str(revcomp) + ";" + str(largeseeds) + ";" + str(flankwidth) + ";" + str(flavor) + ";" + str(score) + "\n")
						exp_outfile.close()
						print "Added " + supermatrix_path

