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

def create_param_str(kmin, revcomp, c, largeseeds, flankwidth, quickdelta):
	param_str = "_kmin" + str(kmin) + "_revcomp" + str(revcomp) + "_c" + str(c) + "_largeseeds" + str(largeseeds) + "_flank" + str(flankwidth) + "_quickdelta" + str(quickdelta) + "_"
	return param_str

def do_experiment(kmin, revcomp, c, largeseeds, flankwidth, quickdelta, sequences_inpath, reference_tree_path):
	param_str = create_param_str(kmin, revcomp, c, largeseeds, flankwidth, quickdelta)
	hakmer_call = "./hakmer-ng-redesign -f " + sequences_inpath + " --kmin " + str(kmin)
	if flankwidth > 0:
		hakmer_call = hakmer_call + " --flankwidth " + str(flankwidth)
	else:
		hakmer_call = hakmer_call + " --dynamic"
	if revcomp:
		hakmer_call = hakmer_call + " --revcomp"
	if largeseeds:
		hakmer_call = hakmer_call + " --largeSeeds"
	if quickdelta:
		hakmer_call = hakmer_call + " --quickDelta"

	supermatrix_path = sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + param_str
	hakmer_call = hakmer_call + " -o " + "hakmer_fasta/" + supermatrix_path

	hakmer_call = hakmer_call + " matrix" + " --minTaxa " + str(c)

	print "calling hakmer-ng supermatrix..."
	print hakmer_call

	check_call(shlex.split(hakmer_call))
	# now we need to convert the supermatrix into a FASTA file
	score = -1
	try:
		# Run the tree inference via raxml-ng
		raxml_call = "./raxml-ng --msa " + "hakmer_fasta/" + supermatrix_path + " --model GTR+G" + " --tree rand{5},pars{5} --force"
		print "calling raxml-ng..."
		check_call(shlex.split(raxml_call))
		sampled_tree_path = "hakmer_fasta/" + supermatrix_path + ".raxml.bestTree"
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
	if not os.path.exists("hakmer_fasta/"):
    		os.makedirs("hakmer_fasta/")

	kmin_values = [8, 16, 32, 64]
	c_values = [4, 8]
	revcomp_values = [False, True]
	flankwidth_values = [-1, 25, 50, 100, 200, 400]
	largeseeds_values = [False, True]
	quickdelta_values = [True, False]

	exp_outpath = sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + "_results.csv"

	if not os.path.exists(exp_outpath):
		exp_outfile = open(exp_outpath, 'w')
		exp_outfile.write("kmin;revComp;c;largeseeds;flanks;quickdelta;norm_rf\n")
		exp_outfile.close()

	for kmin in kmin_values:
		for revcomp in revcomp_values:
			for c in c_values:
				for largeseeds in largeseeds_values:
					for flankwidth in flankwidth_values:
						for quickdelta in quickdelta_values:
							if flankwidth != -1 and quickdelta == True:
								continue
							param_str = create_param_str(kmin, revcomp, c, largeseeds, flankwidth, quickdelta)
							supermatrix_path = "hakmer_fasta/" + sequences_inpath.split("/")[len(sequences_inpath.split("/")) - 1] + param_str
							if os.path.exists(supermatrix_path):
								print "Skipping " + supermatrix_path
								continue
							score = do_experiment(kmin, revcomp, c, largeseeds, flankwidth, quickdelta, sequences_inpath, reference_tree_path)
							exp_outfile = open(exp_outpath, 'a')
							exp_outfile.write(str(kmin) + ";" + str(revcomp) + ";" + str(c) + ";" + str(largeseeds) + ";" + str(flankwidth) + ";" + str(quickdelta) + ";" + str(score) + "\n")
							exp_outfile.close()
							print "Added " + supermatrix_path

