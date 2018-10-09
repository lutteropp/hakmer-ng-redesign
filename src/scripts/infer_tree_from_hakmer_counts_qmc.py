#!usr/env/python

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details at
# http://www.gnu.org/copyleft/gpl.html

import math

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

if __name__ == "__main__":
	import sys
	import subprocess as sp
	import os

	print_all = False
	print_with_val = False
	with_weights = False
	print_counts = False

	lines = open(sys.argv[1]).readlines()
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
		sp.check_call(["./max-cut-tree", "qrtt=tempFile.txt", "weights=on", "otre=" + id_tree_file])
	else:
		sp.check_call(["./max-cut-tree", "qrtt=tempFile.txt", "weights=off", "otre=" + id_tree_file])		

	fix_tree(id_tree_file, id_map, sys.argv[2])

	os.remove("tempFile.txt")
	os.remove(id_tree_file)
