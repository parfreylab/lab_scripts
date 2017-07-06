"""
A script to remove leading zeroes in transposed matrix count files or in the headers of FASTA files, or both simultaneously.

inputs:
1. a transposed matrix count .txt file 
2. a node representative sequences file in .fasta format
3. (optional) path to output folder

outputs:
the transposed matrix count file and/or node representative sequences file with leading zeroes removed

example usage: python rm_lead_zeroes.py -m MATRIX-COUNT.txt -n NODE-REPRESENTATIVES.fasta -o /path/to/output/folder
"""

import argparse
import os
import re

parser = argparse.ArgumentParser(
	description = "This script removes leading zeroes from a transposed MATRIX-COUNT.txt table or from a NODE-REPRESENTATIVES.fasta file, or both simultaneously.")
parser.add_argument(
	"-m",
	"--matrixfile",
	help = "Path to the transposed matrix count .txt file, string.")
parser.add_argument(
	"-n",
	"--noderepfile",
	help = "Path to the node representative sequences .fasta file, string.")
parser.add_argument(
	"-o",
	"--outputpath",
	help = "Path to output directory, string.")

args = parser.parse_args()

matrixfile = args.matrixfile
noderepfile = args.noderepfile
out = args.outputpath

# if the matrix count file is specified
if matrixfile is not None:
	# getting filename of the input file without the extension
	matrixfile_base = os.path.basename(matrixfile)
	matrixfile_base = os.path.splitext(matrixfile_base)[0]
	try: # try to create an outfile in the specified directory
		matrix_path = os.path.join(out, matrixfile_base + ".transposed.txt")
		matrix_out = open(matrix_path, "w")
	except AttributeError: # create an outfile in the current directory
		matrix_out = open((matrixfile_base + ".transposed.txt"), "w")
	# read the input matrix file
	with open(matrixfile, "U") as f:
		for i, line in enumerate(f):
			data = line.strip()
			data_info = data.split("\t") # split each line into a list by tabs
			sample_id = data_info[0] # get sample IDs as the first element
			rest = data_info[1:] # get the rest of the line in a separate list
			if i == 0: # the first line of the file does not require any modification
				matrix_out.write(data + "\n")
			else:
				match = re.match(r"(0+)(\d*)", sample_id) # match based on leading zeroes in the sample IDs
				if match:
					sample_id_format = match.groups()[1] # the sample ID with all the leading zeroes removed
					rest.insert(0, sample_id_format) # insert the formatted sample ID as the first element of the list
					matrix_out.write("\t".join(rest) + "\n")

# if the node representatives file is specified
if noderepfile is not None:
	noderepfile_base = os.path.basename(noderepfile)
	noderepfile_base = os.path.splitext(noderepfile_base)[0]
	try:
		node_path = os.path.join(out, noderepfile_base + ".DOWNSTREAM.fasta")
		node_out = open(node_path, "w")
	except AttributeError:
		node_out = open((noderepfile_base + ".DOWNSTREAM.fasta"), "w")
	with open(noderepfile, "U") as f:
		for i, line in enumerate(f):
			if line.startswith(">"): # if the line starts with '>', it is a header
				line = line.strip()
				data = line.split("|") # split by '|', create a list (>SAMPLEID|size: XXXX)
				match = re.match(r"(>0+)(\d*)", data[0]) # match the first element of the list containing sample IDs
				if match:
					sample_id_format = match.groups()[1] # the sample ID with leading zeroes removed
					node_out.write(">" + sample_id_format + "|" + data[1] + "\n") # write the header line with formatted sample IDs out
			else: # if the line is a sequence, copy and paste!
				node_out.write(line)					


