"""
This script does the following:
1. remove extra information in the headers of FASTA files, and to format the rest for MED (BEFORE MED)
2. remove leading zeroes in transposed matrix count files or in the headers of FASTA files, or both simultaneously (AFTER MED)

inputs:
1. (optional) a trimmed FASTA file
2. (optional) a transposed matrix count .txt file 
3. (optional) a node representative sequences file in .fasta format
4. (optional) path to output folder

outputs:
1. a MED analysis ready FASTA file with headers properly formatted
2. the transposed matrix count file and/or node representative sequences file with leading zeroes removed from MED output

example usage (before MED): python trimheaders.rm_lead_zeroes.py -f seqs.trimmed_filtered_250bp.fna -o /path/to/output/folder
example usage (after MED):  python trimheaders.rm_lead.zeroes.py -m MATRIX-COUNT.txt -n NODE-REPRESENTATIVES.fasta -o /path/to/output/folder
"""

import argparse
import sys
import os
import re
import numpy as np

__author__ = "Kevin Chan"
__email__ = "kevchan1@alumni.ubc.ca"

parser = argparse.ArgumentParser(
	description = "This script removes extraneous information in the headers of FASTA files for MED input, and transposes + removes leading zeroes from MED output.")
parser.add_argument(
	"-f",
	"--trimmed_file",
	help = "Path to trimmed .fna file, string.",)
parser.add_argument(
	"-m",
	"--matrixfile",
	help = "Path to the matrix count .txt file, string.")
parser.add_argument(
	"-n",
	"--noderepfile",
	help = "Path to the node representative sequences .fasta file, string.")
parser.add_argument(
	"-o",
	"--outputpath",
	default = "./",
	help = "Path to output directory, string. Default: current working directory.")

args = parser.parse_args()

# if nothing is specified, display a message
if not len(sys.argv) > 1:
	print "Please enter an argument!"

trimmed_file = args.trimmed_file
matrixfile = args.matrixfile
noderepfile = args.noderepfile
out = args.outputpath

# if the trimmed FASTA file is specified
if trimmed_file is not None:
	# get the basename of the file
	trimmed_file_base = os.path.basename(trimmed_file)
	trimmed_file_base = os.path.splitext(trimmed_file_base)[0]

	# join the output directory with the file to write. by default, output directory is current working directory
	completeName = os.path.join(out, trimmed_file_base + ".MED.fna")
	trim_out = open(completeName, "w")

	with open(trimmed_file, "U") as f:
		for i, line in enumerate(f):
			if line.startswith(">"): 		   # if the line is a read header
				line = line.strip()			   # remove leading and trailing whitespace
				header_info = line.split(" ")  # split the line by spaces
				sample_info = header_info[0]   # keep the sample names
				if sample_info.count('_') > 1: # if the sample names are delimited by underscores
					headertoprint = re.split(r"(\w*_\w*)+_", sample_info) # match based on the last underscore before read number
					readnum = headertoprint[-1]
					del headertoprint[-1] # delete the read number from the list
					headertoprint = "".join(headertoprint) # join all elements of the list 
					headertoprint = headertoprint + "_Read" + readnum # add in the read number and properly format
					trim_out.write(headertoprint + "\n")
				else: # if sample names are not delimited by underscores
					sample_info = sample_info.split('_')
					sample_info = "_Read".join(sample_info)
					trim_out.write(sample_info + "\n")
			else: # if the line contains a sequence
				trim_out.write(line)
	print "successfully trimmed FASTA headers!"

# if the matrix count file is specified
if matrixfile is not None:
	# getting filename of the input file without the extension
	matrixfile_base = os.path.basename(matrixfile)
	matrixfile_base = os.path.splitext(matrixfile_base)[0]

	matrix_path = os.path.join(out, matrixfile_base + ".transposed.txt")
	matrix_out = open(matrix_path, "w")

	to_transpose = [] # initialize a list which will contain line information from matrix-count.txt
	with open(matrixfile, "U") as f:
		for i, line in enumerate(f):
			data = line.strip() # strip forward and backward whitespace
			data_info = data.split("\t") # split each line into a list by tabs
			to_transpose.append(data_info) # append each list (data info) to another list (to_transpose). creates a 2D list for transposition
	to_transpose = np.array(to_transpose).transpose() # change from type list to array, so that transpose can be called
	counter = 0 # counter for iterating through the 2D array
	for i in to_transpose:
		if counter == 0: 
			pass # the first element contains sample info, no leading zeroes so do nothing
		else:
			match = re.match(r"(0+)(\d*)", i[0]) # match based on leading zeroes in first element of each element of array
			if match:
				zeroes_rmed = match.groups()[1]  # get the string of numbers with leading zeroes removed
				to_transpose[counter][0] = zeroes_rmed # replace with leading zeroes removed
		counter += 1
	np.savetxt(matrix_out, to_transpose, fmt = '%s', delimiter = "\t") # save the matrix as a tab-delimited file
	print "successfully transposed and removed leading zeroes from matrix file!"

# if the node representatives file is specified
if noderepfile is not None:
	noderepfile_base = os.path.basename(noderepfile)
	noderepfile_base = os.path.splitext(noderepfile_base)[0]

	node_path = os.path.join(out, noderepfile_base + ".DOWNSTREAM.fasta")
	node_out = open(node_path, "w")
	
	with open(noderepfile, "U") as f:
		for i, line in enumerate(f):
			if line.startswith(">"): # if the line starts with '>', it is a header
				line = line.strip()
				data = line.split("|") # split by '|', create a list (>SAMPLEID|size: XXXX)
				match = re.match(r"(>0+)(\d*)", data[0]) # match the first element of the list containing sample IDs
				if match:
					sample_id_format = match.groups()[1] # the sample ID with leading zeroes removed
					node_out.write(">" + sample_id_format + "\n") # write the header line with formatted sample IDs out
			else: # if the line is a sequence, copy and paste!
				node_out.write(line)			
	print "successfully removed leading zeroes from node rep file!"		
