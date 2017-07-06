"""
A script to remove extra information in the headers of FASTA files, and to format the rest for MED.

inputs:
1. a trimmed FASTA file
2. (optional) path to output folder

output:
a MED analysis ready FASTA file with headers properly formatted

example usage: python trimheaders.add_read.py -f seqs.fna -o /path/to/output/folder
"""

import argparse
import os
import re

parser = argparse.ArgumentParser(
	description = "One input required: trimmed_data.fna")
requiredargs = parser.add_argument_group("required argument")
requiredargs.add_argument(
	"-f",
	"--trimmed_file",
	help = "Path to trimmed .fna file, string.",
	required = True)
parser.add_argument(
	"-o",
	"--outputpath",
	help = "Path to output directory, string. Default: current working directory.")

args = parser.parse_args()

trimmed_file = args.trimmed_file
out = args.outputpath

# get the basename of the file
trimmed_file_base = os.path.basename(trimmed_file)
trimmed_file_base = os.path.splitext(trimmed_file_base)[0]

try: # try to create an outfile in the specified directory
	completeName = os.path.join(out, trimmed_file_base + ".MED.fna")
	trim_out = open(completeName, "w")
except AttributeError: # create an outfile in the current directory
	trim_out = open((trimmed_file_base + ".MED.fna"), "w")

with open(trimmed_file) as f:
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