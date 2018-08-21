#usr/bin/env python

# author: Laura Wegener Parfrey
# date: Dec 2014
# goal: substitute a string for another string. 
# Input is a tab delimited file with new_name \t old_name (common usage is gg_accession \t MEDtype_number) and file to be substituted. 
# usage: python substitute.py names_file.txt file_to_change.txt > out_file.txt
# note that this script is crude, so need to change print line to correct format or modify input names file to contain whole replacement line. 

from sys import argv

#open tab delimited names file
names = open(argv[1], "U")

#open file to be substituted
file_sub = open(argv[2], "U")

# create dict of old and new names
names_dict = {}
for line in names:
	name = line.strip().split('\t')
	names_dict[name[1]] = name[0]
	
# 	to substitute into excel format OTU table
# for line in file_sub:
# 	acc = line.strip().split('\t')
# 	# test for presence of old name
# 	if acc[0] in names_dict:
# 		# replace old name with new name
# 		new_name = names_dict[acc[0]]
# 		# change print line to print in correct format
# 		#updated_line = 
# 		print "%s\t%s" % (new_name, '\t'.join(acc[1:-1]))
# 	else:
# 		print line.strip()
		

# to substitute into fasta header	
from cogent.parse.fasta import MinimalFastaParser
	
for label, seq in MinimalFastaParser(file_sub):
	l = label.strip('>').strip().split(' ')
	curr_label = l[0]
	#test for accession in subtype map
	if curr_label in names_dict:
		# print lable in fasta for this accession
		print ">%s\n%s" % (names_dict[curr_label], seq)
	else: 
		# print out header and seq as is - seq not in rename file
		print ">%s\n%s" % (curr_label, seq)
