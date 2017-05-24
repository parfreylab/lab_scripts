import sys
import re
import argparse
import os.path
import numpy as np


####a script that applies SILVA accessions to matching OTUs, after closed reference OTU picking for a MED output####
#IMPORTANT: this script will collapse nodes that match the same accession and add the read counts sample-wise to produce a new 'node' that represents all reads matching an accession in each sample
#it's important to make sure your closed-reference OTU picking that takes place before running this script is as stringent as you would like it to be, or you will lose a lot of detail/granularity of OTUs present in your original MED output.

##inputs:
###1. MED output MATRIX-COUNTS.txt
###2. Accession/OTU map from closed ref OTU picking output

##output:
###1. MED MATRIX-COUNT file, with inherited accessions

#########
##USAGE##
#########
#python /path/to/inherit_accessions.py -i MATRIX-COUNT.transposed.txt -m NODE-REPRESENTATIVES.DOWNSTREAM_otus.txt -o MATRIX-COUNT.transposed.inherited_accessions.txt

print "run: python inherit_accessions.py -h for help."
print "\n"

#Description of arguments needed
parser = argparse.ArgumentParser(
	description = "Two inputs required: MATRIX-COUNT.transposed.txt NODE-REPRESENTATIVES.DOWNSTREAM_otus.txt") 
requiredargs = parser.add_argument_group("required arguments")
requiredargs.add_argument( #matrix file. user defined.
	"-i",
	"--matrix_count",
	help = "Path to raw MED output file MATRIX-COUNT.txt",
	required = True)
requiredargs.add_argument( #path to otu map. user defined.
	"-m",
	"--otu_map",
	help = "The OTU map produced from closed-reference OTU picking with MED nodes against the SILVA database.",
	required = True)
parser.add_argument( #output file name. user defined. optional.
	"-o",
	"--outputfile",
	help = "Name of output file. Can contain full path to file. Default is in the same directory as the input, labeled MATRIX-COUNT.inherited_accessions.txt")
args = parser.parse_args()


matrix_count = args.matrix_count
otu_mapfile = args.otu_map
outputfile = args.outputfile

#read OTU map
#no header
otu_map = {}
with open(otu_mapfile) as MAP: #open barcode file
        for i, line in enumerate(MAP): #for each line
                data = line.rstrip() #strip whitespace
                all_in_line = data.split('\t') #split on tab
                accession = all_in_line[0] #first element
                OTUs = all_in_line[1:] #second through Nth elements are OTU IDs
		otu_map[accession] = OTUs
MAP.close()


#read matrix file. OTUs in rows, samples in columns
matrix = {}
arraylen = ()
header = []
with open (matrix_count) as MATRIXFILE: #open matrix count file
	for i, line in enumerate(MATRIXFILE): #for each line
                if i == 0:
                        data = line.rstrip()
                        header = data.split('\t')
                else:
	                data = line.rstrip() #strip whitespace
	                all = data.split('\t') #split on tab
	                OTU = all[0] #otu ID is the rowname
	                counts = all[1:] #read counts are the rest of each line
                        counts = map(int, counts)
                        arraylen=len(counts)
	                matrix[OTU] = counts #link counts with OTUs
MATRIXFILE.close()

#now merge OTUs that belong to the same accession
try:
        outputfile2=outputfile[:-4] #remove extension (if it is a 3 letter extension. This file will always be a .txt file so I don't think we need to worry here.)
        completename=os.path.join(outputfile2 + ".inherited_accessions.txt")
        OUTFILE=open(completename, "w")
except AttributeError:
        completename=os.path.join(matrix_count + ".inherited_accessions.txt")
        OUTFILE=open(completename, "w")
accession_counts = {}
otus_merged = {}
for key in sorted(otu_map): #for each accession
        otus = otu_map[key] #.split('\t') #retrieve OTU string from hash, split on tab to make list
        merged = [0] * arraylen #declare empty list to store merged read counts
        merged = map(int, merged)
        merged = np.array(merged)
        for ID in otus: #for each OTU ID belonging to an accession
                otus_merged[ID]=1 #save OTU ID in dictionary. test this hash later to decide which OTUs to include in final OTU map
                counts=np.array(matrix[ID]) #.split('\t') #collect read counts as individual elements of a list
                #tmp=map(add, int(counts), int(merged))
                #tmp=[sum(x) for x in zip(counts,merged)] #store read counts in the 'merged' list, adding the new counts for each OTU we are merging
                tmp = counts + merged
                merged = tmp.tolist()
                #merged = [sum(x) for x in izip_longest(counts, merged, fillvalue=0)]
                #sys.stderr.write(ID + "\t" + '\t'.join(merged) + "\n")
        accession_counts[key]=merged #store merged info in dictionary (hash)
        merged=[] #clear merged list


#and print accessions and unmerged OTUs to the outfile
headerprint = "\t".join(header)
OUTFILE.write(headerprint + "\n")
for key in sorted(accession_counts): #for each accession
        toprint='\t'.join(str(x) for x in accession_counts[key]) #must coerce to str from int (needed int to do math with numpy above)
#        toprint="\t".join(accession_counts[key]) #make string from contents of dictionary[list]
        OUTFILE.write(key + toprint + "\n") #print read counts to file

for key in sorted(matrix): #for OTU IDs
        if key not in otus_merged: #if we didn't merge the OTU into an accession
                toprint='\t'.join(str(x) for x in matrix[key])
                #toprint="\t".join(matrix[key]) #collect read counts into string
                OUTFILE.write(key + toprint + "\n") #print read coutnts to file
                
        else: #otherwise don't do anything
             pass
