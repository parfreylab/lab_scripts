import re
import argparse
import os.path

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
#python /path/to/inherit_accessions.py -i MATRIX-COUNT.txt -m NODE-REPRESENTATIVES.DOWNSTREAM_otus.txt -o MATRIX-COUNT.inherited_accessions.txt

print "run: python inherit_accessions.py -h for help."
print "\n"

#Description of arguments needed
parser = argparse.ArgumentParser(
	description = "Two inputs required: MATRIX-COUNT.txt NODE-REPRESENTATIVES.DOWNSTREAM_otus.txt") 
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
otu_map = args.otu_map
outputfile = args.outputfile

#read OTU map
#no header
otu_map = {}
with open(otu_map) as MAP: #open barcode file
        for i, line in enumerate(MAP): #for each line
                data = line.rstrip() #strip whitespace
                all_in_line = data.split('\t') #split on tab
                accession = all_in_line[0] #first element
                OTUs = all_in_line[1:] #second through Nth elements are OTU IDs
		otu_map[accession] = OTUs
MAP.close()


#read matrix file. OTUs in rows, samples in columns
matrix = {}
with open () as MATRIXFILE: #open matrix count file
	for i, line in enumerate(MATRIXFILE): #for each line
	data = line.rstrip() #strip whitespace
	all = data.split('\t') #split on tab
	OTU = all[0] #otu ID is the rowname
	counts = all[1:] #read counts are the rest of each line
	matrix[OTU] = counts #link counts with OTUs
MATRIXFILE.close()

#now do the operations

#for each accession
for key in sorted(otu_map):
        data=otu_map[key]
        otus=data.split('\t')
        counts=list()
        for ID in otus:
                counts=matrix[ID].split('\t')
                

                
#collect OTUs that belong to it

#merge them by indices? #there might be a cool easy way to do this in python. do some googling before trying to reinvent the wheel
#after merge print to file

#keep track of OTUs merged

#for all other OTUs, just print to file
