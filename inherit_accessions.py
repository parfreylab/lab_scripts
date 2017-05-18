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



