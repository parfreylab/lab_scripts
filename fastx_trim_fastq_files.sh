#!/bin/sh

#this script takes a list of files (specify with wildcard, for example. can include full or relative filepath.) and trims them using fastx toolkit to 250 base pairs
#will work with either fasta or fastq files
#DEPENDS: fastx toolkit

for FILE in "$@"
do
    fastx_trimmer -l 250 -i $FILE -o $FILE.trimmed.fastq
done
