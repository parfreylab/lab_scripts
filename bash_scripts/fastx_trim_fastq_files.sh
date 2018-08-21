#!/bin/sh

#this script takes a list of files (specify with wildcard, for example. can include full or relative filepath.) and trims them using fastx toolkit to 250 base pairs, throwing out reads that are less than 250bp long
#will work with either fasta or fastq files
#DEPENDS: fastx toolkit

for FILE in "$@"
do
    fastx_trimmer -l 250 -i $FILE | fastx_clipper -v -a NNNNNNNNNNNNN -l 250 -o $FILE.trimmed_clipped.fastq >> fastx_trim_clip.log
done
