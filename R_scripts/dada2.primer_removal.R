####PRIMER REMOVAL IN 16S AND 18S SEQUECNES, DADA2 PIPELINE####
#author: Evan Morien
#using this ITS guide for primer removal: https://benjjneb.github.io/dada2/ITS_workflow.html
#last modified: Dec 6th, 2018

####Dependencies and Requirements####
#you must first install cutadapt ( in your system, not in R! ) to do primer removal
#the pipeline requires that you begin with demultiplexed fastq files (one fastq per sample)

####Libraries####
library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)

####Environment Setup####
theme_set(theme_bw())
setwd("/projects/my_project/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "/projects/my_project/data/raw_data/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####fastq Quality Plots####
plotQualityProfile(fnFs[1:20]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
plotQualityProfile(fnRs[1:20])

####prepare file names for filtering####
#this defines names and paths for the filtered fastq files you are about to generate
filtFs <- file.path(path, "dada2_filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "dada2_filtered", paste0(sample.names, "_R_filt.fastq"))

####identify primers####
FWD <- "CYGCGGTAATTCCAGCTC"  ## CHANGE ME to your forward primer sequence
REV <- "CRAAGAYGATYAGATACCRT"  ## CHANGE ME to your reverse primer sequence
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = FALSE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details. it is linked at the top of this guide.

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: sampleID_something_else.fastq
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#adjusting the parameters for the filterAndTrim function are crucial to the success of a dada2 run. truncation length, in particular, will be a strong determinant of the percentage of reads you retain per sample
#the parameters below represent best-practices from several different 18s experiments that we have done, tips from Ramon Masanna's lab, and parameters we have observed in the literature. the truncation length will need to be adjusted for the sequencing technology, of course. this example truncation length is OK for MiSeq data of expected (read: average) quality.
out <- filterAndTrim(cutFs, cutRs, filtFs, filtRs, maxN=0, maxEE=c(6,8), truncQ=c(2, 2), rm.phix=TRUE, compress=TRUE, verbose=TRUE, truncLen=c(240,220), multithread = TRUE, matchIDs=TRUE) #with trunQ filter
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
View(retained)

#from this point forward, use the standard guide for either 18s or 16s data processing with dada2