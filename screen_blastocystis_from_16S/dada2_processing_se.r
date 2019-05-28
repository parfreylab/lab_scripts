#####################################
#Blastocystis screen for 16S V4 data
#####################################
#Dada2 processing pipeline for single end reads

######CHANGE THE VALUES BELOW#######
PATH="/Users/parfreylab/Desktop/lab_member_files/mann" #change this to your working directory
RAW="/Users/parfreylab/Desktop/lab_member_files/mann/raw" #change this to the directory where your raw data (fastq) is stored
FWPRI="GTGYCAGCMGCCGCGGTAA" #change this to your forward primer sequence
RVPRI="GGACTACNVGGGTWTCTAAT" #change this to your reverse primer sequence
CUTAD="/anaconda2/bin/cutadapt" #change this to cutadapt installation path
VSEARCH="/usr/local/bin/vsearch" #change this to vsearch installation path
BACTAX="/Users/parfreylab/Desktop/lab_member_files/mann/ezbiocloud_id_taxonomy.txt" #change this to your bacterial reference db taxonomy file
BACREF="/Users/parfreylab/Desktop/lab_member_files/mann/ezbiocloud_qiime_full.fasta" #change this to your bacterial reference db fasta file
EUKTAX="/Users/parfreylab/Desktop/lab_member_files/mann/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt" #change this to your eukaryotic reference db taxonomy file
EUKREF="/Users/parfreylab/Desktop/lab_member_files/mann /99_SILVA_128_euks_rep_set.fasta" #change this to your eukaryotic reference db fasta file
#####################################

####Libraries####
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd(PATH)

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- RAW
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
pdf("f_qualityplot.pdf")
plotQualityProfile(fnFs[1:3]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
dev.off()

####Primer Removal####
####identify primers####
FWD <- FWPRI
REV <- RVPRI
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

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), #add the index of the sample you'd like to use for this test (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))

#### primer removal ####
cutadapt <- CUTAD
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

# fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, truncLen= 0, minLen = 150,
                     maxN= 0, maxEE= 8, truncQ= 2, rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
#View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)

pdf("error_rates_plot.pdf")
plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaFs[[1]]

getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(dadaFs, getN))
samples_to_keep <- track[,2] > 50

####construct sequence table####
seqtab <- makeSequenceTable(derepFs[samples_to_keep])
dim(seqtab)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), rowSums(seqtab.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nochimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.txt", row.names=FALSE, quote=F, sep="\t")

##if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
# seqtab.nochim <- fread("sequence_table.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
# row.names(seqtab.nochim) <- seqtab.nochim[,1] #set row names
# seqtab.nochim <- seqtab.nochim[,-1] #remove column with the row names in it
# seqtab.nochim <- as.matrix(seqtab.nochim) #cast the object as a matrix

####assign taxonomy####
uniquesToFasta(seqtab.nochim, "rep_set.fa")

##assign taxonomy using vsearch
##first with bacterial database
##if using usearch instead of vsearch comment out the line below
system2("alias", args = paste("uclust=", VSEARCH, sep=""))
system2("assign_taxonomy.py", args = c("-i rep_set.fa", paste("-t", BACTAX, sep=""), paste("-r", BACREF, sep=""), "-o assigntax", "-m uclust"))

##now pull unassigned and those without strong bac assignment (no phylum assignment), reassign these with eukaryotic database
system2("grep", args= "-E 'Unassigned|Bacteria\t' assigntax/rep_set_tax_assignments.txt | awk '{print $1}' > unassigned.ids")
system2("filter_fasta.py", args= "-f rep_set.fa -o unassigned.fa -s unassigned.ids")
system2("assign_taxonomy.py", args = c("-i unassigned.fa", paste("-t", EUKTAX, sep=""), paste("-r", EUKREF, sep=""), "-o assigntax", "-m uclust"))

##any blastocystis?
system2("echo", "Number of Blastocystis reads:")
system2("grep", args = "Blastocystis assigntax/unassigned_tax_assignments.txt")