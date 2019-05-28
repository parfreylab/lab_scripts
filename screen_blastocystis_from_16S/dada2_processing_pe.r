#####################################
#Blastocystis screen for 16S V4 data
#####################################
#Dada2 pipeline for paired end reads

######CHANGE THE VALUES BELOW#######
PATH="/Users/parfreylab/Desktop/lab_member_files/mann/afribiota" #change this to your working directory
RAW="/Users/parfreylab/Desktop/lab_member_files/mann/afribiota/raw" #change this to the directory where your raw data (fastq) is stored
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
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE)) #same for this pattern, for R2 files
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
pdf("f_qualityplot.pdf")
plotQualityProfile(fnFs[1:20]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
dev.off()
pdf("r_qualityplot.pdf")
plotQualityProfile(fnRs[1:20])
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
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), #add the index of the sample you'd like to use for this test (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#### primer removal ####
cutadapt <- CUTAD
cutadapt
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
cutFs <- sort(list.files(path.cut, pattern = "_1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(8,10), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
#View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf("error_rates_plot.pdf")
plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below

getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("length_histogram.pdf")
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
dev.off()

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
sum(seqtab.nochim)/sum(seqtab) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")

##if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
#seqtab.nochim <- fread("sequence_table.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
#row.names(seqtab.nochim) <- seqtab.nochim[,1] #set row names
#seqtab.nochim <- seqtab.nochim[,-1] #remove column with the row names in it
#seqtab.nochim <- as.matrix(seqtab.nochim) #cast the object as a matrix

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