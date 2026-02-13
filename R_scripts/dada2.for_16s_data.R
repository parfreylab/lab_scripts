####DADA2 BEST PRACTICES FOR PROCESSING 16S SEQUENCES####
#author: Evan Morien
#using and modifying this dada2 guide as necessary: https://benjjneb.github.io/dada2/tutorial.html
#last modified: Feb 12th, 2026

####READ FIRST####
#this document is intended as a rough guide for processing 16s metabarcoding data with dada2. it is not meant to present a definitive solution for this kind of work. you will need to adjust parameters according to the dataset you are working with.
#best practices for 16s reads in our lab will involve the analysis of a merged read set
#to run this pipeline on forward reads only, just skip steps that involve the reverse reads, and skip the merging step, modifying code as necessary. should you need it, the code for an R1 only analysis section of the pipeline is present in the lab's 18s dada2 R guide.

#it's necessary to filter out low-abundance/prevalence ASVs after constructing the sequence table in dada2. dada2 developers do not provide any advice on this front, but the process is outlined below in the section after sequence table construction

####Dependencies and Requirements####
#the pipeline requires that you begin with demultiplexed, gzipped fastq files (one fastq per sample)

####IMPORTANT EXPERIMENT-SPECIFIC INFO####
#Unless otherwise noted by the sequencing facility, all datasets should be primer-trimmed. There is code below to primer trim your files. However, if you need additional help there is further documentation and explanation of primer trimming with dada2 in their recommended ITS sequence processing pipeline.

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
setwd("/projects/my_project/16s/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "/projects/my_project/16s/data/raw_data/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) #same for this pattern, for R2 files
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
#if you have hundreds of samples in your dataset, it's easiest to just look at a few (10-50) fastq files
numplot <- 49 #CHANGE ME to the number of samples you want to include in the quality plot ( can be anywhere from 2 to length(fnFs) ). for large experiments a nice 7x7 plot looks pretty good, hence a vaue of 49 in the example
a <- sample(fnFs, numplot) #randomly select a set of N samples
b <- which(fnFs %in% a) #identify the indices of those samples
plotfnFs <- fnFs[b] #use the indices to create two lists of corresponding forward and reverse files to plot
plotfnRs <- fnRs[b]
pdf("quality_plots.dada2.16S.R1s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
plotQualityProfile(plotfnFs) #this plots the quality profiles for each sample
dev.off()
pdf("quality_plots.dada2.16S.R2s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
plotQualityProfile(plotfnRs)
dev.off()

####Primer Removal####
####identify primers####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "CCGYCAATTYMTTTRAGTTT"  ## CHANGE ME to your reverse primer sequence

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

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, trimLeft = c(0,0), trimRight = c(0,0), maxN = 0, multithread = 32, compress = TRUE, matchIDs=TRUE) #right/left trimming here can be useful if there is a very low quaity base visible in every sample in the quality plots, provided it is close enough to the 3' or 5' end of the reads

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
index <- 5 #this is the index of the file we want to check for primers, within the lists "fn*s.filtN", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the sample you'd like to use for this test is used here (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))

####ONLY IF NECESSARY: fixing primer input orientation####
REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details. it is linked at the top of this guide. your primer hits table above should look like the one in the dada2 tutorial example, but with the 16S amplicon there is no read through, so only pay attention to the placement of detected primers in forward orientation (leftmost column)

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 32,# -n 2 required to remove FWD and REV from reads #-j sets no. threads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#sanity check, should report zero for all orientations and read sets
index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE)) #remember to change this so it matches ALL your file names!
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE)) #remember to change this so it matches ALL your file names!

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimLeft = c(0,0), trimRight = c(0,0), minLen = c(110,110),
                     maxN=c(0,0), maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, "retained_reads.filterAndTrim_step.length_var.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf("error_rates.dada2.R1s.length_var.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errF, nominalQ=TRUE) #assess this graph after processing. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered carefully before moving on.
dev.off()
pdf("error_rates.dada2.R2s.length_var.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errR, nominalQ=TRUE) #assess this graph after processing. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered carefully before moving on.
dev.off()

#sometimes, samples with low sequence count don't produce any filtered output, so we have to do something extra here to get rid of their entries in filtFs and filtRs
filtFs <- sort(list.files("raw_data/cutadapt/filtered/", pattern = "_R1_", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!
filtRs <- sort(list.files("raw_data/cutadapt/filtered/", pattern = "_R2_", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!

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
samples_to_keep <- track[,4] > 5 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), "samples_retained.txt", row.names=FALSE, quote=F, sep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))), "samples_removed.txt", row.names=FALSE, quote=F, sep="\n")


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
pdf("length_histogram.merged_reads.length_var.pdf", width = 10, height = 8) # define plot width and height. completely up to user.
plot(x=length.histogram[,1], y=length.histogram[,2],
     xlab = "length (bp)",
     ylab = "frequency (# observed ASVs)",
     main = "ASV Length Histogram") #view length distribution plot
dev.off()
#review after processing as sanity check. most ASVs should be w/in a few bp of the amplicon length. 
#main peak of ASV lengths centres around 355bp. i'm not familiar with the 16S amplicon, but will check this against the published methods paper

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)

otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
removed <- length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
if (removed > 0) {
  otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
} else {
  otus_filt <- otus #nothing to filter
}
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

#Start ASV report
cat("dimensions of unfiltered ASV table:", dim(otus),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# ASVs Removed:", length(intersect(a,b)),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("dimensions of filtered ASV table:", dim(otus_filt),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
cat("dimensions of ASV table after chimera removal", dim(seqtab.nosingletons.nochim),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
cat("proportion of chimeric to non chimeric reads:", sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons),file="ASV_report.txt",sep="\t",append=TRUE)

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

####save output from sequence table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16S_merged.length_var.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.16S.merged.length_var.txt", row.names=FALSE, quote=F, sep="\t")

####assign taxonomy####
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/path/to/taxonomyDBs/silva_for_dada2/v138_by_dada2team/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE) #CHANGE path here based on the reference database and amplicon you are working with

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


#### save sequences and taxonomy table ####
##### replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.length_var.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "16S_ASV_sequences.fasta") #save sequences with new names in fasta format
#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements
row.names(taxa) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtab.nosingletons.nochim) <- ASV.num
row.names(taxa) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.16S.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(taxa),taxa),"taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")


#### OPTIONAL: combine sequence and taxonomy tables into one####
#taxa will be the rows, columns will be samples, followed by each rank of taxonomy assignment, from rank1 (domain-level) to rank7/8 (species-level), followed by accession (if applicable)
#first check if the row names of the taxonomy table match the column headers of the sequence table
length(which(row.names(taxa) %in% colnames(seqtab.nosingletons.nochim)))
dim(taxa)
dim(seqtab.nosingletons.nochim)
#the number of taxa from the last three commands should match

#now ensure that the taxa in the tables are in the same order #this should be true if you haven't reordered one or the other of these matrices inadvertently
row.names(taxa) == colnames(seqtab.nosingletons.nochim) #IMPORTANT: only proceed if this evaluation is true for every element. if it isn't you'll need to re-order your data. I'd suggest sorting both matrices by their rows after transposing the sequence table.

#as long as the ordering of taxa is the same, you can combine like this (note you need to transpose the sequence table so that the taxa are in the rows)
sequence_taxonomy_table <- cbind(t(seqtab.nosingletons.nochim), taxa)
colnames(sequence_taxonomy_table) #the last elements of this list should be "Rank1", "Rank2", etc, followed by "Accession"

#now write to file
write.table(data.frame("row_names"=rownames(sequence_taxonomy_table),sequence_taxonomy_table),"sequence_taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")

####code for getting your data into a phyloseq object####
#### load either .RDS files you may have saved earlier (code not included in tutorial) or .txt formatted sequence table, taxonomy table, and metadata####
#OPTIONAL:loading from .RDS
# phyloseq_object_name <- readRDS(name_of_file_where_you_stored_phyloseq_object.RDS) #modify for whatever file you need to load. NB this works for whatever R objects you like, you can save and load any R object, preserving all formatting, class attributes, etc and load it back in on another computer this way.

####loading sequence/taxonomy tables from tab-delimited .txt ####
library(phyloseq)
library(data.table)
library(tidyverse)
seqtab.nosingletons.nochim <- fread("sequence_table.16S.w_ASV_names.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

taxatable <- fread("taxonomy_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(taxatable) <- taxatable[,1] #set row names
taxatable <- taxatable[,-1] #remove column with the row names in it
taxatable <- as.matrix(taxatable) #cast the object as a matrix

#load sample data (from tab-delimited .txt)
rawmetadata <- fread("metadata.txt", sep="\t", header=T, colClasses = c("sampleID"="character"), data.table=FALSE) #change "sampleID" to the header of the column where your sample names are stored
row.names(rawmetadata) <- rawmetadata[,1] #set row names #change the index to the column with the row names in it
rawmetadata <- rawmetadata[,-1] #remove column with the row names in it

#make a note of which samples are present/absent in the metadata vs the sequence table (samples without a corresponding entry in both of these will not be included in the final phyloseq object, so correct any mislabeled samples now (the matching must be exact)
notinmeta <- setdiff(row.names(seqtab.nosingletons.nochim), row.names(rawmetadata))
notinraw <- setdiff(row.names(rawmetadata), row.names(seqtab.nosingletons.nochim))

#combine into phyloseq object
project_data <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), #taxa_are_rows is a parameter that sets whether phyloseq will look for taxa in the rows or columns of the input matrix. make sure it is set correctly.
                         sample_data(rawmetadata), 
                         tax_table(taxatable))

#### saving phyloseq object ####
# at this point I recommend saving your full unfiltered dataset in phyloseq format as a .RDS, so that you can pick up the analysis from this point easily if you decide to change your downstream filtering criteria later on
saveRDS(project_data, "kelp16S.full_dataset.phyloseq_format.RDS")

#### OPTIONAL but GENERALLY RECOMMENDED UNLESS THERE IS A GOOD REASON NOT TO DO THIS: filtering the complete phyloseq dataset ####
# Remove samples with less than N reads (N=100 in example. adjust per experiment.)
project_data.filt <- prune_samples(sample_sums(project_data) >= 100, project_data) #i generally err on the side of caution for 18s experiments, where low read counts can reflect low diversity/low numbers of 18s organisms in the sample, while not strictly reflecting "sample failure". for 16s, a higher threshold is better. successful samples will have more than 1000 reads, but use a plot of the sorted sample sums to find what looks like an appropriate data-informed cutoff, if possible.

# OPTIONAL: Remove OTUs with less than N total reads. (N = 50 for example. adjust per experiment)
#project_data.filt <- prune_taxa(taxa_sums(project_data) >= 50, project_data) #again, discretion is necessary here. in the dada2 pipeline, this is not necessary becasue a similar filter is applied already in the pipeline. for MED/QIIME, user's discretion.

# Set aside unassigned taxa #with dada2, there may not be any unassigned taxa as dada2's RDP classifier usually ends up classifying everything. You may want to adjust this command to set aside ASVs that have no assignment beyond Rank1 instead of no assignment at Rank1.
project_data.unassigned <- project_data.filt %>%
  subset_taxa(Kingdom == "Unassigned")
#remove unassigned
project_data.filt <- project_data.filt %>%
  subset_taxa(Kingdom != "Unassigned")
# Remove non-prokaryotes
project_data.filt <- project_data.filt %>%
  subset_taxa(Kingdom != "Eukaryota")

# zero out counts of 2 or less from OTU table #not zeroing out these low counts this at this time
# this is a de-noising procedure applicable to both dada2 and MED/QIIME produced datasets which minimizes cross-well contamination
# otu <- as.data.frame(otu_table(project_data.filt))
# otu_table(project_data.filt)[otu <= 2] <- 0

# additional filtering steps, if not already done before loading into R
# 16S ONLY: Remove mitochondrial and chloroplast OTUs
# IMPORTANT: make sure that the filter terms will work with your taxonomy strings, and ranks. mitochondria and chloroplast may be at different ranks in whatever database you are using. check your taxonomy table manually!
project_data <- project_data.filt %>%
  subset_taxa(Family != "Mitochondria") %>%
  subset_taxa(Order != "Chloroplast")

# 18S (and optional for 16S): Remove unwanted clades
# you can chain as many of these subset_taxa calls as you like into this single command using the R pipe (%>%)
# project_data <- project_data %>%
#   subset_taxa(RankX != "UNWANTED_CLADE") %>%
#   subset_taxa(RankY != "UNWANTED_CLADE")

# modify Rank labels in taxa table (check the colnames of the tax_table(project_data) object to see if you want to change them)
colnames(tax_table(project_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#### accessing data within phyloseq objects ####
#this can be a bit tricky, but you just have to treat the component parts of the S4 object as data frames
#it is sometimes necessary to strip the class information out with "unclass()" before viewing the data tables with the View() function
#for viewing
View(as.data.frame(unclass(sample_data(project_data))))
View(as.data.frame(unclass(otu_table(project_data))))
View(as.data.frame(unclass(tax_table(project_data))))

#OPTIONA but usefulL: creating separate data frames for use with other programs
my_sample_data <- as.data.frame(unclass(sample_data(project_data))) #etc

#### save filtered phyloseq as RDS ####
saveRDS(project_data, "kelp16S.filtered_dataset.phyloseq_format.RDS") #this is the ready-to-analyse dataset