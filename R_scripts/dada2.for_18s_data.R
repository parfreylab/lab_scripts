####BEST PRACTICES FOR DADA2 READ PROCESSING WITH 18S DATA####
#author: Evan Morien
#using and modifying this dada2 guide as necessary: https://benjjneb.github.io/dada2/tutorial.html
#last modified: Nov 28th, 2018

####READ FIRST####
#this document is intended as a rough guide for processing 18s metabarcoding data with dada2. it is not meant to present a definitive solution for this kind of work.
#eukaryotic species vary greatly in the length of different regions of the 18s gene, and doing any kind of trimming or truncation will exclude clades that have longer sequenced regions (for example, the V4 region is quite variable in length)
#but it's also sometimes necessary to truncate your reads with dada2 to retain enough reads for analysis, so you will need to balance these two things against each other when using this pipeline.
#best practices in our lab will involve the analysis of both merged and forward read sets separately, to identify clades excluded by any truncation/trimming done prior to merging.
#to run this pipeline on forward reads only, just skip steps that involve the reverse reads, and skip the merging step, modifying code as necessary.

#it's also necessary to filter out low-abundance ASVs after constructing the sequence table in dada2. dada2 developers do not provide any advice on this front, but the process is outlined below in the section after sequence table construction

#our standard procedure with 18s reads processed with dada2 will be to do an analysis with both the merged reads and the R1s only, to make sure we aren't missing any important clades that have V4 regions too long to merge
#the R1-only analysis procedure follows the merged read procedure in this guide

####Dependencies and Requirements####
#the pipeline requires that you begin with demultiplexed fastq files (one fastq per sample)

####Libraries####
library(dada2)
library(phyloseq)
library(ShortRead)
library(Biostrings)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)

####Environment Setup####
theme_set(theme_bw())
setwd("/projects/my_project/18s/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "/projects/my_project/18s/data/raw_data/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####fastq Quality Plots####
plotQualityProfile(fnFs[1:20]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
plotQualityProfile(fnRs[1:20])

####prepare file names for filtering####
#this defines names and paths for the filtered fastq files you are about to generate
filtFs <- file.path(path, "dada2_filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "dada2_filtered", paste0(sample.names, "_R_filt.fastq"))

####filter and trim reads####
#adjusting the parameters for the filterAndTrim function are crucial to the success of a dada2 run. truncation length, in particular, will be a strong determinant of the percentage of reads you retain per sample
#the parameters below represent best-practices from several different 18s experiments that we have done, tips from Ramon Masanna's lab, and parameters have observed in published literature. the truncation length will need to be adjusted for the sequencing technology, of course. this example truncation length is OK for MiSeq data.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(6,8), truncQ=c(2, 2), rm.phix=TRUE, compress=TRUE, verbose=TRUE, truncLen=c(240,220), multithread = TRUE, matchIDs=TRUE) #with trunQ filter
retained_18s <- as.data.frame(out)
retained_18s$percentage_retained <- retained_18s$reads.out/retained_18s$reads.in*100
View(retained_18s)

####learn error rates####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.

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

####merge paired reads####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) #a "subscript out of bounds" error here may indicate that you aren't merging any reads in one or more samples. you can remove samples with low counts from the workflow before the filterAndTrim step (a few steps back), or you can filter samples using the information from the dereplication and sample-inference steps (section just above)
#OPTIONAL: modify command if removing low-sequence samples
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 1)) #how many are there with N reads or fewer? (N=1 in example)

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names
write.table(track, "read_retention.merged_18s.txt", quote=F, row.names=T, col.names=T, sep="\t")

####assign taxonomy####
#note, this takes ages if you have a large dataset. strongly recommend doing on a multi-core machine (zoology cluster, or entamoeba in the lab). another option: saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function (their implementation of RDP).
#for 18s mammalian gut experiments
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/PATH/TO/taxonomy_databases/silva_for_dada2/v128_for_parfreylab/18s/silva_128.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)
#for other 18s experiments
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/PATH/TO/taxonomy_databases/silva_for_dada2/v132_for_parfreylab/18s/silva_132.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
#set column ranks for the taxa table. should be Rank1 through RankN, depending on how many ranks you have. last column is the accession for our parfreylab custom-formatted databases, and species for other publicly available databases.
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Accession") #for parfreylab version of SIVLA 128
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession") #for parfreylab version of SIVLA 132
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7") #for dada2 dev's version of SIVLA 132

####saving data####
#please note you can also save these as plain text files. Those commands will replace these when I have time to make sure I have generalizable code for it.
write.table(taxa, "taxonomy_table.18s.txt", sep="\t", quote=F, row.names=T, col.names=T)
write.table(seqtab.nosingletons.nochim, "sequence_table.18s.txt", sep="\t", quote=F, row.names=T, col.names=T)

#### hand off to PhyloSeq ####
#load sample data
rawmetadata <- read_delim(file = file.path("/projects/my_project/18s/", "mapping_file.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

#set row names as the sample IDs for the metadata
rownames(rawmetadata) <- rawmetadata$SampleID

#record samples absent in either metadata or OTU table
notinmeta <- setdiff(row.names(seqtab.nosingletons.nochim), rawmetadata$SampleID)
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.nosingletons.nochim))

#create phyloseq object from "seqtab.nosingletons.nochim", "rawmetadata", and "taxa"
ps.dada2_join <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(taxa))

#set column names for taxonomy table
colnames(tax_table(ps.dada2_join)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession")

# Remove samples with less than N reads (N=100 in example. adjust per experiment.)
ps.dada2_join <- prune_samples(sample_sums(ps.dada2_join) >= 100, ps.dada2_join)

# Remove OTUs with less than N total reads. (N = 250 for example) 
ps.dada2_join <- prune_taxa(taxa_sums(ps.dada2_join) >= 50, ps.dada2_join)

# Set aside unassigned taxa #with dada2, there may not be any unassigned taxa as dada2's RDP classifier usually ends up classifying everything. You may want to adjust this command to set aside ASVs that have no assignment beyond Rank1 instead of no assignment at Rank1.
ps.dada2_join.unassigned <- ps.dada2_join %>%
  subset_taxa(Rank1 == "Unassigned")
# Remove unassigned taxa
ps.dada2_join <- ps.dada2_join %>%
  subset_taxa(Rank1 != "Unassigned")

# Remove counts of 1 from OTU table
otu <- as.data.frame(otu_table(ps.dada2_join))
otu_table(ps.dada2_join)[otu <= 1] <- 0

#after denoising you can also remove ASVs that are in groups you wish to exclude (i.e. mammalia, embryophyta, etc.)
#to do this, just determine which rank the clade is captured by, and filter like so:
# Remove mitochondrial and chloroplast OTUs
ps.dada2_join <- ps.dada2_join %>%
  subset_taxa(Rank5 != "Mammalia") %>% #you can chain as many of these subset_taxa calls as you like into this single command using the R pipe (%>%)
  subset_taxa(Rank5 != "Embryophyta")

#with your filtered phyloseq object, you are now ready to move on to whatever statistical/diversity analyses you're planning to do. please see other R script guides in the lab github.



####R1-ONLY ANALYSIS GUIDE####
####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "/projects/my_project/18s/data/raw_data/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####fastq Quality Plots####
plotQualityProfile(fnFs[1:20]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.

####prepare file names for filtering####
#this defines names and paths for the filtered fastq files you are about to generate
filtFs <- file.path(path, "dada2_filtered", paste0(sample.names, "_F_filt.fastq"))

####filter and trim reads####
#adjusting the parameters for the filterAndTrim function are crucial to the success of a dada2 run. truncation length, in particular, will be a strong determinant of the percentage of reads you retain per sample
#the parameters below represent best-practices from several different 18s experiments that we have done, tips from Ramon Masanna's lab, and parameters have observed in published literature. the truncation length will need to be adjusted for the sequencing technology, of course. this example truncation length is OK for MiSeq data.
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=c(6), truncQ=c(2), rm.phix=TRUE, compress=TRUE, verbose=TRUE, truncLen=c(240), multithread = TRUE, matchIDs=TRUE) #with trunQ filter
retained_R1s <- as.data.frame(out)
retained_R1s$percentage_retained <- retained_R1s$reads.out/retained_R1s$reads.in*100
View(retained_R1s)

####learn error rates####
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs[[1]]

####construct sequence table####
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 1)) #how many are there with N reads or fewer? (N=1 in example)

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,4]/track[,3]*100, 100-track[,5]/track[,4]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names
write.table(track, "read_retention.merged_18s_R1s.txt", quote=F, row.names=T, col.names=T, sep="\t")

####assign taxonomy####
#note, this takes ages if you have a large dataset. strongly recommend doing on a multi-core machine (zoology cluster, or entamoeba in the lab). another option: saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function (their implementation of RDP).
#for 18s mammalian gut experiments
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/PATH/TO/taxonomy_databases/silva_for_dada2/v128_for_parfreylab/18s/silva_128.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)
#for other 18s experiments
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "/PATH/TO/taxonomy_databases/silva_for_dada2/v132_for_parfreylab/18s/silva_132.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
#set column ranks for the taxa table. should be Rank1 through RankN, depending on how many ranks you have. last column is the accession for our parfreylab custom-formatted databases, and species for other publicly available databases.
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Accession") #for parfreylab version of SIVLA 128
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession") #for parfreylab version of SIVLA 132
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7") #for dada2 dev's version of SIVLA 132

####saving data####
#please note you can also save these as plain text files. Those commands will replace these when I have time to make sure I have generalizable code for it.
write.table(taxa, "taxonomy_table.18s_R1s.txt", sep="\t", quote=F, row.names=T, col.names=T)
write.table(seqtab.nosingletons.nochim, "sequence_table.18s_R1s.txt", sep="\t", quote=F, row.names=T, col.names=T)

#the handoff to phyloseq remains the same for both R1 and merged procedures, see above section.