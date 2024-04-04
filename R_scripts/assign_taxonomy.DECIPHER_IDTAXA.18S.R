#to be done after processing amplicon data with dada2
#last modified: Apr 3rd, 2024
#author: Evan Morien

####libraries####
library(tidyverse)
library(data.table)
library(DECIPHER)
library(Biostrings)
library(ensembleTax) #note, please install this from my fork until the devs fix their issues with PR2's updated taxonomy. use the devtools install method here, but replace "dcat4" with "morien" and set the vignettes parameter to FALSE: https://github.com/morien/ensembleTax

####Inputs####
#sequence table from dada2 (see dada2 processing scripts for any amplicon). ASVs not labeled (i.e. raw sequences as ASV names)
seqtab.nosingletons.nochim <- fread("sequence_table.18S.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix
mode(seqtab.nosingletons.nochim) <- "numeric"

## creating a DNAStringSet object from the ASVs
dna_18S <- DNAStringSet(getSequences(seqtab.nosingletons.nochim))

## loading PR2 DECIPHER train set #this is provided by PR2 developers, see PR2's github page for updated versions
trainingset <- readRDS("/data/taxonomyDBs/DECIPHER_IDTAXA/pr2_version_5.0.0_SSU.decipher.trained.rds")

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_18S <- IdTaxa(dna_18S, trainingset, type="extended", strand = "both", processors = 36) 

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_18S <- colnames(seqtab.nosingletons.nochim)

rubric <- DNAStringSet(getSequences(seqtab.nosingletons.nochim))
# this creates names (ASV1, ASV2, ..., ASVX) for each ASV
snam <- vector(mode = "character", length = length(rubric))
pad <- nchar(length(snam)) #length of zeroes to pad ASV names. we do this so the alphanumeric ASV names are sorted correctly by the idtax2df function
for (i in 1:length(rubric)) {
    snam[i] <- str_pad(as.character(i), width=pad, side="left", pad="0")
    snam[i] <- paste0("ASV", snam[i])
}
names(rubric) <- snam

tax_tab_18S <- idtax2df(tax_info_18S, 
                             db = "pr2", 
                             ranks = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"),
                             boot = 50,
                             rubric = rubric,
                             return.conf = FALSE)


colnames(tax_tab_18S) <- c("ASV_ID", "sequence", "domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")
row.names(tax_tab_18S) <- NULL
#we added leading zeroes earlier to allow correct sorting order by the idtax2df function, now remove them so the ASV IDs match the other IDs in the rest of our pipeline
pattern <- "ASV0*"
replacement <- "ASV"
tax_tab_18S$ASV_ID <- gsub(pattern, replacement, tax_tab_18S$ASV_ID)
write.table(tax_tab_18S, "taxonomy_table.18S.DECIPHER.PR2v5.0.txt", sep = "\t", quote = F, row.names = FALSE)

#### slightly different procedure for SILVA v 138 ####
## loading SILVA v138 DECIPHER train set #this is created and provided by DECIPHER developers, see downloads section of DECIPHER website
load("/data/taxonomyDBs/DECIPHER_IDTAXA/SILVA_SSU_r138_2019.RData")

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_18S <- IdTaxa(dna_18S, trainingSet, strand = "both", processors = 36) 

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV1, ASV2... ASVN)
asv_seqs_18S <- colnames(seqtab.nosingletons.nochim)

asv_headers_18S <- vector(dim(seqtab.nosingletons.nochim)[2], mode = "character")
for (i in 1:dim(seqtab.nosingletons.nochim)[2]) {
  asv_headers_18S[i] <- paste(">ASV", i, sep = "")
}

# creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab_18S <- t(sapply(tax_info_18S, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

colnames(tax_tab_18S) <- ranks
row.names(tax_tab_18S) <- NULL
tax_tab_18S <- data.frame("ASV_ID" = sub(">", "", asv_headers_18S), tax_tab_18S, check.names = FALSE)
write.table(tax_tab_18S, "taxonomy_table.18S.DECIPHER.SILVAv138.txt", sep = "\t", quote = F, row.names = FALSE)


#### formatting final output files ####
#### replace the long ASV names (the actual sequences) with human-readable names ####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.18S.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "18S_ASV_sequences.fasta") #save sequences with new names in fasta format

#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements

#assign new ASV names
colnames(seqtab.nosingletons.nochim) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.18S.merged.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")
