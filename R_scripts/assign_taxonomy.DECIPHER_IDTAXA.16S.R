#to be done after processing amplicon data with dada2
#last modified: Feb 17th, 2026
#author: Evan Morien

#### taxonomy assignment with DECIPHER/IDTAXA ####
#### code for IDTAXA assignment with 16S amplicon data ####
####libraries####
library(tidyverse)
library(data.table)
library(DECIPHER)
library(Biostrings)
library(ensembleTax) #note, please install this from my fork until the devs fix their issues with PR2's updated taxonomy. use the devtools install method here, but replace "dcat4" with "morien" and set the vignettes parameter to FALSE: https://github.com/morien/ensembleTax

#code to install our fork of ensembleTax
#library(devtools)
#devtools::install_github("morien/ensembleTax", build_manual = TRUE, build_vignettes = FALSE, force=TRUE)
#packageVersion("ensembleTax")

####Inputs####
#sequence table from dada2 (see dada2 processing scripts for any amplicon). ASVs not labeled (i.e. raw sequences as ASV names)
seqtab <- fread("sequence_table.16S.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab) <- seqtab[,1] #set row names
seqtab <- seqtab[,-1] #remove column with row names in it
seqtab <- as.matrix(seqtab) #cast the object as a matrix
mode(seqtab) <- "numeric"

## creating a DNAStringSet object from the ASVs
dna_16S <- DNAStringSet(getSequences(seqtab))

#there are three choices for reference database (as of Nov. 2025): PR2 (IDTAXA provided), SILVA (Hakai custom db), or MetaZooGene (Hakai custom DB)
#there is a section of this example code for each, read each carefully if you would like to use it

#### loading PR2 DECIPHER train set #this is provided by PR2 developers, see PR2's github page for updated versions ####
load("~/projects/taxonomyDBs/DECIPHER_IDTAXA/SILVA_SSU_r138.2.RData") #note the name of the object you load in. in this case it is "trainingSet" (see next command) but in other cases it may be different

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_16S <- IdTaxa(dna_16S, trainingSet, type="extended", strand = "both", processors = 8, threshold = 50) 

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_16S <- colnames(seqtab)

# this creates names (ASV1, ASV2, ..., ASVX) for each ASV
snam <- vector(mode = "character", length = length(dna_16S))
pad <- nchar(length(snam)) #length of zeroes to pad ASV names. we do this so the alphanumeric ASV names are sorted correctly by the idtax2df function
for (i in 1:length(dna_16S)) {
  snam[i] <- str_pad(as.character(i), width=pad, side="left", pad="0")
  snam[i] <- paste0("ASV", snam[i])
}
names(dna_16S) <- snam

tax_tab_16S <- idtax2df(tax_info_16S, 
                        db = "silva", 
                        boot = 50,
                        rubric = dna_16S,
                        return.conf = FALSE)

colnames(tax_tab_16S) <- c("ASV_ID", "sequence", "domain", "phylum", "class", "order", "family", "genus")
row.names(tax_tab_16S) <- NULL
#we added leading zeroes earlier to allow correct sorting order by the idtax2df function, now remove them so the ASV IDs match the other IDs in the rest of our pipeline
pattern <- "ASV0*"
replacement <- "ASV"
tax_tab_16S$ASV_ID <- gsub(pattern, replacement, tax_tab_16S$ASV_ID)
write.table(tax_tab_16S, "taxonomy_table.16S.DECIPHER.SILVA_138.2.txt", sep = "\t", quote = F, row.names = FALSE)

# also create corresponding confidence table. this is useful for interpretation, but assignments can be filtered above by passing a threshold to the IdTaxa() function
ranks <- c("domain", "phylum", "class", "order", "family", "genus")
confidences <- t(sapply(tax_info_16S, function (x) {
  m<-match(ranks, x$rank)
  confidence<-x$confidence[m]}
))
colnames(confidences) <- ranks
row.names(confidences) <- NULL
confidences <- data.frame("ASV_ID" = sub(">", "", names(dna_16S)), confidences, check.names = FALSE)
pattern <- "ASV0*"
replacement <- "ASV"
confidences$ASV_ID <- gsub(pattern, replacement, confidences$ASV_ID)
write.table(confidences, "taxonomy_table.16S.DECIPHER.SILVA_138.2.conf.txt", sep = "\t", quote = F, row.names = FALSE)
