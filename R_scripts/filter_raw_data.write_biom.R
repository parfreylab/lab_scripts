#R code (can be generalized, but is not) for filtering and creating a new OTU table from an existing one based on a list of samples
#author: evan morien
#last revision: December 1st, 2017

#load libraries
library(biomformat)
library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the phyloseq package for microbial community analysis
library(phyloseq)

#### Filtering raw data for Rhea ####
#start with the unfiltered raw data, then filter based on her list of samples
setwd("/projects/kelp_seagrass_18s/")

#### Load Required/Optional Files ####
# REQUIRED: Read the raw data in .biom format, and store as a phylo_seq object called "rawdata"
rawdata <- import_biom(file.path("/projects/kelp_seagrass_18s/closed_ref_OTU_picking_similarity100/", "OTU_Table.inherited_accessions.wtaxa.biom"), # file.path() is used for cross-platform compatibility
                       parallel = TRUE,
                       trim_ws = TRUE)

rawmetadata <- read_delim(file = file.path("/projects/kelp_seagrass_18s/metadata/", "seagrass.mapping.all_18s.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE # remove leading and trailing spaces from character string entries
)


project_data <- merge_phyloseq(rawdata, rawmetadata)


#load in rhea's list of samples
rhea_samples <- read.delim(file="metadata/rhea_samples_18s.txt", sep="\t", header=F)
length(intersect(rhea_samples[,1], colnames(otu_table(rawdata))))
length(intersect(rhea_samples[,1], rawmetadata$`#SampleID`))

#436 are in the raw data, but now I'll check to see which ones aren't there
#notinmeta <- setdiff(sample_names(rawdata), rawmetadata$`#SampleID`)
notinraw <- setdiff(rawmetadata$`#SampleID`, sample_names(rawdata))
notinraw.rheasamplelist <- setdiff(rhea_samples$V1, sample_names(rawdata)) #same list of samples not in raw data #these are failed samples

#some water samples are not included in rhea's list, just combine unique entries from her list and the raw metadata, and re-write the list
#new_sample_list <- unique(c(as.character(rawmetadata$`#SampleID`), as.character(rhea_samples$V1)))
#write.table(new_sample_list, "metadata/rhea_samples_18s.txt", sep="\t", quote=F, row.names=F, col.names=F)

#filter phyloseq object by sample names
seagrass_set <- prune_samples(as.character(rhea_samples[,1]), rawdata)

#removing zostera OTUs
seagrass_set <- seagrass_set %>%
  subset_taxa(Rank7 != "__Zostera")
#removing mitochondrial OTUs
seagrass_set <- seagrass_set %>%
  subset_taxa(Rank5 != "__Mitochondria")

#filter OTU table based on counts (no samples with <1000 reads, no OTUs with <250)
seagrass_set <- prune_samples(sample_sums(seagrass_set) >= 1000, seagrass_set)
seagrass_set <- prune_taxa(taxa_sums(seagrass_set) >= 250, seagrass_set)
colnames(tax_table(seagrass_set)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8")

#### saving filtered OTU table ####
#NOTE: for the time being, do not use biomformat to write .biom files, as the files written with this package can't be read in by phyloseq. I will update this file if that changes.
#convert phyloseq object to data frames containing taxonomy and OTU tables
otu_table <- as.data.frame(otu_table(seagrass_set))
taxa_table <- as.data.frame(tax_table(seagrass_set))

#join taxonomies into QIIME-formatted single column
taxa_table <- taxa_table %>% unite("Taxonomy", Rank1, Rank2, Rank3, Rank4, Rank5, Rank6, Rank7, Rank8, sep="; __")

#cbind data frames #they're in the same order #you should double check each time to ensure you don't mis-assign taxonomies at this step
otu_table <- cbind(otu_table, taxa_table$Taxonomy)
names(otu_table)[names(otu_table) == 'taxa_table$Taxonomy'] <- 'Taxonomy' #adjust column name
#write table to disk
write.table(otu_table, "closed_ref_OTU_picking_similarity100/OTU_Table.inherited_accessions.wtaxa.full_seagrass_dataset.2015-2016.txt", col.names=T, row.names=T, quote=F, sep="\t")
#you'll need to add the string "OTU_ID" followed by a single tab to the beginning of the header row after writing this file. do with a simple text editor like nano, vi, gedit, or sublime text. do not use rich text editors or programs like excel that use formatting.

#convert written .txt formatted table to .biom if desired with a command like this:
#biom convert -i closed_ref_OTU_picking_similarity100/OTU_Table.inherited_accessions.wtaxa.full_seagrass_dataset.2015-2016.txt -o closed_ref_OTU_picking_similarity100/OTU_Table.inherited_accessions.wtaxa.full_seagrass_dataset.2015-2016.biom --process-obs-metadata taxonomy --to-json --table-type="OTU table" --header-key taxonomy

