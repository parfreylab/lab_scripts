#basic procedure for loading in a .biom file, processing data with DeSeq2
#NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT IF YOU WISH TO MODIFY IT
#       2. CHANGE ALL GENERALIZED PARAMETERS/VARIABLES (EX: "FACTOR_1") TO MATCH YOUR DATA

#### load libraries ####
#you must install these first if you want to load the data in using phyloseq and process with deseq
library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the data.table package for better metadata manipulation
library(data.table)
# load the deseq2 library
library(DESeq2)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)

#### environment settings ####
#set working directory
setwd("/path/to/working/directory/")

#### load in .biom file, metadata, and (optional) phylogenetic tree ####
# REQUIRED: Read the raw data in .biom format, and store as a phylo_seq object called "rawdata"
rawdata <- import_biom(file.path("/file/path/", "OTU_Table.biom"), # file.path() is used for cross-platform compatibility
                       parallel = TRUE,
                       trim_ws = TRUE) # use multiple processor cores for shorter read times

# REQUIRED: Read the raw metadata in .txt format, and store as a tibble (an improved dataframe) called "rawmetadata"
rawmetadata <- read_delim(file = file.path("/file/path/", "metadata.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the metadata file must be tab delimited and in .txt format
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

# IMPORTANT: ROW NAMES OF THE RAW METADATA FILE MUST MATCH SAMPLE NAMES FROM THE OTU TABLE. IF NOT, THE PHYLOSEQ OBJECT WILL NOT BE MERGED PROPERLY AND NON-MATCHING IDS WILL BE THROWN OUT
# RECOMMENDED: preserve a list of samples that are unique to the OTU table and metadata file to make sure there aren't unexpected differences
# NOTE: change MY_SAMPLE_IDs to match the column where sample IDs are stored in the metadata file. in QIIME compatible metadata files, these will be under `#SampleID`
notinmeta <- setdiff(sample_names(rawdata), rawmetadata$MY_SAMPLE_IDs)
# if there are any samples in "notinmeta", remove these samples from the raw data
allsamples <- sample_names(rawdata)
allsamples <- allsamples[!(allsamples %in% notinmeta)]
biomdata <- prune_samples(allsamples, rawdata)
# Find all samples in raw metadata that are not in raw data
notinraw <- setdiff(rawmetadata$MY_SAMPLE_IDs, sample_names(rawdata))
# Remove these samples from the metadata
metadata <- filter(rawmetadata, !(MY_SAMPLE_IDs %in% notinraw))
# Create vector consisting of ALL removed samples; useful for record-keeping
symmdiff <- c(notinmeta, notinraw)
# Format metadata tibble into the phyloseq metadata format
metadata <- sample_data(metadata)
# Add rownames to the phyloseq metadata which correspond to the OTU IDs
rownames(metadata) <- metadata$MY_SAMPLE_IDs

# OPTIONAL: Read in a phylogenetic tree in .tre format, and store as a tree object called "rawtreedata"
# IMPORTANT: TREE TIPS MUST MATCH OTU IDs FROM TAXA TABLE. THIS MEANS IF YOU ARE USING AN EPA PLACEMENT TREE, YOU MUST REMOVE "QUERY___" from each tree tip label BEFORE IMPORTING if you wish to use the tree in a phyloseq object
# file.path() is used for cross-platform compatibility
rawtreedata <- read_tree(file.path("/file/path", "phylo_tree.tre"))

#IMPORTANT NOTE: at this point you should make sure your sample IDs in the data, metadata, and tree data objects match

#### create phyloseq object with completed metadata, otu table, and tree ####
project_data <- merge_phyloseq(biomdata, metadata, rawtreedata)
#filtering steps, if not already done before loading into R
#filter out samples with less than 1000 reads (arbitrary threshold and generally the minimum, choose your own, use sample_counts() to look at distribution of read counts per sample)
project_data <- prune_samples(sample_sums(project_data) >= 1000, project_data) 
# Remove OTUs with less than N total reads. (N = 250 in example) 
project_data <- prune_taxa(taxa_sums(project_data) >= 250, project_data)
# 16S ONLY: Remove mitochondrial and chloroplast OTUs #IMPORTANT: make sure that the filter terms will work with your taxonomy strings, and ranks
project_data <- project_data %>%
  subset_taxa(Rank5 != "__Mitochondria") %>% 
  subset_taxa(Rank3 != "__Chloroplast")

# 18S (and optional for 16S): Remove unwanted clades 
project_data <- project_data %>%
  subset_taxa(Rank5 != "UNWANTED_HOST_FAMILY") %>% 
  subset_taxa(Rank7 != "UNWANTED_CLADE")

# Preserve unassigned taxa
project_data.unassigned <- project_data %>%
  subset_taxa(Rank1 == "Unassigned") #works as long as your OTUs without a taxonomy assignment are labeled as "Unassigned". adjust accordingly.
# Remove unassigned taxa
project_data <- project_data %>%
  subset_taxa(Rank1 != "Unassigned")
# Remove counts that may represent noise, use a threshold (we are using a threshod of 2 reads for most datasets. be sure to choose the correct value for your own data.)
otu <- as.data.frame(otu_table(project_data)) #get OTU table
otu_table(project_data)[otu <= 2] <- 0 #for entries where the raw abundance of an OTU in a sample is less than N (N=2 in the example), set the raw read count to 0

# OPTIONAL: modify Rank labels in taxa table (check the colnames of the tax_table(project_data) object to see if you want to change them)
colnames(tax_table(project_data)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")

####DESeq2 analysis: taxa correlated with infection and death controlling for timepoint ####
#OPTIONAL/IF NEEDED:assign variables as either factors or numeric
#this can be necesssary for certain types of data, and helpful if you want to control the order of factors in a deseq contrast below
sample_data(project_data)$variable <- as.numeric(sample_data(project_data)$variable)
sample_data(project_data)$variable3 <- factor(sample_data(project_data)$variable3, levels=c("level1", "level2"))
##deseq##
#parameters up front
alpha <- 0.01 #your significance threshold for MULTIPLE TEST CORRECTED pvals
##OPTIONAL: SUBSET DATA BEFORE TEST
#if i want to test only a subset of my data, i can subset it like this. in the example a time series is filtered by whether the factor "experiment_day" is greater than 4
project_data.post_day4 <- prune_samples(sample_data(project_data)$experiment_day > 4, project_data)

## OPTION 1: DESEQ RUN USING LRT MODEL ##
dds.var3 <- phyloseq_to_deseq2(project_data, design = ~ variable1 + variable2 + variable3)
dds.var3 <- DESeq(dds.var3, test = "LRT", fitType = "parametric", reduced=~variable1 + variable2) 
resultsNames(dds.var3) #check results names for the name of the contrast you'd like to examine
res.var3 <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", name="variable3_level2_vs_level1", altHypothesis = "greaterAbs") #do not use "contrast" with the LRT test, it may give you the wrong pvalues #whatever contrast is calcualted, it's done by default with the first category against the others, in alphabetical order, or in the order they appear in the factor levels
# view a summary of the results table with a padj value < 0.01
summary(res.var3, alpha = alpha)

## OPTION 2: DESEQ RUN USING WALD TEST ##
dds.var3 <- phyloseq_to_deseq2(project_data, design = ~ variable1 + variable2 + variable3)
dds.var3 <- DESeq(dds.var3, test = "Wald", fitType = "parametric")
res.var3 <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", contrast=c("variable3", "level2", "level1"), altHypothesis = "greaterAbs") # here we can use "contrast" to select the comparison we are interested in
# view a summary of the results table with a padj value < 0.01
summary(res.var3, alpha = alpha)

#for details on what the difference between the LRT and WALD tests are, see the DeSeq vignette

#filtering the results table #
# reorder the results table by adjusted p-value and remove any "NA" entries
res_p_ordered <- res.var3[order(res.var3$padj, na.last = NA), ]
# filter out any results which have a adjust p-value less than alpha (0.01)
res_p_ordered_filt <- res_p_ordered[which(res_p_ordered$padj < alpha), ]
# show results preserving taxa table
res_p_ordered_wtaxa.var3 <- cbind(as(res_p_ordered_filt, "data.frame"), as(tax_table(project_data)[rownames(res_p_ordered_filt), ], "matrix"))

#### deseq fold change plot ####
#the following code will make a plot that shows you log2 fold change of OTUs colored by their taxonomic levels (you can set what rank you'd like below)
#plot colorspace prep
numcol <- length(unique(res_p_ordered_wtaxa.var3$Rank3))
set.seed(15)
newpal <- qualpal(n = numcol, colorspace = "pretty") # use qualpalr colour palettes for easily distinguishing taxa
# Phylum order
x = tapply(res_p_ordered_wtaxa.var3$log2FoldChange, res_p_ordered_wtaxa.var3$Rank3, function(x) max(x))
x = sort(x, TRUE)
res_p_ordered_wtaxa.var3$Rank3 = factor(as.character(res_p_ordered_wtaxa.var3$Rank3), levels=names(x))
# Genus order
x = tapply(res_p_ordered_wtaxa.var3$log2FoldChange, res_p_ordered_wtaxa.var3$Rank6, function(x) max(x))
x = sort(x, TRUE)
res_p_ordered_wtaxa.var3$Rank6 = factor(as.character(res_p_ordered_wtaxa.var3$Rank6), levels=names(x))
#make the actual plot
pdf("differentially_present_taxa.variable3.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 10, # define plot width and height. completely up to user.
    height = 8)
ggplot(res_p_ordered_wtaxa.var3, aes(x=Rank6, y=log2FoldChange, color=Rank3)) + geom_point(size=6) + 
  scale_color_manual(values = newpal$hex) + #NOTE: to preserve colors between plots, you must use the line below
  scale_color_manual(values = c("clade1" = "#333BFF",  "clade2"="#BBB3FF",  "clade3" = "#CC6600" ,  "etc."="#9633FF")) #NOTE: use this to manually assign color to each clade, or use it in multiple plots to make colors uniform across many plots
  ggtitle("Differentially Abundant Taxa WRT variable3, My Experiment Name") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
dev.off()


#### save results table in plain text format ####
write.table(res_p_ordered_wtaxa.var3, "differentially_present_taxa.variable3.txt", quote=F, row.names=T, col.names=T, sep="\t")

