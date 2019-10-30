#### basic analysis of processed amplicon data with DeSeq2####
#author: Evan Morien
#last modified: October 30th, 2019

#NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT IF YOU WISH TO MODIFY IT PERMANENTLY
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

### ensure you have loaded into R phyloseq object(s) which contain your dataset, including taxonomy strings, sequence table, metadata, and phylogenetic tree (if applicable)####
#instructions for creating a phyloseq object are in our lab_scripts/R_scripts github folder, in the script "load_data_into_phyloseq_object.R"

####DESeq2 analysis preparation####
#OPTIONAL/IF NEEDED:assign variables as either factors or numeric
#VERY IMPORTANT: DESEQ is designed to work for factors that have two levels. if your factor that you want to test has more than two levels, you will need to select two of the levels to do the comparison, or deseq will select them for you.
#this can be necesssary for certain types of data, and helpful if you want to control the order of factors in a deseq contrast below
sample_data(project_data)$variable1 <- factor(sample_data(project_data)$variable1) #not controling levels
sample_data(project_data)$variable3 <- factor(sample_data(project_data)$variable3, levels=c("level1", "level2")) #controling level order here
#IMPORTANT: factors need to be watched carefully, recoding variables as factors can re-arrange data w/in the data frame if you have previously ordered it in a specific way

####DESeq: taxa correlated with a binary factor ("variable3"), controlling for two other factors ("variable1" and "variable2") ####
#parameters up front
alpha <- 0.01 #your significance threshold for MULTIPLE TEST CORRECTED pvals
##OPTIONAL: SUBSET DATA BEFORE TEST
#if i want to test only a subset of my data, i can subset it like this. in the example a time series is filtered by whether the factor "experiment_day" is greater than 4
project_data.numeric_variable.greater_than_4 <- prune_samples(sample_data(project_data)$numeric_variable > 4, project_data)
project_data.factor.only_level1 <- prune_samples(sample_data(project_data)$factor == "level1", project_data)

## OPTION 1: DESEQ RUN USING LRT MODEL ##
dds.var3 <- phyloseq_to_deseq2(project_data, design = ~ variable1 + variable2 + variable3)
dds.var3 <- DESeq(dds.var3, test = "LRT", fitType = "parametric", reduced=~variable1 + variable2) 
resultsNames(dds.var3) #check results names for the name of the contrast you'd like to examine
res.var3 <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", name="variable3_level2_vs_level1", altHypothesis = "greaterAbs") #here "level2" is experimental and "level1" is control or baseline #do not use "contrast" with the LRT test, it may give you the wrong pvalues #whatever contrast is calcualted, it's done by default with the first category against the others, in alphabetical order, or in the order they appear in the factor levels
# view a summary of the results table with a padj value < 0.01
summary(res.var3, alpha = alpha)

## OPTION 2: DESEQ RUN USING WALD TEST ##
dds.var3 <- phyloseq_to_deseq2(project_data, design = ~ variable1 + variable2 + variable3)
dds.var3 <- DESeq(dds.var3, test = "Wald", fitType = "parametric")
res.var3 <- results(dds.var3, cooksCutoff = FALSE, alpha = alpha, pAdjustMethod = "BH", contrast=c("variable3", "level2", "level1"), altHypothesis = "greaterAbs") #here "level2" is experimental and "level1" is control or baseline #here we can use "contrast" to specify the comparison we are interested in
# view a summary of the results table with a padj value < 0.01
summary(res.var3, alpha = alpha)

#for details on what the difference between the LRT and WALD tests are, see the DeSeq vignette

#IMPORTANT: you may get an error like this when trying to run DESeq:
    #estimating size factors
    #Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
    #every gene contains at least one zero, cannot compute log geometric means
#It's then up to you to figure out why you are unable to estimate the size factors for this data, and adjust your factors or apply an alternative method to estimate the size factors.
#You can read more about this online, especially on websites like biostars or on the bioconductor support website.
#Michael Love, one of the main developers of DESeq2 has very extensively and thoroughly answered questions on this error.

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

