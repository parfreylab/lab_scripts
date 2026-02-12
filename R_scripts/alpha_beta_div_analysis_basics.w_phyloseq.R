####basic procedure for completing alpha and beta diversity analyses with a phyloseq object####
#author: Evan Morien
#last modified: Feb 5th, 2026

#IMPORTANT NOTE:
# for both alpha and beta diversity analyses, data should be rarefied.
# although there is an argument to be made against rarefying for beta diversity analyses, we find that with many of our lab experiments, the bias from differentially successful sequencing (i.e. samples with many more reads than others, or many fewer) will prevent the ordination from revealing meaningful biological signals in the data.

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
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)
#load the vegan library
library(vegan)
#load ranacapa library for making rarefaction plots
library(ranacapa)

#### environment settings ####
#set working directory
setwd("/path/to/working/directory/")

####IMPORTANT: first, make sure you have data in phyloseq format####
#if you don't have a phyloseq object ready, please see the script on the lab github that details loading data from various sources into phyloseq (load_data_into_phyloseq_object.R)
#in the examples below, your complete, filtered phyloseq object is called "project_data"

#load phyloseq object from .RDS
project_data <- readRDS("my_phyloseq_object.RDS")

#### plot rarefaction curves ####
plot(sort(sample_sums(project_data))) #looking at sample read counts
sort(sample_sums(project_data))
summary(sample_sums(project_data))

#making a rarefaction plot with ggrare() function
p <- ggrare(project_data, step = 1000, color = "FACTOR_1", se = FALSE) #coloring is optional, factor must be selected
# OPTIONAL: facet_wrap the plot to see differences according to different sample sites, types, etc. FACTOR_1, FACTOR_2, etc represent columns in your metadata, viewable through the sample_data() function from phyloseq package
p + facet_wrap(~FACTOR_1 + FACTOR_2)
#you can use the plot above to judge a rough cutoff for rarefaction. it is also possible to do this with QIIME's alpha rarefaction script if you have a .biom file

#you can use the "which" command like the example below to see which samples you'll lose for a given cutoff value
which(sample_sums(project_data) < 2000)

#### rarefy data ####
set.seed(24) #you must set a numerical seed like this for reproducibility, but keep in mind that if your diversity results differ significantly after changing the seed, then there may be issues with your data. there is no right number to choose, you can use any whole number, just record it for reproducibiity purposes
project_data.rarefied <- rarefy_even_depth(project_data, sample.size = min(sample_sums(project_data)), replace=FALSE)

#### Create Colour Palettes ####
# 1. identify the number of colors you need from the factor you want to plot by
numcol <- length(unique(sample_data(project_data)$FACTOR_1)) #EXAMPLE ONLY, adjust per object being plotted
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(13)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")
# 4. If you want to color only a few things in the plot, try using a more colourblind-friendly palette:
cbPalette <- c("#E69F00", "#56B4E9", "#000000", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#FFFF00", "#999999", "#FF00FF", "#F0E442", "#FFFFFF", "#00FFFF")

#### basic alpha div plot ####
#### step1: calculate alpha diversity and add it as a column in the metadata ####
project_data.chao1 <- estimate_richness(project_data.rarefied, split = TRUE, measures = c("Chao1")) #estimate richness
sample_data(project_data.rarefied)$chao1 <- project_data.chao1$Chao1 #add to metadata (the rows are in the same order already)
sample_data(project_data.rarefied)$chao1 <- as.numeric(sample_data(project_data.rarefied)$chao1)

#### step 2: use calculated alpha diversity to make basic plot with ggplot ####
#this plot lets you customize things a bit more than the plot_richness function, if desired
p <- ggplot(sample_data(project_data.rarefied), aes(x=FACTOR_1, y=chao1, color=FACTOR_3))
p + geom_boxplot() + 
  facet_grid(~ FACTOR_2, drop=TRUE, scales="free", space="free") + #drop, scales, and space are used to ensure that each boxplot bar is equally sized, and that empty levels created by the facet grid are dropped from the final plot
  labs(title="Alpha Diversity, Factor 1 ~ Factor 2 + Factor 3, Chao1", x="Factor1 ~ Factor2", y="Chao1") + 
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### beta diversity (NMDS, PCoA, etc.) ####
#do ordinations #be sure to use rarefied data for beta diversity analyses
set.seed(24) #same deal here as above w/r/t choosing a numeric seed for the ordination. specific number unimportant, always record for reproducibility
NMDS.bray <- ordinate(
  physeq = project_data.rarefied, 
  method = "NMDS", 
  distance = "bray"
) # you can choose different methods and distance metrics, see the ordinate help page for details. this function works with "phyloseq" class objects.

#### making beta div plots ####
#we get more plotting control if we don't use the phyloseq plotting functions for ordination plots, and instead add the results of the ordination to our existing metadata
NMDS <- as.data.frame(sample_data(project_data.rarefied))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$FACTOR_1, NMDS$FACTOR_2),] #if you do create a sorted plotting object, remember to modify the plots below so that you are using it as the input, not the unsorted table

#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
pdf("NMDS.expt_name.FACTOR_1_color.FACTOR_2_shape.pdf"
    , width = 16 # Default is 7
    , height = 9 # Change to 10; make it taller
)
p <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, shape = FACTOR_2, color = FACTOR_1)) # change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=4) + scale_shape_manual(values=1:nlevels(NMDS$FACTOR_2)) +
  labs(title="NMDS Factor1 & Factor2") + 
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#NOTE: beta div ordination plots like the one above are endlessly customizable, if you want to do something like draw lines between points to show how a time series connects, or make the points different sizes based on a factor, or any other advanced plotting technique, there will be examples online describing how to do this with ggplot

####for statistical tests on beta diversity data (permanova, beta dispersion test), please see the lab script permanova.beta_dispersion_test.basics.R ####
